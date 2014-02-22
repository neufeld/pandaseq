/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
#define _POSIX_C_SOURCE 2
#include "config.h"
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#if HAVE_UNAME_SYSCALL
#        include<sys/utsname.h>
#endif
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"
#include "module.h"
#ifdef HAVE_PTHREAD
#        include"pandaseq-mux.h"
#endif

#define MAX_OPT_LIST 53
#define MAX_MODULES 100

int compare_general(
	const void *key,
	const void *item) {
	const panda_tweak_general **tweak = (const panda_tweak_general **) item;
	return *((int *) key) - (*tweak)->flag;
}

int compare_assembler(
	const void *key,
	const void *item) {
	const panda_tweak_assembler **tweak = (const panda_tweak_assembler **) item;
	return *((int *) key) - (*tweak)->flag;
}

#define CLEANUP() for (it = 0; it < *options_used; it++) if(options[it].arg != NULL) free(options[it].arg); *options_used = 0;
#define BSEARCH(item, type) bsearch(item, PANDACONCAT(type, _args), PANDACONCAT(type, _args_length), sizeof(PANDACONCAT(panda_tweak_, type) *), PANDACONCAT(compare_, type))

bool panda_dispatch_args(
	char *const *args,
	int args_length,
	const panda_tweak_assembler *const *const assembler_args,
	size_t assembler_args_length,
	const panda_tweak_general *const *const general_args,
	size_t general_args_length,
	PandaTweakGeneral tweak,
	void *user_data,
	panda_tweak_assembler_opt * options,
	size_t options_length,
	size_t *options_used,
	int *args_unused) {

	int c;
	bool help = false;
	size_t it;
	char optlist[MAX_OPT_LIST + 1];
	char seen_options[MAX_OPT_LIST + 1];
	size_t seen_options_length = 0;
	const panda_tweak_general **general_tweak;
	const panda_tweak_assembler **assembler_tweak;
	size_t opt_it;

#ifdef __GLIBC__
	optind = 0;
#else
	optind = 1;
#endif
	opt_it = 0;
	*options_used = 0;
	MAYBE(args_unused) = 0;

	c = '\0';
	for (it = 0; it < general_args_length; it++) {
		if (general_args[it]->flag <= c) {
			fprintf(stderr, "Internal error: general arguments in a bad order.\n");
			abort();
		}
		c = optlist[opt_it++] = general_args[it]->flag;
		if (general_args[it]->takes_argument != NULL)
			optlist[opt_it++] = ':';
		if (opt_it >= MAX_OPT_LIST) {
			fprintf(stderr, "Internal error: too many options.\n");
			abort();
		}
	}
	c = '\0';
	for (it = 0; it < assembler_args_length; it++) {
		if (assembler_args[it]->flag <= c) {
			fprintf(stderr, "Internal error: assembler arguments in a bad order.\n");
			abort();
		}
		c = optlist[opt_it++] = assembler_args[it]->flag;
		if (assembler_args[it]->takes_argument != NULL)
			optlist[opt_it++] = ':';
		if (opt_it >= MAX_OPT_LIST) {
			fprintf(stderr, "Internal error: too many options.\n");
			abort();
		}
	}
	optlist[opt_it] = '\0';

	/* Process command line arguments. */
	while ((c = getopt(args_length, args, optlist)) != -1) {
		if (c == '?') {
			if (strchr(optlist, optopt) != NULL) {
				if ((general_tweak = BSEARCH(&optopt, general)) != NULL) {
					fprintf(stderr, "Option -%c requires an argument %s.\n", optopt, (*general_tweak)->takes_argument);
				} else if ((assembler_tweak = BSEARCH(&optopt, assembler)) != NULL) {
					fprintf(stderr, "Option -%c requires an argument %s.\n", optopt, (*assembler_tweak)->takes_argument);
				} else {
					fprintf(stderr, "Unhandled command line argument -%c requires an argument. This is a bug.\n", (int) optopt);
				}
			} else if (isprint(optopt)) {
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			} else {
				fprintf(stderr, "Unknown option character `\\x%x'.\n", (unsigned int) optopt);
			}
			CLEANUP();
			return false;
		} else {
			if ((general_tweak = BSEARCH(&c, general)) != NULL) {
				if ((*general_tweak)->takes_argument != NULL || !(*general_tweak)->repeatable) {
					if (memchr(seen_options, c, seen_options_length) != NULL) {
						fprintf(stderr, "Command line argument -%c may not be repeated.\n", (int) c);
						CLEANUP();
						return false;
					}
					if (seen_options_length == MAX_OPT_LIST) {
						fprintf(stderr, "Too many command line arguments.\n");
						CLEANUP();
						return false;
					}
					seen_options[seen_options_length++] = c;
				}

				if (!tweak(user_data, c, (*general_tweak)->takes_argument != NULL ? optarg : NULL)) {
					CLEANUP();
					return false;
				}
			} else if ((assembler_tweak = BSEARCH(&c, assembler)) != NULL) {
				if ((*assembler_tweak)->takes_argument == NULL || !(*assembler_tweak)->repeatable) {
					if (memchr(seen_options, c, seen_options_length) != NULL) {
						fprintf(stderr, "Command line argument -%c may not be repeated.\n", (int) c);
						CLEANUP();
						return false;
					}
					if (seen_options_length == MAX_OPT_LIST) {
						fprintf(stderr, "Too many command line arguments.\n");
						CLEANUP();
						return false;
					}
					seen_options[seen_options_length++] = c;
				}

				it = (*options_used)++;
				if (it >= options_length) {
					fprintf(stderr, "Too many command line arguments.\n");
					CLEANUP();
					return false;
				}
				options[it].tweak = *assembler_tweak;
				if ((*assembler_tweak)->takes_argument != NULL) {
					options[it].arg = malloc(strlen(optarg) + 1);
					memcpy(options[it].arg, optarg, strlen(optarg) + 1);
				} else {
					options[it].arg = NULL;
				}
			} else {
				fprintf(stderr, "Unhandled command line argument -%c. This is a bug.\n", (int) c);
				CLEANUP();
				return false;
			}
		}
	}
	MAYBE(args_unused) = optind;
	return true;
}

#undef CLEANUP

void panda_args_help(
	const char *binary_name,
	const panda_tweak_assembler *const *const assembler_args,
	size_t assembler_args_length,
	const panda_tweak_general *const *const general_args,
	size_t general_args_length) {
	size_t it;
	size_t general_it = 0;
	size_t assembler_it = 0;

	fprintf(stderr, "%s <%s>\nUsage: %s", PACKAGE_STRING, PACKAGE_BUGREPORT, binary_name);
	for (general_it = 0; general_it < general_args_length; general_it++) {
		if (general_args[general_it]->takes_argument != NULL && !general_args[general_it]->optional) {
			fprintf(stderr, " -%c %s", (int) general_args[general_it]->flag, general_args[general_it]->takes_argument);
		}
	}
	general_it = 0;
	while (general_it < general_args_length || assembler_it < assembler_args_length) {
		if (general_it < general_args_length && assembler_it < assembler_args_length && general_args[general_it]->flag == assembler_args[assembler_it]->flag) {
			assembler_it++;
			continue;
		}
		if (general_it < general_args_length && (assembler_it == assembler_args_length || general_args[general_it]->flag < assembler_args[assembler_it]->flag)) {

			if (general_args[general_it]->takes_argument != NULL) {
				if (general_args[general_it]->optional) {
					fprintf(stderr, " [-%c %s]", (int) general_args[general_it]->flag, general_args[general_it]->takes_argument);
				}
			} else {
				fprintf(stderr, " [-%c]", (int) general_args[general_it]->flag);
			}
			general_it++;
		} else {
			if (assembler_args[assembler_it]->takes_argument != NULL) {
				fprintf(stderr, " [-%c %s]", (int) assembler_args[assembler_it]->flag, assembler_args[assembler_it]->takes_argument);
			} else {
				fprintf(stderr, " [-%c]", (int) assembler_args[assembler_it]->flag);
			}
			assembler_it++;
		}
	}
	fprintf(stderr, "\n");
	general_it = 0;
	assembler_it = 0;
	while (general_it < general_args_length || assembler_it < assembler_args_length) {
		if (general_it < general_args_length && assembler_it < assembler_args_length && general_args[general_it]->flag == assembler_args[assembler_it]->flag) {
			assembler_it++;
			continue;
		}
		if (general_it < general_args_length && (assembler_it == assembler_args_length || general_args[general_it]->flag < assembler_args[assembler_it]->flag)) {
			if (general_args[general_it]->takes_argument != NULL) {
				fprintf(stderr, "\t-%c %s\t%s%s\n", (int) general_args[general_it]->flag, general_args[general_it]->takes_argument, general_args[general_it]->help, general_args[general_it]->repeatable ? " (May be repeated.)" : "");
			} else {
				fprintf(stderr, "\t-%c\t%s\n", (int) general_args[general_it]->flag, general_args[general_it]->help);
			}
			general_it++;
		} else {
			if (assembler_args[assembler_it]->takes_argument != NULL) {
				fprintf(stderr, "\t-%c %s\t%s%s\n", (int) assembler_args[assembler_it]->flag, assembler_args[assembler_it]->takes_argument, assembler_args[assembler_it]->help, assembler_args[assembler_it]->repeatable ? " (May be repeated.)" : "");
			} else {
				fprintf(stderr, "\t-%c\t%s\n", (int) assembler_args[assembler_it]->flag, assembler_args[assembler_it]->help);
			}
			assembler_it++;
		}
	}
	fprintf(stderr, "Known algorithms are %s", panda_algorithms[0]->name);
	for (it = 1; it < panda_algorithms_length; it++) {
		fprintf(stderr, ", %s", panda_algorithms[it]->name);
	}
	fprintf(stderr, ".\n");
	module_show_all();
}

struct data {
	bool fastq;
	bool help;
	size_t num_kmers;
	bool version;
	PandaWriter writer_err;
	PandaWriter writer_out;
#ifdef HAVE_PTHREAD
	int threads;
#endif
	PandaTweakGeneral general;
	void *general_data;
};

static const panda_tweak_general logging = {.flag = 'd',.optional = true,.takes_argument = "flags",.help = "Control the logging messages. Capital to enable; small to disable.\n\t\t(R)econstruction detail.\n\t\tSequence (b)uilding information.\n\t\t(F)ile processing.\n\t\t(k)-mer table construction.\n\t\tShow every (m)ismatch.\n\t\tOptional (s)tatistics." };
static const panda_tweak_general kmers = {.flag = 'k',.optional = true,.takes_argument = "kmers",.help = "The number of k-mers in the table." };
static const panda_tweak_general fastq = {.flag = 'F',.optional = true,.takes_argument = NULL,.help = "Output FASTQ instead of FASTA." };
static const panda_tweak_general logfile = {.flag = 'g',.optional = true,.takes_argument = "log.txt",.help = "Output log to a text file." };
static const panda_tweak_general logfile_bz = {.flag = 'G',.optional = true,.takes_argument = "log.txt.bz2",.help = "Output log to a BZip2-compressed text file." };

#		ifdef HAVE_PTHREAD
static const panda_tweak_general threads = {.flag = 'T',.optional = true,.takes_argument = "threads",.help = "Run with a number of parallel threads." };
#		endif
static const panda_tweak_general outputfile = {.flag = 'w',.optional = true,.takes_argument = "output.fasta",.help = "Output seqences to a FASTA (or FASTQ) file." };
static const panda_tweak_general outputfile_bz = {.flag = 'W',.optional = true,.takes_argument = "output.fasta.bz2",.help = "Output seqences to a BZip2-compressed FASTA (or FASTQ) file." };
static const panda_tweak_general help = {.flag = 'h',.optional = true,.takes_argument = NULL,.help = "Show this delightful nonsense." };
static const panda_tweak_general version = {.flag = 'v',.optional = true,.takes_argument = NULL,.help = "Show version and exit." };

static const panda_tweak_general *common_args[] = {
	&fastq,
	&help,
	&kmers,
	&logfile,
	&logfile_bz,
	&logging,
	&outputfile,
	&outputfile_bz,
#		ifdef HAVE_PTHREAD
	&threads,
#		endif
	&version
};

static bool common_tweak_general(
	void *user_data,
	char flag,
	const char *argument) {

	struct data *data = (struct data *) user_data;
	size_t it;
	long int value;

	switch (flag) {
	case 'd':
		for (it = 0; it < strlen(argument); it++) {
			PandaDebug flag = 0;
			switch (tolower(argument[it])) {
			case 'b':
				flag = PANDA_DEBUG_BUILD;
				break;
			case 'f':
				flag = PANDA_DEBUG_FILE;
				break;
			case 's':
				flag = PANDA_DEBUG_STAT;
				break;
			case 'k':
				flag = PANDA_DEBUG_KMER;
				break;
			case 'r':
				flag = PANDA_DEBUG_RECON;
				break;
			case 'm':
				flag = PANDA_DEBUG_MISMATCH;
				break;
			default:
				fprintf(stderr, "Ignoring unknown debug flag `%c'.\n", (int) argument[it]);
				continue;
			}
			if (islower(argument[it])) {
				panda_debug_flags &= ~flag;
			} else {
				panda_debug_flags |= flag;
			}
		}
		return true;
	case 'F':
		data->fastq = true;
		return true;
	case 'g':
	case 'G':
		panda_writer_unref(data->writer_err);
		data->writer_err = panda_writer_open_file(argument, isupper(flag));
		if (data->writer_err == NULL) {
			perror(argument);
			return false;
		}
		return true;
	case 'h':
		data->help = true;
		return true;
	case 'k':
		errno = 0;
		value = strtol(argument, NULL, 10);
		if (errno != 0 || value < 0 || value > PANDA_MAX_LEN) {
			fprintf(stderr, "Bad k-mer list length.\n");
			return false;
		}
		data->num_kmers = (size_t) value;
		return true;
#ifdef HAVE_PTHREAD
	case 'T':
		errno = 0;
		data->threads = (size_t) strtol(argument, NULL, 10);
		if (errno != 0 || data->threads < 1) {
			fprintf(stderr, "Bad number of threads.\n");
			return false;
		}
		return true;
#endif
	case 'w':
	case 'W':
		panda_writer_unref(data->writer_out);
		data->writer_out = panda_writer_open_file(argument, isupper(flag));
		if (data->writer_out == NULL) {
			perror(argument);
			return false;
		}
		return true;
	case 'v':
		data->version = true;
		return true;
	default:
		return data->general(data->general_data, flag, argument);
	}
}

#define BASE_CLEANUP() for (it = 0; it < options_used; it++) if(options[it].arg != NULL) free(options[it].arg); DESTROY_STACK(next); DESTROY_STACK(fail); panda_assembler_unref(assembler); panda_log_proxy_unref(logger); panda_writer_unref(data.writer_out); panda_writer_unref(data.writer_err); free(combined_general_args)
#ifdef HAVE_PTHREAD
#        define CLEANUP() BASE_CLEANUP(); panda_mux_unref(mux)
#else
#        define CLEANUP() BASE_CLEANUP()
#endif

bool panda_parse_args(
	char *const *args,
	int args_length,
	const panda_tweak_assembler *const *const assembler_args,
	size_t assembler_args_length,
	const panda_tweak_general *const *const general_args,
	size_t general_args_length,
	PandaTweakGeneral tweak,
	PandaOpener opener,
	PandaSetup assembler_setup,
	void *user_data,
	PandaAssembler *out_assembler,
	PandaMux *out_mux,
	int *out_threads,
	PandaOutputSeq * output,
	void **output_data,
	PandaDestroy *output_destroy) {

	PandaAssembler assembler = NULL;
	const panda_tweak_general **combined_general_args = NULL;
	size_t combined_general_args_length = 0;
	struct data data;
	PandaLogProxy logger = NULL;
	size_t it;
#ifdef HAVE_PTHREAD
	PandaMux mux = NULL;
#endif
	MANAGED_STACK(PandaFailAlign,
		fail);
	MANAGED_STACK(PandaNextSeq,
		next);
	panda_tweak_assembler_opt options[50];
	size_t options_used;
	int args_unused;
#if HAVE_UNAME_SYSCALL
	struct utsname uname_info;
#endif

	data.general = tweak;
	data.general_data = user_data;
	data.fastq = false;
	data.help = false;
	data.num_kmers = PANDA_DEFAULT_NUM_KMERS;
	data.version = false;
	data.writer_out = panda_writer_new_stdout();
	data.writer_err = panda_writer_new_stderr();
#ifdef HAVE_PTHREAD
	data.threads = panda_get_default_worker_threads();
#endif

	MAYBE(out_mux) = NULL;
	MAYBE(out_assembler) = NULL;
#ifdef HAVE_PTHREAD
	MAYBE(out_threads) = data.threads;
#else
	MAYBE(out_threads) = 1;
#endif
	MAYBE(output) = NULL;
	MAYBE(output_data) = NULL;
	MAYBE(output_destroy) = NULL;

	/* Process command line arguments. */
	panda_tweak_general_append(&combined_general_args, &combined_general_args_length, general_args, general_args_length);
	panda_tweak_general_append(&combined_general_args, &combined_general_args_length, common_args, sizeof(common_args) / sizeof(common_args[0]));
	if (!panda_dispatch_args(args, args_length, assembler_args, assembler_args_length, combined_general_args, combined_general_args_length, common_tweak_general, &data, options, sizeof(options) / sizeof(options[0]), &options_used, &args_unused)) {
		CLEANUP();
		return false;
	}

	if (args_length - args_unused > 1) {
		fprintf(stderr, "Ignoring extra arguments passed.\n");
	}

	logger = panda_log_proxy_new(data.writer_err);
	if (data.version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		CLEANUP();
		return false;
	}
	panda_writer_set_slave(data.writer_out, data.writer_err);

	if (data.help) {
		panda_args_help(args[0], assembler_args, assembler_args_length, combined_general_args, combined_general_args_length);
		CLEANUP();
		return false;
	}
	if ((next = opener(user_data, logger, &fail, &fail_data, &fail_destroy, &next_data, &next_destroy)) == NULL) {
		panda_writer_append(data.writer_err, "Too confused to continue.\nTry -h for help.\n");
		panda_writer_commit(data.writer_err);
		CLEANUP();
		return false;
	}
#ifdef HAVE_PTHREAD
	next = panda_create_async_reader(next, next_data, next_destroy, data.threads, &next_data, &next_destroy);
#endif
	panda_log_proxy_write_str(logger, "INFO\tVER\t" PACKAGE_STRING " <" PACKAGE_BUGREPORT ">");

#if HAVE_UNAME_SYSCALL
	if (uname(&uname_info) == 0) {
		panda_writer_append(data.writer_err, "INFO\tUNAME\t%s %s %s %s\n", uname_info.sysname, uname_info.release, uname_info.version, uname_info.machine);
		panda_writer_commit(data.writer_err);
	} else {
		panda_log_proxy_perror(logger, "uname");
	}
#endif

	for (it = 0; it < args_length; it++) {
		char buf[2048];
		if (snprintf(buf, sizeof(buf), "ARG[%d]\t%s", (int) it, args[it]) < sizeof(buf)) {
			panda_log_proxy_write_str(logger, buf);
		}
	}

#ifdef HAVE_PTHREAD
	mux = panda_mux_new(next, next_data, next_destroy, logger);
	if (mux == NULL) {
		panda_log_proxy_write_str(logger, "ERR\tLIB\tCould not create multiplexer.\n");
		CLEANUP();
		return false;
	}
	assembler = panda_mux_create_assembler_kmer(mux, data.num_kmers);
#else
	assembler = panda_assembler_new_kmer(next, next_data, next_destroy, logger, data.num_kmers);
#endif
	next = NULL;
	next_data = NULL;
	next_destroy = NULL;
	if (assembler == NULL) {
		panda_log_proxy_write_str(logger, "ERR\tLIB\tCould not create assembler.\n");
		CLEANUP();
		return false;
	}
	for (it = 0; it < options_used; it++) {
		char *arg = options[it].arg;
		options[it].arg = NULL;
		if (!(options[it].tweak->setup) (assembler, options[it].tweak->flag, arg)) {
			CLEANUP();
			return false;
		}
	}
	if (assembler_setup != NULL && !assembler_setup(user_data, assembler)) {
		CLEANUP();
		return false;
	}
	if (fail != NULL) {
#if HAVE_PTHREAD
		panda_mux_set_fail_alignment(mux, fail, fail_data, fail_destroy);
#else
		panda_assembler_set_fail_alignment(assembler, fail, fail_data, fail_destroy);
#endif
		fail = NULL;
		fail_data = NULL;
		fail_destroy = NULL;
	}
	if (out_mux) {
#if HAVE_PTHREAD
		*out_mux = mux;
		mux = NULL;
#else
		*out_mux = NULL;
#endif
	}

	if (out_assembler) {
		*out_assembler = assembler;
		assembler = NULL;
	}
#if HAVE_PTHREAD
	MAYBE(out_threads) = data.threads;
#else
	MAYBE(out_threads) = 1;
#endif
	MAYBE(output) = (PandaOutputSeq) (data.fastq ? panda_output_fastq : panda_output_fasta);
	MAYBE(output_data) = data.writer_out;
	MAYBE(output_destroy) = (PandaDestroy) panda_writer_unref;
	data.writer_out = NULL;

	CLEANUP();
	return true;
}

bool panda_parse_key_values(
	const char *str,
	PandaKeyParsed key_parsed,
	void *key_parsed_data) {
	char key[100];
	size_t key_length = 0;
	char value[200];
	size_t value_length = 0;
	bool in_key = true;
	if (str == NULL)
		return true;

	for (; *str != '\0'; str++) {
		if (*str == ',' && !in_key) {
			value[value_length] = '\0';
			if (!key_parsed(key, value, key_parsed_data))
				return false;
			key_length = 0;
			value_length = 0;
			in_key = true;
		}
		if (*str == '=' && in_key) {
			key[key_length] = '\0';
			in_key = false;
		} else if (!in_key) {
			value[value_length++] = *str;
		} else if (isalnum(*str) || *str == '_' || *str == '-') {
			key[key_length++] = *str;
		} else {
			return false;
		}
		if (key_length == sizeof(key) || value_length == sizeof(value))
			return false;
	}
	if (in_key)
		return false;
	if (key_length > 0 && !in_key) {
		value[value_length] = '\0';
		return key_parsed(key, value, key_parsed_data);
	}
	return !in_key;
}
