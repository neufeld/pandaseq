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
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

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

#define CLEANUP() for (it = 0; it < modules_length; it++) panda_module_unref(modules[it]); for (it = 0; it < assembler_args_length; it++) if(opt_assembler_args[it] > (char*)1) free(opt_assembler_args[it]); DESTROY_STACK(next); DESTROY_STACK(fail); if (mux != NULL) panda_mux_unref(mux); if (assembler != NULL) panda_assembler_unref(assembler); panda_log_proxy_unref(logger)
#define BSEARCH(item, type) bsearch(item, PANDACONCAT(type, _args), PANDACONCAT(type, _args_length), sizeof(PANDACONCAT(panda_tweak_, type) *), PANDACONCAT(compare_, type))

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
	MANAGED_STACK(PandaNextSeq,
		next);
	MANAGED_STACK(PandaFailAlign,
		fail);
	PandaAssembler assembler = NULL;
	int c;
	bool fastq = false;
	bool help = false;
	size_t it;
	PandaLogProxy logger = panda_log_proxy_new_stderr();
	size_t modules_length = 0;
	PandaModule modules[MAX_MODULES];
	size_t num_kmers = PANDA_DEFAULT_NUM_KMERS;
	char optlist[MAX_OPT_LIST + 1];
	char *opt_assembler_args[assembler_args_length];
	const panda_tweak_general **general_tweak;
	const panda_tweak_assembler **assembler_tweak;
	size_t opt_it;
	bool version = false;
#ifdef HAVE_PTHREAD
	PandaMux mux = NULL;
	int threads = 1;
#endif

	MAYBE(out_mux) = NULL;
	MAYBE(out_assembler) = NULL;
	MAYBE(out_threads) = threads;
	MAYBE(output) = NULL;
	MAYBE(output_data) = NULL;
	MAYBE(output_destroy) = NULL;

	memset(&opt_assembler_args, 0, assembler_args_length * sizeof(char *));
	optlist[0] = '\0';
	strncat(optlist, "C:d:g:G:hk:vF", MAX_OPT_LIST);
#ifdef HAVE_PTHREAD
	strncat(optlist, "T:", MAX_OPT_LIST);
#endif

	opt_it = strlen(optlist);

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
		switch (c) {
		case 'C':
			if (modules_length == MAX_MODULES || (modules[modules_length] = panda_module_load(optarg)) == NULL) {
				if (modules_length == MAX_MODULES) {
					fprintf(stderr, "Too many modules.\n");
				}
				CLEANUP();
				return false;
			}
			modules_length++;
			break;
		case 'd':
			for (it = 0; it < strlen(optarg); it++) {
				PandaDebug flag = 0;
				switch (tolower(optarg[it])) {
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
					fprintf(stderr, "Ignoring unknown debug flag `%c'.\n", (int) optarg[it]);
					continue;
				}
				if (islower(optarg[it])) {
					panda_debug_flags &= ~flag;
				} else {
					panda_debug_flags |= flag;
				}
			}
			break;
		case 'F':
			fastq = true;
			break;
		case 'g':
		case 'G':
			panda_log_proxy_unref(logger);
			logger = panda_log_proxy_open_file(optarg, isupper(c));
			if (logger == NULL) {
				perror(optarg);
				CLEANUP();
				return false;
			}
			break;
		case 'h':
			help = true;
			break;
		case 'k':
			errno = 0;
			num_kmers = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || num_kmers < 0 || num_kmers > PANDA_MAX_LEN) {
				fprintf(stderr, "Bad k-mer list length.\n");
				CLEANUP();
				return false;
			}
			break;
#ifdef HAVE_PTHREAD
		case 'T':
			errno = 0;
			threads = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || threads < 1) {
				fprintf(stderr, "Bad number of threads.\n");
				CLEANUP();
				return false;
			}
			break;
#endif
		case 'v':
			version = true;
			break;
		case '?':
			if (strchr(optlist, optopt) != NULL) {
				if ((general_tweak = BSEARCH(&optopt, general)) != NULL) {
					fprintf(stderr, "Option -%c requires an argument %s.\n", optopt, (*general_tweak)->takes_argument);
				} else if ((assembler_tweak = BSEARCH(&optopt, assembler)) != NULL) {
					fprintf(stderr, "Option -%c requires an argument.\n", optopt, (*assembler_tweak)->takes_argument);
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
		default:
			if ((general_tweak = BSEARCH(&c, general)) != NULL) {
				if (!tweak(user_data, c, (*general_tweak)->takes_argument != NULL ? optarg : NULL)) {
					CLEANUP();
					return false;
				}
			} else if ((assembler_tweak = BSEARCH(&c, assembler)) != NULL) {
				it = assembler_tweak - assembler_args;
				if ((*assembler_tweak)->takes_argument != NULL) {
					opt_assembler_args[it] = malloc(strlen(optarg) + 1);
					memcpy(opt_assembler_args[it], optarg, strlen(optarg) + 1);
				} else {
					opt_assembler_args[it] = (void *) 1;
				}
			} else {
				fprintf(stderr, "Unhandled command line argument -%c. This is a bug.\n", (int) c);
				CLEANUP();
				return false;
			}
		}
	}

	if (version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		for (it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s %s\n", panda_module_get_name(modules[it]), panda_module_get_version(modules[it]));
		}
		CLEANUP();
		return false;
	}
	if (help || (next = opener(user_data, logger, &fail, &fail_data, &fail_destroy, &next_data, &next_destroy)) == NULL) {
		size_t general_it = 0;
		size_t assembler_it = 0;

		fprintf(stderr, "%s <%s>\nUsage: %s", PACKAGE_STRING, PACKAGE_BUGREPORT, args[0]);
		for (general_it = 0; general_it < general_args_length; general_it++) {
			if (general_args[general_it]->takes_argument != NULL && !general_args[general_it]->optional) {
				fprintf(stderr, " -%c %s", (int) general_args[general_it]->flag, general_args[general_it]->takes_argument);
			}
		}
		fprintf(stderr, " [-C module1 -C module2 ...]");
		fprintf(stderr, " [-d flags]");
		fprintf(stderr, " [-k kmers]");
		fprintf(stderr, " [-F]");
		fprintf(stderr, " [-g log.txt]");
		fprintf(stderr, " [-G log.txt.bz2]");
#		ifdef HAVE_PTHREAD
		fprintf(stderr, " [-T threads]");
#		endif
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
		fprintf(stderr, "\t-C module\tLoad a sequence validation module.\n");
		fprintf(stderr, "\t-d flags\tControl the logging messages. Capital to enable; small to disable.\n");
		fprintf(stderr, "\t\t(R)econstruction detail.\n");
		fprintf(stderr, "\t\tSequence (b)uilding information.\n");
		fprintf(stderr, "\t\t(F)ile processing.\n");
		fprintf(stderr, "\t\t(k)-mer table construction.\n");
		fprintf(stderr, "\t\tShow every (m)ismatch.\n");
		fprintf(stderr, "\t\tOptional (s)tatistics.\n");
		fprintf(stderr, "\t-k kmers\tThe number of k-mers in the table.\n");
		fprintf(stderr, "\t-F\tOutput FASTQ instead of FASTA.\n");
#		ifdef HAVE_PTHREAD
		fprintf(stderr, "\t-T thread\tRun with a number of parallel threads.\n");
#		endif
		general_it = 0;
		assembler_it = 0;
		while (general_it < general_args_length || assembler_it < assembler_args_length) {
			if (general_it < general_args_length && assembler_it < assembler_args_length && general_args[general_it]->flag == assembler_args[assembler_it]->flag) {
				assembler_it++;
				continue;
			}
			if (general_it < general_args_length && (assembler_it == assembler_args_length || general_args[general_it]->flag < assembler_args[assembler_it]->flag)) {
				if (general_args[general_it]->takes_argument != NULL) {
					fprintf(stderr, "\t-%c %s\t%s\n", (int) general_args[general_it]->flag, general_args[general_it]->takes_argument, general_args[general_it]->help);
				} else {
					fprintf(stderr, "\t-%c\t%s\n", (int) general_args[general_it]->flag, general_args[general_it]->help);
				}
				general_it++;
			} else {
				if (assembler_args[assembler_it]->takes_argument != NULL) {
					fprintf(stderr, "\t-%c %s\t%s\n", (int) assembler_args[assembler_it]->flag, assembler_args[assembler_it]->takes_argument, assembler_args[assembler_it]->help);
				} else {
					fprintf(stderr, "\t-%c\t%s\n", (int) assembler_args[assembler_it]->flag, assembler_args[assembler_it]->help);
				}
				assembler_it++;
			}
		}
		for (it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s(%s) %s\n\t%s\n", panda_module_get_name(modules[it]), panda_module_get_description(modules[it]), panda_module_get_version(modules[it]), panda_module_get_usage(modules[it]));
		}
		CLEANUP();
		return false;
	}
	panda_log_proxy_write_str(logger, "INFO\tVER\t" PACKAGE_STRING " <" PACKAGE_BUGREPORT ">\n");

#define BSIZE 2048
	for (it = 0; it < args_length; it++) {
		char buf[BSIZE];
		if (snprintf(buf, BSIZE, "ARG[%d]\t%s", (int) it, args[it]) < BSIZE) {
			panda_log_proxy_write_str(logger, buf);
		}
	}

#ifdef HAVE_PTHREAD
	mux = panda_mux_new(next, next_data, next_destroy, logger);
	if (mux == NULL) {
		fprintf(stderr, "ERR\tLIB\tCould not create multiplexer.\n");
		CLEANUP();
		return false;
	}
	assembler = panda_mux_create_assembler_kmer(mux, num_kmers);

#else
	assembler = panda_assembler_new_kmer(next, next_data, next_destroy, logger, num_kmers);
#endif
	next = NULL;
	next_data = NULL;
	next_destroy = NULL;
	if (assembler == NULL) {
		fprintf(stderr, "ERR\tLIB\tCould not create assembler.\n");
		CLEANUP();
		return false;
	}
	for (it = 0; it < assembler_args_length; it++) {
		char *arg = opt_assembler_args[it];
		opt_assembler_args[it] = NULL;
		if (!assembler_args[it]->setup(assembler, assembler_args[it]->flag, (assembler_args[it]->takes_argument != NULL) ? arg : NULL, arg != NULL)) {
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
	it = panda_assembler_add_modules(assembler, modules, modules_length);
	if (it != modules_length) {
		fprintf(stderr, "Problem with %s(%s) %s\n\t%s\n", panda_module_get_name(modules[it]), panda_module_get_description(modules[it]), panda_module_get_version(modules[it]), panda_module_get_usage(modules[it]));
		CLEANUP();
		panda_assembler_unref(assembler);
#if HAVE_PTHREAD
		panda_mux_unref(mux);
#endif
		return false;
	}
	if (out_mux) {
		*out_mux = mux;
		mux = NULL;
	}

	if (out_assembler) {
		*out_assembler = assembler;
		assembler = NULL;
	}
	MAYBE(out_threads) = threads;
	MAYBE(output) = (PandaOutputSeq) (fastq ? panda_output_fastq : panda_output_fasta);
	MAYBE(output_data) = stdout;
	MAYBE(output_destroy) = NULL;

	CLEANUP();
	return true;
}
