/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2012  Andre Masella

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
#include<ctype.h>
#include<errno.h>
#include<float.h>
#include<limits.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq.h"
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif

#define MAX_MODULES 100
bool fastq = false;
char *forward_primer = NULL;
ssize_t maxlen = 2 * PANDA_MAX_LEN + 1;
size_t minlen = 0;
bool no_n = false;
char *reverse_primer = NULL;
volatile bool some_seqs = false;
time_t starttime;
#ifdef HAVE_PTHREAD
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef HAVE_PTHREAD
#        	define STAT(str, val)	fprintf(stderr, "%p\tSTAT\t" str "\n", (void*) assembler, (val))
#else
#        	define STAT(str, val)	fprintf(stderr, "STAT\t" str "\n", (val))
#endif

static void
printtime(
	PandaAssembler assembler,
	long count,
	time_t starttime) {
	time_t now;
	(void) time(&now);
	STAT("TIME\t%s", ctime(&now));
	STAT("ELAPSED\t%d", (int) (now - starttime));
	STAT("READS\t%ld", count);
}

bool
short_check(
	const panda_result_seq *sequence,
	void *user_data) {
	return sequence->sequence_length >= minlen;
}

bool
long_check(
	const panda_result_seq *sequence,
	void *user_data) {
	size_t length = (size_t) user_data;
	return sequence->sequence_length <= maxlen;
}

static void *
do_assembly(
	PandaAssembler assembler) {
	long count;
	size_t it = 0;
	size_t max;
	const panda_result_seq *result;

	while ((result = panda_assembler_next(assembler)) != NULL) {
		count = panda_assembler_get_count(assembler);
#ifdef HAVE_PTHREAD
		pthread_mutex_lock(&output_mutex);
#endif
		if (count % 1000 == 0) {
			printtime(assembler, count, starttime);
		}
		if (fastq) {
			panda_output_fastq(result, stdout);
		} else {
			panda_output_fasta(result, stdout);
		}
#ifdef HAVE_PTHREAD
		pthread_mutex_unlock(&output_mutex);
#endif
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&output_mutex);
#endif
	count = panda_assembler_get_count(assembler);
	if (count > 0) {
		some_seqs = true;
	}
	printtime(assembler, count, starttime);
	if (forward_primer != NULL)
		STAT("NOFP\t%ld", panda_assembler_get_no_forward_primer_count(assembler));
	if (reverse_primer != NULL)
		STAT("NORP\t%ld", panda_assembler_get_no_reverse_primer_count(assembler));
	STAT("NOALGN\t%ld", panda_assembler_get_failed_alignment_count(assembler));
	STAT("LOWQ\t%ld", panda_assembler_get_low_quality_count(assembler));
	STAT("BADR\t%ld", panda_assembler_get_bad_read_count(assembler));
	STAT("SLOW\t%ld", panda_assembler_get_slow_count(assembler));
	if (no_n)
		STAT("DEGENERATE\t%ld", panda_assembler_get_degenerate_count(assembler));
	panda_assembler_module_stats(assembler);
	STAT("OK\t%ld", panda_assembler_get_ok_count(assembler));

#ifdef HAVE_PTHREAD
	fprintf(stderr, "%p\t", (void *) assembler);
#endif
	fprintf(stderr, "STAT\tOVERLAPS\t%ld", panda_assembler_get_overlap_count(assembler, it));
	max = panda_assembler_get_longest_overlap(assembler);
	for (it = 1; it < max; it++) {
		fprintf(stderr, " %ld", panda_assembler_get_overlap_count(assembler, it));
	}
	fprintf(stderr, "\n");
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&output_mutex);
#endif

	panda_assembler_unref(assembler);
	return NULL;
}

bool
set_primer(
	PandaAssembler assembler,
	void (*set_func) (PandaAssembler,
		panda_nt *,
		size_t),
	const char *str,
	panda_nt (*parse) (char)) {
	panda_nt buffer[PANDA_MAX_LEN];
	size_t it;
	for (it = 0; it < strlen(str); it++) {
		if ((buffer[it] = parse(str[it])) == '\0') {
			return false;
		}
	}
	set_func(assembler, buffer, strlen(str));
	return true;
}

int
main(
	int argc,
	char **argv) {
	PandaAssembler assembler;
	bool primers_after = false;
	bool bzip = false;
	int c;
	size_t foffset = 1;
	char *forward_filename = NULL;
	bool help = false;
	size_t it;
	int minoverlap = 1;
	size_t modules_length = 0;
	PandaModule modules[MAX_MODULES + 2];	// For LONG and SHORT
	const char *no_algn_file_name = NULL;
#ifdef HAVE_PTHREAD
	size_t num_kmers = PANDA_DEFAULT_NUM_KMERS;
#endif
	bool okay;
	const char *optlist = "6aBC:d:f:Fhjl:L:No:p:q:Q:r:t:u:v"
#ifdef HAVE_PTHREAD
		"k:T:"
#endif
		;
	PandaTagging policy = PANDA_TAG_PRESENT;
	double q = 0.36;
	int qualmin = 33;
	char *reverse_filename = NULL;
	size_t roffset = 1;
	double threshold;
	bool version = false;
#ifdef HAVE_PTHREAD
	PandaMux mux;
	int threads = 1;
	pthread_t *thread_list;
#endif
	(void) time(&starttime);
	threshold = 0.6;

	/* Process command line arguments. */
	while ((c = getopt(argc, argv, optlist)) != -1) {
		char *endptr;
		switch (c) {
		case '6':
			qualmin = 64;
			break;
		case 'a':
			primers_after = true;
			break;
		case 'B':
			policy = PANDA_TAG_OPTIONAL;
			break;
		case 'C':
			if (modules_length == MAX_MODULES || (modules[modules_length] = panda_module_load(optarg)) == NULL) {
				if (modules_length == MAX_MODULES) {
					fprintf(stderr, "Too many modules.\n");
				}
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
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
		case 'f':
			forward_filename = optarg;
			break;
		case 'F':
			fastq = 1;
			break;
		case 'h':
			help = true;
			break;
		case 'j':
			bzip = true;
			break;
		case 'l':
			errno = 0;
			minlen = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || minlen < 0 || minlen > 2 * PANDA_MAX_LEN) {
				fprintf(stderr, "Bad minimum length.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'L':
			errno = 0;
			maxlen = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || maxlen < 1 || maxlen > 2 * PANDA_MAX_LEN) {
				fprintf(stderr, "Bad maximum length.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'N':
			no_n = true;
			break;
		case 'o':
			errno = 0;
			minoverlap = strtol(optarg, NULL, 10);
			if (errno != 0 || minoverlap < 1 || minoverlap > PANDA_MAX_LEN) {
				fprintf(stderr, "Bad minimum overlap.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'p':
			errno = 0;
			foffset = strtol(optarg, &endptr, 10);
			if (*endptr != '\0') {
				forward_primer = optarg;
			} else if (errno != 0 || foffset < 1 || foffset > PANDA_MAX_LEN) {
				fprintf(stderr, "Bad forward primer length.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			} else {
				foffset++;
			}
			break;
		case 'q':
			errno = 0;
			roffset = strtol(optarg, &endptr, 10);
			if (*endptr != '\0') {
				reverse_primer = optarg;
			} else if (errno != 0 || roffset < 1 || roffset > PANDA_MAX_LEN) {
				fprintf(stderr, "Bad reverse primer length.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			} else {
				roffset++;
			}
			break;
		case 'Q':
			errno = 0;
			q = strtod(optarg, NULL);
			if (errno != 0 || q < 0 || q > 1) {
				fprintf(stderr, "Bad quality. It should be between 0 and 1.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'r':
			reverse_filename = optarg;
			break;
		case 't':
			errno = 0;
			threshold = strtod(optarg, NULL);
			if (errno != 0 || threshold < 0 || threshold > 1) {
				fprintf(stderr, "Bad threshold. It should be between 0 and 1.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'u':
			no_algn_file_name = optarg;
			break;
		case 'v':
			version = true;
			break;
#ifdef HAVE_PTHREAD
		case 'k':
			errno = 0;
			num_kmers = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || num_kmers < 0 || num_kmers > PANDA_MAX_LEN) {
				fprintf(stderr, "Bad k-mer list length.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
		case 'T':
			errno = 0;
			threads = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || threads < 1) {
				fprintf(stderr, "Bad number of threads.\n");
				for (it = 0; it < modules_length; it++)
					panda_module_unref(modules[it]);
				return 1;
			}
			break;
#endif
		case '?':
			if (strchr(optlist, optopt) != NULL) {
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			} else if (isprint(optopt)) {
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			} else {
				fprintf(stderr, "Unknown option character `\\x%x'.\n", (unsigned int) optopt);
			}
			return 1;
		default:
			abort();
		}
	}

	if (version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		for (it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s %s\n", panda_module_get_name(modules[it]), panda_module_get_version(modules[it]));
			panda_module_unref(modules[it]);
		}
		return 1;
	}
	if (forward_filename == NULL || reverse_filename == NULL || help) {
		fprintf(stderr, "%s <%s>\nUsage: %s ", PACKAGE_STRING, PACKAGE_BUGREPORT, argv[0]);
		fprintf(stderr, "-f forward.fastq ");
		fprintf(stderr, "-r reverse.fastq ");
		fprintf(stderr, "[-6] ");
		fprintf(stderr, "[-a] ");
		fprintf(stderr, "[-B] ");
		fprintf(stderr, "[-C module1 -C module2 ...] ");
		fprintf(stderr, "[-d flags] ");
		fprintf(stderr, "[-F] ");
		fprintf(stderr, "[-j] ");
		fprintf(stderr, "[-L maxlen] ");
		fprintf(stderr, "[-l minlen] ");
		fprintf(stderr, "[-N] ");
		fprintf(stderr, "[-o minoverlap] ");
		fprintf(stderr, "[-p forwardprimer] ");
		fprintf(stderr, "[-q reverseprimer] ");
#		ifdef HAVE_PTHREAD
		fprintf(stderr, "[-T threads] ");
#		endif
		fprintf(stderr, "[-t threshold] ");
		fprintf(stderr, "\n");
		fprintf(stderr, "\t-6\tUse PHRED+64 (CASAVA 1.3-1.7) instead of PHRED+33 (CASAVA 1.8+).\n");
		fprintf(stderr, "\t-a\tStrip the primers after assembly, rather than before.\n");
		fprintf(stderr, "\t-B\tAllow unbarcoded sequences (try this for BADID errors).\n");
		fprintf(stderr, "\t-C module\tLoad a sequence validation module.\n");
		fprintf(stderr, "\t-d flags\tControl the logging messages. Capital to enable; small to disable.\n");
		fprintf(stderr, "\t\t(R)econstruction detail.\n");
		fprintf(stderr, "\t\tSequence (b)uilding information.\n");
		fprintf(stderr, "\t\t(F)ile processing.\n");
		fprintf(stderr, "\t\t(k)-mer table construction.\n");
		fprintf(stderr, "\t\tShow every (m)ismatch.\n");
		fprintf(stderr, "\t\tOptional (s)tatistics.\n");
		fprintf(stderr, "\t-f\tInput FASTQ file containing forward reads.\n");
		fprintf(stderr, "\t-F\tOutput FASTQ instead of FASTA.\n");
		fprintf(stderr, "\t-j\tInput files are bzipped.\n");
		fprintf(stderr, "\t-k kmers\tThe number of k-mers in the table.\n");
		fprintf(stderr, "\t-L maxlen\tMaximum length for a sequence\n");
		fprintf(stderr, "\t-l minlen\tMinimum length for a sequence\n");
		fprintf(stderr, "\t-N\tEliminate all sequences with unknown nucleotides in the output.\n");
		fprintf(stderr, "\t-o minoverlap\tMinimum overlap between forward and reverse reads (default = %d)\n", minoverlap);
		fprintf(stderr, "\t-p\tForward primer sequence or number of bases to be removed.\n");
		fprintf(stderr, "\t-q\tReverse primer sequence or number of bases to be removed.\n");
		fprintf(stderr, "\t-r\tInput FASTQ file containing reverse reads.\n");
#		ifdef HAVE_PTHREAD
		fprintf(stderr, "\t-T thread\tRun with a number of parallel threads.\n");
#		endif
		fprintf(stderr, "\t-t\tThe minimum probability that a sequence must have to match a primer. (default = %e)\n", threshold);
		for (it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s(%s) %s\n\t%s\n", panda_module_get_name(modules[it]), panda_module_get_description(modules[it]), panda_module_get_version(modules[it]), panda_module_get_usage(modules[it]));
			panda_module_unref(modules[it]);
		}
		return 1;
	}
	if (maxlen < minlen || minoverlap > maxlen) {
		fprintf(stderr, "The minimum length (%d), maximum length (%d), and minimum overlap (%d) are not sensible if you think about the sequence reconstruction.\n", (int) minlen, (int) maxlen, minoverlap);
		for (it = 0; it < modules_length; it++)
			panda_module_unref(modules[it]);
		return 1;
	}
	fprintf(stderr, "INFO\tVER\t%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);

	for (it = 0; it < argc; it++) {
		fprintf(stderr, "INFO\tARG[%d]\t%s\n", (int) it, argv[it]);
	}

	modules[modules_length++] = panda_module_new("SHORT", short_check, NULL, NULL, NULL);
	modules[modules_length++] = panda_module_new("LONG", long_check, NULL, NULL, NULL);

#ifdef HAVE_PTHREAD
	mux = bzip ? panda_mux_open_bz2(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin, policy) : panda_mux_open_gz(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin, policy);
	if (mux == NULL) {
		fprintf(stderr, "ERR\tLIB\tCould not create multiplexer.\n");
		for (it = 0; it < modules_length; it++)
			panda_module_unref(modules[it]);
		return 1;
	}
	assembler = panda_mux_create_assembler_kmer(mux, num_kmers);

#else
	assembler = bzip ? panda_assembler_open_bz2(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin, policy) : panda_assembler_open_gz(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin, policy);
#endif
	if (assembler == NULL) {
		fprintf(stderr, "ERR\tLIB\tCould not create assembler.\n");
		for (it = 0; it < modules_length; it++)
			panda_module_unref(modules[it]);
#if HAVE_PTHREAD
		panda_mux_unref(mux);
#endif
		return 1;
	}
	if (no_algn_file_name != NULL) {
		FILE *no_algn_file = fopen(no_algn_file_name, "w");
		if (no_algn_file != NULL) {
#if HAVE_PTHREAD
			panda_mux_set_fail_alignment(mux, (PandaFailAlign) panda_output_fail, no_algn_file, (PandaDestroy) fclose);
#else
			panda_assembler_set_fail_alignment(assembler, (PandaFailAlign) panda_output_fail, no_algn_file, (PandaDestroy) fclose);
#endif
		} else {
			perror(no_algn_file_name);
		}
	}
	okay = true;
	for (it = 0; it < modules_length; it++) {
		bool add = panda_assembler_add_module(assembler, modules[it]);
		if (!add) {
			fprintf(stderr, "Problem with %s(%s) %s\n\t%s\n", panda_module_get_name(modules[it]), panda_module_get_description(modules[it]), panda_module_get_version(modules[it]), panda_module_get_usage(modules[it]));
		}
		okay &= add;
		panda_module_unref(modules[it]);
	}
	if (!okay) {
		panda_assembler_unref(assembler);
#if HAVE_PTHREAD
		panda_mux_unref(mux);
#endif
		return 1;
	}
	modules_length = 0;
	panda_assembler_set_threshold(assembler, threshold);
	panda_assembler_set_minimum_overlap(assembler, minoverlap);
	panda_assembler_set_disallow_degenerates(assembler, no_n);
	panda_assembler_set_primers_after(assembler, primers_after);

	if (forward_primer != NULL) {
		if (!set_primer(assembler, panda_assembler_set_forward_primer, forward_primer, panda_nt_from_ascii)) {
			fprintf(stderr, "ERR\tBADNT\tFPRIMER\n");
			panda_assembler_unref(assembler);
#if HAVE_PTHREAD
			panda_mux_unref(mux);
#endif
			return 1;
		}
	} else {
		panda_assembler_set_forward_trim(assembler, foffset - 1);
	}
	if (reverse_primer != NULL) {
		if (!set_primer(assembler, panda_assembler_set_reverse_primer, reverse_primer, panda_nt_from_ascii_complement)) {
			fprintf(stderr, "ERR\tBADNT\tRPRIMER\n");
			panda_assembler_unref(assembler);
#if HAVE_PTHREAD
			panda_mux_unref(mux);
#endif
			return 1;
		}
	} else {
		panda_assembler_set_reverse_trim(assembler, roffset - 1);
	}

#if HAVE_PTHREAD
	if (threads > 1) {
		thread_list = malloc(sizeof(pthread_t) * (threads - 1));
		for (it = 0; it < threads - 1; it++) {
			PandaAssembler slave_assembler = panda_mux_create_assembler(mux);
			if (slave_assembler == NULL) {
				fprintf(stderr, "ERR\tMUXCREATE\t%d\n", (int) it + 1);
				threads = it + 1;
				break;
			}
			panda_assembler_copy_configuration(slave_assembler, assembler);
			if (pthread_create(&thread_list[it], NULL, (void *(*)(void *)) do_assembly, slave_assembler) != 0) {
				fprintf(stderr, "ERR\tPCREATE\t%d\n", (int) it + 1);
				threads = it + 1;
				break;
			}
		}
	}
	panda_mux_unref(mux);
#endif
	do_assembly(assembler);
#if HAVE_PTHREAD
	if (threads > 1) {
		for (it = 0; it < threads - 1; it++) {
			pthread_join(thread_list[it], NULL);
		}
		free(thread_list);
	}
#endif
	return some_seqs ? 0 : 1;
}
