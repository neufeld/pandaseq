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
#include<ltdl.h>
#include<limits.h>
#include<math.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq.h"

#define STR0(x) #x
#define STR(x) STR0(x)
static void printtime(long count, time_t starttime)
{
	time_t now;
	(void)time(&now);
	fprintf(stderr, "STAT\tTIME\t%s\nSTAT\tELAPSED\t%d\nSTAT\tREADS\t%ld\n",
		ctime(&now), (int)(now - starttime), count);
}

bool set_primer(PandaAssembler assembler, void(*set_func)(PandaAssembler,panda_nt*,size_t), const char *str, panda_nt(*parse)(char)) {
	panda_nt buffer[PANDA_MAX_LEN];
	size_t it;
	for(it = 0; it < strlen(str); it++) {
		if ((buffer[it] = parse(str[it])) == '\0') {
			return false;
		}
	}
	set_func(assembler, buffer, strlen(str));
	return true;
}
int main(int argc, char **argv)
{
	PandaAssembler assembler;
	bool bzip = false;
	int c;
	long count;
	bool fastq = false;
	size_t foffset = 1;
	char *forward_filename = NULL;
	char *forward_primer = NULL;
	bool help = false;
	size_t it;
	long longcount = 0;
	ssize_t maxlen = SSIZE_MAX;
	size_t minlen = 0;
	int minoverlap = 1;
	PandaModule modules[100];
	size_t modules_length = 0;
	bool no_n = false;
	double q = 0.36;
	int qualmin = 33;
	const panda_result_seq *result;
	char *reverse_filename = NULL;
	char *reverse_primer = NULL;
	size_t roffset = 1;
	long shortcount = 0;
	time_t starttime;
	double threshold;
	bool version = false;
	(void)time(&starttime);
	threshold = 0.6;

	if (lt_dlinit() != 0) {
		fprintf(stderr, "ERR\tLTLD\tINIT\n");
		return 1;
	}
	if (lt_dladdsearchdir(STR(PKGLIBDIR)) != 0) {
		fprintf(stderr, "ERR\tLTLD\tPKGDIR\t%s\n", STR(PKGLIBDIR));
		(void)lt_dlexit();
		return 1;
	}

	/* Process command line arguments. */
	while ((c = getopt(argc, argv, "hvjp:q:f:r:t:o:Nl:L:Q:C:6F")) != -1) {
		char *endptr;
		switch (c) {
		case 'h':
			help = true;
			break;
		case 'v':
			version = true;
			break;
		case 'j':
			bzip = true;
			break;
		case 't':
			errno = 0;
			threshold = strtod(optarg, NULL);
			if (errno != 0 || threshold < 0 || threshold > 1) {
				fprintf(stderr,
					"Bad threshold. It should be between 0 and 1.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			}
			break;
		case 'Q':
			errno = 0;
			q = strtod(optarg, NULL);
			if (errno != 0 || q < 0 || q > 1) {
				fprintf(stderr,
					"Bad quality. It should be between 0 and 1.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			}
			break;
		case 'l':
			errno = 0;
			minlen = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || minlen < 0) {
				fprintf(stderr, "Bad minimum length.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			}
			break;
		case 'L':
			errno = 0;
			maxlen = (size_t) strtol(optarg, NULL, 10);
			if (errno != 0 || minlen < 0) {
				fprintf(stderr, "Bad maximum length.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			}
			break;
		case 'o':
			errno = 0;
			minoverlap = strtol(optarg, NULL, 10);
			if (errno != 0 || minoverlap < 1) {
				fprintf(stderr, "Bad minimum overlap.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			}
			break;
		case 'f':
			forward_filename = optarg;
			break;
		case 'r':
			reverse_filename = optarg;
			break;
		case 'N':
			no_n = true;
			break;
		case 'F':
			fastq = 1;
			break;
		case 'p':
			errno = 0;
			foffset = strtol(optarg, &endptr, 10);
			if (*endptr != '\0') {
				forward_primer = optarg;
			} else if (errno != 0 || foffset < 1) {
				fprintf(stderr, "Bad forward primer length.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
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
			} else if (errno != 0 || roffset < 1) {
				fprintf(stderr, "Bad reverse primer length.\n");
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				(void)lt_dlexit();
				return 1;
			} else {
				roffset++;
			}
			break;
		case 'C':
			if ((modules[modules_length++] = panda_module_load(optarg)) == NULL) {
				for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
				lt_dlexit();
				return 1;
			}
			break;
		case '6':
			qualmin = 64;
			break;
		case '?':
			if (optopt == (int)'f' || optopt == (int)'r'
			    || optopt == (int)'p' || optopt == (int)'q'
			    || optopt == (int)'l' || optopt == (int)'L'
			    || optopt == (int)'Q' || optopt == (int)'C') {
				fprintf(stderr,
					"Option -%c requires an argument.\n",
					optopt);
			} else if (isprint(optopt)) {
				fprintf(stderr, "Unknown option `-%c'.\n",
					optopt);
			} else {
				fprintf(stderr,
					"Unknown option character `\\x%x'.\n",
					(unsigned int)optopt);
			}
			lt_dlexit();
			return 1;
		default:
			abort();
		}
	}

	if (version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		for(it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s %s\n", panda_module_get_name(modules[it]), panda_module_get_version(modules[it]));
			panda_module_unref(modules[it]);
		}
		lt_dlexit();
		return 1;
	}
	if (forward_filename == NULL || reverse_filename == NULL || help) {
		fprintf(stderr,
			"%s <%s>\nUsage: %s -f forward.fastq -r reverse.fastq [-j] [-p forwardprimer] [-q reverseprimer] [-t threshold] [-N] [-o minoverlap] [-l minlen] [-L maxlen] [ -C module1 -C module2 ...] [-6] [-F]\n\t-f\tInput FASTQ file containing forward reads.\n\t-r\tInput FASTQ file containing reverse reads.\n\t-j\tInput files are bzipped.\n\t-p\tForward primer sequence or number of bases to be removed.\n\t-q\tReverse primer sequence or number of bases to be removed.\n\t-t\tThe minimum probability that a sequence must have to match a primer. (default = %e)\n\t-N\tEliminate all sequences with unknown nucleotides in the output.\n\t-o minoverlap\tMinimum overlap between forward and reverse reads (default = %d)\n\t-l minlen\tMinimum length for a sequence\n\t-L maxlen\tMaximum length for a sequence\n\t-C module\tLoad a sequence validation module.\n\t-6\tUse PHRED+64 (CASAVA 1.3-1.7) instead of PHRED+33 (CASAVA 1.8+).\n\t-F\tOutput FASTQ instead of FASTA.\n",
			PACKAGE_STRING, PACKAGE_BUGREPORT,
			argv[0], threshold, minoverlap);
			for(it = 0; it < modules_length; it++) {
			fprintf(stderr, "%s(%s) %s\n\t%s\n", panda_module_get_name(modules[it]), panda_module_get_description(modules[it]), panda_module_get_version(modules[it]), panda_module_get_usage(modules[it]));
			panda_module_unref(modules[it]);
		}
		lt_dlexit();
		return 1;
	}
	fprintf(stderr, "INFO\tVER\t%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);

	assembler = bzip ? panda_assembler_open_bz2(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin) : panda_assembler_open_gz(forward_filename, reverse_filename, (PandaLogger) panda_logger_file, stderr, NULL, qualmin);

	if (assembler == NULL) {
		fprintf(stderr, "ERR\tLIB\tCould not create assembler.\n");
		for(it = 0; it < modules_length; it++) panda_module_unref(modules[it]);
		lt_dlexit();
		return 1;
	}
	for(it = 0; it < modules_length; it++) {
		panda_assembler_add_module(assembler, modules[it]);
		panda_module_unref(modules[it]);
	}
	modules_length = 0;
	panda_assembler_set_threshold(assembler, threshold);
	panda_assembler_set_minimum_overlap(assembler, minoverlap);
	panda_assembler_set_disallow_degenerates(assembler, no_n);

	if (forward_primer != NULL) {
		if (!set_primer(assembler, panda_assembler_set_forward_primer, forward_primer, panda_nt_from_ascii)) {
			fprintf(stderr, "ERR\tBADNT\tFPRIMER\n");
			panda_assembler_unref(assembler);
			(void)lt_dlexit();
			return 1;
		}
	} else {
		panda_assembler_set_forward_trim(assembler, foffset - 1);
	}
	if (reverse_primer != NULL) {
		if (!set_primer(assembler, panda_assembler_set_reverse_primer, reverse_primer, panda_nt_from_ascii_complement)) {
			fprintf(stderr, "ERR\tBADNT\tRPRIMER\n");
			panda_assembler_unref(assembler);
			(void)lt_dlexit();
			return 1;
		}
	} else {
		panda_assembler_set_reverse_trim(assembler, roffset - 1);
	}
	while ((result = panda_assembler_next(assembler)) != NULL) {
		count = panda_assembler_get_count(assembler);
		if (count % 1000 == 0) {
			printtime(count, starttime);
		}
		if (result->sequence_length < minlen) {
			fputs("ERR\tSHORT\t%s\n", stderr);
			panda_seqid_print(&result->name, stderr);
			fputc('\n', stderr);
			shortcount++;
		} else if (result->sequence_length > maxlen) {
			fputs("ERR\tSHORT\t%s\n", stderr);
			panda_seqid_print(&result->name, stderr);
			fputc('\n', stderr);
			shortcount++;
		} else if (fastq) {
			panda_output_fastq(result, stdout);
		} else {
			panda_output_fasta(result, stdout);
		}
	}
	count = panda_assembler_get_count(assembler);
	printtime(count, starttime);
	if (forward_primer != NULL)
		fprintf(stderr, "STAT\tNOFP\t%ld\n", panda_assembler_get_no_forward_primer_count(assembler));
	if (reverse_primer != NULL)
		fprintf(stderr, "STAT\tNORP\t%ld\n", panda_assembler_get_no_reverse_primer_count(assembler));
	fprintf(stderr, "STAT\tNOALGN\t%ld\nSTAT\tLOWQ\t%ld\n", panda_assembler_get_failed_alignment_count(assembler),
		panda_assembler_get_low_quality_count(assembler));
	if (no_n)
		fprintf(stderr, "STAT\tDEGENERATE\t%ld\n", panda_assembler_get_degenerate_count(assembler));
	if (minlen > 0)
		fprintf(stderr, "STAT\tSHORT\t%ld\n", shortcount);
	if (maxlen < SSIZE_MAX)
		fprintf(stderr, "STAT\tLONG\t%ld\n", longcount);
	fprintf(stderr, "STAT\tOK\t%ld\n", panda_assembler_get_ok_count(assembler) - shortcount - longcount);

	panda_assembler_unref(assembler);
	(void)lt_dlexit();
	return 0;
}
