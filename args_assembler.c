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
#include "config.h"
#include<ctype.h>
#include<errno.h>
#include<ltdl.h>
#include<stdlib.h>
#include<string.h>
#include "pandaseq.h"

static bool set_algorithm(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	char *extra;
	size_t it;
	(void) flag;
	if (argument == NULL)
		return true;
	extra = strchr(argument, LT_PATHSEP_CHAR);
	if (extra != NULL) {
		*extra = '\0';
		extra++;
	}
	for (it = 0; it < panda_algorithms_length; it++) {
		if (strcmp(argument, panda_algorithms[it]->name) == 0) {
			PandaAlgorithm algo;
			if (panda_algorithms[it]->create == NULL) {
				fprintf(stderr, "It seems that no one wrote the code to use this algorithm properly.\n");
				free(argument);
				return false;
			}
			algo = panda_algorithms[it]->create(extra);
			if (algo == NULL) {
				fprintf(stderr, "Unable to initialise requested algorithm.\n");
				free(argument);
				return false;
			}
			panda_assembler_set_algorithm(assembler, algo);
			panda_algorithm_unref(algo);
			free(argument);
			return true;
		}
	}
	fprintf(stderr, "Unknown algorithm: %s\n", argument);
	free(argument);
	return false;
}

const panda_tweak_assembler panda_stdargs_algorithm = { 'A', "algorithm:parameters", "Select the algorithm to use for overlap selection and scoring.", set_algorithm, false };

static bool set_primers_after(
	PandaAssembler assembler,
	char flag,
	char *argument) {

	(void) flag;
	(void) argument;

	panda_assembler_set_primers_after(assembler, true);
	return true;
}

const panda_tweak_assembler panda_stdargs_primers_after = { 'a', NULL, "Strip the primers after assembly, rather than before.", set_primers_after, false };

bool no_n_check(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {

	(void) logger;
	(void) user_data;

	return sequence->degenerates == 0;
}

static bool add_module(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	PandaModule module;
	(void) flag;

	if (argument == NULL) {
		return true;
	}
	module = panda_module_load(panda_assembler_get_logger(assembler), argument);
	free(argument);
	if (module == NULL) {
		return false;
	}
	panda_assembler_add_module(assembler, module);
	panda_module_unref(module);
	return true;
}

const panda_tweak_assembler panda_stdargs_module = { 'C', "filter", "Load a pluggable filter module.", add_module, true };

static bool set_degenerates(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	PandaModule m = panda_module_new("DEGENERATE", no_n_check, NULL, NULL, NULL);
	(void) flag;
	panda_assembler_add_module(assembler, m);
	panda_module_unref(m);
	if (argument != NULL) {
		free(argument);
	}
	return true;
}

const panda_tweak_assembler panda_stdargs_degenerates = { 'N', NULL, "Eliminate all sequences with unknown nucleotides in the output.", set_degenerates, false };

static bool set_threshold(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	double threshold = 0.6;
	(void) flag;
	if (argument != NULL) {
		errno = 0;
		threshold = strtod(argument, NULL);
		if (errno != 0 || threshold < 0 || threshold > 1) {
			fprintf(stderr, "Bad threshold: %s. It should be between 0 and 1.\n", argument);
			free(argument);
			return false;
		}
		free(argument);
	}
	panda_assembler_set_threshold(assembler, threshold);
	return true;
}

const panda_tweak_assembler panda_stdargs_threshold = { 't', "threshold", "The minimum probability that a sequence must have to assemble and, if used, match a primer.", set_threshold, false };

static bool set_primer(
	PandaAssembler assembler,
	char *argument,
	char *direction,
	void (*set_trim) (PandaAssembler,
		size_t),
	void (*set_func) (PandaAssembler,
		panda_nt *,
		size_t),
	panda_nt (*parse) (char)) {
	if (argument != NULL) {
		char *endptr;
		size_t offset;
		errno = 0;
		offset = strtol(argument, &endptr, 10);
		if (*endptr != '\0') {
			panda_nt buffer[PANDA_MAX_LEN];
			size_t it;
			for (it = 0; it < strlen(argument); it++) {
				if ((buffer[it] = parse(argument[it])) == '\0') {
					fprintf(stderr, "ERR\tBADNT\t%cPRIMER\n", (int) toupper(direction[0]));
					free(argument);
					return false;
				}
			}
			set_func(assembler, buffer, strlen(argument));
		} else if (errno != 0 || offset < 1 || offset > PANDA_MAX_LEN) {
			fprintf(stderr, "Bad %s primer length.\n", direction);
			free(argument);
			return false;
		} else {
			set_trim(assembler, offset);
		}
		free(argument);
	}
	return true;
}

static bool set_primer_group(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	if (flag == 'p') {
		return set_primer(assembler, argument, "forward", panda_assembler_set_forward_trim, panda_assembler_set_forward_primer, panda_nt_from_ascii);
	} else if (flag == 'q') {
		return set_primer(assembler, argument, "reverse", panda_assembler_set_reverse_trim, panda_assembler_set_reverse_primer, panda_nt_from_ascii_complement);
	}
	if (argument != NULL) {
		free(argument);
	}
	return false;
}

const panda_tweak_assembler panda_stdargs_forward_primer = { 'p', "primer", "Forward primer sequence or number of bases to be removed.", set_primer_group, false };
const panda_tweak_assembler panda_stdargs_reverse_primer = { 'q', "primer", "Reverse primer sequence or number of bases to be removed.", set_primer_group, false };

bool short_check(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {
	(void) logger;
	return sequence->sequence_length >= (size_t) user_data;
}

static bool set_short_check(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	long int minlen;
	PandaModule m;
	(void) flag;

	if (argument == NULL) {
		return true;
	}
	errno = 0;
	minlen = (size_t) strtol(argument, NULL, 10);
	if (errno != 0 || minlen < 0 || (size_t) minlen > 2 * PANDA_MAX_LEN) {
		fprintf(stderr, "Bad minimum length.\n");
		free(argument);
		return false;
	}
	m = panda_module_new("SHORT", short_check, NULL, (void *) (size_t) minlen, NULL);
	panda_assembler_add_module(assembler, m);
	panda_module_unref(m);
	free(argument);
	return true;
}

const panda_tweak_assembler panda_stdargs_min_len = { 'l', "length", "Minimum length for a sequence.", set_short_check, false };

bool long_check(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {
	size_t length = (size_t) user_data;
	(void) logger;
	return sequence->sequence_length <= length;
}

static bool set_long_check(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	size_t maxlen;
	PandaModule m;

	(void) flag;
	if (argument == NULL) {
		return true;
	}
	errno = 0;
	maxlen = (size_t) strtol(argument, NULL, 10);
	if (errno != 0 || maxlen < 1 || maxlen > 2 * PANDA_MAX_LEN) {
		fprintf(stderr, "Bad maximum length.\n");
		free(argument);
		return false;
	}

	m = panda_module_new("LONG", long_check, NULL, (void *) maxlen, NULL);
	panda_assembler_add_module(assembler, m);
	panda_module_unref(m);
	free(argument);
	return true;
}

const panda_tweak_assembler panda_stdargs_max_len = { 'L', "length", "Maximum length for a sequence.", set_long_check, false };

static bool set_minimum_overlap(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	int min_overlap;

	(void) flag;
	if (argument == NULL) {
		return true;
	}
	errno = 0;
	min_overlap = strtol(argument, NULL, 10);
	if (errno != 0 || min_overlap < 1 || (size_t) min_overlap > 2 * PANDA_MAX_LEN) {
		fprintf(stderr, "Bad overlap length.\n");
		free(argument);
		return false;
	}

	panda_assembler_set_minimum_overlap(assembler, min_overlap);
	free(argument);
	return true;
}

const panda_tweak_assembler panda_stdargs_min_overlap = { 'o', "length", "Minumum overlap region length for a sequence.", set_minimum_overlap, false };

static bool set_maxiumum_overlap(
	PandaAssembler assembler,
	char flag,
	char *argument) {
	int max_overlap;

	(void) flag;
	if (argument == NULL) {
		return true;
	}
	errno = 0;
	max_overlap = strtol(argument, NULL, 10);
	if (errno != 0 || max_overlap < 0 || (size_t) max_overlap > 2 * PANDA_MAX_LEN) {
		fprintf(stderr, "Bad overlap length.\n");
		free(argument);
		return false;
	}

	panda_assembler_set_maximum_overlap(assembler, max_overlap);
	free(argument);
	return true;
}

const panda_tweak_assembler panda_stdargs_max_overlap = { 'O', "length", "Maximum overlap region length for a sequence. (0 to use read length.)", set_maxiumum_overlap, false };

const panda_tweak_assembler *const panda_stdargs[] = {
	&panda_stdargs_algorithm,
	&panda_stdargs_module,
	&panda_stdargs_max_len,
	&panda_stdargs_degenerates,
	&panda_stdargs_max_overlap,
	&panda_stdargs_primers_after,
	&panda_stdargs_min_len,
	&panda_stdargs_min_overlap,
	&panda_stdargs_forward_primer,
	&panda_stdargs_reverse_primer,
	&panda_stdargs_threshold,
};

const size_t panda_stdargs_length = sizeof(panda_stdargs) / sizeof(panda_tweak_assembler *);
