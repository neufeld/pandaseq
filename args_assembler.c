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
#include<ctype.h>
#include<errno.h>
#include<stdlib.h>
#include<string.h>
#include "config.h"
#include "pandaseq.h"

static bool set_primers_after(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	panda_assembler_set_primers_after(assembler, is_set);
	return true;
}

const panda_tweak_assembler panda_stdargs_primers_after = { 'a', NULL, "Strip the primers after assembly, rather than before.", set_primers_after };

static bool set_degenerates(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	panda_assembler_set_disallow_degenerates(assembler, is_set);
	return true;
}

const panda_tweak_assembler panda_stdargs_degenerates = { 'N', NULL, "Eliminate all sequences with unknown nucleotides in the output.", set_degenerates };

static bool set_threshold(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	double threshold = 0.6;
	if (argument != NULL) {
		errno = 0;
		threshold = strtod(argument, NULL);
		if (errno != 0 || threshold < 0 || threshold > 1) {
			fprintf(stderr, "Bad threshold: %s. It should be between 0 and 1.\n", argument);
			free(argument);
			return false;
		}
	}
	free(argument);
	panda_assembler_set_threshold(assembler, threshold);
	return true;
}

const panda_tweak_assembler panda_stdargs_threshold = { 't', "threshold", "The minimum probability that a sequence must have to match a primer.", set_threshold };

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
	}
	free(argument);
	return true;
}

static bool set_primer_group(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	if (flag == 'p') {
		return set_primer(assembler, argument, "forward", panda_assembler_set_forward_trim, panda_assembler_set_forward_primer, panda_nt_from_ascii_complement);
	} else if (flag == 'q') {
		return set_primer(assembler, argument, "reverse", panda_assembler_set_reverse_trim, panda_assembler_set_reverse_primer, panda_nt_from_ascii_complement);
	}
	free(argument);
	return false;
}

const panda_tweak_assembler panda_stdargs_forward_primer = { 'p', "primer", "Forward primer sequence or number of bases to be removed.", set_primer_group };
const panda_tweak_assembler panda_stdargs_reverse_primer = { 'q', "primer", "Reverse primer sequence or number of bases to be removed.", set_primer_group };

bool short_check(
	const panda_result_seq *sequence,
	void *user_data) {
	return sequence->sequence_length >= (size_t) user_data;
}

static bool set_short_check(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	size_t minlen;
	PandaModule m;

	if (argument == NULL) {
		free(argument);
		return true;
	}
	errno = 0;
	minlen = (size_t) strtol(argument, NULL, 10);
	if (errno != 0 || minlen < 0 || minlen > 2 * PANDA_MAX_LEN) {
		fprintf(stderr, "Bad minimum length.\n");
		free(argument);
		return false;
	}
	m = panda_module_new("SHORT", short_check, NULL, (void *) minlen, NULL);
	panda_assembler_add_module(assembler, m);
	panda_module_unref(m);
	free(argument);
	return true;
}

const panda_tweak_assembler panda_stdargs_min_len = { 'l', "length", "Minimum length for a sequence.", set_short_check };

bool long_check(
	const panda_result_seq *sequence,
	void *user_data) {
	size_t length = (size_t) user_data;
	return sequence->sequence_length <= (size_t) user_data;
}

static bool set_long_check(
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set) {
	size_t maxlen;
	PandaModule m;

	if (argument == NULL) {
		free(argument);
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

const panda_tweak_assembler panda_stdargs_max_len = { 'L', "length", "Maximum length for a sequence.", set_long_check };

const panda_tweak_assembler *const panda_stdargs[] = {
	&panda_stdargs_max_len,
	&panda_stdargs_degenerates,
	&panda_stdargs_primers_after,
	&panda_stdargs_min_len,
	&panda_stdargs_forward_primer,
	&panda_stdargs_reverse_primer,
	&panda_stdargs_threshold,
};

const size_t panda_stdargs_length = sizeof(panda_stdargs) / sizeof(panda_tweak_assembler *);
