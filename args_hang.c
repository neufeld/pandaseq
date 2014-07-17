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
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "config.h"
#include "pandaseq.h"
#include "misc.h"

struct panda_args_hang {
	void *user_data;
	PandaDestroy destroy;
	PandaOpener opener;
	PandaSetup setup;
	PandaTweakGeneral tweak;
	panda_nt *forward;
	size_t forward_length;
	panda_nt *reverse;
	size_t reverse_length;
	bool skip;
	double threshold;
};

PandaArgsHang panda_args_hang_new(
	void *user_data,
	PandaDestroy destroy,
	PandaTweakGeneral tweak,
	PandaOpener opener,
	PandaSetup setup) {
	PandaArgsHang data = malloc(sizeof(struct panda_args_hang));
	data->user_data = user_data;
	data->destroy = destroy;
	data->opener = opener;
	data->setup = setup;
	data->tweak = tweak;
	data->skip = false;
	data->threshold = log(0.6);

	data->forward = calloc(PANDA_MAX_LEN, sizeof(panda_nt));
	data->reverse = calloc(PANDA_MAX_LEN, sizeof(panda_nt));

	return data;
}

void panda_args_hang_free(
	PandaArgsHang data) {
	if (data->destroy != NULL) {
		data->destroy(data->user_data);
	}
	if (data->forward != NULL)
		free(data->forward);
	if (data->reverse != NULL)
		free(data->reverse);
	free(data);
}

static bool set_cutoff_primer(
	panda_nt *array,
	size_t *length,
	const char *argument,
	panda_nt (*parse) (char),
	const char *direction) {
	size_t it;
	*length = strlen(argument);
	if (*length > PANDA_MAX_LEN) {
		fprintf(stderr, "The %s primer given for overhang cut-off is too long.\n", direction);
		return false;
	}
	for (it = 0; it < *length; it++) {
		if ((array[it] = parse(argument[it])) == '\0') {
			fprintf(stderr, "ERR\tBADNT\t%cHANG\n", (int) toupper(direction[0]));
			return false;
		}
	}
	return true;
}

bool panda_args_hang_tweak(
	PandaArgsHang data,
	char flag,
	const char *argument) {
	double threshold;

	switch (flag) {
	case 'P':
		return set_cutoff_primer(data->forward, &data->forward_length, argument, panda_nt_from_ascii, "forward");
	case 'Q':
		return set_cutoff_primer(data->reverse, &data->reverse_length, argument, panda_nt_from_ascii_complement, "reverse");
	case 's':
		data->skip = true;
		return true;
	case 't':
		errno = 0;
		threshold = strtod(argument, NULL);
		if (errno != 0 || threshold < 0 || threshold > 1) {
			fprintf(stderr, "Bad threshold: %s. It should be between 0 and 1.\n", argument);
			return false;
		}
		data->threshold = log(threshold);
		return true;
	default:
		return data->tweak(data->user_data, flag, argument);
	}
}

static const panda_tweak_general hang_forward = { 'P', false, "primer", "The sequence in the forward read overhang cut-off.", false };
static const panda_tweak_general hang_reverse = { 'Q', false, "primer", "The sequence in the reverse read overhang cut-off.", false };
static const panda_tweak_general hang_skip = { 's', true, NULL, "Skip reads that do not contain the primer (otherwise, assembly will be attempted).", false };
static const panda_tweak_general hang_threshold = { 't', true, "threshold", "The minimum probability that a sequence must have to match a primer.", false };

const panda_tweak_general **panda_args_hang_args(
	const panda_tweak_general *const *const general_args,
	size_t general_args_length,
	size_t *length) {
	const panda_tweak_general **array = calloc(4, sizeof(panda_tweak_general *));
	array[0] = &hang_forward;
	array[1] = &hang_reverse;
	array[2] = &hang_skip;
	array[3] = &hang_threshold;
	*length = 4;
	panda_tweak_general_append(&array, length, general_args, general_args_length);
	return array;
}

PandaNextSeq panda_args_hang_opener(
	PandaArgsHang data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy) {
	MANAGED_STACK(PandaNextSeq,
		inner);

	inner = data->opener(data->user_data, logger, fail, fail_data, fail_destroy, &inner_data, &inner_destroy);
	if (inner == NULL) {
		MAYBE(next_data) = NULL;
		MAYBE(next_destroy) = NULL;
		return NULL;
	}

	return panda_trim_overhangs(inner, inner_data, inner_destroy, logger, data->forward, data->forward_length, data->reverse, data->reverse_length, data->skip, data->threshold, next_data, next_destroy);
}

bool panda_args_hang_setup(
	PandaArgsHang data,
	PandaAssembler assembler) {
	panda_assembler_set_threshold(assembler, data->threshold);
	return data->setup(data->user_data, assembler);
}
