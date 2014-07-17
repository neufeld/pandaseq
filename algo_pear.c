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
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pandaseq.h"
#include "algo.h"
#include "prob.h"
#include "table.h"

struct pear {
	double random_base;
};

static double overlap_probability(
	struct pear *data,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	size_t overlap) {
	double probability = 0;
	size_t i;

	for (i = 0; i < overlap; i++) {
		int findex = forward_length + i - overlap;
		int rindex = reverse_length - i - 1;
		if (findex < 0 || rindex < 0 || (size_t) findex >= forward_length || (size_t) rindex >= reverse_length)
			continue;
		panda_nt f = forward[findex].nt;
		panda_nt r = reverse[rindex].nt;
		if (PANDA_NT_IS_N(f) || PANDA_NT_IS_N(r)) {
			probability -= data->random_base;
		} else if ((f & r) != 0) {
			probability += qual_match_pear[PHREDCLAMP(forward[findex].qual)][PHREDCLAMP(forward[rindex].qual)];
		} else {
			probability += qual_mismatch_pear[PHREDCLAMP(forward[findex].qual)][PHREDCLAMP(forward[rindex].qual)];
		}
	}

	return probability;
}

static double match_probability(
	struct pear *data,
	bool match,
	char a,
	char b) {
	(void) data;
	return (match ? qual_match_pear : qual_mismatch_pear)[PHREDCLAMP(a)][PHREDCLAMP(b)];
}

static PandaAlgorithm from_string(
	const char *argument) {
	PandaAlgorithm algo;
	double random_base;
	char *end;

	if (argument == NULL)
		return panda_algorithm_pear_new();
	errno = 0;
	random_base = strtod(argument, &end);
	if (errno == ERANGE || *end != '\0') {
		fprintf(stderr, "Cannot parse value: %s\n", argument);
		return NULL;
	}
	if (random_base < 0 || random_base > 1) {
		fprintf(stderr, "Random base %f is not a probability.\n", random_base);
		return NULL;
	}
	algo = panda_algorithm_pear_new();
	panda_algorithm_pear_set_random_base_log_p(algo, log(random_base));
	return algo;
}

const struct panda_algorithm_class panda_algorithm_pear_class = {
	.data_size = sizeof(struct pear),
	.name = "pear",
	.create = from_string,
	.data_destroy = NULL,
	.overlap_probability = (PandaComputeOverlap) overlap_probability,
	.match_probability = (PandaComputeMatch) match_probability,
	.prob_unpaired = qual_nn_simple_bayesian,
};

PandaAlgorithm panda_algorithm_pear_new(
	void) {
	PandaAlgorithm algo = panda_algorithm_new(&panda_algorithm_pear_class);
	panda_algorithm_pear_set_random_base_log_p(algo, log(0.25));
	return algo;
}

void panda_algorithm_pear_set_random_base_log_p(
	PandaAlgorithm algorithm,
	double log_p) {
	if (panda_algorithm_is_a(algorithm, &panda_algorithm_pear_class)) {
		((struct pear *) panda_algorithm_data(algorithm))->random_base = log_p;
	}
}

double panda_algorithm_pear_get_random_base_log_p(
	PandaAlgorithm algorithm) {
	if (panda_algorithm_is_a(algorithm, &panda_algorithm_pear_class)) {
		return ((struct pear *) panda_algorithm_data(algorithm))->random_base;
	} else {
		return 1;
	}
}
