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
#include <stdlib.h>
#include "pandaseq.h"
#include "algo.h"
#include "prob.h"
#include "table.h"
#include <stdio.h>

struct simple_bayes {
	double q;
	double pmatch;
	double pmismatch;
};

static double overlap_probability(
	struct simple_bayes *data,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	size_t overlap) {
	size_t matches = 0;
	size_t mismatches = 0;
	size_t unknowns = 0;
	size_t i;

	for (i = 0; i < overlap; i++) {
		int findex = forward_length + i - overlap;
		int rindex = reverse_length - i - 1;
		if (findex < 0 || rindex < 0 || findex >= forward_length || rindex >= reverse_length)
			continue;
		panda_nt f = forward[findex].nt;
		panda_nt r = reverse[rindex].nt;
		if (PANDA_NT_IS_N(f) || PANDA_NT_IS_N(r)) {
			unknowns++;
		} else if ((f & r) != 0) {
			matches++;
		} else {
			mismatches++;
		}
	}

	if (overlap >= forward_length && overlap >= reverse_length) {
		return (qual_nn_simple_bayesian * unknowns + matches * data->pmatch + mismatches * data->pmismatch);
	} else {
		return (qual_nn_simple_bayesian * (forward_length + reverse_length - 2 * overlap + unknowns) + matches * data->pmatch + mismatches * data->pmismatch);
	}
}

static double match_probability(
	struct simple_bayes *data,
	bool match,
	char a,
	char b) {
	return (match ? qual_match_simple_bayesian : qual_mismatch_simple_bayesian)[PHREDCLAMP(a)][PHREDCLAMP(b)];
}

static PandaAlgorithm from_string(
	const char *argument) {
	PandaAlgorithm algo;
	double err_estimation;
	char *end;

	if (argument == NULL)
		return panda_algorithm_simple_bayes_new();
	errno = 0;
	err_estimation = strtod(argument, &end);
	if (errno == ERANGE || *end != '\0') {
		fprintf(stderr, "Cannot parse value: %s\n", argument);
		return NULL;
	}
	if (err_estimation < 0 || err_estimation > 1) {
		fprintf(stderr, "Error estimation %f is not a probability.\n", err_estimation);
		return NULL;
	}
	algo = panda_algorithm_simple_bayes_new();
	panda_algorithm_simple_bayes_set_error_estimation(algo, err_estimation);
	return algo;
}

const struct panda_algorithm_class panda_algorithm_simple_bayes_class = {
	.data_size = sizeof(struct simple_bayes),
	.name = "simple_bayesian",
	.create = from_string,
	.data_destroy = NULL,
	.overlap_probability = (PandaComputeOverlap) overlap_probability,
	.match_probability = (PandaComputeMatch) match_probability,
	.prob_unpaired = qual_nn_simple_bayesian,
};

PandaAlgorithm panda_algorithm_simple_bayes_new(
	void) {
	PandaAlgorithm algo = panda_algorithm_new(&panda_algorithm_simple_bayes_class);
	panda_algorithm_simple_bayes_set_error_estimation(algo, 0.36);
	return algo;
}

double panda_algorithm_simple_bayes_get_error_estimation(
	PandaAlgorithm algorithm) {
	if (panda_algorithm_is_a(algorithm, &panda_algorithm_simple_bayes_class)) {
		return ((struct simple_bayes *) panda_algorithm_data(algorithm))->q;
	} else {
		return -1;
	}
}

void panda_algorithm_simple_bayes_set_error_estimation(
	PandaAlgorithm algorithm,
	double q) {
	if (q > 0 && q < 1 && panda_algorithm_is_a(algorithm, &panda_algorithm_simple_bayes_class)) {
		struct simple_bayes *data = panda_algorithm_data(algorithm);
		data->q = q;
		data->pmatch = log(0.25 * (1 - 2 * q + q * q));
		data->pmismatch = log((3 * q - 2 * q * q) / 18.0);
	}
}
