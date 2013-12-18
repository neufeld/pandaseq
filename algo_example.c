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

/*
 * To create a new scoring algorithm:
 * 1. Make a copy of this file with the name of your algorithm.
 * 2. Edit `Makefile.am` and include the new file in the `libpandaseq_la_SOURCES` list.
 * 3. Edit `algo.c` and include your algorithm in the `panda_algorithms` list.
 * 4. Edit `pandaseq-algorithm.h` and add a new section for your algorithm and the weird 3-line definition stanza.
 * 5. Fill in this file, renaming "example" to the name of your algorithm.
 * 6. Compile and test.
 * 7. Edit `pandaseq.1` and include documentation about the parameters.
 * 8. (Optional) create a Vala class in `vapi.in`.
 *
 * Have a look at the existing algorithms to get an idea of how to write these.
 */

struct example {
	/* Create all the parameters your algorithm needs. */
};

static double overlap_probability(
	struct example *data,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	size_t overlap) {

	/* Compute the probability of this overlap being correct, as a log probability and return the value. The overlap region may be longer than either read, so be sure to handle those cases. */
}

static double match_probability(
	struct example *data,
	bool match,
	char a,
	char b) {
	/* Compute the log probability that two bases, of scores `a` and `b` are either matched or mismatched based on `match`. If a calculation can be transformed into a lookup table, it can be precomputed in `mktable.c`. */
}

static PandaAlgorithm from_string(
	const char *argument) {
	PandaAlgorithm algo;

	/* Parse the possibly null command line argument and return a new algorithm. */
}

/* This is the class definition. Just give it a name. */
const struct panda_algorithm_class panda_algorithm_simple_bayes_class = {
	.data_size = sizeof(struct example),
	.name = "example",
	.create = from_string,
	.data_destroy = NULL,
	.overlap_probability = (PandaComputeOverlap) overlap_probability,
	.match_probability = (PandaComputeMatch) match_probability,
	.prob_unpaired = qual_nn_simple_bayesian,
};

/* The constructor for your algorithm. It needs to call the super constructor with your class, then you may initialise your variables. Also, include the definition in pandaseq-algorithm.h. */
PandaAlgorithm panda_algorithm_example_new(
	void) {
	PandaAlgorithm algo = panda_algorithm_new(&panda_algorithm_example_class);
	/* Set default parameters here. Use the setters defined below. */
	return algo;
}

/* Create getters and setters for all the parameters of the algorithm. Create a pair for each parameter needed. Also, include the definition in pandaseq-algorithm.h. */
double panda_algorithm_example_get_parameter(
	PandaAlgorithm algorithm) {
	if (panda_algorithm_is_a(algorithm, &panda_algorithm_example_bayes_class)) {
		return ((struct example *) panda_algorithm_data(algorithm))->parameter;
	} else {
		return -1;
	}
}

void panda_algorithm_example_set_parameterrror_estimation(
	PandaAlgorithm algorithm,
	double parameter) {
	if (panda_algorithm_is_a(algorithm, &panda_algorithm_example_class)) {
		struct example *data = panda_algorithm_data(algorithm);
		data->parameter = parameter;
	}
}
