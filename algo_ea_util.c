/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella

    Based on the work of Expression Analysis / Erik Aronesty:
      https://code.google.com/p/ea-utils/wiki/FastqJoin

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
#include <stdio.h>
#include <stdlib.h>
#include "pandaseq.h"
#include "prob.h"
#include "table.h"

static double overlap_probability(
	void *data,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	size_t overlap) {
	size_t mismatches = 0;
	size_t i;
	size_t real_overlap = 0;

	(void) data;

	for (i = 0; i < overlap; i++) {
		int findex = forward_length + i - overlap;
		int rindex = reverse_length - i - 1;
		if (findex < 0 || rindex < 0 || (size_t) findex >= forward_length || (size_t) rindex >= reverse_length)
			continue;
		panda_nt f = forward[findex].nt;
		panda_nt r = reverse[rindex].nt;
		if (PANDA_NT_IS_N(f) || PANDA_NT_IS_N(r) || (f & r) == 0) {
			mismatches++;
		}
		real_overlap++;
	}

	return (((double) mismatches) * mismatches + 1) / real_overlap;
}

static double match_probability(
	void *data,
	bool match,
	char a,
	char b) {
	int score = (a > b) ? PHREDCLAMP(a) : PHREDCLAMP(b);
	(void) data;
	(void) match;
	return qual_score[score];
}

static PandaAlgorithm from_string(
	const char *argument) {

	if (argument != NULL) {
		fprintf(stderr, "No arguments allowed: %s\n", argument);
		return NULL;
	}
	return panda_algorithm_ea_util_new();
}

const struct panda_algorithm_class panda_algorithm_ea_util_class = {
	.data_size = 0,
	.name = "ea_util",
	.create = from_string,
	.data_destroy = NULL,
	.overlap_probability = (PandaComputeOverlap) overlap_probability,
	.match_probability = (PandaComputeMatch) match_probability,
	.prob_unpaired = qual_nn_simple_bayesian,
};

PandaAlgorithm panda_algorithm_ea_util_new(
	void) {
	return panda_algorithm_new(&panda_algorithm_ea_util_class);
}
