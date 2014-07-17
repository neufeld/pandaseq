/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella
     Copyright (C) 2013       Qiong Wang

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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pandaseq.h"
#include "algo.h"
#include "prob.h"
#include "table.h"

static double match_probability(
	void *data,
	bool match,
	char a,
	char b) {
	(void) data;
	if (match) {
		char max = (a >= b) ? a : b;
		return qual_score[PHREDCLAMP(max)];
	} else {
		return qual_mismatch_assembled_rdp_mle[PHREDCLAMP(a)][PHREDCLAMP(b)];
	}
}

static double overlap_probability(
	void *data,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	size_t overlap) {

	double probability = 0;
	size_t i;

	(void) data;
	for (i = 0; i < overlap; i++) {
		int findex = forward_length + i - overlap;
		int rindex = reverse_length - i - 1;
		if (findex < 0 || rindex < 0 || (size_t) findex > forward_length || (size_t) rindex > reverse_length)
			continue;
		panda_nt f = forward[findex].nt;
		panda_nt r = reverse[rindex].nt;
		char fqual = PHREDCLAMP(forward[findex].qual);
		char rqual = PHREDCLAMP(reverse[rindex].qual);
		bool ismatch = ((f & r) != 0);

		/* when two bases match, the assumption that the forward and reverse bases are from independent observations doesn't work with the MiSeq mock community data we tested. Instead, the higher score of the two raw base q scores is close to the predicated error rate */
		if (ismatch) {
			probability += qual_match_simple_bayesian[(int) fqual][(int) rqual] - qual_nn_simple_bayesian;
		} else {
			probability += qual_mismatch_rdp_mle[(int) fqual][(int) rqual] - qual_nn_simple_bayesian;
		}
	}
	return probability;
}

static PandaAlgorithm from_string(
	const char *argument) {

	if (argument == NULL || strlen(argument) == 0)
		return panda_algorithm_rdp_mle_new();
	return NULL;
}

const struct panda_algorithm_class panda_algorithm_rdp_mle_class = {
	.data_size = 0,
	.name = "rdp_mle",
	.create = from_string,
	.data_destroy = NULL,
	.overlap_probability = (PandaComputeOverlap) overlap_probability,
	.match_probability = (PandaComputeMatch) match_probability,
	.prob_unpaired = qual_nn_simple_bayesian,
};

PandaAlgorithm panda_algorithm_rdp_mle_new(
	void) {
	PandaAlgorithm algo = panda_algorithm_new(&panda_algorithm_rdp_mle_class);
	return algo;
}
