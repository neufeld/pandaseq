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
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include "pandaseq.h"
#include "prob.h"
#include "table.h"

#ifndef M_LN2
#        define M_LN2 0.69314718055994530942
#endif

#define CIRC(index, len) (((index) + (len)) % (len))

/* Compute 1-exp(p) See <http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf> */
double panda_log1mexp(
	double p) {
	return (p > M_LN2) ? log1p(-exp(-p)) : log(-expm1(-p));
}

typedef void (
	*base_score) (
	const void *data,
	panda_nt *base,
	double *prob,
	double *notprob);

static size_t computeoffset(
	double threshold,
	double penalty,
	bool reverse,
	const unsigned char *seq,
	size_t seq_length,
	size_t size,
	base_score score,
	const panda_nt *primer,
	size_t primerlen) {
	/* Circular buffer of probabilities of primer alignment indexed by the offset. */
	double probabilities[primerlen];
	double bestpr = exp(primerlen * threshold);
	size_t bestindex = 0;
	size_t index;
	if (primerlen > seq_length) {
		return 0;
	}

	for (index = 0; index < primerlen; index++) {
		probabilities[index] = -INFINITY;
	}

	for (index = 0; index < seq_length; index++) {
		ptrdiff_t x;
		double last_pr = exp(probabilities[CIRC(index, primerlen)] / (index + 1)) - index * penalty;
		/* The last bucket in the buffer holds the probability of a complete alignment. If it so better than we have seen previously, store it. */
		if (last_pr > bestpr) {
			bestpr = last_pr;
			bestindex = index + 1;
		}
		probabilities[CIRC(index, primerlen)] = 0;
		for (x = (ptrdiff_t) (primerlen > index ? index : primerlen - 1); x >= 0; x--) {
			if (!PANDA_NT_IS_N(primer[x])) {
				panda_nt nt;
				double p;
				double notp;
				score(&seq[size * (reverse ? (seq_length - index - 1) : index)], &nt, &p, &notp);
				probabilities[CIRC(index - x, primerlen)] += ((nt & primer[x]) != 0) ? p : notp;
			}
		}
	}
	return bestindex;
}

void qual_base_score(
	const void *data,
	panda_nt *base,
	double *prob,
	double *notprob) {
	int phred = PHREDCLAMP(((panda_qual *) data)->qual);
	*base = ((panda_qual *) data)->nt;
	*prob = qual_score[phred];
	*notprob = qual_score_err[phred];
}

size_t panda_compute_offset_qual(
	double threshold,
	double penalty,
	bool reverse,
	const panda_qual *haystack,
	size_t haystack_length,
	const panda_nt *needle,
	size_t needle_length) {
	return computeoffset(threshold, penalty, reverse, (const unsigned char *) haystack, haystack_length, sizeof(panda_qual), qual_base_score, needle, needle_length);
}

void result_base_score(
	const void *data,
	panda_nt *base,
	double *prob,
	double *notprob) {
	*base = ((panda_result *) data)->nt;
	*prob = ((panda_result *) data)->p;
	*notprob = panda_log1mexp(*prob);
}

size_t panda_compute_offset_result(
	double threshold,
	double penalty,
	bool reverse,
	const panda_result *haystack,
	size_t haystack_length,
	const panda_nt *needle,
	size_t needle_length) {
	return computeoffset(threshold, penalty, reverse, (const unsigned char *) haystack, haystack_length, sizeof(panda_result), result_base_score, needle, needle_length);
}
