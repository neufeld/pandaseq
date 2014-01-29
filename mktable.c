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

#include <float.h>
#include <math.h>
#include "pandaseq-tablebuilder.h"

static double match(
	double p,
	double q,
	void *data) {
	return (1 - p) * (1 - q) + p * q / 3;
}

static double mismatch(
	double p,
	double q,
	void *data) {
	return (1 - p) * q / 3 + (1 - q) * p / 3 + 2 * p * q / 9;
}

static double match_pear(
	double p,
	double q,
	void *data) {
	return (1 - (1 - q) * p / 3 - (1 - p) * q / 3 - 2 * (1 - p) * (1 - q) / 9);
}

static double mismatch_pear(
	double p,
	double q,
	void *data) {
	return (1 - p) * q / 3 + (1 - q) * p / 3 + p * q / 2;
}

static double score(
	double p,
	void *data) {
	if (p == 1) {
		return -2;
	}
	return log(1.0 - p);
}

static double score_err(
	double p,
	void *data) {
	return log(p);
}

double mismatch_rdp(
	double p,
	double q,
	void *data) {
	return ((1 - p) * q / 3 + (1 - q) * p / 3 + 2 * p * q / 9);
}

double mismatch_rdp_assembled(
	double p,
	double q,
	void *data) {

	double min = p <= q ? p : q;
	return (min - p * q / 3.0) / (p + q - 4.0 / 3.0 * p * q);
}

int main(
	int argc,
	char **argv) {
	PandaTBld t_bld;

	t_bld = panda_tbld_open("table");
	if (t_bld == NULL) {
		return 1;
	}

	panda_tbld_constant(t_bld, "qual_nn_simple_bayesian", log(0.25));
	panda_tbld_matrix_prob(t_bld, "qual_match_simple_bayesian", match, NULL, true);
	panda_tbld_matrix_prob(t_bld, "qual_mismatch_simple_bayesian", mismatch, NULL, true);
	panda_tbld_matrix_prob(t_bld, "qual_match_pear", match_pear, NULL, true);
	panda_tbld_matrix_prob(t_bld, "qual_mismatch_pear", mismatch_pear, NULL, true);
	panda_tbld_matrix_prob(t_bld, "qual_mismatch_rdp_mle", mismatch_rdp, NULL, true);
	panda_tbld_matrix_prob(t_bld, "qual_mismatch_assembled_rdp_mle", mismatch_rdp_assembled, NULL, true);
	panda_tbld_array_prob(t_bld, "qual_score", score, NULL, false);
	panda_tbld_array_prob(t_bld, "qual_score_err", score_err, NULL, false);

	panda_tbld_free(t_bld);
	return 0;
}
