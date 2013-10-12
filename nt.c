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
#include "config.h"
#include <math.h>
#include <stdlib.h>
#include "pandaseq.h"
#include "nt.h"
#include "prob.h"
#include "table.h"

static char ntchar[16] = { 'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' };

panda_nt iupac_forward[32] = {
	/*@ */ PANDA_NT_Z,
	 /*A*/ PANDA_NT_A,
	 /*B*/ PANDA_NT_C | PANDA_NT_G | PANDA_NT_T,
	 /*C*/ PANDA_NT_C,
	 /*D*/ PANDA_NT_A | PANDA_NT_G | PANDA_NT_T,
	 /*E*/ PANDA_NT_Z,
	 /*F*/ PANDA_NT_Z,
	 /*G*/ PANDA_NT_G,
	 /*H*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_T,
	 /*I*/ PANDA_NT_Z,
	 /*J*/ PANDA_NT_Z,
	 /*K*/ PANDA_NT_G | PANDA_NT_T,
	 /*L*/ PANDA_NT_Z,
	 /*M*/ PANDA_NT_A | PANDA_NT_C,
	 /*N*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T,
	 /*O*/ PANDA_NT_Z,
	 /*P*/ PANDA_NT_Z,
	 /*Q*/ PANDA_NT_Z,
	 /*R*/ PANDA_NT_A | PANDA_NT_G,
	 /*S*/ PANDA_NT_C | PANDA_NT_G,
	 /*T*/ PANDA_NT_T,
	 /*U*/ PANDA_NT_T,
	 /*V*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G,
	 /*W*/ PANDA_NT_A | PANDA_NT_T,
	 /*X*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T,
	 /*Y*/ PANDA_NT_C | PANDA_NT_T,
	 /*Z*/ PANDA_NT_Z,
	/*[ */ PANDA_NT_Z,
	/*\ */ PANDA_NT_Z,
	/*] */ PANDA_NT_Z,
	/*^ */ PANDA_NT_Z,
	/*_*/ PANDA_NT_Z
};

panda_nt iupac_reverse[32] = {
/* @ */ PANDA_NT_Z,
	 /*A*/ PANDA_NT_T,
	 /*B*/ PANDA_NT_G | PANDA_NT_C | PANDA_NT_A,
	 /*C*/ PANDA_NT_G,
	 /*D*/ PANDA_NT_T | PANDA_NT_C | PANDA_NT_A,
	 /*E*/ PANDA_NT_Z,
	 /*F*/ PANDA_NT_Z,
	 /*G*/ PANDA_NT_C,
	 /*H*/ PANDA_NT_T | PANDA_NT_G | PANDA_NT_A,
	 /*I*/ PANDA_NT_Z,
	 /*J*/ PANDA_NT_Z,
	 /*K*/ PANDA_NT_C | PANDA_NT_A,
	 /*L*/ PANDA_NT_Z,
	 /*M*/ PANDA_NT_T | PANDA_NT_G,
	 /*N*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T,
	 /*O*/ PANDA_NT_Z,
	 /*P*/ PANDA_NT_Z,
	 /*Q*/ PANDA_NT_Z,
	 /*R*/ PANDA_NT_T | PANDA_NT_C,
	 /*S*/ PANDA_NT_G | PANDA_NT_C,
	 /*T*/ PANDA_NT_A,
	 /*U*/ PANDA_NT_A,
	 /*V*/ PANDA_NT_T | PANDA_NT_G | PANDA_NT_C,
	 /*W*/ PANDA_NT_T | PANDA_NT_A,
	 /*X*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T,
	 /*Y*/ PANDA_NT_G | PANDA_NT_A,
	 /*Z*/ PANDA_NT_Z,
	/*[ */ PANDA_NT_Z,
	/*\ */ PANDA_NT_Z,
	/*] */ PANDA_NT_Z,
	/*^ */ PANDA_NT_Z,
	/*_*/ PANDA_NT_Z
};

double panda_quality_probability(
	const panda_qual *q) {
	return exp(panda_quality_log_probability(q));
}

double panda_quality_log_probability(
	const panda_qual *q) {
	return qual_score[PHREDCLAMP(q->qual)];
}

double panda_quality_compare(
	const panda_qual *a,
	const panda_qual *b) {
	return ((a->nt & b->nt) != '\0' ? qual_match : qual_mismatch)[PHREDCLAMP(a->qual)][PHREDCLAMP(b->qual)];
}

static int find_double(
	const double *key,
	const double *value) {
	if (*key == *value)
		return 0;
	return *key < *value ? -1 : 1;
}

char panda_result_phred(
	const panda_result *r) {
	const double *item;

	if (r->p <= -2)
		return 1;

	item = bsearch(&r->p, qual_score, PHREDMAX, sizeof(double), (int (*)(const void *, const void *)) find_double);
	if (item == NULL)
		return 1;
	return item - qual_score;
}

panda_nt panda_nt_from_ascii(
	char c) {
	return iupac_forward[(int) c & 0x1F];
}

panda_nt panda_nt_from_ascii_complement(
	char c) {
	return iupac_reverse[(int) c & 0x1F];
}

char panda_nt_to_ascii(
	panda_nt val) {
	if (val < PANDA_NT_Z || val > (panda_nt) 15) {
		return 'N';
	}
	return ntchar[(int) (val)];
}
