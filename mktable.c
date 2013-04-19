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

#include <stdio.h>
#include <math.h>
#include "prob.h"

static void buildmatrix(
	FILE *header,
	FILE *source,
	char *name,
	double (*formula) (double,
		double)) {
	int i, j;
	fprintf(header, "extern double %s[][%d];\n", name, PHREDMAX + 1);
	fprintf(source, "double %s[][%d] = {\n", name, PHREDMAX + 1);
	for (i = 0; i <= PHREDMAX; i++) {
		if (i > 0) {
			fprintf(source, ", \n");
		}
		fprintf(source, "\t{");
		for (j = 0; j <= PHREDMAX; j++) {
			if (j > 0) {
				fprintf(source, ",");
			}

			fprintf(source, " %g", log(formula(PROBABILITY(i), PROBABILITY(j))));

		}
		fprintf(source, "}");
	}
	fprintf(source, "};\n");

}

static void buildlist(
	FILE *header,
	FILE *source,
	char *name,
	double (*formula) (double)) {
	int i;
	fprintf(header, "extern double %s[%d];\n", name, PHREDMAX + 1);
	fprintf(source, "double %s[%d] = {\n", name, PHREDMAX + 1);
	for (i = 0; i <= PHREDMAX; i++) {
		if (i > 0) {
			fprintf(source, ",");
		}
		fprintf(source, " %g", (double) formula(PROBABILITY(i)));
	}
	fprintf(source, "};\n");
}

static double match(
	double p,
	double q) {
	return (1 - p) * (1 - q) + p * q / 3;
}

static double mismatch(
	double p,
	double q) {
	return (1 - p) * q / 3 + (1 - q) * p / 3 + 2 * p * q / 9;
}

static double nmatch(
	double p) {
	if (p == 1) {
		return -2;
	}
	return log(p / 2 + 0.25);
}

static double score(
	double p) {
	if (p == 1) {
		return -2;
	}
	return log(1.0 - p);
}

static double score_err(
	double p) {
	return log(p);
}

int main(
	int argc,
	char **argv) {
	FILE *header;
	FILE *source = fopen("table.c", "w");
	if (source == NULL) {
		perror("table.c");
		return 1;
	}
	header = fopen("table.h", "w");
	if (header == NULL) {
		perror("table.h");
		return 1;
	}

	fprintf(header, "#ifndef _TABLE_H\n#define _TABLE_H\n#define qual_nn %f\n", log(0.25));
	buildmatrix(header, source, "qual_match", match);
	buildmatrix(header, source, "qual_mismatch", mismatch);
	buildlist(header, source, "qual_nmatch", nmatch);
	buildlist(header, source, "qual_score", score);
	buildlist(header, source, "qual_score_err", score_err);
	fprintf(header, "#endif\n");

	(void) fclose(source);
	(void) fclose(header);
	return 0;
}
