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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pandaseq-tablebuilder.h"
#include "prob.h"

struct panda_tbld {
	FILE *header;
	FILE *source;
};

PandaTBld panda_tbld_open(
	const char *base_name) {
	char buffer[1024];
	PandaTBld t_bld;
	FILE *header;
	FILE *source;
	size_t it;
	if (strlen(base_name) > 1000)
		return NULL;

	snprintf(buffer, 1024, "%s.c", base_name);
	source = fopen(buffer, "w");
	if (source == NULL) {
		perror(buffer);
		return NULL;
	}
	snprintf(buffer, 1024, "%s.h", base_name);
	header = fopen(buffer, "w");
	if (header == NULL) {
		perror(buffer);
		return NULL;
	}
	t_bld = malloc(sizeof(struct panda_tbld));
	t_bld->source = source;
	t_bld->header = header;
	for (it = 0; it < strlen(base_name); it++) {
		buffer[it] = isalnum(base_name[it]) ? toupper(base_name[it]) : '_';
	}
	buffer[it] = '\0';
	fprintf(header, "#ifndef _%s_H\n#define _%s_H\n", buffer, buffer);
	return t_bld;
}

void panda_tbld_free(
	PandaTBld t_bld) {
	fprintf(t_bld->header, "#endif\n");
	(void) fclose(t_bld->source);
	(void) fclose(t_bld->header);
	free(t_bld);
}

void panda_tbld_array(
	PandaTBld t_bld,
	const char *name,
	PandaArrayFormula formula,
	void *formula_context,
	size_t max) {
	size_t i;
	fprintf(t_bld->header, "extern const double %s[%zd];\n", name, max);
	fprintf(t_bld->source, "const double %s[%zd] = {\n", name, max);
	for (i = 0; i < max; i++) {
		if (i > 0) {
			fprintf(t_bld->source, ",");
		}
		fprintf(t_bld->source, " %g", formula(i, formula_context));
	}
	fprintf(t_bld->source, "};\n");
}

struct array_prob {
	PandaArrayProbFormula formula;
	void *formula_context;
	bool log_output;
};

static double array_prob_formula(
	size_t x,
	void *_data) {
	struct array_prob *data = (struct array_prob *) _data;
	double result = data->formula(PROBABILITY(x), data->formula_context);
	return data->log_output ? log(result) : result;
}

void panda_tbld_array_prob(
	PandaTBld t_bld,
	const char *name,
	PandaArrayProbFormula formula,
	void *formula_context,
	bool log_output) {
	struct array_prob context;

	context.formula = formula;
	context.formula_context = formula_context;
	context.log_output = log_output;

	panda_tbld_array(t_bld, name, array_prob_formula, &context, PHREDMAX + 1);
}

void panda_tbld_constant(
	PandaTBld t_bld,
	const char *name,
	double value) {
	fprintf(t_bld->header, "#define %s %g\n", name, value);
}

void panda_tbld_matrix(
	PandaTBld t_bld,
	const char *name,
	PandaMatrixFormula formula,
	void *formula_context,
	size_t x_max,
	size_t y_max) {
	size_t i, j;
	fprintf(t_bld->header, "extern const double %s[][%zd];\n", name, y_max);
	fprintf(t_bld->source, "const double %s[][%zd] = {\n", name, y_max);
	for (i = 0; i < x_max; i++) {
		if (i > 0) {
			fprintf(t_bld->source, ", \n");
		}
		fprintf(t_bld->source, "\t{");
		for (j = 0; j < y_max; j++) {
			if (j > 0) {
				fprintf(t_bld->source, ",");
			}

			fprintf(t_bld->source, " %g", formula(i, j, formula_context));

		}
		fprintf(t_bld->source, "}");
	}
	fprintf(t_bld->source, "};\n");
}

struct matrix_prob {
	PandaMatrixProbFormula formula;
	void *formula_context;
	bool log_output;
};

static double matrix_prob_formula(
	size_t x,
	size_t y,
	void *_data) {
	struct matrix_prob *data = (struct matrix_prob *) _data;
	double result = data->formula(PROBABILITY(x), PROBABILITY(y), data->formula_context);
	return data->log_output ? log(result) : result;
}

void panda_tbld_matrix_prob(
	PandaTBld t_bld,
	const char *name,
	PandaMatrixProbFormula formula,
	void *formula_context,
	bool log_output) {
	struct matrix_prob context;

	context.formula = formula;
	context.formula_context = formula_context;
	context.log_output = log_output;

	panda_tbld_matrix(t_bld, name, matrix_prob_formula, &context, PHREDMAX + 1, PHREDMAX + 1);
}
