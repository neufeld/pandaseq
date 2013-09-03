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

#ifndef _PANDASEQ_TBLD_H
#        define _PANDASEQ_TBLD_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <stdarg.h>
#        include <stdio.h>
#        include <stdbool.h>
EXTERN_C_BEGIN typedef struct panda_tbld *PandaTBld;
/**
 * A formula to compute over a range [0, max) of type size_t.
 * @context: (closure): the provided user context.
 */
typedef double (
	*PandaArrayFormula) (
	size_t i,
	void *context);
/**
 * A formula to compute over a range P_i : i in [0, PHRED_MAX], where P_i is a PHRED score converted to a probability.
 * @context: (closure): the provided user context.
 */
typedef double (
	*PandaArrayProbFormula) (
	double p_i,
	void *context);
/**
 * A formula to compute over a matrix of size_t pairs.
 * @context: (closure): the provided user context.
 * @see PandaArrayFormula
 */
typedef double (
	*PandaMatrixFormula) (
	size_t x,
	size_t y,
	void *context);
/**
 * A formula to compute over a matrix of probability pairs.
 * @context: (closure): the provided user context.
 * @see PandaArrayProbFormula
 */
typedef double (
	*PandaMatrixProbFormula) (
	double p_x,
	double p_y,
	void *context);

/* === Constructors === */
/**
 * Create a new table written to a file.
 *
 * The file name and header symbols will be inferred from the provided name.
 * Returns: (allow-none): the builder, unless an error occurred.
 */
PandaTBld panda_tbld_open(
	const char *base_name);

/* === Methods === */
void panda_tbld_free(
	PandaTBld t_bld);
/**
 * Write an array filling it with the provided formula.
 * @name: the C symbol for the formula.
 * @formula: (closure formula_context): the formula to compute.
 * @max: the array length of the output.
 */
void panda_tbld_array(
	PandaTBld t_bld,
	const char *name,
	PandaArrayFormula formula,
	void *formula_context,
	size_t max);
/**
 * Write an array filling it with the provided formula, computed over PHRED probabilities.
 * @name: the C symbol for the formula.
 * @formula: (closure formula_context): the formula to compute.
 * @log_output: write the logarithm of the output, rather than the output.
 */
void panda_tbld_array_prob(
	PandaTBld t_bld,
	const char *name,
	PandaArrayProbFormula formula,
	void *formula_context,
	bool log_output);
/**
 * Write `#define` constant.
 * @name: the C symbol for the constant.
 * @value: the value to write
 */
void panda_tbld_constant(
	PandaTBld t_bld,
	const char *name,
	double value);
/**
 * Write an array of arrays filling it with the provided formula.
 * @name: the C symbol for the formula.
 * @formula: (closure formula_context): the formula to compute.
 * @x_max: the outer array length of the output.
 * @y_max: the inner array length of the output.
 */
void panda_tbld_matrix(
	PandaTBld t_bld,
	const char *name,
	PandaMatrixFormula formula,
	void *formula_context,
	size_t x_max,
	size_t y_max);
/**
 * Write an array of arrays filling it with the provided formula, computed over PHRED probabilities.
 * @name: the C symbol for the formula.
 * @formula: (closure formula_context): the formula to compute.
 * @log_output: write the logarithm of the output, rather than the output.
 */
void panda_tbld_matrix_prob(
	PandaTBld t_bld,
	const char *name,
	PandaMatrixProbFormula formula,
	void *formula_context,
	bool log_output);
EXTERN_C_END
#endif
