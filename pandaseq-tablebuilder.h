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
typedef double (
	*PandaArrayFormula) (
	size_t i,
	void *context);
typedef double (
	*PandaArrayProbFormula) (
	double p_i,
	void *context);
typedef double (
	*PandaMatrixFormula) (
	size_t x,
	size_t y,
	void *context);
typedef double (
	*PandaMatrixProbFormula) (
	double p_x,
	double p_y,
	void *context);

/* === Constructors === */
PandaTBld panda_tbld_open(
	const char *base_name);

/* === Methods === */
void panda_tbld_free(
	PandaTBld t_bld);
void panda_tbld_array(
	PandaTBld t_bld,
	char *name,
	PandaArrayFormula formula,
	void *formula_context,
	size_t max);
void panda_tbld_array_prob(
	PandaTBld t_bld,
	char *name,
	PandaArrayProbFormula formula,
	void *formula_context,
	bool log_output);
void panda_tbld_constant(
	PandaTBld t_bld,
	char *name,
	double value);
void panda_tbld_matrix(
	PandaTBld t_bld,
	char *name,
	PandaMatrixFormula formula,
	void *formula_context,
	size_t x_max,
	size_t y_max);
void panda_tbld_matrix_prob(
	PandaTBld t_bld,
	char *name,
	PandaMatrixProbFormula formula,
	void *formula_context,
	bool log_output);
EXTERN_C_END
#endif
