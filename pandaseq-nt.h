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

#ifndef _PANDASEQ_NT_H
#        define _PANDASEQ_NT_H
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
EXTERN_C_BEGIN
/**
 * Nothing (invalid nucleotide)
 */
#        define PANDA_NT_Z ((panda_nt)0)
/**
 * Adenine
 */
#        define PANDA_NT_A ((panda_nt)1)
/**
 * Cytosine
 */
#        define PANDA_NT_C ((panda_nt)2)
/**
 * Guanine
 */
#        define PANDA_NT_G ((panda_nt)4)
/**
 * Thyamine
 */
#        define PANDA_NT_T ((panda_nt)8)
/**
 * Is nucleotide degenerate?
 */
#        define PANDA_NT_IS_DEGN(v) (((((unsigned int)(v)) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf) != 1)
/**
 * Is nucleotide all possible values?
 */
#        define PANDA_NT_IS_N(n) ((n) == (panda_nt)0x0F)
/**
 * Get the nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii(
	char c);
/**
 * Get the complement nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii_complement(
	char c);
/**
 * Get the complementary nucleotide.
 */
panda_nt panda_nt_complement(
	panda_nt nt);
/**
 * Convert a nucleotide to an IUPAC representation
 */
char panda_nt_to_ascii(
	panda_nt val);
EXTERN_C_END
#endif
