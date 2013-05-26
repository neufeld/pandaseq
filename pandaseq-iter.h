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

#ifndef _PANDASEQ_ITER_H
#        define _PANDASEQ_ITER_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/* === Constructors === */
/**
 * Create an iterator over a sequence of nucleotides.
 * @seq: (array length=seq_length) (scope container): the sequence to iterate over, and its length. This sequence must not be freed during the life of the iterator.
 * @reverse: true to iterate from the end of the sequence rather than the beginning
 * @k: the length of the output words. This must range between 1 and 4 * sizeof(size_t). Any other values will be converted to the standard k-mer length of 8.
 */
PandaIter panda_iterate_nt(
	panda_nt *seq,
	size_t seq_length,
	bool reverse,
	int k);

/**
 * Iterate over quality-annotated sequence.
 * @see panda_iterate_nt
 */
PandaIter panda_iterate_qual(
	panda_qual *seq,
	size_t seq_length,
	bool reverse,
	int k);

/**
 * Iterate over probability-annotated sequence.
 * @see panda_iterate_nt
 */
PandaIter panda_iterate_result(
	panda_result *seq,
	size_t seq_length,
	bool reverse,
	int k);

/* === Methods === */

/**
 * Copy an iterator to a new one, preserving its current state.
 */
PandaIter panda_iter_dup(
	PandaIter iter);

/**
 * Destroy an iterator.
 */
void panda_iter_free(
	PandaIter iter);

/**
 * Advance to the next position in the sequence.
 * Returns: (allow-none) (transfer none): if null, there are no more k-mers in the sequence
 */
const panda_kmer *panda_iter_next(
	PandaIter iter);

/**
 * Set an iterator back to the beginning of the sequence.
 */
void panda_iter_reset(
	PandaIter iter);

/* === Getters and Setters === */

/**
 * Get the number of useful bits in the output.
 *
 * This is the maximum value of panda_kmer.kmer for this iterator.
 */
size_t panda_iter_bits(
	PandaIter iter);

/**
 * Get the k-mer length for the iterator.
 */
int panda_iter_k(
	PandaIter iter);

EXTERN_C_END
#endif
