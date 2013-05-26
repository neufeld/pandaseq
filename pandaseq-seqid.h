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

#ifndef _PANDASEQ_SEQID_H
#        define _PANDASEQ_SEQID_H
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
 * Reset a sequnce identifier.
 * @id: (out caller-allocates): The structure to clear.
 */
void panda_seqid_clear(
	panda_seq_identifier *id);

/**
 * Parse an Illumina header
 *
 * @id: (out caller-allocates): The structure to fill with the parse result.
 * Returns: The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
 */
int panda_seqid_parse(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy);

/**
 * Parse the Illumina header
 *
 * @id: (out caller-allocates): The structure to fill with the parse result.
 * Returns: The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
 * @old: (out): Whether the sequence is from CASAVA 1.3-1.5 or not.
 * @end_ptr: (out) (transfer none): The point in the input where parsing stopped. If parsing was successful, this will be the end of the string.
 * @see panda_seqid_parse
 */
int panda_seqid_parse_fail(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy,
	bool *old,
	const char **end_ptr);

/* === Methods === */

/**
 * Order two Illumina headers
 */
int panda_seqid_compare(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two);

/**
 * Compare two Illumina headers
 */
bool panda_seqid_equal(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two);

/**
 * Write an Illumina header for a sequence identifier to a file
 */
void panda_seqid_print(
	const panda_seq_identifier *id,
	FILE *file);

/**
 * Create an Illumina header for a sequence identifier
 * @id: (allow-none): The identifer to be formatted
 * Returns: (transfer none): Subsequent calls will obliterate the previously returned string.
 */
const char *panda_seqid_str(
	const panda_seq_identifier *id);

/**
 * Write the Illumina header to a printf-like function
 * @xprintf: (closure x): The callback to accept the input.
 */
void panda_seqid_xprint(
	const panda_seq_identifier *id,
	PandaPrintf xprintf,
	void *x);
EXTERN_C_END
#endif
