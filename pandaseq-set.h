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

#ifndef _PANDASEQ_SET_H
#        define _PANDASEQ_SET_H
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
/* === Constructor === */
/**
 * Create a new, empty set.
 */
PandaSet panda_idset_new(
	void);

/* === Methods === */

/**
 * Add a sequence identifier to a set.
 */
void panda_idset_add(
	PandaSet set,
	const panda_seq_identifier *id);

/**
 * Parse a sequence identifier and add it to the set.
 * @id: the text id to parse
 * @old: (out): Whether the sequence is from CASAVA 1.3-1.5 or not.
 * @end_ptr: (out) (transfer none): The point in the input where parsing stopped. If parsing was successful, this will be the end of the string.
 * Returns: true on success
 * @see panda_seqid_parse_fail
 */
bool panda_idset_add_str(
	PandaSet set,
	const char *id,
	PandaTagging policy,
	bool *old,
	const char **end_ptr);

/**
 * Check if a sequence identifier has been added to the set.
 */
bool panda_idset_contains(
	PandaSet set,
	const panda_seq_identifier *id);

/**
 * Increase the reference count on a set.
 *
 * This is thread-safe.
 */
PandaSet panda_idset_ref(
	PandaSet set);

/**
 * Decrease the reference count on a set.
 *
 * This is thread-safe.
 * @set: (transfer full): the set to be released.
 */
void panda_idset_unref(
	PandaSet set);
EXTERN_C_END
#endif
