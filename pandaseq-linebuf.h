/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2013  Andre Masella

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

#ifndef _PANDASEQ_LINEBUF_H
#        define _PANDASEQ_LINEBUF_H
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
 * Create a new line reader from a buffer reading source.
 * @read: (closure read_data) (scope notified): the function to do reading.
 */
PandaLineBuf panda_linebuf_new(
	PandaBufferRead read,
	void *read_data,
	PandaDestroy read_destroy);

/* === Methods === */
/**
 * Destroy the line buffer.
 */
void panda_linebuf_free(
	PandaLineBuf linebuf);
/**
 * Read the next line.
 * Returns: (transfer none) (allow-none): the next line in the file. This is only valid until the next call.
 */
const char *panda_linebuf_next(
	PandaLineBuf linebuf);
EXTERN_C_END
#endif
