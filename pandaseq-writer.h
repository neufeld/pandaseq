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

#ifndef _PANDASEQ_WRITER_H
#        define _PANDASEQ_WRITER_H
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
 * Create a new writer, backed by some target.
 */
PandaWriter panda_writer_new(
	PandaBufferWrite write,
	void *write_data,
	PandaDestroy write_destroy);
/**
 * Create a new writer, backed by an open file.
 * @file: (transfer full): the open file.
 */
PandaWriter panda_writer_new_file(
	FILE *file);
/**
 * Create a new writer, backed by standard error.
 */
PandaWriter panda_writer_new_stderr(
	void);
/**
 * Create a new writer, backed by standard output.
 */
PandaWriter panda_writer_new_stdout(
	void);

/**
 * Open a file for writing.
 * @filename: The file to write.
 * @bzip: Write BZipped text rather than plain text.
 * Returns: (allow-none): A writer.
 */
PandaWriter panda_writer_open_file(
	const char *filename,
	bool bzip);

/* === Methods === */
/**
 * Write a printf-like formatted string to the output.
 */
void panda_writer_append(
	PandaWriter writer,
	const char *format,
	...);
/**
 * Write a single character to the output.
 */
void panda_writer_append_c(
	PandaWriter writer,
	char c);
/**
 * Write a sequence identifier to the output.
 */
void panda_writer_append_id(
	PandaWriter writer,
	const panda_seq_identifier *id);
/**
 * Write a printf-like formatted string to the output.
 */
void panda_writer_append_v(
	PandaWriter writer,
	const char *format,
	va_list va);
/**
 * End the current transaction and start another.
 *
 * This will consider the appending done so far to this writer to be a unit
 * that can be passed to the output when necessary.
 */
void panda_writer_commit(
	PandaWriter writer);

/**
 * A target write to be commited at the same time as this one.
 *
 * A writer may have a slave writer that will recieve a commit whenever panda_writer_commit is called on this writer.
 */
void panda_writer_set_slave(
	PandaWriter writer,
	PandaWriter slave);
PandaWriter panda_writer_get_slave(
	PandaWriter writer);

/**
 * Force writing all buffered data to the output.
 *
 * This requires getting a lock and happens automatically under normal circumstances.
 */
void panda_writer_flush(
	PandaWriter writer);

/**
 * Increase the reference count on a writer.
 *
 * This is thread-safe.
 */
PandaWriter panda_writer_ref(
	PandaWriter writer);

/**
 * Decrease the reference count on a writer.
 *
 * This is thread-safe.
 * @set: (transfer full): the set to be released.
 */
void panda_writer_unref(
	PandaWriter writer);
EXTERN_C_END
#endif
