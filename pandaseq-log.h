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

#ifndef _PANDASEQ_LOG_H
#        define _PANDASEQ_LOG_H
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
/* === Constructors === */
/**
 * Write the log to a transactional writer.
 */
PandaLogProxy panda_log_proxy_new(
	PandaWriter writer);

/**
 * Create a new proxy to standard error.
 */
PandaLogProxy panda_log_proxy_new_stderr(
	);

/**
 * Write the log to an open file.
 */
PandaLogProxy panda_log_proxy_new_file(
	FILE *file);

/**
 * Open a file for writing error messages.
 * @filename: The file to write.
 * @bzip: Write BZipped text rather than plain text.
 * Returns: (allow-none): A logger proxy.
 */
PandaLogProxy panda_log_proxy_open_file(
	const char *filename,
	bool bzip);

/* === Getters and Setters === */

PandaWriter panda_log_proxy_get_writer(
	PandaLogProxy proxy);

/* === Methods === */

/**
 * Increase the reference count on a proxy.
 */
PandaLogProxy panda_log_proxy_ref(
	PandaLogProxy proxy);

/**
 * Decrease the reference count on a proxy.
 * @proxy: (transfer full): the proxy to release.
 */
void panda_log_proxy_unref(
	PandaLogProxy proxy);

/**
 * Writes an error message to the log with the same semantics as the POSIX perror function.
 */
void panda_log_proxy_perror(
	PandaLogProxy proxy,
	const char *prefix);

/**
 * Print a message to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const char *message);

/**
 * Print the overlap histogram of an assember to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_write_overlap(
	PandaLogProxy proxy,
	PandaAssembler assembler);

/**
 * Put a printf-like message in the log.
 */
void panda_log_proxy_write_f(
	PandaLogProxy proxy,
	const char *format,
	...);

/**
 * Print a string to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_write_str(
	PandaLogProxy proxy,
	const char *str);

/**
 * Print a double with a STAT header to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_stat_double(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	double value);

/**
 * Print a double with a STAT header to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_stat_long(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	long value);

/**
 * Print a size_t with a STAT header to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_stat_size_t(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	size_t value);

/**
 * Print a string with a STAT header to the log.
 *
 * This method is thread-safe.
 */
void panda_log_proxy_stat_str(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	const char *value);

EXTERN_C_END
#endif
