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
/**
 * Write errors and information to a file.
 */
bool panda_logger_file(
	PandaCode code,
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const char *message,
	FILE *file);

/* === Constructors === */

/**
 * Create a new proxy with a callback.
 */
PandaLogProxy panda_log_proxy_new(
	PandaLogger log,
	void *log_data,
	PandaDestroy log_destroy);

/**
 * Create a new proxy to standard error.
 */
PandaLogProxy panda_log_proxy_new_stderr(
	);

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

bool panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const char *message);
EXTERN_C_END
#endif
