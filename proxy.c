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
#include "config.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

struct panda_log_proxy {
	size_t refcnt;
#ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#endif
	PandaWriter writer;
};

PandaLogProxy panda_log_proxy_new(
	PandaWriter writer) {
	PandaLogProxy proxy = malloc(sizeof(struct panda_log_proxy));
	proxy->refcnt = 1;
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&proxy->mutex, NULL);
#endif
	proxy->writer = panda_writer_ref(writer);
	return proxy;
}

PandaLogProxy panda_log_proxy_new_stderr(
	) {
	PandaWriter writer = panda_writer_new_stderr();
	PandaLogProxy proxy = panda_log_proxy_new(writer);
	panda_writer_unref(writer);
	return proxy;
}

PandaLogProxy panda_log_proxy_new_file(
	FILE *file) {
	PandaWriter writer = panda_writer_new_file(file);
	PandaLogProxy proxy = panda_log_proxy_new(writer);
	panda_writer_unref(writer);
	return proxy;
}

PandaLogProxy panda_log_proxy_open_file(
	const char *filename,
	bool bzip) {
	PandaWriter writer = panda_writer_open_file(filename, bzip);
	PandaLogProxy proxy = (writer == NULL) ? NULL : panda_log_proxy_new(writer);
	panda_writer_unref(writer);
	return proxy;
}

PandaLogProxy panda_log_proxy_ref(
	PandaLogProxy proxy) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	proxy->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
	return proxy;
}

void panda_log_proxy_unref(
	PandaLogProxy proxy) {
	size_t count;
	if (proxy == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	count = --(proxy->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
	if (count == 0) {
		panda_writer_unref(proxy->writer);
#ifdef HAVE_PTHREAD
		pthread_mutex_destroy(&proxy->mutex);
#endif
		free(proxy);
	}
}

void panda_log_proxy_perror(
	PandaLogProxy proxy,
	const char *prefix) {
	const char *message = strerror(errno);
	if (prefix == NULL) {
		panda_writer_append(proxy->writer, "%s\n", message);
	} else {
		panda_writer_append(proxy->writer, "%s: %s\n", prefix, message);
	}
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const char *message) {
	const char *name;

	name = panda_assembler_get_name(assembler);
	panda_writer_append(proxy->writer, "%s%s%s%s", name == NULL ? "" : name, name == NULL ? "" : "\t", panda_code_str(code), id == NULL ? "" : "\t");
	if (id != NULL) {
		panda_writer_append_id(proxy->writer, id);
	}
	if (message != NULL) {
		panda_writer_append(proxy->writer, "\t%s\n", message);
	} else {
		panda_writer_append(proxy->writer, "\n");
	}
	if (code == PANDA_CODE_ID_PARSE_FAILURE && message != NULL) {
		panda_writer_append(proxy->writer, "* * * * * Something is wrong with this ID. If tags are absent, try passing the -B option.\n* * * * * Consult `pandaseq-checkid \"%s\"` to get an idea of the problem..\n", message);
	} else if (code == PANDA_CODE_PHRED_OFFSET) {
		panda_writer_append(proxy->writer, "* * * * * Using the default PHRED+33 offset, but no sequences had quality data under PHRED+64.\n* * * * * This is probably not what you want. Consult the manual about the -6 option.\n");
	} else if (code == PANDA_CODE_READ_TOO_LONG) {
		panda_writer_append(proxy->writer, "* * * * * The input reads are longer than this version of PANDAseq can handle. Currently %zd nucleotides.\n", PANDA_MAX_LEN);
	}
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_write_str(
	PandaLogProxy proxy,
	const char *str) {
	panda_writer_append(proxy->writer, "%s\n", str);
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_write_f(
	PandaLogProxy proxy,
	const char *format,
	...) {
	va_list va;
	va_start(va, format);
	panda_writer_append_v(proxy->writer, format, va);
	va_end(va);
	panda_writer_commit(proxy->writer);
}

static void write_assembler_name(
	PandaLogProxy proxy,
	PandaAssembler assembler) {
	const char *name;
	name = panda_assembler_get_name(assembler);

	if (name != NULL) {
		panda_writer_append(proxy->writer, "%s\t", name);
	}
}

void panda_log_proxy_write_overlap(
	PandaLogProxy proxy,
	PandaAssembler assembler) {
	size_t it = 0;
	size_t max;

	write_assembler_name(proxy, assembler);

	panda_writer_append(proxy->writer, "STAT\tOVERLAPS\t%ld", panda_assembler_get_overlap_count(assembler, it));
	max = panda_assembler_get_longest_overlap(assembler);
	for (it = 1; it <= max; it++) {
		panda_writer_append(proxy->writer, " %ld", panda_assembler_get_overlap_count(assembler, it));
	}
	panda_writer_append_c(proxy->writer, '\n');
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_stat_double(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	double value) {
	write_assembler_name(proxy, assembler);
	panda_writer_append(proxy->writer, "STAT\t%s\t%f\n", name, value);
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_stat_long(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	long value) {
	write_assembler_name(proxy, assembler);
	panda_writer_append(proxy->writer, "STAT\t%s\t%ld\n", name, value);
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_stat_size_t(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	size_t value) {
	write_assembler_name(proxy, assembler);
	panda_writer_append(proxy->writer, "STAT\t%s\t%zd\n", name, value);
	panda_writer_commit(proxy->writer);
}

void panda_log_proxy_stat_str(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	const char *value) {
	write_assembler_name(proxy, assembler);
	panda_writer_append(proxy->writer, "STAT\t%s\t%s\n", name, value);
	panda_writer_commit(proxy->writer);
}

PandaWriter panda_log_proxy_get_writer(
	PandaLogProxy proxy) {
	return proxy->writer;
}
