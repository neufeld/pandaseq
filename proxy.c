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
#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

struct panda_log_proxy {
	size_t refcnt;
	 MANAGED_MEMBER(
		PandaPrintf,
		printf);
#ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#endif
};

PandaLogProxy panda_log_proxy_new(
	PandaPrintf printf,
	void *printf_data,
	PandaDestroy printf_destroy) {
	PandaLogProxy proxy = malloc(sizeof(struct panda_log_proxy));
	proxy->refcnt = 1;
	proxy->printf = printf;
	proxy->printf_data = printf_data;
	proxy->printf_destroy = printf_destroy;
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&proxy->mutex, NULL);
#endif
	return proxy;
}

PandaLogProxy panda_log_proxy_new_stderr(
	) {
	return panda_log_proxy_new((PandaPrintf) fprintf, stderr, NULL);
}

PandaLogProxy panda_log_proxy_new_file(
	FILE *file) {
	return panda_log_proxy_new((PandaPrintf) fprintf, file, (PandaDestroy) fclose);
}

#define BSIZE 2048
void bzprintf(
	BZFILE * file,
	const char *format,
	...) {
	char buf[BSIZE];
	int len;
	va_list va;
	va_start(va, format);
	len = vsnprintf(buf, BSIZE, format, va);
	va_end(va);
	(void) BZ2_bzwrite(file, buf, len);
}

PandaLogProxy panda_log_proxy_open_file(
	const char *filename,
	bool bzip) {
	if (bzip) {
		BZFILE *file;
		file = BZ2_bzopen(filename, "w");
		return (file == NULL) ? NULL : panda_log_proxy_new((PandaPrintf) bzprintf, file, (PandaDestroy) BZ2_bzclose);
	} else {
		FILE *file;
		file = fopen(filename, "w");
		return (file == NULL) ? NULL : panda_log_proxy_new_file(file);
	}
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
		DESTROY_MEMBER(proxy, printf);
		pthread_mutex_destroy(&proxy->mutex);
		free(proxy);
	}
}

void panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const char *message) {
	const char *name;
	if (proxy->printf == NULL)
		return;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif

	name = panda_assembler_get_name(assembler);
	proxy->printf(proxy->printf_data, "%s%s%s%s", name == NULL ? "" : name, name == NULL ? "" : "\t", panda_code_str(code), id == NULL ? "" : "\t");
	if (id != NULL) {
		panda_seqid_xprint(id, proxy->printf, proxy->printf_data);
	}
	if (message != NULL) {
		proxy->printf(proxy->printf_data, "\t%s\n", message);
	} else {
		proxy->printf(proxy->printf_data, "\n");
	}
	if (code == PANDA_CODE_ID_PARSE_FAILURE && message != NULL) {
		proxy->printf(proxy->printf_data, "* * * * * Something is wrong with this ID. If tags are absent, try passing the -B option.\n* * * * * Consult `pandaseq-checkid \"%s\"` to get an idea of the problem..\n", message);
	} else if (code == PANDA_CODE_PHRED_OFFSET) {
		proxy->printf(proxy->printf_data, "* * * * * Using the default PHRED+33 offset, but no sequences had quality data under PHRED+64.\n* * * * * This is probably not what you want. Consult the manual about the -6 option.\n");
	} else if (code == PANDA_CODE_READ_TOO_LONG) {
		proxy->printf(proxy->printf_data, "* * * * * The input reads are longer than this version of PANDAseq can handle. Currently %zd nucleotides.\n", PANDA_MAX_LEN);
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}

void panda_log_proxy_write_str(
	PandaLogProxy proxy,
	const char *str) {
	if (proxy->printf == NULL)
		return;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	proxy->printf(proxy->printf_data, "%s\n", str);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}

void panda_log_proxy_write_overlap(
	PandaLogProxy proxy,
	PandaAssembler assembler) {
	char buffer[BSIZE];
	size_t it = 0;
	int len;
	size_t max;
	const char *name;
	if (proxy->printf == NULL)
		return;

	name = panda_assembler_get_name(assembler);

	len = name == NULL ? snprintf(buffer, BSIZE, "%s\t", name) : 0;

	len += snprintf(buffer + len, BSIZE - len, "STAT\tOVERLAPS\t%ld", panda_assembler_get_overlap_count(assembler, it));
	max = panda_assembler_get_longest_overlap(assembler);
	for (it = 1; it <= max; it++) {
		len += snprintf(buffer + len, BSIZE - len, " %ld", panda_assembler_get_overlap_count(assembler, it));
	}
	if (len < BSIZE) {
		panda_log_proxy_write_str(proxy, buffer);
	}
}

void panda_log_proxy_stat_double(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	double value) {
	const char *asm_name;
	if (proxy->printf == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	asm_name = panda_assembler_get_name(assembler);
	if (asm_name == NULL) {
		proxy->printf(proxy->printf_data, "STAT\t%s\t%f\n", panda_assembler_get_name(assembler), name, value);
	} else {
		proxy->printf(proxy->printf_data, "%s\tSTAT\t%s\t%f\n", asm_name, name, value);
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}

void panda_log_proxy_stat_long(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	long value) {
	const char *asm_name;
	if (proxy->printf == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	asm_name = panda_assembler_get_name(assembler);
	if (asm_name == NULL) {
		proxy->printf(proxy->printf_data, "STAT\t%s\t%ld\n", panda_assembler_get_name(assembler), name, value);
	} else {
		proxy->printf(proxy->printf_data, "%s\tSTAT\t%s\t%ld\n", asm_name, name, value);
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}

void panda_log_proxy_stat_size_t(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	size_t value) {
	const char *asm_name;
	if (proxy->printf == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	asm_name = panda_assembler_get_name(assembler);
	if (asm_name == NULL) {
		proxy->printf(proxy->printf_data, "STAT\t%s\t%zd\n", panda_assembler_get_name(assembler), name, value);
	} else {
		proxy->printf(proxy->printf_data, "%s\tSTAT\t%s\t%zd\n", asm_name, name, value);
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}

void panda_log_proxy_stat_str(
	PandaLogProxy proxy,
	PandaAssembler assembler,
	const char *name,
	const char *value) {
	const char *asm_name;
	if (proxy->printf == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	asm_name = panda_assembler_get_name(assembler);
	if (asm_name == NULL) {
		proxy->printf(proxy->printf_data, "STAT\t%s\t%s\n", panda_assembler_get_name(assembler), name, value);
	} else {
		proxy->printf(proxy->printf_data, "%s\tSTAT\t%s\t%s\n", asm_name, name, value);
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
}
