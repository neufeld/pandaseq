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
#include <stdlib.h>
#include <zlib.h>
#if HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

static bool buff_read_gz(
	char *buf,
	size_t buf_len,
	size_t *read,
	void *data) {
	gzFile *file = (gzFile *) data;
	int code;
	code = gzread(file, buf, buf_len);
	if (code < 1) {
		*read = 0;
		return gzeof(file);
	}
	*read = code;
	return true;
}

PandaNextSeq panda_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy) {
	gzFile *forward_file;
	gzFile *reverse_file;
	*user_data = NULL;
	*destroy = NULL;

	forward_file = gzopen(forward, "r");
	if (forward_file == NULL) {
		panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, forward);
		return NULL;
	}
	reverse_file = gzopen(reverse, "r");
	if (reverse_file == NULL) {
		panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, reverse);
		gzclose(forward_file);
		return NULL;
	}
	return panda_create_fastq_reader(buff_read_gz, forward_file, (PandaDestroy) gzclose, buff_read_gz, reverse_file, (PandaDestroy) gzclose, logger, qualmin, policy, user_data, destroy);
}

PandaAssembler panda_assembler_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_gz(forward, reverse, logger, qualmin, policy, &next_data, &next_destroy)) == NULL) {
		return NULL;
	}
	return panda_assembler_new(next, next_data, next_destroy, logger);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_gz(forward, reverse, logger, qualmin, policy, &next_data, &next_destroy)) == NULL) {
		return NULL;
	}

	return panda_mux_new(next, next_data, next_destroy, logger);
}
#endif

static bool buff_read_bz2(
	char *buf,
	size_t buf_len,
	size_t *read,
	void *data) {
	BZFILE *file = (BZFILE *) data;
	int bzerror;
	*read = BZ2_bzRead(&bzerror, file, buf, buf_len);
	return bzerror == BZ_OK || bzerror == BZ_STREAM_END;
}

PandaNextSeq panda_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy) {
	BZFILE *forward_file;
	BZFILE *reverse_file;
	forward_file = BZ2_bzopen(forward, "r");
	if (forward_file == NULL) {
		panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, forward);
		return NULL;
	}
	reverse_file = BZ2_bzopen(reverse, "r");
	if (reverse_file == NULL) {
		panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, reverse);
		BZ2_bzclose(forward_file);
		return NULL;
	}
	return panda_create_fastq_reader(buff_read_bz2, forward_file, BZ2_bzclose, buff_read_bz2, reverse_file, BZ2_bzclose, logger, qualmin, policy, user_data, destroy);
}

PandaAssembler panda_assembler_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_bz2(forward, reverse, logger, qualmin, policy, &next_data, &next_destroy)) == NULL) {
		return NULL;
	}

	return panda_assembler_new(next, next_data, next_destroy, logger);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_bz2(forward, reverse, logger, qualmin, policy, &next_data, &next_destroy)) == NULL) {

		return NULL;
	}
	return panda_mux_new(next, next_data, next_destroy, logger);
}
#endif
