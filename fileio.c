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
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <zlib.h>
#if HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"
#ifdef HAVE_PTHREAD
#        include"pandaseq-mux.h"
#endif

static bool buff_read_gz(
	char *buf,
	size_t buf_len,
	size_t *read,
	void *data) {
	gzFile file = (gzFile) data;
	int code;
	code = gzread(file, buf, buf_len);
	if (code < 1) {
		*read = 0;
		return gzeof(file);
	}
	*read = code;
	return true;
}

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

PandaBufferRead panda_open_buffer(
	const char *file_name,
	PandaLogProxy logger,
	void **user_data,
	PandaDestroy *destroy) {
	char buffer[2];
	int fd;
	*user_data = NULL;
	*destroy = NULL;

	fd = open(file_name, O_RDONLY);
	if (fd < 0 || read(fd, &buffer, 2) != 2 || lseek(fd, 0, SEEK_SET) != 0) {
		panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, file_name);
		return NULL;
	}
	if (buffer[0] == 'B' && buffer[1] == 'Z') {
		BZFILE *bz_file;
		bz_file = BZ2_bzdopen(fd, "r");
		if (bz_file == NULL) {
			panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, file_name);
			close(fd);
			return NULL;
		}
		*user_data = bz_file;
		*destroy = BZ2_bzclose;
		return buff_read_bz2;
	} else {
		gzFile gz_file;
		gz_file = gzdopen(fd, "r");
		if (gz_file == NULL) {
			panda_log_proxy_write(logger, PANDA_CODE_NO_FILE, NULL, NULL, file_name);
			close(fd);
			return NULL;
		}
		*user_data = gz_file;
		*destroy = (PandaDestroy) gzclose;
		return buff_read_gz;
	}
}

PandaNextSeq panda_open_fastq(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	const char *index,
	void **user_data,
	PandaDestroy *destroy) {
	MANAGED_STACK(PandaBufferRead,
		forward_file);
	MANAGED_STACK(PandaBufferRead,
		reverse_file);
	MANAGED_STACK(PandaBufferRead,
		index_file);

	*user_data = NULL;
	*destroy = NULL;

	forward_file = panda_open_buffer(forward, logger, &forward_file_data, &forward_file_destroy);
	if (forward_file == NULL) {
		return NULL;
	}

	reverse_file = panda_open_buffer(reverse, logger, &reverse_file_data, &reverse_file_destroy);
	if (reverse_file == NULL) {
		DESTROY_STACK(forward_file);
		return NULL;
	}
	index_file = index == NULL ? NULL : panda_open_buffer(index, logger, &index_file_data, &index_file_destroy);
	if (index != NULL && index_file == NULL) {
		DESTROY_STACK(forward_file);
		DESTROY_STACK(reverse_file);
		return NULL;
	}

	return panda_create_fastq_reader(forward_file, forward_file_data, forward_file_destroy, reverse_file, reverse_file_data, reverse_file_destroy, logger, qualmin, policy, index_file, index_file_data, index_file_destroy, user_data, destroy);
}

PandaAssembler panda_assembler_open_fastq(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_fastq(forward, reverse, logger, qualmin, policy, NULL, &next_data, &next_destroy)) == NULL) {
		return NULL;
	}

	return panda_assembler_new(next, next_data, next_destroy, logger);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_fastq(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	PandaNextSeq next;
	void *next_data;
	PandaDestroy next_destroy;
	if ((next = panda_open_fastq(forward, reverse, logger, qualmin, policy, NULL, &next_data, &next_destroy)) == NULL) {

		return NULL;
	}
	return panda_mux_new(next, next_data, next_destroy, logger);
}
#endif
