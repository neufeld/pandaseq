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
#include <string.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

struct panda_writer {
	size_t refcnt;
	 MANAGED_MEMBER(
		PandaBufferWrite,
		write);
	PandaWriter commit_slave;
#ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
	pthread_key_t buffers;
	struct write_buffer *buffer_list;
#endif
};

#ifdef HAVE_PTHREAD
struct write_buffer {
	char uncommitted[2048];
	size_t uncommitted_length;
	char committed[20480];
	size_t committed_length;
	PandaWriter owner;
	struct write_buffer *next;
};

static struct write_buffer *get_write_buffer(
	PandaWriter writer) {
	struct write_buffer *data = pthread_getspecific(writer->buffers);
	if (data == NULL) {
		data = malloc(sizeof(struct write_buffer));
		data->uncommitted_length = 0;
		data->committed_length = 0;
		data->owner = writer;
		pthread_setspecific(writer->buffers, data);
		pthread_mutex_lock(&writer->mutex);
		data->next = writer->buffer_list;
		writer->buffer_list = data;
		pthread_mutex_unlock(&writer->mutex);
	}
	return data;
}

static void flush_buffer(
	PandaWriter writer,
	struct write_buffer *data) {
	pthread_mutex_lock(&writer->mutex);
	data->owner->write(data->committed, data->committed_length, data->owner->write_data);
	data->owner->write(data->uncommitted, data->uncommitted_length, data->owner->write_data);
	data->uncommitted_length = 0;
	data->committed_length = 0;
	pthread_mutex_unlock(&writer->mutex);
}
#endif

PandaWriter panda_writer_new(
	PandaBufferWrite write,
	void *write_data,
	PandaDestroy write_destroy) {
	PandaWriter writer = malloc(sizeof(struct panda_writer));
	writer->refcnt = 1;
	writer->commit_slave = NULL;
	writer->write = write;
	writer->write_data = write_data;
	writer->write_destroy = write_destroy;
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&writer->mutex, NULL);
	pthread_key_create(&writer->buffers, NULL);
	writer->buffer_list = NULL;
#endif
	return writer;
}

static void file_write(
	const char *buffer,
	size_t buffer_length,
	FILE *file) {
	fwrite(buffer, 1, buffer_length, file);
}

static void bzip_write(
	const char *buffer,
	size_t buffer_length,
	BZFILE *file) {
	int bzerror;
	BZ2_bzWrite(&bzerror, file, (void *) buffer, buffer_length);
}

PandaWriter panda_writer_new_stderr(
	void) {
	return panda_writer_new((PandaBufferWrite) file_write, stderr, NULL);
}

PandaWriter panda_writer_new_stdout(
	void) {
	return panda_writer_new((PandaBufferWrite) file_write, stdout, NULL);
}

PandaWriter panda_writer_new_file(
	FILE *file) {
	return panda_writer_new((PandaBufferWrite) file_write, file, (PandaDestroy) fclose);
}

static void null_write(
	const char *buffer,
	size_t buffer_length,
	void *data) {
}

PandaWriter panda_writer_new_null(
	) {
	return panda_writer_new(null_write, NULL, NULL);
}

PandaWriter panda_writer_open_file(
	const char *filename,
	bool bzip) {
	if (bzip) {
		BZFILE *file;
		file = BZ2_bzopen(filename, "w");
		return (file == NULL) ? NULL : panda_writer_new((PandaBufferWrite) bzip_write, file, (PandaDestroy) BZ2_bzclose);
	} else {
		FILE *file;
		file = fopen(filename, "w");
		return (file == NULL) ? NULL : panda_writer_new_file(file);
	}
}

PandaWriter panda_writer_ref(
	PandaWriter writer) {
#ifdef HAVE_PTHREAD
	flush_buffer(writer, get_write_buffer(writer));
	pthread_mutex_lock(&writer->mutex);
#endif
	writer->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&writer->mutex);
#endif
	return writer;
}

void panda_writer_unref(
	PandaWriter writer) {
	size_t count;
	if (writer == NULL)
		return;
#ifdef HAVE_PTHREAD
	flush_buffer(writer, get_write_buffer(writer));
	pthread_mutex_lock(&writer->mutex);
#endif
	count = --(writer->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&writer->mutex);
#endif
	if (count == 0) {
#ifdef HAVE_PTHREAD
		struct write_buffer *data;

		pthread_key_delete(writer->buffers);
		pthread_mutex_destroy(&writer->mutex);

		data = writer->buffer_list;
		while (data != NULL) {
			struct write_buffer *temp = data->next;
			writer->write(data->committed, data->committed_length, data->owner->write_data);
			writer->write(data->uncommitted, data->uncommitted_length, data->owner->write_data);
			free(data);
			data = temp;
		}
		if (writer->commit_slave != NULL)
			panda_writer_unref(writer->commit_slave);
#endif
		DESTROY_MEMBER(writer, write);
		free(writer);
	}
}

void panda_writer_append(
	PandaWriter writer,
	const char *format,
	...) {
	va_list va;
	va_start(va, format);
	panda_writer_append_v(writer, format, va);
	va_end(va);
}

void panda_writer_append_c(
	PandaWriter writer,
	char c) {
#ifdef HAVE_PTHREAD
	struct write_buffer *data = get_write_buffer(writer);
	if (data->uncommitted_length < sizeof(data->uncommitted)) {
		data->uncommitted[data->uncommitted_length] = c;
		data->uncommitted_length++;
	}
#else
	writer->write(&c, 1, writer->write_data);
#endif
}

void panda_writer_append_id(
	PandaWriter writer,
	const panda_seq_identifier *id) {
	panda_seqid_xprint(id, (PandaPrintf) panda_writer_append, writer);
}

void panda_writer_append_v(
	PandaWriter writer,
	const char *format,
	va_list va) {
#ifdef HAVE_PTHREAD
	struct write_buffer *data = get_write_buffer(writer);
	data->uncommitted_length += vsnprintf(data->uncommitted + data->uncommitted_length, sizeof(data->uncommitted) - data->uncommitted_length, format, va);
#else
	char buffer[2048];
	size_t buffer_length;
	buffer_length = vsnprintf(buffer, sizeof(buffer), format, va);
	writer->write(buffer, buffer_length, writer->write_data);
#endif
}

void panda_writer_commit(
	PandaWriter writer) {
#ifdef HAVE_PTHREAD
	struct write_buffer *data = get_write_buffer(writer);
	if (sizeof(data->committed) - data->committed_length < data->uncommitted_length) {
		flush_buffer(writer, data);
	} else {
		memcpy(data->committed + data->committed_length, data->uncommitted, data->uncommitted_length);
		data->committed_length += data->uncommitted_length;
		data->uncommitted_length = 0;
	}
	if (writer->commit_slave != NULL) {
		panda_writer_commit(writer->commit_slave);
	}
#endif
}

void panda_writer_flush(
	PandaWriter writer) {
#ifdef HAVE_PTHREAD
	flush_buffer(writer, get_write_buffer(writer));
#endif
}

PandaWriter panda_writer_get_slave(
	PandaWriter writer) {
	return writer->commit_slave;
}

void panda_writer_set_slave(
	PandaWriter writer,
	PandaWriter slave) {
	PandaWriter check;
	for (check = slave; check != NULL; check = check->commit_slave) {
		if (check->commit_slave == writer) {
			return;
		}
	}
	if (writer->commit_slave != NULL) {
		panda_writer_unref(writer->commit_slave);
	}
	writer->commit_slave = panda_writer_ref(slave);
}
