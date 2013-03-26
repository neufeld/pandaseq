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
#include <stdlib.h>
#include "pandaseq.h"
#include "buffer.h"
#include "misc.h"
#include "nt.h"
#include "prob.h"

struct stream_data {
	MANAGED_MEMBER(
		PandaNextChar,
		next);
};

struct fastq_data {
	struct stream_data forward;
	struct stream_data reverse;
	PandaLogger logger;
	void *logger_data;
	unsigned char qualmin;
	panda_qual forward_seq[PANDA_MAX_LEN];
	size_t forward_seq_length;
	panda_qual reverse_seq[PANDA_MAX_LEN];
	size_t reverse_seq_length;
	PandaTagging policy;
	bool seen_under_64;
	bool non_empty;
};

static bool
read_line(
	char *buffer,
	size_t max_len,
	struct stream_data *stream,
	size_t *len) {
	int pos = 0;
	if (stream->next == NULL)
		return false;
	while (pos < max_len) {
		int v = stream->next(stream->next_data);
		if (v == EOF) {
			DESTROY_MEMBER(stream, next);
			buffer[pos] = '\0';
			*len = pos;
			return pos > 0;
		}
		if (v == '\r' || v == '\n') {
			if (pos != 0) {
				buffer[pos] = '\0';
				return true;
			}
		} else {
			buffer[pos++] = v;
		}
	}
	return false;
}

#define LOG(flag, code) do { if(panda_debug_flags & flag) data->logger((code), id, NULL, data->logger_data); } while(0)
#define LOGV(flag, code, fmt, ...) do { if(panda_debug_flags & flag) { snprintf(static_buffer(), BUFFER_SIZE, fmt, __VA_ARGS__); data->logger((code), id, static_buffer(), data->logger_data); }} while(0)
#define TOINDEX(val) (((int)(val)) < data->qualmin ? 0 : ((((int)(val)) > data->qualmin + PHREDMAX ? PHREDMAX : (int)(val)) - data->qualmin))
static bool
read_seq(
	panda_seq_identifier *id,
	panda_qual *buffer,
	size_t max_len,
	struct stream_data *stream,
	char *table,
	struct fastq_data *data,
	size_t *length) {
	int pos = 0;
	int qpos = 0;
	int v = EOF;
	if (stream->next == NULL)
		return false;
	while (pos < max_len) {
		v = stream->next(stream->next_data);
		if (v == EOF) {
			DESTROY_MEMBER(stream, next);
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
			return false;
		}
		if (v == '\r' || v == '\n') {
			if (pos != 0) {
				break;
			}
		} else {
			if ((buffer[pos++].nt = table[v & 0x1F]) == '\0') {
				LOGV(PANDA_DEBUG_FILE, PANDA_CODE_BAD_NT, "%c@%d", v, pos);
				return false;
			}
		}
	}
	for (; v == '\n' || v == '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	if (v != '+') {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PARSE_FAILURE);
		DESTROY_MEMBER(stream, next);
		return false;
	}
	for (; v != '\n' && v != '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	for (; v == '\n' || v == '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	for (; v != '\n' && v != '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PARSE_FAILURE);
			DESTROY_MEMBER(stream, next);
			return false;
		}
		if (v < 64) {
			data->seen_under_64 = true;
		}
		buffer[qpos++].qual = TOINDEX(v);
	}

	if (pos == 0) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_NO_DATA);
		return false;
	}
	if (qpos != pos) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_NO_QUALITY_INFO);
		return false;
	}
	*length = pos;
	data->non_empty = true;
	return true;
}

#undef TOINDEX

static bool
stream_next_seq(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	struct fastq_data *data) {
	panda_seq_identifier rid;
	char buffer[BUFFER_SIZE];
	size_t fsize;
	size_t rsize;
	int fdir;
	int rdir;

	*forward = NULL;
	*reverse = NULL;

	if (!read_line(buffer, BUFFER_SIZE, &data->forward, &fsize)) {
		return false;
	}
	if ((fdir = panda_seqid_parse(id, &buffer[1], data->policy)) == 0) {
		LOGV(PANDA_DEBUG_FILE, PANDA_CODE_ID_PARSE_FAILURE, "%s", buffer);
		return false;
	}
	if (!read_line(buffer, BUFFER_SIZE, &data->reverse, &rsize)) {
		return false;
	}
	if ((rdir = panda_seqid_parse(&rid, &buffer[1], data->policy)) == 0) {
		LOGV(PANDA_DEBUG_FILE, PANDA_CODE_ID_PARSE_FAILURE, "%s", buffer);
		return false;
	}
	if (!panda_seqid_equal(id, &rid)) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_NOT_PAIRED);
		return false;
	}
	if (!read_seq(id, data->forward_seq, PANDA_MAX_LEN, &data->forward, iupac_forward, data, forward_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		return false;
	}
	if (!read_seq(id, data->reverse_seq, PANDA_MAX_LEN, &data->reverse, iupac_reverse, data, reverse_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		return false;
	}
	*forward = data->forward_seq;
	*reverse = data->reverse_seq;
	return true;
}

static void
stream_destroy(
	struct fastq_data *data) {
	if (data->non_empty && !data->seen_under_64 && data->qualmin < 64) {
		/* Used in the LOG macro. */
		panda_seq_identifier *id = NULL;
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PHRED_OFFSET);
	}
	DESTROY_MEMBER(&data->forward, next);
	DESTROY_MEMBER(&data->reverse, next);
	free(data);
}

PandaNextSeq
panda_create_fastq_reader(
	PandaNextChar forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaNextChar reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogger logger,
	void *logger_data,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy) {
	struct fastq_data *data;
	data = malloc(sizeof(struct fastq_data));
	data->forward.next = forward;
	data->forward.next_data = forward_data;
	data->forward.next_destroy = forward_destroy;
	data->reverse.next = reverse;
	data->reverse.next_data = reverse_data;
	data->reverse.next_destroy = reverse_destroy;
	data->logger = logger;
	data->logger_data = logger_data;
	data->qualmin = qualmin;
	data->policy = policy;
	data->seen_under_64 = false;
	data->non_empty = false;
	*user_data = data;
	*destroy = (PandaDestroy) stream_destroy;
	return (PandaNextSeq) stream_next_seq;
}
