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

struct fastq_data {
	PandaLineBuf forward;
	PandaLineBuf reverse;
	PandaLogProxy logger;
	unsigned char qualmin;
	panda_qual forward_seq[MAX_LEN];
	size_t forward_seq_length;
	panda_qual reverse_seq[MAX_LEN];
	size_t reverse_seq_length;
	PandaTagging policy;
	bool seen_under_64;
	bool non_empty;
};

#define LOG(flag, code) do { if(panda_debug_flags & flag) panda_log_proxy_write(data->logger, (code), NULL, id, NULL); } while(0)
#define LOGV(flag, code, fmt, ...) do { if(panda_debug_flags & flag) { snprintf(static_buffer(), BUFFER_SIZE, fmt, __VA_ARGS__); panda_log_proxy_write(data->logger, (code), NULL, id, static_buffer()); }} while(0)
#define TOINDEX(val) (((int)(val)) < data->qualmin ? 0 : ((((int)(val)) > data->qualmin + PHREDMAX ? PHREDMAX : (int)(val)) - data->qualmin))
static bool read_seq(
	panda_seq_identifier *id,
	panda_qual *buffer,
	size_t max_len,
	PandaLineBuf linebuf,
	char *table,
	struct fastq_data *data,
	size_t *length) {
	const char *input;
	size_t pos = 0;
	size_t qpos = 0;
	input = panda_linebuf_next(linebuf);
	if (input == NULL) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
		return false;
	}
	for (; pos < max_len; input++) {
		if (*input == '\0') {
			if (pos != 0) {
				break;
			}
		} else {
			if ((buffer[pos++].nt = table[*input & 0x1F]) == '\0') {
				LOGV(PANDA_DEBUG_FILE, PANDA_CODE_BAD_NT, "%c@%zd", *input, pos);
				return false;
			}
		}
	}
	input = panda_linebuf_next(linebuf);
	if (input == NULL) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
		return false;
	}
	if (*input != '+') {
		/* Check if we have more sequence... */
		if ((table[*input & 0x1F]) != '\0') {
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_READ_TOO_LONG);
		} else {
			/* Or just junk... */
			LOG(PANDA_DEBUG_FILE, PANDA_CODE_PARSE_FAILURE);
		}
		return false;
	}
	input = panda_linebuf_next(linebuf);
	if (input == NULL) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PREMATURE_EOF);
		return false;
	}
	for (; *input != '\0'; input++) {
		if (*input < 64) {
			data->seen_under_64 = true;
		}
		buffer[qpos++].qual = TOINDEX(*input);
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

static bool stream_next_seq(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	struct fastq_data *data) {
	panda_seq_identifier rid;
	PandaIdFmt format;
	const char *line;
	int fdir;
	int rdir;

	*forward = NULL;
	*reverse = NULL;
	*forward_length = 0;
	*reverse_length = 0;

	if ((line = panda_linebuf_next(data->forward)) == NULL) {
		return false;
	}
	if ((fdir = panda_seqid_parse_fail(id, line + 1, data->policy, &format, NULL)) == 0) {
		LOGV(PANDA_DEBUG_FILE, PANDA_CODE_ID_PARSE_FAILURE, "%s", line + 1);
		return false;
	}
	if ((line = panda_linebuf_next(data->reverse)) == NULL) {
		return false;
	}
	if ((rdir = panda_seqid_parse(&rid, line + 1, data->policy)) == 0) {
		LOGV(PANDA_DEBUG_FILE, PANDA_CODE_ID_PARSE_FAILURE, "%s", line + 1);
		return false;
	}
	if (!panda_seqid_equal(id, &rid) || (panda_idfmt_has_direction(format) && rdir == fdir)) {
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_NOT_PAIRED);
		return false;
	}
	if (format == PANDA_IDFMT_CASAVA_1_7) {
		/* We know that CASAVA 1.7+ is always PHRED+33, so supress the warning. */
		data->seen_under_64 = true;
	}
	if (!read_seq(id, data->forward_seq, MAX_LEN, data->forward, iupac_forward, data, forward_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		return false;
	}
	if (!read_seq(id, data->reverse_seq, MAX_LEN, data->reverse, iupac_reverse, data, reverse_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		return false;
	}
	*forward = data->forward_seq;
	*reverse = data->reverse_seq;
	return true;
}

static void stream_destroy(
	struct fastq_data *data) {
	if (data->non_empty && !data->seen_under_64 && data->qualmin < 64) {
		/* Used in the LOG macro. */
		panda_seq_identifier *id = NULL;
		LOG(PANDA_DEBUG_FILE, PANDA_CODE_PHRED_OFFSET);
	}
	panda_linebuf_free(data->forward);
	panda_linebuf_free(data->reverse);
	panda_log_proxy_unref(data->logger);
	free(data);
}

PandaNextSeq panda_create_fastq_reader(
	PandaBufferRead forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaBufferRead reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy) {
	struct fastq_data *data;
	data = malloc(sizeof(struct fastq_data));
	data->forward = panda_linebuf_new(forward, forward_data, forward_destroy);
	data->reverse = panda_linebuf_new(reverse, reverse_data, reverse_destroy);
	data->logger = panda_log_proxy_ref(logger);
	data->qualmin = qualmin;
	data->policy = policy;
	data->seen_under_64 = false;
	data->non_empty = false;
	*user_data = data;
	*destroy = (PandaDestroy) stream_destroy;
	return (PandaNextSeq) stream_next_seq;
}
