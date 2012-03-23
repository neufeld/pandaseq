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

#include<bzlib.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<zlib.h>
#include "pandaseq.h"
#include "assembler.h"
#include "parser.h"
#include "prob.h"
#include "table.h"

#define BUFSIZE 1024
#define PARSE_CHUNK if (*input == '\0') return 0; for(;*input != '\0' && *input != ':' && *input != '#' && *input != '/' && *input != ' '; input++)
#define PARSE_INT do { value = 0; PARSE_CHUNK { if (*input >= '0' && *input <= '9') { value = 10*value + (*input - '0'); } else { return 0; } } } while(0)

static char iupac_forward[32] =
    { /* @ */ 0, /*A*/ PANDA_NT_A, /*B*/ PANDA_NT_C | PANDA_NT_G | PANDA_NT_T, /*C*/ PANDA_NT_C,
	 /*D*/ PANDA_NT_A | PANDA_NT_G | PANDA_NT_T, /*E*/ 0, /*F*/ 0, /*G*/ PANDA_NT_G,
	 /*H*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_T,
	 /*I*/ 0, /*J*/ 0, /*K*/ PANDA_NT_G | PANDA_NT_T, /*L*/ 0, /*M*/ PANDA_NT_A | PANDA_NT_C,
	 /*N*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T, /*O*/ 0, /*P*/ 0, /*Q*/ 0,
	 /*R*/ PANDA_NT_A | PANDA_NT_G,
	 /*S*/ PANDA_NT_C | PANDA_NT_G, /*T*/ PANDA_NT_T, /*U*/ PANDA_NT_T, /*V*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G,
	 /*W*/ PANDA_NT_A | PANDA_NT_T,
	 /*X*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T, /*Y*/ PANDA_NT_C | PANDA_NT_T, /*Z*/ 0, /*[ */ 0, /*\ */ 0,	/*] */
	0, /*^ */ 0, /*_*/ 0
};

static char iupac_reverse[32] =
    { /*@ */ 0, /*A*/ PANDA_NT_T, /*B*/ PANDA_NT_G | PANDA_NT_C | PANDA_NT_A, /*C*/ PANDA_NT_G,
	 /*D*/ PANDA_NT_T | PANDA_NT_C | PANDA_NT_A, /*E*/ 0, /*F*/ 0, /*G*/ PANDA_NT_C,
	 /*H*/ PANDA_NT_T | PANDA_NT_G | PANDA_NT_A,
	 /*I*/ 0, /*J*/ 0, /*K*/ PANDA_NT_C | PANDA_NT_A, /*L*/ 0, /*M*/ PANDA_NT_T | PANDA_NT_G,
	 /*N*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T, /*O*/ 0, /*P*/ 0, /*Q*/ 0,
	 /*R*/ PANDA_NT_T | PANDA_NT_C,
	 /*S*/ PANDA_NT_G | PANDA_NT_C, /*T*/ PANDA_NT_A, /*U*/ PANDA_NT_A, /*V*/ PANDA_NT_T | PANDA_NT_G | PANDA_NT_C,
	 /*W*/ PANDA_NT_T | PANDA_NT_A,
	 /*X*/ PANDA_NT_A | PANDA_NT_C | PANDA_NT_G | PANDA_NT_T, /*Y*/ PANDA_NT_G | PANDA_NT_A, /*Z*/ 0, /*[ */ 0, /*\ */ 0,	/*] */
	0, /*^ */ 0, /*_*/ 0
};
double panda_quality_probability(const panda_qual *q) {
	return exp(panda_quality_log_probability(q));
}
double panda_quality_log_probability(const panda_qual *q) {
	int index = q->qual;
	if (index < 0) {
		index = 0;
	} else if (index > PHREDMAX) {
		index = PHREDMAX;
	}
	return qual_score[index];
}

PandaAssembler panda_assembler_open_gz(char *forward, char *reverse, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy, unsigned char qualmin)
{
	gzFile *forward_file;
	gzFile *reverse_file;
	forward_file = gzopen(forward, "r");
	if (forward_file == NULL) {
		logger(logger_data, PANDA_CODE_NO_FILE, forward);
		return NULL;
	}
	reverse_file = gzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(logger_data, PANDA_CODE_NO_FILE, reverse);
		gzclose(forward_file);
		return NULL;
	}

	return panda_assembler_new_fastq_reader(gzgetc, forward_file, (PandaDestroy) gzclose, gzgetc, reverse_file, (PandaDestroy) gzclose, logger, logger_data, logger_destroy, qualmin);
}

static int bzgetc(BZFILE *file) {
	char c;
	if (BZ2_bzread(file, &c, 1) != 1) {
		return EOF;
	}
	return c;
}

PandaAssembler panda_assembler_open_bz2(char *forward, char *reverse, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy, unsigned char qualmin)
{
	BZFILE *forward_file;
	BZFILE *reverse_file;
	forward_file = BZ2_bzopen(forward, "r");
	if (forward_file == NULL) {
		logger(logger_data, PANDA_CODE_NO_FILE, forward);
		return NULL;
	}
	reverse_file = BZ2_bzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(logger_data, PANDA_CODE_NO_FILE, reverse);
		BZ2_bzclose(forward_file);
		return NULL;
	}
	return panda_assembler_new_fastq_reader(bzgetc, forward_file, BZ2_bzclose, bzgetc, reverse_file, BZ2_bzclose, logger, logger_data, logger_destroy, qualmin);
}

struct stream_data {
	MANAGED_MEMBER(PandaNextChar, next);
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
};

int panda_seqid_parse(panda_seq_identifier *id, char *input)
{
	char *dest;
	int value;
	if (strchr(input, '#') != NULL) {
		/* Old CASAVA 1.4-1.6 format */
		id->run = 0;
		id->flowcell[0] = '\0';
		dest = id->instrument;
		PARSE_CHUNK {
			*dest++ = (*input);
		}
		input++;
		*dest = '\0';
		PARSE_INT;
		input++;
		id->lane = value;
		PARSE_INT;
		input++;
		id->tile = value;
		PARSE_INT;
		input++;
		id->x = value;
		PARSE_INT;
		input++;
		id->y = value;
		dest = id->tag;
		*dest = '\0';
		PARSE_CHUNK {
			if (dest >= &id->tag[PANDA_TAG_LEN])
				return 0;
			*dest++ = (*input);
			*dest = '\0';
		}
		if (id->tag[0] == '\0') {
			return 0;
		}
		input++;
		PARSE_INT;
		input++;
		return value;
	} else {
		/* New CASAVA 1.7+ format */
		int mate;
		dest = id->instrument;
		PARSE_CHUNK {
			*dest++ = (*input);
		}
		input++;
		*dest = '\0';
		PARSE_INT;
		input++;
		id->run = value;
		dest = id->flowcell;
		PARSE_CHUNK {
			*dest++ = (*input);
		}
		input++;
		*dest = '\0';
		PARSE_INT;
		input++;
		id->lane = value;
		PARSE_INT;
		input++;
		id->tile = value;
		PARSE_INT;
		input++;
		id->x = value;
		PARSE_INT;
		input++;
		id->y = value;
		PARSE_INT;
		input++;
		mate = value;
		PARSE_CHUNK;
		input++;
		/* filtered */
		PARSE_INT;
		input++;
		/* control bits */
		dest = id->tag;
		*dest = '\0';
		PARSE_CHUNK {
			if (dest >= &id->tag[PANDA_TAG_LEN])
				return 0;
			*dest++ = (*input);
			*dest = '\0';
		}
		if (id->tag[0] == '\0') {
			return 0;
		}
		return mate;
	}
}

static bool read_line(char* buffer, size_t max_len, struct stream_data *stream, size_t *len) {
	int pos = 0;
	if (stream->next == NULL)
		return false;
	while(pos < max_len) {
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

#define TOINDEX(val) (((int)(val)) < data->qualmin ? 0 : ((((int)(val)) > data->qualmin + PHREDMAX ? PHREDMAX : (int)(val)) - data->qualmin))
static bool read_seq(panda_seq_identifier *id, panda_qual* buffer, size_t max_len, struct stream_data *stream, char *table, struct fastq_data *data, size_t *length) {
	int pos = 0;
	int qpos = 0;
	int v = EOF;
	if (stream->next == NULL) 
		return false;
	while(pos < max_len) {
		v = stream->next(stream->next_data);
		if (v == EOF) {
			DESTROY_MEMBER(stream, next);
			LOG(data, PANDA_CODE_PREMATURE_EOF, NULL);
			return false;
		}
		if (v == '\r' || v == '\n') {
			if (pos != 0) {
				break;
			}
		} else {
			if ((buffer[pos++].nt = table[v & 0x1F]) == '\0') {
				LOG(data, PANDA_CODE_BAD_NT, id, v, pos);
				return false;
			}
		}
	}
	for(; v == '\n' || v == '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(data, PANDA_CODE_PREMATURE_EOF, id);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	if (v != '+') {
			LOG(data, PANDA_CODE_PARSE_FAILURE, id);
			DESTROY_MEMBER(stream, next);
			return false;
	}
	for(; v != '\n' && v != '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(data, PANDA_CODE_PREMATURE_EOF, id);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	for(; v == '\n' || v == '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(data, PANDA_CODE_PREMATURE_EOF, id);
			DESTROY_MEMBER(stream, next);
			return false;
		}
	}
	for(; v != '\n' && v != '\r'; v = stream->next(stream->next_data)) {
		if (v == EOF) {
			LOG(data, PANDA_CODE_PARSE_FAILURE, id);
			DESTROY_MEMBER(stream, next);
			return false;
		}
		buffer[qpos++].qual = TOINDEX(v);
	}
	
	if (pos == 0) {
		LOG(data, PANDA_CODE_NO_DATA, id);
		return false;
	}
	if (qpos != pos) {
		LOG(data, PANDA_CODE_NO_QUALITY_INFO, id);
		return false;
	}
	*length = pos;
	return true;
}
#undef TOINDEX

static bool stream_next_seq(panda_seq_identifier *id, panda_qual **forward, size_t *forward_length, panda_qual **reverse, size_t *reverse_length, struct fastq_data *data) {
	panda_seq_identifier rid;
	char buffer[BUFSIZE];
	size_t fsize;
	size_t rsize;
	int fdir;
	int rdir;

	*forward = NULL;
	*reverse = NULL;

	if (!read_line(buffer, BUFSIZE, &data->forward, &fsize)) {
		return false;
	}
	if ((fdir = panda_seqid_parse(id, &buffer[1])) == 0) {
		LOG(data, PANDA_CODE_ID_PARSE_FAILURE, buffer);
		return false;
	}	
	if (!read_line(buffer, BUFSIZE, &data->reverse, &rsize)) {
		return false;
	}
	if ((rdir = panda_seqid_parse(&rid, &buffer[1])) == 0) {
		LOG(data, PANDA_CODE_ID_PARSE_FAILURE, buffer);
		return false;
	}
	if (!panda_seqid_equal(id, &rid)) {
		LOG(data, PANDA_CODE_NOT_PAIRED, id, &rid);
		return false;
	}
	if (!read_seq(id, data->forward_seq, PANDA_MAX_LEN, &data->forward, iupac_forward, data, forward_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		LOG(data, PANDA_CODE_PARSE_FAILURE, id);
		return false;
	}
	if (!read_seq(id, data->reverse_seq, PANDA_MAX_LEN, &data->reverse, iupac_reverse, data, reverse_length)) {
		*forward_length = 0;
		*reverse_length = 0;
		LOG(data, PANDA_CODE_PARSE_FAILURE, id);
		return false;
	}
	*forward = data->forward_seq;
	*reverse = data->reverse_seq;
	return true;
}

static void stream_destroy(struct fastq_data *data) {
	DESTROY_MEMBER(&data->forward, next);
	DESTROY_MEMBER(&data->reverse, next);
	free(data);
}

PandaNextSeq panda_assembler_create_fastq_reader(PandaNextChar forward, void *forward_data, PandaDestroy forward_destroy, PandaNextChar reverse, void *reverse_data, PandaDestroy reverse_destroy, PandaLogger logger, void *logger_data, unsigned char qualmin, void **user_data, PandaDestroy *destroy) {
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
	*user_data = data;
	*destroy = (PandaDestroy) stream_destroy;
	return (PandaNextSeq) stream_next_seq;
}
PandaAssembler panda_assembler_new_fastq_reader(PandaNextChar forward, void *forward_data, PandaDestroy forward_destroy, PandaNextChar reverse, void *reverse_data, PandaDestroy reverse_destroy, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy, unsigned char qualmin) {
	void* user_data;
	PandaDestroy destroy;
	PandaNextSeq next;
	next = panda_assembler_create_fastq_reader(forward, forward_data, forward_destroy, reverse, reverse_data, reverse_destroy, logger, logger_data, qualmin, &user_data, &destroy);
	return panda_assembler_new(next, user_data, destroy, logger, logger_data, logger_destroy);
}

void panda_seqid_print(const panda_seq_identifier *id, FILE *file) {
	panda_seqid_xprint(id, (PandaPrintf) fprintf, file);
}

void panda_seqid_xprint(const panda_seq_identifier *id, PandaPrintf xprintf, void *x)
{
	if (id == NULL)
		return;
	xprintf(x, "%s:%d:%s:%d:%d:%d:%d:%s", id->instrument,
	       id->run, id->flowcell, id->lane, id->tile, id->x, id->y,
	       id->tag);
}

const char *panda_seqid_str(const panda_seq_identifier *id)
{
	static char buffer[1024];
	if (id == NULL)
		return NULL;
	panda_seqid_xprint(id, (PandaPrintf)sprintf, buffer);
	return buffer;
}

bool panda_seqid_equal(const panda_seq_identifier *one, const panda_seq_identifier *two)
{
	return one->run == two->run && one->lane == two->lane
	    && one->tile == two->tile && one->x == two->x && one->y == two->y
	    && strcmp(one->instrument, two->instrument) == 0
	    && strcmp(one->flowcell, two->flowcell) == 0
	    && strcmp(one->tag, two->tag) == 0;
	;
}

static char ntchar[16] =
    { 'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B',
	'N'
};

panda_nt panda_nt_from_ascii(char c) {
	return iupac_forward[c & 0x1F];
}
	
panda_nt panda_nt_from_ascii_complement(char c) {
	return iupac_reverse[c & 0x1F];
}
	
char panda_nt_to_ascii(panda_nt val) {
	if (val < 0 || val > 15) {
		return 'N';
	}
	return ntchar[(int)(val)];
}

void panda_assembler_set_forward_primer(PandaAssembler assembler, panda_nt *sequence, size_t length) {
	size_t it;
	if (length < PANDA_MAX_LEN) {
		for(it = 0; it < length; it++) {
			assembler->forward_primer[it] = sequence[it];
		}
		assembler->forward_primer_length = length;
		assembler->forward_trim = 0;
	}
}

void panda_assembler_set_reverse_primer(PandaAssembler assembler, panda_nt *sequence, size_t length) {
	size_t it;
	if (length < PANDA_MAX_LEN) {
		for(it = 0; it < length; it++) {
			assembler->reverse_primer[it] = sequence[it];
		}
		assembler->reverse_primer_length = length;
		assembler->reverse_trim = 0;
	}
}

panda_nt *panda_assembler_get_reverse_primer(PandaAssembler assembler, size_t *length) {
	*length = assembler->reverse_primer_length;
	return *length == 0 ? NULL : assembler->reverse_primer;
}

panda_nt *panda_assembler_get_forward_primer(PandaAssembler assembler, size_t *length) {
	*length = assembler->forward_primer_length;
	return *length == 0 ? NULL : assembler->forward_primer;
}

size_t panda_assembler_get_forward_trim(PandaAssembler assembler) {
	return assembler->forward_trim;
}

void panda_assembler_set_forward_trim(PandaAssembler assembler, size_t trim) {
	assembler->forward_trim = trim;
	assembler->forward_primer_length = 0;
}

size_t panda_assembler_get_reverse_trim(PandaAssembler assembler) {
	return assembler->reverse_trim;
}

void panda_assembler_set_reverse_trim(PandaAssembler assembler, size_t trim) {
	assembler->reverse_trim = trim;
	assembler->reverse_primer_length = 0;
}
