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
#include <string.h>
#include "pandaseq.h"
#include "buffer.h"

const char *panda_idfmt_str(
	PandaIdFmt format) {
	switch (format) {
	case PANDA_IDFMT_EBI_SRA:
		return "EBI Short Read Archive";
	case PANDA_IDFMT_SRA:
		return "NCBI Short Read Archive";
	case PANDA_IDFMT_CASAVA_1_4:
		return "CASAVA 1.4-1.6";
	case PANDA_IDFMT_CASAVA_1_7:
		return "CASAVA 1.7+";
	default:
		return "unknown";
	}
}

bool panda_idfmt_has_direction(
	PandaIdFmt format) {
	return format != PANDA_IDFMT_EBI_SRA && format != PANDA_IDFMT_SRA;
}

void panda_seqid_clear(
	panda_seq_identifier *id) {
	id->instrument[0] = '\0';
	id->run[0] = '\0';
	id->flowcell[0] = '\0';
	id->lane = 0;
	id->tile = 0;
	id->x = 0;
	id->y = 0;
	id->tag[0] = '\0';
}

#define STRCMP(x) do { if (result == 0) result = strncmp(one->x, two->x, sizeof(((panda_seq_identifier*)0)->x)); } while (0)
int panda_seqid_compare(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two) {
	int result = 0;
	STRCMP(instrument);
	STRCMP(run);
	STRCMP(flowcell);
	if (result == 0)
		result = one->lane - two->lane;
	if (result == 0)
		result = one->tile - two->tile;
	if (result == 0)
		result = one->x - two->x;
	if (result == 0)
		result = one->y - two->y;
	STRCMP(tag);
	return result;
}

#undef STRCMP

void panda_seqid_copy(
	const panda_seq_identifier *src,
	panda_seq_identifier *dest) {
	dest->lane = src->lane;
	dest->tile = src->tile;
	dest->x = src->x;
	dest->y = src->y;
	strncpy(dest->instrument, src->instrument, sizeof(src->instrument));
	strncpy(dest->run, src->run, sizeof(src->run));
	strncpy(dest->flowcell, src->flowcell, sizeof(src->flowcell));
	strncpy(dest->tag, src->tag, PANDA_TAG_LEN);
}

#define STRCMP(x) (strncmp(one->x, two->x, sizeof(((panda_seq_identifier*)0)->x)) == 0)
bool panda_seqid_equal(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two) {
	return one->lane == two->lane && one->tile == two->tile && one->x == two->x && one->y == two->y && STRCMP(instrument) && STRCMP(run) && STRCMP(flowcell) && STRCMP(tag);
}

#undef STRCMP

void panda_seqid_print(
	const panda_seq_identifier *id,
	FILE *file) {
	panda_seqid_xprint(id, (PandaPrintf) fprintf, file);
}

void panda_seqid_xprint(
	const panda_seq_identifier *id,
	PandaPrintf xprintf,
	void *x) {
	if (id == NULL)
		return;
	xprintf(x, "%s:%s:%s:%d:%d:%d:%d:%s", id->instrument, id->run, id->flowcell, id->lane, id->tile, id->x, id->y, id->tag);
}

const char *panda_seqid_str(
	const panda_seq_identifier *id) {
	char *buffer = seqid_buffer();
	if (id == NULL)
		return NULL;
	panda_seqid_xprint(id, (PandaPrintf) bufferprintf, buffer);
	return buffer;
}

int panda_seqid_parse(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy) {
	const char *endptr;
	PandaIdFmt detected_format;
	return panda_seqid_parse_fail(id, input, policy, &detected_format, &endptr);
}

#define PARSE_CHUNK_MAYBE for(;**endptr != '\0' && **endptr != ':' && **endptr != '#' && **endptr != '/' && **endptr != ' '; (*endptr)++)
#define PARSE_CHUNK if (**endptr == '\0') return 0; PARSE_CHUNK_MAYBE
#define PARSE_INT do { value = 0; PARSE_CHUNK { if (**endptr >= '0' && **endptr <= '9') { value = 10*value + (int)(**endptr - '0'); } else { return 0; } } } while(0)
#define PARSE_SRA_INT do { value = 0; for(;**endptr != '\0' && **endptr != '.' && **endptr != ' '; (*endptr)++){ if (**endptr >= '0' && **endptr <= '9') { value = 10*value + (int)(**endptr - '0'); } else { return 0; } } } while(0)
#define PARSE_STR(target) do { dest = target; PARSE_CHUNK { if ((size_t) (dest - target) > sizeof(target)) return 0; *dest++ = (**endptr); } *dest = '\0'; } while(0);

int panda_seqid_parse_fail(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy,
	PandaIdFmt *detected_format,
	const char **endptr) {
	char *dest;
	const char *temp;
	int value;

	if (endptr == NULL)
		endptr = &temp;
	*endptr = input;

	if (strlen(input) > 3 && (input[0] == 'E' || input[0] == 'S') && input[1] == 'R' && input[2] == 'R') {
		/* Short Read Archive. Most of the data is too mangled to get out. */
		if (detected_format != NULL)
			*detected_format = (input[0] == 'S') ? PANDA_IDFMT_SRA : PANDA_IDFMT_EBI_SRA;
		*endptr += 3;
		panda_seqid_clear(id);
		PARSE_SRA_INT;
		(*endptr)++;
		sprintf(id->instrument, "%cRR%d", (int) input[0], value);
		PARSE_SRA_INT;
		(*endptr)++;
		id->lane = value;
		(*endptr)++;
		/* endptr still contains stuff, but it's just mangled SRA version of the
		 * Illumina header, so ignore it. */
		return 1;
	} else if (strchr(input, '/') != NULL) {
		/* Old CASAVA 1.4-1.6 format */
		if (detected_format != NULL)
			*detected_format = PANDA_IDFMT_CASAVA_1_4;
		id->run[0] = '\0';
		id->flowcell[0] = '\0';
		PARSE_STR(id->instrument);
		(*endptr)++;
		PARSE_INT;
		(*endptr)++;
		id->lane = value;
		PARSE_INT;
		(*endptr)++;
		id->tile = value;
		PARSE_INT;
		(*endptr)++;
		id->x = value;
		PARSE_INT;
		(*endptr)++;
		id->y = value;
		dest = id->tag;
		id->tag[0] = '\0';
		if (*(*endptr - 1) == '#') {
			PARSE_CHUNK_MAYBE {
				if (dest >= &id->tag[PANDA_TAG_LEN])
					return 0;
				*dest++ = (**endptr);
				*dest = '\0';
			}
			(*endptr)++;
		}
		if (policy != PANDA_TAG_OPTIONAL && policy != ((id->tag[0] == '\0') ? PANDA_TAG_ABSENT : PANDA_TAG_PRESENT)) {
			return 0;
		}
		PARSE_INT;
		(*endptr)++;
		return value;
	} else {
		/* New CASAVA 1.7+ format */
		int mate;
		if (detected_format != NULL)
			*detected_format = PANDA_IDFMT_CASAVA_1_7;
		PARSE_STR(id->instrument);
		(*endptr)++;
		PARSE_STR(id->run);
		(*endptr)++;
		PARSE_STR(id->flowcell);
		(*endptr)++;
		PARSE_INT;
		(*endptr)++;
		id->lane = value;
		PARSE_INT;
		(*endptr)++;
		id->tile = value;
		PARSE_INT;
		(*endptr)++;
		id->x = value;
		PARSE_INT;
		(*endptr)++;
		id->y = value;
		PARSE_INT;
		(*endptr)++;
		mate = value;
		PARSE_CHUNK;
		(*endptr)++;
		/* filtered */
		PARSE_INT;
		(*endptr)++;
		/* control bits */
		dest = id->tag;
		id->tag[0] = '\0';
		PARSE_CHUNK_MAYBE {
			if (dest >= &id->tag[PANDA_TAG_LEN])
				return 0;
			*dest++ = (**endptr);
			*dest = '\0';
		}
		if (policy != PANDA_TAG_OPTIONAL && policy != ((id->tag[0] == '\0') ? PANDA_TAG_ABSENT : PANDA_TAG_PRESENT)) {
			return 0;
		}
		return mate;
	}
}
