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
	xprintf(x, "%s:%d:%s:%d:%d:%d:%d:%s", id->instrument, id->run, id->flowcell, id->lane, id->tile, id->x, id->y, id->tag);
}

const char *panda_seqid_str(
	const panda_seq_identifier *id) {
	char *buffer = static_buffer();
	if (id == NULL)
		return NULL;
	panda_seqid_xprint(id, (PandaPrintf) bufferprintf, buffer);
	return buffer;
}

bool panda_seqid_equal(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two) {
	return one->run == two->run && one->lane == two->lane && one->tile == two->tile && one->x == two->x && one->y == two->y && strcmp(one->instrument, two->instrument) == 0 && strcmp(one->flowcell, two->flowcell) == 0 && strcmp(one->tag, two->tag) == 0;
	;
}

int panda_seqid_compare(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two) {
	int result;
	result = strcmp(one->instrument, two->instrument);
	if (result == 0)
		result = one->run - two->run;
	if (result == 0)
		result = strcmp(one->flowcell, two->flowcell);
	if (result == 0)
		result = one->lane - two->lane;
	if (result == 0)
		result = one->tile - two->tile;
	if (result == 0)
		result = one->x - two->x;
	if (result == 0)
		result = one->y - two->y;
	if (result == 0)
		result = strcmp(one->tag, two->tag);
	return result;
}

#define PARSE_CHUNK_MAYBE for(;**endptr != '\0' && **endptr != ':' && **endptr != '#' && **endptr != '/' && **endptr != ' '; (*endptr)++)
#define PARSE_CHUNK if (**endptr == '\0') return 0; PARSE_CHUNK_MAYBE
#define PARSE_INT do { value = 0; PARSE_CHUNK { if (**endptr >= '0' && **endptr <= '9') { value = 10*value + (int)(**endptr - '0'); } else { return 0; } } } while(0)

void panda_seqid_clear(
	panda_seq_identifier *id) {
	id->instrument[0] = '\0';
	id->run = 0;
	id->flowcell[0] = '\0';
	id->lane = 0;
	id->tile = 0;
	id->x = 0;
	id->y = 0;
	id->tag[0] = '\0';
}

int panda_seqid_parse(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy) {
	const char *endptr;
	bool old;
	return panda_seqid_parse_fail(id, input, policy, &old, &endptr);
}

int panda_seqid_parse_fail(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy,
	bool *old,
	const char **endptr) {
	char *dest;
	bool has_tag;
	const char *temp;
	int value;

	if (endptr == NULL)
		endptr = &temp;

	*endptr = input;
	if (strchr(input, '/') != NULL) {
		/* Old CASAVA 1.4-1.6 format */
		if (old != NULL)
			*old = true;
		id->run = 0;
		id->flowcell[0] = '\0';
		dest = id->instrument;
		PARSE_CHUNK {
			if (dest - id->instrument > sizeof(id->instrument))
				return 0;
			*dest++ = (**endptr);
		}
		(*endptr)++;
		*dest = '\0';
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
		has_tag = id->tag[0] != '\0';
		if (!has_tag && policy == PANDA_TAG_PRESENT || has_tag && policy == PANDA_TAG_ABSENT) {
			return 0;
		}
		PARSE_INT;
		(*endptr)++;
		return value;
	} else {
		/* New CASAVA 1.7+ format */
		int mate;
		if (old != NULL)
			*old = false;
		dest = id->instrument;
		PARSE_CHUNK {
			if (dest - id->instrument > sizeof(id->instrument))
				return 0;
			*dest++ = (**endptr);
		}
		(*endptr)++;
		*dest = '\0';
		PARSE_INT;
		(*endptr)++;
		id->run = value;
		dest = id->flowcell;
		PARSE_CHUNK {
			if (dest - id->flowcell > sizeof(id->flowcell))
				return 0;
			*dest++ = (**endptr);
		}
		(*endptr)++;
		*dest = '\0';
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
		has_tag = id->tag[0] != '\0';
		if (!has_tag && policy == PANDA_TAG_PRESENT || has_tag && policy == PANDA_TAG_ABSENT) {
			return 0;
		}
		return mate;
	}
}
