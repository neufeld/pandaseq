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

void panda_seqid_print(const panda_seq_identifier *id, FILE *file)
{
	panda_seqid_xprint(id, (PandaPrintf) fprintf, file);
}

void panda_seqid_xprint(const panda_seq_identifier *id, PandaPrintf xprintf,
			void *x)
{
	if (id == NULL)
		return;
	xprintf(x, "%s:%d:%s:%d:%d:%d:%d:%s", id->instrument,
		id->run, id->flowcell, id->lane, id->tile, id->x, id->y,
		id->tag);
}

const char *panda_seqid_str(const panda_seq_identifier *id)
{
	char *buffer = static_buffer();
	if (id == NULL)
		return NULL;
	panda_seqid_xprint(id, (PandaPrintf) bufferprintf, buffer);
	return buffer;
}

bool panda_seqid_equal(const panda_seq_identifier *one,
		       const panda_seq_identifier *two)
{
	return one->run == two->run && one->lane == two->lane
	    && one->tile == two->tile && one->x == two->x && one->y == two->y
	    && strcmp(one->instrument, two->instrument) == 0
	    && strcmp(one->flowcell, two->flowcell) == 0
	    && strcmp(one->tag, two->tag) == 0;
	;
}

#define PARSE_CHUNK if (*input == '\0') return 0; for(;*input != '\0' && *input != ':' && *input != '#' && *input != '/' && *input != ' '; input++)
#define PARSE_INT do { value = 0; PARSE_CHUNK { if (*input >= '0' && *input <= '9') { value = 10*value + (*input - '0'); } else { return 0; } } } while(0)

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
