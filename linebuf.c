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
#include <string.h>
#include "pandaseq.h"
#include "misc.h"

struct panda_linebuf {
	char data[10 * MAX_LEN];
	size_t data_length;
	size_t offset;
	 MANAGED_MEMBER(
		PandaBufferRead,
		read);
};

PandaLineBuf panda_linebuf_new(
	PandaBufferRead read,
	void *read_data,
	PandaDestroy read_destroy) {
	PandaLineBuf buffer;
	if (read == NULL)
		return NULL;
	buffer = malloc(sizeof(struct panda_linebuf));
	buffer->data_length = 0;
	buffer->offset = 0;
	buffer->read = read;
	buffer->read_data = read_data;
	buffer->read_destroy = read_destroy;
	return buffer;
}

void panda_linebuf_free(
	PandaLineBuf linebuf) {
	if (linebuf == NULL)
		return;
	DESTROY_MEMBER(linebuf, read);
	free(linebuf);
}

const char *panda_linebuf_next(
	PandaLineBuf linebuf) {
	char *end;
	if (linebuf->offset > 0) {
		memmove(linebuf->data, linebuf->data + linebuf->offset, linebuf->data_length - linebuf->offset);
		linebuf->data_length -= linebuf->offset;
		linebuf->offset = 0;
	}

	while ((end = memchr(linebuf->data, '\n', linebuf->data_length)) == NULL && linebuf->data_length < sizeof(linebuf->data)) {
		size_t new_bytes = 0;
		if (!linebuf->read(linebuf->data + linebuf->data_length, sizeof(linebuf->data) - linebuf->data_length, &new_bytes, linebuf->read_data)) {
			return NULL;
		}
		if (new_bytes == 0) {
			end = linebuf->data + linebuf->data_length + 1;
			break;
		}
		linebuf->data_length += new_bytes;
	}
	if (end == NULL || linebuf->data_length == 0 || *end == '\0') {
		return NULL;
	}

	/* White out any carriage returns if we get DOS-formatted files. */
	if (end != linebuf->data && end[-1] == '\r') {
		end[-1] = '\0';
	}

	*end = '\0';
	linebuf->offset = end - linebuf->data + 1;
	return linebuf->data;
}
