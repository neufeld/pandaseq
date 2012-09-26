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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "buffer.h"

PandaDebug panda_debug_flags = PANDA_DEBUG_DEFAULT;

#if HAVE_PTHREAD

#        define BUFFER(name, type, size) pthread_key_t PANDACONCAT(name, _key);
#        include "buffer.list"
#        undef BUFFER

__attribute__ ((constructor))
static void
lib_init(
	void) {
#        define BUFFER(name, type, size) pthread_key_create(&PANDACONCAT(name, _key), &free);
#        include "buffer.list"
#        undef BUFFER
}

__attribute__ (())
void
lib_destroy(
	void) {
#        define BUFFER(name, type, size) pthread_key_delete(PANDACONCAT(name, _key));
#        include "buffer.list"
#        undef BUFFER
}

static void *
get_buffer(
	pthread_key_t key,
	size_t size) {
	void *buffer;
	if ((buffer = pthread_getspecific(key)) == NULL) {
		buffer = malloc(size);
		pthread_setspecific(key, buffer);
	}
	return buffer;
}

#        define BUFFER(name, type, size) type *PANDACONCAT(name, _buffer)(void) { return get_buffer(PANDACONCAT(name, _key), sizeof(type) * size); }
#        include "buffer.list"
#        undef BUFFER
#else
#        define BUFFER(name, type, size) static type PANDACONCAT(name, buffer)[size]; type *PANDACONCAT(name, buffer)(void) { return PANDACONCAT(name, buffer); }
#        include "buffer.list"
#        undef BUFFER
#endif

void
bufferprintf(
	char *buffer,
	char *fmt,
	...) {
	va_list va;
	va_start(va, fmt);
	vsnprintf(buffer, BUFFER_SIZE, fmt, va);
	va_end(va);
}
