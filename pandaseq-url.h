/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2013  Andre Masella

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
#ifndef _PANDASEQ_URL_H
#        define _PANDASEQ_URL_H
#        include <pandaseq.h>
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
EXTERN_C_BEGIN
/**
 * Decompress BZipped data.
 * @source: (closure source_data) (scope notified) (allow-none): the underlying stream to decompress.
 * @verbosity: the BZip logging level.
 * Returns:(closure data) (scope notified) (allow-none): the read function.
 */
PandaBufferRead panda_bz_decompress(
	PandaBufferRead source,
	void *source_data,
	PandaDestroy source_destroy,
	int verbosity,
	void **user_data,
	PandaDestroy *destroy);

/**
 * Open a URL and read the sequence.
 * @url: the URL, as understood by cURL.
 * Returns:(closure data) (scope notified) (allow-none): the read function.
 */
PandaBufferRead panda_open_url(
	const char *url,
	PandaLogProxy logger,
	void **data,
	PandaDestroy *destroy);

/**
 * Increment the reference count on the cURL library.
 *
 * Since cURL needs to be initialised, PANDAseq will do this automatically when a URL is opened and automatically call the matching clean up when all readers have been disposed.
 *
 * If the program wishes to use cURL, it should call this method to increment the reference count on PANDAseq's internal counter, such that it will not clean up the cURL library while in use.
 *
 * Returns: whether the library was successfully initialised.
 */
bool panda_curl_ref(
	void);

/**
 * Decrement the reference count on the cURL library.
 */
void panda_curl_unref(
	void);

EXTERN_C_END
#endif
