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

#ifndef _PANDASEQ_MUX_H
#        define _PANDASEQ_MUX_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/* === Constructors === */
/**
 * Create a new multiplexed data source from a sequence callback.
 *
 * The interface will guarantee that only one call will be made at a time to the data source or the logger. However, the interface makes no guarantees in which thread the call will be made. Furthermore, the logger may be call multiple times by different assembly processes (i.e., the logging messages from different sequences may be interleaved).
 * @next: (closure next_data) (scope notified): the next sequence handler
 * @logger: the logger callback
 */
PandaMux panda_mux_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogProxy logger);

/**
 * Create a new multiplexed reader for given to FASTQ streams.
 * @see panda_create_fastq_reader
 */
PandaMux panda_mux_new_fastq_reader(
	PandaNextChar forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaNextChar reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/**
 * Open a pair of bzipped files for multi-threaded assembly.
 *
 * @see panda_assembler_open_bz2
 */
PandaMux panda_mux_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/**
 * Open a pair of gzipped files for multi-threaded assembled.
 * @see panda_assembler_open_gz
 */
PandaMux panda_mux_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/* === Methods === */

/**
 * Create a new assembler using the multiplexer as it sequence source.
 *
 * The new assembler will draw sequences from the original source in a thread-safe way. Each assembler is not thread-safe. This means that, to use the interface correctly, one creates a sequence source, wraps it in a multiplexer, then creates an assembler for every thread. Each assembler should be accessed in only one thread. It may be advisable to create a single assembler and set its configuration, then copy the settings to subsequently created assemblers.
 * @see panda_assembler_copy_configuration
 */
PandaAssembler panda_mux_create_assembler(
	PandaMux mux);

/**
 * Create a new assembler using the multiplexer as it sequence source with a custom k-mer table size.
 * @see panda_mux_create_assembler
 * @see panda_assembler_new_kmer
 */
PandaAssembler panda_mux_create_assembler_kmer(
	PandaMux mux,
	size_t num_kmers);

/**
 * Increase the reference count on a multiplexer.
 */
PandaMux panda_mux_ref(
	PandaMux mux);

/**
 * Attached a callback for every sequence that fails to have an overlap.
 *
 * This will be called when a sequence fails to have an overlap computed. This does not include sequences that are missing primers or sequences that are assembled and discarded by modules.
 *
 * Concurrency will be handled by the mulitplexer; all calls to this function will be serialised.
 *
 * @handler: (closure handler_data) (scope notified): the callback for a failed pair
 */

void panda_mux_set_fail_alignment(
	PandaMux mux,
	PandaFailAlign handler,
	void *handler_data,
	PandaDestroy handler_destroy);

/**
 * Decrease the reference count on a multiplexer.
 * @mux: (transfer full): the mux to be released.
 */
void panda_mux_unref(
	PandaMux mux);

/**
 * Get the number of assemblers created so far.
 */
size_t panda_mux_get_child_count(
	PandaMux mux);

/**
 * The logging proxy used by this mux
 * Returns: (transfer none): the proxy
 */
PandaLogProxy panda_mux_get_loggger(
	PandaMux mux);

EXTERN_C_END
#endif
