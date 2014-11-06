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

#ifndef _PANDASEQ_ASSEMBLER_H
#        define _PANDASEQ_ASSEMBLER_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <stdarg.h>
#        include <stdio.h>
#        include <stdbool.h>
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/* === Costructors === */
/**
 * Create a new assembler from a sequence source.
 *
 * @next: (closure next_data) (scope notified) (allow-none): the function to call to get the next sequence. The assembler does not manage the memory of the returned arrays, but assume it may use them until the next call of next(next_data) or next_destroy(next_data). When the assembler is destroyed, it will call next_destroy(next_data). If null, only panda_assembler_assemble may be used and not panda_assembler_next.
 * @logger: the function to call to report information to the user
 */
PandaAssembler panda_assembler_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogProxy logger);

/**
 * The default number of locations in the k-mer look up table.
 *
 * When attempting to align the sequences, the assembler will store the location of every k-mer in a table. If the same k-mer is present multiple times, only the first ones will be store until the table is full. If the sequences are highly repetitive, lost positions can prevent good alignments.
 */
#        define PANDA_DEFAULT_NUM_KMERS 2

/**
 * Create a new assembler from a sequence source with a custom k-mer table size.
 *
 * @next: (closure next_data) (scope notified) (allow-none): the function to call to get the next sequence. The assembler does not manage the memory of the returned arrays, but assume it may use them until the next call of next(next_data) or next_destroy(next_data). When the assembler is destroyed, it will call next_destroy(next_data). If null, only panda_assembler_assemble may be used and not panda_assembler_next.
 * @logger: the proxy to call to report information to the user
 * @num_kmers: the number of sequence locations for a particular k-mer. The default is PANDA_DEFAULT_NUM_KMERS. This should be small (no more than 10), or the k-mer table will be extremely large.
 * @see panda_assembler_new
 */
PandaAssembler panda_assembler_new_kmer(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogProxy logger,
	size_t num_kmers);

/**
 * Create a new assembler for given to FASTQ streams.
 * @forward: (closure forward_data) (scope notified): the functions to provide the stream of forward characters.
 * @reverse: (closure reverse_data) (scope notified): the same for the reverse sequence.
 * @qualmin: the quality to subtract from the incoming file (usually 33 or 64, depending on CASAVA version)
 * @policy: method to handle unbarcoded sequences
 * @logger: the logging function to use during assembly.
 * @see panda_create_fastq_reader
 */
PandaAssembler panda_assembler_new_fastq_reader(
	PandaBufferRead forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaBufferRead reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/**
 * Open a pair of FASTQ files for assembly.
 *
 * File may be plain text, or compressed with gzip or bzip2.
 *
 * @logger: The proxy for error logging.
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 */
PandaAssembler panda_assembler_open(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/* === Methods === */

/**
 * Add a module to this assembly process.
 *
 * Sequences will be checked using this module.
 * Returns: true if the module was successfully initialised and added. If the
 * module's command line arguments are not processed correctly, this will fail.
 */
bool panda_assembler_add_module(
	PandaAssembler assembler,
	PandaModule module);

/**
 * Add a collection of modules to this assembly process.
 * @modules: (array length=modules_length): the modules to add
 * Returns: the index of the last successfully added module
 */
size_t panda_assembler_add_modules(
	PandaAssembler assembler,
	PandaModule *modules,
	size_t modules_length);

/**
 * Assemble a single sequence pair not drawn from the sequence stream.
 *
 * This works exactly like panda_assembler_next, but instead of asking the PandaNextSeq for the data, it expects this information to be provided.
 * @id: the sequence identifier for this read pair
 * @forward: (array length=forward_length): the forward read
 * @reverse: (array length=reverse_length): the reverse read
 */
const panda_result_seq *panda_assembler_assemble(
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length);

/**
 * Clone the configuration of one assembler to another.
 *
 * This does not affect the sequence source or logging. Only the primers, trimming, and modules. The recorded statistics are separate.
 */
void panda_assembler_copy_configuration(
	PandaAssembler dest,
	PandaAssembler src);

/**
 * Review all the modules associated with an assembler.
 *
 * @assembler: the assembler who owns the modules.
 * @callback: (closure callback_data): the callback to invoke for each module.
 * Returns: true if callback has seen each module
 */
bool panda_assembler_foreach_module(
	PandaAssembler assembler,
	PandaModuleCallback callback,
	void *callback_data);

/**
 * Log the number of sequences rejected by each module.
 */
void panda_assembler_module_stats(
	PandaAssembler assembler);

/**
 * Assemble the next sequence from the input
 *
 * Returns: (transfer none) (allow-none): The next successfully assembled sequence sequences until one is assembled successfully or no more sequences are available from the input stream, after which it will return null.
 * The returned sequence becomes invalid after the next call or after calling panda_assembler_unref.
 */
const panda_result_seq *panda_assembler_next(
	PandaAssembler assembler);

/**
 * Increase the reference count on an assembler.
 *
 * This is thread-safe.
 */
PandaAssembler panda_assembler_ref(
	PandaAssembler assembler);
/**
 * Decrease the reference count on an assembler.
 *
 * This is thread-safe.
 * @assembler: (transfer full): The assembler to destroy.
 */
void panda_assembler_unref(
	PandaAssembler assembler);

/* === Getters and Setters === */

/**
 * The scoring algorithm used.
 */
PandaAlgorithm panda_assembler_get_algorithm(
	PandaAssembler assembler);
void panda_assembler_set_algorithm(
	PandaAssembler assembler,
	PandaAlgorithm algorithm);

/**
 * The number of sequences rejected because the reads are unsatisfactory in some way.
 */
long panda_assembler_get_bad_read_count(
	PandaAssembler assembler);

/**
 * The number of sequences processed so far.
 */
long panda_assembler_get_count(
	PandaAssembler assembler);

/**
 * Attached a callback for every sequence that fails to have an overlap.
 *
 * This will be called when a sequence fails to have an overlap computed. This does not include sequences that are missing primers or sequences that are assembled and discarded by modules.
 *
 * @handler: (closure handler_data) (scope notified): the callback for a failed pair
 */
void panda_assembler_set_fail_alignment(
	PandaAssembler assembler,
	PandaFailAlign handler,
	void *handler_data,
	PandaDestroy handler_destroy);

/**
 * The number of sequences rejected because the overlap could not be determined.
 */
long panda_assembler_get_failed_alignment_count(
	PandaAssembler assembler);

/**
 * The forward primer sequence to be stripped
 * 
 * This is mutually exclusive with forward_trim
 * Returns: (array length=length) (transfer none) (allow-none): If no primer sequence is set, this will return null and set length to 0.
 */
panda_nt *panda_assembler_get_forward_primer(
	PandaAssembler assembler,
	size_t *length);

/**
 * The forward primer sequence to be stripped
 * @sequence: (array length_length) (allow-none): The primer sequene.
 */
void panda_assembler_set_forward_primer(
	PandaAssembler assembler,
	panda_nt *sequence,
	size_t length);

/**
 * The amount of forward sequence to strip
 * 
 * This is mutually exclusive with forward_primer
 */
size_t panda_assembler_get_forward_trim(
	PandaAssembler assembler);
void panda_assembler_set_forward_trim(
	PandaAssembler assembler,
	size_t trim);

/**
 * Report the longest overlap assembled so far.
 */
size_t panda_assembler_get_longest_overlap(
	PandaAssembler assembler);

/**
 * The number of sequences rejected because the quality score is too low.
 */
long panda_assembler_get_low_quality_count(
	PandaAssembler assembler);

/**
 * The minimum overlap two sequences must have to be accepted. It must be greater than one.
 */
int panda_assembler_get_minimum_overlap(
	PandaAssembler assembler);
void panda_assembler_set_minimum_overlap(
	PandaAssembler assembler,
	int overlap);

/**
 * The maximum overlap two sequences must have to be accepted.
 */
int panda_assembler_get_maximum_overlap(
	PandaAssembler assembler);
void panda_assembler_set_maximum_overlap(
	PandaAssembler assembler,
	int overlap);

/**
 * The name assoicated with this assembler.
 *
 * This is shown in logging output.
 */
const char *panda_assembler_get_name(
	PandaAssembler assembler);
void panda_assembler_set_name(
	PandaAssembler assembler,
	const char *name);

/**
 * The number of sequences rejected because the forward primer could not be aligned.
 */
long panda_assembler_get_no_forward_primer_count(
	PandaAssembler assembler);
/**
 * The number of sequences rejected because the reverse primer could not be aligned.
 */
long panda_assembler_get_no_reverse_primer_count(
	PandaAssembler assembler);

/**
 * The size of the k-mer table in this assembler.
 */
size_t panda_assembler_get_num_kmer(
	PandaAssembler assembler);

/**
 * The number of sequences accepted.
 */
long panda_assembler_get_ok_count(
	PandaAssembler assembler);

/**
 * Report the number of sequences that have been assembled with the overlap specified.
 */
long panda_assembler_get_overlap_count(
	PandaAssembler assembler,
	size_t overlap);

/**
 * Whether to strip the primers before or after assembly.
 *
 * Stripping before is faster, but means that there can be no primer sequence in the opposite read.
 */
bool panda_assembler_get_primers_after(
	PandaAssembler assembler);
void panda_assembler_set_primers_after(
	PandaAssembler assembler,
	bool after);

/**
 * The reverse primer sequence to be stripped
 * 
 * This is mutually exclusive with reverse_trim
 * Returns: (array length=length) (transfer none) (allow-none): If no primer sequence is set, this will return null and set length to 0.
 */
panda_nt *panda_assembler_get_reverse_primer(
	PandaAssembler assembler,
	size_t *length);
/**
 * The reverse primer sequence to be stripped
 * @sequence: (array length_length) (allow-none): The primer sequene.
 */
void panda_assembler_set_reverse_primer(
	PandaAssembler assembler,
	panda_nt *sequence,
	size_t length);
/**
 * The amount of reverse sequence to strip
 * 
 * This is mutually exclusive with reverse_primer
 */
size_t panda_assembler_get_reverse_trim(
	PandaAssembler assembler);
void panda_assembler_set_reverse_trim(
	PandaAssembler assembler,
	size_t trim);

/**
 * The numer of sequences where all possible overlaps had to be examined, instead of a quick hashing.
 */
long panda_assembler_get_slow_count(
	PandaAssembler assembler);

/**
 * The minimum quality threshold to have an assembly accepted. Must be between 0 and 1, exclusive.
 */
double panda_assembler_get_threshold(
	PandaAssembler assembler);
void panda_assembler_set_threshold(
	PandaAssembler assembler,
	double threshold);

/**
 * The logging proxy used by this assembler
 * Returns: (transfer none): the proxy
 */
PandaLogProxy panda_assembler_get_logger(
	PandaAssembler assembler);

/**
 * The penalty for moving the primer down the window.
 */
double panda_assembler_get_primer_penalty(
	PandaAssembler assembler);
void panda_assembler_set_primer_penalty(
	PandaAssembler assembler,
	double threshold);

EXTERN_C_END
#endif
