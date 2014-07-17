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

#ifndef _PANDASEQ_H
#        define _PANDASEQ_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/*
 * While this is not a GLib-based file, it generally follows GLib-style
 * conventions, particularly in the documentation.
 * 
 * Every pointer is assumed to be not-null unless specified otherwise by
 * (allow-none).
 *
 * All return types and pointers are must be freed or unreferenced by the
 * caller unless marked as (transfer none).
 *
 * Arrays are generally followed by a length. If the array is being returned,
 * the array length is an out parameter.
 *
 * The major library interfaces are refence counted. When all copies are
 * unreferenced, using the appropriate function, it will be garbage collected.
 * If it needs to be stored in multiple places, it can simply be referenced.
 *
 * PANDAseq makes heavy use of closures which involve a function pointer and a
 * void pointer containing opaque data passed to the function pointer. If the
 * opaque data must persist beyond the life of the call, another function
 * pointer is required to clean up the data. If neither the closure nor the
 * clean up is necessary, these may be null.
 *
 * If compiled against pthreads, PANDAseq is relatively thread-safe. All
 * functions that change reference counts are may be called at any time from
 * any thread. All other functions are have no guaranteeds unless explicity
 * stated.
 *
 * See https://live.gnome.org/GObjectIntrospection/Annotations for more
 * information.
 */
#        define PANDA_API 3
/**
 * Find the best offset of a small sequence in a large sequence.
 * @threshold: the minimum log probability to match
 * @reverse: if false, scan the sequence from start to finish, else, scan in the opposite direction
 * @haystack: (array length=haystack_length): the sequence to be searched
 * @needle: (array length=needle_length): the sequence for which to look
 * Returns: 0 if the sequence is not found, or one more than the offset
 */
size_t panda_compute_offset_qual(
	double threshold,
	bool reverse,
	const panda_qual *haystack,
	size_t haystack_length,
	const panda_nt *needle,
	size_t needle_length);

/**
 * Find the best offset of a small sequence in a large sequence.
 * @threshold: the minimum log probability to match
 * @reverse: if false, scan the sequence from start to finish, else, scan in the opposite direction
 * @haystack: (array length=haystack_length): the sequence to be searched
 * @needle: (array length=needle_length): the sequence for which to look
 * Returns: 0 if the sequence is not found, or one more than the offset
 */
size_t panda_compute_offset_result(
	double threshold,
	bool reverse,
	const panda_result *haystack,
	size_t haystack_length,
	const panda_nt *needle,
	size_t needle_length);

/**
 * Create an object to read sequences from two character streams of FASTQ data
 *
 * @forward: (closure forward_data) (scope notified): the functions to provide the stream of forward characters. Every time a new data is required, forward(forward_data) is called. When the stream has retrned EOF or the assembler is deallocated, forward_destroy(forward_data) is called.
 * @reverse: (closure reverse_data) (scope notified): the same for the reverse
 * sequence.
 * @qualmin: the quality to subtract from the incoming file (usually 33 or 64, depending on CASAVA versi
n)
 * @policy: method to handle unbarcoded sequences
 * @logger: the logging function to use during assembly.
 * Returns: (closure user_data) (scope notified): The function to call.
 */
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
	PandaDestroy *destroy);

/**
 * Compare sequences assembled by two different assemblers.
 *
 * @reader: (closure reader_data): the source of the sequences.
 * @control: (closure control_data): the control assembly process.
 * @experiment: (closure experiment_data): the experiment assembly process.
 * @suppress_quality_diffs: consider nucleotides that have different quality scores to be identical.
 */
bool panda_diff(
	PandaNextSeq reader,
	void *reader_data,
	PandaAssemble control,
	void *control_data,
	PandaAssemble experiment,
	void *experiment_data,
	bool suppress_quality_diffs);

/**
 * Wraps an existing stream of reads and clips off reads that have the too-long overlap problem.
 * @inner:(closure inner_data) (scope notified): the stream to wrap.
 * @forward:(array length=forward_length): the sequence to trim from the forward read.
 * @reverse:(array length=reverse_length): the sequence to trim from the reverse read.
 * @skip: whether to try to assemble sequences that don't contain a trim sequence.
 * @threshold: a log probability threshold for cut-off alignment
 * Returns: (closure next_data) (scope notified): the sequence stream
 */
PandaNextSeq panda_trim_overhangs(
	PandaNextSeq inner,
	void *inner_data,
	PandaDestroy inner_destroy,
	PandaLogProxy logger,
	panda_nt *forward,
	size_t forward_length,
	panda_nt *reverse,
	size_t reverse_length,
	bool skip,
	double threshold,
	void **next_data,
	PandaDestroy *next_destroy);

/**
 * Compute log(1 - exp(p)) efficiently.
 *
 * See [[http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf|Mächler, 2012]].
 */
double panda_log1mexp(
	double p);

/**
 * The current version string
 */
char const *panda_version(
	void);

/* === Flags === */

/** Usual output about assembly. */
#        define PANDA_DEBUG_BUILD ((PandaDebug) 1)
/** Input processing-related errors. */
#        define PANDA_DEBUG_FILE ((PandaDebug) 2)
/** Extra statistics. */
#        define PANDA_DEBUG_STAT ((PandaDebug) 4)
/** Information about building the k-mer table (long and boring). */
#        define PANDA_DEBUG_KMER ((PandaDebug) 8)
/** Excruciating detail about the reconstruction. */
#        define PANDA_DEBUG_RECON ((PandaDebug) 16)
/** Bucket loads of data about mistatches. */
#        define PANDA_DEBUG_MISMATCH ((PandaDebug) 32)
#        define PANDA_DEBUG_DEFAULT (PANDA_DEBUG_BUILD | PANDA_DEBUG_FILE | PANDA_DEBUG_STAT)

/**
 * The current flags used be the assember to report errors. Some errors are always reported.
 */
PANDA_EXTERN PandaDebug panda_debug_flags;

/* === Methods (for things from elsewhere) === */

/**
 * The output-file representation of an error code.
 * Returns: (transfer none): The string representation
 */
char const *panda_code_str(
	PandaCode code);

/**
 * Convert the PHRED quality score to a log probability.
 */
double panda_quality_log_probability(
	const panda_qual *q);

/**
 * Convert the PHRED quality score to a probability.
 */
double panda_quality_probability(
	const panda_qual *q);

/**
 * Convert the probability to a PHRED quality score.
 */
char panda_result_phred(
	const panda_result *r);

/* === I/O Methods === */

/**
 * Open a file that might be uncompressed or compressed with gzip or bzip2.
 *
 * @file_name: the file to open
 * @logger: the logger to write data
 * Returns: (scope notified) (closure user_data): the buffer read function to use.
 */
PandaBufferRead panda_open_buffer(
	const char *file_name,
	PandaLogProxy logger,
	void **user_data,
	PandaDestroy *destroy);

/**
 * Open a pair of FASTQ files for reading.
 *
 * Files may be uncompressed, or compressed with gzip or bzip2.
 *
 * @forward: the forward filename
 * @reverse: the reverse filename
 * @logger: the logging function to use during assembly.
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 * Returns: (closure user_data) (scope notified): The function to call.
 */
PandaNextSeq panda_open_fastq(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy);

/**
 * Write an unassembled sequence to a FASTA file as a concatenated pair.
 */
void panda_output_fail(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	PandaWriter writer);

/**
 * Write an unassembled sequence to a FASTQ file as a concatenated pair.
 */
void panda_output_fail_qual(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	PandaWriter writer);

/**
 * Write an assembly to a FASTA file.
 */
bool panda_output_fasta(
	const panda_result_seq *sequence,
	PandaWriter writer);
/**
 * Write an assembly to a FASTQ file.
 */
bool panda_output_fastq(
	const panda_result_seq *sequence,
	PandaWriter writer);

PandaNextSeq panda_create_async_reader(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	size_t length,
	void **user_data,
	PandaDestroy *destroy);

	/**
	 * The first base in the result sequence in the overlap.
	 */
#        define PANDA_RESULT_OVERLAP_OFFSET(result) ((result)->forward_length - (result)->forward_offset - (result)->overlap)
	/**
	 * The first base in the forward sequence in the overlap.
	 */
#        define PANDA_RESULT_OVERLAP_FORWARD_OFFSET(result) ((result)->forward_length - (result)->overlap)
	/**
	 * The first base in the reverse sequence in the overlap.
	 */
#        define PANDA_RESULT_OVERLAP_REVERSE_OFFSET(result) ((result)->result_length - (result)->overlap)

/* === Convenience macro is for Vala === */
#        define PANDA_FAIL(file, append, user_data, destroy) (*(user_data) = fopen(file, append ? "a" : "w"), *(destroy) = fclose, *(user_data) == NULL ? NULL : (PandaFailAlign) panda_output_fail)
#        define PANDACONCATE(x,y) x ## y
#        define PANDACONCAT(x,y) PANDACONCATE(x, y)

/* === Everything else === */

#        include<pandaseq-algorithm.h>
#        include<pandaseq-args.h>
#        include<pandaseq-assembler.h>
#        include<pandaseq-iter.h>
#        include<pandaseq-linebuf.h>
#        include<pandaseq-log.h>
#        include<pandaseq-module.h>
#        include<pandaseq-nt.h>
#        include<pandaseq-seqid.h>
#        include<pandaseq-set.h>
#        include<pandaseq-writer.h>
EXTERN_C_END
#endif
