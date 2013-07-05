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

#ifndef _PANDASEQ_COMMON_H
#        define _PANDASEQ_COMMON_H
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
EXTERN_C_BEGIN
/**
 * Maximum length of a sequence
 */
#        define PANDA_MAX_LEN (panda_max_len())
#        define PANDA_TAG_LEN 50
extern size_t panda_max_len(
	void);

/* === Objects === */

/**
 * The manager for an assembly
 */
typedef struct panda_assembler *PandaAssembler;

/**
 * The standard argument handler for a pair of FASTQ files from Illumina.
 */
typedef struct panda_args_fastq *PandaArgsFastq;

/**
 * The standard argument handler for overhanging read pair trimmer.
 */
typedef struct panda_args_hang *PandaArgsHang;

/**
 * Iterate over a sequence presenting all k-mers without Ns or other denegerate bases.
 *
 * Iterators are not reference counted.
 */
typedef struct panda_iter *PandaIter;

/**
 * Logging proxy object.
 */
typedef struct panda_log_proxy *PandaLogProxy;

/**
 * Sequence validity checker
 */
typedef struct panda_module *PandaModule;

/**
 * A threading-safe wrapper to allow multiple assemblers to share a single data source.
 */
typedef struct panda_mux *PandaMux;

/**
 * A set of sequence identifiers against which to match.
 */
typedef struct panda_idset *PandaSet;

/* === Enum and Flags === */

/**
 * Codes used for logging conditions during the assembly.
 *
 * Some of these are errors and some are informational.
 */
typedef enum {
	PANDA_CODE_BAD_NT,
	PANDA_CODE_BEST_OVERLAP,
	PANDA_CODE_BUILD_FORWARD,
	PANDA_CODE_BUILD_OVERLAP,
	PANDA_CODE_BUILD_REVERSE,
	PANDA_CODE_FORWARD_KMER,
	PANDA_CODE_ID_PARSE_FAILURE,
	PANDA_CODE_INSUFFICIENT_KMER_TABLE,
	PANDA_CODE_LOST_KMER,
	PANDA_CODE_LOW_QUALITY_REJECT,
	PANDA_CODE_MISMATCHED_BASE,
	PANDA_CODE_MOD_INFO,
	PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH,
	PANDA_CODE_NO_DATA,
	PANDA_CODE_NO_FILE,
	PANDA_CODE_NO_FORWARD_PRIMER,
	PANDA_CODE_NO_QUALITY_INFO,
	PANDA_CODE_NO_REVERSE_PRIMER,
	PANDA_CODE_NOT_PAIRED,
	PANDA_CODE_OVERLAP_POSSIBILITY,
	PANDA_CODE_PARSE_FAILURE,
	PANDA_CODE_PREMATURE_EOF,
	PANDA_CODE_READ_TOO_LONG,
	PANDA_CODE_RECONSTRUCTION_PARAM,
	PANDA_CODE_REJECT_STAT,
	PANDA_CODE_REVERSE_KMER,
	PANDA_CODE_SEQUENCE_TOO_LONG,
	PANDA_CODE_PHRED_OFFSET,
} PandaCode;

/**
 * Decide what kinds of messages are passed to the logger.
 */
typedef unsigned int PandaDebug;

/**
 * The policy for Illumina tags/barcodes in sequence names.
 */
typedef enum {
	/**
	 * The parsing should return an error if the sequence does not have a tag.
	 */
	PANDA_TAG_PRESENT,
	/**
	 * The parsing should return an error if the sequence has a tag.
	 */
	PANDA_TAG_ABSENT,
	/**
	 * The parsing should not care if the sequence a tag.
	 */
	PANDA_TAG_OPTIONAL,
} PandaTagging;

/**
 * A single nucleotide
 */
typedef char panda_nt;

/* === Structures === */

/**
 * A k-mer and its position in the original sequence.
 */
typedef struct {
	size_t kmer;
	size_t posn;
} panda_kmer;

/**
 * A single nucleotide with quality information
 */
typedef struct {
	/**
	 * The nucleotide
	 */
	panda_nt nt;
	/**
	 * The quality score as a PHRED score
	 */
	char qual;
} panda_qual;

typedef struct {
	/**
	 * The nucleotide
	 */
	panda_nt nt;
	/**
	 * The quality score as a log probability
	 */
	double p;
} panda_result;

/**
 * Illumina sequence information from the FASTQ header
 */
typedef struct {
	char instrument[100];
	int run;
	char flowcell[100];
	int lane;
	int tile;
	int x;
	int y;
	char tag[PANDA_TAG_LEN];
} panda_seq_identifier;

/**
 * Describes a command line option.
 */
typedef struct panda_tweak_general {
	/**
	 * The command line option.
	 */
	char flag;
	/**
	 * Whether the flag needs to be specified.
	 *
	 * This is used in the help output only.
	 */
	bool optional;
	/**
	 * The name of the argument as it appears in the help. If null, the argument is assumed to be boolean.
	 */
	const char *takes_argument;
	/**
	 * The help information to display to the user.
	 */
	const char *help;
} panda_tweak_general;

/**
 * A reconstructed sequence with meta information
 */
typedef struct {
	/**
	 * Calculated quality score as the geometric mean of the product of the Illumina quality scores of the included bases.
	 *
	 * It will always be between 0 and 1.
	 */
	double quality;
	/**
	 * Number of uncalled bases in the sequence.
	 */
	size_t degenerates;
	/**
	 * The sequence identification information
	 */
	panda_seq_identifier name;
	/**
	 * The reconstructed sequence with quality information
	 */
	panda_result *sequence;
	size_t sequence_length;
	/**
	 * The original forward sequence
	 */
	panda_qual *forward;
	size_t forward_length;
	/**
	 * The original reverse sequence
	 */
	panda_qual *reverse;
	size_t reverse_length;

	/**
	 * The number of nucleotides clipped from the forward sequence
	 */
	size_t forward_offset;
	/**
	 * The number of nucleotides clipped from the reverse sequence
	 */
	size_t reverse_offset;
	/**
	 * The number of mismatches in the overlap region.
	 */
	size_t overlap_mismatches;
	/**
	 * The number of overlaps that were examined to determine the one finally used.
	 */
	size_t overlaps_examined;
	/**
	 * The overlap chosen.
	 */
	size_t overlap;
} panda_result_seq;

/* === Function Pointers === */

/**
 * Check a sequence after reconstruction for validity.
 */
typedef bool (
	*PandaCheck) (
	const panda_result_seq *sequence,
	void *user_data);

/**
 * Free user data
 *
 * Any method which takes user data with a function pointer will call a destroy function when the user data is not longer needed such that the memory can be freed, if necessary. A destroy function may always be null, in which case, the memory managment is the responsibility of the caller.
 */
typedef void (
	*PandaDestroy) (
	void *user_data);

/**
 * Handle a failed alignment
 *
 * This is called when an assembler fails to align a sequence because it can't compute a reasonable overlap.
 * @assembler: The assembler that made the attempt
 * @id: the sequence id of the failed pair
 * @forward: (array length=forward_length): the forward read
 * @reverse: (array length=reverse_length): the reverse read
 * @user_data: (closure): context data
 */
typedef void (
	*PandaFailAlign) (
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	void *user_data);

/**
 * Get the next character from a FASTQ file or EOF.
 *
 * For assembly from an alternate source of data, this function returns the next character in the stream.
 */
typedef int (
	*PandaNextChar) (
	void *user_data);

/**
 * A callback for iterating over the current modules.
 * @assembler: the assembler which is being queried
 * @module: the module selected
 * @rejected: the number of sequences rejected by this module in the context of the current assembler.
 * @data: (closure): some user context data provided
 * Returns: true to continue iterating, false to stop
 */
typedef bool (
	*PandaModuleCallback) (
	PandaAssembler assembler,
	PandaModule module,
	size_t rejected,
	void *data);

/**
 * Get the next sequence pair.
 *
 * For assembly from a non-FASTQ text source, this function can provide the next sequence. The function must provide the sequences and metadata for assembly by modifing the values of its parameters.
 * @id: (out caller-allocates): the identifier information for the sequence pair
 * @forward: (array length=forward_length) (allow-none): the location of the parsed sequence data of the forward read. This memory is not managed by the assembler.
 * @reverse: (array length=reverse_length) (allow-none): the location of the parsed sequence data of the reverse read. This memory is not managed by the assembler.
 * Returns: true if there is a sequence available. All the parameters must be set correctly. If false, no more sequences will be read and the values in the parameters are ignored.
 */
typedef bool (
	*PandaNextSeq) (
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	void *user_data);

/**
 * Write a finsihed sequence to an appropriate place.
 * @sequence: the sequence from assembly
 * @user_data: (closure): the context provided
 */
typedef bool (
	*PandaOutputSeq) (
	const panda_result_seq *sequence,
	void *user_data);

/**
 * Check a sequence before reconstruction for validity.
 * @forward: (array length=forward_length): The forward read.
 * @reverse: (array length=reverse_length): The reverse read.
 */
typedef bool (
	*PandaPreCheck) (
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	void *user_data);

/**
 * Printf-like function for output.
 *
 * @data: (closure): The context for the error logging.
 * @format: (printf-like): The format string for use by printf.
 */
typedef void (
	*PandaPrintf) (
	void *data,
	const char *format,
	...);

/**
 * Create a sequence reader after argument parsing.
 * 
 * This returns a sequence source and a failure handler so that assembly can proceed.
 * @logger: The logging proxy to use, if needed.
 * @fail: (closure fail_data) (transfer full) (allow-none) (out callee-allocates) (scope notified): the handler for any sequences which do not align, if desired.
 * @user_data:(closure): the context
 * Returns: (scope notified) (allow-none): the sequence source, or null to indicate a failure
 */
typedef PandaNextSeq (
	*PandaOpener) (
	void *user_data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy);

/**
 * Perform any modifications to the assembler after creation.
 */
typedef bool (
	*PandaSetup) (
	void *user_data,
	PandaAssembler assembler);

/**
 * Process a command-line flag specified by the user.
 * @assembler: the assembler to which to make the adjustments
 * @flag: the command line flag specified
 * @argument: (transfer full) (allow-none): the command line argument, or null if not set.
 * @is_set: whether the command line argument is given in the case of boolean flags.
 * Returns: whether the flag was parsed successfully
 */
typedef bool (
	*PandaTweakAssembler) (
	PandaAssembler assembler,
	char flag,
	char *argument,
	bool is_set);

/**
 * Process a command-line flag specified by the user.
 * @user_data: (closure): the context
 * @flag: the command line flag specified
 * @argument: the option passed specified with the flag, if requested.
 * Returns: whether the flag was parsed succesfully
 */
typedef bool (
	*PandaTweakGeneral) (
	void *user_data,
	char flag,
	const char *argument);

/* === Round 2 Structures === */

/**
 * Describes a command line option that can be applied to an assembler.
 */
typedef struct panda_tweak_assembler {
	/**
	 * The command line option.
	 */
	char flag;
	/**
	 *  The name of the argument as it appears in the help. If null, the argument is assumed to be boolean.
	 */
	const char *takes_argument;
	/**
	 * The description of the option.
	 */
	const char *help;
	/**
	 * The callback to make the appropriate changes to the assembler.
	 */
	PandaTweakAssembler setup;
} panda_tweak_assembler;
EXTERN_C_END
#endif
