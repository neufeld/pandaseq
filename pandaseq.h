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
#        include <stdarg.h>
#        include <stdio.h>
#        include <stdbool.h>
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
/**
 * Maximum length of a sequence
 */
#        define PANDA_MAX_LEN (panda_max_len())
extern size_t panda_max_len(
	void);
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
 * The current version string
 */
const char *panda_version(
	void);
/**
 * The current module API version of the running library
 */
int panda_api_version(
	void);

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
 * The output-file representation of an error code.
 * Returns: (transfer none): The string representation
 */
char const *const panda_code_str(
	PandaCode code);

/**
 * Decide what kinds of messages are passed to the logger.
 */
typedef unsigned int PandaDebug;

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
extern PandaDebug panda_debug_flags;

/**
 * A single nucleotide
 */
typedef char panda_nt;
/**
 * Nothing (invalid nucleotide)
 */
#        define PANDA_NT_Z ((panda_nt)0)
/**
 * Adenine
 */
#        define PANDA_NT_A ((panda_nt)1)
/**
 * Cytosine
 */
#        define PANDA_NT_C ((panda_nt)2)
/**
 * Guanine
 */
#        define PANDA_NT_G ((panda_nt)4)
/**
 * Thyamine
 */
#        define PANDA_NT_T ((panda_nt)8)
/**
 * Is nucleotide degenerate?
 */
#        define PANDA_NT_IS_DEGN(v) (((((unsigned int)(v)) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf) != 1)
/**
 * Is nucleotide all possible values?
 */
#        define PANDA_NT_IS_N(n) ((n) == (panda_nt)0x0F)
/**
 * Get the nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii(
	char c);
/**
 * Get the complement nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii_complement(
	char c);
/**
 * Convert a nucleotide to an IUPAC representation
 */
char panda_nt_to_ascii(
	panda_nt val);

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
/**
 * Convert the PHRED quality score to a probability.
 */
double panda_quality_probability(
	const panda_qual *q);
/**
 * Convert the PHRED quality score to a log probability.
 */
double panda_quality_log_probability(
	const panda_qual *q);

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

#        define PANDA_TAG_LEN 50

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
 * Write the Illumina header to a printf-like function
 * @xprintf: (closure x): The callback to accept the input.
 */
void panda_seqid_xprint(
	const panda_seq_identifier *id,
	PandaPrintf xprintf,
	void *x);
/**
 * Write an Illumina header for a sequence identifier to a file
 */
void panda_seqid_print(
	const panda_seq_identifier *id,
	FILE *file);
/**
 * Create an Illumina header for a sequence identifier
 * @id: (allow-none): The identifer to be formatted
 * Returns: (transfer none): Subsequent calls will obliterate the previously returned string.
 */
const char *panda_seqid_str(
	const panda_seq_identifier *id);
/**
 * Compare two Illumina headers
 */
bool panda_seqid_equal(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two);

/**
 * Order two Illumina headers
 */
int panda_seqid_compare(
	const panda_seq_identifier *one,
	const panda_seq_identifier *two);

/**
 * Reset a sequnce identifier.
 * @id: (out caller-allocates): The structure to clear.
 */
void panda_seqid_clear(
	panda_seq_identifier *id);

/**
 * Parse an Illumina header
 *
 * @id: (out caller-allocates): The structure to fill with the parse result.
 * Returns: The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
 */
int panda_seqid_parse(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy);

/**
 * Parse the Illumina header
 *
 * @id: (out caller-allocates): The structure to fill with the parse result.
 * Returns: The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
 * @old: (out): Whether the sequence is from CASAVA 1.3-1.5 or not.
 * @end_ptr: (out) (transfer none): The point in the input where parsing stopped. If parsing was successful, this will be the end of the string.
 * @see panda_seqid_parse
 */
int panda_seqid_parse_fail(
	panda_seq_identifier *id,
	const char *input,
	PandaTagging policy,
	bool *old,
	const char **end_ptr);

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

/**
 * Sequence validity checker
 */
typedef struct panda_module *PandaModule;

/**
 * Check a sequence after reconstruction for validity.
 */
typedef bool (
	*PandaCheck) (
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
 * Free user data
 *
 * Any method which takes user data with a function pointer will call a destroy function when the user data is not longer needed such that the memory can be freed, if necessary. A destroy function may always be null, in which case, the memory managment is the responsibility of the caller.
 */
typedef void (
	*PandaDestroy) (
	void *user_data);
/**
 * Log an error/event
 *
 * @message: (allow-none): The error message, for the user.
 * Returns: If the function returns false, assembly will be halted.
 * @see PandaCode
 */
typedef bool (
	*PandaLogger) (
	PandaCode code,
	panda_seq_identifier *id,
	const char *message,
	void *user_data);

/**
 * Logging proxy object.
 */
typedef struct panda_log_proxy *PandaLogProxy;

/**
 * Create a new proxy with a callback.
 */
PandaLogProxy panda_log_proxy_new(
	PandaLogger log,
	void *log_data,
	PandaDestroy log_destroy);

/**
 * Create a new proxy to standard error.
 */
PandaLogProxy panda_log_proxy_new_stderr(
	);

/**
 * Increase the reference count on a proxy.
 */
PandaLogProxy panda_log_proxy_ref(
	PandaLogProxy proxy);
/**
 * Decrease the reference count on a proxy.
 * @module: (transfer full): the proxy to release.
 */
void panda_log_proxy_unref(
	PandaLogProxy proxy);

bool panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	panda_seq_identifier *id,
	const char *message);

/**
 * Create a module given sequence checking parameters.
 *
 * @name: the name of the module, for user interaction
 * @check: (closure user_data): the function to be run after assembly
 * @precheck: (closure user_data): a function to be run before assembly
 * @user_data: (transfer full): the context data for the functions. The user is responsible for managing the memory associated with user_data, but the cleanup function will always be called to do so.
 * @cleanup: (closure user_data): a function to be called when this module is garbage collected
 */
PandaModule panda_module_new(
	const char *name,
	PandaCheck check,
	PandaPreCheck precheck,
	void *user_data,
	PandaDestroy cleanup);
/**
 * Load a module from a string containg the module name and arguments.
 *
 * @path: the name or path to a module separated by LT_PATHSEP_CHAR and any arguments to the initialisation function of that module
 */
PandaModule panda_module_load(
	const char *path);

/**
 * Increase the reference count on a module.
 */
PandaModule panda_module_ref(
	PandaModule module);
/**
 * Decrease the reference count on a module.
 * @module: (transfer full): the module to release.
 */
void panda_module_unref(
	PandaModule module);
/**
 * Get the name of a module.
 *
 * Returns: (transfer none): the module's name
 */
const char *panda_module_get_name(
	PandaModule module);
/**
 * Get the description of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The description help text.
 */
const char *panda_module_get_description(
	PandaModule module);
/**
 * Get the usage information (i.e., help text) of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_usage(
	PandaModule module);
/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_version(
	PandaModule module);
/**
 * Get the arguments passed on loading of a module of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_args(
	PandaModule module);
/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules. Modules constructed by panda_module_new will always return PANDA_API.
 */
int panda_module_get_api(
	PandaModule module);

/**
 * The manager for an assembly
 */
typedef struct panda_assembler *PandaAssembler;
/**
 * Open a pair of gzipped (or uncompressed files) for assembly.
 * @logger: The proxy for error logging.
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 */
PandaAssembler panda_assembler_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);
/**
 * Open a pair of bzipped files for assembly.
 *
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 */
PandaAssembler panda_assembler_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy);

/**
 * Get the next character from a FASTQ file or EOF.
 *
 * For assembly from an alternate source of data, this function returns the next character in the stream.
 */
typedef int (
	*PandaNextChar) (
	void *user_data);
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
 * Create an object to read sequences from two character streams of FASTQ data
 *
 * @forward: (closure forward_data) (scope notified): the functions to provide the stream of forward characters. Every time a new character is required, forward(forward_data) is called. When the stream has returned EOF or the assembler is deallocated, forward_destroy(forward_data) is called.
 * @reverse: (closure reverse_data) (scope notified): the same for the reverse sequence.
 * @qualmin: the quality to subtract from the incoming file (usually 33 or 64, depending on CASAVA version)
 * @policy: method to handle unbarcoded sequences
 * @logger: the logging function to use during assembly.
 * Returns: (closure user_data) (scope notified): The function to call.
 */
PandaNextSeq panda_create_fastq_reader(
	PandaNextChar forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaNextChar reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy);
/**
 * Open a pair of gzipped (or uncompressed files).
 *
 * @forward: the forward filename
 * @reverse: the reverse filename
 * @logger: the logging function to use during assembly.
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 * Returns: (closure user_data) (scope notified): The function to call.
 */
PandaNextSeq panda_open_gz(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy);
/**
 * Open a pair of bzipped files.
 *
 * @forward: the forward filename
 * @reverse: the reverse filename
 * @logger: the logging function to use during assembly.
 * @qualmin: the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 * Returns: (closure user_data) (scope notified): The function to call.
 */
PandaNextSeq panda_open_bz2(
	const char *forward,
	const char *reverse,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy,
	void **user_data,
	PandaDestroy *destroy);
/**
 * Create a new assembler for given to FASTQ streams.
 * @forward: (closure forward_data) (scope notified): the functions to provide the stream of forward characters. Every time a new character is required, forward(forward_data) is called. When the stream has returned EOF or the assembler is deallocated, forward_destroy(forward_data) is called.
 * @reverse: (closure reverse_data) (scope notified): the same for the reverse sequence.
 * @qualmin: the quality to subtract from the incoming file (usually 33 or 64, depending on CASAVA version)
 * @policy: method to handle unbarcoded sequences
 * @logger: the logging function to use during assembly.
 * @see panda_create_fastq_reader
 */
PandaAssembler panda_assembler_new_fastq_reader(
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
 * Clone the configuration of one assembler to another.
 *
 * This does not affect the sequence source or logging. Only the primers, trimming, and modules. The recorded statistics are separate.
 */
void panda_assembler_copy_configuration(
	PandaAssembler dest,
	PandaAssembler src);
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
 * The minimum error estimation in the sequence data (epsilon)
 */
double panda_assembler_get_error_estimation(
	PandaAssembler assembler);
void panda_assembler_set_error_estimation(
	PandaAssembler assembler,
	double q);

/**
 * The minimum overlap two sequences must have to be accepted. It must be greater than one.
 */
int panda_assembler_get_minimum_overlap(
	PandaAssembler assembler);
void panda_assembler_set_minimum_overlap(
	PandaAssembler assembler,
	int overlap);

/**
 * The minimum quality threshold to have an assembly accepted. Must be between 0 and 1, exclusive.
 */
double panda_assembler_get_threshold(
	PandaAssembler assembler);
void panda_assembler_set_threshold(
	PandaAssembler assembler,
	double threshold);

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
 * The number of sequences accepted.
 */
long panda_assembler_get_ok_count(
	PandaAssembler assembler);
/**
 * The number of sequences rejected because the quality score is too low.
 */
long panda_assembler_get_low_quality_count(
	PandaAssembler assembler);
/**
 * The number of sequences rejected because the reads are unsatisfactory in some way.
 */
long panda_assembler_get_bad_read_count(
	PandaAssembler assembler);
/**
 * The numer of sequences where all possible overlaps had to be examined, instead of a quick hashing.
 */
long panda_assembler_get_slow_count(
	PandaAssembler assembler);
/**
 * The number of sequences rejected because they contain degenerate (N) bases.
 */
long panda_assembler_get_degenerate_count(
	PandaAssembler assembler);
/**
 * The number of sequences rejected because the overlap could not be determined.
 */
long panda_assembler_get_failed_alignment_count(
	PandaAssembler assembler);
/**
 * The number of sequences processed so far.
 */
long panda_assembler_get_count(
	PandaAssembler assembler);
/**
 * Reject sequences with degenerate (N) bases.
 */
bool panda_assembler_get_disallow_degenerates(
	PandaAssembler assembler);
void panda_assembler_set_disallow_degenerates(
	PandaAssembler assembler,
	bool allow);
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
 * Log the number of sequences rejected by each module.
 */
void panda_assembler_module_stats(
	PandaAssembler assembler);

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
 * Assemble the next sequence from the input
 *
 * Returns: (transfer none) (allow-none): The next successfully assembled sequence sequences until one is assembled successfully or no more sequences are available from the input stream, after which it will return null.
 * The returned sequence becomes invalid after the next call or after calling panda_assembler_unref.
 */
const panda_result_seq *panda_assembler_next(
	PandaAssembler assembler);
/**
 * Assemble a single sequence pair not drawn from the sequence stream.
 *
 * This works exactly like panda_assembler_next, but instead of asking the PandaSeqNext for the data, it expects this information to be provided.
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
 * Report the number of sequences that have been assembled with the overlap specified.
 */
long panda_assembler_get_overlap_count(
	PandaAssembler assembler,
	size_t overlap);

/**
 * Report the longest overlap assembled so far.
 */
size_t panda_assembler_get_longest_overlap(
	PandaAssembler assembler);

/**
 * Write an assembly to a FASTA file.
 */
bool panda_output_fasta(
	const panda_result_seq *sequence,
	FILE *file);
/**
 * Write an assembly to a FASTQ file.
 */
bool panda_output_fastq(
	const panda_result_seq *sequence,
	FILE *file);
/**
 * Write errors and information to a file.
 */
bool panda_logger_file(
	PandaCode code,
	panda_seq_identifier *id,
	const char *message,
	FILE *file);

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
	FILE *file);

/**
 * A threading-safe wrapper to allow multiple assemblers to share a single data source.
 */
typedef struct panda_mux *PandaMux;
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
 * Increase the reference count on a multiplexer.
 */
PandaMux panda_mux_ref(
	PandaMux mux);
/**
 * Decrease the reference count on a multiplexer.
 * @mux: (transfer full): the mux to be released.
 */
void panda_mux_unref(
	PandaMux mux);
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
 * Open a pair of gzipped files for multi-threaded assembled.
 * @see panda_assembler_open_gz
 */
PandaMux panda_mux_open_gz(
	const char *forward,
	const char *reverse,
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
	panda_qual *haystack,
	size_t haystack_length,
	panda_nt *needle,
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
	panda_result *haystack,
	size_t haystack_length,
	panda_nt *needle,
	size_t needle_length);

/**
 * Compute log(1 - exp(p)) efficiently.
 *
 * See [[http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf|Mächler, 2012]].
 */
double panda_log1mexp(
	double p);

/**
 * A set of sequence identifiers against which to match.
 */
typedef struct panda_idset *PandaSet;

/**
 * Create a new, empty set.
 */
PandaSet panda_idset_new(
	void);
/**
 * Increase the reference count on a set.
 *
 * This is thread-safe.
 */
PandaSet panda_idset_ref(
	PandaSet set);
/**
 * Decrease the reference count on a set.
 *
 * This is thread-safe.
 * @set: (transfer full): the mux to be released.
 */
void panda_idset_unref(
	PandaSet set);

/**
 * Add a sequence identifier to a set.
 */
void panda_idset_add(
	PandaSet set,
	const panda_seq_identifier *id);
/**
 * Parse a sequence identifier and add it to the set.
 * @id: the text id to parse
 * @old: (out): Whether the sequence is from CASAVA 1.3-1.5 or not.
 * @end_ptr: (out) (transfer none): The point in the input where parsing stopped. If parsing was successful, this will be the end of the string.
 * Returns: true on success
 * @see panda_seqid_parse_fail
 */
bool panda_idset_add_str(
	PandaSet set,
	const char *id,
	PandaTagging policy,
	bool *old,
	const char **end_ptr);
/**
 * Check if a sequence identifier has been added to the set.
 */
bool panda_idset_contains(
	PandaSet set,
	const panda_seq_identifier *id);

/**
 * A k-mer and its position in the original sequence.
 */
typedef struct {
	size_t kmer;
	size_t posn;
} panda_kmer;
/**
 * Iterate over a sequence presenting all k-mers without Ns or other denegerate bases.
 *
 * Iterators are not reference counted.
 */
typedef struct panda_iter *PandaIter;
/**
 * Destroy an iterator.
 */
void panda_iter_free(
	PandaIter iter);
/**
 * Copy an iterator to a new one, preserving its current state.
 */
PandaIter panda_iter_dup(
	PandaIter iter);
/**
 * Set an iterator back to the beginning of the sequence.
 */
void panda_iter_reset(
	PandaIter iter);
/**
 * Get the k-mer length for the iterator.
 */
int panda_iter_k(
	PandaIter iter);
/**
 * Get the number of useful bits in the output.
 *
 * This is the maximum value of panda_kmer.kmer for this iterator.
 */
size_t panda_iter_bits(
	PandaIter iter);
/**
 * Advance to the next position in the sequence.
 * Returns: (allow-none) (transfer none): if null, there are no more k-mers in the sequence
 */
const panda_kmer *panda_iter_next(
	PandaIter iter);
/**
 * Create an iterator over a sequence of nucleotides.
 * @seq: (array length=seq_length) (scope container): the sequence to iterate over, and its length. This sequence must not be freed during the life of the iterator.
 * @reverse: true to iterate from the end of the sequence rather than the beginning
 * @k: the length of the output words. This must range between 1 and 4 * sizeof(size_t). Any other values will be converted to the standard k-mer length of 8.
 */
PandaIter panda_iterate_nt(
	panda_nt *seq,
	size_t seq_length,
	bool reverse,
	int k);
/**
 * Iterate over quality-annotated sequence.
 * @see panda_iterate_nt
 */
PandaIter panda_iterate_qual(
	panda_qual *seq,
	size_t seq_length,
	bool reverse,
	int k);
/**
 * Iterate over probability-annotated sequence.
 * @see panda_iterate_nt
 */
PandaIter panda_iterate_result(
	panda_result *seq,
	size_t seq_length,
	bool reverse,
	int k);

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
 * Spawn threads and assemble sequences.
 * @threads: the number of threads to spawn
 * @assembler: (transfer full): the main assembler to use. If multiple threads are to be used, the configuration of this assembler will be copied to all the slave assemblers.
 * @mux: (transfer full) (allow-none): the multiplexer to use. If null, no threads will be created. The provided assembler must be a product of this multiplexer.
 * @output: (closure output_data) (scope notified): the function that will write assembled sequences to where they belong.
 */
int panda_run_pool(
	int threads,
	PandaAssembler assembler,
	PandaMux mux,
	PandaOutputSeq output,
	void *output_data,
	PandaDestroy output_destroy);

/*
 * Convenience macro is for Vala
 */
#        define PANDA_FAIL(file, append, user_data, destroy) (*(user_data) = fopen(file, append ? "a" : "w"), *(destroy) = fclose, *(user_data) == NULL ? NULL : (PandaFailAlign) panda_output_fail)

#        define PANDA_API 2
#        define PANDACONCATE(x,y) x ## y
#        define PANDACONCAT(x,y) PANDACONCATE(x, y)

/* Convenience macros for creating modules */
#        ifdef PANDASEQ_MODULE
#                define PRECHECK bool PANDACONCAT(PANDASEQ_MODULE,_LTX_precheck) (const panda_seq_identifier *id, panda_qual *forward, size_t forward_length, panda_qual *reverse, size_t reverse_length)
#                define CHECK bool PANDACONCAT(PANDASEQ_MODULE,_LTX_check) (const panda_result_seq *sequence)
#                define INIT bool PANDACONCAT(PANDASEQ_MODULE,_LTX_init)(const char *args)
#                define CLEANUP void PANDACONCAT(PANDASEQ_MODULE,_LTX_destroy)(void)
#                define HELP(desc, usage) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_desc) = desc; const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_usage) = usage
#                define VER_INFO(version) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_version) = version
#        endif
EXTERN_C_END
#endif
