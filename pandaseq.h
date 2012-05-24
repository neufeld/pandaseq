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
#define _PANDASEQ_H
#include <stdarg.h>
#include <stdio.h>
#include <stdbool.h>

/**
 * Maximum length of a sequence
 */
#define PANDA_MAX_LEN 300

typedef void (*PandaPrintf)(void *data, const char *format, ...);

const char *panda_version(void);

/**
 * Codes used for logging conditions during the assembly.
 *
 * Some of these are errors and some are informational.
 */
typedef enum {
	PANDA_CODE_API_VERSION/*(int api_version)*/,
	PANDA_CODE_BAD_NT/*(panda_seq_identifier *id, int character)*/,
	PANDA_CODE_BEST_OVERLAP/*(panda_seq_identifier *id, int overlap)*/,
	PANDA_CODE_BUILD_FORWARD/*(panda_seq_identifier *id, int index_in_assembly, int index_in_reverse, panda_result *nt)*/,
	PANDA_CODE_BUILD_OVERLAP/*(panda_seq_identifier *id, int index_in_assembly, int index_in_forward, int index_in_reverse, panda_result *nt, panda_nt *forward, panda_nt *reverse)*/,
	PANDA_CODE_BUILD_REVERSE/*(panda_seq_identifier *id, int index_in_assembly, int index_in_reverse, panda_result *nt)*/,
	PANDA_CODE_FORWARD_KMER/*(panda_seq_identifier *id, unsigned int kmer, ssize_t position)*/,
	PANDA_CODE_ID_PARSE_FAILURE/* (const char *id) */,
	PANDA_CODE_INSUFFICIENT_KMER_TABLE/*(panda_seq_identifer *id)*/,
	PANDA_CODE_LOST_KMER/*(panda_seq_identifier *id, unsigned int kmer, ssize_t position)*/,
	PANDA_CODE_LOW_QUALITY_REJECT/*(panda_seq_identifier *id, double quality, double threshold)*/,
	PANDA_CODE_MISMATCHED_BASE/*(panda_seq_identifier *id, int index_in_forward, int index_in_reverse, panda_qual *forward, panda_qual *reverse)*/,
	PANDA_CODE_MOD_INFO/*(PandaModule module)*/,
	PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH/*(panda_seq_identifier *id)*/,
	PANDA_CODE_NO_DATA/*(panda_seq_identifier *id)*/,
	PANDA_CODE_NO_FILE/*(const char *filename)*/,
	PANDA_CODE_NO_FORWARD_PRIMER/*(panda_seq_identifier *id)*/,
	PANDA_CODE_NO_QUALITY_INFO/*(panda_seq_identifier *id)*/,
	PANDA_CODE_NO_REVERSE_PRIMER/*(panda_seq_identifier *id)*/,
	PANDA_CODE_NOT_PAIRED/*(panda_seq_identifier *foward, panda_seq_identifier *reverse)*/,
	PANDA_CODE_OVERLAP_POSSIBILITY/*(panda_seq_identifier *id, size_t overlap, size_t matches, size_t mismatches, size_t unknowns, double probability)*/,
	PANDA_CODE_PARSE_FAILURE/*(panda_seq_identifier *id)*/,
	PANDA_CODE_PREMATURE_EOF/*(panda_seq_identifier *id)*/,
	PANDA_CODE_RECONSTRUCTION_PARAM/*(panda_seq_identifier *id, int overlap, int forward_unpaired, int reverse_unpaired)*/,
	PANDA_CODE_REJECT_STAT/*(PandaModule module, long rejected_sequence_count)*/,
	PANDA_CODE_REVERSE_KMER/*(panda_seq_identifier *id, unsigned int kmer, ssize_t position)*/,
	PANDA_CODE_SEQUENCE_TOO_LONG/*(panda_seq_identifier *id)*/,
} PandaCode;

/**
 * A single nucleotide
 */
typedef char panda_nt;
/**
 * Adenine
 */
#define PANDA_NT_A ((panda_nt)1)
/**
 * Cytosine
 */
#define PANDA_NT_C ((panda_nt)2)
/**
 * Guanine
 */
#define PANDA_NT_G ((panda_nt)4)
/**
 * Thyamine
 */
#define PANDA_NT_T ((panda_nt)8)
/**
 * Is nucleotide degenerate?
 */
#define PANDA_NT_IS_DEGN(v) (((((unsigned int)(v)) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf) != 1)
/**
 * Is nucleotide all possible values?
 */
#define PANDA_NT_IS_N(n) ((n) == (panda_nt)0x0F)
/**
 * Get the nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii(char c);
/**
 * Get the complement nucleotide code for an ASCII character in IUPAC
 */
panda_nt panda_nt_from_ascii_complement(char c);
/**
 * Convert a nucleotide to an IUPAC representation
 */
char panda_nt_to_ascii(panda_nt val);

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
double panda_quality_probability(const panda_qual *q);
/**
 * Convert the PHRED quality score to a log probability.
 */
double panda_quality_log_probability(const panda_qual *q);


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
#define PANDA_TAG_LEN 50
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
 */
void panda_seqid_xprint(const panda_seq_identifier *id, PandaPrintf xprintf, void *x);
/**
 * Write an Illumina header for a sequence identifer to a file
 */
void panda_seqid_print(const panda_seq_identifier *id, FILE *file);
/**
 * Create an Illumina header for a sequence identifier
 *
 * The return string must not be freed and subsequent calls will obliterate the previously returned string.
 */
/*@dependent@*/const char *panda_seqid_str(const panda_seq_identifier *id);
/**
 * Compare two Illumina headers
 */
bool panda_seqid_equal(const panda_seq_identifier *one, const panda_seq_identifier *two);
/**
 * Parse an Illumina header
 *
 * Fills `id` with the parse result. The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
 */
int panda_seqid_parse(panda_seq_identifier *id, char *input);

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
	panda_result sequence[PANDA_MAX_LEN];
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
} panda_result_seq;

/**
 * Sequence validity checker
 */
typedef /*@refcounted@*/ struct panda_module* PandaModule;

/**
 * Check a sequence after reconstruction for validity.
 */
typedef bool (*PandaCheck)(const panda_result_seq *sequence, void *user_data);
/**
 * Check a sequence before reconstruction for validity.
 */
typedef bool (*PandaPreCheck)(const panda_seq_identifier *id, const panda_qual *forward, size_t forward_length, const panda_qual *reverse, size_t reverse_length, void *user_data);
/**
 * Free user data
 *
 * Any method which takes user data with a function pointer will call a destroy function when the user data is not longer needed such that the memory can be freed, if necessary. A destroy function may always be null, in which case, the memory managment is the responsibility of the caller.
 */
typedef void (*PandaDestroy) (void *user_data);
/**
 * Log an error/event
 *
 * If the function returns false, assembly will be halted.
 * The variadic arguments will provide context based on the particular code passed.
 * @see PandaCode
 */
typedef bool (*PandaLogger)(void *user_data, PandaCode code, ...);

/**
 * Create a module given sequence checking parameters.
 *
 * @param name the name of the module, for user interaction
 * @param check the check function, which must not be null
 * @param precheck an optional check to be done before the module
 * @param cleanup an optional function to be called when this module is garbage collected
 * The user is responsible for managing the memory associated with user_data, but the cleanup function will always be called.
 */
/*@notnull@*/PandaModule panda_module_new(/*@notnull@*/char *name, /*@notnull@*/PandaCheck check, /*@null@*/PandaPreCheck precheck, /*@null@*/void *user_data, /*@null@*/PandaDestroy cleanup);
/**
 * Load a module from a string containg the module name and arguments.
 *
 * @param path the name or path to a module separated by LT_PATHSEP_CHAR and any arguments to the initialisation function of that module
 */
/*@null@*/PandaModule panda_module_load(/*@notnull@*/char *path);

/**
 * Increase the reference count on a module.
 */
/*@notnull@*/PandaModule panda_module_ref(/*@notnull@*/PandaModule module);
/**
 * Decrease the reference count on a module.
 */
void panda_module_unref(/*@notnull@*/PandaModule module);
/**
 * Get the name of a module.
 *
 * The string returned must NOT be freed.
 */
/*@dependent@*/const char *panda_module_get_name(/*@notnull@*/PandaModule module);
/**
 * Get the description of a module.
 *
 * This is only appropriate for loaded modules.
 * Possibly null. The string returned must NOT be freed.
 */
/*@dependent@*/const char *panda_module_get_description(/*@notnull@*/PandaModule module);
/**
 * Get the usage information (i.e., help text) of a module.
 *
 * This is only appropriate for loaded modules.
 * Possibly null. The string returned must NOT be freed.
 */
/*@dependent@*/const char *panda_module_get_usage(/*@notnull@*/PandaModule module);
/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules.
 * Possibly null. The string returned must NOT be freed.
 */
/*@dependent@*/const char *panda_module_get_version(/*@notnull@*/PandaModule module);
/**
 * Get the arguments passed on loading of a module of a module.
 *
 * This is only appropriate for loaded modules.
 * Possibly null. The string returned must NOT be freed.
 */
/*@dependent@*/const char *panda_module_get_args(/*@notnull@*/PandaModule module);
/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules. Modules constructed by panda_module_new will always return PANDA_API.
 */
int panda_module_get_api(/*@notnull@*/PandaModule module);

/**
 * The manager for an assembly
 */
typedef /*@refcounted@*/ struct panda_assembler *PandaAssembler;
/**
 * Open a pair of gzipped (or uncompressed files) for assembly.
 *
 * @param qualmin the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 */
/*@null@*/PandaAssembler panda_assembler_open_gz(/*@notnull@*/char *forward, /*@notnull@*/char *reverse, /*@notnull@*/PandaLogger logger, /*@null@*/void *logger_data, /*@null@*/PandaDestroy logger_destroy, unsigned char qualmin);
/**
 * Open a pair of bzipped for assembly.
 *
 * @param qualmin the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
 */
/*@null@*/PandaAssembler panda_assembler_open_bz2(/*@notnull@*/char *forward, /*@notnull@*/char *reverse, /*@notnull@*/PandaLogger logger, /*@null@*/void *logger_data, /*@null@*/PandaDestroy logger_destroy, unsigned char qualmin);

/**
 * Get the next character from a FASTQ file or EOF.
 *
 * For assembly from an alternate source of data, this function returns the next character in the stream.
 */
typedef int (*PandaNextChar)(/*@null@*/void* user_data);
/**
 * Get the next sequence pair.
 *
 * For assembly from a non-FASTQ text source, this function can provide the next sequence. The function must provide the sequences and metadata for assembly by modifing the values of its parameters.
 * @param id the identifier information for the sequence pair
 * @param forward the location of the parsed sequence data of the forward read. This memory is not managed by the assembler.
 * @param forward_length the number of nucleotides in the forward read
 * @param reverse the location of the parsed sequence data of the reverse read. This memory is not managed by the assembler.
 * @param reverse_length the number of nucleotides in the reverse read
 */
typedef bool (*PandaNextSeq)(/*@notnull@ @out@*/panda_seq_identifier *id, /*@notnull@ @out@*/panda_qual **forward, /*@notnull@ @out@*/size_t *forward_length, /*@notnull@ @out@*/panda_qual **reverse, /*@notnull@ @out@*/size_t *reverse_length, /*@null@*/void *user_data);
/**
 * Create an object to read sequences from two character streams of FASTQ data
 *
 * @param forward, forward_data, forward_destroy the functions to provide the stream of forward characters. Every time a new character is required, forward(forward_data) is called. When the stream has returned EOF or the assembler is deallocated, forward_destroy(forward_data) is called.
 * @param reverse, reverse_data, reverse_destroy the same for the reverse sequence.
 * @param logger, logger_data the logging function to use during assembly. The logging function will not be memory managed.
 * @param user_data where to store the user_data for this function
 * @param destroy where to store the destroy function for the user data
 */
PandaNextSeq panda_assembler_create_fastq_reader(/*@notnull@*/PandaNextChar forward, /*@null@*/void *forward_data, /*@null@*/PandaDestroy forward_destroy, /*@notnull@*/PandaNextChar reverse, /*@null@*/void *reverse_data, /*@null@*/PandaDestroy reverse_destroy, /*@notnull@*/PandaLogger logger, /*@null@*/void *logger_data, unsigned char qualmin, /*@notnull@ @out@*/void **user_data, /*@notnull@ @out@*/PandaDestroy *destroy);
/**
 * Create a new assembler for given to FASTQ streams.
 * @see panda_assembler_create_fastq_reader
 */
/*@notnull@*/PandaAssembler panda_assembler_new_fastq_reader(/*@notnull@*/PandaNextChar forward, /*@null@*/void *forward_data, /*@null@*/PandaDestroy forward_destroy, /*@notnull@*/PandaNextChar reverse, /*@null@*/void *reverse_data, /*@null@*/PandaDestroy reverse_destroy, /*@notnull@*/PandaLogger logger, /*@null@*/void *logger_data, /*@null@*/PandaDestroy logger_destroy, unsigned char qualmin);
/**
 * Create a new assembler from a sequence source.
 *
 * @param next, next_data, next_destroy the function to call to get the next sequence. The assembler does not manage the memory of the returned arrays, but assume it may use them until the next call of next(next_data) or next_destroy(next_data). When the assembler is destroy, it will call next_destroy(next_data). If null, only panda_assembler_assemble may be used and not panda_assembler_next.
 * @param logger, logger_data, logger_destroy the function to call to report information to the user
 */
/*@notnull@*/PandaAssembler panda_assembler_new(/*@null@*/PandaNextSeq next, /*@null@*/void *next_data, /*@null@*/PandaDestroy next_destroy, /*@notnull@*/PandaLogger logger, /*@null@*/void *logger_data, /*@null@*/PandaDestroy logger_destroy);

/**
 * Increase the reference count on an assembler.
 *
 * This is thread-safe.
 */
PandaAssembler panda_assembler_ref(/*@notnull@*/PandaAssembler assembler);
/**
 * Decrease the reference count on an assembler.
 *
 * This is thread-safe.
 */
void panda_assembler_unref(/*@notnull@*/PandaAssembler assembler);
/**
 * Add a module to this assembly process.
 *
 * Sequences will be checked using this module.
 */
void panda_assembler_add_module(/*@notnull@*/PandaAssembler assembler, /*@notnull@*/PandaModule module);

/**
 * The minimum error estimation in the sequence data (epsilon)
 */
double panda_assembler_get_error_estimation(/*@notnull@*/PandaAssembler assembler);
void panda_assembler_set_error_estimation(/*@notnull@*/PandaAssembler assembler, double q);

/**
 * The minimum overlap two sequences must have to be accepted. It must be greater than one.
 */
int panda_assembler_get_minimum_overlap(/*@notnull@*/PandaAssembler assembler);
void panda_assembler_set_minimum_overlap(/*@notnull@*/PandaAssembler assembler, int overlap);

/**
 * The minimum quality threshold to have an assembly accepted. Must be between 0 and 1, exclusive.
 */
double panda_assembler_get_threshold(/*@notnull@*/PandaAssembler assembler);
void panda_assembler_set_threshold(/*@notnull@*/PandaAssembler assembler, double threshold);

/**
 * The number of sequences rejected because the forward primer could not be aligned.
 */
long panda_assembler_get_no_forward_primer_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences rejected because the reverse primer could not be aligned.
 */
long panda_assembler_get_no_reverse_primer_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences accepted.
 */
long panda_assembler_get_ok_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences rejected because the quality score is too low.
 */
long panda_assembler_get_low_quality_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences rejected because they contain degenerate (N) bases.
 */
long panda_assembler_get_degenerate_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences rejected because the overlap could not be determined.
 */
long panda_assembler_get_failed_alignment_count(/*@notnull@*/PandaAssembler assembler);
/**
 * The number of sequences processed so far.
 */
long panda_assembler_get_count(/*@notnull@*/PandaAssembler assembler);
/**
 * Reject sequences with degenerate (N) bases.
 */
bool panda_assembler_get_disallow_degenerates(/*@notnull@*/PandaAssembler assembler);
void panda_assembler_set_disallow_degenerates(/*@notnull@*/PandaAssembler assembler, bool allow);
/**
 * The forward primer sequence to be stripped
 * 
 * This is mutually exclusive with forward_trim
 */
/*@null@*/panda_nt *panda_assembler_get_forward_primer(/*@notnull@*/PandaAssembler assembler, /*@notnull@*/size_t *length);
void panda_assembler_set_forward_primer(/*@notnull@*/PandaAssembler assembler, panda_nt *sequence, size_t length);
/**
 * The reverse primer sequence to be stripped
 * 
 * This is mutually exclusive with reverse_trim
 */
/*@null@*/panda_nt *panda_assembler_get_reverse_primer(/*@notnull@*/PandaAssembler assembler, /*@notnull@*/size_t *length);
void panda_assembler_set_reverse_primer(/*@notnull@*/PandaAssembler assembler, panda_nt *sequence, size_t length);
/**
 * The amount of forward sequence to strip
 * 
 * This is mutually exclusive with forward_primer
 */
size_t panda_assembler_get_forward_trim(PandaAssembler assembler);
void panda_assembler_set_forward_trim(PandaAssembler assembler, size_t trim);
/**
 * The amount of reverse sequence to strip
 * 
 * This is mutually exclusive with reverse_primer
 */
size_t panda_assembler_get_reverse_trim(PandaAssembler assembler);
void panda_assembler_set_reverse_trim(PandaAssembler assembler, size_t trim);

/**
 * Log the number of sequences rejected by each module.
 */
void panda_assembler_module_stats(/*@notnull@*/PandaAssembler assembler);
/**
 * Assemble the next sequence from the input
 *
 * This function will process sequences until one is assembled successfully or no more sequences are available from the input stream, after which it will return null.
 * The returned sequence becomes invalid after the next call or after calling panda_assembler_unref.
 */
/*@null@*/const panda_result_seq *panda_assembler_next(/*@notnull@*/PandaAssembler assembler);
/**
 * Assemble a single sequence pair not drawn from the sequence stream.
 *
 * This works exactly like panda_assembler_next, but instead of asking the PandaSeqNext for the data, it expects this information to be provided.
 */
/*@null@*/const panda_result_seq *panda_assembler_assemble(/*@notnull@*/PandaAssembler assembler, /*@notnull@*/panda_seq_identifier *id, /*@notnull@*/const panda_qual *forward, size_t forward_length, /*@notnull@*/const panda_qual *reverse, size_t reverse_length);
/**
 * Write an assembly to a FASTA file.
 */
bool panda_output_fasta(/*@notnull@*/const panda_result_seq *sequence, /*@notnull@*/FILE *file);
/**
 * Write an assembly to a FASTQ file.
 */
bool panda_output_fastq(/*@notnull@*/const panda_result_seq *sequence, /*@notnull@*/FILE *file);
/**
 * Write errors and information to a file.
 */
bool panda_logger_file(/*@notnull@*/FILE *file, PandaCode code, ...);

/**
 * Writes an error message using the supplied printf-like function.
 */
bool panda_logger_v(PandaPrintf xprintf, void *x, PandaCode code, va_list va);
/*
 * Convenience macro is for Vala
 */
#define PANDA_LOGGER(file, user_data, destroy) (*user_data = file, *destroy = NULL, (PandaLogger) panda_logger_file)

/* Convenience macros for creating modules */
#define PANDA_API 1
#define PANDACONCATE(x,y) x ## y
#define PANDACONCAT(x,y) PANDACONCATE(x, y)
#define PRECHECK bool PANDACONCAT(PANDASEQ_MODULE,_LTX_precheck) (const panda_seq_identifier *id, panda_qual *forward, size_t forward_length, panda_qual *reverse, size_t reverse_length)
#define CHECK int PANDACONCAT(PANDASEQ_MODULE,_LTX_api) = PANDA_API; bool PANDACONCAT(PANDASEQ_MODULE,_LTX_check) (const panda_result_seq *sequence)
#define INIT bool PANDACONCAT(PANDASEQ_MODULE,_LTX_init)(const char *args)
#define CLEANUP void PANDACONCAT(PANDASEQ_MODULE,_LTX_destroy)(void)
#define HELP(desc, usage) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_desc) = desc; const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_usage) = usage
#define VER_INFO(version) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_version) = version
#endif
