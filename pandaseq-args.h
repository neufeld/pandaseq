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

#ifndef _PANDASEQ_ARGS_H
#        define _PANDASEQ_ARGS_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/**
 * Parse command line arguments to in order to construct assemblers.
 *
 * This is meant to serve as a general framework for parsing command line arguments with maximum code reuse. There are two kinds of arguments: assembler-only and general. Assembler-only arguments need no context (i.e., they on modify the assembler based on their argument). General arguments might do this or they might be involved in selecting the sequence source.
 *
 * Lists of command line arguments are passed in and parsed. The opener is then called to open the data source and provide a sequence source. Then, an assembler and multiplexer will be created. All the assembler-only arguments and any additional setup are applied to the assembler. Any needed modules are loaded. Finally, the assembler and multiplexer are output, for use by the thread pool code.
 * @args:(array length=args_length): the strings from the command line
 * @assembler_args:(array length=assembler_args_length): descriptors of all the assembler-only command line arguments
 * @general_args:(array length=general_args_length): descriptors of all the command line arguments userstood by the callbacks
 * @tweak:(closure user_data): a callback for every command line argument matching a general argument
 * @opener:(closure user_data): a callback to open the sequence source
 * @assembler_setup:(closure user_data) (allow-none): a callback to configure the assembler
 * @out_assembler:(out callee-allocates) (transfer full): the assembler constructed after argument parsing
 * @out_mux:(out callee-allocates) (transfer full): the multiplexer constructed after argument parsing
 * @out_threads:(out caller-allocates): the number of threads the user wishes to use
 * Returns: whether command line parsing was successful and the output parameters have been populated
 */
bool panda_parse_args(
	char *const *args,
	int args_length,
	const panda_tweak_assembler *const *const assembler_args,
	size_t assembler_args_length,
	const panda_tweak_general *const *const general_args,
	size_t general_args_length,
	PandaTweakGeneral tweak,
	PandaOpener opener,
	PandaSetup assembler_setup,
	void *user_data,
	PandaAssembler *out_assembler,
	PandaMux *out_mux,
	int *out_threads,
	PandaOutputSeq * output,
	void **output_data,
	PandaDestroy *output_destroy);

/**
 * Spawn threads and assemble sequences.
 * @threads: the number of threads to spawn
 * @assembler: (transfer full): the main assembler to use. If multiple threads are to be used, the configuration of this assembler will be copied to all the slave assemblers.
 * @mux: (transfer full) (allow-none): the multiplexer to use. If null, no threads will be created. The provided assembler must be a product of this multiplexer.
 * @output: (closure output_data) (scope notified): the function that will write assembled sequences to where they belong.
 */
bool panda_run_pool(
	int threads,
	PandaAssembler assembler,
	PandaMux mux,
	PandaOutputSeq output,
	void *output_data,
	PandaDestroy output_destroy);

/**
 * Construct a new list by combining existing lists.
 * @array:(allow-none)(array length=array_length): The storage location of the array. If the array is initially empty, this may be null. This must be from reallocable storage.
 * @additions:(array length_additions_length): The array whose items to copy.
 */
void panda_tweak_assembler_append(
	const panda_tweak_assembler ***array,
	size_t *length,
	const panda_tweak_assembler *const *additions,
	size_t additions_length);

/**
 * Sort a list of arguments into flag order.
 */
void panda_tweak_assembler_sort(
	const panda_tweak_assembler **array,
	size_t length);

/**
 * Construct a new list by combining existing lists.
 * @array:(allow-none)(array length=array_length): The storage location of the array. If the array is initially empty, this may be null. This must be from reallocable storage.
 * @additions:(array length_additions_length): The array whose items to copy.
 */
void panda_tweak_general_append(
	const panda_tweak_general ***array,
	size_t *length,
	const panda_tweak_general *const *additions,
	size_t additions_length);

/**
 * Sort a list of arguments into flag order.
 */
void panda_tweak_general_sort(
	const panda_tweak_general **array,
	size_t length);

/* === Standard Arguments === */

/**
 * The standard list of assembler-only arguments for PANDAseq binaries.
 */
PANDA_EXTERN const panda_tweak_assembler *const panda_stdargs[];
PANDA_EXTERN const size_t panda_stdargs_length;
/**
 * The strip primers after switch (-a).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_primers_after;
/**
 * The minimum length filter switch (-l).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_min_len;
/**
 * The maximum length filter switch (-L).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_max_len;
/**
 * The no N's switch (-N).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_degenerates;
/**
 * The minimum overlap switch (-o).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_min_overlap;
/**
 * The forward primer filter switch (-p).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_forward_primer;
/**
 * The reverse primer filter switch (-q).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_reverse_primer;
/**
 * The threshold filter switch (-t).
 */
PANDA_EXTERN const panda_tweak_assembler panda_stdargs_threshold;

/* === FASTQ file arguments === */

/**
 * Command line arguments for a pair of FASTQ files from Illumina.
 */
PANDA_EXTERN const panda_tweak_general *const panda_args_fastq_args[];
PANDA_EXTERN const size_t panda_args_fastq_args_length;

/**
 * Create a new argument handler.
 */
PandaArgsFastq panda_args_fastq_new(
	);

/**
 * Cleanup the argument handler.
 */
void panda_args_fastq_free(
	PandaArgsFastq data);

/**
 * Initialise the sequence stream for the FASTQ argument handler.
 */
PandaNextSeq panda_args_fastq_opener(
	PandaArgsFastq data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy);

/**
 * Do additional assembly setup for the FASTQ argument handler.
 */
bool panda_args_fastq_setup(
	PandaArgsFastq data,
	PandaAssembler assembler);

/**
 * Process the command line arguments for the FASTQ argument handler.
 */
bool panda_args_fastq_tweak(
	PandaArgsFastq data,
	char flag,
	const char *argument);

/* === Overhanging wrapper arguments === */

/**
 * Command line arguments for overhanging read pair trimmer.
 */
const panda_tweak_general **panda_args_hang_args(
	const panda_tweak_general *const *const general_args,
	size_t general_args_length,
	size_t *length);

/**
 * Create a new argument handler.
 */
PandaArgsHang panda_args_hang_new(
	void *user_data,
	PandaDestroy destroy,
	PandaTweakGeneral tweak,
	PandaOpener opener,
	PandaSetup setup);

/**
 * Cleanup the argument handler.
 */
void panda_args_hang_free(
	PandaArgsHang data);

/**
 * Initialise the sequence stream.
 */
PandaNextSeq panda_args_hang_opener(
	PandaArgsHang data,
	PandaLogProxy logger,
	PandaFailAlign *fail,
	void **fail_data,
	PandaDestroy *fail_destroy,
	void **next_data,
	PandaDestroy *next_destroy);

/**
 * Do additional assembly setup.
 */
bool panda_args_hang_setup(
	PandaArgsHang data,
	PandaAssembler assembler);

/**
 * Process the command line arguments.
 */
bool panda_args_hang_tweak(
	PandaArgsHang data,
	char flag,
	const char *argument);
EXTERN_C_END
#endif
