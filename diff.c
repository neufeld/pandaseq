/* PANDAseq-Diff -- Check differences between two versions of PANDAseq.
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
#include "config.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pandaseq.h"
#include "misc.h"

bool panda_diff(
	PandaNextSeq reader,
	void *reader_data,
	PandaAssemble control,
	void *control_data,
	PandaAssemble experiment,
	void *experiment_data,
	bool suppress_quality_diffs) {

	size_t length_diffs = 0;
	size_t nt_diffs = 0;
	size_t diffs_better_score = 0;
	size_t diffs_worse_score = 0;
	size_t gained = 0;
	size_t lost = 0;
	size_t total = 0;

	const panda_qual *forward;
	size_t forward_length;
	const panda_qual *reverse;
	size_t reverse_length;
	panda_seq_identifier id;

	while (reader(&id, &forward, &forward_length, &reverse, &reverse_length, reader_data)) {
		const panda_result_seq *control_result = control(control_data, &id, forward, forward_length, reverse, reverse_length);
		const panda_result_seq *experiment_result = experiment(experiment_data, &id, forward, forward_length, reverse, reverse_length);
		total++;
		if (control_result == NULL && experiment_result == NULL) {
			/* Both fail. That's a match. */
			continue;
		} else if (control_result == NULL || experiment_result == NULL) {
			panda_seqid_print(&id, stdout);
			if (control_result == NULL) {
				gained++;
				fprintf(stdout, " has been gained.\n");
			} else {
				lost++;
				fprintf(stdout, " has been lost.\n");
			}
			continue;
		}

		if (control_result->quality < experiment_result->quality) {
			diffs_better_score++;
		} else if (control_result->quality > experiment_result->quality) {
			diffs_worse_score++;
		}

		if (experiment_result->sequence_length != control_result->sequence_length) {
			length_diffs++;
			panda_seqid_print(&id, stdout);
			fprintf(stdout, " differ in length %zd → %zd.\n", control_result->sequence_length, experiment_result->sequence_length);
		} else {
			bool nt_diff = false;
			size_t it;
			for (it = 0; it < experiment_result->sequence_length; it++) {
				if (control_result->sequence[it].nt != experiment_result->sequence[it].nt) {
					panda_seqid_print(&id, stdout);
					fprintf(stdout, " differ at nucleotide %zd, %c → %c.\n", it, panda_nt_to_ascii(control_result->sequence[it].nt), panda_nt_to_ascii(experiment_result->sequence[it].nt));
					nt_diff = true;
				} else if (control_result->sequence[it].p != experiment_result->sequence[it].p && !suppress_quality_diffs) {
					panda_seqid_print(&id, stdout);
					fprintf(stdout, " differ at nucleotide %zd (%c), quality %f → %f.\n", it, (int) panda_nt_to_ascii(control_result->sequence[it].nt), exp(control_result->sequence[it].p), exp(experiment_result->sequence[it].p));
					nt_diff = true;
				}
			}
			if (nt_diff) {
				nt_diffs++;
			}
		}
	}
	fprintf(stdout, "%zd sequences compared.\n%zd scored better.\n%zd scored worse.\n%zd changed (%zd length changed, %zd sequence changed).\n%zd gained.\n%zd lost.\n", total, diffs_better_score, diffs_worse_score, nt_diffs + length_diffs, length_diffs, nt_diffs, gained, lost);
	return (total == 0 || length_diffs > 0 || nt_diffs > 0 || gained > 0 || lost > 0);
}

struct data {
	bool help;
	size_t num_kmers;
	PandaTweakGeneral general;
	void *general_data;
	bool verbose;
};

static const panda_tweak_general kmers = {.flag = 'k',.optional = true,.takes_argument = "kmers",.help = "The number of k-mers in the table.", false };

static const panda_tweak_general help = {.flag = 'h',.optional = true,.takes_argument = NULL,.help = "Show this delightful nonsense.", false };
static const panda_tweak_general verbose = {.flag = 'v',.optional = true,.takes_argument = NULL,.help = "Be more verbose and show differences in per-base quality.", false };

static const panda_tweak_general *common_args[] = {
	&help,
	&kmers,
	&verbose
};

static const panda_tweak_general *unique_args[] = {
	&kmers,
};

static bool common_tweak_general(
	void *user_data,
	char flag,
	const char *argument) {

	struct data *data = (struct data *) user_data;
	long int value;

	switch (flag) {
	case 'h':
		data->help = true;
		return true;
	case 'k':
		errno = 0;
		value = strtol(argument, NULL, 10);
		if (errno != 0 || value < 0 || (size_t) value > PANDA_MAX_LEN) {
			fprintf(stderr, "Bad k-mer list length.\n");
			return false;
		}
		data->num_kmers = (size_t) value;
		return true;
	case 'v':
		data->verbose = true;
		return true;
	default:
		return data->general(data->general_data, flag, argument);
	}
}

#define CLEANUP() for (it = 0; it < options_used; it++) if(options[it].arg != NULL) free(options[it].arg)
static PandaAssembler prep_assembler(
	const char *name,
	PandaLogProxy logger,
	char *const **args,
	int *args_length,
	const panda_tweak_assembler *const *const assembler_args,
	size_t assembler_args_length,
	panda_tweak_assembler_opt * common_options,
	size_t common_options_length,
	size_t num_kmers,
	struct data *data,
	PandaSetup assembler_setup) {

	panda_tweak_assembler_opt options[50];
	size_t options_used;
	PandaAssembler assembler = NULL;
	int args_unused;
	size_t it;

	data->num_kmers = 0;

	if (!panda_dispatch_args(*args, *args_length, assembler_args, assembler_args_length, unique_args, sizeof(unique_args) / sizeof(unique_args[0]), common_tweak_general, data, options, sizeof(options) / sizeof(options[0]), &options_used, &args_unused)) {
		CLEANUP();
		return NULL;
	}

	assembler = panda_assembler_new_kmer(NULL, NULL, NULL, logger, (data->num_kmers > 0) ? data->num_kmers : num_kmers);
	if (assembler == NULL) {
		CLEANUP();
		return NULL;
	}
	for (it = 0; it < common_options_length; it++) {
		char *arg = malloc(strlen(common_options[it].arg) + 1);
		memcpy(arg, common_options[it].arg, strlen(common_options[it].arg) + 1);
		fprintf(stderr, "Applying common flag -%c %s to %s assembler.\n", common_options[it].tweak->flag, common_options[it].arg, name);
		if (!(common_options[it].tweak->setup) (assembler, common_options[it].tweak->flag, arg)) {
			CLEANUP();
			panda_assembler_unref(assembler);
			return NULL;
		}
	}
	for (it = 0; it < options_used; it++) {
		char *arg = options[it].arg;
		fprintf(stderr, "Applying flag -%c %s to %s assembler.\n", options[it].tweak->flag, options[it].arg, name);
		options[it].arg = NULL;
		if (!(options[it].tweak->setup) (assembler, options[it].tweak->flag, arg)) {
			CLEANUP();
			panda_assembler_unref(assembler);
			return NULL;
		}
	}
	if (assembler_setup != NULL && !assembler_setup(data->general_data, assembler)) {
		CLEANUP();
		panda_assembler_unref(assembler);
		return NULL;
	}

	*args_length -= args_unused - 1;
	*args += args_unused - 1;

	CLEANUP();
	return assembler;
}

#undef CLEANUP

#define CLEANUP() for (it = 0; it < options_used; it++) if(options[it].arg != NULL) free(options[it].arg); panda_log_proxy_unref(logger); panda_writer_unref(writer); free(combined_general_args)

bool panda_diff_parse_args(
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
	PandaAssembler *out_control_assembler,
	PandaAssembler *out_experimental_assembler,
	PandaNextSeq *next,
	void **next_data,
	PandaDestroy *next_destroy,
	bool *out_suppress_quality_diffs) {

	PandaAssembler assembler;
	int args_unused;
	const panda_tweak_general **combined_general_args = NULL;
	size_t combined_general_args_length = 0;
	struct data data;
	PandaLogProxy logger = NULL;
	size_t it;
	size_t num_kmers;
	MANAGED_STACK(PandaFailAlign,
		fail);
	panda_tweak_assembler_opt options[50];
	size_t options_used;
	PandaWriter writer = panda_writer_new_null();

	data.general = tweak;
	data.general_data = user_data;
	data.help = false;
	data.num_kmers = 0;
	data.verbose = false;

	MAYBE(out_control_assembler) = NULL;
	MAYBE(out_experimental_assembler) = NULL;
	MAYBE(out_suppress_quality_diffs) = false;
	*next = NULL;
	*next_data = NULL;
	*next_destroy = NULL;

	/* Process command line arguments. */
	panda_tweak_general_append(&combined_general_args, &combined_general_args_length, general_args, general_args_length);
	panda_tweak_general_append(&combined_general_args, &combined_general_args_length, common_args, sizeof(common_args) / sizeof(common_args[0]));
	if (!panda_dispatch_args(args, args_length, assembler_args, assembler_args_length, combined_general_args, combined_general_args_length, common_tweak_general, &data, options, sizeof(options) / sizeof(options[0]), &options_used, &args_unused)) {
		CLEANUP();
		return false;
	}

	DESTROY_STACK(fail);
	args_length -= args_unused - 1;
	args += args_unused - 1;

	logger = panda_log_proxy_new(writer);
	if (data.help || (*next = opener(user_data, logger, &fail, &fail_data, &fail_destroy, next_data, next_destroy)) == NULL) {
		panda_args_help(args[0], assembler_args, assembler_args_length, combined_general_args, combined_general_args_length);
		fprintf(stderr, "\n\nRepeat arguments in blocks: %s common arguments -- control arguments -- experimental arguments\n", args[0]);
		CLEANUP();
		return false;
	}

	num_kmers = data.num_kmers;
	if (num_kmers == 0) {
		num_kmers = PANDA_DEFAULT_NUM_KMERS;
	}
	assembler = prep_assembler("control", logger, &args, &args_length, assembler_args, assembler_args_length, options, options_used, num_kmers, &data, assembler_setup);
	if (assembler == NULL) {
		CLEANUP();
		return false;
	}
	MAYBE(out_control_assembler) = assembler;
	assembler = prep_assembler("experimental", logger, &args, &args_length, assembler_args, assembler_args_length, options, options_used, num_kmers, &data, assembler_setup);
	if (assembler == NULL) {
		CLEANUP();
		return false;
	}
	MAYBE(out_experimental_assembler) = assembler;
	MAYBE(out_suppress_quality_diffs) = !data.verbose;
	CLEANUP();
	return true;
}
