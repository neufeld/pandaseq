/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2013  Andre Masella

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
#define _POSIX_C_SOURCE 2
#include "config.h"
#include<stdlib.h>
#include<time.h>
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq.h"
#include "assembler.h"
#include "misc.h"

struct shared_info {
	MANAGED_MEMBER(
		PandaOutputSeq,
		output);
	volatile bool some_seqs;
	time_t starttime;
#ifdef HAVE_PTHREAD
	pthread_mutex_t output_mutex;
	pthread_mutex_t stderr_mutex;
#endif
};

struct thread_info {
	struct shared_info *shared;
	pthread_t thread;
	PandaAssembler assembler;
	size_t index;
};

#define STAT(name, type, val)	PANDACONCAT(panda_log_proxy_stat_, type)(info->assembler->logger, info->assembler, name, (val))

static void printtime(
	struct thread_info *info,
	long count) {
	time_t now;
	(void) time(&now);
	STAT("TIME", str, ctime(&now));
	STAT("ELAPSED", long,
		 (int) (now - info->shared->starttime));
	STAT("READS", long,
		count);
}

static void *do_assembly(
	struct thread_info *info) {
	long count;
	const panda_result_seq *result;

	while ((result = panda_assembler_next(info->assembler)) != NULL) {
		count = panda_assembler_get_count(info->assembler);
		if (count % 1000 == 0) {
#ifdef HAVE_PTHREAD
			pthread_mutex_lock(&info->shared->stderr_mutex);
#endif
			printtime(info, count);
#ifdef HAVE_PTHREAD
			pthread_mutex_unlock(&info->shared->stderr_mutex);
#endif
		}
#ifdef HAVE_PTHREAD
		pthread_mutex_lock(&info->shared->output_mutex);
#endif
		info->shared->output(result, info->shared->output_data);
#ifdef HAVE_PTHREAD
		pthread_mutex_unlock(&info->shared->output_mutex);
#endif
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&info->shared->stderr_mutex);
#endif
	count = panda_assembler_get_count(info->assembler);
	if (count > 0) {
		info->shared->some_seqs = true;
	}
	printtime(info, count);
	if (panda_assembler_get_forward_primer(info->assembler, NULL) != NULL)
		STAT("NOFP", long,
			panda_assembler_get_no_forward_primer_count(info->assembler));
	if (panda_assembler_get_reverse_primer(info->assembler, NULL) != NULL)
		STAT("NORP", long,
			panda_assembler_get_no_reverse_primer_count(info->assembler));
	STAT("NOALGN", long,
		panda_assembler_get_failed_alignment_count(info->assembler));
	STAT("LOWQ", long,
		panda_assembler_get_low_quality_count(info->assembler));
	STAT("BADR", long,
		panda_assembler_get_bad_read_count(info->assembler));
	STAT("SLOW", long,
		panda_assembler_get_slow_count(info->assembler));
	if (panda_assembler_get_disallow_degenerates(info->assembler))
		STAT("DEGENERATE", long,
			panda_assembler_get_degenerate_count(info->assembler));
	panda_assembler_module_stats(info->assembler);
	STAT("OK", long,
		panda_assembler_get_ok_count(info->assembler));

	panda_log_proxy_write_overlap(info->assembler->logger, info->assembler);

	panda_assembler_unref(info->assembler);
	return NULL;
}

bool panda_run_pool(
	int threads,
	PandaAssembler assembler,
	PandaMux mux,
	PandaOutputSeq output,
	void *output_data,
	PandaDestroy output_destroy) {
	size_t it;
	struct thread_info self;
	struct shared_info shared_info;
	struct thread_info *thread_list;

	if (assembler == NULL)
		return false;

	shared_info.some_seqs = false;
	shared_info.output = output;
	shared_info.output_data = output_data;
	shared_info.output_destroy = output_destroy;
	(void) time(&shared_info.starttime);

#if HAVE_PTHREAD
	pthread_mutex_init(&shared_info.output_mutex, NULL);
	pthread_mutex_init(&shared_info.stderr_mutex, NULL);
	if (threads > 1 && mux != NULL) {

		thread_list = calloc(sizeof(struct thread_info), threads - 1);
		for (it = 0; it < threads - 1; it++) {
			thread_list[it].assembler = panda_mux_create_assembler(mux);
			thread_list[it].index = it + 1;
			thread_list[it].shared = &shared_info;
			if (thread_list[it].assembler == NULL) {
				fprintf(stderr, "ERR\tMUXCREATE\t%d\n", (int) it + 1);
				threads = it + 1;
				break;
			}
			panda_assembler_copy_configuration(thread_list[it].assembler, assembler);
			if (pthread_create(&thread_list[it].thread, NULL, (void *(*)(void *)) do_assembly, &thread_list[it]) != 0) {
				fprintf(stderr, "ERR\tPCREATE\t%d\n", (int) it + 1);
				threads = it + 1;
				break;
			}
		}
	}
	panda_mux_unref(mux);
#endif
	self.shared = &shared_info;
	self.index = 0;
	self.assembler = assembler;
	do_assembly(&self);
#if HAVE_PTHREAD
	if (threads > 1 && mux != NULL) {
		for (it = 0; it < threads - 1; it++) {
			pthread_join(thread_list[it].thread, NULL);
		}
		free(thread_list);
	}
	pthread_mutex_destroy(&shared_info.output_mutex);
	pthread_mutex_destroy(&shared_info.stderr_mutex);
#endif
	DESTROY_MEMBER(&shared_info, output);
	return shared_info.some_seqs;
}
