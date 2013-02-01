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
#include "config.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pandaseq.h"
#include "assembler.h"
#include "module.h"

PandaAssembler
panda_assembler_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogger logger,
	void *logger_data,
	PandaDestroy logger_destroy) {
	return panda_assembler_new_kmer(next, next_data, next_destroy, logger, logger_data, logger_destroy, PANDA_DEFAULT_NUM_KMERS);
}

PandaAssembler
panda_assembler_new_kmer(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogger logger,
	void *logger_data,
	PandaDestroy logger_destroy,
	size_t num_kmers) {
	PandaAssembler assembler = malloc(sizeof(struct panda_assembler));
	if (assembler == NULL) {
		if (next_destroy != NULL) {
			next_destroy(next_data);
		}
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		return NULL;
	}
	assembler->refcnt = 1;
	assembler->next = next;
	assembler->next_data = next_data;
	assembler->next_destroy = next_destroy;
	assembler->logger = logger;
	assembler->logger_data = logger_data;
	assembler->logger_destroy = logger_destroy;
	assembler->noalgn = NULL;
	assembler->noalgn_data = NULL;
	assembler->noalgn_destroy = NULL;
	assembler->rejected = NULL;
	assembler->modules = NULL;
	assembler->modules_length = 0;
	assembler->modules_size = 0;
	assembler->result.forward = NULL;
	assembler->forward_primer_length = 0;
	assembler->result.reverse = NULL;
	assembler->reverse_primer_length = 0;
	assembler->forward_trim = 0;
	assembler->reverse_trim = 0;
	assembler->nofpcount = 0;
	assembler->norpcount = 0;
	assembler->okcount = 0;
	assembler->lowqcount = 0;
	assembler->degencount = 0;
	assembler->noalgncount = 0;
	assembler->badreadcount = 0;
	assembler->slowcount = 0;
	assembler->count = 0;
	assembler->no_n = false;
	assembler->post_primers = false;
	assembler->num_kmers = num_kmers;
	assert(1 << (8 * sizeof(seqindex)) > PANDA_MAX_LEN);
	assembler->kmerseen = malloc(KMERSEEN_SIZE(num_kmers));
	if (assembler->kmerseen == NULL) {
		if (next_destroy != NULL) {
			next_destroy(next_data);
		}
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		free(assembler);
		return NULL;
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&assembler->mutex, NULL);
#endif
	memset(assembler->kmerseen, 0, KMERSEEN_SIZE(num_kmers));
	panda_assembler_set_error_estimation(assembler, 0.36);
	panda_assembler_set_threshold(assembler, 0.6);
	panda_assembler_set_minimum_overlap(assembler, 2);
	return assembler;
}

void
panda_assembler_copy_configuration(
	PandaAssembler dest,
	PandaAssembler src) {
	int it;
	for (it = 0; it < src->modules_length; it++) {
		panda_assembler_add_module(dest, src->modules[it]);
	}
	panda_assembler_set_forward_primer(dest, src->forward_primer, src->forward_primer_length);
	panda_assembler_set_reverse_primer(dest, src->reverse_primer, src->reverse_primer_length);
	dest->forward_trim = src->forward_trim;
	dest->reverse_trim = src->reverse_trim;
	dest->threshold = src->threshold;
	dest->minoverlap = src->minoverlap;
	dest->no_n = src->no_n;
	dest->post_primers = src->post_primers;
}

int
panda_assembler_get_minimum_overlap(
	PandaAssembler assembler) {
	return assembler->minoverlap;
}

void
panda_assembler_set_minimum_overlap(
	PandaAssembler assembler,
	int overlap) {
	if (overlap > 1 && overlap < PANDA_MAX_LEN) {
		assembler->minoverlap = overlap;
	}
}

double
panda_assembler_get_threshold(
	PandaAssembler assembler) {
	return exp(assembler->threshold);
}

void
panda_assembler_set_threshold(
	PandaAssembler assembler,
	double threshold) {
	if (threshold > 0 && threshold < 1) {
		assembler->threshold = log(threshold);
	}
}

double
panda_assembler_get_error_estimation(
	PandaAssembler assembler) {
	return assembler->q;
}

void
panda_assembler_set_error_estimation(
	PandaAssembler assembler,
	double q) {
	if (q > 0 && q < 1) {
#ifdef HAVE_PTHREAD
		pthread_mutex_lock(&assembler->mutex);
#endif
		assembler->q = q;
		assembler->pmatch = log(0.25 * (1 - 2 * q + q * q));
		assembler->pmismatch = log((3 * q - 2 * q * q) / 18.0);
#ifdef HAVE_PTHREAD
		pthread_mutex_unlock(&assembler->mutex);
#endif
	}
}

long
panda_assembler_get_no_forward_primer_count(
	PandaAssembler assembler) {
	return assembler->nofpcount;
}

long
panda_assembler_get_no_reverse_primer_count(
	PandaAssembler assembler) {
	return assembler->norpcount;
}

long
panda_assembler_get_ok_count(
	PandaAssembler assembler) {
	return assembler->okcount;
}

long
panda_assembler_get_low_quality_count(
	PandaAssembler assembler) {
	return assembler->lowqcount;
}

long
panda_assembler_get_bad_read_count(
	PandaAssembler assembler) {
	return assembler->badreadcount;
}

long
panda_assembler_get_slow_count(
	PandaAssembler assembler) {
	return assembler->slowcount;
}

long
panda_assembler_get_degenerate_count(
	PandaAssembler assembler) {
	return assembler->degencount;
}

long
panda_assembler_get_failed_alignment_count(
	PandaAssembler assembler) {
	return assembler->noalgncount;
}

long
panda_assembler_get_count(
	PandaAssembler assembler) {
	return assembler->count;
}

bool
panda_assembler_get_disallow_degenerates(
	PandaAssembler assembler) {
	return assembler->no_n;
}

void
panda_assembler_set_disallow_degenerates(
	PandaAssembler assembler,
	bool allow) {
	assembler->no_n = allow;
}

PandaAssembler
panda_assembler_ref(
	PandaAssembler assembler) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&assembler->mutex);
#endif
	assembler->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&assembler->mutex);
#endif
	return assembler;
}

void
panda_assembler_unref(
	PandaAssembler assembler) {
	int count;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&assembler->mutex);
#endif
	count = --(assembler->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&assembler->mutex);
#endif
	if (count == 0) {
#ifdef HAVE_PTHREAD
		pthread_mutex_destroy(&assembler->mutex);
#endif
		free(assembler->kmerseen);
		module_destroy(assembler);
		DESTROY_MEMBER(assembler, next);
		DESTROY_MEMBER(assembler, logger);
		DESTROY_MEMBER(assembler, noalgn);
		free(assembler);
	}
}

PandaAssembler
panda_assembler_new_fastq_reader(
	PandaNextChar forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaNextChar reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogger logger,
	void *logger_data,
	PandaDestroy logger_destroy,
	unsigned char qualmin,
	PandaTagging policy) {
	void *user_data;
	PandaDestroy destroy;
	PandaNextSeq next;
	next = panda_create_fastq_reader(forward, forward_data, forward_destroy, reverse, reverse_data, reverse_destroy, logger, logger_data, qualmin, policy, &user_data, &destroy);
	return panda_assembler_new(next, user_data, destroy, logger, logger_data, logger_destroy);
}

void
panda_assembler_set_forward_primer(
	PandaAssembler assembler,
	panda_nt *sequence,
	size_t length) {
	size_t it;
	if (length < PANDA_MAX_LEN) {
		for (it = 0; it < length; it++) {
			assembler->forward_primer[it] = sequence[it];
		}
		assembler->forward_primer_length = length;
		assembler->forward_trim = 0;
	}
}

void
panda_assembler_set_reverse_primer(
	PandaAssembler assembler,
	panda_nt *sequence,
	size_t length) {
	size_t it;
	if (length < PANDA_MAX_LEN) {
		for (it = 0; it < length; it++) {
			assembler->reverse_primer[it] = sequence[it];
		}
		assembler->reverse_primer_length = length;
		assembler->reverse_trim = 0;
	}
}

panda_nt *
panda_assembler_get_reverse_primer(
	PandaAssembler assembler,
	size_t *length) {
	*length = assembler->reverse_primer_length;
	return *length == 0 ? NULL : assembler->reverse_primer;
}

panda_nt *
panda_assembler_get_forward_primer(
	PandaAssembler assembler,
	size_t *length) {
	*length = assembler->forward_primer_length;
	return *length == 0 ? NULL : assembler->forward_primer;
}

size_t
panda_assembler_get_forward_trim(
	PandaAssembler assembler) {
	return assembler->forward_trim;
}

void
panda_assembler_set_forward_trim(
	PandaAssembler assembler,
	size_t trim) {
	assembler->forward_trim = trim;
	assembler->forward_primer_length = 0;
}

bool
panda_assembler_get_primers_after(
	PandaAssembler assembler) {
	return assembler->post_primers;
}

void
panda_assembler_set_primers_after(
	PandaAssembler assembler,
	bool after) {
	assembler->post_primers = after;
}

size_t
panda_assembler_get_reverse_trim(
	PandaAssembler assembler) {
	return assembler->reverse_trim;
}

void
panda_assembler_set_reverse_trim(
	PandaAssembler assembler,
	size_t trim) {
	assembler->reverse_trim = trim;
	assembler->reverse_primer_length = 0;
}

void
panda_assembler_set_fail_alignment(
	PandaAssembler assembler,
	PandaFailAlign handler,
	void *handler_data,
	PandaDestroy handler_destroy) {
	DESTROY_MEMBER(assembler, noalgn);
	assembler->noalgn = handler;
	assembler->noalgn_data = handler_data;
	assembler->noalgn_destroy = handler_destroy;
}
