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
#define __USE_UNIX98 1
#define _XOPEN_SOURCE 500
#include "config.h"
#if HAVE_PTHREAD
#        include <pthread.h>
#        include <stdlib.h>
#        include <string.h>
#        include "pandaseq.h"
#        include "pandaseq-mux.h"
#        include "assembler.h"
#        include "buffer.h"

struct panda_mux {
	pthread_mutex_t mutex;
	pthread_mutex_t next_mutex;
	pthread_rwlock_t noalgn_rwlock;
	PandaLogProxy logger;
	 MANAGED_MEMBER(
		PandaNextSeq,
		next);
	 MANAGED_MEMBER(
		PandaFailAlign,
		noalgn);
	volatile size_t refcnt;
	volatile size_t child_count;
};

PandaMux panda_mux_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogProxy logger) {
	PandaMux mux = malloc(sizeof(struct panda_mux));
	mux->refcnt = 1;
	mux->next = next;
	mux->next_data = next_data;
	mux->next_destroy = next_destroy;
	mux->logger = panda_log_proxy_ref(logger);
	mux->noalgn = NULL;
	mux->noalgn_data = NULL;
	mux->noalgn_destroy = NULL;
	mux->child_count = 0;
	pthread_mutex_init(&mux->mutex, NULL);
	pthread_mutex_init(&mux->next_mutex, NULL);
	pthread_rwlock_init(&mux->noalgn_rwlock, NULL);
	return mux;
}

PandaMux panda_mux_new_fastq_reader(
	PandaBufferRead forward,
	void *forward_data,
	PandaDestroy forward_destroy,
	PandaBufferRead reverse,
	void *reverse_data,
	PandaDestroy reverse_destroy,
	PandaLogProxy logger,
	unsigned char qualmin,
	PandaTagging policy) {
	void *user_data;
	PandaDestroy destroy;
	PandaNextSeq next;
	next = panda_create_fastq_reader(forward, forward_data, forward_destroy, reverse, reverse_data, reverse_destroy, logger, qualmin, policy, NULL, NULL, NULL, &user_data, &destroy);
	return panda_mux_new(next, user_data, destroy, logger);
}

size_t child_count(
	PandaMux mux) {
	size_t count;
	pthread_mutex_lock(&mux->mutex);
	count = mux->child_count++;
	pthread_mutex_unlock(&mux->mutex);
	return count;
}

PandaMux panda_mux_ref(
	PandaMux mux) {
	pthread_mutex_lock(&mux->mutex);
	mux->refcnt++;
	pthread_mutex_unlock(&mux->mutex);
	return mux;
}

void panda_mux_unref(
	PandaMux mux) {
	size_t count;
	if (mux == NULL)
		return;
	pthread_mutex_lock(&mux->mutex);
	count = --(mux->refcnt);
	pthread_mutex_unlock(&mux->mutex);
	if (count == 0) {
		pthread_mutex_destroy(&mux->mutex);

		panda_log_proxy_unref(mux->logger);

		pthread_mutex_lock(&mux->next_mutex);
		DESTROY_MEMBER(mux, next);
		pthread_mutex_unlock(&mux->next_mutex);
		pthread_mutex_destroy(&mux->next_mutex);

		pthread_rwlock_wrlock(&mux->noalgn_rwlock);
		DESTROY_MEMBER(mux, noalgn);
		pthread_rwlock_unlock(&mux->noalgn_rwlock);
		pthread_rwlock_destroy(&mux->noalgn_rwlock);
		free(mux);
	}
}

struct mux_data {
	PandaMux mux;
	panda_qual forward[MAX_LEN];
	panda_qual reverse[MAX_LEN];
};

static bool mux_next(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	struct mux_data *data) {
	bool result;
	const panda_qual *common_forward;
	const panda_qual *common_reverse;
	pthread_mutex_lock(&data->mux->next_mutex);

	result = data->mux->next(id, &common_forward, forward_length, &common_reverse, reverse_length, data->mux->next_data);
	if (common_forward == NULL || *forward_length == 0) {
		*forward = NULL;
		*forward_length = 0;
	} else {
		*forward = data->forward;
		memcpy(*forward, common_forward, sizeof(panda_qual) * *forward_length);
	}
	if (common_reverse == NULL || *reverse_length == 0) {
		*reverse = NULL;
		*reverse_length = 0;
	} else {
		*reverse = data->reverse;
		memcpy(*reverse, common_reverse, sizeof(panda_qual) * *reverse_length);
	}
	pthread_mutex_unlock(&data->mux->next_mutex);
	return result;
}

void mux_free(
	struct mux_data *data) {
	panda_mux_unref(data->mux);
	free(data);
}

void mux_fail_algn(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	PandaMux mux) {
	if (mux->noalgn == NULL)
		return;
	pthread_rwlock_rdlock(&mux->noalgn_rwlock);
	/* We repeat this check to solve a potential race condition. Locking the mutex is expensive, so if the handler is not set, we bail out. However, it is possible that someone is concurrently setting the handler via panda_mux_set_fail_alignment, so we have to check again once we have an exclusive lock. */
	if (mux->noalgn != NULL) {
		mux->noalgn(assembler, id, forward, forward_length, reverse, reverse_length, mux->noalgn_data);
	}
	pthread_rwlock_unlock(&mux->noalgn_rwlock);
}

PandaAssembler panda_mux_create_assembler(
	PandaMux mux) {
	return panda_mux_create_assembler_kmer(mux, PANDA_DEFAULT_NUM_KMERS);
}

PandaAssembler panda_mux_create_assembler_kmer(
	PandaMux mux,
	size_t num_kmers) {
	PandaAssembler assembler;
	struct mux_data *data = malloc(sizeof(struct mux_data));
	data->mux = panda_mux_ref(mux);
	assembler = panda_assembler_new_kmer((PandaNextSeq) mux_next, data, (PandaDestroy) mux_free, mux->logger, num_kmers);
	if (assembler != NULL) {
		char buffer[MAX_LEN];
		panda_assembler_set_fail_alignment(assembler, (PandaFailAlign) mux_fail_algn, mux, NULL);
		sprintf(buffer, "%p:%zd", (void *) mux, child_count(mux));
		panda_assembler_set_name(assembler, buffer);
	}
	return assembler;
}

size_t panda_mux_get_child_count(
	PandaMux mux) {
	return mux->child_count;
}

PandaLogProxy panda_mux_get_loggger(
	PandaMux mux) {
	return mux->logger;
}

void panda_mux_set_fail_alignment(
	PandaMux mux,
	PandaFailAlign handler,
	void *handler_data,
	PandaDestroy handler_destroy) {
	pthread_rwlock_wrlock(&mux->noalgn_rwlock);
	DESTROY_MEMBER(mux, noalgn);
	mux->noalgn_data = handler_data;
	mux->noalgn_destroy = handler_destroy;
	mux->noalgn = handler;
	pthread_rwlock_unlock(&mux->noalgn_rwlock);
}

#endif
