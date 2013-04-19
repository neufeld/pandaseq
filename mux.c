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
#if HAVE_PTHREAD
#        include <stdlib.h>
#        include <string.h>
#        include "pandaseq.h"
#        include "assembler.h"
#        include "buffer.h"

struct panda_mux {
	pthread_mutex_t mutex;
	pthread_mutex_t next_mutex;
	pthread_mutex_t logger_mutex;
	pthread_mutex_t noalgn_mutex;
	 MANAGED_MEMBER(
		PandaNextSeq,
		next);
	 MANAGED_MEMBER(
		PandaLogger,
		logger);
	 MANAGED_MEMBER(
		PandaFailAlign,
		noalgn);
	volatile size_t refcnt;
};

PandaMux panda_mux_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogger logger,
	void *logger_data,
	PandaDestroy logger_destroy) {
	PandaMux mux = malloc(sizeof(struct panda_mux));
	mux->refcnt = 1;
	mux->next = next;
	mux->next_data = next_data;
	mux->next_destroy = next_destroy;
	mux->logger = logger;
	mux->logger_data = logger_data;
	mux->logger_destroy = logger_destroy;
	mux->noalgn = NULL;
	mux->noalgn_data = NULL;
	mux->noalgn_destroy = NULL;
	pthread_mutex_init(&mux->mutex, NULL);
	pthread_mutex_init(&mux->next_mutex, NULL);
	pthread_mutex_init(&mux->logger_mutex, NULL);
	pthread_mutex_init(&mux->noalgn_mutex, NULL);
	return mux;
}

PandaMux panda_mux_new_fastq_reader(
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
	next = panda_create_fastq_reader(forward, forward_data, forward_destroy, reverse, reverse_data, reverse_destroy, (PandaLogger) logger, logger_data, qualmin, policy, &user_data, &destroy);
	return panda_mux_new(next, user_data, destroy, logger, logger_data, logger_destroy);
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
	int count;
	pthread_mutex_lock(&mux->mutex);
	count = --(mux->refcnt);
	pthread_mutex_unlock(&mux->mutex);
	if (count == 0) {
		pthread_mutex_destroy(&mux->mutex);

		pthread_mutex_lock(&mux->next_mutex);
		DESTROY_MEMBER(mux, next);
		pthread_mutex_unlock(&mux->next_mutex);
		pthread_mutex_destroy(&mux->next_mutex);

		pthread_mutex_lock(&mux->logger_mutex);
		DESTROY_MEMBER(mux, logger);
		pthread_mutex_unlock(&mux->logger_mutex);
		pthread_mutex_destroy(&mux->logger_mutex);

		pthread_mutex_lock(&mux->noalgn_mutex);
		DESTROY_MEMBER(mux, noalgn);
		pthread_mutex_unlock(&mux->noalgn_mutex);
		pthread_mutex_destroy(&mux->noalgn_mutex);
		free(mux);
	}
}

static bool mux_logger(
	PandaCode code,
	panda_seq_identifier *id,
	char *msg,
	PandaMux mux) {
	bool ret;
	pthread_mutex_lock(&mux->logger_mutex);
	ret = (mux->logger) (code, id, msg, mux->logger_data);
	pthread_mutex_unlock(&mux->logger_mutex);
	return ret;
}

static bool mux_next(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	PandaMux mux) {
	bool result;
	panda_qual *common_forward;
	panda_qual *common_reverse;
	pthread_mutex_lock(&mux->next_mutex);

	result = mux->next(id, &common_forward, forward_length, &common_reverse, reverse_length, mux->next_data);
	if (common_forward == NULL) {
		*forward = NULL;
		*forward_length = 0;
	} else {
		*forward = forward_buffer();
		memcpy(*forward, common_forward, sizeof(panda_qual) * *forward_length);
	}
	if (common_reverse == NULL) {
		*reverse = NULL;
		*reverse_length = 0;
	} else {
		*reverse = reverse_buffer();
		memcpy(*reverse, common_reverse, sizeof(panda_qual) * *reverse_length);
	}
	pthread_mutex_unlock(&mux->next_mutex);
	return result;
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
	pthread_mutex_lock(&mux->noalgn_mutex);
	/* We repeat this check to solve a potential race condition. Locking the mutex is expensive, so if the handler is not set, we bail out. However, it is possible that someone is concurrently setting the handler via panda_mux_set_fail_alignment, so we have to check again once we have an exclusive lock. */
	if (mux->noalgn != NULL) {
		mux->noalgn(assembler, id, forward, forward_length, reverse, reverse_length, mux->noalgn_data);
	}
	pthread_mutex_unlock(&mux->noalgn_mutex);
}

PandaAssembler panda_mux_create_assembler(
	PandaMux mux) {
	return panda_mux_create_assembler_kmer(mux, PANDA_DEFAULT_NUM_KMERS);
}

PandaAssembler panda_mux_create_assembler_kmer(
	PandaMux mux,
	size_t num_kmers) {
	PandaAssembler assembler;
	pthread_mutex_lock(&mux->mutex);
	mux->refcnt += 2;
	pthread_mutex_unlock(&mux->mutex);
	assembler = panda_assembler_new_kmer((PandaNextSeq) mux_next, mux, (PandaDestroy) panda_mux_unref, (PandaLogger) mux_logger, mux, (PandaDestroy) panda_mux_unref, num_kmers);
	if (assembler != NULL) {
		panda_assembler_set_fail_alignment(assembler, (PandaFailAlign)
			mux_fail_algn, mux, NULL);
	}
	return assembler;
}

void panda_mux_set_fail_alignment(
	PandaMux mux,
	PandaFailAlign handler,
	void *handler_data,
	PandaDestroy handler_destroy) {
	pthread_mutex_lock(&mux->noalgn_mutex);
	DESTROY_MEMBER(mux, noalgn);
	mux->noalgn_data = handler_data;
	mux->noalgn_destroy = handler_destroy;
	mux->noalgn = handler;
	pthread_mutex_unlock(&mux->noalgn_mutex);
}
#endif
