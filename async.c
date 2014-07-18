/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
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
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "pandaseq.h"
#include "misc.h"

#ifdef HAVE_PTHREAD
struct seq_data {
	panda_seq_identifier id;
	panda_qual forward[MAX_LEN];
	size_t forward_length;
	panda_qual reverse[MAX_LEN];
	size_t reverse_length;
	struct seq_data *next;
};

struct async_data {
	MANAGED_MEMBER(
		PandaNextSeq,
		next);
	pthread_t reader;
	pthread_cond_t is_ready;
	pthread_cond_t has_free;
	volatile bool done;

	pthread_key_t in_flight;
	pthread_mutex_t free_mutex;
	pthread_mutex_t ready_mutex;
	struct seq_data *ready;
	struct seq_data *free;
};

static bool async_next_seq(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	struct async_data *data) {
	struct seq_data *seq = NULL;

	*forward = NULL;
	*reverse = NULL;
	*forward_length = 0;
	*reverse_length = 0;

	seq = pthread_getspecific(data->in_flight);
	if (seq != NULL) {
		pthread_mutex_lock(&data->free_mutex);
		seq->next = data->free;
		data->free = seq;
		pthread_setspecific(data->in_flight, NULL);
		pthread_cond_signal(&data->has_free);
		pthread_mutex_unlock(&data->free_mutex);
	}

	pthread_mutex_lock(&data->ready_mutex);
	do {
		seq = data->ready;
		if (seq != NULL) {
			data->ready = seq->next;
			pthread_mutex_unlock(&data->ready_mutex);

			pthread_setspecific(data->in_flight, seq);
			*forward = seq->forward;
			*forward_length = seq->forward_length;
			*reverse = seq->reverse;
			*reverse_length = seq->reverse_length;
			*id = seq->id;
			seq->next = NULL;
			return true;
		}
		if (data->done) {
			pthread_mutex_unlock(&data->ready_mutex);
			return false;
		}
	} while (pthread_cond_wait(&data->is_ready, &data->ready_mutex) == 0);
	pthread_mutex_unlock(&data->ready_mutex);
	return false;
}

static void *async_thread(
	struct async_data *data) {
	while (true) {
		struct seq_data *seq;

		pthread_mutex_lock(&data->free_mutex);
		while ((seq = data->free) == NULL) {
			if (data->done || pthread_cond_wait(&data->has_free, &data->free_mutex) != 0 || data->done) {
				data->done = true;
				pthread_mutex_unlock(&data->free_mutex);
				pthread_mutex_lock(&data->ready_mutex);
				pthread_cond_broadcast(&data->is_ready);
				pthread_mutex_unlock(&data->ready_mutex);
				pthread_exit(NULL);
			}
		}
		data->free = NULL;
		pthread_mutex_unlock(&data->free_mutex);
		while (seq != NULL) {
			const panda_qual *forward;
			const panda_qual *reverse;
			if (data->next(&seq->id, &forward, &seq->forward_length, &reverse, &seq->reverse_length, data->next_data)) {
				struct seq_data *next;
				memcpy(seq->forward, forward, seq->forward_length * sizeof(panda_qual));
				memcpy(seq->reverse, reverse, seq->reverse_length * sizeof(panda_qual));
				next = seq->next;
				pthread_mutex_lock(&data->ready_mutex);
				seq->next = data->ready;
				data->ready = seq;
				pthread_cond_broadcast(&data->is_ready);
				pthread_mutex_unlock(&data->ready_mutex);
				seq = next;
			} else {
				struct seq_data *next;
				data->done = true;
				for (; seq != NULL; seq = next) {
					next = seq->next;
					free(seq);
				}
				pthread_mutex_lock(&data->ready_mutex);
				pthread_cond_broadcast(&data->is_ready);
				pthread_mutex_unlock(&data->ready_mutex);
				seq = next;
				pthread_exit(NULL);
			}
		}
	}
}

static void async_destroy(
	struct async_data *data) {
	struct seq_data *seq;
	struct seq_data *temp;
	data->done = true;

	pthread_mutex_lock(&data->free_mutex);
	pthread_cond_broadcast(&data->has_free);
	pthread_mutex_unlock(&data->free_mutex);
	pthread_join(data->reader, NULL);
	pthread_cond_destroy(&data->is_ready);
	pthread_cond_destroy(&data->has_free);
	pthread_mutex_destroy(&data->free_mutex);
	pthread_mutex_destroy(&data->ready_mutex);
	pthread_key_delete(data->in_flight);

	DESTROY_MEMBER(data, next);
	for (seq = data->ready; seq != NULL; seq = (temp = seq->next, free(seq), temp)) ;
	for (seq = data->free; seq != NULL; seq = (temp = seq->next, free(seq), temp)) ;
	free(data);
}

PandaNextSeq panda_create_async_reader(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	size_t length,
	void **user_data,
	PandaDestroy *destroy) {
	struct async_data *data;

	if (length < 2) {
		*user_data = next_data;
		*destroy = next_destroy;
		return next;
	}

	data = malloc(sizeof(struct async_data));
	data->done = false;
	data->next = next;
	data->next_data = next_data;
	data->next_destroy = next_destroy;

	pthread_cond_init(&data->is_ready, NULL);
	pthread_cond_init(&data->has_free, NULL);
	pthread_mutex_init(&data->free_mutex, NULL);
	pthread_mutex_init(&data->ready_mutex, NULL);
	pthread_key_create(&data->in_flight, free);

	data->ready = NULL;
	data->free = NULL;
	length *= 4;
	for (; length > 0; length--) {
		struct seq_data *seq = malloc(sizeof(struct seq_data));
		seq->next = data->free;
		data->free = seq;
	}

	pthread_create(&data->reader, NULL, (void *(*)(void *)) &async_thread, data);

	*user_data = data;
	*destroy = (PandaDestroy) async_destroy;
	return (PandaNextSeq) async_next_seq;
}
#else
PandaNextSeq panda_create_async_reader(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	size_t length,
	void **user_data,
	PandaDestroy *destroy) {

	(void) length;

	*user_data = next_data;
	*destroy = next_destroy;
	return next;
}
#endif
