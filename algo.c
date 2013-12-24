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
#include "config.h"
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "algo.h"

double panda_algorithm_quality_compare(
	PandaAlgorithm algorithm,
	const panda_qual *a,
	const panda_qual *b) {
	return algorithm->clazz->match_probability(&algorithm->end, (a->nt & b->nt) != '\0', a->qual, b->qual);
}

void *panda_algorithm_data(
	PandaAlgorithm algo) {
	return &algo->end;
}

PandaAlgorithmClass panda_algorithm_class(
	PandaAlgorithm algo) {
	return algo->clazz;
}

bool panda_algorithm_is_a(
	PandaAlgorithm algo,
	PandaAlgorithmClass clazz) {
	return algo != NULL && algo->clazz == clazz;
}

PandaAlgorithm panda_algorithm_ref(
	PandaAlgorithm algo) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&algo->mutex);
#endif
	algo->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&algo->mutex);
#endif
	return algo;
}

void panda_algorithm_unref(
	PandaAlgorithm algo) {
	size_t count;
	if (algo == NULL)
		return;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&algo->mutex);
#endif
	count = --(algo->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&algo->mutex);
#endif
	if (count == 0) {
#ifdef HAVE_PTHREAD
		pthread_mutex_destroy(&algo->mutex);
#endif
		if (algo->clazz->data_destroy != NULL) {
			algo->clazz->data_destroy(&algo->end);
		}
		free(algo);
	}
}

PandaAlgorithm panda_algorithm_new(
	PandaAlgorithmClass clazz) {
	PandaAlgorithm instance = malloc(sizeof(struct panda_algorithm) + clazz->data_size);
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&instance->mutex, NULL);
#endif
	instance->refcnt = 1;
	instance->clazz = clazz;
	return instance;
}

const PandaAlgorithmClass panda_algorithms[] = {
	&panda_algorithm_simple_bayes_class,
	&panda_algorithm_pear_class,
};

const size_t panda_algorithms_length = sizeof(panda_algorithms) / sizeof(PandaAlgorithmClass);
