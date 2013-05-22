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
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

struct panda_log_proxy {
	size_t refcnt;
	 MANAGED_MEMBER(
		PandaLogger,
		log);
#ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#endif
};

PandaLogProxy panda_log_proxy_new(
	PandaLogger log,
	void *log_data,
	PandaDestroy log_destroy) {
	PandaLogProxy proxy = malloc(sizeof(struct panda_log_proxy));
	proxy->refcnt = 1;
	proxy->log = log;
	proxy->log_data = log_data;
	proxy->log_destroy = log_destroy;
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&proxy->mutex, NULL);
#endif
	return proxy;
}

PandaLogProxy panda_log_proxy_new_stderr(
	) {
	return panda_log_proxy_new((PandaLogger) panda_logger_file, stderr, NULL);
}

PandaLogProxy panda_log_proxy_ref(
	PandaLogProxy proxy) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	proxy->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
	return proxy;
}

void panda_log_proxy_unref(
	PandaLogProxy proxy) {
	size_t count;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	count = --(proxy->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
	if (count == 0) {
		DESTROY_MEMBER(proxy, log);
		pthread_mutex_destroy(&proxy->mutex);
		free(proxy);
	}
}

bool panda_log_proxy_write(
	PandaLogProxy proxy,
	PandaCode code,
	panda_seq_identifier *id,
	const char *message) {
	bool result;
	if (proxy->log == NULL)
		return true;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&proxy->mutex);
#endif
	result = proxy->log(code, id, message, proxy->log_data);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&proxy->mutex);
#endif
	return result;
}
