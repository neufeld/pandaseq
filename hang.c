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
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "config.h"
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

struct hang_data {
	MANAGED_MEMBER(
		PandaNextSeq,
		next);
	PandaLogProxy logger;
	panda_nt forward[MAX_LEN];
	size_t forward_length;
	panda_nt reverse[MAX_LEN];
	size_t reverse_length;
	bool skip;
	double threshold;
};

bool hang_next(
	panda_seq_identifier *id,
	panda_qual **forward,
	size_t *forward_length,
	panda_qual **reverse,
	size_t *reverse_length,
	void *user_data) {
	struct hang_data *data = (struct hang_data *) user_data;
	while (data->next(id, forward, forward_length, reverse, reverse_length, data->next_data)) {
		size_t offset;
		if (data->forward_length > 0) {
			offset = panda_compute_offset_qual(data->threshold, true, *forward, *forward_length, data->forward, data->forward_length);
			if (offset == 0) {
				panda_log_proxy_write(data->logger, PANDA_CODE_NO_FORWARD_PRIMER, NULL, id, "OVERHANGING REJECT");
				if (!data->skip)
					continue;
			} else {
				*forward_length -= offset - 1;
			}
		}
		if (data->reverse_length > 0) {
			offset = panda_compute_offset_qual(data->threshold, true, *reverse, *reverse_length, data->reverse, data->reverse_length);
			if (offset == 0) {
				panda_log_proxy_write(data->logger, PANDA_CODE_NO_REVERSE_PRIMER, NULL, id, "OVERHANGING REJECT");
				if (!data->skip)
					continue;
			} else {
				*reverse_length -= offset - 1;
			}
		}
		return true;
	}
	return false;
}

void hang_free(
	void *user_data) {
	struct hang_data *hang_data = (struct hang_data *) user_data;
	DESTROY_MEMBER(hang_data, next);
	panda_log_proxy_unref(hang_data->logger);
	free(hang_data);
}

PandaNextSeq panda_trim_overhangs(
	PandaNextSeq inner,
	void *inner_data,
	PandaDestroy inner_destroy,
	PandaLogProxy logger,
	panda_nt *forward,
	size_t forward_length,
	panda_nt *reverse,
	size_t reverse_length,
	bool skip,
	double threshold,
	void **next_data,
	PandaDestroy *next_destroy) {
	struct hang_data *data = malloc(sizeof(struct hang_data));
	size_t it;

	data->next = inner;
	data->next_data = inner_data;
	data->next_destroy = inner_destroy;
	data->skip = skip;
	data->threshold = threshold;
	for (it = 0; it < forward_length; it++)
		data->forward[forward_length - it - 1] = forward[it];
	for (it = 0; it < reverse_length; it++)
		data->reverse[reverse_length - it - 1] = reverse[it];
	data->forward_length = forward_length;
	data->reverse_length = reverse_length;
	data->logger = panda_log_proxy_ref(logger);

	*next_data = data;
	*next_destroy = hang_free;
	return hang_next;
}
