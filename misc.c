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
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "pandaseq.h"
#include "misc.h"

const char const *panda_version(
	void) {
	return PACKAGE_STRING;
}

int panda_api_version(
	void) {
	return PANDA_API;
}

size_t panda_max_len(
	void) {
	return MAX_LEN;
}

enum iter_type {
	ITER_QUAL,
	ITER_NT,
	ITER_RESULT
};

struct panda_iter {
	enum iter_type type;
	bool reverse;
	int k;
	kmer_it it;

	panda_qual *qual;
	size_t qual_length;

	panda_nt *nt;
	size_t nt_length;

	panda_result *result;
	size_t result_length;

	panda_kmer output;
};

void panda_iter_free(
	PandaIter iter) {
	free(iter);
}

PandaIter panda_iter_dup(
	PandaIter iter) {
	PandaIter new_iter = malloc(sizeof(struct panda_iter));
	memcpy(new_iter, iter, sizeof(struct panda_iter));
	return new_iter;
}

void panda_iter_reset(
	PandaIter iter) {
	if (iter->reverse) {
		switch (iter->type) {
		case ITER_QUAL:
			iter->it.posn = iter->qual_length;
			break;
		case ITER_NT:
			iter->it.posn = iter->nt_length;
			break;
		case ITER_RESULT:
			iter->it.posn = iter->result_length;
			break;
		}
	} else {
		iter->it.posn = -1;
	}
	iter->it.bad = iter->k;
}

int panda_iter_k(
	PandaIter iter) {
	return iter->k;
}

size_t panda_iter_bits(
	PandaIter iter) {
	return iter->k * 2;
}

#define RETURN_KMER return iter->output.kmer = iter->it.kmer, iter->output.posn = iter->it.posn, &iter->output
const panda_kmer *panda_iter_next(
	PandaIter iter) {
	switch (iter->type) {
	case ITER_QUAL:
		if (iter->reverse) {
			iter->it.posn--;
			_FOREACH_KMER(iter->it, iter->qual,.nt, iter->it.posn, iter->it.bad, >=0, --, iter->k) {
				RETURN_KMER;
			}
		} else {
			iter->it.posn++;
			_FOREACH_KMER(iter->it, iter->qual,.nt, iter->it.posn, iter->it.bad, <iter->qual_length, ++, iter->k) {
				RETURN_KMER;
			}
		}
		return NULL;
	case ITER_NT:
		if (iter->reverse) {
			iter->it.posn--;
			_FOREACH_KMER(iter->it, iter->nt,, iter->it.posn, iter->it.bad, >=0, --, iter->k) {
				RETURN_KMER;
			}
		} else {
			iter->it.posn++;
			_FOREACH_KMER(iter->it, iter->nt,, iter->it.posn, iter->it.bad, <iter->nt_length, ++, iter->k) {
				RETURN_KMER;
			}
		}
		return NULL;
	case ITER_RESULT:
		if (iter->reverse) {
			iter->it.posn--;
			_FOREACH_KMER(iter->it, iter->result,.nt, iter->it.posn, iter->it.bad, >=0, --, iter->k) {
				RETURN_KMER;
			}
		} else {
			iter->it.posn++;
			_FOREACH_KMER(iter->it, iter->result,.nt, iter->it.posn, iter->it.bad, <iter->result_length, ++, iter->k) {
				RETURN_KMER;
			}
		}
		return NULL;
	}
}

static PandaIter iter_new(
	enum iter_type type,
	bool reverse,
	int k) {
	PandaIter iter = malloc(sizeof(struct panda_iter));
	iter->type = type;
	iter->reverse = reverse;
	if (k < 1) {
		iter->k = KMER_LEN;
	} else {
		iter->k = (k < sizeof(size_t) * 4) ? k : (sizeof(size_t) * 4);
	}
	return iter;
}

PandaIter panda_iterate_qual(
	panda_qual *seq,
	size_t seq_length,
	bool reverse,
	int k) {
	PandaIter iter = iter_new(ITER_QUAL, reverse, k);
	iter->qual = seq;
	iter->qual_length = seq_length;
	panda_iter_reset(iter);
	return iter;
}

PandaIter panda_iterate_nt(
	panda_nt *seq,
	size_t seq_length,
	bool reverse,
	int k) {
	PandaIter iter = iter_new(ITER_NT, reverse, k);
	iter->nt = seq;
	iter->nt_length = seq_length;
	panda_iter_reset(iter);
	return iter;
}

PandaIter panda_iterate_result(
	panda_result *seq,
	size_t seq_length,
	bool reverse,
	int k) {
	PandaIter iter = iter_new(ITER_RESULT, reverse, k);
	iter->result = seq;
	iter->result_length = seq_length;
	panda_iter_reset(iter);
	return iter;
}
