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
#include<assert.h>
#include<limits.h>
#include<math.h>
#include<stdbool.h>
#include<stddef.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq.h"
#include "assembler.h"
#include "parser.h"
#include "prob.h"
#include "plugin.h"
#include "table.h"
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#define KMER_LEN 8
typedef unsigned int bitstype;
#define FOR_BITS_IN_LIST(bits,index) for (index = 0; index < bits##_size; index++) if ((bits)[index / sizeof(bitstype) / 8] & (1 << (index % (8 * sizeof(bitstype)))))
#define BIT_LIST_SET(bits,index) (bits)[(index) / sizeof(bitstype) / 8] |= (1 << ((index) % (8 * sizeof(bitstype))));
#define BIT_LIST_GET(bits,index) ((bits)[(index) / sizeof(bitstype) / 8] & (1 << ((index) % (8 * sizeof(bitstype)))))
#define BITS_INIT(bits,size) bitstype bits[(size) / 8 + 1]; size_t bits##_size = (size); memset(&bits, 0, (size) / 8 + 1)
#define ALL_BITS_IF_NONE(bits,size) do { int _all = 0; int _bitctr; for (_bitctr = 0; _bitctr < ((size) / 8 + 1); _bitctr++) { _all |= (bits)[_bitctr]; } if (_all != 0) { memset(&(bits), 0xFF, (size) / 8 + 1); }} while (0)

typedef struct {
	unsigned int kmer;
	ssize_t posn;
	ssize_t bad;
} kmer_it;
#define _FOREACH_KMER(iterator,sequence,start,check,step) for ((iterator).posn = (start), (iterator).bad = KMER_LEN; (iterator).posn check; (iterator).posn step, (iterator).kmer = (((iterator).kmer << 2) | ((sequence)[(iterator).posn].nt == PANDA_NT_T ? 3 : (sequence)[(iterator).posn].nt == PANDA_NT_G ? 2 : (sequence)[(iterator).posn].nt == PANDA_NT_C ? 1 : 0)) & ((1 << (2 * KMER_LEN)) - 1)) if (PANDA_NT_IS_N((sequence)[(iterator).posn].nt)) { (iterator).bad = KMER_LEN; } else if ((iterator).bad > 0) { (iterator).bad--; } else
#define FOREACH_KMER(iterator,sequence) _FOREACH_KMER(iterator,sequence, 0, < sequence ## _length, ++)
#define FOREACH_KMER_REVERSE(iterator,sequence) _FOREACH_KMER(iterator,sequence, sequence ## _length - 1, >= 0, --)
#define KMER(kmerit) ((kmerit).kmer)
#define KMER_POSITION(kmerit) ((kmerit).posn)
#define NUM_KMERS 2
#define KMERSEEN_SIZE (sizeof(seqindex) * NUM_KMERS * (1 << (2 * KMER_LEN)))
#define VEEZ(x) ((x) < 0 ? 0 : (x))
#define WEDGEZ(x) ((x) > 0 ? 0 : (x))
PandaAssembler panda_assembler_new(PandaNextSeq next, void *next_data, PandaDestroy next_destroy, PandaLogger logger, void *logger_data, PandaDestroy logger_destroy)
{
	PandaAssembler assembler = malloc(sizeof(struct panda_assembler));
	if (assembler == NULL)
		return NULL;
	assembler->refcnt = 1;
	assembler->next = next;
	assembler->next_data = next_data;
	assembler->next_destroy = next_destroy;
	assembler->logger = logger;
	assembler->logger_data = logger_data;
	assembler->logger_destroy = logger_destroy;
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
	assembler->count = 0;
	assembler->no_n = false;
	assembler->kmerseen = malloc(KMERSEEN_SIZE);
#ifdef HAVE_PTHREAD
	pthread_mutex_init(&assembler->mutex, NULL);
#endif
	memset(assembler->kmerseen, 0, KMERSEEN_SIZE);
	panda_assembler_set_error_estimation(assembler, 0.36);
	panda_assembler_set_threshold(assembler, 0.6);
	panda_assembler_set_minimum_overlap(assembler, 1);
	return assembler;
}

int panda_assembler_get_minimum_overlap(PandaAssembler assembler) {
	return assembler->minoverlap;
}
void panda_assembler_set_minimum_overlap(PandaAssembler assembler, int overlap) {
	if (overlap > 1 && overlap < PANDA_MAX_LEN) {
		assembler->minoverlap = overlap;
	}
}

double panda_assembler_get_threshold(PandaAssembler assembler) {
	return exp(assembler->threshold);
}

void panda_assembler_set_threshold(PandaAssembler assembler, double threshold) {
	if (threshold > 0 && threshold < 1) {
		assembler->threshold = log(threshold);
	}
}

double panda_assembler_get_error_estimation(PandaAssembler assembler) {
	return assembler->q;
}
void panda_assembler_set_error_estimation(PandaAssembler assembler, double q) {
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

long panda_assembler_get_no_forward_primer_count(PandaAssembler assembler) {
	return assembler->nofpcount;
}
long panda_assembler_get_no_reverse_primer_count(PandaAssembler assembler) {
	return assembler->norpcount;
}
long panda_assembler_get_ok_count(PandaAssembler assembler) {
	return assembler->okcount;
}
long panda_assembler_get_low_quality_count(PandaAssembler assembler) {
	return assembler->lowqcount;
}
long panda_assembler_get_degenerate_count(PandaAssembler assembler) {
	return assembler->degencount;
}
long panda_assembler_get_failed_alignment_count(PandaAssembler assembler) {
	return assembler->noalgncount;
}
long panda_assembler_get_count(PandaAssembler assembler) {
	return assembler->count;
}
bool panda_assembler_get_disallow_degenerates(PandaAssembler assembler) {
	return assembler->no_n;
}
void panda_assembler_set_disallow_degenerates(PandaAssembler assembler, bool allow) {
	assembler->no_n = allow;
}

PandaAssembler panda_assembler_ref(PandaAssembler assembler) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&assembler->mutex);
#endif
	assembler->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&assembler->mutex);
#endif
	return assembler;
}

void panda_assembler_unref(PandaAssembler assembler) {
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
		free(assembler);
	}
}

#define CIRC(index, len) (((index) + (len)) % (len))
static double qualscore(panda_nt nt1, char qual1, panda_nt nt2, char qual2)
{
	if (PANDA_NT_IS_N(nt1)) {
		if (PANDA_NT_IS_N(nt2)) {
			return qual_nn;
		}
		return qual_nmatch[qual2];
	}
	if (PANDA_NT_IS_N(nt2)) {
		return qual_nmatch[qual1];
	}
	return (nt1 ==
		nt2 ? qual_match :
		qual_mismatch)[qual1][qual2];
}

/* Find the offset of a primer in a sequence and return the offset (i.e., the start of the useful sequence. */
static size_t computeoffset(PandaAssembler assembler, panda_qual *seq, size_t seq_length, panda_nt *primer, size_t primerlen)
{
	/* Circular buffer of probabilities of primer alignment indexed by the offset. */
	double probabilities[primerlen];
	double bestpr = primerlen * assembler->threshold;
	size_t bestindex = 0;
	size_t index;

	if (primerlen > seq_length) {
		return 0;
	}

	for (index = 0; index < primerlen; index++) {
		probabilities[index] = -INFINITY;
	}

	for (index = 0; index < seq_length; index++) {
		ssize_t x;
		for (x = (ssize_t) (primerlen > index ? index : primerlen - 1);
		     x >= 0; x--) {
			probabilities[CIRC(index - x, primerlen)] +=
			    qualscore(seq[index].nt, seq[index].qual,
				      primer[x], (char)PHREDMAX);
		}
		/* The last bucket in the buffer holds the probability of a complete alignment. If it so better than we have seen previously, store it. */
		if (probabilities[CIRC(index, primerlen)] > bestpr) {
			bestpr = probabilities[CIRC(index, primerlen)];
			bestindex = index + 1;
		}
		probabilities[CIRC(index, primerlen)] = 0;
	}
	return bestindex;
}

/* Try to align forward and reverse reads and return the quality of the aligned sequence and the sequence itself. */
static bool
align(PandaAssembler assembler, panda_result_seq* result, int maxresult)
{
	ssize_t i, j;
	ssize_t df, dr;
	/* For determining overlap. */
	size_t maxoverlap =
	    result->forward_length < result->reverse_length ? result->forward_length : result->reverse_length;
	double bestprobability = qual_nn * (result->forward_length + result->reverse_length);
	int bestoverlap = -1;
	size_t overlap;
	size_t counter;
	kmer_it it;

	/* For computing new sequence. */
	double fquality = 0;
	double oquality = 0;
	double rquality = 0;
	ssize_t len;
	BITS_INIT(posn, maxoverlap - assembler->minoverlap + 1);

	if (result->forward_length >= (1 << (8 * sizeof(seqindex)))) {
		LOG(assembler, PANDA_CODE_INSUFFICIENT_KMER_TABLE, &result->name);
		return false;
	}

	/* Scan forward sequence building k-mers and appending the position to kmerseen[k] */
	FOREACH_KMER(it, result->forward) {
#ifdef DEBUG
		LOG(assembler, PANDA_CODE_FORWARD_KMER, &result->name, KMER(it), KMER_POSITION(it));
#endif
		for (j = 0;
		     j < NUM_KMERS && assembler->kmerseen[(KMER(it) << 1) + j] != 0; j++) ;
		if (j == NUM_KMERS) {
			/* If we run out of storage, we lose k-mers. */
			LOG(assembler, PANDA_CODE_LOST_KMER, &result->name, KMER(it), KMER_POSITION(it));
		} else {
			assembler->kmerseen[(KMER(it) << 1) + j] = KMER_POSITION(it);
		}
	}

	/* Scan reverse sequence building k-mers. For each position in the forward sequence for this kmer (i.e., kmerseen[k]), flag that we should check the corresponding overlap. */
	FOREACH_KMER_REVERSE(it, result->reverse) {
#ifdef DEBUG
		LOG(assembler, PANDA_CODE_REVERSE_KMER, &result->name, KMER(it), KMER_POSITION(it));
#endif
		for (j = 0;
		     j < NUM_KMERS && assembler->kmerseen[(KMER(it) << 1) + j] != (seqindex)0; j++) {
			int index =
			    result->forward_length + result->reverse_length -
			    KMER_POSITION(it) - assembler->kmerseen[(KMER(it) << 1) + j] -
			    assembler->minoverlap - 1;

			if (index >= 0) {
				BIT_LIST_SET(posn, index);
			}
		}
	}

	/* Reset kmerseen */
	FOREACH_KMER(it, result->forward) {
		for (j = 0; j < NUM_KMERS; j++)
			assembler->kmerseen[(KMER(it) << 1) + j] = 0;
	}

	ALL_BITS_IF_NONE(posn, maxoverlap - assembler->minoverlap + 1);

	/* Compute the quality of the overlapping region for the various overlaps and pick the best one. */
	FOR_BITS_IN_LIST(posn, counter) {
		size_t matches = 0;
		size_t mismatches = 0;
		size_t unknowns = 0;
		double probability;
		overlap = counter + assembler->minoverlap;

		for (i = 0; i < overlap; i++) {
			int findex = result->forward_length + i - overlap;
			int rindex = result->reverse_length - i - 1;
			panda_nt f = result->forward[findex].nt;
			panda_nt r = result->reverse[rindex].nt;
			if (PANDA_NT_IS_N(f) || PANDA_NT_IS_N(r)) {
				unknowns++;
			} else if ((f & r) != 0) {
				matches++;
			} else {
				mismatches++;
			}
		}

		probability =
		    (qual_nn *
		     (result->forward_length + result->reverse_length - 2 * overlap +
		      unknowns) + matches * assembler->pmatch + mismatches * assembler->pmismatch);

#ifdef DEBUG
		LOG(assembler, PANDA_CODE_OVERLAP_POSSIBILITY, overlap, matches, mismatches, unknowns, probability);
#endif
		if (probability > bestprobability) {
			bestprobability = probability;
			bestoverlap = overlap;
		}
	}

	LOG(assembler, PANDA_CODE_BEST_OVERLAP, &result->name, bestoverlap);

	if (bestoverlap == -1) {
		return false;
	}

	/* Compute the correct alignment and the quality score of the entire sequence. */
	len =
	    result->forward_length - (ssize_t) result->forward_offset - bestoverlap + result->reverse_length -
	    (ssize_t) result->reverse_offset + 1;
	if (len <= 0) {
		LOG(assembler, PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH, &result->name);
		return false;
	}
	if (len > maxresult) {
		LOG(assembler, PANDA_CODE_SEQUENCE_TOO_LONG, &result->name);
		return false;
	}
	result->sequence_length = len - 1;
	result->degenerates = 0;

	df = (ssize_t) result->forward_length - (ssize_t) result->forward_offset - bestoverlap;
	dr = (ssize_t) result->reverse_length - (ssize_t) result->reverse_offset - bestoverlap;
	/* Copy the unpaired forward sequence. */
#ifdef DEBUG
	LOG(assembler, PANDA_CODE_RECONSTRUCTION_PARAM, (int)bestoverlap, (int)df, (int)dr);
#endif
	for (i = 0; i < VEEZ(df); i++) {
		int findex = i + result->forward_offset;
		panda_nt fbits = result->forward[findex].nt;
		double q = qual_score[result->forward[findex].qual];
		result->sequence[i].nt = fbits;
		result->sequence[i].p = q;
		if (PANDA_NT_IS_DEGN(fbits)) {
			result->degenerates++;
		}
		fquality += q;
#ifdef DEBUG
		LOG(assembler, PANDA_CODE_BUILD_FORWARD, &result->name, i, findex, &result->sequence[i]);
#endif
	}

	/* Mask out the B-cliff at the end of sequences */
	for (i = result->forward_length - 1;
	     i > 0 && result->forward[i].qual == (char)2; i--) {
		result->forward[i].qual = '\0';
	}
	for (i = result->reverse_length - 1;
	     i > 0 && result->reverse[i].qual == (char)2; i--) {
		result->reverse[i].qual = '\0';
	}
	/* Copy the paired sequence adjusting the probabilities based on the quality information from both sequences. */
	for (i = 0; i < bestoverlap + WEDGEZ(df) + WEDGEZ(dr); i++) {
		int index = VEEZ(df) + i;
		int findex = result->forward_offset + VEEZ(df) + i;
		int rindex = result->reverse_length - i - 1 + WEDGEZ(df);
		bool ismatch =
		    (result->reverse[rindex].nt & result->forward[findex].nt) != '\0';
		double fpr;
		double rpr;
		double q;
		char nt;

		if (index < 0 || findex < 0 || rindex < 0
		    || findex >= result->forward_length || rindex >= result->reverse_length)
			continue;

		fpr = result->forward[findex].qual == '\0' ? qual_nn :
		    qual_score[result->forward[findex].qual];
		rpr = result->reverse[rindex].qual == '\0' ? qual_nn :
		    qual_score[result->reverse[rindex].qual];

		if (!ismatch) {
			LOG(assembler, PANDA_CODE_MISMATCHED_BASE, findex, rindex, &result->forward[findex],
				&result->reverse[rindex]);
		}

		if (result->forward[findex].qual == '\0'
		    && result->reverse[rindex].qual == '\0') {
			q = qual_nn;
		} else if (result->forward[findex].qual == '\0') {
			q = ismatch ? rpr : qual_nn;
		} else if (result->reverse[rindex].qual == '\0') {
			q = ismatch ? fpr : qual_nn;
		} else {
			q = (ismatch ? qual_match :
			     qual_mismatch)[result->forward[findex].qual][result->reverse[rindex].qual];
		}

		if (ismatch) {
			nt = (result->reverse[rindex].nt & result->forward[findex].nt);
		} else {
			if (result->forward[rindex].qual < result->reverse[rindex].qual) {
				nt = result->reverse[rindex].nt;
			} else {
				nt = result->forward[findex].nt;
			}
		}
		result->sequence[index].nt = nt;
		result->sequence[index].p = q;
		if (PANDA_NT_IS_DEGN(nt)) {
			result->degenerates++;
		}
		oquality += q;
#ifdef DEBUG
		LOG(assembler, PANDA_CODE_BUILD_OVERLAP, &result->name, index, findex, rindex, &result->sequence[index], &result->forward[findex], &result->reverse[rindex]);
#endif
	}

	/* Copy the unpaired reverse sequence. */
	for (i = 0; i < VEEZ(dr); i++) {
		int index = df + bestoverlap + i;
		int rindex = result->reverse_length - bestoverlap - i - 1;
		panda_nt rbits = result->reverse[rindex].nt;
		double q = qual_score[result->reverse[rindex].qual];
		rquality += q;
		result->sequence[index].nt = rbits;
		result->sequence[index].p = q;
		if (PANDA_NT_IS_DEGN(rbits)) {
			result->degenerates++;
		}
#ifdef DEBUG
		LOG(assembler, PANDA_CODE_BUILD_REVERSE, &result->name, index, rindex, &result->sequence[index]);
#endif
	}
	result->quality = (fquality + rquality + oquality) / len;

	return true;
}

bool assemble_seq(PandaAssembler assembler) {
	assembler->count++;
	if (!module_precheckseq(assembler, &assembler->result.name, assembler->result.forward, assembler->result.forward_length, assembler->result.reverse, assembler->result.reverse_length)) {
		return false;
	}
	if (assembler->forward_primer_length > 0) {
		assembler->result.forward_offset = computeoffset(assembler, assembler->result.forward, assembler->result.forward_length, assembler->forward_primer, assembler->forward_primer_length);
		if (assembler->result.forward_offset == 0) {
			LOG(assembler, PANDA_CODE_NO_FORWARD_PRIMER, &assembler->result.name);
			assembler->nofpcount++;
			return false;
		}
		assembler->result.forward_offset--;
	} else {
		assembler->result.forward_offset = assembler->forward_trim;
	}
	if (assembler->reverse_primer_length > 0) {
		assembler->result.reverse_offset = computeoffset(assembler, assembler->result.reverse, assembler->result.reverse_length, assembler->reverse_primer, assembler->reverse_primer_length);
		if (assembler->result.reverse_offset == 0) {
			LOG(assembler, PANDA_CODE_NO_REVERSE_PRIMER, &assembler->result.name);
			assembler->norpcount++;
			return false;
		}
		assembler->result.reverse_offset--;
	} else {
		assembler->result.reverse_offset = assembler->reverse_trim;
	}
	if (!align(assembler, &assembler->result, PANDA_MAX_LEN)) {
		return false;
	}
	if (assembler->result.quality < assembler->threshold) {
		assembler->lowqcount++;
		LOG(assembler, PANDA_CODE_LOW_QUALITY_REJECT, exp(assembler->result.quality), exp(assembler->threshold));
		return false;
	}
	if (module_checkseq(assembler, &assembler->result)) {
		assembler->okcount++;
		return true;
	}
}

const panda_result_seq *panda_assembler_next(PandaAssembler assembler) {
	if (assembler->next == NULL) {
		return NULL;
	}
	while(true) {
		if (!assembler->next(&assembler->result.name, &assembler->result.forward, &assembler->result.forward_length, &assembler->result.reverse, &assembler->result.reverse_length, assembler->next_data)) {
			return NULL;
		}
		assert(assembler->result.forward_length <= PANDA_MAX_LEN);
		assert(assembler->result.reverse_length <= PANDA_MAX_LEN);
		if (assemble_seq(assembler)) {
			return &assembler->result;
		}
	}
	return NULL;
}

const panda_result_seq *panda_assembler_assemble(PandaAssembler assembler, panda_seq_identifier *id, const panda_qual *forward, size_t forward_length, const panda_qual *reverse, size_t reverse_length) {
	size_t it;
	assert(forward_length <= PANDA_MAX_LEN);
	assert(reverse_length <= PANDA_MAX_LEN);
	assembler->result.name = *id;
	assembler->result.forward_length = forward_length;
	assembler->result.reverse_length = reverse_length;
	for (it = 0; it < forward_length; it++)
		assembler->result.forward[it] = forward[it];
	for (it = 0; it < reverse_length; it++)
		assembler->result.reverse[it] = reverse[it];
	return assemble_seq(assembler) ? &assembler->result : NULL;
}
