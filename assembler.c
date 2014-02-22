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
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "pandaseq.h"
#include "algo.h"
#include "assembler.h"
#include "buffer.h"
#include "misc.h"
#include "module.h"
#include "prob.h"
#include "table.h"

#define LOG(flag, code) do { if(panda_debug_flags & flag) panda_log_proxy_write(assembler->logger, (code), assembler, &assembler->result.name, NULL); } while(0)
#define LOGV(flag, code, fmt, ...) do { if(panda_debug_flags & flag) { snprintf(static_buffer(), BUFFER_SIZE, fmt, __VA_ARGS__); panda_log_proxy_write(assembler->logger, (code), assembler, &assembler->result.name, static_buffer()); }} while(0)

typedef unsigned int bitstype;
#define FOR_BITS_IN_LIST(bits,index) for (index = 0; index < bits##_size; index++) if ((bits)[index / sizeof(bitstype) / 8] & (1 << (index % (8 * sizeof(bitstype)))))
#define BIT_LIST_SET(bits,index) (bits)[(index) / sizeof(bitstype) / 8] |= (1 << ((index) % (8 * sizeof(bitstype))));
#define BIT_LIST_GET(bits,index) ((bits)[(index) / sizeof(bitstype) / 8] & (1 << ((index) % (8 * sizeof(bitstype)))))
#define BITS_INIT(bits,size) bitstype bits[(size) / 8 / sizeof(bitstype) + 1]; size_t bits##_size = (size); memset(&bits, 0, ((size) / 8 / sizeof(bitstype) + 1) * sizeof(bitstype))
#define ALL_BITS_IF_NONE(bits) do { bitstype _all = 0; size_t _bitctr; for (_bitctr = 0; _bitctr < ((bits##_size) / 8 / sizeof(bitstype) + 1); _bitctr++) { _all |= (bits)[_bitctr]; } if (_all == 0) { memset(&bits, 0xFF, (bits##_size / 8 / sizeof(bitstype) + 1) * sizeof(bitstype)); }} while (0)

#define VEEZ(x) ((x) < 0 ? 0 : (x))
#define WEDGEZ(x) ((x) > 0 ? 0 : (x))

/* Try to align forward and reverse reads and return the quality of the aligned sequence and the sequence itself. */
static bool align(
	PandaAssembler assembler,
	panda_result_seq *result) {
	ssize_t i, j;
	ssize_t df, dr;
	/* Cache all algorithm information. */
	double qual_nn = assembler->algo->clazz->prob_unpaired;
	void *algo_data = panda_algorithm_data(assembler->algo);
	PandaComputeOverlap overlap_probability = assembler->algo->clazz->overlap_probability;
	PandaComputeMatch match_probability = assembler->algo->clazz->match_probability;
	/* For determining overlap. */
	size_t maxoverlap = result->forward_length + result->reverse_length - assembler->minoverlap - result->forward_offset - result->reverse_offset - 1;
	double bestprobability = qual_nn * (result->forward_length + result->reverse_length);
	ssize_t bestoverlap = -1;
	size_t overlap;
	size_t counter;
	kmer_it it;
	size_t unmasked_forward_length;
	size_t unmasked_reverse_length;

	/* For computing new sequence. */
	double fquality = 0;
	double oquality = 0;
	double rquality = 0;
	ssize_t len;

	if (maxoverlap > assembler->maxoverlap) {
		maxoverlap = assembler->maxoverlap;
	}

	BITS_INIT(posn, assembler->minoverlap <= maxoverlap ? (maxoverlap - assembler->minoverlap + 1) : 1);

	if (result->forward_length >= (1 << (8 * sizeof(seqindex)))) {
		LOG(PANDA_DEBUG_BUILD, PANDA_CODE_INSUFFICIENT_KMER_TABLE);
		return false;
	}

	/* Scan forward sequence building k-mers and appending the position to kmerseen[k] */
	FOREACH_KMER(it, result->forward,.nt) {
		LOGV(PANDA_DEBUG_KMER, PANDA_CODE_FORWARD_KMER, "%zd@%zu", KMER(it), KMER_POSITION(it));
		for (j = 0; j < assembler->num_kmers && assembler->kmerseen[(KMER(it) << 1) + j] != 0; j++) ;
		if (j == assembler->num_kmers) {
			/* If we run out of storage, we lose k-mers. */
			LOGV(PANDA_DEBUG_BUILD, PANDA_CODE_LOST_KMER, "%zd@%zu", KMER(it), KMER_POSITION(it));
		} else {
			assembler->kmerseen[(KMER(it) * assembler->num_kmers) + j] = KMER_POSITION(it);
		}
	}

	/* Scan reverse sequence building k-mers. For each position in the forward sequence for this kmer (i.e., kmerseen[k]), flag that we should check the corresponding overlap. */
	FOREACH_KMER_REVERSE(it, result->reverse,.nt) {
		LOGV(PANDA_DEBUG_KMER, PANDA_CODE_REVERSE_KMER, "%zd@%zu", KMER(it), KMER_POSITION(it));
		for (j = 0; j < assembler->num_kmers && assembler->kmerseen[(KMER(it) * assembler->num_kmers) + j] != (seqindex) 0; j++) {
			int index = result->forward_length + result->reverse_length - KMER_POSITION(it) - assembler->kmerseen[(KMER(it) * assembler->num_kmers) + j] - assembler->minoverlap - 1;
			if (index >= 0) {
				BIT_LIST_SET(posn, index);
			}
		}
	}

	/* Reset kmerseen */
	FOREACH_KMER(it, result->forward,.nt) {
		for (j = 0; j < assembler->num_kmers; j++)
			assembler->kmerseen[(KMER(it) * assembler->num_kmers) + j] = 0;
	}

	ALL_BITS_IF_NONE(posn);

	result->overlaps_examined = 0;
	/* Compute the quality of the overlapping region for the various overlaps and pick the best one. */
	FOR_BITS_IN_LIST(posn, counter) {
		double probability;
		size_t overlap = counter + assembler->minoverlap;
		probability = overlap_probability(algo_data, result->forward, result->forward_length, result->reverse, result->reverse_length, overlap);

		LOGV(PANDA_DEBUG_RECON, PANDA_CODE_OVERLAP_POSSIBILITY, "overlap = %zd probability = %f", overlap, probability);
		if (probability > bestprobability && overlap >= assembler->minoverlap) {
			bestprobability = probability;
			bestoverlap = overlap;
		}
		result->overlaps_examined++;
	}

	if (result->overlaps_examined == maxoverlap - assembler->minoverlap + 1) {
		assembler->slowcount++;
	}

	LOGV(PANDA_DEBUG_BUILD, PANDA_CODE_BEST_OVERLAP, "%zd", bestoverlap);

	if (bestoverlap == -1) {
		return false;
	}

	/* Compute the correct alignment and the quality score of the entire sequence. */
	len = result->forward_length - (ssize_t) result->forward_offset - bestoverlap + result->reverse_length - (ssize_t) result->reverse_offset + 1;
	if (len <= 0) {
		LOG(PANDA_DEBUG_BUILD, PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH);
		return false;
	}
	if (len > 2 * PANDA_MAX_LEN) {
		LOG(PANDA_DEBUG_BUILD, PANDA_CODE_SEQUENCE_TOO_LONG);
		return false;
	}
	result->sequence_length = len - 1;
	result->degenerates = 0;

	df = (ssize_t) result->forward_length - (ssize_t) result->forward_offset - bestoverlap;
	dr = (ssize_t) result->reverse_length - (ssize_t) result->reverse_offset - bestoverlap;
	/* Copy the unpaired forward sequence. */
	LOGV(PANDA_DEBUG_RECON, PANDA_CODE_RECONSTRUCTION_PARAM, "bestoverlap = %zd, dforward = %zd, dreverse = %zd, len = %zd", bestoverlap, df, dr, len);
	for (i = 0; i < VEEZ(df); i++) {
		int findex = i + result->forward_offset;
		panda_nt fbits = result->forward[findex].nt;
		double q = qual_score[PHREDCLAMP(result->forward[findex].qual)];
		result->sequence[i].nt = fbits;
		result->sequence[i].p = q;
		if (PANDA_NT_IS_DEGN(fbits)) {
			result->degenerates++;
		}
		fquality += q;
		LOGV(PANDA_DEBUG_RECON, PANDA_CODE_BUILD_FORWARD, "S[%zd] = F[%d] = %c", i, findex, panda_nt_to_ascii(result->sequence[i].nt));
	}

	/* Mask out the B-cliff at the end of sequences */
	for (unmasked_forward_length = result->forward_length; unmasked_forward_length > 0 && result->forward[unmasked_forward_length - 1].qual == (char) 2; unmasked_forward_length--) ;
	for (unmasked_reverse_length = result->reverse_length; unmasked_reverse_length > 0 && result->reverse[unmasked_reverse_length - 1].qual == (char) 2; unmasked_reverse_length--) ;

	/* Copy the paired sequence adjusting the probabilities based on the quality information from both sequences. */
	result->overlap_mismatches = 0;
	for (i = 0; i < bestoverlap + WEDGEZ(df) + WEDGEZ(dr); i++) {
		int index = VEEZ(df) + i;
		int findex = result->forward_offset + VEEZ(df) + i;
		int rindex = result->reverse_length - i - 1 + WEDGEZ(df);
		bool ismatch = (result->reverse[rindex].nt & result->forward[findex].nt) != '\0';
		double fpr;
		double rpr;
		double q;
		char nt;

		if (index < 0 || findex < 0 || rindex < 0 || findex >= result->forward_length || rindex >= result->reverse_length)
			continue;

		fpr = findex >= unmasked_forward_length ? qual_nn : qual_score[PHREDCLAMP(result->forward[findex].qual)];
		rpr = rindex >= unmasked_reverse_length ? qual_nn : qual_score[PHREDCLAMP(result->reverse[rindex].qual)];

		if (!ismatch) {
			LOGV(PANDA_DEBUG_MISMATCH, PANDA_CODE_MISMATCHED_BASE, "(F[%d] = %c) != (R[%d] = %c)", findex, panda_nt_to_ascii(result->forward[findex].nt), rindex, panda_nt_to_ascii(result->reverse[rindex].nt));
			result->overlap_mismatches++;
		}

		if (findex >= unmasked_forward_length && rindex >= unmasked_reverse_length) {
			q = qual_nn;
		} else if (findex >= unmasked_forward_length) {
			q = rpr;
		} else if (rindex >= unmasked_reverse_length) {
			q = fpr;
		} else {
			q = match_probability(algo_data, ismatch, result->forward[findex].qual, result->reverse[rindex].qual);
		}

		if (ismatch) {
			nt = (result->reverse[rindex].nt & result->forward[findex].nt);
		} else {
			if (result->forward[findex].qual < result->reverse[rindex].qual) {
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
		LOGV(PANDA_DEBUG_RECON, PANDA_CODE_BUILD_OVERLAP, "S[%d] = %c, F[%d] = %c, R[%d] = %c", index, panda_nt_to_ascii(nt), findex, panda_nt_to_ascii(result->forward[findex].nt), rindex, panda_nt_to_ascii(result->reverse[rindex].nt));
	}

	/* Copy the unpaired reverse sequence. */
	for (i = 0; i < VEEZ(dr); i++) {
		int index = df + bestoverlap + i;
		int rindex = result->reverse_length - bestoverlap - i - 1;
		panda_nt rbits = result->reverse[rindex].nt;
		double q = qual_score[PHREDCLAMP(result->reverse[rindex].qual)];
		rquality += q;
		result->sequence[index].nt = rbits;
		result->sequence[index].p = q;
		if (PANDA_NT_IS_DEGN(rbits)) {
			result->degenerates++;
		}
		LOGV(PANDA_DEBUG_RECON, PANDA_CODE_BUILD_REVERSE, "S[%d] = R[%d] = %c", index, rindex, panda_nt_to_ascii(result->sequence[index].nt));
	}
	result->quality = (fquality + rquality + oquality) / len;

	result->overlap = bestoverlap;
	result->estimated_overlap_probability = bestprobability;

	return true;
}

bool assemble_seq(
	PandaAssembler assembler) {
	assembler->count++;
	if (assembler->result.forward_length < 2 || assembler->result.reverse_length < 2) {
		assembler->badreadcount++;
		return false;
	}
	if (!module_precheckseq(assembler, &assembler->result.name, assembler->result.forward, assembler->result.forward_length, assembler->result.reverse, assembler->result.reverse_length)) {
		return false;
	}
	if (!assembler->post_primers) {
		if (assembler->forward_primer_length > 0) {
			assembler->result.forward_offset = panda_compute_offset_qual(assembler->threshold, false, assembler->result.forward, assembler->result.forward_length, assembler->forward_primer, assembler->forward_primer_length);
			if (assembler->result.forward_offset == 0) {
				LOG(PANDA_DEBUG_STAT, PANDA_CODE_NO_FORWARD_PRIMER);
				assembler->nofpcount++;
				return false;
			}
			assembler->result.forward_offset--;
		} else {
			assembler->result.forward_offset = assembler->forward_trim;
		}
		if (assembler->reverse_primer_length > 0) {
			assembler->result.reverse_offset = panda_compute_offset_qual(assembler->threshold, false, assembler->result.reverse, assembler->result.reverse_length, assembler->reverse_primer, assembler->reverse_primer_length);
			if (assembler->result.reverse_offset == 0) {
				LOG(PANDA_DEBUG_STAT, PANDA_CODE_NO_REVERSE_PRIMER);
				assembler->norpcount++;
				return false;
			}
			assembler->result.reverse_offset--;
		} else {
			assembler->result.reverse_offset = assembler->reverse_trim;
		}
	} else {
		assembler->result.forward_offset = 0;
		assembler->result.reverse_offset = 0;
	}
	if (((assembler->result.forward_length < assembler->result.reverse_length) ? assembler->result.forward_length : assembler->result.reverse_length) < assembler->minoverlap) {
		assembler->badreadcount++;
		return false;
	}
	if (!align(assembler, &assembler->result)) {
		if (assembler->noalgn != NULL) {
			assembler->noalgn(assembler, &assembler->result.name, assembler->result.forward, assembler->result.forward_length, assembler->result.reverse, assembler->result.reverse_length, assembler->noalgn_data);
		}
		assembler->noalgncount++;
		return false;
	}
	if (assembler->post_primers) {
		size_t it;
		if (assembler->forward_primer_length > 0) {
			assembler->result.forward_offset = panda_compute_offset_result(assembler->threshold, false, assembler->result.sequence, assembler->result.sequence_length, assembler->forward_primer, assembler->forward_primer_length);
			if (assembler->result.forward_offset == 0) {
				LOG(PANDA_DEBUG_STAT, PANDA_CODE_NO_FORWARD_PRIMER);
				assembler->nofpcount++;
				return false;
			}
			assembler->result.forward_offset--;
		} else {
			assembler->result.forward_offset = assembler->forward_trim;
		}
		if (assembler->reverse_primer_length > 0) {
			assembler->result.reverse_offset = panda_compute_offset_result(assembler->threshold, true, assembler->result.sequence, assembler->result.sequence_length, assembler->reverse_primer, assembler->reverse_primer_length);
			if (assembler->result.reverse_offset == 0) {
				LOG(PANDA_DEBUG_STAT, PANDA_CODE_NO_REVERSE_PRIMER);
				assembler->norpcount++;
				return false;
			}
			assembler->result.reverse_offset--;
		} else {
			assembler->result.reverse_offset = assembler->reverse_trim;
		}
		if (assembler->result.sequence_length <= assembler->result.forward_offset + assembler->result.reverse_offset) {
			LOG(PANDA_DEBUG_STAT, PANDA_CODE_NO_FORWARD_PRIMER);
			assembler->nofpcount++;
			return false;
		}
		assembler->result.sequence_length -= assembler->result.forward_offset + assembler->result.reverse_offset;
		for (it = 0; it < assembler->result.sequence_length; it++) {
			assembler->result.sequence[it] = assembler->result.sequence[it + assembler->result.forward_offset];
		}
	}
	if (assembler->result.quality < assembler->threshold) {
		assembler->lowqcount++;
		LOGV(PANDA_DEBUG_STAT, PANDA_CODE_LOW_QUALITY_REJECT, "%f < %f", exp(assembler->result.quality), exp(assembler->threshold));
		return false;
	}
	if (module_checkseq(assembler, &assembler->result)) {
		assembler->okcount++;
		assembler->overlapcount[assembler->result.overlap]++;
		if (assembler->longest_overlap < assembler->result.overlap) {
			assembler->longest_overlap = assembler->result.overlap;
		}
		return true;
	}
	return false;
}

const panda_result_seq *panda_assembler_next(
	PandaAssembler assembler) {
	if (assembler->next == NULL) {
		return NULL;
	}
	while (true) {
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

const panda_result_seq *panda_assembler_assemble(
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length) {
	size_t it;
	assert(forward_length <= PANDA_MAX_LEN);
	assert(reverse_length <= PANDA_MAX_LEN);
	assembler->result.name = *id;
	assembler->result.forward_length = forward_length;
	assembler->result.reverse_length = reverse_length;
	assembler->result.forward = forward;
	assembler->result.reverse = reverse;
	return assemble_seq(assembler) ? &assembler->result : NULL;
}
