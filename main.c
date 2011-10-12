/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011  Andre Masella

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
#include<bzlib.h>
#include<ctype.h>
#include<errno.h>
#include<fcntl.h>
#include<float.h>
#include<ltdl.h>
#include<limits.h>
#include<math.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>
#include<time.h>
#include<unistd.h>
#include<zlib.h>
#include "config.h"
#include "kseq.h"
#include "parser.h"
#include "prob.h"
#include "plugin.h"
#include "pandaseq.h"
#include "table.h"

#define RESULT_SIZE 4096
#define CIRC(index, len) (((index) + (len)) % (len))
double threshold;
int minoverlap = 1;
double pmatch;
double pmismatch;

#define STR0(x) #x
#define STR(x) STR0(x)
/* Function pointers for file I/O such that we can deal with compressed files. */
void *(*fileopen) (char *, char *) = (void *(*)(char *, char *))gzopen;
int (*fileread) (void *, void *, int) = (int (*)(void *, void *, int))gzread;
int (*fileclose) (void *) = (int (*)(void *))gzclose;

#define NT_A ((char)1)
#define NT_C ((char)2)
#define NT_G ((char)4)
#define NT_T ((char)8)
#define IS_N_NT(n) ((n) == (char)0x0F)

char iupac_forward[32] =
    { /* @ */ 0, /*A*/ NT_A, /*B*/ NT_C | NT_G | NT_T, /*C*/ NT_C,
	 /*D*/ NT_A | NT_G | NT_T, /*E*/ 0, /*F*/ 0, /*G*/ NT_G,
	 /*H*/ NT_A | NT_C | NT_T,
	 /*I*/ 0, /*J*/ 0, /*K*/ NT_G | NT_T, /*L*/ 0, /*M*/ NT_A | NT_C,
	 /*N*/ NT_A | NT_C | NT_G | NT_T, /*O*/ 0, /*P*/ 0, /*Q*/ 0,
	 /*R*/ NT_A | NT_G,
	 /*S*/ NT_C | NT_G, /*T*/ NT_T, /*U*/ NT_T, /*V*/ NT_A | NT_C | NT_G,
	 /*W*/ NT_A | NT_T,
	 /*X*/ NT_A | NT_C | NT_G | NT_T, /*Y*/ NT_C | NT_T, /*Z*/ 0, /*[ */ 0, /*\ */ 0,	/*] */
	0, /*^ */ 0, /*_*/ 0
};

char iupac_reverse[32] =
    { /*@ */ 0, /*A*/ NT_T, /*B*/ NT_G | NT_C | NT_A, /*C*/ NT_G,
	 /*D*/ NT_T | NT_C | NT_A, /*E*/ 0, /*F*/ 0, /*G*/ NT_C,
	 /*H*/ NT_T | NT_G | NT_A,
	 /*I*/ 0, /*J*/ 0, /*K*/ NT_C | NT_A, /*L*/ 0, /*M*/ NT_T | NT_G,
	 /*N*/ NT_A | NT_C | NT_G | NT_T, /*O*/ 0, /*P*/ 0, /*Q*/ 0,
	 /*R*/ NT_T | NT_C,
	 /*S*/ NT_G | NT_C, /*T*/ NT_A, /*U*/ NT_A, /*V*/ NT_T | NT_G | NT_C,
	 /*W*/ NT_T | NT_A,
	 /*X*/ NT_A | NT_C | NT_G | NT_T, /*Y*/ NT_G | NT_A, /*Z*/ 0, /*[ */ 0, /*\ */ 0,	/*] */
	0, /*^ */ 0, /*_*/ 0
};

#define NTCHAR(val) (ntchar[(int)(val)])
char ntchar[16] =
    { 'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B',
	'N'
};

/* Compatibility function to make BZ2_bzRead look like gzread. */
int bzread(BZFILE * file, void *buf, int len)
{
	int bzerror = BZ_OK;
	int retval = BZ2_bzRead(&bzerror, file, buf, len);
	if (bzerror == BZ_OK || bzerror == BZ_STREAM_END) {
		return retval;
	} else {
		fprintf(stderr, "ERR\tBZIP\t%d\n", bzerror);
		return -1;
	}
}

int qualmin = 33;

KSEQ_INIT(void *, fileread)
#define TOINDEX(val) (((int)(val)) < qualmin ? 0 : ((((int)(val)) > qualmin + PHREDMAX ? PHREDMAX : (int)(val)) - qualmin))
double qualscore(char nt1, char qual1, char nt2, char qual2)
{
	if (IS_N_NT(nt1)) {
		if (IS_N_NT(nt2)) {
			return qual_nn;
		}
		return qual_nmatch[TOINDEX(qual2)];
	}
	if (IS_N_NT(nt2)) {
		return qual_nmatch[TOINDEX(qual1)];
	}
	return (nt1 ==
		nt2 ? qual_match :
		qual_mismatch)[TOINDEX(qual1)][TOINDEX(qual2)];
}

#define BITS_SET(v) ((((unsigned int)(v)) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf)
/* Find the offset of a primer in a sequence and return the offset (i.e., the start of the useful sequence. */
size_t computeoffset(kseq_t * seq, char *primer, size_t primerlen)
{
	/* Circular buffer of probabilities of primer alignment indexed by the offset. */
	double probabilities[primerlen];
	double bestpr = primerlen * threshold;
	size_t bestindex = 0;
	size_t index;

	if (primerlen > seq->seq.l) {
		return 0;
	}

	for (index = 0; index < primerlen; index++) {
		probabilities[index] = -INFINITY;
	}

	for (index = 0; index < seq->seq.l; index++) {
		ssize_t x;
		for (x = (ssize_t) (primerlen > index ? index : primerlen - 1);
		     x >= 0; x--) {
			probabilities[CIRC(index - x, primerlen)] +=
			    qualscore(seq->seq.s[index], seq->qual.s[index],
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

#define KMER_LEN 8
typedef unsigned int bitstype;
typedef unsigned char seqindex;
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
#define _FOREACH_KMER(iterator,sequence,start,check,step) for ((iterator).posn = (start), (iterator).bad = KMER_LEN; (iterator).posn check; (iterator).posn step, (iterator).kmer = (((iterator).kmer << 2) | ((sequence)->seq.s[(iterator).posn] == NT_T ? 3 : (sequence)->seq.s[(iterator).posn] == NT_G ? 2 : (sequence)->seq.s[(iterator).posn] == NT_C ? 1 : 0)) & ((1 << (2 * KMER_LEN)) - 1)) if (IS_N_NT((sequence)->seq.s[(iterator).posn])) { (iterator).bad = KMER_LEN; } else if ((iterator).bad > 0) { (iterator).bad--; } else
#define FOREACH_KMER(iterator,sequence) _FOREACH_KMER(iterator,sequence, 0, < (sequence)->seq.l, ++)
#define FOREACH_KMER_REVERSE(iterator,sequence) _FOREACH_KMER(iterator,sequence, (sequence)->seq.l - 1, >= 0, --)
#define KMER(kmerit) ((kmerit).kmer)
#define KMER_POSITION(kmerit) ((kmerit).posn)
#define NUM_KMERS 2
#define KMERSEEN_SIZE (sizeof(seqindex) * NUM_KMERS * (1 << (2 * KMER_LEN)))
#define VEEZ(x) ((x) < 0 ? 0 : (x))
#define WEDGEZ(x) ((x) > 0 ? 0 : (x))
seqindex *kmerseen = NULL;
/* Try to align forward and reverse reads and return the quality of the aligned sequence and the sequence itself. */
int
align(resultseq * result, int maxresult, kseq_t * forward, size_t foffset,
      kseq_t * reverse, size_t roffset)
{
	ssize_t i, j;
	ssize_t df, dr;
	/* For determining overlap. */
	size_t maxoverlap =
	    forward->seq.l < reverse->seq.l ? reverse->seq.l : reverse->seq.l;
	double bestprobability = qual_nn * (forward->seq.l + reverse->seq.l);
	int bestoverlap = -1;
	size_t overlap;
	size_t counter;
	kmer_it it;

	/* For computing new sequence. */
	double fquality = 0;
	double oquality = 0;
	double rquality = 0;
	ssize_t len;
	BITS_INIT(posn, maxoverlap - minoverlap + 1);

	if (forward->seq.l >= (1 << (8 * sizeof(seqindex)))) {
		fprintf(stderr, "ERR\tKLNG\n");
		return 0;
	}

	/* Scan forward sequence building k-mers and appending the position to kmerseen[k] */
	FOREACH_KMER(it, forward) {
#ifdef DEBUG
		fprintf(stderr, "DBG\tFMER\t%x @ %d\n", KMER(it),
			KMER_POSITION(it));
#endif
		for (j = 0;
		     j < NUM_KMERS && kmerseen[(KMER(it) << 1) + j] != 0; j++) ;
		if (j == NUM_KMERS) {
			/* If we run out of storage, we lose k-mers. */
			fprintf(stderr, "DBG\tFML\t%x\n", KMER(it));
		} else {
			kmerseen[(KMER(it) << 1) + j] = KMER_POSITION(it);
		}
	}

	/* Scan reverse sequence building k-mers. For each position in the forward sequence for this kmer (i.e., kmerseen[k]), flag that we should check the corresponding overlap. */
	FOREACH_KMER_REVERSE(it, reverse) {
#ifdef DEBUG
		fprintf(stderr, "DBG\tRMER\t%x @ %d\n", KMER(it),
			KMER_POSITION(it));
#endif
		for (j = 0;
		     j < NUM_KMERS && kmerseen[(KMER(it) << 1) + j] != (seqindex)0; j++) {
			int index =
			    forward->seq.l + reverse->seq.l -
			    KMER_POSITION(it) - kmerseen[(KMER(it) << 1) + j] -
			    minoverlap - 1;

			if (index >= 0) {
				BIT_LIST_SET(posn, index);
			} else {
				fprintf(stderr, "DBG\tIOR\t%d\n", index);
			}
		}
	}

	/* Reset kmerseen */
	FOREACH_KMER(it, forward) {
		for (j = 0; j < NUM_KMERS; j++)
			kmerseen[(KMER(it) << 1) + j] = 0;
	}

	ALL_BITS_IF_NONE(posn, maxoverlap - minoverlap + 1);

	/* Compute the quality of the overlapping region for the various overlaps and pick the best one. */
	FOR_BITS_IN_LIST(posn, counter) {
		size_t matches = 0;
		size_t mismatches = 0;
		size_t unknowns = 0;
		double probability;
		overlap = counter + minoverlap;

		for (i = 0; i < overlap; i++) {
			int findex = forward->seq.l + i - overlap;
			int rindex = reverse->seq.l - i - 1;
			char f = forward->seq.s[findex];
			char r = reverse->seq.s[rindex];
			if (IS_N_NT(f) || IS_N_NT(r)) {
				unknowns++;
			} else if ((f & r) != 0) {
				matches++;
			} else {
				mismatches++;
			}
		}

		probability =
		    (qual_nn *
		     (forward->seq.l + reverse->seq.l - 2 * overlap +
		      unknowns) + matches * pmatch + mismatches * pmismatch);

#ifdef DEBUG
		fprintf(stderr,
			"INFO\tOLD\t%d mat = %d, mismat = %d, unk = %d, prob = %f\n",
			overlap, matches, mismatches, unknowns, probability);
#endif
		if (probability > bestprobability) {
			bestprobability = probability;
			bestoverlap = overlap;
		}
	}

	fprintf(stderr, "INFO\tBESTOLP\t%d\n", bestoverlap);

	if (bestoverlap == -1) {
		return 0;
	}

	/* Compute the correct alignment and the quality score of the entire sequence. */
	len =
	    forward->seq.l - (ssize_t) foffset - bestoverlap + reverse->seq.l -
	    (ssize_t) roffset + 1;
	if (len <= 0) {
		fprintf(stderr, "ERR\tNEGS\n");
		return 0;
	}
	if (len > maxresult) {
		fprintf(stderr, "ERR\tOOM\n");
		return 0;
	}
	result->len = len;
	result->degenerates = 0;
	result->sequence[len - 1] = '\0';
	result->scores = (double *)&result->sequence[len];

	df = (ssize_t) forward->seq.l - (ssize_t) foffset - bestoverlap;
	dr = (ssize_t) reverse->seq.l - (ssize_t) roffset - bestoverlap;
	/* Copy the unpaired forward sequence. */
#ifdef DEBUG
	fprintf(stderr,
		"INFO\tRECR\tflen = %d rlen = %d overlap = %d df = %d dr = %d\n",
		(int)forward->seq.l, (int)reverse->seq.l, (int)bestoverlap,
		(int)df, (int)dr);
#endif
	for (i = 0; i < VEEZ(df); i++) {
		int findex = i + foffset;
		char fbits = forward->seq.s[findex];
		char f = NTCHAR(fbits);
		double q = qual_score[TOINDEX(forward->qual.s[findex])];
		result->sequence[i] = f;
		result->scores[i] = q;
		if (BITS_SET(fbits) != 1) {
			result->degenerates++;
		}
		fquality += q;
#ifdef DEBUG
		fprintf(stderr, "INFO\tBUILD\tr[%u] = %c[%u] | %e\n",
			(unsigned int)i, (int)f, (unsigned int)findex, exp(q));
#endif
	}

	/* Mask out the B-cliff at the end of sequences */
	for (i = forward->qual.l - 1;
	     i > 0 && forward->qual.s[i] == (char)(qualmin + 2); i--) {
		forward->qual.s[i] = '\0';
	}
	for (i = reverse->qual.l - 1;
	     i > 0 && reverse->qual.s[i] == (char)(qualmin + 2); i--) {
		reverse->qual.s[i] = '\0';
	}
	/* Copy the paired sequence adjusting the probabilities based on the quality information from both sequences. */
	for (i = 0; i < bestoverlap + WEDGEZ(df) + WEDGEZ(dr); i++) {
		int index = VEEZ(df) + i;
		int findex = foffset + VEEZ(df) + i;
		int rindex = reverse->seq.l - i - 1 + WEDGEZ(df);
		bool ismatch =
		    (reverse->seq.s[rindex] & forward->seq.s[findex]) != '\0';
		double fpr;
		double rpr;
		double q;
		char nt;

		if (index < 0 || findex < 0 || rindex < 0
		    || findex >= forward->seq.l || rindex >= reverse->seq.l)
			continue;

		fpr = forward->qual.s[findex] == '\0' ? qual_nn :
		    qual_score[TOINDEX(forward->qual.s[findex])];
		rpr = reverse->qual.s[rindex] == '\0' ? qual_nn :
		    qual_score[TOINDEX(reverse->qual.s[rindex])];

		if (!ismatch) {
			fprintf(stderr, "INFO\tMISM\t%d %d\n",
				(int)forward->qual.s[findex],
				(int)reverse->qual.s[rindex]);
		}

		if (forward->qual.s[findex] == '\0'
		    && reverse->qual.s[rindex] == '\0') {
			q = qual_nn;
		} else if (forward->qual.s[findex] == '\0') {
			q = ismatch ? rpr : qual_nn;
		} else if (reverse->qual.s[rindex] == '\0') {
			q = ismatch ? fpr : qual_nn;
		} else {
			q = (ismatch ? qual_match :
			     qual_mismatch)[TOINDEX(forward->
						    qual.s[rindex])][TOINDEX
								     (reverse->
								      qual.s
								      [rindex])];
		}

		if (ismatch) {
			nt = (reverse->seq.s[rindex] & forward->seq.s[findex]);
		} else {
			if (forward->qual.s[rindex] < reverse->qual.s[rindex]) {
				nt = reverse->seq.s[rindex];
			} else {
				nt = forward->seq.s[findex];
			}
		}
		result->sequence[index] = NTCHAR(nt);
		result->scores[index] = q;
		result->scores[index] = q;
		if (BITS_SET(nt) != 1) {
			result->degenerates++;
		}
		oquality += q;
#ifdef DEBUG
		fprintf(stderr,
			"INFO\tBUILD\tr[%d] = %c | %e given %c[%d] | %e and %c[%d] | %e\n",
			index, result->sequence[index], exp(q),
			(int)NTCHAR(forward->seq.s[findex]), findex, exp(fpr),
			(int)NTCHAR(reverse->seq.s[rindex]), rindex, exp(rpr));
#endif
	}

	/* Copy the unpaired reverse sequence. */
	for (i = 0; i < VEEZ(dr); i++) {
		int index = df + bestoverlap + i;
		int rindex = reverse->seq.l - bestoverlap - i - 1;
		char rbits = reverse->seq.s[rindex];
		double q = qual_score[TOINDEX(reverse->qual.s[rindex])];
		rquality += q;
		result->sequence[index] = NTCHAR(rbits);
		result->scores[index] = q;
		if (BITS_SET(rbits) != 1) {
			result->degenerates++;
		}
#ifdef DEBUG
		fprintf(stderr, "INFO\tBUILD\tr[%d] = %c[%d] | %e\n", index,
			(int)result->sequence[index], rindex, exp(q));
#endif
	}
	result->quality = (fquality + rquality + oquality) / len;

	return 1;
}

/* Check that the a sequence looks like a valid sequence and convert it to an encoded represetnation. */
size_t cleanseq(char *seq, char *name, char *table)
{
	size_t i = 0;
	while (seq[i] != '\0') {
		if ((seq[i] = table[((int)seq[i]) & 0x1F]) == '\0') {
			fprintf(stderr, "ERR\tBADNT\tposition = %u\t%s\n",
				(unsigned int)i, name);
			return 0;
		}
		i++;
	}
	return i;
}

void printtime(long count, time_t starttime)
{
	time_t now;
	time(&now);
	fprintf(stderr, "STAT\tTIME\t%s\nSTAT\tELAPSED\t%d\nSTAT\tREADS\t%ld\n",
		ctime(&now), (int)(now - starttime), count);
}

int main(int argc, char **argv)
{
	int help = 0;
	int version = 0;
	int c;
	int bzip = 0;
	int fastq = 0;
	char *forward_filename = NULL;
	char *reverse_filename = NULL;
	char *forward_primer = NULL;
	char *reverse_primer = NULL;
	void *fforward;
	void *freverse;
	kseq_t *fseq, *rseq;
	int flen, rlen;
	size_t fprimerlen = 0;
	size_t rprimerlen = 0;
	size_t minlen = 0;
	ssize_t maxlen = SSIZE_MAX;
	size_t foffset = 1;
	size_t roffset = 1;
	time_t starttime;
	long count = 0;
	long nofpcount = 0;
	long norpcount = 0;
	long lowqcount = 0;
	long shortcount = 0;
	long longcount = 0;
	long degencount = 0;
	long noalgncount = 0;
	long okcount = 0;
	double q = 0.36;
	int no_n = 0;
	resultseq *complete = malloc(RESULT_SIZE);
	int maxresult =
	    (RESULT_SIZE - sizeof(resultseq) - 1) / (sizeof(char) +
						     sizeof(double));
	time(&starttime);
	threshold = log(0.6);

	if (lt_dlinit() != 0) {
		fprintf(stderr, "ERR\tLTLD\tINIT\n");
		return 1;
	}
	if (lt_dladdsearchdir(STR(PKGLIBDIR)) != 0) {
		fprintf(stderr, "ERR\tLTLD\tPKGDIR\t%s\n", STR(PKGLIBDIR));
		return 1;
	}

	/* Process command line arguments. */
	while ((c = getopt(argc, argv, "hvjp:q:f:r:t:o:Nl:L:Q:C:6F")) != -1) {
		char *endptr;
		switch (c) {
		case 'h':
			help = 1;
			break;
		case 'v':
			version = 1;
			break;
		case 'j':
			fileopen = (void *(*)(char *, char *))BZ2_bzopen;
			fileread = (int (*)(void *, void *, int))bzread;
			fileclose = (int (*)(void *))BZ2_bzclose;
			bzip = 1;
			break;
		case 't':
			errno = 0;
			threshold = strtod(optarg, NULL);
			if (errno != 0 || threshold < 0 || threshold > 1) {
				fprintf(stderr,
					"Bad threshold. It should be between 0 and 1.\n");
				return 1;
			}
			threshold = log(threshold);
			break;
		case 'Q':
			errno = 0;
			q = strtod(optarg, NULL);
			if (errno != 0 || q < 0 || q > 1) {
				fprintf(stderr,
					"Bad quality. It should be between 0 and 1.\n");
				return 1;
			}
			threshold = log(threshold);
			break;
		case 'l':
			errno = 0;
			minlen = strtol(optarg, NULL, 10);
			if (errno != 0 || minlen < 0) {
				fprintf(stderr, "Bad minimum length.\n");
				return 1;
			}
			break;
		case 'L':
			errno = 0;
			maxlen = strtol(optarg, NULL, 10);
			if (errno != 0 || minlen < 0) {
				fprintf(stderr, "Bad maximum length.\n");
				return 1;
			}
			break;
		case 'o':
			errno = 0;
			minoverlap = strtol(optarg, NULL, 10);
			if (errno != 0 || minoverlap < 1) {
				fprintf(stderr, "Bad minimum overlap.\n");
				return 1;
			}
			break;
		case 'f':
			forward_filename = optarg;
			break;
		case 'r':
			reverse_filename = optarg;
			break;
		case 'N':
			no_n = 1;
			break;
		case 'F':
			fastq = 1;
			break;
		case 'p':
			errno = 0;
			foffset = strtol(optarg, &endptr, 10);
			if (*endptr != '\0') {
				forward_primer = optarg;
			} else if (errno != 0 || foffset < 1) {
				fprintf(stderr, "Bad forward primer length.\n");
				return 1;
			} else {
				foffset++;
			}
			break;
		case 'q':
			errno = 0;
			roffset = strtol(optarg, &endptr, 10);
			if (*endptr != '\0') {
				reverse_primer = optarg;
			} else if (errno != 0 || roffset < 1) {
				fprintf(stderr, "Bad reverse primer length.\n");
				return 1;
			} else {
				roffset++;
			}
			break;
		case 'C':
			if (!module_load(optarg)) {
				module_cleanup();
				return 1;
			}
			break;
		case '6':
			qualmin = 64;
			break;
		case '?':
			if (optopt == (int)'f' || optopt == (int)'r'
			    || optopt == (int)'p' || optopt == (int)'q'
			    || optopt == (int)'l' || optopt == (int)'L'
			    || optopt == (int)'Q' || optopt == (int)'C') {
				fprintf(stderr,
					"Option -%c requires an argument.\n",
					optopt);
			} else if (isprint(optopt)) {
				fprintf(stderr, "Unknown option `-%c'.\n",
					optopt);
			} else {
				fprintf(stderr,
					"Unknown option character `\\x%x'.\n",
					(unsigned int)optopt);
			}
			return 1;
		default:
			abort();
		}
	}

	if (version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		module_version();
		return 1;
	}
	if (forward_filename == NULL || reverse_filename == NULL || help) {
		fprintf(stderr,
			"%s <%s>\nUsage: %s -f forward.fastq -r reverse.fastq [-j] [-p forwardprimer] [-q reverseprimer] [-t threshold] [-N] [-o minoverlap] [-l minlen] [-L maxlen] [ -C module1 -C module2 ...] [-6] [-F]\n\t-f\tInput FASTQ file containing forward reads.\n\t-r\tInput FASTQ file containing reverse reads.\n\t-j\tInput files are bzipped.\n\t-p\tForward primer sequence or number of bases to be removed.\n\t-q\tReverse primer sequence or number of bases to be removed.\n\t-t\tThe minimum probability that a sequence must have to match a primer. (default = %e)\n\t-N\tEliminate all sequences with unknown nucleotides in the output.\n\t-o minoverlap\tMinimum overlap between forward and reverse reads (default = %d)\n\t-l minlen\tMinimum length for a sequence\n\t-L maxlen\tMaximum length for a sequence\n\t-C module\tLoad a sequence validation module.\n\t-6\tUse PHRED+64 (CASAVA 1.3-1.7) instead of PHRED+33 (CASAVA 1.8+).\n\t-F\tOutput FASTQ instead of FASTA.\n",
			PACKAGE_STRING, PACKAGE_BUGREPORT,
			argv[0], exp(threshold), minoverlap);
			module_help();
		return 1;
	}
	fprintf(stderr, "INFO\tVER\t%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
	if (!module_init()) {
		return 1;
	}
	pmatch = log(0.25 * (1 - 2 * q + q * q));
	pmismatch = log((3 * q - 2 * q * q) / 18.0);

	/* Open files and initialise FASTQ reader. */
	fforward = fileopen(forward_filename, "r");
	if (fforward == NULL) {
		perror(forward_filename);
		return 1;
	}

	freverse = fileopen(reverse_filename, "r");
	if (freverse == NULL) {
		perror(reverse_filename);
		return 1;
	}

	if (forward_primer != NULL) {
		fprimerlen =
		    cleanseq(forward_primer, "forward primer", iupac_forward);
		if (fprimerlen == 0) {
			return 1;
		}
	}
	if (reverse_primer != NULL) {
		rprimerlen =
		    cleanseq(reverse_primer, "reverse primer", iupac_reverse);
		if (rprimerlen == 0) {
			return 1;
		}
	}
	fseq = kseq_init(fforward);
	rseq = kseq_init(freverse);

	kmerseen = malloc(KMERSEEN_SIZE);
	memset(kmerseen, 0, KMERSEEN_SIZE);
	/* Process sequences... */
	while ((flen = kseq_read(fseq)) >= 0 && (rlen = kseq_read(rseq)) >= 0) {
		seqidentifier fid;
		seqidentifier rid;
		int fmate;
		int rmate;

		if (count % 1000 == 0) {
			printtime(count, starttime);
		}
		count++;

		/* Ensure reads are paired. */
		fmate = seqid_parse(&fid, fseq->name.s);
		rmate = seqid_parse(&rid, rseq->name.s);

		if (fmate == 0) {
			fprintf(stderr,
				"ERR\tSTOP\tCould not parse: \"%s\"\n", fseq->name.s);
			break;
		}
		if (rmate == 0) {
			fprintf(stderr,
				"ERR\tSTOP\tCould not parse: \"%s\"\n", rseq->name.s);
			break;
		}
		if (fmate == rmate) {
			fprintf(stderr,
				"ERR\tSTOP\tMate-pairs do not match: \"%s\" and \"%s\"\n",
				fseq->name.s, rseq->name.s);
			break;
		}
		if (!seqid_equal(&fid, &rid)) {
			fprintf(stderr,
				"ERR\tSTOP\tNames do not match: \"%s\" and \"%s\"\n",
				fseq->name.s, rseq->name.s);
			break;
		}

		/* Check there is sufficient quality information. */
		if (fseq->qual.l != fseq->seq.l) {
			fprintf(stderr, "ERR\tNOFQ\t%s\n", fseq->name.s);
			continue;
		}
		if (rseq->qual.l != rseq->seq.l) {
			fprintf(stderr, "ERR\tNORQ\t%s\n", fseq->name.s);
			continue;
		}
		if (!module_precheckseq(&fid, fseq->seq.s, rseq->seq.s)) {
			continue;
		}

		/* Check and encode the sequences. */
		if (cleanseq(fseq->seq.s, fseq->name.s, iupac_forward) == 0
		    || cleanseq(rseq->seq.s, fseq->name.s,
				iupac_reverse) == 0) {
			continue;
		}

		/* Try to find primer offsets. */
		if (forward_primer != NULL) {
			foffset =
			    computeoffset(fseq, forward_primer, fprimerlen);
			if (foffset == 0) {
				fprintf(stderr, "ERR\tNOFP\t%s\n",
					fseq->name.s);
				nofpcount++;
				continue;
			}
		}
		if (reverse_primer != NULL) {
			roffset =
			    computeoffset(rseq, reverse_primer, rprimerlen);
			if (roffset == 0) {
				fprintf(stderr, "ERR\tNORP\t%s\n",
					fseq->name.s);
				norpcount++;
				continue;
			}
		}

		/* Align reads and emit sequence if good enough. */
		if (!align
		    (complete, maxresult, fseq, foffset - 1, rseq,
		     roffset - 1)) {
			fprintf(stderr, "ERR\tNOALGN\t%s\n", fseq->name.s);
			noalgncount++;
		} else {
			complete->name = &fid;
			if (complete->quality < threshold) {
				fprintf(stderr, "ERR\tLOWQ\t%e\t%s\n",
					exp(complete->quality), fseq->name.s);
				lowqcount++;
			} else if (no_n && complete->degenerates > 0) {
				fprintf(stderr, "ERR\tDEGN\t%s\n",
					fseq->name.s);
				degencount++;
			} else if (complete->len < minlen) {
				fprintf(stderr, "ERR\tSHORT\t%s\n",
					fseq->name.s);
				shortcount++;
			} else if (complete->len > maxlen) {
				fprintf(stderr, "ERR\tLONG\t%s\n",
					fseq->name.s);
				longcount++;
			} else if (module_checkseq(complete)) {
				fprintf(stderr, "OK\t%e\t%s\n",
					exp(complete->quality), fseq->name.s);
				if (fastq) {
					int it;
					seqid_print(complete->name, '@');
					printf("%s\n", complete->sequence);
					seqid_print(complete->name, '+');
					for (it = 0; it < complete->len - 1;
					     it++) {
						putchar((int)
							(33 -
							 10 * log(1 -
								  exp(complete->
								      scores
								      [it])) /
							 M_LN10));
					}
					putchar('\n');
				} else {
					seqid_print(complete->name, '>');
					printf("%s\n", complete->sequence);
				}
				okcount++;
			}
		}
	}
	printtime(count, starttime);
	if (forward_primer != NULL)
		fprintf(stderr, "STAT\tNOFP\t%ld\n", nofpcount);
	if (reverse_primer != NULL)
		fprintf(stderr, "STAT\tNORP\t%ld\n", norpcount);
	fprintf(stderr, "STAT\tNOALGN\t%ld\nSTAT\tLOWQ\t%ld\n", noalgncount,
		lowqcount);
	if (no_n)
		fprintf(stderr, "STAT\tDEGENERATE\t%ld\n", degencount);
	if (minlen > 0)
		fprintf(stderr, "STAT\tSHORT\t%ld\n", shortcount);
	if (maxlen < SSIZE_MAX)
		fprintf(stderr, "STAT\tLONG\t%ld\n", longcount);
	fprintf(stderr, "STAT\tOK\t%ld\n", okcount);

	/* Clean up. */
	free(kmerseen);
	kseq_destroy(fseq);
	kseq_destroy(rseq);
	if (fileclose(fforward) != Z_OK && bzip == 0) {
		perror(forward_filename);
	}
	if (fileclose(freverse) != Z_OK && bzip == 0) {
		perror(reverse_filename);
	}
	module_cleanup();
	lt_dlexit();
	free(complete);
	return 0;
}
