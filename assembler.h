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

#ifndef ASM_H
#        define ASM_H
#        include "config.h"
#        include "pandaseq.h"
#        include "misc.h"
#        ifdef HAVE_PTHREAD
#                include <pthread.h>
#        endif

typedef unsigned short seqindex;
#        define KMER_LEN 8
#        define KMERSEEN_SIZE(num_kmers) (sizeof(seqindex) * (num_kmers) * (1 << (2 * KMER_LEN)))

struct panda_assembler {
	volatile size_t refcnt;

	 MANAGED_MEMBER(
		PandaNextSeq,
		next);
	 MANAGED_MEMBER(
		PandaLogger,
		logger);
	 MANAGED_MEMBER(
		PandaFailAlign,
		noalgn);

	size_t *rejected;
	PandaModule *modules;
	size_t modules_length;
	size_t modules_size;

	double threshold;
	int minoverlap;
	double q;
	double pmatch;
	double pmismatch;

	seqindex *kmerseen;
	size_t num_kmers;

	panda_result_seq result;

	size_t forward_primer_length;
	size_t reverse_primer_length;
	size_t forward_trim;
	size_t reverse_trim;

	long nofpcount;
	long norpcount;
	long okcount;
	long lowqcount;
	long degencount;
	long noalgncount;
	long badreadcount;
	long slowcount;
	long count;
	bool no_n;
	bool post_primers;
#        ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#        endif
	panda_nt forward_primer[MAX_LEN];
	panda_nt reverse_primer[MAX_LEN];
	panda_result result_seq[2 * MAX_LEN];
	long overlapcount[MAX_LEN];
	size_t longest_overlap;
};

#endif
