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

#ifndef PANDAASM_H
#define PANDAASM_H
#include "pandaseq.h"
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#define LOG(assembler, ...) (assembler)->logger((assembler)->logger_data, __VA_ARGS__);
#define DESTROY_MEMBER(self, name) if ((self)->name ## _destroy != NULL && (self)->name != NULL) { (self)->name ## _destroy((self)->name ## _data); } (self)->name = NULL; (self)->name ## _data = NULL; (self)->name ## _destroy = NULL
#define MANAGED_MEMBER(type, name) type name; void * name ## _data; PandaDestroy name ## _destroy
#define free0(val) if ((val) != NULL) free(val); (val) = NULL
typedef unsigned char seqindex;

struct panda_assembler {
	volatile size_t refcnt;

	MANAGED_MEMBER(PandaNextSeq, next);
	MANAGED_MEMBER(PandaLogger, logger);

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
	long count;
	bool no_n;
#ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#endif
	panda_nt forward_primer[PANDA_MAX_LEN];
	panda_nt reverse_primer[PANDA_MAX_LEN];
};

#endif
