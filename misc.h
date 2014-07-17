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
#ifndef MISC_H
#        define MISC_H

#        include <unistd.h>
#        include "pandaseq.h"

#        define DESTROY_MEMBER(self, name) if ((self)->name ## _destroy != NULL && (self)->name != NULL) { (self)->name ## _destroy((self)->name ## _data); } (self)->name = NULL; (self)->name ## _data = NULL; (self)->name ## _destroy = NULL
#        define DESTROY_STACK(name) if (name ## _destroy != NULL && name != NULL) { name ## _destroy(name ## _data); } name = NULL; name ## _data = NULL; name ## _destroy = NULL
#        define MANAGED_MEMBER(type, name) type name; void * name ## _data; PandaDestroy name ## _destroy
#        define MANAGED_STACK(type, name) type name = NULL; void * name ## _data = NULL; PandaDestroy name ## _destroy = NULL
#        define MAYBE(x) if (x != NULL) *x
#        define free0(val) if ((val) != NULL) free(val); (val) = NULL

typedef unsigned short seqindex;
#        define KMER_LEN ((size_t) 8)
#        define KMERSEEN_SIZE(num_kmers) (sizeof(seqindex) * (num_kmers) * (1 << (2 * KMER_LEN)))

typedef struct {
	size_t kmer;
	ssize_t posn;
	ssize_t bad;
} kmer_it;
#        define _FOREACH_KMER(iterator,sequence,suffix,start,badstart,check,step,badreset) for ((iterator).posn = (start), (iterator).kmer = 0, (iterator).bad = badstart; (iterator).posn check; (iterator).posn step) if ((iterator).kmer = (((iterator).kmer << 2) | ((sequence)[(iterator).posn]suffix == PANDA_NT_T ? 3 : (sequence)[(iterator).posn]suffix == PANDA_NT_G ? 2 : (sequence)[(iterator).posn]suffix == PANDA_NT_C ? 1 : 0)) & ((1 << (2 * badreset)) - 1), PANDA_NT_IS_N((sequence)[(iterator).posn]suffix)) { (iterator).bad = badreset; } else if ((iterator).bad > 0) { (iterator).bad--; } else
#        define FOREACH_KMER(iterator,sequence,suffix) _FOREACH_KMER(iterator, sequence, suffix, 0, KMER_LEN, < (ssize_t)sequence ## _length, ++, KMER_LEN)
#        define FOREACH_KMER_REVERSE(iterator,sequence,suffix) _FOREACH_KMER(iterator, sequence, suffix, sequence ## _length - 1, KMER_LEN, >= 0, --, KMER_LEN)
#        define KMER(kmerit) ((kmerit).kmer)
#        define KMER_POSITION(kmerit) ((kmerit).posn)

#endif
