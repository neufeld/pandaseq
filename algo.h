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

#ifndef ALGO_H
#        define ALGO_H
#        include "config.h"
#        include "pandaseq.h"
#        include "misc.h"
#        ifdef HAVE_PTHREAD
#                include <pthread.h>
#        endif

struct panda_algorithm {
	const struct panda_algorithm_class *clazz;
	volatile size_t refcnt;
#        ifdef HAVE_PTHREAD
	pthread_mutex_t mutex;
#        endif
	void *end;
};

#endif
