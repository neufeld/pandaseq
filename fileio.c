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
#include <bzlib.h>
#include <stdlib.h>
#include <zlib.h>
#if HAVE_PTHREAD
#include <pthread.h>
#endif
#include "pandaseq.h"

PandaAssembler panda_assembler_open_gz(char *forward, char *reverse,
				       PandaLogger logger, void *logger_data,
				       PandaDestroy logger_destroy,
				       unsigned char qualmin, PandaTagging policy)
{
	gzFile *forward_file;
	gzFile *reverse_file;

	forward_file = gzopen(forward, "r");
	if (forward_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, forward, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		return NULL;
	}
	reverse_file = gzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, reverse, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		gzclose(forward_file);
		return NULL;
	}

	return panda_assembler_new_fastq_reader(gzgetc, forward_file,
						(PandaDestroy) gzclose, gzgetc,
						reverse_file,
						(PandaDestroy) gzclose, logger,
						logger_data, logger_destroy,
						qualmin, policy);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_gz(char *forward, char *reverse, PandaLogger logger,
			   void *logger_data, PandaDestroy logger_destroy,
			   unsigned char qualmin, PandaTagging policy)
{
	gzFile *forward_file;
	gzFile *reverse_file;
	forward_file = gzopen(forward, "r");
	if (forward_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, forward, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		return NULL;
	}
	reverse_file = gzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, reverse, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		gzclose(forward_file);
		return NULL;
	}

	return panda_mux_new_fastq_reader(gzgetc, forward_file,
					  (PandaDestroy) gzclose, gzgetc,
					  reverse_file, (PandaDestroy) gzclose,
					  (PandaLogger) logger, logger_data,
					  (PandaDestroy) logger_destroy,
					  qualmin, policy);
}
#endif

static int bzgetc(BZFILE * file)
{
	char c;
	if (BZ2_bzread(file, &c, 1) != 1) {
		return EOF;
	}
	return c;
}

PandaAssembler panda_assembler_open_bz2(char *forward, char *reverse,
					PandaLogger logger, void *logger_data,
					PandaDestroy logger_destroy,
					unsigned char qualmin, PandaTagging policy)
{
	BZFILE *forward_file;
	BZFILE *reverse_file;
	forward_file = BZ2_bzopen(forward, "r");
	if (forward_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, forward, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		return NULL;
	}
	reverse_file = BZ2_bzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, reverse, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		BZ2_bzclose(forward_file);
		return NULL;
	}
	return panda_assembler_new_fastq_reader(bzgetc, forward_file,
						BZ2_bzclose, bzgetc,
						reverse_file, BZ2_bzclose,
						logger, logger_data,
						logger_destroy, qualmin, policy);
}

#ifdef HAVE_PTHREAD
PandaMux panda_mux_open_bz2(char *forward, char *reverse, PandaLogger logger,
			    void *logger_data, PandaDestroy logger_destroy,
			    unsigned char qualmin, PandaTagging policy)
{
	BZFILE *forward_file;
	BZFILE *reverse_file;
	forward_file = BZ2_bzopen(forward, "r");
	if (forward_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, forward, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		return NULL;
	}
	reverse_file = BZ2_bzopen(reverse, "r");
	if (reverse_file == NULL) {
		logger(PANDA_CODE_NO_FILE, NULL, reverse, logger_data);
		if (logger_destroy != NULL) {
			logger_destroy(logger_data);
		}
		BZ2_bzclose(forward_file);
		return NULL;
	}
	return panda_mux_new_fastq_reader(bzgetc, forward_file, BZ2_bzclose,
					  bzgetc, reverse_file, BZ2_bzclose,
					  (PandaLogger) logger, logger_data,
					  (PandaDestroy) logger_destroy,
					  qualmin, policy);
}
#endif
