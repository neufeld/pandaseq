/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2013  Andre Masella

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
#include <ltdl.h>
#include <pandaseq.h>
#define CLEANUP() if (handle != NULL) lt_dlclose(handle); (void)lt_dlexit();

PandaAssembler panda_assembler_new_from_file(
	const char *file_name,
	PandaLogProxy logger,
	const panda_result_seq *(*assemble) (PandaAssembler assemble,
		panda_seq_identifier *id,
		const panda_qual *forward,
		size_t forward_length,
		const panda_qual *reverse,
		size_t reverse_length)) {
	lt_dlhandle handle = NULL;
	lt_dladvise advise;
	PandaAssembler (
		*new_func) (
		PandaNextSeq next,
		void *next_data,
		PandaDestroy next_destroy,
		PandaLogProxy logger);

	if (lt_dlinit() != 0) {
		CLEANUP();
		return NULL;
	}
	lt_dladvise_init(&advise);
	lt_dladvise_ext(&advise);
	lt_dladvise_local(&advise);
	handle = lt_dlopenadvise(file_name, advise);
	lt_dladvise_destroy(&advise);
	if (handle == NULL) {
		CLEANUP();
		return NULL;
	}

	*(void **) (&new_func) = lt_dlsym(handle, "panda_assembler_new");
	*(void **) assemble = lt_dlsym(handle, "panda_assembler_assemble");
	if (new_func == NULL) {
		CLEANUP();
		return NULL;
	}
	return new_func(NULL, NULL, NULL, logger);
}
