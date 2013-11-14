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
#include "test_wrapper.h"

extern PandaAssembler xxpanda_assembler_new(
	PandaNextSeq next,
	void *next_data,
	PandaDestroy next_destroy,
	PandaLogProxy logger);

extern void xxpanda_assembler_unref(
	PandaAssembler assembler);

extern const panda_result_seq *xxpanda_assembler_assemble(
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length);

Assemble panda_assembler_new_from_file(
	PandaLogProxy logger,
	void ** data
, PandaDestroy *destroy) {

	*data = xxpanda_assembler_new(NULL, NULL, NULL, logger);
	printf("New Assembler: %p\n", *data);
	*destroy = (PandaDestroy)xxpanda_assembler_unref;
	return xxpanda_assembler_assemble;
}

