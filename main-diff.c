/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2014  Andre Masella

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
#include "pandaseq.h"
#include "misc.h"

int main(
	int argc,
	char **argv) {
	MANAGED_STACK(PandaNextSeq,
		next);
	PandaAssembler control_assembler;
	PandaAssembler experimental_assembler;
	PandaArgsFastq data = panda_args_fastq_new();
	bool result;
	bool suppress_quality_diffs;

	if (!panda_diff_parse_args(argv, argc, panda_stdargs, panda_stdargs_length, panda_args_fastq_args, panda_args_fastq_args_length, (PandaTweakGeneral) panda_args_fastq_tweak, (PandaOpener) panda_args_fastq_opener, (PandaSetup) panda_args_fastq_setup, data, &control_assembler, &experimental_assembler, &next, &next_data, &next_destroy, &suppress_quality_diffs)) {
		panda_args_fastq_free(data);
		DESTROY_STACK(next);
		return 1;
	}
	result = panda_diff(next, next_data, (PandaAssemble) panda_assembler_assemble, control_assembler, (PandaAssemble) panda_assembler_assemble, experimental_assembler, suppress_quality_diffs);
	DESTROY_STACK(next);
	panda_assembler_unref(control_assembler);
	panda_assembler_unref(experimental_assembler);
	panda_args_fastq_free(data);
	return result ? 0 : 1;
}
