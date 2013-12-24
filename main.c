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
#include "pandaseq.h"
#include "misc.h"

int main(
	int argc,
	char **argv) {
	MANAGED_STACK(PandaOutputSeq, output);
	PandaAssembler assembler;
	PandaArgsFastq data = panda_args_fastq_new();
	PandaMux mux;
	bool result;
	int threads;

	if (!panda_parse_args(argv, argc, panda_stdargs, panda_stdargs_length, panda_args_fastq_args, panda_args_fastq_args_length, (PandaTweakGeneral) panda_args_fastq_tweak, (PandaOpener) panda_args_fastq_opener, (PandaSetup) panda_args_fastq_setup, data, &assembler, &mux, &threads, &output, &output_data, &output_destroy)) {
		panda_args_fastq_free(data);
		return 1;
	}
	result = panda_run_pool(threads, assembler, mux, output, output_data, output_destroy);
	panda_args_fastq_free(data);
	return result ? 0 : 1;
}
