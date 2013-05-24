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
#define _POSIX_C_SOURCE 2
#include<stdio.h>
#include<stdlib.h>
#include "config.h"
#ifdef HAVE_PTHREAD
#        include<pthread.h>
#endif
#include "pandaseq.h"
#include "misc.h"

int main(
	int argc,
	char **argv) {
	MANAGED_STACK(PandaOutputSeq, output);
	PandaAssembler assembler;
	PandaArgsHang data = panda_args_hang_new(panda_args_fastq_new(), (PandaDestroy) panda_args_fastq_free, (PandaTweakGeneral) panda_args_fastq_tweak, (PandaOpener) panda_args_fastq_opener, (PandaSetup) panda_args_fastq_setup);
	const panda_tweak_general **general_args;
	size_t general_args_length;
	PandaMux mux;
	bool result;
	int threads;

	general_args = panda_args_hang_args(panda_args_fastq_args, panda_args_fastq_args_length, &general_args_length);

	if (!panda_parse_args(argv, argc, panda_stdargs, panda_stdargs_length, general_args, general_args_length, (PandaTweakGeneral) panda_args_hang_tweak, (PandaOpener) panda_args_hang_opener, (PandaSetup) panda_args_hang_setup, data, &assembler, &mux, &threads, &output, &output_data, &output_destroy)) {
		free(general_args);
		panda_args_hang_free(data);
		return 1;
	}
	free(general_args);
	result = panda_run_pool(threads, assembler, mux, output, output_data, output_destroy);
	panda_args_hang_free(data);
	return result ? 0 : 1;
}
