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
#ifndef PLUGIN_H
#        define PLUGIN_H
#        include "pandaseq.h"

extern bool module_checkseq(
	PandaAssembler assembler,
	panda_result_seq *sequence);
extern bool module_precheckseq(
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length);
extern void module_help(
	PandaAssembler assembler);
extern void module_version(
	PandaAssembler assembler);
extern bool module_init(
	PandaAssembler assembler);
extern void module_cleanup(
	PandaAssembler assembler);
extern void module_destroy(
	PandaAssembler assembler);
#endif
