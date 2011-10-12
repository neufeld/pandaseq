/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011  Andre Masella

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
#ifndef _PLUGIN_H
#define _PLUGIN_H
#include "pandaseq.h"

extern int module_checkseq(resultseq * sequence);
extern int module_precheckseq(seqidentifier *id, char const *forward, char const *reverse);
extern void module_help();
extern void module_version();
extern void module_cleanup();
extern int module_load(char *path);
extern int module_init();
#endif
