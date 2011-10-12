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

#ifndef _RESULTSEQ_H
#define _RESULTSEQ_H

/* Illumina sequence information from the FASTQ header */
typedef struct {
	char instrument[100];
	int run;
	char flowcell[100];
	int lane;
	int tile;
	int x;
	int y;
	char tag[6];
} seqidentifier;

/* Structure containing a sequence */
typedef struct {
	/* Length of the final sequence. */
	size_t len;
	/* Calculated quality score as the geometric mean of the product of the Illumina quality scores of the included bases. It will always be between 0 and 1. */
	double quality;
	/* Number of uncalled bases in the sequence. */
	size_t degenerates;
	/* Name of the sequence as provided in the input except lacking the forward/reverse identifier (i.e., /1 or /2). */
	seqidentifier *name;
	/* An array of log probabilities for each base the assembled sequence. */
	double *scores;
	/* The sequence, in uppercase ASCII. */
	char sequence[];
} resultseq;

/* Convenience macros for creating modules */
#define PANDA_API 1
#define PANDACONCATE(x,y) x ## y
#define PANDACONCAT(x,y) PANDACONCATE(x, y)
#define PRECHECK int PANDACONCAT(PANDASEQ_MODULE,_LTX_precheck) (const seqidentifier *id, char const *forward, char const *reverse)
#define CHECK int PANDACONCAT(PANDASEQ_MODULE,_LTX_api) = PANDA_API; int PANDACONCAT(PANDASEQ_MODULE,_LTX_check) (const resultseq * sequence)
#define INIT int PANDACONCAT(PANDASEQ_MODULE,_LTX_init)(const char *args)
#define CLEANUP void PANDACONCAT(PANDASEQ_MODULE,_LTX_destroy)()
#define HELP(desc, usage) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_desc) = desc; const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_usage) = usage
#define VER_INFO(version) const char *PANDACONCAT(PANDASEQ_MODULE,_LTX_version) = version;
#endif
