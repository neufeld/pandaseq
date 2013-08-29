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
#include<ctype.h>
#include<errno.h>
#include<float.h>
#include<limits.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq.h"
int main(
	int argc,
	char **argv) {
	int c;
	bool help = false;
	panda_seq_identifier id;
	int it;
	const char *optlist = "hv";
	bool version = false;

	while ((c = getopt(argc, argv, optlist)) != -1) {
		char *endptr;
		switch (c) {
		case 'h':
			help = true;
			break;
		case 'v':
			version = true;
			break;
		case '?':
			if (strchr(optlist, optopt) != NULL) {
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			} else if (isprint(optopt)) {
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			} else {
				fprintf(stderr, "Unknown option character `\\x%x'.\n", (unsigned int) optopt);
			}
			return 1;
		default:
			abort();
		}
	}

	if (version) {
		fprintf(stderr, "%s <%s>\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
		return 1;
	}
	if (argc < 2 || help) {
		fprintf(stderr, "%s <%s>\nUsage: %s \"seq header\" ...\nCheck is the sequence header is recognised by PANDAseq.\n", PACKAGE_STRING, PACKAGE_BUGREPORT, argv[0]);
		return 1;
	}

	for (it = 1; it < argc; it++) {
		const char *endptr;
		int dir;
		PandaIdFmt detected_format;
		panda_seqid_clear(&id);
		dir = panda_seqid_parse_fail(&id, argv[it] + (argv[it][0] == '@' ? 1 : 0), PANDA_TAG_OPTIONAL, &detected_format, &endptr);
		if (dir == 0) {
			int count;
			printf("%s\n", argv[it]);
			for (count = endptr - argv[it] - 2; count > 0; count--) {
				putchar(' ');
			}
			printf("^\n\tBAD\n");
		} else {
			panda_seqid_print(&id, stdout);
			printf("\n\tGOOD\n\tdirection = %s\n\thastag = %s\n", dir == 1 ? "forward" : "reverse", id.tag[0] == '\0' ? "no" : "yes");
		}
		printf("\tinstrument = \"%s\"\n\trun = %d\n\tflowcell = \"%s\"\n\tlane = %d\n\ttile = %d\n\tx = %d\n\ty = %d\n\ttag = \"%s\"\n\tgenerator = %s\n", id.instrument, id.run, id.flowcell, id.lane, id.tile, id.x, id.y, id.tag, panda_idfmt_str(detected_format));
	}
	return 0;
}
