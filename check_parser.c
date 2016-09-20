#define _POSIX_C_SOURCE 2
#include<ctype.h>
#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include "config.h"
#include "pandaseq.h"

typedef struct {
	const char *str;
	int dir;
	PandaIdFmt format;
	panda_seq_identifier id;
} test_case;

bool check(
	const test_case * test) {
	PandaIdFmt detected_format;
	panda_seq_identifier id;
	panda_seqid_clear(&id);
	int dir = panda_seqid_parse_fail(&id, test->str, PANDA_TAG_OPTIONAL, &detected_format, NULL);
	return panda_seqid_equal(&id, &test->id) && dir == test->dir && detected_format == test->format;
}

const test_case checks[] = {
	{"M01271:10:000000000-A3WGH:1:1101:18786:6175 1:N:0:1", 1, PANDA_IDFMT_CASAVA_1_7, {"M01271", "10", "000000000-A3WGH", 1, 1101, 18786, 6175, "1"}},
	{"ILLUMINA-BE9C3F:29:FC:3:1:2462:1120 1:N:0:GCTATA", 1, PANDA_IDFMT_CASAVA_1_7, {"ILLUMINA-BE9C3F", "29", "FC", 3, 1, 2462, 1120, "GCTATA"}},
	{"M00958:47:000000000-A3GH7:1:1101:15028:1512 2:N:0:3", 2, PANDA_IDFMT_CASAVA_1_7, {"M00958", "47", "000000000-A3GH7", 1, 1101, 15028, 1512, "3"}},
	{"1468:1:1:12675:1118#ATCACGA/1", 1, PANDA_IDFMT_CASAVA_1_4, {"1468", "", "", 1, 1, 12675, 1118, "ATCACGA"}},
	{"1468:1:1:12675:1118#ATCACGA/2", 2, PANDA_IDFMT_CASAVA_1_4, {"1468", "", "", 1, 1, 12675, 1118, "ATCACGA"}},
	{"MISEQ03:18:000000000-A1REG:1:1101:14774:1712#GATAGTGCCAC/1", 1, PANDA_IDFMT_CASAVA_CONVERTED, {"MISEQ03", "18", "000000000-A1REG", 1, 1101, 14774, 1712, "GATAGTGCCAC"}}
};

int main(
	) {
	int exit_code = 0;
	for (size_t it = 0; it < sizeof(checks) / sizeof(*checks); it++) {
		if (!check(checks + it)) {
			fprintf(stderr, "FAILED: %s\n", checks[it].str);
			exit_code = 1;
		}
	}
	return exit_code;
}
