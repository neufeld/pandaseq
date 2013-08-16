#define _POSIX_C_SOURCE 2
#include<errno.h>
#include<stdio.h>
#include<stdlib.h>
#include<pandaseq-plugin.h>

HELP("Filter out sequences that have mismatches in the overlap region.", "completely_miss_the_point:mismatches");

VER_INFO("1.0");

int mismatches;

CHECK {
	return sequence->overlap_mismatches <= mismatches;
}

INIT {
	if (args == NULL || *args == '\0') {
		fprintf(stderr, "Please supply the maximum allowed mismatches.\n");
		return false;
	}
	errno = 0;
	mismatches = (size_t) strtol(args, NULL, 10);
	if (errno != 0 || mismatches < 0 || mismatches > PANDA_MAX_LEN) {
		fprintf(stderr, "Bad maximum allowed mismatches.\n");
		return false;
	}
	return true;
}
