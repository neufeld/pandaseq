#include<stdio.h>
#include<stdlib.h>
#include<pandaseq-plugin.h>

HELP("Ensure the minimum score of all the output bases is above a certain PHRED value.", "min_phred:value");

VER_INFO("1.0");

static char min_score;

CHECK {
	size_t it;
	for (it = 0; it < sequence->sequence_length; it++) {
		if (panda_result_phred(&sequence->sequence[it]) < min_score) {
			return false;
		}
	}
	return true;
}

INIT {
	long int value;
	char *endptr;
	if (args == NULL || *args == '\0') {
		fprintf(stderr, "Need a number for a PHRED score.\n", args);
		return false;
	}
	
	value = strtol(args, &endptr, 10);
	if (endptr != NULL && *endptr != '\0' || value < 0 || value > 127) {
		fprintf(stderr, "PHRED score must be a number between 0 and 127.\n", args);
		return false;
	}
	min_score = value;
	return true;
}
