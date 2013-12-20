#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pandaseq-plugin.h>

HELP("Ensure the minimum read Q score of an assembled read is above certain value.", "min_readqscore:value");

VER_INFO("1.0");

static int min_score;

CHECK {
	size_t it;
	/* sum of the error probabilities of the assembled bases */

	double errorsum = 0;
	for (it = 0; it < sequence->sequence_length; it++) {
		errorsum += -sequence->sequence[it].p;
	}
	int readqscore = (int) floor(-10.0 * log10(errorsum / sequence->sequence_length));
	if (readqscore < min_score) {
		return false;
	}
	return true;
}

INIT {
	long int value;
	char *endptr;
	if (args == NULL || *args == '\0') {
		panda_log_proxy_write_str(logger, "Need a number for read Q score.\n");
		return false;
	}

	value = strtol(args, &endptr, 10);
	if (endptr != NULL && *endptr != '\0' || value < 0) {
		panda_log_proxy_write_str(logger, "Read Q score must be a number greater than 0.\n");
		return false;
	}
	min_score = value;
	return true;
}
