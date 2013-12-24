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
		panda_log_proxy_write_str(logger, "Need a number for a PHRED score.\n");
		return false;
	}

	value = strtol(args, &endptr, 10);
	if (endptr != NULL && *endptr != '\0' || value < 0 || value > 127) {
		panda_log_proxy_write_str(logger, "PHRED score must be a number between 0 and 127.\n");
		return false;
	}
	min_score = value;
	return true;
}
