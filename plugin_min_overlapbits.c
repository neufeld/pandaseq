#include<errno.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pandaseq-plugin.h>
#include"table.h"

#ifndef M_LN2
#        define M_LN2 0.69314718055994530942
#endif

HELP("Check the number of bits saved (Cole 2013).", "min_overlapbits:15");

VER_INFO("1.0");

static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {
	(void) logger;
	return *(double *) user_data <= sequence->estimated_overlap_probability;
}

OPEN {
	double bits_saved = 15 * M_LN2;	// change from bits to nats
	char *remainder = NULL;
	double orig_value;

	(void) precheck;

	if (args != NULL) {
		errno = 0;
		orig_value = strtod(args, &remainder);
		bits_saved = orig_value * M_LN2;	// change from bits to nats

		if (errno != 0) {
			panda_log_proxy_write_str(logger, "bits_saved");
			return false;
		} else if (*remainder != '\0') {
			panda_log_proxy_write_f(logger, "bits_saved: trailing garbage: %s\n", remainder);
			return false;
		}
		if (bits_saved < 0) {
			panda_log_proxy_write_f(logger, "Value %f out of range for bits saved cut-off.", orig_value);
			return false;
		}
	}
	*check = check_func;
	*user_data = PANDA_STRUCT_DUP(&bits_saved);
	*destroy = free;
	return true;
}
