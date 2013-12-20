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

double bits_saved = 15 * M_LN2;	// change from bits to nats

HELP("Check the number of bits saved (Cole 2013).", "min_overlapbits:15");

VER_INFO("1.0");

CHECK {
	return bits_saved <= sequence->estimated_overlap_probability;
}

INIT {
	char *remainder = NULL;
	if (args == NULL) {
		return true;
	}
	errno = 0;
	double orig_value = strtod(args, &remainder);
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
	return true;
}
