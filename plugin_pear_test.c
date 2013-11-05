#define _XOPEN_SOURCE 500
#include<errno.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<pandaseq-plugin.h>

double alpha = 1;
double beta = -1;
double cutoff = 0.01;

HELP("This is a sample module that does nothing", "pear_test:alpha=1.0,beta=-1.0,cutoff=0.01");

VER_INFO("1.0");

CHECK {
	double product = 1;
	size_t i;
	double oes = alpha * (sequence->overlap - sequence->overlap_mismatches) + beta * sequence->overlap_mismatches;
	for (i = sequence->overlap; i < sequence->forward_length && i < sequence->reverse_length; i++) {
		double sum = 0;
		size_t k;
		size_t l_i = ceil((oes - beta * i) / (alpha - beta)) - 1;
		for (k = 0; k < l_i; k++) {
			double i_choose_k = lgamma(i + 1) - lgamma(k + 1) - lgamma(i - k + 1);
			sum += exp(i_choose_k + k * log(0.25) + (i - k) * log(0.75));
		}
		product *= sum;
	}
	return cutoff > 1 - product * product;
}

enum {
	OPT_ALPHA = 0,
	OPT_BETA,
	OPT_CUTOFF
};

char *const token[] = {
	[OPT_ALPHA] = "alpha",
	[OPT_BETA] = "beta",
	[OPT_CUTOFF] = "cutoff",
	NULL
};

bool parse_argument(
	char *value,
	char *arg_name,
	double *output) {
	char *remainder;
	errno = 0;
	*output = strtod(value, &remainder);
	if (errno != 0) {
		perror(arg_name);
		return false;
	} else if (*remainder != '\0') {
		fprintf(stderr, "%s: trailing garbage: %s\n", arg_name, remainder);
		return false;
	}
	return true;
}

INIT {
	char *value;
	if (args == NULL)
		return true;
	switch (getsubopt(&args, token, &value)) {
	case OPT_ALPHA:
		if (!parse_argument(value, token[OPT_ALPHA], &alpha))
			return false;
		break;
	case OPT_BETA:
		if (!parse_argument(value, token[OPT_BETA], &beta))
			return false;
		break;
	case OPT_CUTOFF:
		if (!parse_argument(value, token[OPT_CUTOFF], &cutoff))
			return false;
		if (cutoff < 0 || cutoff > 1) {
			fprintf(stderr, "Value %f out of range for p-value cut-off.", cutoff);
			return false;
		}
		break;
	default:
		fprintf(stderr, "No match found for token: /%s/\n", value);
		return false;
	}
	return true;
}
