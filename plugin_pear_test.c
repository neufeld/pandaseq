#include<errno.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pandaseq-plugin.h>

double alpha = 1;
double beta = -1;
double cutoff = 0.01;

HELP("Use the statistical test from PEAR (Zhang 2013)", "pear_test:alpha=1.0,beta=-1.0,cutoff=0.01");

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

struct {
	const char *name;
	double *holder;
} const token[] = {
	{.name = "alpha",.holder = &alpha},
	{.name = "beta",.holder = &beta},
	{.name = "cutoff",.holder = &cutoff},
	{NULL}
};

bool parse_argument(
	PandaLogProxy logger,
	const char *value,
	const char *arg_name,
	double *output) {
	char *remainder = NULL;
	errno = 0;
	*output = strtod(value, &remainder);
	if (errno != 0) {
		perror(arg_name);
		return false;
	} else if (*remainder != '\0') {
		panda_log_proxy_write_f(logger, "%s: trailing garbage: %s\n", arg_name, remainder);
		return false;
	}
	return true;
}

static bool key_processor(
	const char *key,
	const char *value,
	void *data) {
	size_t it;
	for (it = 0; token[it].name != NULL; it++) {
		if (strcmp(key, token[it].name) == 0) {
			return parse_argument((PandaLogProxy) data, value, token[it].name, token[it].holder);
		}
	}
	panda_log_proxy_write_f((PandaLogProxy) data, "Unknown setting: /%s/\n", key);
}

INIT {
	char *value;
	if (!panda_parse_key_values(args, key_processor, logger))
		return false;
	if (cutoff < 0 || cutoff > 1) {
		panda_log_proxy_write_f(logger, "Value %f out of range for p-value cut-off.", cutoff);
		return false;
	}
	return true;
}
