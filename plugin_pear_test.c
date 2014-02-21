#include<errno.h>
#include<math.h>
#include<stddef.h>
#include<stdlib.h>
#include<string.h>
#include<pandaseq-plugin.h>

struct data {
	double alpha;
	double beta;
	double cutoff;
};

HELP("Use the statistical test from PEAR (Zhang 2013)", "pear_test:alpha=1.0,beta=-1.0,cutoff=0.01");

VER_INFO("1.0");

static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {

	struct data *data = (struct data *) user_data;
	double product = 1;
	size_t i;
	double oes = data->alpha * (sequence->overlap - sequence->overlap_mismatches) + data->beta * sequence->overlap_mismatches;
	for (i = sequence->overlap; i < sequence->forward_length && i < sequence->reverse_length; i++) {
		double sum = 0;
		size_t k;
		size_t l_i = ceil((oes - data->beta * i) / (data->alpha - data->beta)) - 1;
		for (k = 0; k < l_i; k++) {
			double i_choose_k = lgamma(i + 1) - lgamma(k + 1) - lgamma(i - k + 1);
			sum += exp(i_choose_k + k * log(0.25) + (i - k) * log(0.75));
		}
		product *= sum;
	}
	return data->cutoff > 1 - product * product;
}

struct {
	const char *name;
	size_t holder;
} const token[] = {
	{.name = "alpha",.holder = offsetof(struct data, alpha)},
	{.name = "beta",.holder = offsetof(struct data, beta)},
	{.name = "cutoff",.holder = offsetof(struct data, cutoff)},
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
		panda_log_proxy_perror(logger, arg_name);
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
			return parse_argument((PandaLogProxy) data, value, token[it].name, (double *) ((char *) data + token[it].holder));
		}
	}
	panda_log_proxy_write_f((PandaLogProxy) data, "Unknown setting: /%s/\n", key);
	return false;
}

OPEN {
	char *value;
	struct data data;

	data.alpha = 1;
	data.beta = -1;
	data.cutoff = 0.01;

	if (!panda_parse_key_values(args, key_processor, logger))
		return false;
	if (data.cutoff < 0 || data.cutoff > 1) {
		panda_log_proxy_write_f(logger, "Value %f out of range for p-value cut-off.", data.cutoff);
		return false;
	}
	*check = check_func;
	*user_data = PANDA_STRUCT_DUP(&data);
	*destroy = free;
	return true;
}
