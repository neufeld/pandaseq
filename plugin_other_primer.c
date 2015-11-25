#include<pandaseq-plugin.h>
#include<math.h>

HELP("Remove reads with another primer. Use f for forward, r for reverse.", "other_primer:[fr]:NNNNN");

VER_INFO("1.0");

struct data {
	size_t primer_length;
	bool forward;
	panda_nt primer[1];
};

static bool precheck_func(
	PandaLogProxy logger,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	struct data *data) {

	return panda_compute_offset_qual(log(0.9), 0.01, !data->forward, data->forward ? forward : reverse, data->forward ? forward_length : reverse_length, data->primer, data->primer_length) == 0;
}

OPEN {
	struct data *data;
	bool forward;
	size_t it;

	if (args == NULL || *args == '\0') {
		return false;
	}
	if (*args == 'f' || *args == 'p') {
		forward = true;
	} else if (*args == 'r' || *args == 'q') {
		forward = false;
	} else {
		panda_log_proxy_write_f(logger, "ERR\tOTHER_PRIMER\tINIT\tExpected f or r, but got %c.\n", (int) *args);
		return false;
	}
	args++;
	if (*args != ':') {
		panda_log_proxy_write_f(logger, "ERR\tOTHER_PRIMER\tINIT\tExpected :, but got %c.\n", (int) *args);
		return false;
	}
	args++;
	if (*args == '\0') {
		panda_log_proxy_write_f(logger, "ERR\tOTHER_PRIMER\tINIT\tPrimer cannot be empty.\n");
		return false;
	}

	data = malloc(sizeof(struct data) + sizeof(panda_nt) * strlen(args));
	data->forward = forward;
	data->primer_length = strlen(args);
	for (it = 0; it < data->primer_length; it++) {
		if ((data->primer[it] = (forward ? panda_nt_from_ascii : panda_nt_from_ascii_complement) (args[it])) == '\0') {
			panda_log_proxy_write_f(logger, "ERR\tOTHER_PRIMER\tBADNT\t%c\n", (int) args[it]);
			free(data);
			return false;
		}
	}

	*precheck = (PandaPreCheck) precheck_func;
	*user_data = data;
	*destroy = (PandaDestroy) free;
	return true;
}
