#include<pandaseq-plugin.h>

HELP("Include only sequences in the one before the provided sequence", "before:sequenceid");

VER_INFO("1.0");

struct data {
	panda_seq_identifier marker_id;
	bool state;
};

static bool precheck_func(
	PandaLogProxy logger,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	void *user_data) {

	struct data *data = (struct data *) user_data;

	if (panda_seqid_equal(&data->marker_id, id)) {
		data->state = true;
	}
	return !data->state;
}

OPEN {
	struct data data;

	if (args == NULL) {
		panda_log_proxy_write_str(logger, "ERR\tBEFORE\tNO ID\n");
		return false;
	}

	if (panda_seqid_parse(&data.marker_id, args[0] == '@' ? (args + 1) : args, PANDA_TAG_OPTIONAL) == 0) {
		panda_log_proxy_write_f(logger, "ERR\tBEFORE\tBAD\t%s\n", args);
		return false;
	} else {
		data.state = false;
		*precheck = precheck_func;
		*destroy = free;
		*user_data = PANDA_STRUCT_DUP(&data);
		return true;
	}
}
