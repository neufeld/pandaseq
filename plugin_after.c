#include<pandaseq-plugin.h>

HELP("Include only sequences in the one after the provided sequence", "after:sequenceid");

VER_INFO("1.0");

panda_seq_identifier marker_id;
bool state = false;

PRECHECK {
	if (panda_seqid_equal(&marker_id, id)) {
		state = true;
	}
	return state;
}

INIT {
	if (args == NULL) {
		panda_log_proxy_write_str(logger, "ERR\tAFTER\tNO ID\n");
		return false;
	}

	if (panda_seqid_parse(&marker_id, args[0] == '@' ? (args + 1) : args, PANDA_TAG_OPTIONAL) == 0) {
		panda_log_proxy_write_f(logger, "ERR\tAFTER\tBAD\t%s\n", args);
		return false;
	} else {
		return true;
	}
}
