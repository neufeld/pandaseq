#include<stdio.h>
#include<stdlib.h>
#include<pandaseq.h>

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
		fprintf(stderr, "ERR\tAFTER\tNO ID\n");
		return false;
	}

	if (panda_seqid_parse(&marker_id, args[0] == '@' ? (args + 1) : args, PANDA_TAG_OPTIONAL) == 0) {
		fprintf(stderr, "ERR\tAFTER\tBAD\t%s\n", args);
		return false;
	} else {
		return true;
	}
}
