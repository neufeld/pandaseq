#include<stdlib.h>
#include<string.h>
#include<pandaseq-plugin.h>

HELP("Filter out any sequences without a valid index tag.", "validtag:TAG1:TAG2:TAG3");
VER_INFO("1.0");

struct data {
	char **tags;
	char *tag_data;
	int numtags;
	int taglen;
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

	int it;
	const char *tag = id->tag;
	if (tag == NULL)
		return false;

	for (it = 0; it < data->numtags; it++) {
		if (strncmp(tag, data->tags[it], data->taglen) == 0) {
			return true;
		}
	}
	return false;
}

static void destroy_func(
	struct data *data) {
	free(data->tags);
	free(data->tag_data);
	free(data);
}

OPEN {
	struct data data;
	const char *it = args;
	char *wit;
	char **currtag;

	data.numtags = 1;
	data.taglen = 0;

	if (args == NULL) {
		panda_log_proxy_write_f(logger, "ERR\tVALTAG\tNOTAGS\n");
		return false;
	}
	while (*it != '\0' && *it != ':') {
		data.taglen++;
		it++;
	}
	if (data.taglen == 0) {
		panda_log_proxy_write_f(logger, "ERR\tVALTAG\tNOTAGS\n");
		return false;
	}

	if (*it != '\0') {
		it++;
		while (*it != '\0') {
			int currtaglen = 0;
			while (*it != ':' && *it != '\0') {
				it++;
				currtaglen++;
			}
			data.numtags++;
			if (*it == ':')
				it++;
			if (currtaglen != data.taglen) {
				panda_log_proxy_write_f(logger, "ERR\tVALTAG\tBADTLEN\t%d != %d %s\n", currtaglen, data.taglen, it - currtaglen);
				return false;
			}
		}
	}

	data.tags = malloc(sizeof(char *) * data.numtags);
	data.tag_data = malloc(strlen(args) + 1);
	memcpy(data.tag_data, args, strlen(args) + 1);
	currtag = data.tags;
	wit = data.tag_data;
	*currtag++ = wit;
	while (*wit != '\0') {
		if (*wit == ':') {
			*wit = '\0';
			*currtag++ = ++wit;
		}
		wit++;
	}

	*precheck = precheck_func;
	*user_data = PANDA_STRUCT_DUP(&data);
	*destroy = (PandaDestroy) destroy_func;
	return true;
}
