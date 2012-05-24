#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pandaseq.h>

HELP("Filter out any sequences without a valid index tag.", "validtag:TAG1:TAG2:TAG3");
VER_INFO("1.0");

char **tags;
char *tag_data;
int numtags = 1;
int taglen = 0;

PRECHECK {
	int it;
	const char *tag = id->tag;
	if (tag == NULL)
		return false;

	for (it = 0; it < numtags; it++) {
		if (strncmp(tag, tags[it], taglen) == 0) {
			return true;
		}
	}
	return false;
}

CHECK {
	return true;
}

INIT {
	const char *it = args;
	char *wit;
	char **currtag;
	if (args == NULL) {
		fprintf(stderr, "ERR\tVALTAG\tNOTAGS\n");
		return false;
	}
	while (*it != '\0' && *it != ':') {
		taglen++;
		it++;
	}
	if (taglen == 0) {
		fprintf(stderr, "ERR\tVALTAG\tNOTAGS\n");
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
			numtags++;
			if (*it == ':')
				it++;
			if (currtaglen != taglen) {
				fprintf(stderr, "ERR\tVALTAG\tBADTLEN\t%d != %d %s\n", currtaglen, taglen, it - currtaglen);
				return false;
			}
		}
	}

	tags = malloc(sizeof(char *) * numtags);
	tag_data = malloc(strlen(args));
	memcpy(tag_data, args, strlen(args));
	currtag = tags;
	wit = tag_data;
	*currtag++ = wit;
	while (*wit != '\0') {
		if (*wit == ':') {
			*wit = '\0';
			*currtag++ = ++wit;
		}
		wit++;
	}

	return true;
}

CLEANUP {
	free(tags);
	free(tag_data);
	return;
}
