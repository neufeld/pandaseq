#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pandaseq.h>

HELP("Filter out any sequences without a valid index tag.", "validtag:TAG1:TAG2:TAG3");
VER_INFO("1.0");

char **tags;
int numtags = 1;
int taglen = 0;

PRECHECK {
	int it;
	char *tag = id->tag;
	if (tag == NULL)
		return 0;

	for (it = 0; it < numtags; it++) {
		if (strncmp(tag, tags[it], taglen) == 0) {
			return 1;
		}
	}
	return 0;
}

CHECK {
	return 1;
}

INIT {
	char *it = args;
	char **currtag;
	if (args == NULL) {
		fprintf(stderr, "ERR\tVALTAG\tNOTAGS\n");
		return 0;
	}
	while (*it != '\0' && *it != ':') {
		taglen++;
		it++;
	}
	if (taglen == 0) {
		fprintf(stderr, "ERR\tVALTAG\tNOTAGS\n");
		return 0;
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
				fprintf(stderr, "ERR\tVALTAG\tBADTLEN\n");
				return 0;
			}
		}
	}

	tags = malloc(sizeof(char *) * numtags);
	currtag = tags;
	it = args;
	*currtag++ = it;
	while (*it != '\0') {
		if (*it == ':') {
			*it = '\0';
			*currtag++ = ++it;
		}
		it++;
	}

	return 1;
}

CLEANUP {
	free(tags);
	return;
}
