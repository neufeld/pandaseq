#include<stdio.h>
#include<pandaseq-plugin.h>

HELP("Filters sequences based on the contents of a file of ids, one sequence ID per line.", "filter:file");

VER_INFO("1.0");

static bool precheck_func(
	PandaLogProxy logger,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	void *user_data) {
	return panda_idset_contains((PandaSet) user_data, id);
}

OPEN {
	char buffer[1024];
	bool close = false;
	FILE *file;
	PandaSet set;

	if (args == NULL || *args == '\0') {
		file = stdin;
	} else {
		file = fopen(args, "r");
		if (file == NULL) {
			panda_log_proxy_perror(logger, args);
			return false;
		}
	}
	set = panda_idset_new();
	while (fgets(buffer, sizeof(buffer), file) != NULL) {
		int it;
		for (it = 0; buffer[it] != '\n'; it++) ;
		buffer[it] = '\0';

		if (!panda_idset_add_str(set, buffer[0] == '@' ? (buffer + 1) : buffer, PANDA_TAG_OPTIONAL, NULL, NULL)) {
			panda_log_proxy_write_f(logger, "ERR\tFILTER\tBAD\t%s\n", buffer);
			if (close)
				fclose(file);
			return false;
		}
	}
	if (ferror(file)) {
		panda_log_proxy_perror(logger, args);
		if (close)
			fclose(file);
		return false;
	}
	if (close)
		fclose(file);
	*precheck = precheck_func;
	*user_data = set;
	*destroy = (PandaDestroy) panda_idset_unref;
	return true;
}
