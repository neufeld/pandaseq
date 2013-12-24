#include<stdio.h>
#include<pandaseq-plugin.h>

HELP("Filters sequences based on the contents of a file of ids, one sequence ID per line.", "filter[:file]");

VER_INFO("1.0");

#define BUFF_SIZE 1024
char buffer[BUFF_SIZE];
PandaSet set;

PRECHECK {
	return panda_idset_contains(set, id);
}

INIT {
	FILE *file;
	bool close = false;

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
	while (fgets(buffer, BUFF_SIZE, file) != NULL) {
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
	return true;
}

CLEANUP {
	panda_idset_unref(set);
	return;
}
