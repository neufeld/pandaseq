#include<pandaseq-plugin.h>

HELP("Drops empty (zero-length) output sequences.", "empty");
VER_INFO("1.0");

static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void* data) {
	(void) logger;
	(void) data;

	return sequence->sequence_length > 0;
}

OPEN {
	*precheck = NULL;
	*check = (PandaCheck) check_func;
	*destroy = NULL;
	*user_data = NULL;

	if (args != NULL && *args != '\0') {
		panda_log_proxy_write_f(logger, "No arguments allowed to empty filter.");
		return false;
	}
	return true;
}
