#include<pandaseq-plugin.h>
#include<stdlib.h>

HELP("Produce statistics on the overlaps examined. Somewhat interesting to see the efficiency of the k-mer table.", "overlap_stat");

VER_INFO("1.0");

struct data {
	PandaWriter writer;
	size_t counts[];
};

static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {

	struct data *data = (struct data *) user_data;

	(void) logger;

	if (sequence->overlaps_examined > 0) {
		data->counts[sequence->overlaps_examined - 1]++;
	}
	return true;
}

static void cleanup(
	struct data *data) {
	size_t it;
	size_t max;

	for (max = PANDA_MAX_LEN - 1; data->counts[max] == 0; max--) ;

	panda_writer_append(data->writer, "STAT\tEXAMINED");
	for (it = 0; it <= max; it++) {
		panda_writer_append(data->writer, " %d", data->counts[it]);
	}
	panda_writer_append_c(data->writer, '\n');
	panda_writer_commit(data->writer);
	panda_writer_unref(data->writer);
	free(data);
}

OPEN {
	struct data *data;
	size_t it;

	if (args != NULL && args[0] != '\0') {
		panda_log_proxy_write_str(logger, "ERR\tOVERLAPSTAT\n");
		return false;
	}
	data = malloc(sizeof(struct data) + sizeof(size_t) * PANDA_MAX_LEN);
	for (it = 0; it < PANDA_MAX_LEN; it++) {
		data->counts[it] = 0;
	}
	data->writer = panda_writer_ref(panda_log_proxy_get_writer(logger));
	*precheck = NULL;
	*check = check_func;
	*destroy = (PandaDestroy) cleanup;
	*user_data = data;
	return true;
}
