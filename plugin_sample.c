#include<pandaseq-plugin.h>

/* Adding a validation plugin:
 *
 * To create a validation plugin, copy this file.
 * In your module, you can provide 3 functions, as shown below.
 *
 * When you are ready to compile, use ``pandaxs yourmodule.c'' to compile.
 * You can also link against other libraries
 * (e.g., ``pandaxs yourmodule.c -lglib -I/usr/include/glib-2.0'')
 * If you do this as root, the module will be installed into /usr/lib/pandaseq
 * or a similar location.
 */

/*
 * Provide a description and usage information
 */

HELP("This is a sample module that does nothing", "sample:args");

/*
 * Provide version information
 */
VER_INFO("1.0");

/*
 * If data is needed, store it in a structure.
 */

struct data {
	int some_value;
};

/* Given a sequence, determine if the sequence is valid.
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @sequence: View pandaseq-common.h for more information about panda_result_seq.
 * Returns: true if the sequence should be kept.
 * At least one of this function or precheck is required.
 */
static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	struct data *data) {
	panda_log_proxy_write_str(logger, "INFO\tSAMPLE\tCHECK\n");
	return true;
}

/* Given the forward and reverse reads, determine if the sequence is worth assembling.
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @id: View pandaseq-common.h for more information about panda_seqid.
 * @forward: (array length=forward_length): The forward read, of type panda_qual.
 * @reverse: (array length=reverse_length): The reverse read, of type panda_qual.
 * Returns: true if the reads should be assembled.
 * At least one of this function or check is required.
 */
static bool precheck_func(
	PandaLogProxy logger,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	struct data *data) {
	panda_log_proxy_write_f(logger, "INFO\tSAMPLE\tPRECHECK\n");
	return true;
}

/* Called once upon completion to perform any needed cleanup.
 * This function is optional.
 */
static void destroy_func(
	struct data *data) {
	free(data);
}

/* Called once to initialise the module upon loading. Arguments can be provided
 * to the module upon loading.
 * (e.g., "-C /usr/lib/pandaseq/mynewmodule.so:foo=bar", then args = "foo=bar")
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @args: the argument string provided.
 * Returns: false if there is a failure to initialise.
 */
OPEN {
	struct data data;
	panda_log_proxy_write_f(logger, "INFO\tSAMPLE\tINIT\t%s\n", args);

	*check = (PandaCheck) check_func;
	*precheck = (PandaPreCheck) precheck_func;
	*user_data = PANDA_STRUCT_DUP(&data);
	*destroy = (PandaDestroy) destroy_func;
	return true;
}
