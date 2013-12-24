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

/* Given a sequence, determine if the sequence is valid.
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @sequence: View pandaseq-common.h for more information about panda_result_seq.
 * Returns: true if the sequence should be kept.
 * At least one of this function or PRECHECK is required.
 */
CHECK {
	panda_log_proxy_write_str(logger, "INFO\tSAMPLE\tCHECK\n");
	return true;
}

/* Given the forward and reverse reads, determine if the sequence is worth assembling.
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @id: View pandaseq-common.h for more information about panda_seqid.
 * @forward: (array length=forward_length): The forward read, of type panda_qual.
 * @reverse: (array length=reverse_length): The reverse read, of type panda_qual.
 * Returns: true if the reads should be assembled.
 * At least one of this function or CHECK is required.
 */
PRECHECK {
	panda_log_proxy_write_f(logger, "INFO\tSAMPLE\tPRECHECK\n");
	return true;
}

/* Called once to initialise the module upon loading. Arguments can be provided
 * to the module upon loading.
 * (e.g., "-C /usr/lib/pandaseq/mynewmodule.so:foo=bar", then args = "foo=bar")
 * @logger: View pandaseq-log.h for more information about PandaLogProxy.
 * @args: the argument string provided.
 * Returns: false if there is a failure to initialise.
 * This function is optional.
 */
INIT {
	panda_log_proxy_write_f(logger, "INFO\tSAMPLE\tINIT\t%s\n", args);
	return false;
}

/* Called once upon completion to perform any needed cleanup.
 * This function is optional.
 */
CLEANUP {
	return;
}
