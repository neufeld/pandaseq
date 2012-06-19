#include<stdio.h>
#include<stdlib.h>
#include<pandaseq.h>

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

HELP("This is a sample module that does nothing", "No configurable options.");

/*
 * Provide version information
 */
VER_INFO("1.0");

/* Given a sequence, determine if the sequence is valid.
 * Arguments: resultseq* sequence
 *                       View pandaseq.h for more information about resultseq.
 * Return true if the sequence should be kept.
 * This function is required. */
CHECK
{
	fprintf(stderr, "INFO\tSAMPLE\tCHECK\n");
	return true;
}

/* Given the forward and reverse reads, determine if the sequence is worth assembling.
 * Arguments: seqidentifier* id
 *                       View pandaseq.h for more information about id.
 *            char* forward
 *                       The forward sequence.
 *            char* reverse
 *                       The reverse sequence.
 * Return true if the reads should be assembled.
 * This function is optional. */
PRECHECK
{
	fprintf(stderr, "INFO\tSAMPLE\tPRECHECK\n");
	return true;
}

/* Called once to initialise the module upon loading. Arguments can be provided
 * to the module upon loading.
 * (e.g., "-C /usr/lib/pandaseq/mynewmodule.so:foo=bar", then args = "foo=bar")
 * Arguments: char* args
 * Returns false if there is a failure to initialise.
 * This function is optional. */
INIT
{
	fprintf(stderr, "INFO\tSAMPLE\tINIT\t%s\n", args);
	return false;
}

/* Called once upon completion to perform any needed cleanup.
 * This function is optional. */
CLEANUP
{
	fprintf(stderr, "INFO\tSAMPLE\tDESTROY\n");
	return;
}
