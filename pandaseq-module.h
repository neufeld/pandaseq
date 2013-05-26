/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2012  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef _PANDASEQ_MODULE_H
#        define _PANDASEQ_MODULE_H
#        ifdef __cplusplus
#                define EXTERN_C_BEGIN  extern "C" {
#                define EXTERN_C_END    }
#        else
#                define EXTERN_C_BEGIN
#                define EXTERN_C_END
#        endif
#        include <pandaseq-common.h>
EXTERN_C_BEGIN
/**
 * The current module API version of the running library
 */
int panda_api_version(
	void);

/* === Constructors === */

/**
 * Create a module given sequence checking parameters.
 *
 * @name: the name of the module, for user interaction
 * @check: (closure user_data): the function to be run after assembly
 * @precheck: (closure user_data): a function to be run before assembly
 * @user_data: (transfer full): the context data for the functions. The user is responsible for managing the memory associated with user_data, but the cleanup function will always be called to do so.
 * @cleanup: (closure user_data): a function to be called when this module is garbage collected
 */
PandaModule panda_module_new(
	const char *name,
	PandaCheck check,
	PandaPreCheck precheck,
	void *user_data,
	PandaDestroy cleanup);

/**
 * Load a module from a string containg the module name and arguments.
 *
 * @path: the name or path to a module separated by LT_PATHSEP_CHAR and any arguments to the initialisation function of that module
 */
PandaModule panda_module_load(
	const char *path);

/* === Methods === */

/**
 * Increase the reference count on a module.
 */
PandaModule panda_module_ref(
	PandaModule module);

/**
 * Decrease the reference count on a module.
 * @module: (transfer full): the module to release.
 */
void panda_module_unref(
	PandaModule module);

/* === Getter and Setters === */

/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules. Modules constructed by panda_module_new will always return PANDA_API.
 */
int panda_module_get_api(
	PandaModule module);

/**
 * Get the arguments passed on loading of a module of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_args(
	PandaModule module);

/**
 * Get the description of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The description help text.
 */
const char *panda_module_get_description(
	PandaModule module);

/**
 * Get the name of a module.
 *
 * Returns: (transfer none): the module's name
 */
const char *panda_module_get_name(
	PandaModule module);

/**
 * Get the usage information (i.e., help text) of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_usage(
	PandaModule module);

/**
 * Get the version of a module.
 *
 * This is only appropriate for loaded modules.
 * Returns: (transfer none) (allow-none): The usage help text.
 */
const char *panda_module_get_version(
	PandaModule module);
EXTERN_C_END
#endif
