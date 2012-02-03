/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011  Andre Masella

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
#include<ltdl.h>
#include<stdio.h>
#include<stdlib.h>
#include"pandaseq.h"
#include"plugin.h"

typedef struct module_t module_t;
struct module_t {
	lt_dlhandle handle;
	int (*check) (resultseq * sequence);
	int (*precheck) (seqidentifier *id, char const *forward, char const *reverse);
	void (*destroy) ();
	char *args;
	int (*init)(char* args);
	int api;
	char **version;
	long rejected;
	module_t *next;
};

module_t *modules = NULL;

int module_checkseq(resultseq * sequence)
{
	module_t *current = modules;
	while (current != NULL) {
		if (!current->check(sequence)) {
			current->rejected++;
			return 0;
		}
		current = current->next;
	}
	return 1;
}

int module_precheckseq(seqidentifier *id, char const *forward, char const *reverse)
{
	module_t *current = modules;
	while (current != NULL) {
		if (current->precheck != NULL && !current->precheck(id, forward, reverse)) {
			current->rejected++;
			return 0;
		}
		current = current->next;
	}
	return 1;
}

void module_help()
{
	module_t *current = modules;
	while (current != NULL) {
		const lt_dlinfo *info = lt_dlgetinfo(current->handle);
		char **desc = lt_dlsym(current->handle, "desc");
		char **usage = lt_dlsym(current->handle, "usage");
		fprintf(stderr, "%s: %s\n", info == NULL ? "unknown" : info->name, desc == NULL ? "no description" : *desc);
		if (usage)
			fprintf(stderr, "\tUsage: %s\n", *usage);
		if (info)
			fprintf(stderr, "\tFile: %s\n", info->filename);
		if (current->version)
			fprintf(stderr, "\tVersion: %s\n", *(current->version));
		current = current->next;
	}
}

void module_version()
{
	module_t *current = modules;
	fprintf(stderr, "API: %d\n", PANDA_API);
	while (current != NULL) {
		const lt_dlinfo *info = lt_dlgetinfo(current->handle);
		char **desc = lt_dlsym(current->handle, "desc");
		fprintf(stderr, "%s\n\tFile: %s\n\tVersion: %s\n\tAPI: %d\n", info == NULL ? "unknown" : info->name, info == NULL ? "unknown" : info->filename, current->version == NULL ? "no version" : *(current->version), current->api);
		current = current->next;
	}
}

int module_init()
{
	module_t *current = modules;
	while (current != NULL) {
		const lt_dlinfo *info = lt_dlgetinfo(current->handle);
		fprintf(stderr, "INFO\tMOD\t%s(%s:%d)\t%s\n", info == NULL ? "unknown" : info->name, current->version == NULL ? "?" : *(current->version), current->api, current->args);
		if (current->init != NULL && !current->init(current->args)) {
			return 0;
		}
		current = current->next;
	}
	return 1;

}


void module_cleanup()
{
	while (modules != NULL) {
		module_t *next = modules->next;
		if (modules->destroy != NULL) {
			modules->destroy();
		}
		const lt_dlinfo *info = lt_dlgetinfo(modules->handle);
		fprintf(stderr, "STAT\t%s\t%ld\n", info == NULL ? "unknown" : info->name, modules->rejected);
		lt_dlclose(modules->handle);
		free(modules);
		modules = next;
	}
}

int module_load(char *path)
{
	module_t *m;
	lt_dlhandle handle;
	int (*init) (char *args);
	int (*check) (resultseq * sequence);
	int *api;
	char **version;
	char *args = path;
	while (*args != '\0' && *args != LT_PATHSEP_CHAR) {
		args++;
	}
	if (*args == '\0') {
		args = NULL;
	} else {
		*args = '\0';
		args++;
	}

	handle = lt_dlopenext(path);
	if (handle == NULL) {
		fprintf(stderr, "Could not open module %s: %s\n", path,
			lt_dlerror());
		return 0;
	}

	api = lt_dlsym(handle, "api");
	if (api == NULL || *api > PANDA_API) {
		lt_dlclose(handle);
		fprintf(stderr, "Invalid API in %s. Are you sure this module was compiled for this version of PANDAseq?\n", path);
		return 0;
	}

	*(void **)(&check) = lt_dlsym(handle, "check");
	if (check == NULL) {
		lt_dlclose(handle);
		fprintf(stderr, "Could not find check function in %s\n", path);
		return 0;
	}

	m = malloc(sizeof(module_t));
	m->next = modules;
	m->handle = handle;
	m->check = check;
	m->args = args;
	m->init = init;
	m->api = *api;
	m->rejected = 0;
	*(void **)(&m->init) = lt_dlsym(handle, "init");
	*(void **)(&m->precheck) = lt_dlsym(handle, "precheck");
	*(void **)(&m->destroy) = lt_dlsym(handle, "destroy");
	m->version = lt_dlsym(handle, "version");
	modules = m;
	return 1;
}
