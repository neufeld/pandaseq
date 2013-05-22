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
#include "config.h"
#include <ltdl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif
#include "pandaseq.h"
#include "assembler.h"
#include "buffer.h"

#define STR0(x) #x
#define STR(x) STR0(x)
#define LOGV(code, fmt, ...) snprintf(static_buffer(), BUFFER_SIZE, fmt, __VA_ARGS__); panda_log_proxy_write(assembler->logger, (code), NULL, static_buffer());

#ifdef HAVE_PTHREAD
/* All modules share a single mutex to control reference counts */
static pthread_mutex_t ref_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

struct panda_module {
	volatile size_t refcnt;

	bool (
		*init) (
		char *args);
	PandaCheck check;
	PandaPreCheck precheck;
	PandaDestroy destroy;
	void *user_data;

	lt_dlhandle handle;
	char *name;
	char *args;

	int api;
	char **version;
};

bool module_checkseq(
	PandaAssembler assembler,
	panda_result_seq *sequence) {
	int it;
	for (it = 0; it < assembler->modules_length; it++) {
		PandaModule module = assembler->modules[it];
		if (module->check != NULL && !module->check(sequence, module->user_data)) {
			assembler->rejected[it]++;
			return false;
		}
	}
	return true;
}

bool module_precheckseq(
	PandaAssembler assembler,
	panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length) {
	int it;
	for (it = 0; it < assembler->modules_length; it++) {
		PandaModule module = assembler->modules[it];
		if (module->precheck != NULL && !module->precheck(id, forward, forward_length, reverse, reverse_length, module->user_data)) {
			assembler->rejected[it]++;
			return false;
		}
	}
	return true;
}

void module_init(
	PandaAssembler assembler) {
	int it;
	for (it = 0; it < assembler->modules_length; it++) {
		PandaModule module = assembler->modules[it];
		if (module->handle != NULL) {
			const lt_dlinfo *info = lt_dlgetinfo(module->handle);
			LOGV(PANDA_CODE_MOD_INFO, "%s(%s:%d)\t%s", info == NULL ? "unknown" : info->name, module->version == NULL ? "?" : *(module->version), module->api, module->args);
			assembler->rejected[it] = 0;
		}
	}
}

void panda_assembler_module_stats(
	PandaAssembler assembler) {
	size_t it;
	for (it = 0; it < assembler->modules_length; it++) {
		if (assembler->rejected[it] > 0) {
			LOGV(PANDA_CODE_REJECT_STAT, "%s\t%ld", assembler->modules[it]->name, assembler->rejected[it]);
		}
	}
}

bool panda_assembler_foreach_module(
	PandaAssembler assembler,
	PandaModuleCallback callback,
	void *data) {
	int it;
	for (it = 0; it < assembler->modules_length; it++) {
		if (!callback(assembler, assembler->modules[it], assembler->rejected[it], data)) {
			return false;
		}
	}
	return true;
}

void module_destroy(
	PandaAssembler assembler) {
	int it;
	if (assembler->rejected)
		free(assembler->rejected);
	for (it = 0; it < assembler->modules_length; it++) {
		panda_module_unref(assembler->modules[it]);
	}
	assembler->modules_length = 0;
	free(assembler->modules);
}

static volatile int ltdl_count = 0;
static bool ref_ltdl(
	void) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	if (ltdl_count == 0) {
		const char *path;
		if (lt_dlinit() != 0) {
#ifdef HAVE_PTHREAD
			pthread_mutex_unlock(&ref_lock);
#endif
			return false;
		}
		path = lt_dlgetsearchpath();
		if (path == NULL || strstr(path, STR(PKGLIBDIR)) == NULL) {
			if (lt_dladdsearchdir(STR(PKGLIBDIR)) != 0) {
				(void) lt_dlexit();
				return false;
			}
		}
	}
	ltdl_count++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	return true;
}

static void unref_ltdl(
	void) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	if (--ltdl_count == 0) {
		lt_dlexit();
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
}

static const char path_sep_string[] = { LT_PATHSEP_CHAR, '\0' };

PandaModule panda_module_load(
	const char *path) {
	PandaModule m;
	lt_dlhandle handle;
	PandaCheck check;
	PandaPreCheck precheck;
	size_t name_length;
	int *api;
	char **version;
	char *name;

	name_length = strcspn(path, path_sep_string);
	name = malloc(name_length + 1);
	memcpy(name, path, name_length);
	name[name_length] = '\0';

	if (!ref_ltdl()) {
		free(name);
		return NULL;
	}

	handle = lt_dlopenext(name);
	if (handle == NULL) {
		fprintf(stderr, "Could not open module %s: %s\n", name, lt_dlerror());
		free(name);
		unref_ltdl();
		return NULL;
	}

	api = lt_dlsym(handle, "api");
	if (api == NULL || *api != PANDA_API) {
		lt_dlclose(handle);
		fprintf(stderr, "Invalid API in %s. Are you sure this module was compiled for this version of PANDAseq?\n", name);
		free(name);
		unref_ltdl();
		return NULL;
	}

	*(void **) (&check) = lt_dlsym(handle, "check");
	*(void **) (&precheck) = lt_dlsym(handle, "precheck");
	if (check == NULL && precheck == NULL) {
		lt_dlclose(handle);
		fprintf(stderr, "Could not find check or precheck function in %s\n", name);
		free(name);
		unref_ltdl();
		return NULL;
	}

	m = malloc(sizeof(struct panda_module));
	m->api = *api;
	if (path[name_length] == LT_PATHSEP_CHAR) {
		size_t len = strlen(path + name_length + 1) + 1;
		m->args = malloc(len);
		memcpy(m->args, path + name_length + 1, len);
	} else {
		m->args = NULL;
	}
	m->check = check;
	m->precheck = precheck;
	m->handle = handle;
	m->name = name;
	m->refcnt = 1;
	m->user_data = NULL;
	m->version = lt_dlsym(handle, "version");
	*(void **) (&m->destroy) = lt_dlsym(handle, "destroy");
	*(void **) (&m->init) = lt_dlsym(handle, "init");

	return m;
}

PandaModule panda_module_new(
	const char *name,
	PandaCheck check,
	PandaPreCheck precheck,
	void *user_data,
	PandaDestroy cleanup) {
	PandaModule m;
	if (check == NULL && precheck == NULL)
		return NULL;
	m = malloc(sizeof(struct panda_module));
	m->api = PANDA_API;
	m->args = NULL;
	m->check = check;
	m->destroy = cleanup;
	m->handle = NULL;
	m->init = NULL;
	m->name = malloc(strlen(name) + 1);
	memcpy(m->name, name, strlen(name) + 1);
	m->precheck = precheck;
	m->refcnt = 1;
	m->user_data = user_data;
	m->version = NULL;
	return m;
}

PandaModule panda_module_ref(
	PandaModule module) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	module->refcnt++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	return module;
}

void panda_module_unref(
	PandaModule module) {
	size_t count;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	count = --(module->refcnt);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	if (count == 0) {
		if (module->destroy != NULL)
			module->destroy(module->user_data);
		if (module->name != NULL)
			free(module->name);
		if (module->args != NULL)
			free(module->args);
		if (module->handle != NULL)
			lt_dlclose(module->handle);
		free(module);
		unref_ltdl();
	}
}

const char *panda_module_get_name(
	PandaModule module) {
	return module->name;
}

const char *panda_module_get_args(
	PandaModule module) {
	return module->args;
}

int panda_module_get_api(
	PandaModule module) {
	return module->api;
}

const char *panda_module_get_version(
	PandaModule module) {
	return module->version == NULL ? NULL : *module->version;
}

const char *panda_module_get_description(
	PandaModule module) {
	char **val;
	const lt_dlinfo *info;
	if (module->handle == NULL)
		return NULL;
	info = lt_dlgetinfo(module->handle);
	val = lt_dlsym(module->handle, "desc");
	return val == NULL ? NULL : *val;
}

const char *panda_module_get_usage(
	PandaModule module) {
	char **val;
	const lt_dlinfo *info;
	if (module->handle == NULL)
		return NULL;
	info = lt_dlgetinfo(module->handle);
	val = lt_dlsym(module->handle, "usage");
	return val == NULL ? NULL : *val;
}

bool panda_assembler_add_module(
	PandaAssembler assembler,
	PandaModule module) {
	bool init = true;
	if (module == NULL) {
		return false;
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	if (module->init != NULL) {
		init = module->init(module->args);
		if (init)
			module->init = NULL;
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	if (!init)
		return false;
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&assembler->mutex);
#endif
	if (assembler->modules_length == assembler->modules_size) {
		if (assembler->modules_size == 0) {
			assembler->modules_size = 8;
		} else {
			assembler->modules_size *= 2;
		}
		assembler->modules = realloc(assembler->modules, assembler->modules_size * sizeof(PandaModule));
		assembler->rejected = realloc(assembler->rejected, assembler->modules_size * sizeof(long));
	}
	assembler->rejected[assembler->modules_length] = 0;
	assembler->modules[assembler->modules_length++] = panda_module_ref(module);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&assembler->mutex);
#endif
	return true;
}

size_t panda_assembler_add_modules(
	PandaAssembler assembler,
	PandaModule *modules,
	size_t modules_length) {
	size_t it;
	for (it = 0; it < modules_length; it++) {
		if (!panda_assembler_add_module(assembler, modules[it])) {
			return it;
		}
	}
	return it;
}
