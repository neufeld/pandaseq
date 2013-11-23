#define _XOPEN_SOURCE 500

#include "config.h"
#include <curl/curl.h>
#include <errno.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#ifdef _WIN32
#        include <windows.h>
typedef void *ucontext_t;
#else
#        include <sys/resource.h>
#        include <ucontext.h>
#endif
#ifdef HAVE_PTHREAD
#        include <pthread.h>
#endif

#include "pandaseq.h"
#include "pandaseq-url.h"

static volatile int curl_count = 0;
#ifdef HAVE_PTHREAD
static pthread_mutex_t ref_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

static bool ref_curl(
	void) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	if (curl_count == 0) {
		const char *path;
		if (curl_global_init(CURL_GLOBAL_ALL) != CURLE_OK) {
#ifdef HAVE_PTHREAD
			pthread_mutex_unlock(&ref_lock);
#endif
			return false;
		}
	}
	curl_count++;
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
	return true;
}

static void unref_curl(
	void) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&ref_lock);
#endif
	if (--curl_count == 0) {
		curl_global_cleanup();
	}
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&ref_lock);
#endif
}

struct curl_data {
	char *url;
	PandaLogProxy logger;
	bool exited;
	bool please_exit;
	ucontext_t read_context;
	ucontext_t curl_context;
	CURL *curl_handle;
	char *curl_memory;
	size_t curl_memory_length;
	char *stack;
};

union data_ints {
	struct curl_data *data;
	int ints[2];
};

static bool read_curl(
	char *buffer,
	size_t buffer_length,
	size_t *read,
	struct curl_data *data) {
	bool buffer_is_small;

	if (data->curl_memory_length == 0) {
		if (data->exited) {
			*read = 0;
			return false;
		}
#ifdef _WIN32
		SwitchToFiber(data->curl_context);
#else
		if (swapcontext(&data->read_context, &data->curl_context) != 0) {
			panda_log_proxy_write_f(data->logger, "%s: swapcontext to curl: %s", data->url, strerror(errno));
			*read = 0;
			return false;
		}
#endif
		if (data->curl_memory_length == 0) {
			*read = 0;
			return true;
		}
	}
	buffer_is_small = buffer_length < data->curl_memory_length;
	*read = buffer_is_small ? buffer_length : data->curl_memory_length;
	memcpy(buffer, data->curl_memory, *read);
	if (buffer_is_small) {
		data->curl_memory += buffer_length;
		data->curl_memory_length -= buffer_length;
	} else {
		data->curl_memory_length = 0;
	}
	return true;
}

static size_t data_ready(
	void *contents,
	size_t size,
	size_t nmemb,
	struct curl_data *data) {
	data->curl_memory = contents;
	data->curl_memory_length = size * nmemb;
#ifdef _WIN32
	SwitchToFiber(data->read_context);
#else
	if (swapcontext(&data->curl_context, &data->read_context) != 0) {
		panda_log_proxy_write_f(data->logger, "%s: swapcontext to read: %s", url, strerror(errno));
		return 0;
	}
#endif
	if (data->please_exit) {
		return 0;
	}
	return size * nmemb;
}

#ifdef _WIN32
static void start_curl(
	struct curl_data *data) {
#else
static void start_curl(
	int i0,
	int i1) {
	union data_ints passed;
	struct curl_data *data;
#endif
	CURLcode res;
#ifndef _WIN32
	passed.ints[0] = i0;
	passed.ints[1] = i1;
	data = passed.data;
#endif
	res = curl_easy_perform(data->curl_handle);

	if (res != CURLE_OK && (res != CURLE_WRITE_ERROR || !data->please_exit)) {
		panda_log_proxy_write_f(data->logger, "%s: %s", data->url, curl_easy_strerror(res));
	}

	curl_easy_cleanup(data->curl_handle);
	data->exited = true;
}

void destroy_curl(
	struct curl_data *data) {
	if (!data->exited) {
		data->please_exit = true;
#ifdef _WIN32
		SwitchToFiber(data->curl_context);
#else
		if (swapcontext(&data->read_context, &data->curl_context) != 0) {
			panda_log_proxy_write_f(data->logger, "%s: swapcontext: %s", data->url, strerror(errno));
		}
#endif
	}
	if (!data->exited)
		curl_easy_cleanup(data->curl_handle);
	free(data->stack);
	free(data->url);
	panda_log_proxy_unref(data->logger);
	free(data);
	unref_curl();
}

#define SET_OPT(opt, value)  if ((res = curl_easy_setopt(curl_handle, (opt), (value))) != CURLE_OK) { panda_log_proxy_write_f(logger, "%s: %s", data->url, curl_easy_strerror(res)); curl_easy_cleanup(curl_handle); unref_curl(); if (data != NULL) free(data); return NULL; }

PandaBufferRead panda_open_url(
	const char *url,
	PandaLogProxy logger,
	void **out_data,
	PandaDestroy *destroy) {
	struct curl_data *data = NULL;
	CURL *curl_handle;
	CURLcode res;
#ifndef _WIN32
	union data_ints passed;
	struct rlimit stack_size;
	if (getrlimit(RLIMIT_STACK, &stack_size) != 0) {
		panda_log_proxy_write_f(data->logger, "%s: getrlimit(RLIMIT_STACK): %s", url, strerror(errno));
		return NULL;
	}
#endif

	if (!ref_curl())
		return NULL;
	curl_handle = curl_easy_init();
	if (curl_handle == NULL) {
		return NULL;
	}
	data = malloc(sizeof(struct curl_data));
	SET_OPT(CURLOPT_URL, url);
	SET_OPT(CURLOPT_WRITEFUNCTION, (curl_write_callback) data_ready);
	SET_OPT(CURLOPT_WRITEDATA, data);
	SET_OPT(CURLOPT_USERAGENT, PACKAGE_STRING);

	data->curl_handle = curl_handle;
	data->url = malloc(strlen(url) + 1);
	memcpy(data->url, url, strlen(url) + 1);
	data->logger = panda_log_proxy_ref(logger);
	data->exited = false;
	data->please_exit = false;
	data->curl_memory_length = 0;
#ifdef _WIN32
	data->read_context = ConvertThreadToFiber(data);
	if (data->read_context == NULL) {
		panda_log_proxy_write_f(logger, "%s: ConvertThreadToFiber error (%d)", url, GetLastError());
		return NULL;
	}

	data->curl_context = CreateFiber(0, (LPFIBER_START_ROUTINE) start_curl, data);

	if (data->curl_context == NULL) {
		panda_log_proxy_write_f(logger, "%s: CreateFiber error (%d)\n", url, GetLastError());
		return NULL;
	}
#else
	passed.data = data;
	getcontext(&data->curl_context);
	data->stack = malloc(stack_size.rlim_cur);
	data->curl_context.uc_stack.ss_sp = data->stack;
	data->curl_context.uc_stack.ss_size = stack_size.rlim_cur;
	data->curl_context.uc_link = &data->read_context;
	makecontext(&data->curl_context, (void (*)(void)) start_curl, 2, passed.ints[0], passed.ints[1]);
#endif
	*out_data = data;
	*destroy = (PandaDestroy) destroy_curl;
	return (PandaBufferRead) read_curl;
}
