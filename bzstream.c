#include "config.h"
#include <bzlib.h>
#include <stdlib.h>
#include "pandaseq.h"
#include "misc.h"

struct bz_stream_data {
	bz_stream strm;
	size_t buffer_length;
	 MANAGED_MEMBER(
		PandaBufferRead,
		source);
	char buffer[];
};

static bool read_stream(
	char *buffer,
	size_t buffer_length,
	size_t *read,
	struct bz_stream_data *data) {

	int ret;
	while (true) {
		if (data->strm.avail_in == 0) {
			size_t avail_in;
			data->strm.next_in = data->buffer;
			if (!data->source(data->buffer, data->buffer_length, &avail_in, data->source_data)) {
				*read = 0;
				return false;
			}
			if (avail_in == 0) {
				*read = 0;
				return true;
			}
			data->strm.avail_in = avail_in;
		}
		data->strm.next_out = buffer;
		data->strm.avail_out = buffer_length;
		ret = BZ2_bzDecompress(&data->strm);

		if (ret == BZ_OK || ret == BZ_STREAM_END) {
			*read = buffer_length - data->strm.avail_out;
			if (*read > 0) {
				return true;
			}
		} else {
			*read = 0;
			return false;
		}
	}
}

static void destroy_stream(
	struct bz_stream_data *data) {
	BZ2_bzDecompressEnd(&data->strm);
	DESTROY_MEMBER(data, source);
	free(data);
}

#define BUFFER_SIZE 1024
PandaBufferRead panda_bz_decompress(
	PandaBufferRead source,
	void *source_data,
	PandaDestroy source_destroy,
	int verbosity,
	void **user_data,
	PandaDestroy *destroy) {

	struct bz_stream_data *data;

	if (source == NULL) {
		*user_data = NULL;
		*destroy = NULL;
		return NULL;
	}

	data = malloc(sizeof(struct bz_stream_data) + BUFFER_SIZE);
	data->buffer_length = BUFFER_SIZE;

	data->strm.bzalloc = NULL;
	data->strm.bzfree = NULL;
	data->strm.opaque = NULL;
	if (BZ2_bzDecompressInit(&data->strm, verbosity, 0) == BZ_OK) {
		data->strm.avail_in = 0;

		data->source = source;
		data->source_data = source_data;
		data->source_destroy = source_destroy;

		*user_data = data;
		*destroy = (PandaDestroy) destroy_stream;
		return (PandaBufferRead) read_stream;
	} else {
		free(data);
		source_destroy(source_data);
		*user_data = NULL;
		*destroy = NULL;
		return NULL;
	}
}
