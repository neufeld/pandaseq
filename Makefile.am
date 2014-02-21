NULL =
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq pandaseq-checkid pandaseq-diff pandaseq-hang
lib_LTLIBRARIES = libpandaseq.la
bin_SCRIPTS = pandaxs
library_includedir=$(includedir)/$(LIB_NAME)
library_include_HEADERS = \
	pandaseq.h \
	pandaseq-algorithm.h \
	pandaseq-args.h \
	pandaseq-assembler.h \
	pandaseq-common.h \
	pandaseq-iter.h \
	pandaseq-linebuf.h \
	pandaseq-log.h \
	pandaseq-module.h \
	pandaseq-mux.h \
	pandaseq-nt.h \
	pandaseq-plugin.h \
	pandaseq-seqid.h \
	pandaseq-set.h \
	pandaseq-writer.h \
	panda_api.c \
	$(NULL)
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(LIB_NAME).pc
vapidir = $(datadir)/vala/vapi
dist_vapi_DATA = $(LIB_NAME).vapi
man1_MANS = pandaseq.1 pandaxs.1 pandaseq-checkid.1 pandaseq-diff.1 pandaseq-hang.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README plugin_sample.c

mktable$(EXEEXT): mktable.c tablebuilder.c
	@CC_FOR_BUILD@ -std=c99 -pedantic -o $@ $^ -lm

pandaseq_CPPFLAGS = $(MAIN_CPPFLAGS) -std=c99 -pedantic -g
pandaseq_CFLAGS = $(MAIN_CFLAGS)
pandaseq_LDFLAGS = $(MAIN_LIBS)
pandaseq_SOURCES = main.c
pandaseq_LDADD = libpandaseq.la
pandaseq_checkid_CPPFLAGS = $(MAIN_CPPFLAGS) -std=c99 -pedantic -g
pandaseq_checkid_CFLAGS = $(MAIN_CFLAGS)
pandaseq_checkid_LDFLAGS = $(MAIN_LIBS)
pandaseq_checkid_SOURCES = main-parse.c
pandaseq_checkid_LDADD = libpandaseq.la
pandaseq_diff_CPPFLAGS = $(MAIN_CPPFLAGS) -std=c99 -pedantic -g
pandaseq_diff_CFLAGS = $(MAIN_CFLAGS)
pandaseq_diff_LDFLAGS = $(MAIN_LIBS)
pandaseq_diff_SOURCES = main-diff.c
pandaseq_diff_LDADD = libpandaseq.la
pandaseq_hang_CPPFLAGS = $(MAIN_CPPFLAGS) -std=c99 -pedantic -g
pandaseq_hang_CFLAGS = $(MAIN_CFLAGS)
pandaseq_hang_LDFLAGS = $(MAIN_LIBS)
pandaseq_hang_SOURCES = main-hang.c
pandaseq_hang_LDADD = libpandaseq.la
libpandaseq_la_CPPFLAGS = \
	$(PANDA_CPPFLAGS) \
	$(PTHREAD_CFLAGS) \
	-DPKGLIBDIR=$(pkglibdir)$(LIB_MAJOR) \
	-DPANDA_LIB_COMPILING \
	-std=c99 \
	-g \
	$(NULL)
libpandaseq_la_LDFLAGS = \
	$(PANDA_LIBS) \
	$(PTHREAD_LIBS) \
	-export-symbols-regex '^panda_' \
	-version-info $(LIB_VER) \
	$(NO_UNDEFINED) \
	$(NULL)
libpandaseq_la_SOURCES = \
	algo.c \
	algo_pear.c \
	algo_simple_bayes.c \
	algo_rdp_mle.c \
	args.c \
	args_array.c \
	args_assembler.c \
	args_fastq.c \
	args_hang.c \
	assembler.c \
	assembler_support.c \
	async.c \
	buffer.c \
	diff.c \
	fastq.c \
	fileio.c \
	hang.c \
	idset.c \
	iter.c \
	linebuf.c \
	misc.c \
	module.c \
	mux.c \
	nt.c \
	offset.c \
	output.c \
	proxy.c \
	pool.c \
	seqid.c \
	table.c \
	writer.c \
	$(NULL)
if IS_WINDOWS
libpandaseq_la_SOURCES += lib.rc
.rc.lo:
	$(LIBTOOL) --tag=RC --mode=compile $(RC) $< -o $@
endif
pkglibxdir = $(pkglibdir)$(LIB_MAJOR)
pkglibx_LTLIBRARIES = \
	after.la \
	before.la \
	completely_miss_the_point.la \
	filter.la \
	min_phred.la \
	min_overlapbits.la \
	min_readqscore.la \
	pear_test.la \
	validtag.la \
	$(NULL)
after_la_CPPFLAGS = -std=c99
after_la_CFLAGS = -DPANDASEQ_MODULE=after -g
after_la_SOURCES = plugin_after.c panda_api.c
after_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
after_la_LIBADD = libpandaseq.la
before_la_CPPFLAGS = -std=c99
before_la_CFLAGS = -DPANDASEQ_MODULE=before -g
before_la_SOURCES = plugin_before.c panda_api.c
before_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
before_la_LIBADD = libpandaseq.la
completely_miss_the_point_la_CPPFLAGS = -std=c99
completely_miss_the_point_la_CFLAGS = -DPANDASEQ_MODULE=completely_miss_the_point -g
completely_miss_the_point_la_SOURCES = plugin_completely_miss_the_point.c panda_api.c
completely_miss_the_point_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
completely_miss_the_point_la_LIBADD = libpandaseq.la
filter_la_CPPFLAGS = -std=c99
filter_la_CFLAGS = -DPANDASEQ_MODULE=filter -g
filter_la_SOURCES = plugin_filter.c panda_api.c
filter_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
filter_la_LIBADD = libpandaseq.la
min_phred_la_CPPFLAGS = -std=c99
min_phred_la_CFLAGS = -DPANDASEQ_MODULE=min_phred -g
min_phred_la_SOURCES = plugin_min_phred.c panda_api.c
min_phred_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
min_phred_la_LIBADD = libpandaseq.la
min_overlapbits_la_CPPFLAGS = -std=c99
min_overlapbits_la_CFLAGS = -DPANDASEQ_MODULE=min_overlapbits -g
min_overlapbits_la_SOURCES = plugin_min_overlapbits.c panda_api.c
min_overlapbits_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
min_overlapbits_la_LIBADD = libpandaseq.la
min_readqscore_la_CPPFLAGS = -std=c99
min_readqscore_la_CFLAGS = -DPANDASEQ_MODULE=min_readqscore -g
min_readqscore_la_SOURCES = plugin_min_readqscore.c panda_api.c
min_readqscore_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
min_readqscore_la_LIBADD = libpandaseq.la
pear_test_la_CPPFLAGS = -std=c99
pear_test_la_CFLAGS = -DPANDASEQ_MODULE=pear_test -g
pear_test_la_SOURCES = plugin_pear_test.c panda_api.c
pear_test_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS) -lm
pear_test_la_LIBADD = libpandaseq.la
validtag_la_CPPFLAGS = -std=c99
validtag_la_CFLAGS = -DPANDASEQ_MODULE=validtag -g
validtag_la_SOURCES = plugin_validtag.c panda_api.c
validtag_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
validtag_la_LIBADD = libpandaseq.la
BUILT_SOURCES = table.c table.h
CLEANFILES = \
	table.c \
	table.h \
	$(bin_SCRIPTS) \
	mktable$(EXEEXT) \
	$(NULL)
if LIBCURL
lib_LTLIBRARIES += libpandaseq-url.la
library_include_HEADERS += pandaseq-url.h
pkgconfig_DATA += $(LIB_NAME)-url.pc
dist_vapi_DATA += $(LIB_NAME)-url.vapi
libpandaseq_url_la_CPPFLAGS = \
	$(PANDA_CPPFLAGS) \
	$(PTHREAD_CFLAGS) \
	$(LIBCURL_CPPFLAGS) \
	-I. \
	-g \
	-std=c99 \
	$(NULL)
libpandaseq_url_la_LIBADD = libpandaseq.la
libpandaseq_url_la_LDFLAGS = \
	$(PANDA_LIBS) \
	$(PTHREAD_LIBS) \
	$(LIBCURL) \
	-export-symbols-regex '^panda_' \
	-version-info $(LIB_URL_VER) \
	$(NO_UNDEFINED) \
	$(NULL)
libpandaseq_url_la_SOURCES = \
	bzstream.c \
	curl_reader.c \
	$(NULL)
endif

table.c table.h: mktable$(EXEEXT)
	./mktable$(EXEEXT)

