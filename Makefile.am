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
	pandaseq-nt.h \
	pandaseq-plugin.h \
	pandaseq-seqid.h \
	pandaseq-set.h \
	pandaseq-writer.h \
	panda_api.c \
	$(NULL)
if PTHREAD
library_include_HEADERS += pandaseq-mux.h
endif
EXTRA_DIST = \
	algo.h \
	assembler.h \
	buffer.h \
	buffer.list \
	config.h \
	misc.h \
	mktable.c \
	module.h \
	nt.h \
	pandaseq-tablebuilder.h \
	prob.h \
	README.md \
	tablebuilder.c \
	table.h \
	$(doc_DATA) \
	$(man1_MANS) \
	$(NULL)
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(LIB_NAME).pc
vapidir = $(datadir)/vala/vapi
dist_vapi_DATA = $(LIB_NAME).vapi
man1_MANS = pandaseq.1 pandaxs.1 pandaseq-checkid.1 pandaseq-diff.1 pandaseq-hang.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README plugin_sample.c
TESTS = \
	./check_parser \
	$(NULL)
check_PROGRAMS = \
	check_parser \
	$(NULL)

mktable$(EXEEXT): mktable.c tablebuilder.c
	@CC_FOR_BUILD@ $(COMMON_CPPFLAGS) -o $@ $^ -lm

COMMON_CPPFLAGS = \
	-std=c99 \
	-g \
  -Wall -Wextra -Wformat \
	$(NULL)

check_parser_CPPFLAGS = $(COMMON_CPPFLAGS)
check_parser_SOURCES = check_parser.c
check_parser_LDADD = libpandaseq.la
pandaseq_CPPFLAGS = $(COMMON_CPPFLAGS)
pandaseq_SOURCES = main.c
pandaseq_LDADD = libpandaseq.la
pandaseq_checkid_CPPFLAGS = $(COMMON_CPPFLAGS)
pandaseq_checkid_SOURCES = main-parse.c
pandaseq_checkid_LDADD = libpandaseq.la
pandaseq_diff_CPPFLAGS = $(COMMON_CPPFLAGS)
pandaseq_diff_SOURCES = main-diff.c
pandaseq_diff_LDADD = libpandaseq.la
pandaseq_hang_CPPFLAGS = $(COMMON_CPPFLAGS)
pandaseq_hang_SOURCES = main-hang.c
pandaseq_hang_LDADD = libpandaseq.la
libpandaseq_la_CPPFLAGS = \
	$(BZ_CFLAGS) \
	$(LTDL_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	$(Z_CFLAGS) \
	-DPKGLIBDIR=$(pkglibdir)$(LIB_MAJOR) \
	-DPANDA_LIB_COMPILING \
	$(COMMON_CPPFLAGS) \
	$(NULL)
libpandaseq_la_LDFLAGS = \
	$(BZ_LIBS) \
	$(LIBM) \
	$(LTDL_LIBS) \
	$(PTHREAD_LIBS) \
	$(Z_LIBS) \
	-export-symbols-regex '^panda_' \
	-version-info $(LIB_VER) \
	-no-undefined \
	$(NULL)
libpandaseq_la_SOURCES = \
	algo.c \
	algo_ea_util.c \
	algo_flash.c \
	algo_pear.c \
	algo_rdp_mle.c \
	algo_simple_bayes.c \
	algo_stitch.c \
	algo_uparse.c \
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
	nt.c \
	offset.c \
	output.c \
	proxy.c \
	pool.c \
	seqid.c \
	table.c \
	writer.c \
	$(NULL)
if PTHREAD
libpandaseq_la_SOURCES += mux.c
endif

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
	empty.la \
	filter.la \
	min_phred.la \
	min_overlapbits.la \
	other_primer.la \
	overlap_stat.la \
	pear_test.la \
	validtag.la \
	$(NULL)
after_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=after \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
after_la_SOURCES = plugin_after.c panda_api.c
after_la_LDFLAGS = -module $(MODULE_LDFLAGS)
after_la_LIBADD = libpandaseq.la
before_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=before \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
before_la_SOURCES = plugin_before.c panda_api.c
before_la_LDFLAGS = -module $(MODULE_LDFLAGS)
before_la_LIBADD = libpandaseq.la
completely_miss_the_point_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=completely_miss_the_point \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
completely_miss_the_point_la_SOURCES = plugin_completely_miss_the_point.c panda_api.c
completely_miss_the_point_la_LDFLAGS = -module $(MODULE_LDFLAGS)
completely_miss_the_point_la_LIBADD = libpandaseq.la
empty_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=empty \
	$(COMMON_CPPFLAGS)
empty_la_SOURCES = plugin_empty.c panda_api.c
empty_la_LDFLAGS = -module $(MODULE_LDFLAGS)
empty_la_LIBADD = libpandaseq.la
filter_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=filter \
	$(COMMON_CPPFLAGS)
filter_la_SOURCES = plugin_filter.c panda_api.c
filter_la_LDFLAGS = -module $(MODULE_LDFLAGS)
filter_la_LIBADD = libpandaseq.la
min_phred_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=min_phred \
	$(COMMON_CPPFLAGS)
min_phred_la_SOURCES = plugin_min_phred.c panda_api.c
min_phred_la_LDFLAGS = -module $(MODULE_LDFLAGS)
min_phred_la_LIBADD = libpandaseq.la
min_overlapbits_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=min_overlapbits \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
min_overlapbits_la_SOURCES = plugin_min_overlapbits.c panda_api.c
min_overlapbits_la_LDFLAGS = -module $(MODULE_LDFLAGS)
min_overlapbits_la_LIBADD = libpandaseq.la
other_primer_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=other_primer \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
other_primer_la_SOURCES = plugin_other_primer.c panda_api.c
other_primer_la_LDFLAGS = -module $(MODULE_LDFLAGS) $(LIBM)
other_primer_la_LIBADD = libpandaseq.la
overlap_stat_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=overlap_stat \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
overlap_stat_la_SOURCES = plugin_overlap_stat.c panda_api.c
overlap_stat_la_LDFLAGS = -module $(MODULE_LDFLAGS) $(LIBM)
overlap_stat_la_LIBADD = libpandaseq.la
pear_test_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=pear_test \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
pear_test_la_SOURCES = plugin_pear_test.c panda_api.c
pear_test_la_LDFLAGS = -module $(MODULE_LDFLAGS) $(LIBM)
pear_test_la_LIBADD = libpandaseq.la
validtag_la_CPPFLAGS = \
	-DPANDASEQ_MODULE=validtag \
	$(MODULE_CFLAGS) \
	$(COMMON_CPPFLAGS)
validtag_la_SOURCES = plugin_validtag.c panda_api.c
validtag_la_LDFLAGS = -module $(MODULE_LDFLAGS)
validtag_la_LIBADD = libpandaseq.la
BUILT_SOURCES = table.c table.h
CLEANFILES = \
	table.c \
	table.h \
	$(bin_SCRIPTS) \
	mktable$(EXEEXT) \
	*.[ch]~ \
	$(NULL)
if LIBCURL
lib_LTLIBRARIES += libpandaseq-url.la
library_include_HEADERS += pandaseq-url.h
pkgconfig_DATA += $(LIB_NAME)-url.pc
dist_vapi_DATA += $(LIB_NAME)-url.vapi
libpandaseq_url_la_CPPFLAGS = \
	$(BZ_CFLAGS) \
	$(CURL_CFLAGS) \
	$(PTHREAD_CFLAGS) \
	-I. \
	-g \
	-std=c99 \
	$(NULL)
libpandaseq_url_la_LIBADD = libpandaseq.la
libpandaseq_url_la_LDFLAGS = \
	$(BZ_LIBS) \
	$(CURL_LIBS) \
	$(PTHREAD_LIBS) \
	-export-symbols-regex '^panda_' \
	-version-info $(LIB_URL_VER) \
	-no-undefined \
	$(NULL)
libpandaseq_url_la_SOURCES = \
	bzstream.c \
	curl_reader.c \
	$(NULL)
endif

table.c table.h: mktable$(EXEEXT)
	./mktable$(EXEEXT)

