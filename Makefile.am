NULL =
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq pandaseq-checkid pandaseq-hang
lib_LTLIBRARIES = libpandaseq.la
bin_SCRIPTS = pandaxs
library_includedir=$(includedir)/$(LIB_NAME)
library_include_HEADERS = \
	pandaseq.h \
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
man1_MANS = pandaseq.1 pandaxs.1 pandaseq-checkid.1 pandaseq-hang.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README plugin_sample.c

mktable$(EXEEXT): mktable.c
	@CC_FOR_BUILD@ -std=c99 -pedantic -o $@ $< -lm

pandaseq_CPPFLAGS = $(MAIN_CPPFLAGS) -std=gnu99 -pedantic -g
pandaseq_CFLAGS = $(MAIN_CFLAGS)
pandaseq_LDFLAGS = $(MAIN_LIBS)
pandaseq_SOURCES = main.c
pandaseq_LDADD = libpandaseq.la
pandaseq_checkid_CPPFLAGS = $(MAIN_CPPFLAGS) -std=gnu99 -pedantic -g
pandaseq_checkid_CFLAGS = $(MAIN_CFLAGS)
pandaseq_checkid_LDFLAGS = $(MAIN_LIBS)
pandaseq_checkid_SOURCES = main-parse.c
pandaseq_checkid_LDADD = libpandaseq.la
pandaseq_hang_CPPFLAGS = $(MAIN_CPPFLAGS) -std=gnu99 -pedantic -g
pandaseq_hang_CFLAGS = $(MAIN_CFLAGS)
pandaseq_hang_LDFLAGS = $(MAIN_LIBS)
pandaseq_hang_SOURCES = main-hang.c
pandaseq_hang_LDADD = libpandaseq.la
libpandaseq_la_CPPFLAGS = \
	$(PANDA_CPPFLAGS) \
	$(PTHREAD_CFLAGS) \
	-DPKGLIBDIR=$(pkglibdir)$(LIB_MAJOR) \
	-DPANDA_LIB_COMPILING \
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
	args.c \
	args_array.c \
	args_assembler.c \
	args_fastq.c \
	args_hang.c \
	assembler.c \
	assembler_support.c \
	async.c \
	buffer.c \
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
	output.c \
	proxy.c \
	pool.c \
	seqid.c \
	table.c \
	writer.c \
	$(NULL)
pkglibxdir = $(pkglibdir)$(LIB_MAJOR)
pkglibx_LTLIBRARIES = \
	after.la \
	before.la \
	completely_miss_the_point.la \
	filter.la \
	validtag.la \
	$(NULL)
after_la_CFLAGS = -DPANDASEQ_MODULE=after -g
after_la_SOURCES = plugin_after.c panda_api.c
after_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
after_la_LIBADD = libpandaseq.la
before_la_CFLAGS = -DPANDASEQ_MODULE=before -g
before_la_SOURCES = plugin_before.c panda_api.c
before_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
before_la_LIBADD = libpandaseq.la
completely_miss_the_point_la_CFLAGS = -DPANDASEQ_MODULE=completely_miss_the_point -g
completely_miss_the_point_la_SOURCES = plugin_completely_miss_the_point.c panda_api.c
completely_miss_the_point_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
completely_miss_the_point_la_LIBADD = libpandaseq.la
filter_la_CFLAGS = -DPANDASEQ_MODULE=filter -g
filter_la_SOURCES = plugin_filter.c panda_api.c
filter_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED) $(MODULE_LIBS)
filter_la_LIBADD = libpandaseq.la
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

table.c table.h: mktable$(EXEEXT)
	./mktable$(EXEEXT)

