NULL =
ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq
lib_LTLIBRARIES = libpandaseq.la
noinst_PROGRAMS = mktable
bin_SCRIPTS = pandaxs
library_includedir=$(includedir)/$(LIB_NAME)
library_include_HEADERS = pandaseq.h
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(LIB_NAME).pc
vapidir = $(datadir)/vala/vapi
dist_vapi_DATA = $(LIB_NAME).vapi
man1_MANS = pandaseq.1 pandaxs.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README plugins/sample.c
mktable_CFLAGS = -std=gnu99 -pedantic
mktable_SOURCES = mktable.c
pandaseq_CPPFLAGS = $(MAIN_CPPFLAGS) -std=gnu99 -pedantic -DPKGLIBDIR=$(pkglibdir) $(MAIN_CFLAGS)
pandaseq_LDFLAGS = $(MAIN_LDFLAGS) libpandaseq.la
pandaseq_SOURCES = main.c
libpandaseq_la_CPPFLAGS = $(PANDA_CPPFLAGS) $(PTHREAD_CFLAGS)
libpandaseq_la_LDFLAGS = $(PANDA_LDFLAGS) $(PTHREAD_LIBS) -export-symbols-regex '^panda_'
libpandaseq_la_SOURCES = \
	assembler.c \
	assembler_support.c \
	buffer.c \
	fastq.c \
	fileio.c \
	misc.c \
	module.c \
	mux.c \
	nt.c \
	output.c \
	seqid.c \
	table.c \
	$(NULL)
pkglib_LTLIBRARIES = \
	completely_miss_the_point.la \
	validtag.la \
	$(NULL)
validtag_la_CFLAGS = -DPANDASEQ_MODULE=validtag
validtag_la_SOURCES = plugins/validtag.c
validtag_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED)
completely_miss_the_point_la_CFLAGS = -DPANDASEQ_MODULE=completely_miss_the_point
completely_miss_the_point_la_SOURCES = plugins/completely_miss_the_point.c
completely_miss_the_point_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED)
BUILT_SOURCES = table.c table.h
CLEANFILES = table.c table.h $(bin_SCRIPTS)

table.c table.h: mktable
	./mktable
