ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = pandaseq
noinst_PROGRAMS = mktable
bin_SCRIPTS = pandaxs
library_includedir=$(includedir)
library_include_HEADERS = pandaseq.h
man1_MANS = pandaseq.1 pandaxs.1
docdir = $(datadir)/doc/@PACKAGE@
doc_DATA = README plugins/sample.c
mktable_CFLAGS = -std=gnu99 -pedantic
mktable_SOURCES = mktable.c
pandaseq_CPPFLAGS = $(PANDA_CPPFLAGS)
pandaseq_CFLAGS = -std=gnu99 -pedantic -DPKGLIBDIR=$(pkglibdir) $(PANDA_CFLAGS)
pandaseq_LDFLAGS = $(PANDA_LDFLAGS)
pandaseq_SOURCES = main.c table.c plugin.c parser.c
pkglib_LTLIBRARIES = validtag.la
validtag_la_CFLAGS = -DPANDASEQ_MODULE=validtag
validtag_la_SOURCES = plugins/validtag.c
validtag_la_LDFLAGS = -module -shared -avoid-version $(NO_UNDEFINED)
BUILT_SOURCES = table.c table.h
CLEANFILES = table.c table.h $(bin_SCRIPTS)

table.c table.h: mktable
	./mktable
