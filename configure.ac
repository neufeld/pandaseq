AC_INIT([pandaseq], [2.9], [andre@masella.name])
AC_PREREQ([2.60])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug], [build with debugging output [default=no].]), [DEBUG="$enableval"], [DEBUG="no"])
AC_ARG_ENABLE(threads, AC_HELP_STRING([--disable-threads], [disable thread support (default is autodetect)]))
AC_ARG_WITH(max-len, AC_HELP_STRING([--with-max-len=LEN], [sets the maximum read length]), [MAX_LEN="$withval"], [MAX_LEN="450"])

if ! test "$MAX_LEN" -eq "$MAX_LEN" 2> /dev/null
then
	AC_MSG_ERROR([Bad read length $MAX_LEN])
fi
AC_DEFINE_UNQUOTED([MAX_LEN], $MAX_LEN, [Maximum read length])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR
AM_PROG_CC_C_O
AC_PROG_LIBTOOL
LT_INIT
LT_LIB_M
LT_PROG_RC
AX_CHECK_CFLAGS([-pedantic])
AC_HEADER_STDC

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

if test "x$enable_debug" = xyes; then
	AC_DEFINE(DEBUG, 1, [Enable debugging code])
fi

AC_CANONICAL_BUILD
AC_CANONICAL_HOST

AC_MSG_CHECKING([(non-cross) C compiler for building tools])
if test "$cross_compiling" = "yes"; then
  CC_FOR_BUILD="${CC_FOR_BUILD-gcc}"
else
  CC_FOR_BUILD="${CC_FOR_BUILD-$CC}"
fi
AC_MSG_RESULT([$CC_FOR_BUILD])
AC_ARG_VAR(CC_FOR_BUILD,[C compiler for building tools])

AC_MSG_CHECKING([for module flags])
is_windows=false
MODULE_CFLAGS="$CFLAGS"
MODULE_LDFLAGS="$LDFLAGS -no-undefined -avoid-version -shared"
MODULE_LIBTOOL=libtool
case "$host" in
*-mingw*)
	MODULE_LDFLAGS="$MODULE_LDFLAGS -Wl,-static"
	is_windows=true
	;;
*darwin*|*Darwin*)
	MODULE_LDFLAGS="$MODULE_LDFLAGS -flat_namespace"
	MODULE_LIBTOOL=glibtool
	;;
esac
AC_MSG_RESULT([$MODULE_LDFLAGS])
AC_SUBST(MODULE_LIBTOOL)
AC_SUBST(MODULE_CFLAGS)
AC_SUBST(MODULE_LDFLAGS)

if test "$enable_threads" != no; then
	ACX_PTHREAD
else
	acx_pthread_ok=no
fi
AM_CONDITIONAL([PTHREAD], [test x$acx_pthread_ok = xyes])

AG_CHECK_UNAME_SYSCALL
AC_CHECK_HEADERS_ONCE([sys/param.h])
AC_CHECK_HEADERS([sys/sysctl.h], [], [],
[[#if HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif
]])

PKG_CHECK_MODULES(Z, [ zlib ])
PKG_CHECK_MODULES(CURL, [ libcurl ], [have_curl=true], [have_curl=false])
AM_CONDITIONAL([LIBCURL], [test x$have_curl = xtrue])
LEGACY_CHECK_MODULES(BZ, [bzlib.h], [bz2], [BZ2_bzDecompressInit], [], [], [AC_MSG_ERROR([*** bzip2 required, install bzip2 library])])
LEGACY_CHECK_MODULES(LTDL, [ltdl.h], [ltdl], [lt_dlinit], [], [], [AC_MSG_ERROR([*** ltld required, install libtool library. ])])

AM_CONDITIONAL([IS_WINDOWS], [test x$is_windows = xtrue])

LIB_NAME=pandaseq-$(echo $PACKAGE_VERSION | sed -e 's/\..*$//g')
AC_SUBST(LIB_NAME)

VERSION_MAJOR=$(echo $PACKAGE_VERSION | cut -f 1 -d . )
VERSION_MINOR=$(echo $PACKAGE_VERSION | cut -f 2 -d . )
AC_DEFINE_UNQUOTED(VERSION_MAJOR, ${VERSION_MAJOR}, [Package major version number.])
AC_DEFINE_UNQUOTED(VERSION_MINOR, ${VERSION_MINOR}, [Package minor version number.])

LIB_NAME=pandaseq-$(echo $PACKAGE_VERSION | sed -e 's/\..*$//g')
AC_SUBST(LIB_NAME)

# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html#Updating-version-info
LIB_MAJOR=7
LIB_MINOR=0
LIB_VER=${LIB_MAJOR}:${LIB_MINOR}:0
LIB_URL_VER=0:0:0
AC_DEFINE_UNQUOTED(LIB_MAJOR, ${LIB_MAJOR}, [Library major version number.])
AC_DEFINE_UNQUOTED(LIB_MINOR, ${LIB_MINOR}, [Library minor version number.])
AC_SUBST(LIB_MAJOR)
AC_SUBST(LIB_VER)
AC_SUBST(LIB_URL_VER)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([pandaseq.spec])
AC_CONFIG_FILES(${LIB_NAME}.pc:pc.in, [], [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES(${LIB_NAME}.vapi:vapi.in, [], [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES(${LIB_NAME}-url.deps:deps-url.in, [], [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES(${LIB_NAME}-url.pc:pc-url.in, [], [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES(${LIB_NAME}-url.vapi:vapi-url.in, [], [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES([pandaxs:pandaxs.in], [chmod +x pandaxs])
AC_CONFIG_FILES([build-macos-pkg:build-macos-pkg.in], [chmod +x build-macos-pkg])
AC_OUTPUT
