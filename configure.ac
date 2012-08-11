AC_INIT([pandaseq], [2.0], [andre@masella.name])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_ARG_ENABLE(debug, [  --enable-debug          build with debugging output [default=no].], [DEBUG="$enableval"], [DEBUG="no"])
AC_PROG_LIBTOOL
AX_CHECK_CFLAGS([-std=c99 -pedantic])
AC_HEADER_STDC
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m],[pow])
ACX_PTHREAD

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

if test "x$enable_debug" = xyes; then
	AC_DEFINE(DEBUG, 1, [Enable debugging code])
fi

MODULE_CPPFLAGS="$CPPFLAGS"
MODULE_CFLAGS="$CFLAGS"
MODULE_LDFLAGS="$LDFLAGS"
AC_MSG_CHECKING([for module flags])
case "$host" in
*-mingw*)
	NO_UNDEFINED='-no-undefined'
	MODULE_LDFLAGS="$LDFLAGS -Wl,-static"
	MOD_LIBTOOL=libtool
	;;
*openbsd*|*freebsd*)
	NO_UNDEFINED=''
	MOD_LIBTOOL=libtool
	;;
*darwin*|*Darwin*)
	os_win32=no
	NO_UNDEFINED='-no-undefined'
	MOD_LIBTOOL=glibtool
	;;
*)
	os_win32=no
	NO_UNDEFINED='-no-undefined'
	MOD_LIBTOOL=libtool
	;;
esac
AC_MSG_RESULT([$NO_UNDEFINED])
AC_SUBST(NO_UNDEFINED)
AC_SUBST(MOD_LIBTOOL)

AC_CHECK_LIB(ltdl, lt_dlinit,,AC_MSG_ERROR(Cannot link with libtool libs! (libltdl)))
MAIN_CPPFLAGS="$CPPFLAGS"
MAIN_CFLAGS="$CFLAGS"
MAIN_LDFLAGS="$LDFLAGS"
AC_SUBST(MAIN_CPPFLAGS)
AC_SUBST(MAIN_CFLAGS)
AC_SUBST(MAIN_LDFLAGS)

AC_CHECK_HEADER([zlib.h], [], [AC_MSG_ERROR([*** zlib.h is required, install zlib header files])])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([*** gzopen is required, install zlib library files])])
AC_CHECK_HEADER([bzlib.h], [], [AC_MSG_ERROR([*** bzlib.h is required, install bzip2 header files])])
AC_CHECK_LIB([bz2], [BZ2_bzDecompressInit], [], [AC_MSG_ERROR([*** BZ2_bzDecompressInit is required, install bzip2 library files])])
PANDA_CPPFLAGS="$CPPFLAGS"
PANDA_CFLAGS="$CFLAGS"
PANDA_CFLAGS="$LDFLAGS"

AC_SUBST(PANDA_CPPFLAGS)
AC_SUBST(PANDA_CFLAGS)
AC_SUBST(PANDA_LDFLAGS)
CPPFLAGS="$MODULE_CPPFLAGS"
CFLAGS="$MODULE_CFLAGS"
LDFLAGS="$MODULE_LDFLAGS"

LIB_NAME=pandaseq-2
AC_SUBST(LIB_NAME)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([${LIB_NAME}.pc:${LIB_NAME}.pc.in], , [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES([pandaxs:pandaxs.in], [chmod +x pandaxs])
AC_OUTPUT
