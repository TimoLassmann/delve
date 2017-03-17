AC_INIT([delve2], [m4_esyscmd_s([git describe --always])])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_CONFIG_FILES([Makefile src/Makefile])

AC_CHECK_HEADERS([stdlib.h])

LIB_TLDEVELDIR="libtldevel"
AC_SUBST(LIB_TLDEVELDIR)

# Call hand's ./configure script recursively.
CFLAGS="-O2 -funroll-loops -Wall" 
ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"

AC_CONFIG_SUBDIRS([libtldevel])

AC_OUTPUT
