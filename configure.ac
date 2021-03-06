AC_INIT([delve2], [m4_esyscmd_s([git describe --always])])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_PROG_INSTALL


tldevel_check_f=$ac_abs_confdir/tldevel/README

AC_CHECK_FILE([$tldevel_check_f], [],[have_tldevel=no])


if test "x${have_tldevel}" = xno; then
   echo "------------------------------------------"
   echo " Unable to find pthreads on this system.  "
   echo " Building a single-threaded version.      "
   echo "------------------------------------------"
   AC_MSG_ERROR([ASFSAF])
fi



# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])

htslib_support=ok
AC_CHECK_LIB([hts], [main],[], [htslib_support=missing])


if test $htslib_support != ok; then
  AC_MSG_ERROR([htslib not found

This packages uses the htslib library. Please install and try again!


https://github.com/samtools/htslib


FAILED.  This error must be resolved in order to build this package successfully.])
fi

AC_CHECK_HEADERS([stdlib.h string.h])

LIB_TLDEVELDIR="tldevel"
AC_SUBST(LIB_TLDEVELDIR)

AC_TYPE_UINT8_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T

AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([arc4random arc4random_uniform])


AC_ARG_ENABLE(debugging,[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging="no")

case $enable_debugging in
yes)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(DEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(DEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(DEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac

if test "$enable_debugging" != "no"; then
        AC_DEFINE(DEBUG,1,[Defines debugging .])
        CFLAGS="-ggdb -Wall -m64 -std=c99 $MEM $SIGNAL32"
else
        CFLAGS="-O2 -funroll-loops -Wall -std=c99 $MEM $SIGNAL32"
        DEBUG=0
fi

ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"



AC_ARG_ENABLE(valgrind-tests,
    [--enable-valgrind-tests run self tests under valgrind],
    opt_valgrind_tests="$enableval", opt_valgrind_tests="no")
    

OPTS="-q --trace-children=no --error-exitcode=1 --leak-check=yes"
  # Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
        AC_DEFINE(DEBUG,1,[Defines debugging .])   
        CFLAGS="-pg -ggdb -pedantic -Wall -m64 -std=gnu11"
	AC_CHECK_PROGS(VALGRIND, valgrind)
	VALGRIND="$VALGRIND $OPTS"
else
	VALGRIND=
fi
AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)



ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"


AC_CONFIG_FILES([Makefile
                 src/Makefile])
		 
AC_CONFIG_SUBDIRS([tldevel])

AC_OUTPUT
