# Build the library in the hand subdirectory first.


LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 



MYLIBDIRS = -L${top_builddir}/${LIB_TLDEVELDIR} 

MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -ltldevel -lhtsglue -lhts


AM_CPPFLAGS=${MYINCDIRS} 

bin_PROGRAMS = arm delve

arm_SOURCES = arm.c
arm_LDADD = ${MYLIBDIRS}


delve_SOURCES = \
delve.h \
delve.c

delve_LDADD =  ${MYLIBDIRS}

