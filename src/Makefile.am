# Build the library in the hand subdirectory first.


LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

MYLIBDIRS = ../libtldevel/libhtsglue.a \
	../libtldevel/librtree.a \
	../libtldevel/libthrpool.a \
	../libtldevel/libtldevel.a \
	../libtldevel/libpwrite.a

MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -lpthread -lhts  


AM_CPPFLAGS=${MYINCDIRS} 

bin_PROGRAMS =  delve bedsim

bedsim_SOURCES = \
bedsim.c 
bedsim_LDADD =  ${MYLIBDIRS}

delve_SOURCES = \
delve.h \
delve.c
delve_LDADD =  ${MYLIBDIRS}

