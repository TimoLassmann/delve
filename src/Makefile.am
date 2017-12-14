# Build the library in the hand subdirectory first.


LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

MYLIBDIRS = ../tldevel/libhtsglue.a \
	../tldevel/librtree.a \
	../tldevel/libthrpool.a \
	../tldevel/libtldevel.a \
	../tldevel/libpwrite.a

MYINCDIRS = -I${top_builddir}/${LIB_TLDEVELDIR} \
            -I${top_srcdir}/${LIB_TLDEVELDIR}

LIBS = -lm -lpthread -lhts  


AM_CPPFLAGS=${MYINCDIRS} 

bin_PROGRAMS =  delve bedsim run_delve dseed

run_delve_SOURCES = \
run_delve.c
run_delve_LDADD =  ${MYLIBDIRS}

dseed_SOURCES = \
seed.c
dseed_LDADD =  ${MYLIBDIRS}



bedsim_SOURCES = \
bedsim.c 
bedsim_LDADD =  ${MYLIBDIRS}

delve_SOURCES = \
delve_struct.h \
delve.h \
delve.c \
shared_data.h \
shared_data.c \
genome_sequences.h \
genome_sequences.c \
dyn_prog.h \
dyn_prog.c \
thread_data.h \
thread_data.c \
hmm.h \
hmm.c \
genome_rtree.h \
genome_rtree.c \
genome_priors.h \
genome_priors.c \
io.h \
io.c

delve_LDADD =  ${MYLIBDIRS}
