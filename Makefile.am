EXTRA_DIST = makeTAGS.sh autogen.sh
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@ 

# top_builddir = @top_builddir@

SUBDIRS = libtldevel src



TAGS:
	-rm -f ${top_srcdir}/TAGS
	make tags-append	

tags-append:
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/src/*.c
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/src/*.h
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/libtldevel/*.h
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/libtldevel/*.c
