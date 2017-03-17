
EXTRA_DIST = htslib-1.4.tar

all:  htslib

htslib:
	tar -xvf $(top_srcdir)/libinc/htslib-1.4
	cd htslib-1.4 && ./configure &&  $(MAKE)


clean-local:
	rm -rf htslib-4/


