include Makefile.auto

.PHONY: all gsl lib

all:
	cd $(DIR_GSL) && ./configure && $(MAKE)
	env DIR_SRC_GSL=../$(DIR_GSL) CFLAGS_EXTRA="-DLIBEXTENSION_R=1" $(MAKE) -e -C $(DIR_LIB)
	mv $(DIR_LIB)/*.so ./

distclean:
	$(MAKE) -C $(DIR_LIB) distclean
	$(MAKE) -C $(DIR_GSL) distclean || true
	$(RM) *.so
	
