## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha

SUBDIRS = gsl utils sys err const block complex eigen vector matrix permutation sort linalg ieee-utils cblas blas specfunc rng randist statistics histogram cdf poly test

SUBLIBS = block/libgslblock.la linalg/libgsllinalg.la complex/libgslcomplex.la eigen/libgsleigen.la blas/libgslblas.la err/libgslerr.la histogram/libgslhistogram.la ieee-utils/libgslieeeutils.la matrix/libgslmatrix.la permutation/libgslpermutation.la randist/libgslrandist.la rng/libgslrng.la sort/libgslsort.la specfunc/libgslspecfunc.la statistics/libgslstatistics.la sys/libgslsys.la utils/libutils.la vector/libgslvector.la cdf/libgslcdf.la poly/libgslpoly.la test/libgsltest.la

pkginclude_HEADERS = gsl_math.h gsl_pow_int.h gsl_nan.h gsl_machine.h gsl_mode.h gsl_precision.h gsl_types.h gsl_version.h gsl_minmax.h gsl_inline.h

bin_SCRIPTS = gsl-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= gsl.pc

CLEANFILES = gsl.pc gsl-config

EXTRA_DIST = autogen.sh gsl-config.in gsl.pc.in configure.ac THANKS BUGS gsl.spec.in gsl.m4 test_gsl_histogram.sh pkgconfig.test

lib_LTLIBRARIES = libgsl.la
libgsl_la_SOURCES = version.c
libgsl_la_LIBADD = $(SUBLIBS)
libgsl_la_LDFLAGS = -version-info $(GSL_LT_VERSION)
noinst_HEADERS = templates_on.h templates_off.h build.h

MINGW32_HOST = @MINGW32_HOST@
if MINGW32_HOST
libgsl_la_LIBADD += cblas/libgslcblas.la
libgsl_la_LDFLAGS += -no-undefined
endif

m4datadir = $(datadir)/aclocal
m4data_DATA = gsl.m4

bin_PROGRAMS = gsl-randist gsl-histogram

gsl_randist_SOURCES = gsl-randist.c
gsl_randist_LDADD = libgsl.la cblas/libgslcblas.la

gsl_histogram_SOURCES = gsl-histogram.c
gsl_histogram_LDADD = libgsl.la cblas/libgslcblas.la

check_SCRIPTS = test_gsl_histogram.sh pkgconfig.test
TESTS = test_gsl_histogram.sh pkgconfig.test

#bin_PROGRAMS = main dummy
#dummy_SOURCES = version.c
#dummy_LDADD = $(SUBLIBS)
#main_SOURCES = version.c env.c
#main_LDADD = libgsl.la

edit = $(SED) \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@GSL_CFLAGS[@]|$(GSL_CFLAGS)|g' \
	-e 's|@GSL_LIBM[@]|$(GSL_LIBM)|g' \
	-e 's|@GSL_LIBS[@]|$(GSL_LIBS)|g' \
	-e 's|@LIBS[@]|$(LIBS)|g' \
	-e 's|@VERSION[@]|$(VERSION)|g'

gsl-config gsl.pc: Makefile 
	@rm -f $@ $@.tmp
	@$(edit) '$(srcdir)/$@.in' >>$@.tmp
	@mv $@.tmp $@
	@echo creating $@

gsl-config: $(srcdir)/gsl-config.in
gsl.pc: $(srcdir)/gsl.pc.in

