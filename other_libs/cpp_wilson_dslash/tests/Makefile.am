INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

AM_CXXFLAGS=$(INCFLAGS) @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @SSEXX_AUX_LDFLAGS@ @QMT_LDFLAGS@
LDADD = -ldslash @SSEXX_AUX_LIBS@ @LIBS@ @QMT_LIBS@

check_PROGRAMS = t_dslash t_dslash_3d
# t_clover

#if BUILD_SITE_TESTS
#check_PROGRAMS += t_site_ops
#endif

check_PROGRAMS += time_dslash time_dslash_3d
# time_clover
TESTS=
if RUN_TESTS
TESTS += $(check_PROGRAMS)
endif


FRAMEWORK_SRCS =  unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	reunit.h \
	reunit.cc

FRAMEWORK_3D_SRCS = unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	dslash_3d_w.h \
	dslash_3d_w.cc \
	reunit.h \
	reunit.cc

t_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testDslashFull.h \
	testDslashFull.cc \
	t_dslash.cc

#t_clover_SOURCES = $(FRAMEWORK_SRCS) \
#	testClover.h \
#	testClover.cc \
#	t_clover.cc

t_dslash_3d_SOURCES = $(FRAMEWORK_3D_SRCS) \
	testDslash3D.h \
	testDslash3D.cc \
	t_dslash_3d.cc

time_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	timeDslash.h \
	timeDslash.cc \
	time_dslash.cc


time_dslash_3d_SOURCES = $(FRAMEWORK_3D_SRCS) \
	timeDslash3D.h \
	timeDslash3D.cc \
	time_dslash_3d.cc

#time_clover_SOURCES = $(FRAMEWORK_SRCS) \
#	timeClover.h \
#	timeClover.cc \
#	time_clover.cc


#if SCALAR_SINGLE
#t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
#	t_site_ops_scalar.cc \
#	testSiteDslash.h \
#	sse32.h \
#	testSiteDslash32bit.cc
#endif

#if SCALAR_DOUBLE
#t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
#	t_site_ops_scalar.cc \
#	testSiteDslash.h \
#	sse64.h \
#	testSiteDslash64bit.cc
#endif


#if PARSCALAR_SINGLE
#
#check_PROGRAMS += time_decomp
#
#time_decomp_SOURCES = $(FRAMEWORK_SRCS) \
#	timeDecomp.h \
#	timeDecomp.cc \
#	time_decomp.cc

#t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
#	t_site_ops_parscalar.cc \
#	sse32.h \
#	testDecomp.h \
#	testDecompHvv.h \
#	testMvvRecons.h \
#	testRecons.h \
#	testDecomp32bit.cc \
#	testDecompHvv32bit.cc \
#	testMvvRecons32bit.cc \
#	testRecons32bit.cc
#endif

#if PARSCALAR_DOUBLE
#t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
#	t_site_ops_parscalar.cc \
#	sse64.h \
#	testDecomp.h \
#	testDecompHvv.h \
#	testMvvRecons.h \
#	testRecons.h \
#	testDecomp64bit.cc \
#	testDecompHvv64bit.cc \
#	testMvvRecons64bit.cc \
#	testRecons64bit.cc
#endif


t_dslash_DEPENDENCIES = build_lib
# t_clover_DEPENDENCIES = build_lib
#t_dslash_3d_DEPENDENCIES = build_lib
#t_site_ops_DEPENDENCIES = build_lib
#time_dslash_DEPENDENCIES = build_lib
#time_dslash_3d_DEPENDENCIES = build_lib
DEPENDENCIES = build_lib
.PHONY: build_lib
	cd @top_builddir@/lib; $(MAKE)
