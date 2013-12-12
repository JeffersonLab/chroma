#
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

#
# Include Path
#
INCFLAGS=-I$(TOPDIR)/mainprogs/main -I$(TOPDIR)/lib -I$(BUILDDIR)/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @QDPXX_CXXFLAGS@ @BAGEL_WILSON_DSLASH_CXXFLAGS@ @GMP_CXXFLAGS@ @QOP_MG_CXXFLAGS@ @CXXFLAGS@
AM_LDFLAGS = -L$(BUILDDIR)/lib  @QDPXX_LDFLAGS@ @BAGEL_WILSON_DSLASH_LDFLAGS@ @GMP_LDFLAGS@ @LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @BAGEL_WILSON_DSLASH_LIBS@ @GMP_LIBS@  @LIBS@ 

if BUILD_MDWF
AM_CXXFLAGS += @MDWF_CXXFLAGS@
AM_LDFLAGS  += @MDWF_LDFLAGS@
LDADD     += @MDWF_LIBS@
endif


if BUILD_QUDA
AM_CXXFLAGS += @QUDA_CXXFLAGS@
AM_LDFLAGS += @QUDA_LDFLAGS@
LDADD	+= @QUDA_LIBS@
endif

if BUILD_QOP_MG
AM_LDFLAGS += -L@top_builddir@/other_libs/wilsonmg/lib @QOP_MG_LDFLAGS@
LDADD += @QOP_MG_LIBS@
endif

if BUILD_ASQTAD_LEVEL3_INVERTER_WRAPPER

AM_CXXFLAGS += -I@QDPC_DIR@/qmp/include
AM_CXXFLAGS += -I@QDPC_DIR@/qio/include
AM_CXXFLAGS += -I@QDPC_DIR@/qla/include
AM_CXXFLAGS += -I@QDPC_DIR@/qdp/include
AM_CXXFLAGS += -I@QDPC_DIR@/qopqdp/include

LDADD  += -L@QDPC_DIR@/qopqdp/lib -lqopqdp
LDADD  += -L@QDPC_DIR@/qdp/lib
LDADD  +=   -lqdp_d3  -lqdp_df3 -lqdp_f3 -lqdp_df -lqdp_d -lqdp_f -lqdp_int -lqdp_common
LDADD += -L@QDPC_DIR@/qmp/lib -lqmp 
LDADD += -L@QDPC_DIR@/qio/lib  -lqio -llime

LDADD  += -L@QDPC_DIR@/qla/lib
LDADD  += -lqla_c99 -lqla_cmath -lqla_d3 -lqla_d -lqla_df3 -lqla_df -lqla_dq3 -lqla_dq -lqla_f3 -lqla_f -lqla_int -lqla_q3 -lqla_q -lqla_random

endif

AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@QDP_LAPACK_DIR@/lib
LDADD += -lqdp-lapack

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
endif

if BUILD_CPP_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@CPP_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@CPP_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@CPP_DSLASH_DIR@/lib
LDADD += -ldslash @QDPXX_LIBS@
endif

if BUILD_CG_DWF
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@CG_DWF_DIR@
AM_CXXFLAGS += -I$(BUILDDIR)/other_libs/@CG_DWF_DIR@
AM_LDFLAGS  += -L$(BUILDDIR)/other_libs/@CG_DWF_DIR@
LDADD += -lcg-dwf
endif

if BUILD_BAGEL_CLOVER_APPLY
AM_CXXFLAGS += @BAGEL_CLOVER_CXXFLAGS@
AM_LDFLAGS += @BAGEL_CLOVER_LDFLAGS@
LDADD += @BAGEL_CLOVER_LIBS@
endif

if BUILD_QMT
AM_CXXFLAGS += @QMT_CXXFLAGS@
AM_LDFLAGS  += @QMT_LDFLAGS@
LDADD += @QMT_LIBS@
endif


#
# Local Headers
#
HDRS =

## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
#noinst_PROGRAMS = propagator
##
## BJ (21/02/05)
## Moved these to be bin_PROGRAMS so that they are installed. Easier
## than for users to have to always make them by hand
##
bin_PROGRAMS = chroma  \
        purgaug \
        cfgtransf \
        hmc spectrum_s

# Wilson specific programs
check_PROGRAMS = collect_propcomp qpropgfix qproptrev qpropqio qproptransf wallformfac \
        make_meson_ops make_baryon_ops baryon_ops_red_dil


# Staggered specific programs
#check_PROGRAMS += 



## These are programs still under development. They won't compile in
## the check_PROGRAMS line, so we put them here to keep from breaking
## make check
EXTRA_PROGRAMS = qpropadd


# Wilson specific programs
# Temporarily moved here while fixing codes.
EXTRA_PROGRAMS += propagator_comp multi_propagator multi_propagator_comp \
       collect_multi_propcomp 

# Staggered specific programs
EXTRA_PROGRAMS += spectrum_s


#
# The program and its dependencies
#
chroma_SOURCES= chroma.cc
propagator_comp_SOURCES= propagator_comp.cc
multi_propagator_comp_SOURCES = multi_propagator_comp.cc
collect_propcomp_SOURCES = collect_propcomp.cc
collect_multi_propcomp_SOURCES = collect_multi_propcomp.cc
spectrum_s_SOURCES= spectrum_s.cc
hmc_SOURCES= hmc.cc
purgaug_SOURCES= purgaug.cc
wallformfac_SOURCES= wallformfac.cc
qproptransf_SOURCES= qproptransf.cc
qpropqio_SOURCES= qpropqio.cc
qpropadd_SOURCES= qpropadd.cc
qpropgfix_SOURCES= qpropgfix.cc
qproptrev_SOURCES= qproptrev.cc
cfgtransf_SOURCES= cfgtransf.cc

make_baryon_ops_SOURCES= make_baryon_ops.cc
make_meson_ops_SOURCES= make_meson_ops.cc
baryon_ops_red_dil_SOURCES= baryon_ops_red_dil.cc

#mainparam_io_SOURCES= mainparam_io.cc

#
# The latter rule will always try to rebuild libchroma.a when you 
# try to compile example

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target

DEPENDENCIES = build_chroma_libs rebuild_other_libs
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

.PHONY: build_chroma_libs
build_chroma_libs:
	cd $(BUILDDIR)/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd $(BUILDDIR)/other_libs ; $(MAKE)
