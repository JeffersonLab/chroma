#
# VPATH
#
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@
VPATH = @top_srcdir@/mainprogs/main

#
# Include Path
#
INCFLAGS=-I$(TOPDIR)/mainprogs/main -I$(TOPDIR)/lib -I$(BUILDDIR)/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS = -L$(BUILDDIR)/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @LIBS@

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
endif

if BUILD_SSE_DWF_CG
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DWF_CG_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DWF_CG_DIR@/lib
LDADD += -lsse_dwf_cg 
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
bin_PROGRAMS = propagator \
        collect_propcomp \
        seqsource \
	make_source \
        ExampleBuildingBlocks hypsmear cfgtransf sink_smearing qpropgfix \
        qproptrev apesmear stoutsmear hmc hypsmear3d hypsmear4d

# Wilson specific programs
bin_PROGRAMS += spectrum_w qqq_w bar3ptfn qproptransf wallformfac mres

# Staggered specific programs
#check_PROGRAMS += 



## These are programs still under development. They won't compile in
## the check_PROGRAMS line, so we put them here to keep from breaking
## make check
EXTRA_PROGRAMS = qpropadd 

# Wilson specific programs
# Temporarily moved here while fixing codes.
EXTRA_PROGRAMS += propagator_comp multi_propagator multi_propagator_comp \
       collect_multi_propcomp propagator_hack

# Staggered specific programs
#EXTRA_PROGRAMS += 


#
# The program and its dependencies
#
propagator_SOURCES= propagator.cc
propagator_DEPENDENCIES= build_chroma_libs rebuild_other_libs

propagator_comp_SOURCES= propagator_comp.cc
propagator_comp_DEPENDENCIES= build_chroma_libs rebuild_other_libs

propagator_hack_SOURCES= propagator_hack.cc
propagator_hack_DEPENDENCIES= build_chroma_libs rebuild_other_libs

seqsource_SOURCES= seqsource.cc
seqsource_DEPENDENCIES= build_chroma_libs rebuild_other_libs

multi_propagator_SOURCES= multi_propagator.cc
multi_propagator_DEPENDENCIES= build_chroma_libs rebuild_other_libs

multi_propagator_comp_SOURCES = multi_propagator_comp.cc
multi_propagator_comp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

collect_propcomp_SOURCES = collect_propcomp.cc
collect_propcomp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

collect_multi_propcomp_SOURCES = collect_multi_propcomp.cc
collect_multi_propcomp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

spectrum_w_SOURCES= spectrum_w.cc
spectrum_w_DEPENDENCIES= build_chroma_libs rebuild_other_libs

hmc_SOURCES= hmc.cc
hmc_DEPENDENCIES= build_chroma_libs rebuild_other_libs

mres_SOURCES= mres.cc
mres_DEPENDENCIES= build_chroma_libs rebuild_other_libs

qqq_w_SOURCES= qqq_w.cc
qqq_w_DEPENDENCIES= build_chroma_libs rebuild_other_libs

bar3ptfn_SOURCES= bar3ptfn.cc
bar3ptfn_DEPENDENCIES= build_chroma_libs rebuild_other_libs

wallformfac_SOURCES= wallformfac.cc
wallformfac_DEPENDENCIES= build_chroma_libs rebuild_other_libs

make_source_SOURCES= make_source.cc
make_source_DEPENDENCIES= build_chroma_libs rebuild_other_libs

apesmear_SOURCES= apesmear.cc
apesmear_DEPENDENCIES= build_chroma_libs rebuild_other_libs

hypsmear_SOURCES= hypsmear.cc
hypsmear_DEPENDENCIES= build_chroma_libs rebuild_other_libs

hypsmear3d_SOURCES= hypsmear3d.cc
hypsmear3d_DEPENDENCIES= build_chroma_libs rebuild_other_libs

hypsmear4d_SOURCES= hypsmear4d.cc
hypsmear4d_DEPENDENCIES= build_chroma_libs rebuild_other_libs

stoutsmear_SOURCES= stoutsmear.cc
stoutsmear_DEPENDENCIES= build_chroma_libs rebuild_other_libs

qproptransf_SOURCES= qproptransf.cc
qproptransf_DEPENDENCIES= build_chroma_libs rebuild_other_libs

qpropadd_SOURCES= qpropadd.cc
qpropadd_DEPENDENCIES= build_chroma_libs rebuild_other_libs

qpropgfix_SOURCES= qpropgfix.cc
qpropgfix_DEPENDENCIES= build_chroma_libs rebuild_other_libs

qproptrev_SOURCES= qproptrev.cc
qproptrev_DEPENDENCIES= build_chroma_libs rebuild_other_libs

ExampleBuildingBlocks_SOURCES= ExampleBuildingBlocks.cc
ExampleBuildingBlocks_DEPENDENCIES= build_chroma_libs rebuild_other_libs

cfgtransf_SOURCES= cfgtransf.cc
cfgtransf_DEPENDENCIES= build_chroma_libs rebuild_other_libs

sink_smearing_SOURCES= sink_smearing.cc
sink_smearing_DEPENDENCIES= build_chroma_libs rebuild_other_libs

#mainparam_io_SOURCES= mainparam_io.cc
#sink_smearing_DEPENDENCIES= build_chroma_libs rebuild_other_libs

#
# The latter rule will always try to rebuild libchroma.a when you 
# try to compile example

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd $(BUILDDIR)/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd $(BUILDDIR)/other_libs ; $(MAKE)
