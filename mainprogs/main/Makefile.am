#
# VPATH
#
TOPDIR=@top_srcdir@
VPATH = @top_srcdir@/mainprogs/main
VPATH += @top_srcdir@/include

#
# Include Path
#
INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include -I${TOPDIR}/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS = -L@top_builddir@/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @LIBS@

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L@top_builddir@/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
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
check_PROGRAMS = propagator spectrum_w qqq_w bar3ptfn make_source ExampleBuildingBlocks

#
# The program and its dependencies
#
propagator_SOURCES= propagator.cc
propagator_DEPENDENCIES= build_chroma_libs

spectrum_w_SOURCES= spectrum_w.cc
spectrum_w_DEPENDENCIES= build_chroma_libs

qqq_w_SOURCES= qqq_w.cc
qqq_w_DEPENDENCIES= build_chroma_libs

bar3ptfn_SOURCES= bar3ptfn.cc
bar3ptfn_DEPENDENCIES= build_chroma_libs

make_source_SOURCES= make_source.cc
make_source_DEPENDENCIES= build_chroma_libs

ExampleBuildingBlocks_SOURCES= ExampleBuildingBlocks.cc
ExampleBuildingBlocks_DEPENDENCIES= build_chroma_libs


#
# The latter rule will always try to rebuild libchroma.a when you 
# try to compile example

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd @top_builddir@/lib ; $(MAKE)
