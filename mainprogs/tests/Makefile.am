#
# VPATH
#
TOPDIR=@top_srcdir@
VPATH = $(TOPDIR)/mainprogs/tests
VPATH += $(TOPDIR)/include

#
# Include Path
#
INCFLAGS=-I$(TOPDIR)/include -I@top_builddir@/include -I$(TOPDIR)/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS = -L@top_builddir@/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @LIBS@

#
# Local Headers
#
HDRS =

#
# add the programs to build in here
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
##noinst_PROGRAMS = t_mesplq t_dslashm t_io
check_PROGRAMS = t_mesplq t_dslashm t_io t_lwldslash

DEPENDENCIES = build_chroma_libs

#
# The programs and their dependencies
#
t_mesplq_SOURCES = t_mesplq.cc
# t_mesplq_DEPENDENCIES = build_chroma_libs

t_dslashm_SOURCES = t_dslashm.cc
# t_dslashm_DEPENDENCIES = build_chroma_libs

t_io_SOURCES = t_io.cc
t_lwldslash_SOURCES = t_lwldslash.cc

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd @top_builddir@/lib ; $(MAKE)
