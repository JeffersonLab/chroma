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
check_PROGRAMS = t_mesplq t_dslashm t_io t_lwldslash t_formfac t_mesons_w t_conslinop

## DEPENDENCIES = build_chroma_libs

#
# The programs and their dependencies
#
## GTF NOTE: It's getting kind of tedious listing a t_progname_DEPENDENCIES
##   for each target program.  Why doesn't it work to just replace them all
##   with a single 'check_DEPENDENCIES = build_chroma_libs'?  Must rtfm.
##   Also, DEPENDENCIES seems to be a documented feature only
##   in automake-1.7x but seems to do the right thing for automake-1.6x.
##   Worth investigating.
##
## check_DEPENDENCIES = build_chroma_libs

t_mesplq_SOURCES = t_mesplq.cc
t_mesplq_DEPENDENCIES = build_chroma_libs

t_dslashm_SOURCES = t_dslashm.cc
t_dslashm_DEPENDENCIES = build_chroma_libs

t_io_SOURCES = t_io.cc
t_io_DEPENDENCIES = build_chroma_libs

t_lwldslash_SOURCES = t_lwldslash.cc
t_lwldslash_DEPENDENCIES = build_chroma_libs

t_formfac_SOURCES = t_formfac.cc
t_formfac_DEPENDENCIES = build_chroma_libs

t_mesons_w_SOURCES = t_mesons_w.cc
t_mesons_w_DEPENDENCIES = build_chroma_libs

t_conslinop_SOURCES = t_conslinop.cc
t_conslinop_DEPENDENCIES = build_chroma_libs

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd @top_builddir@/lib ; $(MAKE)
