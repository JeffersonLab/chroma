#
# VPATH
#
VPATH = @top_srcdir@/mainprogs/tests
VPATH += @top_srcdir@/include

#
# Include Path
#
INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS = -L@top_builddir@/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD = -lszin @QDPXX_LIBS@ @LIBS@

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
noinst_PROGRAMS = t_mesplq t_dslashm t_io
##check_PROGRAMS = t_mesplq t_dslashm t_io

#
# The programs and their dependencies
#
t_mesplq_SOURCES = t_mesplq.cc
t_mesplq_DEPENDENCIES = build_szin_libs

t_dslashm_SOURCES = t_dslashm.cc
t_dslashm_DEPENDENCIES = build_szin_libs

t_io_SOURCES = t_io.cc

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_szin_libs
build_szin_libs:
	cd @top_builddir@/lib ; $(MAKE)
