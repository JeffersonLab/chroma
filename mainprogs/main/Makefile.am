# Boilerplate stuff
VPATH = @top_srcdir@/mainprogs/tests
VPATH += @top_srcdir@/include

#
# Include Path
#
INCFLAGS = -I@top_srcdir@/include 


#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDP_CXXFLAGS@ @QMP_CXXFLAGS@
AM_LDFLAGS = -L@top_builddir@/lib @QDP_LDFLAGS@ @QMP_LDFLAGS@ \
	      @QMP_COMMS_LDFLAGS@ @LDFLAGS@

LDADD = -lszin @QDP_LIBS@ @QMP_LIBS@ @QMP_COMMS_LIBS@ @LIBS@

#
# add programs in here
# eg:  
#
# noinst_PROGRAMS = propagator
# propagator_SOURCES = propagator.cc
# propagator_DEPENDENCIES = build_szin_libs
#
# The latter rule will always try to rebuild libszin.a when you 
# try to compile propagator

# Always try to rebuild this, just in case it needs updating
.PHONY: build_szin_libs
build_szin_libs:
	cd @top_builddir@/lib ; gmake
