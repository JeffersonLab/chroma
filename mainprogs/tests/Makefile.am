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

noinst_PROGRAMS = t_mesplq t_dslashm t_io

t_mesplq_SOURCES = t_mesplq.cc
t_mesplq_DEPENDENCIES = build_szin_libs

t_dslashm_SOURCES = t_dslashm.cc
t_dslashm_DEPENDENCIES = build_szin_libs

t_io_SOURCES = t_io.cc

#
# Always try to rebuild this, just in case it needs updating
.PHONY: build_szin_libs
build_szin_libs:
	cd @top_builddir@/lib ; $(MAKE)
