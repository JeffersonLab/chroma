# $Id: Makefile.am,v 1.1 2002-12-16 06:53:33 edwards Exp $
#
# Path for all files used in building the library
# This file was automatically generated

TOPDIR=@top_srcdir@
VPATH := $(TOPDIR)
VPATH += $(TOPDIR)/lib
VPATH += $(TOPDIR)/lib/actions
VPATH += $(TOPDIR)/lib/actions/construct
VPATH += $(TOPDIR)/lib/actions/dirac
VPATH += $(TOPDIR)/lib/actions/dslash
VPATH += $(TOPDIR)/lib/actions/gauge
VPATH += $(TOPDIR)/lib/actions/invert
VPATH += $(TOPDIR)/lib/actions/linop
VPATH += $(TOPDIR)/lib/actions/qprop
VPATH += $(TOPDIR)/lib/geom
VPATH += $(TOPDIR)/lib/info
VPATH += $(TOPDIR)/lib/io
VPATH += $(TOPDIR)/lib/meas
VPATH += $(TOPDIR)/lib/meas/eig
VPATH += $(TOPDIR)/lib/meas/gfix
VPATH += $(TOPDIR)/lib/meas/glue
VPATH += $(TOPDIR)/lib/meas/hadron
VPATH += $(TOPDIR)/lib/meas/pbp
VPATH += $(TOPDIR)/lib/meas/schrfun
VPATH += $(TOPDIR)/lib/meas/smear
VPATH += $(TOPDIR)/lib/update
VPATH += $(TOPDIR)/lib/update/heatbath
VPATH += $(TOPDIR)/lib/update/molecdyn
VPATH += $(TOPDIR)/lib/update/molecdyn/dsdu
VPATH += $(TOPDIR)/lib/update/molecdyn/est
VPATH += $(TOPDIR)/lib/update/molecdyn/mese
VPATH += $(TOPDIR)/lib/update/molecdyn/noise
VPATH += $(TOPDIR)/lib/update/molecdyn/step
VPATH += $(TOPDIR)/lib/util
VPATH += $(TOPDIR)/lib/util/comm
VPATH += $(TOPDIR)/lib/util/ferm
VPATH += $(TOPDIR)/lib/util/ft
VPATH += $(TOPDIR)/lib/util/gauge
VPATH += $(TOPDIR)/lib/util/misc
VPATH += $(TOPDIR)/lib/util/random
VPATH += $(TOPDIR)/include
VPATH += $(TOPDIR)/mainprogs
VPATH += $(TOPDIR)/mainprogs/main
VPATH += $(TOPDIR)/mainprogs/tests
VPATH += $(TOPDIR)/root

#
# Include Path
#
INCFLAGS = -I$(TOPDIR)/include
#
# Any other compiler flags
# 
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDP_CXXFLAGS@ @QMP_CXXFLAGS@

#
# Various header files. These may be out of date -- I just copied 
# them from Roberts Previous Makefile
#
HDRS =  word.h inner.h reality.h primitive.h outer.h qdptype.h 
	defs.h special.h subset.h params.h multi.h random.h layout.h proto.h \
	scalar_specific.h \
	geom.h

#
# Object files to add if building the Wilson Builds...
# This is OK, since Ns is determined at configure time
#
WILSON_SRCS = dslashm_w.cc baryon_w.cc mesons_w.cc

WILSON_OBJS = $(WILSON_SRCS:.cc=.o)
#
# This is the list of files to be used in Staggered Builds
# This is OK since Ns is determined at configure time 
#
STAGGERED_SRCS = 

STAGGERED_OBJS = $(STAGGERED_SRCS:.cc=.o)

lib_LIBRARIES = libszin.a


#
# These are the permanent source files
# 
nodist_libszin_a_SOURCES = $(HDRS) mesplq.cc reunit.cc geom.cc

#
# This may need to be substituted depending on arch but
# for example on the QCDOC
# 
libszin_a_LIBADD = stand_main.o
libszin_a_DEPENDENCIES = stand_main.o

#
# An AM_CONDITIONAL if statement to decide which fermion type
#
if BUILD_WILSON_LIBS
libszin_a_LIBADD += $(WILSON_OBJS)
libszin_a_DEPENDENCIES += $(WILSON_OBJS)
else 
libszin_a_LIBADD += $(STAGGERED_OBJS)
libszin_a_DEPENDENCIES += $(STAGGERED_OBJS)
endif
