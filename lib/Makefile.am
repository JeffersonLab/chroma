# $Id: Makefile.am,v 1.21 2003-04-09 17:13:56 edwards Exp $
#
# Path for all files used in building the library
# This file was automatically generated

TOPDIR=@top_srcdir@
VPATH := $(TOPDIR)
VPATH += $(TOPDIR)/lib
VPATH += $(TOPDIR)/lib/actions
VPATH += $(TOPDIR)/lib/actions/boson
VPATH += $(TOPDIR)/lib/actions/boson/operator
VPATH += $(TOPDIR)/lib/actions/ferm
VPATH += $(TOPDIR)/lib/actions/ferm/fermacts
VPATH += $(TOPDIR)/lib/actions/ferm/invert
VPATH += $(TOPDIR)/lib/actions/ferm/linop
VPATH += $(TOPDIR)/lib/actions/ferm/qprop
VPATH += $(TOPDIR)/lib/actions/gauge
VPATH += $(TOPDIR)/lib/geom
VPATH += $(TOPDIR)/lib/info
VPATH += $(TOPDIR)/lib/io
VPATH += ${TOPDIR}/lib/io/ioqprop
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
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@

lib_LIBRARIES = libchroma.a

#
# Various header files. 
#
HDRS =  chroma.h chromabase.h chromainc.h primitives.h defs.h common_declarations.h geom.h

#
# Fermion independent sources
#
GENERIC_SRCS = $(HDRS) gaus_smear.cc taproj.cc hotst.cc expm12.cc \
               klein_gord.cc mesplq.cc reunit.cc sftmom.cc srcfil.cc \
               readszin.cc su3proj.cc sun_proj.cc hyp_smear.cc \
               su2extract.cc sunfill.cc invcg2.cc

#
# Object files to add if building the Wilson Builds...
# This is OK, since Ns is determined at configure time
#
WILSON_SRCS = baryon_w.cc mesons_w.cc lwldslash_w.cc\
	      formfac_w.cc sink_smear2_w.cc barcomp_w.cc transf_w.cc \
              unprec_wilson_fermact_w.cc unprec_wilson_linop_w.cc unprec_qprop_w.cc \
              readszinqprop_w.cc 
#              zolotarev4d_fermact_w.cc zolotarev4d_linop_w.cc lovlapms_w.cc

# dslashm_w.cc 

#
# This is the list of files to be used in Staggered Builds
# This is OK since Ns is determined at configure time 
#
STAGGERED_SRCS = 

#
# An AM_CONDITIONAL if statement to decide which fermion type
#
if BUILD_WILSON_LIBS
  libchroma_a_SOURCES = $(GENERIC_SRCS) $(WILSON_SRCS)
  EXTRA_libchroma_a_SOURCES = $(STAGGERED_SRCS)
endif

if BUILD_STAGGERED_LIBS
  libchroma_a_SOURCES = $(GENERIC_SRCS) $(STAGGERED_SRCS)
  EXTRA_libchroma_a_SOURCES = $(WILSON_SRCS)
endif
