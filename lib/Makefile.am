# $Id: Makefile.am,v 1.80 2004-01-05 16:51:57 ikuro Exp $
#
# Path for all files used in building the library
# This file was automatically generated

TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@
VPATH := $(TOPDIR)
VPATH += $(TOPDIR)/lib
VPATH += $(TOPDIR)/lib/actions
VPATH += $(TOPDIR)/lib/actions/boson
VPATH += $(TOPDIR)/lib/actions/boson/operator
VPATH += $(TOPDIR)/lib/actions/ferm
VPATH += $(TOPDIR)/lib/actions/ferm/dsdu
VPATH += $(TOPDIR)/lib/actions/ferm/fermacts
VPATH += $(TOPDIR)/lib/actions/ferm/invert
VPATH += $(TOPDIR)/lib/actions/ferm/linop
VPATH += $(TOPDIR)/lib/actions/ferm/qprop
VPATH += $(TOPDIR)/lib/actions/gauge
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
VPATH += $(TOPDIR)/lib/meas/sources
VPATH += $(TOPDIR)/lib/update
VPATH += $(TOPDIR)/lib/update/heatbath
VPATH += $(TOPDIR)/lib/update/molecdyn
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

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/include
endif

lib_LIBRARIES = libchroma.a

#
# Various header files. 
#
HDRS =  chroma.h chromabase.h chromainc.h primitives.h common_declarations.h \
        geom.h

#
# Fermion independent sources
#
libchroma_a_SOURCES = rgauge.cc taproj.cc hotst.cc expm12.cc \
               mesplq.cc polylp.cc reunit.cc unit_check.cc sftmom.cc srcfil.cc \
               szin_io.cc readszin.cc writeszin.cc su3proj.cc sun_proj.cc \
               gaus_smear.cc hyp_smear.cc ape_smear.cc \
               klein_gord.cc laplacian.cc displacement.cc \
               su2extract.cc sunfill.cc invcg1.cc invcg2.cc invcg2_array.cc qprop_io.cc \
               milc_io.cc readmilc.cc writemilc.cc \
               axgauge.cc coulgauge.cc grelax.cc polar_dec.cc rot_colvec.cc \
	       transf.cc

#
# Object files to add if building the Wilson Builds...
# This is OK, since Ns is determined at configure time
#
if BUILD_WILSON_LIBS
libchroma_a_SOURCES += baryon_w.cc mesons_w.cc curcor2_w.cc lwldslash_w.cc formfac_w.cc \
              mesonseqsrc_w.cc barseqsrc_w.cc seqpiontest_w.cc \
              sink_smear2_w.cc barcomp_w.cc ldwfdslash_w.cc walfil_w.cc \
              fermact_qprop.cc prec_fermact_qprop.cc \
              unprec_wilson_fermact_w.cc unprec_wilson_linop_w.cc \
              prec_wilson_fermact_w.cc prec_wilson_linop_w.cc \
              unprec_dwf_fermact_w.cc unprec_dwf_linop_w.cc \
              dwffld_w.cc \
              fermact_qprop_array.cc prec_fermact_qprop_array.cc \
              unprec_dwf_fermact_base_w.cc unprec_dwf_fermact_qprop_array_w.cc \
              unprec_dwf_fermact_array_w.cc unprec_dwf_linop_array_w.cc \
              unprec_ovdwf_fermact_array_w.cc unprec_ovdwf_linop_array_w.cc \
              unprec_ovext_fermact_array_w.cc unprec_ovext_linop_array_w.cc \
              prec_dwf_fermact_base_w.cc prec_dwf_fermact_qprop_array_w.cc \
              prec_dwf_fermact_array_w.cc prec_dwf_linop_array_w.cc \
              readszinqprop_w.cc readszinferm_w.cc writeszinqprop_w.cc \
	      dwf_quarkprop4_w.cc quarkprop4_w.cc \
              overlapbu_linop_w.cc \
              p_src_w.cc d_src_w.cc D_j_w.cc DjDk_w.cc BuildingBlocks_w.cc \
	      zolotarev_coeffs.cc overlap_fermact_base_w.cc lovlapms_w.cc  zolotarev4d_fermact_bj_w.cc
#	      wall_qprop_w.cc
# zolotarev4d_linop_w.cc 
# dslashm_w.cc 
endif

#
# Hooks for optimized code
#
if BUILD_SSE_WILSON_DSLASH
libchroma_a_SOURCES += lwldslash_w_sse.cc
endif


#
# This is the list of files to be used in Staggered Builds
# This is OK since Ns is determined at configure time 
#
if BUILD_STAGGERED_LIBS
libchroma_a_SOURCES += mesphas_s.cc mesphas_follana_s.cc pions_follana_s.cc \
	follana_io_s.cc fermact_qprop.cc walfil_s.cc \
	asqtad_linop_s.cc asq_dsl_s.cc fat7_links_s.cc naik_term_s.cc  \
	asqtad_fermact_s.cc asqtad_mdagm_s.cc staggered_qprop.cc stag_phases_s.cc
endif


# Need this header as part of libchroma_a_sources
nodist_libchroma_a_SOURCES = $(HDRS)
