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

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I@top_builddir@/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L@top_builddir@/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
endif

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
check_PROGRAMS  = t_mesplq t_io t_mesons_w t_conslinop t_hypsmear \
                  t_overbu t_ape_smear t_precact

## Moved t_dslashm and t_lwldslash to EXTRA_PROGRAMS because they currently
## do not compile.  This way, they don't break 'make check' but the targets
## 'make t_dslashm' and 'make t_lwldslash' are still there for further
## development

EXTRA_PROGRAMS  = t_dslashm t_lwldslash t_spprod t_lwldslash_sse t_follana_io_s t_formfac  t_follana_pion_s t_dwflinop t_precdwf t_ovlap_bj t_propagator_s t_propagator_w t_gfix t_ovlap5d_bj t_ritz t_ritz_KS t_ritz5d_KS  t_fermion_loop_w t_propagator_fuzz_s t_disc_loop_s t_ov_pbp

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
t_mesplq_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ape_smear_SOURCES = t_ape_smear.cc
t_ape_smear_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_hypsmear_SOURCES = t_hypsmear.cc
t_hypsmear_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_dslashm_SOURCES = t_dslashm.cc
t_dslashm_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_io_SOURCES = t_io.cc
t_io_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_lwldslash_SOURCES = t_lwldslash.cc
t_lwldslash_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_dwflinop_SOURCES = t_dwflinop.cc
t_dwflinop_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_precact_SOURCES = t_precact.cc
t_precact_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_precdwf_SOURCES = t_precdwf.cc
t_precdwf_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_formfac_SOURCES = t_formfac.cc
t_formfac_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_mesons_w_SOURCES = t_mesons_w.cc
t_mesons_w_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_conslinop_SOURCES = t_conslinop.cc
t_conslinop_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_overbu_SOURCES = t_overbu.cc
t_overbu_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_spprod_SOURCES = t_spprod.cc
t_spprod_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_follana_io_s_SOURCES = t_follana_io_s.cc
t_follana_io_s_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_follana_pion_s_SOURCES = t_follana_pion_s.cc
t_follana_pion_s_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_lwldslash_sse_SOURCES = t_lwldslash_sse.cc
t_lwldslash_sse_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ovlap_bj_SOURCES = t_ovlap_bj.cc
t_ovlap_bj_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ovlap5d_bj_SOURCES = t_ovlap5d_bj.cc
t_ovlap5d_bj_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_propagator_s_SOURCES = t_propagator_s.cc
t_propagator_s_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_disc_loop_s_SOURCES = t_disc_loop_s.cc
t_disc_loop_s_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_propagator_fuzz_s_SOURCES = t_propagator_fuzz_s.cc
t_propagator_fuzz_s_DEPENDENCIES = build_chroma_libs rebuild_other_libs


t_propagator_w_SOURCES = t_propagator_w.cc
t_propagator_w_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_fermion_loop_w_SOURCES = t_fermion_loop_w.cc
t_fermion_loop_w_DEPENDENCIES = build_chroma_libs rebuild_other_libs


t_gfix_SOURCES = t_gfix.cc
t_gfix_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ritz_SOURCES = t_ritz.cc
t_ritz_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ritz_KS_SOURCES = t_ritz_KS.cc
t_ritz_KS_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ritz5d_KS_SOURCES = t_ritz5d_KS.cc
t_ritz5d_KS_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ov_pbp_SOURCES = t_ov_pbp.cc
t_ov_pbp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd @top_builddir@/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd @top_builddir@/other_libs ; $(MAKE)
