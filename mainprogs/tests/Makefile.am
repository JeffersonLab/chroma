#
# VPATH
#
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@
VPATH = $(TOPDIR)/mainprogs/tests

#
# Include Path
#
INCFLAGS=-I$(TOPDIR)/lib -I$(BUILDDIR)/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @CXXFLAGS@ @QDPXX_CXXFLAGS@
AM_LDFLAGS = -L$(BUILDDIR)/lib @LDFLAGS@ @QDPXX_LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @LIBS@

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
endif

if BUILD_SSE_DWF_CG
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DWF_CG_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DWF_CG_DIR@/lib
LDADD += -lsse_dwf_cg 
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

EXTRA_PROGRAMS  = t_dslashm t_lwldslash t_spprod t_lwldslash_sse t_follana_io_s t_formfac  t_follana_pion_s t_dwflinop t_precdwf t_ovlap_bj t_propagator_s t_propagator_w t_gfix t_ovlap5d_bj t_ritz t_ritz_KS t_ritz5d_KS  t_fermion_loop_w t_propagator_fuzz_s t_disc_loop_s t_ov_pbp t_lwldslash_pab t_propagator_nrqcd t_dwflocality t_invert3_precwilson t_invert4_precwilson t_read_eigen t_fuzwilp t_ovlap_double_pass t_g5eps_bj t_sumr t_msumr t_invrelcg t_bicgstab t_wilslp t_rel_gmresr t_invborici t_hamsys t_hamsys_ferm t_precact_sse t_hmc_pg t_neflinop t_prec_contfrac t_mres_4d t_prec_nef

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

t_fuzwilp_SOURCES = t_fuzwilp.cc
t_fuzwilp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_wilslp_SOURCES = t_wilslp.cc
t_wilslp_DEPENDENCIES = build_chroma_libs rebuild_other_libs

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

t_lwldslash_pab_SOURCES = t_lwldslash_pab.cc
t_lwldslash_pab_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ovlap_bj_SOURCES = t_ovlap_bj.cc
t_ovlap_bj_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_ovlap_double_pass_SOURCES = t_ovlap_double_pass.cc
t_ovlap_double_pass_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_g5eps_bj_SOURCES = t_g5eps_bj.cc
t_g5eps_bj_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_sumr_SOURCES = t_sumr.cc
t_sumr_DEPENDENCIES= build_chroma_libs rebuild_other_libs

t_rel_gmresr_SOURCES = t_rel_gmresr.cc
t_rel_gmresr_DEPENDENCIES= build_chroma_libs rebuild_other_libs

t_msumr_SOURCES = t_msumr.cc
t_msumr_DEPENDENCIES= build_chroma_libs rebuild_other_libs

t_invrelcg_SOURCES = t_invrelcg.cc
t_invrelcg_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_neflinop_SOURCES = t_neflinop.cc
t_neflinop_DEPENDENCIES = build_chroma_libs rebuild_other_libs

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


t_propagator_nrqcd_SOURCES = t_propagator_nrqcd.cc
t_propagator_nrqcd_DEPENDENCIES = build_chroma_libs rebuild_other_libs


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

t_dwflocality_SOURCES = t_dwflocality.cc
t_dwflocality_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_invert3_precwilson_SOURCES = t_invert3_precwilson.cc
t_invert3_precwilson_DEPENDENCIES =  build_chroma_libs rebuild_other_libs

t_invert4_precwilson_SOURCES = t_invert4_precwilson.cc invcg2_timing_hacks_2.h invcg2_timing_hacks_2.cc invcg2_prec_wilson.h invcg2_prec_wilson.cc
t_invert4_precwilson_DEPENDENCIES =  build_chroma_libs rebuild_other_libs

t_read_eigen_SOURCES = t_read_eigen.cc
t_read_eigen_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_bicgstab_SOURCES = t_bicgstab.cc
t_bicgstab_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_invborici_SOURCES = t_invborici.cc
t_invborici_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_hamsys_SOURCES = t_hamsys.cc
t_hamsys_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_hamsys_ferm_SOURCES = t_hamsys_ferm.cc
t_hamsys_ferm_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_hmc_pg_SOURCES = t_hmc_pg.cc
t_hmc_pg_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_precact_sse_SOURCES = t_precact_sse.cc
t_precact_sse_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_prec_contfrac_SOURCES = t_prec_contfrac.cc
t_prec_contfrac_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_mres_4d_SOURCES = t_mres_4d.cc
t_mres_4d_DEPENDENCIES = build_chroma_libs rebuild_other_libs

t_prec_nef_SOURCES = t_prec_nef.cc
t_prec_nef_DEPENDENCIES = build_chroma_libs rebuild_other_libs

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_chroma_libs
build_chroma_libs:
	cd $(BUILDDIR)/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd $(BUILDDIR)/other_libs ; $(MAKE)

# xmldiff targets
.PHONY: xcheck
xcheck:
	perl ../../scripts/run_xmldiff.pl
