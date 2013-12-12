#
TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

#
# Include Path
#
INCFLAGS=-I$(TOPDIR)/lib -I$(BUILDDIR)/lib

#
# Any other compiler flags
#
AM_CXXFLAGS = $(INCFLAGS) @QDPXX_CXXFLAGS@ @BAGEL_WILSON_DSLASH_CXXFLAGS@ @GMP_CXXFLAGS@ @QOP_MG_CXXFLAGS@ @CXXFLAGS@
AM_LDFLAGS = -L$(BUILDDIR)/lib @QDPXX_LDFLAGS@ @BAGEL_WILSON_DSLASH_LDFLAGS@ @GMP_LDFLAGS@ @LDFLAGS@
LDADD = -lchroma @QDPXX_LIBS@ @BAGEL_WILSON_DSLASH_LIBS@ @GMP_LIBS@  @LIBS@

AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@QDP_LAPACK_DIR@/lib
LDADD += -lqdp-lapack

if BUILD_SSE_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@SSE_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@SSE_DSLASH_DIR@/lib
LDADD += -llevel3 @QDPXX_LIBS@
endif

if BUILD_CPP_WILSON_DSLASH
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@CPP_DSLASH_DIR@/include -I$(BUILDDIR)/other_libs/@CPP_DSLASH_DIR@/include
AM_LDFLAGS += -L$(BUILDDIR)/other_libs/@CPP_DSLASH_DIR@/lib
LDADD += -ldslash @QDPXX_LIBS@
endif

if BUILD_CG_DWF
AM_CXXFLAGS += -I$(TOPDIR)/other_libs/@CG_DWF_DIR@
AM_CXXFLAGS += -I$(BUILDDIR)/other_libs/@CG_DWF_DIR@
AM_LDFLAGS  += -L$(BUILDDIR)/other_libs/@CG_DWF_DIR@
LDADD += -lcg-dwf
endif

if BUILD_MDWF
AM_CXXFLAGS += @MDWF_CXXFLAGS@
AM_LDFLAGS  += @MDWF_LDFLAGS@
LDADD     += @MDWF_LIBS@
endif

if BUILD_QUDA
AM_CXXFLAGS += @QUDA_CXXFLAGS@
AM_LDFLAGS += @QUDA_LDFLAGS@
LDADD	+= @QUDA_LIBS@
endif

if BUILD_QOP_MG
AM_LDFLAGS += -L@top_builddir@/other_libs/wilsonmg/lib @QOP_MG_LDFLAGS@
LDADD += @QOP_MG_LIBS@
endif

if BUILD_BAGEL_CLOVER_APPLY
AM_CXXFLAGS += @BAGEL_CLOVER_CXXFLAGS@
AM_LDFLAGS += @BAGEL_CLOVER_LDFLAGS@
LDADD += @BAGEL_CLOVER_LIBS@
endif

if BUILD_QMT
AM_CXXFLAGS += @QMT_CXXFLAGS@
AM_LDFLAGS  += @QMT_LDFLAGS@
LDADD += @QMT_LIBS@
endif

#
# Local Headers
#
HDRS =

## Production tests:
bin_PROGRAMS = t_mesplq t_lwldslash_sse t_lwldslash_pab t_ritz_KS t_lwldslash_array t_leapfrog t_lwldslash_new t_minvert t_meas_wilson_flow 

#
# add the programs to build in here
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.
check_PROGRAMS  = t_io t_mesons_w  t_conslinop t_hypsmear \
    t_ape_smear t_dwf4d t_propagator_s t_disc_loop_s \
    t_remez t_ritz t_dwflocality t_precact_4d t_precact_5d \
    t_gauge_force t_stout_state t_aniso_gaugeact t_temp_prec t_meas_wilson_flow_loop

if BUILD_QUDA
check_PROGRAMS += t_quda_tprec t_minvert_quda
endif


#
## Move programs under development to EXTRA_PROGRAMS if they currently
## do not compile.  This way, they don't break 'make check' but the targets
## are still there for further development
#
EXTRA_PROGRAMS  = t_dslashm t_lwldslash t_spprod \
    t_follana_io_s t_formfac  t_follana_pion_s t_dwflinop t_precdwf \
    t_ovlap_bj t_propagator_s t_propagator_w t_gfix t_ovlap5d_bj \
    t_ritz5d_KS  t_fermion_loop_w t_propagator_fuzz_s t_propagator_fuzz_baryon_s\
    t_disc_loop_s t_ov_pbp t_propagator_nrqcd \
    t_invert3_precwilson t_invert4_precwilson t_read_eigen \
    t_fuzwilp t_ovlap_double_pass t_g5eps_bj t_sumr t_msumr t_invrelcg t_bicgstab \
    t_wilslp t_rel_gmresr t_invborici t_hamsys t_hamsys_ferm \
    t_precact_sse t_preccfz_opt \
    t_hmc_pg t_neflinop t_prec_contfrac t_mres_4d t_prec_nef t_su3 \
    t_overbu t_lower_tests t_seqsource \
    t_unprec_wilson_force t_prec_wilson_force t_unprec_twoflav_wilson_monomial \
    t_prec_twoflav_wilson_monomial t_hamiltonian t_hmc \
    t_stagg_baryon \
    t_propagator_twisted \
    t_circular_buffer \
    t_precnef \
    t_preccfz \
    t_aniso_sym_force \
    t_monomial_force \
    t_clover \
    t_db \
    t_solver_accum \
    t_eigcginv

#
# The programs and their dependencies
#
t_su3_SOURCES = t_su3.cc
t_mesplq_SOURCES = t_mesplq.cc
t_seqsource_SOURCES = t_seqsource.cc
t_remez_SOURCES = t_remez.cc
t_ape_smear_SOURCES = t_ape_smear.cc
t_lower_tests_SOURCES = t_lower_tests.cc
t_fuzwilp_SOURCES = t_fuzwilp.cc
t_wilslp_SOURCES = t_wilslp.cc
t_hypsmear_SOURCES = t_hypsmear.cc
t_dslashm_SOURCES = t_dslashm.cc
t_io_SOURCES = t_io.cc
t_lwldslash_SOURCES = t_lwldslash.cc
t_dwflinop_SOURCES = t_dwflinop.cc
t_dwf4d_SOURCES = t_dwf4d.cc
t_precact_4d_SOURCES = t_precact_4d.cc
t_precact_5d_SOURCES = t_precact_5d.cc
t_precdwf_SOURCES = t_precdwf.cc
t_formfac_SOURCES = t_formfac.cc
t_mesons_w_SOURCES = t_mesons_w.cc
t_conslinop_SOURCES = t_conslinop.cc
t_overbu_SOURCES = t_overbu.cc
t_spprod_SOURCES = t_spprod.cc
t_follana_io_s_SOURCES = t_follana_io_s.cc
t_follana_pion_s_SOURCES = t_follana_pion_s.cc
t_lwldslash_sse_SOURCES = t_lwldslash_sse.cc
t_lwldslash_pab_SOURCES = t_lwldslash_pab.cc
t_lwldslash_new_SOURCES = t_lwldslash_new.cc
t_ovlap_bj_SOURCES = t_ovlap_bj.cc
t_ovlap_double_pass_SOURCES = t_ovlap_double_pass.cc
t_g5eps_bj_SOURCES = t_g5eps_bj.cc
t_sumr_SOURCES = t_sumr.cc
t_rel_gmresr_SOURCES = t_rel_gmresr.cc
t_msumr_SOURCES = t_msumr.cc
t_invrelcg_SOURCES = t_invrelcg.cc
t_neflinop_SOURCES = t_neflinop.cc
t_ovlap5d_bj_SOURCES = t_ovlap5d_bj.cc

t_propagator_s_SOURCES = t_propagator_s.cc
t_disc_loop_s_SOURCES = t_disc_loop_s.cc
t_propagator_fuzz_s_SOURCES = t_propagator_fuzz_s.cc
t_stagg_baryon_SOURCES = t_stagg_baryon.cc

t_propagator_fuzz_baryon_s_SOURCES = t_propagator_fuzz_baryon_s.cc
t_propagator_w_SOURCES = t_propagator_w.cc

t_propagator_twisted_SOURCES = t_propagator_twisted.cc
t_propagator_nrqcd_SOURCES = t_propagator_nrqcd.cc
t_fermion_loop_w_SOURCES = t_fermion_loop_w.cc

t_gfix_SOURCES = t_gfix.cc
t_ritz_SOURCES = t_ritz.cc
t_ritz_KS_SOURCES = t_ritz_KS.cc
t_ritz5d_KS_SOURCES = t_ritz5d_KS.cc
t_ov_pbp_SOURCES = t_ov_pbp.cc
t_dwflocality_SOURCES = t_dwflocality.cc
t_invert3_precwilson_SOURCES = t_invert3_precwilson.cc
t_invert4_precwilson_SOURCES = t_invert4_precwilson.cc
# invcg2_timing_hacks_2.h invcg2_timing_hacks_2.cc invcg2_prec_wilson.h invcg2_prec_wilson.cc

t_read_eigen_SOURCES = t_read_eigen.cc
t_bicgstab_SOURCES = t_bicgstab.cc
t_invborici_SOURCES = t_invborici.cc
t_hamsys_SOURCES = t_hamsys.cc
t_hamsys_ferm_SOURCES = t_hamsys_ferm.cc
t_hmc_pg_SOURCES = t_hmc_pg.cc
t_precact_sse_SOURCES = t_precact_sse.cc
t_preccfz_opt_SOURCES = t_preccfz_opt.cc
t_prec_contfrac_SOURCES = t_prec_contfrac.cc
t_mres_4d_SOURCES = t_mres_4d.cc
t_prec_nef_SOURCES = t_prec_nef.cc
t_unprec_wilson_force_SOURCES = t_unprec_wilson_force.cc
t_prec_wilson_force_SOURCES = t_prec_wilson_force.cc
t_unprec_twoflav_wilson_monomial_SOURCES = t_unprec_twoflav_wilson_monomial.cc
t_prec_twoflav_wilson_monomial_SOURCES = t_prec_twoflav_wilson_monomial.cc
t_hamiltonian_SOURCES = t_hamiltonian.cc
t_leapfrog_SOURCES = t_leapfrog.cc
t_hmc_SOURCES = t_hmc.cc
t_gauge_force_SOURCES = t_gauge_force.cc
t_circular_buffer_SOURCES = t_circular_buffer.cc
t_precnef_SOURCES = t_precnef.cc
t_preccfz_SOURCES = t_preccfz.cc
t_lwldslash_array_SOURCES = t_lwldslash_array.cc
t_stout_state_SOURCES = t_stout_state.cc
t_aniso_gaugeact_SOURCES = t_aniso_gaugeact.cc
t_temp_prec_SOURCES = t_temp_prec.cc
t_aniso_sym_force_SOURCES = t_aniso_sym_force.cc
t_monomial_force_SOURCES = t_monomial_force.cc
t_clover_SOURCES = t_clover.cc
t_db_SOURCES = t_db.cc
t_solver_accum_SOURCES = t_solver_accum.cc
t_eigcginv_SOURCES = t_eigcginv.cc

t_meas_wilson_flow_SOURCES  = t_meas_wilson_flow.cc
t_meas_wilson_flow_loop_SOURCES = t_meas_wilson_flow_loop.cc

t_minvert_SOURCES = t_minvert.cc
if BUILD_QUDA
t_quda_tprec_SOURCES = t_quda_tprec.cc
t_minvert_quda_SOURCES = t_minvert_quda.cc
endif

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target

DEPENDENCIES = build_chroma_libs rebuild_other_libs
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

.PHONY: build_chroma_libs
build_chroma_libs:
	cd $(BUILDDIR)/lib ; $(MAKE)

.PHONY: rebuild_other_libs
rebuild_other_libs:
	cd $(BUILDDIR)/other_libs ; $(MAKE)

