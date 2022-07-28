PROJECT(ChromaTests VERSION 3.44.1 LANGUAGES CXX)

get_directory_property(SubBuild PARENT_DIRECTORY)

if(NOT SubBuild)
  find_package(Chroma REQUIRED) 
endif()

set(ChromaSingleExecs
	t_mesplq
	t_ritz_KS 
	t_lwldslash_array 
	t_leapfrog 
	t_lwldslash_new 
	t_minvert 
	t_meas_wilson_flow
	t_io t_mesons_w
	t_conslinop
	t_hypsmear
    t_ape_smear
    t_dwf4d 
    t_propagator_s 
    t_disc_loop_s
    t_ritz
    t_dwflocality
    t_precact_4d
    t_precact_5d
    t_gauge_force
    t_stout_state
    t_aniso_gaugeact
    t_temp_prec
    t_multishift_quda
    t_invmdagm_quda
)

if( GMP_FOUND )
  list(APPEND ChromaSingleExecs t_remez )
endif()

############################
# IS THIS STILL RELEVANT
############################
#EXTRA_PROGRAMS  = t_dslashm t_lwldslash t_spprod \
#    t_follana_io_s t_formfac  t_follana_pion_s t_dwflinop t_precdwf \
#    t_ovlap_bj t_propagator_s t_propagator_w t_gfix t_ovlap5d_bj \
#    t_ritz5d_KS  t_fermion_loop_w t_propagator_fuzz_s t_propagator_fuzz_baryon_s\
#    t_disc_loop_s t_ov_pbp t_propagator_nrqcd \
#    t_invert3_precwilson t_invert4_precwilson t_read_eigen \
#    t_fuzwilp t_ovlap_double_pass t_g5eps_bj t_sumr t_msumr t_invrelcg t_bicgstab \
#    t_wilslp t_rel_gmresr t_invborici t_hamsys t_hamsys_ferm \
#    t_precact_sse t_preccfz_opt \
#    t_hmc_pg t_neflinop t_prec_contfrac t_mres_4d t_prec_nef t_su3 \
#    t_overbu t_lower_tests t_seqsource \
#    t_unprec_wilson_force t_prec_wilson_force t_unprec_twoflav_wilson_monomial \
#    t_prec_twoflav_wilson_monomial t_hamiltonian t_hmc \
#    t_stagg_baryon \
#    t_propagator_twisted \
#    t_circular_buffer \
#    t_precnef \
#    t_preccfz \
#    t_aniso_sym_force \
#    t_monomial_force \
#    t_clover \
#    t_db \
#    t_solver_accum \
#    t_eigcginv

#### Conditionals 
#### Add Later
####	
#if BUILD_QUDA
#check_PROGRAMS += t_quda_tprec t_minvert_quda
#endif
#if BUILD_QPHIX
#check_PROGRAMS += t_minvert_qphix
#endif
#if BUILD_MGPROTO
#check_PROGRAMS += t_mgproto
#endif
#if BUILD_GTEST
#check_PROGRAMS += t_inv_fgmres_dr  t_symm_prec
	
	
#t_inv_fgmres_dr_SOURCES = t_inv_fgmres_dr.cc chroma_gtest_env.h \
#	fgmres_dr_tests.cc
	
#t_symm_prec_SOURCES = t_symm_prec.cc chroma_gtest_env.h \
#	symm_prec_xml.h symm_prec_tests.cc
#endif

# Single source guys
foreach( EXEC ${ChromaSingleExecs})
  add_executable( ${EXEC} ${EXEC}.cc )

  if (QDP_IS_QDPJIT)
    if(QDP_PROP_OPT)
      set_target_properties( ${EXEC} PROPERTIES CXX_STANDARD 20 )
    endif()
    
    if( QDP_BACKEND_L0 )
      find_package(PkgConfig REQUIRED)
      pkg_search_module(PKG_LLVMSPIRV REQUIRED LLVMSPIRVLib)
      target_link_libraries(${EXEC} PUBLIC ${PKG_LLVMSPIRV_LIBRARIES})
      target_link_options(${EXEC} PUBLIC ${PKG_LLVMSPIRV_LDFLAGS})
    endif()

  endif()
  
  target_link_libraries( ${EXEC} PUBLIC Chroma::chromalib )
  target_include_directories( ${EXEC} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()

# The GTest ones
add_library(GTest ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/gtest/gtest-all.cc)
target_include_directories(GTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib)
set_target_properties(GTest PROPERTIES CXX_STANDARD 11 CXX_EXTENSIONS OFF)

add_executable(t_symm_prec 
	t_symm_prec.cc 
	chroma_gtest_env.h 
	symm_prec_xml.h 
	symm_prec_tests.cc
	)

if (QDP_IS_QDPJIT)
  if(QDP_PROP_OPT)
    set_target_properties( t_symm_prec PROPERTIES CXX_STANDARD 20 )
  endif()
    
  if( QDP_BACKEND_L0 )
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PKG_LLVMSPIRV REQUIRED LLVMSPIRVLib)
    target_link_libraries(t_symm_prec PUBLIC ${PKG_LLVMSPIRV_LIBRARIES})
    target_link_options(t_symm_prec PUBLIC ${PKG_LLVMSPIRV_LDFLAGS})
  endif()
endif()

target_link_libraries( t_symm_prec PUBLIC Chroma::chromalib GTest)
target_include_directories(t_symm_prec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


if ( Chroma_ENABLE_LAPACK )
add_executable(t_inv_fgmres_dr 
	t_inv_fgmres_dr.cc 
	chroma_gtest_env.h
	fgmres_dr_tests.cc
	)

if (QDP_IS_QDPJIT)
  if(QDP_PROP_OPT)
    set_target_properties( t_inv_fgmres_dr PROPERTIES CXX_STANDARD 20 )
  endif()
    
  if( QDP_BACKEND_L0 )
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PKG_LLVMSPIRV REQUIRED LLVMSPIRVLib)
    target_link_libraries(t_inv_fgmres_dr PUBLIC ${PKG_LLVMSPIRV_LIBRARIES})
    target_link_options(t_inv_fgmres_dr PUBLIC ${PKG_LLVMSPIRV_LDFLAGS})
  endif()
endif()

target_link_libraries( t_inv_fgmres_dr PUBLIC Chroma::chromalib GTest)
target_include_directories(t_inv_fgmres_dr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()
