cmake_minimum_required(VERSION 3.17)
project(Chroma VERSION 3.44.1 LANGUAGES CXX C)
include(FindOpenMP)

option(Chroma_ENABLE_CPP_WILSON_DSLASH "Use CPP Dslash for efficiecy" OFF)
option(Chroma_ENABLE_SSE2 "Use SSE2 Optimizations" OFF)
option(Chroma_ENABLE_SSE3 "Use SSE3 Optimizations" OFF)
option(Chroma_ENABLE_SANITIZERS "Enable memory and address sanitizers" OFF)
#option(Chroma_ENABLE_CG_SOLVER_RESTART "Enable A single Restart in the CG Solvers" OFF)
#option(Chroma_ENABLE_MDWF "Use MDWF Solver From MIT" OFF)
#koption(Chroma_ENABLE_LLVM_WILSON_DSLASH "Use LLVM Based Wilson Dslash (QDP-JIT builds only)" OFF)
#option(Chroma_ENABLE_STATIC_PACKED_GAUGE "Use static space for packed gauges" OFF)
#option(Chroma_ENABLE_SSED_CLOVER "Use SSE Double Prec for Clover Op" OFF)
option(Chroma_ENABLE_JIT_CLOVER "Use JIT-ed Clover Operator (QDP-JIT builds only" OFF)
#option(Chroma_ENABLE_QUDA_DEVICE_INTERFACE "Use Device Interfae with QUDA (QDP-JIT Builds only)" OFF)
#option(Chroma_ENABLE_JIT_CONTRACTION_KERNELS "Enable JIT-ed Contraction kernels (QDP-JIT Builds only)" OFF)
#option(Chroma_ENABLE_OPT_EIGCG "Enable EigCG Solver" OFF)
#option(Chroma_ENABLE_OPT_CFZ_LINOP "Enable Optimized Continued Fraction operator" OFF)
#option(Chroma_ENABLE_GENERIC_BICGSTAB_KERNELS "Enable generic-C BicGStab kernels" OFF)
#option(Chroma_ENABLE_SSE_BICGSTAB_KERNELS "Enable SSE optimized BiCGStab kernels" OFF)
option(Chroma_ENABLE_GTEST "Enable Google Test based tests" OFF)
option(Chroma_ENABLE_QUDA  "Enable the QUDA Library for GPUs" OFF)
#option(Chroma_ENABLE_QUDA_DEBUG "Enable Extra debugging for QUDA" OFF)
#option(Chroma_ENABLE_QUDA_MG_DUMP "Enable to dump on failed QUDA solves" OFF)
option(Chroma_ENABLE_QPHIX  "Enable QPhiX Solvers" OFF)
option(Chroma_ENABLE_MGPROTO "Enable MG-Proto" OFF)
option(Chroma_ENABLE_LAPACK "Enable Lapack Binding in QDPLapack" OFF)
option(Chroma_ENABLE_OPENMP "Enable OpenMP" OFF)
option(Chroma_ENABLE_SUPERBBLAS  "Enable tasks using Superbblas" OFF)
option(Chroma_ENABLE_PRIMME  "Build with PRIMME" OFF)
set(BLAS_SUFFIX_LINK_FLAGS "" CACHE STRING "Extra linking flags placed after BLAS")
separate_arguments(Chroma_extra_link_flags UNIX_COMMAND ${BLAS_SUFFIX_LINK_FLAGS})
message(STATUS "Chroma extra links: ${Chroma_extra_link_flags}")

set(Chroma_QPHIX_SOLVER_ARCH "avx" CACHE STRING "Architecture of QPhiX Architecture")
set_property(CACHE Chroma_QPHIX_SOLVER_ARCH PROPERTY
	STRINGS sse avx avx2 avx512 mic
)

set(Chroma_QPHIX_SOALEN "" CACHE STRING "SOALen for QPhiX")
set(Chroma_QPHIX_INNER_SOALEN "" CACHE STRING "SOALEN for QPhiX Inner Solver")

set(Chroma_QPHIX_INNER_PRECISION "float" CACHE STRING "Float Type for Inner QPhiX Solever")
set_property(CACHE Chroma_QPHIX_INNER_PRECISION PROPERTY
	STRINGS half float double
)

option(Chroma_QPHIX_COMPRESS12 "Set to ON to enable 12-compression" OFF)

option(Chroma_QPHIX_DSLASH "Set to ON to enable QPhiX Dslash. QDP-JIT builds with matchin SOALENs only" OFF)
option(Chroma_USE_MGPROTO "Use the MGProto Multigrid solver" OFF)

set(SUPERBBLAS_DIR "" CACHE PATH "path to superbblas installation")
set(PRIMME_DIR "" CACHE PATH "path to PRIMME installation")

##################
# End of Options
##################

##################
# Configuration
##################

if( Chroma_ENABLE_JIT_CLOVER )
 set(BUILD_JIT_CLOVER_TERM 1)
endif()


# SSE3 automatically sets SSE2
if( Chroma_ENABLE_SSE3 )
  set(Chroma_ENABLE_SSE2)
endif()

##########################
# Check if we can find GMP
##########################
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

find_package(GMP)
if( GMP_FOUND ) 
   set(BUILD_GMP_REMEZ 1) # For the header file
endif()
# Create an imported target

##################
# LibGMP and Remez
##################

###################
# OpenMP
###################
if(Chroma_ENABLE_OPENMP)
  if( NOT OpenMP_FOUND )
  	find_package(OpenMP REQUIRED)
  endif()
  
  if( NOT Threads_FOUND )
  	find_package(Threads REQUIRED)
  endif()
endif()

###################
# 1st: QDPXX
###################
find_package(QDPXX REQUIRED)

if( Chroma_ENABLE_SANITIZERS )
       include(cmake/CheckSanitizeOpts.cmake)
       check_sanitizer_options( "${Chroma_ENABLE_SANITIZERS}" Chroma_SANITIZER_OPTS )
       message(STATUS "Chroma: Setting Sanitizer options: ${Chroma_SANITIZER_OPTS}")
       set(CPPWilsonDslash_ENABLE_SANITIZERS ${Chroma_ENABLE_SANITIZERS} )
       set(QDPLapack_ENABLE_SANITIZERS ${Chroma_ENABLE_SANITIZERS} )
endif()

##########################
# 2nd: Deal with QPhiX
##########################
if( Chroma_ENABLE_QPHIX )
#  
#  Deal with QPhix
#  find_package(QPhiX)
endif()


##########################
# 3rd: Deal with MGProto
#########################
if( Chroma_ENABLE_MGPROTO)
#
# Deal with MG Proto
# find_package(MGPROTO)
#
endif()

###########################
# 4th: Deal with QUDA
###########################
if( Chroma_ENABLE_QUDA )
  message(STATUS "Enabling QUDA")
  set(BUILD_QUDA 1)
  find_package(QUDA REQUIRED)
  if( QUDA_TARGET_CUDA )
    set(CHROMA_TARGET_CUDA 1)
  elseif( QUDA_TARGET_HIP )
    set(CHROMA_TARGET_HIP 1)
  endif()
  if (QDP_IS_QDPJIT)
     message(STATUS "QDP is QDPJIT so enabling QDPJIT Interfaces")
     set(BUILD_QUDA_DEVIFACE_CLOVER 1)
     set(BUILD_QUDA_DEVIFACE_GAUGE 1)
     set(BUILD_QUDA_DEVIFACE_SPINOR 1)
  endif()
endif()

###########################
# 5th: Deal with MDWF
###########################
if( Chroma_ENABLE_MDWF )
# 
# Gonna need to write a find module for this
# 
endif()

####################################################
# 4th: Configure and link against CPPWilsonDslash
####################################################
if(Chroma_ENABLE_CPP_WILSON_DSLASH)
  set(BUILD_CPP_WILSON_DSLASH 1)  # For include file
  # Set Cache Options for CPP Wilson Dslash
  # Other options can be set with -DCPPWilsonDslash_<OPTION>=ON etc
  set(CPPWilsonDslash_ENABLE_QDPXX ON)
  add_subdirectory(other_libs/cpp_wilson_dslash)
endif()

#####################################
# 5th: Setup QDPLapack options
####################################
if( Chroma_ENABLE_LAPACK )
  set(BUILD_LAPACK 1)
  set(QDPLapack_BINDING "lapack")
endif()
add_subdirectory(other_libs/qdp-lapack)
 
###########################
# 6th: Deal with superbblas
###########################
if( Chroma_ENABLE_SUPERBBLAS )
  message(STATUS "Enabling superbblas")

  include(FindBLAS)
  find_package(BLAS REQUIRED)
  
  if( QDP_BACKEND_ROCM AND NOT HIP_CPP_CONFIG )
     find_package(HIP REQUIRED)
     execute_process(
                 COMMAND ${HIP_HIPCONFIG_EXECUTABLE} -C
                 OUTPUT_VARIABLE _hip_cpp_config_raw
                 OUTPUT_STRIP_TRAILING_WHITESPACE
                 )
     separate_arguments(_hip_cpp_config UNIX_COMMAND ${_hip_cpp_config_raw})
     message(STATUS "HIP cpp_config ${_hip_cpp_config}")
     set(HIP_CPP_CONFIG ${_hip_cpp_config} CACHE STRING "HIP cpp_config as computed by hipconfig")
     mark_as_advanced(HIP_CPP_CONFIG)
  endif()

  if( QDP_BACKEND_ROCM )
     find_package(hipBLAS REQUIRED)
  endif()
  if( Chroma_ENABLE_CUDA )
     find_package(CUDAToolkit REQUIRED)
  endif()

  set(BUILD_SB 1)

  # Deal with PRIMME
  if( Chroma_ENABLE_PRIMME )
    message(STATUS "Enabling PRIMME")
    set(BUILD_PRIMME 1)
  endif()
endif()


#####################################
# Configure chroma_config_internal.h
#####################################
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lib/chroma_config_internal.h.cmake.in
			   lib/chroma_config_internal.h)
			
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/chroma_config_internal.h DESTINATION include)

# Now we are ready to build in lib
# Heaven help us. 

add_subdirectory(lib)
add_subdirectory(mainprogs/main)
add_subdirectory(mainprogs/tests)

#install the FindGMP module so we can find GMP in depenendices
# if needed
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGMP.cmake
		DESTINATION lib/cmake/Chroma )
		  
#install the export description of the targets
install(EXPORT ChromaTargets
        FILE ChromaTargets.cmake
        NAMESPACE Chroma::
        DESTINATION lib/cmake/Chroma)

#Write a basic Verson file
include(CMakePackageConfigHelpers)
configure_package_config_file(ChromaConfig.cmake.in ChromaConfig.cmake
                              INSTALL_DESTINATION lib/cmake/Chroma)

write_basic_package_version_file(
	ChromaConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ChromaConfigVersion.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/ChromaConfig.cmake
              DESTINATION lib/cmake/Chroma)
