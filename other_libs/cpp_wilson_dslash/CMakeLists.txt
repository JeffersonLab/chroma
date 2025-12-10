# Boilerplate 
cmake_minimum_required(VERSION 3.17)
project(CPPWilsonDslash VERSION 1.0.0 LANGUAGES CXX C)


if( CPPWilsonDslash_ENABLE_QDPXX )
  message(STATUS "CPPWilsonDslash: Using QDP++")
  find_package(QDPXX REQUIRED)
endif()

set(CPPWilsonDslash_PARALLEL_ARCH ${QDP_PARALLEL_ARCH} CACHE STRING 
    "The parallel arch: values are scalar and parscalar. Default is from QDPXX if QDPXX is enabled")
 
set_property(CACHE CPPWilsonDslash_PARALLEL_ARCH PROPERTY
		STRINGS scalar parscalar)
# --------------------------------------------------
# Cache file options
#     Things one can set in CCMake or CMakeGUI 
#     or on the command line
# --------------------------------------------------
option(CPPWilsonDslash_ENABLE_OPENMP "Use OpenMP Threads" OFF)
option(CPPWilsonDslash_ENABLE_SSE2 "USE SSE2" OFF )

option(CPPWilsonDslash_ENABLE_NOCOMMS "Disable Communications" OFF)
option(CPPWilsonDslash_ENABLE_NOCOMPUTE "Disable Computation"  OFF)

option(CPPWilsonDslash_ENABLE_SANITIZERS "Enable memory and undefined behaviour Sanitizers" OFF)

# Configuration based on the options: Parallel arch
# Convert user entry to lower so we only need to match
string(TOLOWER ${CPPWilsonDslash_PARALLEL_ARCH} dslashParallelArch)
if( dslashParallelArch STREQUAL "scalar")
  message(STATUS "CPPWilsonDslash: using scalar arch")
  set(CPP_DSLASH_SCALAR 1)      # Macro in Config internal.h
elseif( dslashParallelArch STREQUAL "parscalar")
  message(STATUS "CPPWilsonDslash: using parscalar arch")
  set(CPP_DSLASH_PARSCALAR 1)   # Macro in  Config internal.h
else() 
  message(FATAL_ERROR "CPPWilsonDslash: Unknown value ${CPP_DSLASH_PARALLEL_ARCH}. Specify scalar or parscalar or use QDPXX")
endif()


# Whether to use SSE2
if( CPPWilsonDslash_ENABLE_SSE2 ) 
  message(STATUS "Setting Dslash SSE2")
  set(DSLASH_USE_SSE2 1)    # Convert to value used by header
endif()

# Whether to use OpenMP
if( CPPWilsonDslash_ENABLE_OPENMP ) 
  set(DSLASH_USE_OMP_THREADS 1)
  message(STATUS "Enabling OpenMP Threads in CPPWilsonDslash")
   
  # OpenMP may already have been found, e.g. as part of QDPXX
  if( NOT OpenMP_CXX_FOUND )
	find_package(OpenMP REQUIRED)
  	find_package(Threads REQUIRED)
  endif()
endif()

# Whether to enable nocomute and nocomms 
# options. Do these still do anything useful
if( CPPWilsonDslash_NOCOMPUTE )
  message(STATUS "Disabling Compute")
  set(SSEDSLASH_4D_NOCOMPUTE 1)
endif()

if( CPPWilsonDslash_NOCOMMS )
  message(STATUS "Disabling Comms")
  set(SSEDSLASH_4D_NOCOMMS 1)
endif()

# Write the config_internal.h file
configure_file(include/dslash_config_internal.h.cmake.in 
			   include/dslash_config_internal.h)

# Install it? (Or make it a public header later?
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/dslash_config_internal.h DESTINATION include)



# Deal with Sanitizer
if( CPPWilsonDslash_ENABLE_SANITIZERS )
  include(cmake/CheckSanitizeOpts.cmake)
  check_sanitizer_options( "${CPPWilsonDslash_ENABLE_SANITIZERS}" CPPWilsonDslash_SANITIZER_OPTS )
  message(STATUS "CPPWilsonDslash: Setting Sanitizer options: ${CPPWilsonDslash_SANITIZER_OPTS}")
endif()
	
add_subdirectory(lib)
add_subdirectory(tests)

#install the headers 
install(DIRECTORY include DESTINATION . 
	    FILES_MATCHING PATTERN "*.h")
	    
	   # Make the config file
include(CMakePackageConfigHelpers)
configure_package_config_file(Config.cmake.in CPPWilsonDslashConfig.cmake 
							  INSTALL_DESTINATION lib/cmake/CPPWilsonDslash)

# Make the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    CPPWilsonDslashConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CPPWilsonDslashConfigVersion.cmake 
			  ${CMAKE_CURRENT_BINARY_DIR}/CPPWilsonDslashConfig.cmake 
			  DESTINATION lib/cmake/CPPWilsonDslash)
	    
