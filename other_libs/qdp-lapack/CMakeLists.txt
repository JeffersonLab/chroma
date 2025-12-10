cmake_minimum_required(VERSION 3.17)
PROJECT(QDPLapack VERSION 1.0.0 LANGUAGES CXX C)

if(NOT QDPLapack_BINDING)
	set(QDPLapack_BINDING "dummy")
endif()
	
set(QDPLapack_BINDING ${QDPLapack_BINDING} CACHE STRING "Lapack Binding. One of dummy or lapack (default is dummy")


option(QDPLapack_BLAS_CDOT "Enable Blas CDOT" OFF)
option(QDPLapack_ENABLE_SANITIZERS "Enable address and undefined behaviour sanitizers" OFF)

string(TOLOWER ${QDPLapack_BINDING} qdpLapackBinding)

# Check Binding 
if( qdpLapackBinding STREQUAL "lapack")
  enable_language(Fortran)
  include(FindLAPACK)
  message(STATUS "QDPLapack: Building LAPACK Bindings")
  set(FORTRAN_LAPACK_BINDING 1)  # Header file substitution
  find_package(LAPACK REQUIRED)
else() 
  message(STATUS "QDPLapack: Building DUMMY Bindings")
  set(DUMMY_LAPACK_BINDING 1)    # Header file substitution
endif()

if(QDPLapack_BLAS_CDOT)
  message(STATUS "QDPLapack: Using BLAS CDOT")
  set(USE_BLAS_CDOT 1)   # Header file substitution
endif()

# QDP is required
message(STATUS "QDPLapack: Looking for QDP++") 
find_package(QDPXX REQUIRED)

# If QDPXX is parallel in order to pass link test QMP must have been specified
# So I do not need to check for QMP specifically it should alread be on the inclulde 
# and link lines, at best I need to set the -D use QMP falg in the configured header
if(QDP_ARCH_PARSCALAR) 
  message(STATUS "QDPLapack: QDP++ is Parscalar")
  set(USE_QMP 1)
endif()

# Set up configurtion file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/qdp-lapack_config_internal.h.cmake.in
			   ${CMAKE_CURRENT_BINARY_DIR}/include/qdp-lapack_config_internal.h)
			   

# Install genenerated headers
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/qdp-lapack_config_internal.h DESTINATION include)

# Install headers in include
install(DIRECTORY include DESTINATION . FILES_MATCHING PATTERN "*.h")

# Deal with Sanitizer
if( QDPLapack_ENABLE_SANITIZERS )
  include(cmake/CheckSanitizeOpts.cmake)
  check_sanitizer_options( "${QDPLapack_ENABLE_SANITIZERS}" QDPLapack_SANITIZER_OPTS )
  message(STATUS "QDPLapack: Setting Sanitizer options: ${QDPLapack_SANITIZER_OPTS}")
endif()

add_subdirectory(lib)

 # Make the config file
include(CMakePackageConfigHelpers)
configure_package_config_file(Config.cmake.in QDPLapackConfig.cmake
                              INSTALL_DESTINATION lib/cmake/QDPLapack)

# Make the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    QDPLapackConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/QDPLapackConfigVersion.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/QDPLapackConfig.cmake
              DESTINATION lib/cmake/QDPLapack)
		  
  
