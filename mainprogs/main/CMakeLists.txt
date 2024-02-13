cmake_minimum_required(VERSION 3.17)
PROJECT(ChromaMain VERSION 3.44.1 LANGUAGES CXX C)

get_directory_property(SubBuild PARENT_DIRECTORY)

if(NOT SubBuild)
        find_package(Chroma REQUIRED) 
endif()

set( ChromaMain_EXECUTABLES 
	chroma
	purgaug
	cfgtransf
	const_hmc
	hmc
	spectrum_s
    collect_propcomp 
    qpropgfix qproptrev 
    qpropqio  
    qproptransf 
    wallformfac 
)

# These didn't compile
#
#  multi_propagator_comp
#  collect_multi_propcomp
#  qpropadd
#  propagator_comp 
#
foreach( EXEC ${ChromaMain_EXECUTABLES})
  add_executable( ${EXEC} ${EXEC}.cc )

  if (QDP_IS_QDPJIT)
    set_target_properties( ${EXEC} PROPERTIES CXX_STANDARD 20 )

    if( QDP_BACKEND_L0 )
      find_package(PkgConfig REQUIRED)
      pkg_search_module(PKG_LLVMSPIRV REQUIRED LLVMSPIRVLib)
      target_link_libraries(${EXEC} PUBLIC ${PKG_LLVMSPIRV_LIBRARIES})
      target_link_options(${EXEC} PUBLIC ${PKG_LLVMSPIRV_LDFLAGS})
    endif()

  endif()
  
  target_link_libraries( ${EXEC} PUBLIC Chroma::chromalib )
endforeach()

INSTALL( TARGETS hmc const_hmc chroma purgaug cfgtransf 
		 EXPORT ChromaTargets
		 ARCHIVE DESTINATION lib
 		 LIBRARY  DESTINATION lib
  		 RUNTIME DESTINATION bin
  		 INCLUDES DESTINATION include
)
