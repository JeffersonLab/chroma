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
    set_target_properties( ${EXEC} PROPERTIES CXX_STANDARD 14 )
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
