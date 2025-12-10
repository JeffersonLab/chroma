add_library(qdp_lapack)

if( qdpLapackBinding STREQUAL "dummy")
  target_sources(qdp_lapack PRIVATE
  	lapack_stub.cc
  	IncrEigpcg_stub.c 
  	IncrEigbicg_stub.c
  )
else()
  target_sources(qdp_lapack PRIVATE
  	fortran_lapack.cc 
  	Complex.c
  	IncrEigpcg.c  
  	numerical.c  
  	restart_X.c 
  	eigpcg.c   
  	ortho.c 
  	wtime.c 
  	cdotcsub.f  
  	bicgstab.c 
  	biortho.c 
  	conjugate.c 
  	eigbicg.c 
  	G_eval.c 
  	IncrEigbicg.c 
  	sort.c 
  	zdotcsub.f	
  )
  target_link_libraries(qdp_lapack PUBLIC LAPACK::LAPACK)
endif()


if (QDP_IS_QDPJIT)
  set_target_properties( qdp_lapack PROPERTIES CXX_STANDARD 20 )
endif()


target_link_libraries(qdp_lapack PUBLIC QDPXX::qdp)
if( QDPLapack_SANITIZER_OPTS )  
	target_compile_options(qdp_lapack PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${QDPLapack_SANITIZER_OPTS}>)
  target_compile_options(qdp_lapack PRIVATE $<$<COMPILE_LANGUAGE:C>:${QDPLapack_SANITIZER_OPTS}>)
	target_link_options(qdp_lapack PUBLIC ${QDPLapack_SANITIZER_OPTS})
endif()

target_include_directories(qdp_lapack PUBLIC
							$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
							$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include>
							$<INSTALL_INTERFACE:include>
						   )
						   
install(TARGETS qdp_lapack EXPORT QDPLapackTargets
	    ARCHIVE DESTINATION lib
  		LIBRARY  DESTINATION lib
  		RUNTIME DESTINATION bin
	    INCLUDES DESTINATION include
)

#install the export description of the targets
install(EXPORT QDPLapackTargets
        FILE QDPLapackTargets.cmake
        NAMESPACE QDPLapack::
        DESTINATION lib/cmake/QDPLapack)

add_library(QDPLapack::qdp_lapack ALIAS qdp_lapack)
