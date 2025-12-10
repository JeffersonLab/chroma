cmake_minimum_required(VERSION 3.17)
PROJECT(CPPWilsonDslashTests VERSION 1.0.0 LANGUAGES CXX C)

get_directory_property(SubBuild PARENT_DIRECTORY)

if(NOT SubBuild)
        find_package(CPPWilsonDslash REQUIRED)
endif()

add_executable(t_dslash)
add_executable(t_dslash_3d)
add_executable(time_dslash)
add_executable(time_dslash_3d)

add_library(test_utils 
	unittest.h
	dslashm_w.h
	dslashm_w.cc
	reunit.h
	reunit.cc
)
target_link_libraries(test_utils PUBLIC CPPWilsonDslash::dslash)

add_library(test_utils_3d
	unittest.h
	dslashm_w.h
	dslashm_w.cc
	dslash_3d_w.h
	dslash_3d_w.cc
	reunit.h
	reunit.cc
)
target_link_libraries(test_utils_3d PUBLIC CPPWilsonDslash::dslash)

target_sources(t_dslash PUBLIC
	testDslashFull.h 
	testDslashFull.cc
	t_dslash.cc
)
target_link_libraries(t_dslash PUBLIC test_utils)

target_sources(t_dslash_3d PUBLIC 
	testDslash3D.h
	testDslash3D.cc
	t_dslash_3d.cc
)
target_link_libraries(t_dslash_3d PUBLIC test_utils_3d)

target_sources(time_dslash PUBLIC 
	timeDslash.h
	timeDslash.cc
	time_dslash.cc
)
target_link_libraries(time_dslash PUBLIC test_utils)


target_sources(time_dslash_3d PUBLIC
	timeDslash3D.h
	timeDslash3D.cc
	time_dslash_3d.cc
)
target_link_libraries(time_dslash_3d PUBLIC test_utils_3d)
