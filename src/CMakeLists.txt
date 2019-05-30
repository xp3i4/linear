# ===========================================================================
# Linear: Mapping reads efficiently
# CMakeLists.txt.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (linear)
message (STATUS "Configuring linear/src")

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package (SeqAn REQUIRED)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
#set (CMAKE_BUILD_TYPE Release)
set (CMAKE_BUILD_TYPE Debug)

add_library(linear_base STATIC base.cpp)
add_library(linear_cords STATIC cords.cpp)
add_library(linear_shape_extend STATIC shape_extend.cpp)
add_library(linear_index_util STATIC index_util.cpp)
add_library(linear_args_parser STATIC args_parser.cpp)
add_library(linear_pmpfinder STATIC pmpfinder.cpp)
add_library(linear_gap STATIC gap.cpp)
add_library(linear_fio STATIC f_io.cpp)
add_library(linear_align_interface STATIC align_interface.cpp)

add_executable (linear mapper.cpp) 
target_link_libraries (linear 
                       linear_align_interface
                       linear_args_parser 
                       linear_pmpfinder 
                       linear_gap 
                       linear_fio 
                       linear_index_util 
                       linear_shape_extend
                       linear_cords 
                       linear_base
                      ${SEQAN_LIBRARIES})

add_executable (test_case test_case.cpp)
target_link_libraries (test_case  
                       linear_shape_extend
                       linear_base 
                       linear_args_parser
                      ${SEQAN_LIBRARIES}
    )
install (TARGETS linear RUNTIME DESTINATION bin)

