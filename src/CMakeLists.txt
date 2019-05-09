# ===========================================================================
# Linear: Mapping reads efficiently
# CMakeLists.txt.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (linear)
message (STATUS "Configuring linear/src")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
#find_package(CXX11 REQUIRED)
#if (NOT CXX11_FOUND)
#    message (STATUS " C++11 compiler needed to build linear")
#    return ()
#endif (NOT CXX11_FOUND)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message (STATUS "OpenMP needed to build linear")
endif()

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")


add_library(linear_base STATIC base.cpp)
add_library(linear_index_util STATIC index_util.cpp)
add_library(linear_args_parser STATIC args_parser.cpp)
add_library(linear_pmpfinder STATIC pmpfinder.cpp)
add_library(linear_index_util2 STATIC index_util2.cpp)
add_library(linear_gap STATIC gap.cpp)
add_library(linear_fio STATIC f_io.cpp)
add_library(linear_align STATIC align_interface.cpp)

add_executable (linear mapper.cpp) 
target_link_libraries (linear ${SEQAN_LIBRARIES} linear_base linear_index_util linear_index_util2 linear_args_parser linear_pmpfinder linear_gap linear_fio linear_align)
add_executable (test_case test_case.cpp)
target_link_libraries (test_case ${SEQAN_LIBRARIES} linear_base linear_index_util linear_args_parser linear_pmpfinder linear_index_util2 linear_gap linear_fio linear_align)

# Add CXX flags found by find_package (SeqAn).
SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
#SET(COMP_PROF_FLAGS "-ftime-report")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${COMP_PROF_FLAGS}")
set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_BUILD_TYPE Debug)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install linear in ${PREFIX}/bin directory
install (TARGETS linear RUNTIME DESTINATION "${CMAKE_BINARY_DIR}/bin")



