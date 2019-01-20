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

# Update the list of file names below if you add source files to your application.
add_executable (linear mapper.cpp mapper.h base.h index_extend.h shape_extend.h mapparm.h gap.h align_interface.h f_io.h)
add_executable (testCase testCase.cpp mapper.h base.h index_extend.h shape_extend.h mapparm.h)
add_library(linear_fio STATIC f_io.cpp)
# Add dependencies found by find_package (SeqAn).
target_link_libraries (linear ${SEQAN_LIBRARIES} linear_fio)
target_link_libraries (testCase ${SEQAN_LIBRARIES})


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

#set (CMAKE_BUILD_TYPE Release)
set (CMAKE_BUILD_TYPE Debug)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install linear in ${PREFIX}/bin directory
install (TARGETS linear RUNTIME DESTINATION "${CMAKE_BINARY_DIR}/bin")



