# ===========================================================================
# Linear
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (linear CXX)

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include/
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
endif ()

# ----------------------------------------------------------------------------
# Add linear targets
# ----------------------------------------------------------------------------

add_subdirectory(src)

