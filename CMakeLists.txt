# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/pacMapper/CMakeLists.txt
#
# CMakeLists.txt file for INDEX_BENCH.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_pacMapper)
message (STATUS "Configuring apps/pacMapper")

set (SEQAN_APP_VERSION "1.1.1")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
find_package(CXX11 REQUIRED)
if (NOT CXX11_FOUND)
    message (STATUS " C++11 compiler needed to build pacMapper")
    return ()
endif (NOT CXX11_FOUND)

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
add_executable (pacMapper pacMapper.cpp base.h index_extend.h shape_extend.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (pacMapper ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (pacMapper)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install pacMapper in ${PREFIX}/bin directory
install (TARGETS pacMapper
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/pacMapper
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/small.fasta
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (pacMapper)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqan_tcoffee in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} pacMapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:pacMapper")
  set (CPACK_PACKAGE_NAME "pacMapper")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "pacMapper - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Chenxu Pan <chenxu.pan@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (pacMapper "pacMapper")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:pacMapper")
