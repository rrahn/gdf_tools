# ===========================================================================
#                  gdf_tools
# ===========================================================================

# ----------------------------------------------------------------------------
# Load SeqAn build system
# ----------------------------------------------------------------------------

project (seqan)
cmake_minimum_required (VERSION 2.8.2)
set (CMAKE_INCLUDE_PATH
${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include 
${CMAKE_CURRENT_SOURCE_DIR}/include/dev
${CMAKE_INCLUDE_PATH}
)
set (CMAKE_MODULE_PATH
${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake/FindTBB
${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
${CMAKE_MODULE_PATH})
include (SeqAnContribs)
set (SEQAN_USE_SEQAN_BUILD_SYSTEM TRUE CACHE INTERNAL "Use SeqAn build system." FORCE)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/" CACHE INTERNAL "Root source directory." FORCE)
# include(include/seqan/CMakeLists.txt)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)
include (SeqAnCtdSetup)

# ----------------------------------------------------------------------------
# Load Jst bench code.
# ----------------------------------------------------------------------------

add_subdirectory(src)