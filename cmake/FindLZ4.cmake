# Finds liblz4.
#
# This module defines:
# LZ4_FOUND
# LZ4_INCLUDE_DIR
# LZ4_LIBRARY
#

find_path(LZ4_INCLUDE_DIR NAMES lz4.h)
find_library(LZ4_LIBRARY NAMES lz4)

# We require LZ4_compress_default() which was added in v1.7.0
if( LZ4_LIBRARY )
  include( CheckCSourceRuns )
  set( CMAKE_REQUIRED_INCLUDES ${LZ4_INCLUDE_DIR} )
  set (CMAKE_REQUIRED_LIBRARIES ${LZ4_LIBRARY} )

  set( CMAKE_REQUIRED_INCLUDES )
  set( CMAKE_REQUIRED_LIBRARIES )
endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    LZ4 DEFAULT_MSG
    LZ4_LIBRARY LZ4_INCLUDE_DIR )

if( LZ4_FOUND )
  message(STATUS "FindLZ4: Found both LZ4 headers and library")
else( LZ4_FOUND )
  if( LZ4_FIND_REQUIRED )
    message( FATAL_ERROR "FindLZ4: Could not find LZ4 headers or library" )
  endif( LZ4_FIND_REQUIRED )
endif( LZ4_FOUND )

mark_as_advanced(LZ4_INCLUDE_DIR LZ4_LIBRARY)