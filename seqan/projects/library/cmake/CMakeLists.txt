cmake_minimum_required (VERSION 2.6) 
project (Seqan CXX)

################################################################################
# Set path variables
################################################################################

if (NOT SEQAN_LIBRARY)
	get_filename_component (SEQAN_LIBRARY_ABS .. ABSOLUTE)
	set (SEQAN_LIBRARY ${SEQAN_LIBRARY_ABS} CACHE FILEPATH "Path to SeqAn's parent folder (contains apps,demos,seqan)")
endif (NOT SEQAN_LIBRARY)

set (SEQAN_TARGET "Seqan")

################################################################################
# Library sources
################################################################################

# Find all library headers
get_filename_component (SEQAN_BASE_ABS ${SEQAN_LIBRARY}/seqan ABSOLUTE)
file (GLOB SEQAN_HEADERS ${SEQAN_BASE_ABS}/[A-z]*/[A-z]*.h)
list (SORT SEQAN_HEADERS)

#foreach (SEQAN_HEADER ${SEQAN_HEADERS})
#	file (RELATIVE_PATH SEQAN_HEADER_REL ${SEQAN_BASE_ABS} ${SEQAN_HEADER})
#	get_filename_component (SEQAN_MODULE ${SEQAN_HEADER_REL} PATH)
#	get_filename_component (FILENAME ${SEQAN_HEADER_REL} NAME)
# file(APPEND filename ${FILENAME} "\n")
# file(APPEND filename ${SEQAN_MODULE}_generated_forwards.h "\n")
# if (${FILENAME} STREQUAL ${SEQAN_MODULE}_generated_forwards.h)
#   list(REMOVE_ITEM SEQAN_HEADERS ${SEQAN_HEADER})
#   list(APPEND SEQAN_FORWARDS ${SEQAN_HEADER})
# endif (${FILENAME} STREQUAL ${SEQAN_MODULE}_generated_forwards.h)
#endforeach (SEQAN_HEADER)

# Generate forwards for g++ version 4.0 and greater
if (CMAKE_COMPILER_IS_GNUCXX)
  set (SEQAN_MAKE_FORWARDS COMMAND python ${SEQAN_LIBRARY}/../../misc/build_forwards.py ${SEQAN_BASE_ABS})
endif (CMAKE_COMPILER_IS_GNUCXX)

# For Non-GUI Makefiles there is a bug in cmake
# the SOURCE command in add_custom_target is not
# recognized there
if (NOT ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles")
  set (SEQAN_GUI_SOURCES SOURCES ${SEQAN_HEADERS})
endif (NOT ${CMAKE_GENERATOR} STREQUAL "Unix Makefiles")

add_custom_target (${SEQAN_TARGET}
  ${SEQAN_MAKE_FORWARDS}
  DEPENDS ${SEQAN_HEADERS}
  ${SEQAN_GUI_SOURCES}
)

# Group library headers into modules
foreach (SEQAN_HEADER ${SEQAN_HEADERS})
	file (RELATIVE_PATH SEQAN_HEADER_REL ${SEQAN_BASE_ABS} ${SEQAN_HEADER})
	get_filename_component (SEQAN_MODULE ${SEQAN_HEADER_REL} PATH)
	source_group (${SEQAN_MODULE} FILES ${SEQAN_HEADER})
endforeach (SEQAN_HEADER)

################################################################################
# Set compiler flags
################################################################################

# Make sure the compiler can find include files from the SeqAn library. 
include_directories (${SEQAN_LIBRARY})

if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_WARNING_LEVEL 4)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
	add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Include sub-projects
################################################################################

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

add_subdirectory (apps)
add_subdirectory (demos)
add_subdirectory (tests)
