cmake_minimum_required (VERSION 2.6) 
project (Seqan_Demos CXX)

# For the demos, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE "-DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0")
  set (CMAKE_CXX_FLAGS_DEBUG "-DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0")
endif (CMAKE_COMPILER_IS_GNUCXX)

#set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)
if (MSVC)
	add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Demos
################################################################################

file(GLOB SEQAN_DEMOS ${SEQAN_LIBRARY}/demos/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	add_executable (${DEMO} ${DEMOFILE})
	add_dependencies (${DEMO} ${SEQAN_TARGET})
	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		target_link_libraries (${DEMO} rt)
	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
endforeach (DEMOFILE)
