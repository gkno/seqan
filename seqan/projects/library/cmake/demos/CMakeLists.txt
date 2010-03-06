cmake_minimum_required (VERSION 2.6) 
project (Seqan_Demos)

################################################################################
# For the demos, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
################################################################################

set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0")

################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)

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
