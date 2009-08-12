cmake_minimum_required (VERSION 2.6) 
project (Seqan_Demos CXX)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

add_subdirectory (../seqan ../seqan)

################################################################################
# Set path variables
################################################################################

set (SEQAN_BASE ../..)

################################################################################
# Demos
################################################################################

file(GLOB SEQAN_DEMOS ${SEQAN_BASE}/demos/[A-z]*.cpp)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	add_executable (${DEMO} ${DEMOFILE})
endforeach (DEMOFILE)
