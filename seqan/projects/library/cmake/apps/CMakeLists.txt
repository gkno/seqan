cmake_minimum_required (VERSION 2.6) 
project (Seqan_Apps CXX)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)
set (TBB_BASE ${SEQAN_BASE}/extra/tbb)

add_subdirectory (../seqan seqan)

################################################################################
# Set path variables
################################################################################

# Make sure the compiler can find include files from the SeqAn and TBB library.
include_directories (${SEQAN_BASE} ${TBB_BASE}/include)
link_directories (${SEQAN_BASE}/lib)
if (MSVC)
        add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

# Make sure the compiler can find object files from the TBB library.
link_directories (${SEQAN_BASE}/lib)

if (MSVC)
        add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Applications
################################################################################

# Find application header files
file (GLOB RAZERS_HEADERS ${SEQAN_BASE}/apps/razers/[A-z]*.h)
file (GLOB RAZERS2_HEADERS ${SEQAN_BASE}/apps/razers2/[A-z]*.h)
file (GLOB SEQAN_TCOFFEE_HEADERS ${SEQAN_BASE}/apps/seqan_tcoffee/[A-z]*.h)
list (SORT RAZERS_HEADERS)
list (SORT SEQAN_TCOFFEE_HEADERS)

# Make sure the compiler can find object files from the TBB library. 
link_directories (${SEQAN_BASE}/lib)

add_executable (dfi ${SEQAN_BASE}/apps/dfi/dfi.cpp)
add_executable (paramChooser ${SEQAN_BASE}/apps/razers/paramChooser.cpp ${RAZERS_HEADERS})
add_executable (razers ${SEQAN_BASE}/apps/razers/razers.cpp ${RAZERS_HEADERS})
add_executable (razers2 ${SEQAN_BASE}/apps/razers2/razers.cpp ${RAZERS2_HEADERS})
add_executable (seqan_tcoffee ${SEQAN_BASE}/apps/seqan_tcoffee/seqan_tcoffee.cpp ${SEQAN_TCOFFEE_HEADERS})
add_executable (seqcons ${SEQAN_BASE}/apps/seqcons/seqcons.cpp)
add_executable (pair_align ${SEQAN_BASE}/apps/pair_align/pair_align.cpp)
add_executable (tree_recon ${SEQAN_BASE}/apps/tree_recon/tree_recon.cpp)
#add_executable (transquant ${SEQAN_BASE}/apps/transquant/transquant.cpp)

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	target_link_libraries (razers rt)
	target_link_libraries (paramChooser rt)
	target_link_libraries (seqan_tcoffee rt)
	target_link_libraries (seqcons rt)
	#target_link_libraries (razers tbb_debug pthread)
	#target_link_libraries (razers tbb_debug)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
