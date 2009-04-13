cmake_minimum_required (VERSION 2.6) 
project (razers)

# Set path variables
set (SEQAN_BASE ../..)
set (TBB_BASE ${SEQAN_BASE}/extra/tbb)

set (CMAKE_BUILD_TYPE Debug)

# Make sure the compiler can find include files from the SeqAn and TBB library. 
include_directories (${SEQAN_BASE} ${TBB_BASE}/include)

# Make sure the compiler can find object files from the TBB library. 
link_directories (${TBB_BASE}/projects/library/lib)

add_executable (razers razers.cpp)
#add_executable (paramChooser paramChooser.cpp)

#target_link_libraries (paramChooser rt) 
target_link_libraries (razers tbb pthread) 

set_source_files_properties(razers.cpp PROPERTIES COMPILE_FLAGS "-W -Wall")
