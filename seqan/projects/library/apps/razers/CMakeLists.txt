#make_minimum_required (VERSION 2.6) 
project (razers)

# Set path variables
set (SEQAN_BASE ../..)
set (SEQAN_BASE_PLATFORM ${SEQAN_BASE}/../../platforms/gcc)
set (TBB_BASE ${SEQAN_BASE}/extra/tbb)

# Make sure the compiler can find include files from the SeqAn and TBB library. 
include_directories (${SEQAN_BASE} ${SEQAN_BASE_PLATFORM} ${TBB_BASE}/include)

# Make sure the compiler can find object files from the TBB library. 
link_directories (${TBB_BASE}/build/macos_em64t_gcc_cc4.0.1_os10.5.5_debug)

add_executable (razers razers.cpp)
#add_executable (paramChooser paramChooser.cpp)

#target_link_libraries (paramChooser rt) 
target_link_libraries (razers tbb) 
