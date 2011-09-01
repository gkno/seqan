cmake_minimum_required (VERSION 2.6)

# ===========================================================================
# Iff the file INST exists in the source root
# turn on SeqAn instrumentation for experimental data collection.
# ===========================================================================

if(EXISTS ${CMAKE_SOURCE_DIR}/INST)
  set(SEQAN_INSTRUMENTATION YES)
else(EXISTS ${CMAKE_SOURCE_DIR}/INST)
  set(SEQAN_INSTRUMENTATION NO)
endif(EXISTS ${CMAKE_SOURCE_DIR}/INST)


# ===========================================================================
# Warn People About the Perils Of In-Source Builds
# ===========================================================================

if (EXISTS "$ENV{PWD}/CMakeLists.txt" OR EXISTS "$ENV{CD}/CMakeLists.txt")
  if (EXISTS "$ENV{PWD}/I_WANT_IN_SOURCE_BUILDS" OR EXISTS "$ENV{CD}/I_WANT_IN_SOURCE_BUILDS")
    set(IN_SOURCE_BUILD YES)
  else (EXISTS "$ENV{PWD}/I_WANT_IN_SOURCE_BUILDS" OR EXISTS "$ENV{CD}/I_WANT_IN_SOURCE_BUILDS")
    message("CAUTION: In-source builds are not allowed by default!")
    message("")
    message("(If you want to use Eclipse and have read the 'Getting Started,")
    message("with Eclipse Tutorial' then see the instructions below.)")
    message("")
    message("Instead of doing in source builds, go to into the directory build")
    message("and call CMake from there.")
    message("We recommend you to create an additional directory level, this")
    message("allows you to have multiple build configurations (e.g. Debug,")
    message("Release) or project types.")
    message("")
    message("  cd build")
    message("  mkdir Debug")
    message("  cd Debug")
    message("  cmake ../.. -DCMAKE_BUILD_TYPE=Debug")
    message("")
    message("NOTE: CMake has created some files in this directory that you have")
    message("to remove before executing the above commands.")
    message("")
    message("  [Linux, Mac Os X]  rm -rf CMakeCache.txt CMakeFiles")
    message("  [Windows]          del CMakeCache.txt ; rmdir /S CMakeFiles")
    message("")
    message("If you want to use Eclipse then you might want to have in-source")
    message("builds because of limitations of Eclipse CDT4.  In this case")
    message("Create a file in this directory called I_WANT_IN_SOURCE_BUILDS:")
    message("")
    message("  [Linux, Mac Os X]  touch I_WANT_IN_SOURCE_BUILDS")
    message("  [Windows]          type nul >> I_WANT_IN_SOURCE_BUILDS")
    message("")
    message(FATAL_ERROR "Please follow the instructions above.")
  endif (EXISTS "$ENV{PWD}/I_WANT_IN_SOURCE_BUILDS" OR EXISTS "$ENV{CD}/I_WANT_IN_SOURCE_BUILDS")
endif (EXISTS "$ENV{PWD}/CMakeLists.txt" OR EXISTS "$ENV{CD}/CMakeLists.txt")

cmake_minimum_required (VERSION 2.6)
project (seqan)

# ===========================================================================
# Setup Variables
# ===========================================================================

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Debug CACHE STRING
         "Choose the type of build, options are: Debug Release RelDebug ..."
         FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# ===========================================================================
# Setup Modules and Find Packages.
# ===========================================================================

set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/FindTBB
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/cuda
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake
    ${CMAKE_MODULE_PATH})

include(CTest)

find_package (PythonInterp)
find_package (SeqAn)

# ===========================================================================
# Prerequisites
# ===========================================================================

seqan_find_dependencies()  # Find SeqAn optional dependencies.
seqan_setup_global()

# ===========================================================================
# Include Subdirectories
# ===========================================================================

add_subdirectory(core)
add_subdirectory(extras)
add_subdirectory(sandbox)
