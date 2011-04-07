# ===========================================================================
# FindSeqAn.cmake -- Utilities for building SeqAn applications.
# ===========================================================================

# ---------------------------------------------------------------------------
# Macro seqan_setup_global ()
# ---------------------------------------------------------------------------

# Global setup for SeqAn.
#
# This consists of setting the warning level and suppressing some warnings.

macro (seqan_setup_global)
    # This is used for calling anything in util.
    set (SEQAN_ROOT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL
         "Used to get the directory to util, for example." FORCE)

    # -----------------------------------------------------------------------
    # GCC Setup
    # -----------------------------------------------------------------------
    if (CMAKE_COMPILER_IS_GNUCXX)
        # For the GCC, enable warnings.
        set (CMAKE_CXX_WARNING_LEVEL 4)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wno-long-long -fstrict-aliasing -Wstrict-aliasing")

        # Determine GCC version.
        # message("Determining GCC version.")
        EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
                     ARGS --version
                     OUTPUT_VARIABLE _GCC_VERSION)
        STRING(REGEX REPLACE ".* ([0-9])\\.([0-9])\\.([0-9]) .*" "\\1\\2\\3"
               _GCC_VERSION ${_GCC_VERSION})
        # message("  GCC version is ${_GCC_VERSION}")

        # Add -Wno-longlong if the GCC version is < 4.0.0.  Add -pedantic flag
        # but disable warnings for variadic macros with GCC >= 4.0.0.  Earlier
        # versions warn because of anonymous variadic macros in pedantic mode
        # but do not have a flag to disable these warnings.
        if (400 GREATER _GCC_VERSION)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long")
        else (400 GREATER _GCC_VERSION)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wno-variadic-macros")
        endif (400 GREATER _GCC_VERSION)

        # Force GCC to keep the frame pointer when debugging is enabled.
        # This is mainly important for 64 bit but does not get into the way
        # on 32 bit either at minimal performance impact.
        set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
        set (CMAKE_CXX_FLAGS_RELDEBUG "${CMAKE_CXX_FLAGS_RELEASE} -g")# -fno-omit-frame-pointer")
  
        # Pass CXX flags to flags.
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_CXX_FLAGS_=\"${CMAKE_CXX_FLAGS}\"")
    endif (CMAKE_COMPILER_IS_GNUCXX)

    # -----------------------------------------------------------------------
    # Visual Studio Setup
    # -----------------------------------------------------------------------
    if (MSVC)
        # Warning level 3 for MSVC is disabled for now to see how much really bad warnings there are.
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3  /wd4996 -D_CRT_SECURE_NO_WARNINGS")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2  /wd4996 -D_CRT_SECURE_NO_WARNINGS")
        # Disable Microsoft C++ language extensions.
        # TODO(holtgrew): Re-disable again, Microsoft's header do not compile with this option! Look for a workaround.
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Za")
    endif (MSVC)
endmacro (seqan_setup_global)

# ---------------------------------------------------------------------------
# Macro seqan_setup_includes ()
# ---------------------------------------------------------------------------

# Setup an "include" directory, i.e. add it to the include directories via
# add_includes() and register forward building for everything inside a
# subdirectory "seqan".  Also registers it in variable SEQAN_LIBRARY_TARGETS.

function (seqan_setup_includes REL_PATH TARGET_NAME)
    set (PATH ${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH})
    file (GLOB HEADERS_TMP ${PATH}/seqan/[A-z]*/[A-z]*.h)
    file (GLOB SUPER_HEADERS ${PATH}/seqan/[A-z]*.h)

    set (SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH} CACHE INTERNAL "asdf" FORCE)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${REL_PATH})

    set (SEQAN_LIBRARY_TARGETS ${TARGET_NAME} CACHE INTERNAL "asdf" FORCE)

    foreach (HEADER ${HEADERS_TMP})
        if (NOT ${HEADER} MATCHES ".*generated.*")
            list (APPEND HEADERS ${HEADER})
        else (NOT ${HEADER} MATCHES ".*generated.*")
            list (APPEND FORWARDS ${HEADER})
        endif (NOT ${HEADER} MATCHES ".*generated.*")
    endforeach (HEADER ${HEADERS})
    
    # Sort headers and forwards.
    if (HEADERS)
        list (SORT HEADERS)
    endif (HEADERS)
    if (FORWARDS)
        list (SORT FORWARDS)
    endif (FORWARDS)

    # ---------------------------------------------------------------------------
    # Forwards Generation.
    # ---------------------------------------------------------------------------
    if (CMAKE_COMPILER_IS_GNUCXX)
        file (GLOB BASE_CONTENTS RELATIVE ${PATH}/seqan ${PATH}/seqan/[A-z]*)
        foreach (ENTRY ${BASE_CONTENTS})
            if (IS_DIRECTORY ${PATH}/seqan/${ENTRY})
                list (APPEND MODULES ${ENTRY})
            endif (IS_DIRECTORY ${PATH}/seqan/${ENTRY})
        endforeach (ENTRY ${SEQAN_BASE_CONTENTS})
        
        if (MODULES)
            list (SORT MODULES)
            list (REMOVE_DUPLICATES MODULES)
        
            # Build a list of generated forwards headers.  Goes into SEQAN_FORWARDS.
            foreach (MODULE ${MODULES})
                list (APPEND FORWARDS ${PATH}/seqan/${MODULE}/${MODULE}_generated_forwards.h)
            endforeach (MODULE ${MODULES})

            # Now tell CMake that the forward headers can be generated with
            # build_forwards.py
            add_custom_command (
                OUTPUT ${FORWARDS}
                COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/util/bin/build_forwards.py ${PATH}/seqan all
                DEPENDS ${HEADERS})
            #message(STATUS "add_custom_command (
            #    OUTPUT ${FORWARDS}
            #    COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_ROOT_SOURCE_DIR}/util/bin/build_forwards.py ${PATH}/seqan all
            #    DEPENDS ${HEADERS})")
        endif (MODULES)
    endif (CMAKE_COMPILER_IS_GNUCXX)

    # ---------------------------------------------------------------------------
    # GUI Setup Stuff.
    # ---------------------------------------------------------------------------
    get_filename_component (BASE_ABS ${PATH}/seqan ABSOLUTE)
    
    # CMake bug workaround: For Non-GUI generators there is a bug in cmake.
    # The SOURCE command in add_custom_target is not recognized there.
    set (NONGUI_GENERATORS "Unix Makefiles" "MinGW Makefiles")
    list (FIND NONGUI_GENERATORS ${CMAKE_GENERATOR} FOUND)
    if (FOUND EQUAL -1)
        set (GUI_SOURCES SOURCES ${HEADERS} ${SUPER_HEADERS})
    endif (FOUND EQUAL -1)

    # Add SeqAn Pseudo Target for GUIs.
    #
    # This target contains all headers, forwards and umbrella headers.
    add_custom_target(
        ${TARGET_NAME}
        DEPENDS ${HEADERS}
                ${FORWARDS}
                ${GUI_SOURCES}
    )
    # message("add_custom_target(${TARGET_NAME} DEPENDS ${HEADERS} ${FORWARDS} ${GUI_SOURCES})")
    
    # Group library headers into modules.  The CMake documentation says this
    # is mostly (only?) used for Visual Studio Projects.
    foreach (HEADER ${HEADERS})
        file (RELATIVE_PATH HEADER_REL ${BASE_ABS} ${HEADER})
        get_filename_component (MODULE ${HEADER_REL} PATH)
        source_group (${MODULE} FILES ${HEADER})
        # message("source_group(${MODULE} FILES ${HEADER})")    
    endforeach (HEADER ${HEADERS})
endfunction (seqan_setup_includes)

# ---------------------------------------------------------------------------
# Macro seqan_make_seqan_available ()
# ---------------------------------------------------------------------------

# Register a seqan extension that has been previously defined available.
# TODO(holtgrew): Is order really important at all?

function (seqan_make_seqan_available TARGET_NAME)
   set (SEQAN_LIBRARY_TARGETS "${SEQAN_LIBRARY_TARGETS} ${TARGET_NAME}" CACHE INTERNAL "asdf" FORCE)
   include_directories(${SEQAN_INCLUDE_DIR_FOR_${TARGET_NAME}})
endfunction (seqan_make_seqan_available TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_tests ()
# ---------------------------------------------------------------------------

# Switch to testing mode.  This function should be called in the
# CMakeLists.txt in the tests directories before including subdirectories.
#
# The following will happen:
#
# * Setting definitions SEQAN_ENABLE_DEBUG=1 and SEQAN_ENABLE_TESTING=1.
# * If the ${MODEL} variable is NightlyCoverage or ExperimentalCoverage,
#   and the compiler is GCC C++ then symbols for test coverate are added.

macro (seqan_setup_tests TEST_TARGET)
    # Setup flags for tests.
    add_definitions(-DSEQAN_ENABLE_DEBUG=1)
    add_definitions(-DSEQAN_ENABLE_TESTING=1)
    
    # Add a target for the tests.
    add_custom_target(${TEST_TARGET})
    # Create a CMake variable for storing the current test target.
    set (SEQAN_CURRENT_TEST_TARGET ${TEST_TARGET} CACHE INTERNAL
         "Test target, communicated to seqan_add_test_executable" FORCE)
    # message (STATUS "TEST_TARGET " ${TEST_TARGET})
    # message (STATUS "SEQAN_CURRENT_TEST_TARGET " ${SEQAN_CURRENT_TEST_TARGET})
    
    # Conditionally enable coverage mode.
    if (MODEL STREQUAL "NightlyCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
        endif (CMAKE_COMPILER_IS_GNUCXX)
    endif (MODEL STREQUAL "NightlyCoverage")
    if (MODEL STREQUAL "ExperimentalCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
        endif (CMAKE_COMPILER_IS_GNUCXX)
    endif (MODEL STREQUAL "ExperimentalCoverage")
endmacro (seqan_setup_tests)

# ---------------------------------------------------------------------------
# Macro seqan_find_dependencies ()
# ---------------------------------------------------------------------------

# Try to find all dependencies using the find() function.
#
# Currently the libraries SeqAn can use to extend its functionality are:
#
#  * zlib
#  * bzlib
#  * OpenMP
#  * CUDA

macro (seqan_find_dependencies)
    find_package(ZLIB QUIET)
    find_package(BZip2 QUIET)
    find_package(OpenMP QUIET)
    find_package (CUDA)
endmacro (seqan_find_dependencies)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_subdirectories ()
# ---------------------------------------------------------------------------

# This macro calls add_subdirectory() for all subdirectories below
# ${SEQAN_CURRENT_SOURCE_DIR} if they contain a CMakeLists.txt.
#
# Example:
#
#   seqan_add_all_subdirectories(seqan_tests)

macro (seqan_add_all_subdirectories)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/*)
    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            # message(STATUS "Going into ${ENTRY}")
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory(${ENTRY})
            endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
        endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_add_all_subdirectories)

# ---------------------------------------------------------------------------
# Macro seqan_project_just_include_all_subdirs ()
# ---------------------------------------------------------------------------

# Call this macro as the only command in a CMakeLists.txt file to include
# all subdirectories.
#
# Args:
#   PROJECT_NAME  String, the project name of this subdirectory.
#
# Example:
#
#   # Only line in file tests/CMakeLists.txt
#   seqan_project_just_include_all_subdirs(seqan_tests)

macro (seqan_project_just_include_all_subdirs PROJECT_NAME)
    cmake_minimum_required (VERSION 2.6)
    project(${PROJECT_NAME})
    seqan_add_all_subdirectories()
endmacro (seqan_project_just_include_all_subdirs)

# ---------------------------------------------------------------------------
# Macro seqan_add_executable (TARGET_NAME source1.cpp source2.[cpp|h] ...)
#       seqan_add_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn executable from the given source files.  If no such files
# are given then all files in the current directory will be used.

macro (seqan_add_executable TARGET_NAME)
    # TODO(holtgrew): Use all files in directory if ${ARGC} == 0
    # message(STATUS "add_executable (${ARGV})")
    add_executable (${ARGV})
    
    # Link against librt on Linux.
  	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    	target_link_libraries (${TARGET_NAME} rt)
  	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  	
  	# Dependencies on all registered seqan extensions.
  	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})
    # message(STATUS "add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})")
	
  	# Link against zlib and bzlib if found.
  	if (ZLIB_FOUND)
  	    include_directories (${ZLIB_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
  	endif (ZLIB_FOUND)
  	if (BZIP2_FOUND)
  	    include_directories (${BZIP2_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
  	endif (BZIP2_FOUND)
endmacro (seqan_add_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_test_executable (TARGET_NAME source1.cpp source2.[cpp|h])
#       seqan_add_test_executable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn executable from the given source files.  If no such files
# are given then all files in the current directory will be used.
#
# Also the test will be registered with add_test and depend on the current
# test target, i.e. from closest seqan_setup_tests(TEST_TARGET_NAME) call.

macro (seqan_add_test_executable TARGET_NAME)
    # TODO(holtgrew): Use all files in directory if ${ARGC} == 0
    # message(STATUS "add_executable (${ARGV})")
    add_executable (${ARGV})
    add_test(test_${ARGV0} ${ARGV0})
                     
    add_dependencies(${SEQAN_CURRENT_TEST_TARGET} ${ARGV0})
    
    # Link against librt on Linux.
  	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  		  target_link_libraries (${ARGV0} rt)
  	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

  	# Dependencies on all registered seqan extensions.
    # message("add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})")
  	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})

  	# Link against zlib and bzlib if found.
  	if (ZLIB_FOUND)
  	    include_directories (${ZLIB_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
  	endif (ZLIB_FOUND)
  	if (BZIP2_FOUND)
  	    include_directories (${BZIP2_INCLUDE_DIR})
  	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
  	endif (BZIP2_FOUND)
endmacro (seqan_add_test_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Function seqan_add_cuda_testexecutable (TARGET_NAME source1.cu source2.[cu|cpp|h])
#          seqan_add_cuda_testexecutable (TARGET_NAME)
# ---------------------------------------------------------------------------

# Create a SeqAn CUDA executable from the given source files.  If no such
# files are given then all files in the current directory will be used.
#
# Also the test will be registered with add_test and depend on the current
# test target, i.e. from closest seqan_setup_tests(TEST_TARGET_NAME) call.
#
# Only adds target if CUDA_FOUND is set to true.

function (seqan_add_cuda_executable TARGET_NAME)
    if (CUDA_FOUND)
        # TODO(holtgrew): Use all files in directory if ${ARGC} == 0

        # -------------------------------------------------------------------
        # Set CUDA variables
        # -------------------------------------------------------------------
        set (CUDA_PROPAGATE_HOST_FLAGS OFF)
        set (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
        #set (CUDA_NVCC_FLAGS "-pg")
        list (APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS "cu")
        cuda_include_directories(${CUDA_CUT_INCLUDE_DIR})
        #string (REGEX REPLACE "\\-(W( |all)|pedantic)" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})
        string (REGEX REPLACE "\\-pedantic" "" CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})

        # -------------------------------------------------------------------
        # Go on normally as in apps.
        # -------------------------------------------------------------------
        # message(STATUS "cuda_add_executable (${ARGV})")
        cuda_add_executable (${ARGV})
    
        # Link against librt on Linux.
      	if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        	target_link_libraries (${TARGET_NAME} rt)
      	endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  	
      	# Dependencies on all registered seqan extensions.
      	add_dependencies(${ARGV0} ${SEQAN_LIBRARY_TARGETS})
	
      	# Link against zlib and bzlib if found.
      	if (ZLIB_FOUND)
      	    include_directories (${ZLIB_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${ZLIB_LIBRARIES})
      	endif (ZLIB_FOUND)
      	if (BZIP2_FOUND)
      	    include_directories (${BZIP2_INCLUDE_DIR})
      	    target_link_libraries (${TARGET_NAME} ${BZIP2_LIBRARIES})
      	endif (BZIP2_FOUND)
    endif (CUDA_FOUND)
endfunction (seqan_add_cuda_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_executables ([TARGET])
# ---------------------------------------------------------------------------

# This macro calls seqan_add_executable() for all all .cpp files in the
# current directory.  If the optional TARGET parameter is given, a target
# with this name is created and the executable targets will depend on this
# target.
#
# Example:
#
#   seqan_add_all_executables()
#   seqan_add_all_executables(depend_on_this)

macro (seqan_add_all_executables)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    if (${ARGC} GREATER 0)
        if (TARGET ${ARGV0})  # Add target only if it does not exist yet.
        else (TARGET ${ARGV0})
            add_custom_target(${ARGV0})
        endif (TARGET ${ARGV0})
    endif (${ARGC} GREATER 0)
    foreach (ENTRY ${ENTRIES})
        get_filename_component(BIN_NAME ${ENTRY} NAME_WE)
        seqan_add_executable(${BIN_NAME} ${ENTRY})
        if (${ARGC} GREATER 0)
            add_dependencies(${ARGV0} ${BIN_NAME})
        endif (${ARGC} GREATER 0)
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_add_all_executables)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_cuda_executables ([TARGET])
# ---------------------------------------------------------------------------

# This macro calls seqan_add_cuda_executable() for all all .cu files in the
# current directory.  If the optional TARGET parameter is given, a target
# with this name is created and the executable targets will depend on this
# target.
#
# Only adds target if CUDA_FOUND is set to true.
#
# Example:
#
#   seqan_add_all_cuda_executables()
#   seqan_add_all_cuda_executables(depend_on_this)

macro (seqan_add_all_cuda_executables)
    if (CUDA_FOUND)
        file (GLOB ENTRIES
              RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
              ${CMAKE_CURRENT_SOURCE_DIR}/*.cu)
        if (${ARGC} GREATER 0)
            if (TARGET ${ARGV0})  # Add target only if it does not exist yet.
            else (TARGET ${ARGV0})
                add_custom_target(${ARGV0})
            endif (TARGET ${ARGV0})
        endif (${ARGC} GREATER 0)
        foreach (ENTRY ${ENTRIES})
            get_filename_component(BIN_NAME ${ENTRY} NAME_WE)
            # message("seqan_add_cuda_executable(${BIN_NAME} ${ENTRY})")
            seqan_add_cuda_executable(${BIN_NAME} ${ENTRY})
            if (${ARGC} GREATER 0)
                add_dependencies(${ARGV0} ${BIN_NAME})
            endif (${ARGC} GREATER 0)
        endforeach (ENTRY ${ENTRIES})
    endif (CUDA_FOUND)
endmacro (seqan_add_all_cuda_executables)