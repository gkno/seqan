cmake_minimum_required (VERSION 2.6)
project (Seqan_Tests CXX)

#set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

# For the tests, we switch SeqAn into testing mode.
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=1")
endif (CMAKE_COMPILER_IS_GNUCXX)


################################################################################
# Set path variables
################################################################################

include_directories (${SEQAN_LIBRARY})
link_directories (${SEQAN_LIBRARY}/lib)
if (MSVC)
        add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Tests
################################################################################

#-----------------------------------------------------------------------------
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked.
INCLUDE(CTest)
ENABLE_TESTING()
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

get_filename_component (SEQAN_TESTS_ABS ${SEQAN_LIBRARY}/../tests ABSOLUTE)
#if (EXISTS SEQAN_TESTS_ABS)
	file(GLOB SEQAN_TESTS ${SEQAN_TESTS_ABS}/[A-z]*)
	foreach (TEST_DIR ${SEQAN_TESTS})
		file (RELATIVE_PATH TEST_DIR_REL ${SEQAN_TESTS_ABS} ${TEST_DIR})
		file (GLOB TEST_SOURCES ${TEST_DIR}/[A-z]*.h ${TEST_DIR}/test_[A-z]*.cpp)
		list (SORT TEST_SOURCES)
		if (TEST_SOURCES)
			add_executable (test_${TEST_DIR_REL} ${TEST_SOURCES})
			add_dependencies (test_${TEST_DIR_REL} ${SEQAN_TARGET})
			add_test (test_${TEST_DIR_REL} test_${TEST_DIR_REL})
			if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
				target_link_libraries (test_${TEST_DIR_REL} rt)
			endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
		endif (TEST_SOURCES)
	endforeach (TEST_DIR)
#endif (EXISTS SEQAN_TESTS_ABS)
