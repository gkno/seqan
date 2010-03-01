cmake_minimum_required (VERSION 2.6) 
project (Seqan_Apps CXX)

# For the apps, use assertions in debug mode, no assertions in
# release mode and never go into testing mode.
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_RELEASE "-DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0")
  set (CMAKE_CXX_FLAGS_DEBUG "-DSEQAN_ENABLE_DEBUG=1 -DSEQAN_ENABLE_TESTING=0")
endif (CMAKE_COMPILER_IS_GNUCXX)

#set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

# Use CMake's package location functionality to find the directories of
# TBB and Boost.
find_package (TBB)
find_package (Boost)

################################################################################
# Set path variables
################################################################################

# Make sure the compiler can find include files from the SeqAn, Boost and TBB library.
include_directories (${SEQAN_LIBRARY}
                     ${TBB_INCLUDE_DIRS}
                     ${Boost_INCLUDE_DIRS})
link_directories (${SEQAN_LIBRARY}/lib
                  ${TBB_LIBRARY_DIRS})
if (MSVC)
        add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

# Make sure the compiler can find object files from the TBB library.
link_directories (${SEQAN_LIBRARY}/lib)

if (MSVC)
        add_definitions (/W2 /wd4996 /D "_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

################################################################################
# Applications
################################################################################

# TODO(holtgrew): Disable apps if dependencies are missing.  This would mean removing the automatic generation of targets but this invites errors anyway.

# Find targets automatically
get_filename_component (SEQAN_APPS_ABS ${SEQAN_LIBRARY}/apps ABSOLUTE)
file (GLOB SEQAN_APPS ${SEQAN_APPS_ABS}/[A-z]*)

file (GLOB INSEGT_HEADERS ${SEQAN_LIBRARY}/apps/insegt/[A-z]*.h)
file (GLOB RAZERS_HEADERS ${SEQAN_LIBRARY}/apps/razers/[A-z]*.h)
file (GLOB RAZERS2_HEADERS ${SEQAN_LIBRARY}/apps/razers2/[A-z]*.h)
file (GLOB RAZERS3_HEADERS ${SEQAN_LIBRARY}/apps/razers3/[A-z]*.h)
file (GLOB SEQAN_TCOFFEE_HEADERS ${SEQAN_LIBRARY}/apps/seqan_tcoffee/[A-z]*.h)

# Manually add predefined targets
list (APPEND SEQAN_APPS ${SEQAN_LIBRARY}/apps/paramChooser)

# Manually set predefined sets of source files for some targets
set (APP_razers_SOURCES ${SEQAN_LIBRARY}/apps/razers/razers.cpp ${RAZERS_HEADERS})
set (APP_razers2_SOURCES ${SEQAN_LIBRARY}/apps/razers2/razers.cpp ${RAZERS2_HEADERS})
set (APP_razers3_SOURCES ${SEQAN_LIBRARY}/apps/razers3/razers.cpp ${RAZERS3_HEADERS})
set (APP_micro_razers_SOURCES ${SEQAN_LIBRARY}/apps/micro_razers/micro_razers.cpp ${RAZERS_HEADERS})
set (APP_paramChooser_SOURCES ${SEQAN_LIBRARY}/apps/razers/paramChooser.cpp ${RAZERS_HEADERS})
set (APP_insegt_SOURCES ${SEQAN_LIBRARY}/apps/insegt/insegt.cpp ${INSEGT_HEADERS})
if (Boost_FOUND)
  set (APP_transquant_SOURCES ${SEQAN_LIBRARY}/apps/transquant/transquant.cpp)
else (Boost_FOUND)
  message (WARNING "Boost not found, not building transquant!")
endif (Boost_FOUND)
set (APP_rectifier_SOURCES ${SEQAN_LIBRARY}/apps/transquant/rectifier.cpp)

# Make sure the compiler can find object files from the TBB library. 
link_directories (${SEQAN_LIBRARY}/lib)

foreach (APP_DIR ${SEQAN_APPS})
	file (RELATIVE_PATH APP_DIR_REL ${SEQAN_APPS_ABS} ${APP_DIR})
	
	# if there is no predefined set of files for this target
	# find the source files automatically
	if (NOT APP_${APP_DIR_REL}_SOURCES)
		file (GLOB APP_${APP_DIR_REL}_SOURCES ${APP_DIR}/[A-z]*.h ${APP_DIR}/[A-z]*.cpp)
	endif (NOT APP_${APP_DIR_REL}_SOURCES)
	
	# add target, dependencies and linked libraries
	if (APP_${APP_DIR_REL}_SOURCES)
		list (SORT APP_${APP_DIR_REL}_SOURCES)
		add_executable (${APP_DIR_REL} ${APP_${APP_DIR_REL}_SOURCES})
		add_dependencies (${APP_DIR_REL} ${SEQAN_TARGET})
		if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
			target_link_libraries (${APP_DIR_REL} rt)
		endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	endif (APP_${APP_DIR_REL}_SOURCES)
endforeach (APP_DIR)
