cmake_minimum_required (VERSION 2.6)
project (core_apps_razers)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(razers razers.cpp
			    razers.h
			    outputFormat.h
			    paramChooser.h
			    razers_matepairs.h
			    razers_spliced.h
			    readSimulator.h)

# Install pair_align in /bin directory
install(TARGETS razers
        DESTINATION bin)

# Run app tests if python was found
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_razers COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  message(STATUS "add_test(NAME app_test_pair_align COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})")
endif(PYTHONINTERP_FOUND)