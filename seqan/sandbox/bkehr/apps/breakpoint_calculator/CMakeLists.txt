cmake_minimum_required (VERSION 2.6)
project (sandbox_bkehr_apps_breakpoint_calculator)

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_breakpoint_calculator COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

# Lemon, in case we compile on Linux.
include_directories($ENV{HOME}/workspace/LEMON/include)
# Lemon, in case we compile on Windows.
include_directories(Z:/workspace/LEMON/include)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(breakpoint_calculator breakpoint_calculator.cpp
                                           breakpoint_calculator.h
                                           parse_alignment.h
                                           breakpoint_counts.h)
