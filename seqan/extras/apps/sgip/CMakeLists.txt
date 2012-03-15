cmake_minimum_required (VERSION 2.6)
project (extras_apps_sgip)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(sgip sgip.cpp sgip.h sgip_base.h sgip_output.h)

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_sgip COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)