cmake_minimum_required (VERSION 2.6)
project (core_apps_stellar)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(stellar stellar.cpp
			     stellar.h
                             stellar_extension.h
			     stellar_output.h
			     stellar_types.h)
