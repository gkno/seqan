cmake_minimum_required (VERSION 2.6)
project (core_apps_splazers)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(splazers splazers.cpp
			      outputFormat.h
			      paramChooser.h
			      razers.h
			      razers_matepairs.h
			      razers_parallel.h
			      razers_spliced.h
			      readSimulator.h)
