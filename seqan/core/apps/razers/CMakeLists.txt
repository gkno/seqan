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
