cmake_minimum_required (VERSION 2.6)
project (core_apps_mason)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(mason mason.cpp
			   mason.h
			   simulate_454.h
			   simulate_454_base_calling.h
			   simulate_illumina.h
			   simulate_sanger.h
			   store_config.h
			   util.h)