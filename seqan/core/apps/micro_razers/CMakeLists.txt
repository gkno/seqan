cmake_minimum_required (VERSION 2.6)
project (core_apps_micro_razers)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(micro_razers micro_razers.cpp outputFormat.h)
