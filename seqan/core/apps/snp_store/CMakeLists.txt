cmake_minimum_required (VERSION 2.6)
project (core_apps_snp_store)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../extras/include)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(snp_store snp_store.cpp snp_store.h)

install(TARGETS snp_store
        DESTINATION bin)
