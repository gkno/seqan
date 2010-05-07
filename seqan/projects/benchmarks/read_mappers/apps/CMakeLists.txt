CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
ENABLE_TESTING()

# Full debug and warnings.
set(CMAKE_CXX_FLAGS "-g -W -Wall")
#set(CMAKE_CXX_FLAGS "-g -W -Wall -O3")

# Make sure the compiler can find Seqan.
INCLUDE_DIRECTORIES(../../../library)

# ------------------------------------------------------------
# Program: wit_builder
# ------------------------------------------------------------

# Add wit_builder executable.
ADD_EXECUTABLE(wit_builder wit_builder.cpp)

# Make sure the fixtures validate.  Then, we can simply compare the
# output of our program with the result fixtures.
ADD_TEST(fixtures_should_validate_01 ./validate_wit.py fixtures/01_result.wit)
ADD_TEST(fixtures_should_validate_02 ./validate_wit.py fixtures/02_result.wit)
ADD_TEST(fixtures_should_validate_03 ./validate_wit.py fixtures/03_result.wit)
ADD_TEST(fixtures_should_validate_04 ./validate_wit.py fixtures/04_result.wit)
ADD_TEST(fixtures_should_validate_05 ./validate_wit.py fixtures/05_result.wit)
ADD_TEST(fixtures_should_validate_06 ./validate_wit.py fixtures/06_result.wit)

# Run the program on the fixtures and compare the result to the
# fixture'd result with diff.
ADD_TEST(run_on_fixtures_01_edit ./run_and_diff.sh 0 edit fixtures/01_seq.fasta fixtures/01_golden.sam fixtures/01_result.wit fixtures/01_result.json)
ADD_TEST(run_on_fixtures_02_edit ./run_and_diff.sh 5 edit fixtures/02_seq.fasta fixtures/02_golden.sam fixtures/02_result.wit fixtures/02_result.json)
ADD_TEST(run_on_fixtures_03_edit ./run_and_diff.sh 0 edit fixtures/03_seq.fasta fixtures/03_golden.sam fixtures/03_result.wit fixtures/03_result.json)
ADD_TEST(run_on_fixtures_04_edit ./run_and_diff.sh 5 edit fixtures/04_seq.fasta fixtures/04_golden.sam fixtures/04_result.wit fixtures/04_result.json)
ADD_TEST(run_on_fixtures_05_edit ./run_and_diff.sh 0 edit fixtures/05_seq.fasta fixtures/05_golden.sam fixtures/05_result.wit fixtures/05_result.json)
ADD_TEST(run_on_fixtures_06_edit ./run_and_diff.sh 8 edit fixtures/06_seq.fasta fixtures/06_golden.sam fixtures/06_result.wit fixtures/06_result.json)
ADD_TEST(run_on_fixtures_07_edit ./run_and_diff.sh 4 edit fixtures/07_seq.fasta fixtures/07_golden.sam fixtures/07_result.wit fixtures/07_result.json)

# With hamming distance, the results should be the same as the edit
# distance for the hamming distance.
ADD_TEST(run_on_fixtures_09_hamming ./run_and_diff.sh 0 hamming fixtures/09_seq.fasta fixtures/09_golden.sam fixtures/09_result.wit fixtures/09_result.json)
ADD_TEST(run_on_fixtures_10_hamming ./run_and_diff.sh 5 hamming fixtures/10_seq.fasta fixtures/10_golden.sam fixtures/10_result.wit fixtures/10_result.json)
ADD_TEST(run_on_fixtures_11_hamming ./run_and_diff.sh 0 hamming fixtures/11_seq.fasta fixtures/11_golden.sam fixtures/11_result.wit fixtures/11_result.json)
ADD_TEST(run_on_fixtures_12_hamming ./run_and_diff.sh 5 hamming fixtures/12_seq.fasta fixtures/12_golden.sam fixtures/12_result.wit fixtures/12_result.json)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_CXX_FLAGS_=\"${CMAKE_CXX_FLAGS}\"")

# ------------------------------------------------------------
# Program: verify_sam
# ------------------------------------------------------------
ADD_EXECUTABLE(verify_sam verify_sam.cpp)

# ------------------------------------------------------------
# Program: do_search
# ------------------------------------------------------------
ADD_EXECUTABLE(do_search do_search.cpp)

ADD_EXECUTABLE(reweight_wit reweight_wit.cpp)

# ------------------------------------------------------------
# Program: compare_sam_wit
# ------------------------------------------------------------

# Add compare_sam_wit executable.
ADD_EXECUTABLE(compare_sam_wit compare_sam_wit.cpp)
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (compare_sam_wit rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

ADD_EXECUTABLE(sam_info sam_info.cpp)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=1")
ADD_EXECUTABLE(test_dp_quality test_dp_quality.cpp)
ADD_TEST(run_test_dp_quality test_dp_quality)

ADD_EXECUTABLE(fastq_info fastq_info.cpp)

ADD_EXECUTABLE(sam_quality_scale sam_quality_scale.cpp)

# ------------------------------------------------------------
# Program: test_wit_builder
# ------------------------------------------------------------
# XXX Uncomment the following line in production.
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_TESTING=1")
ADD_EXECUTABLE(test_curve_smoothing test_curve_smoothing.cpp)
ADD_EXECUTABLE(test_find_myers_ukkonen_reads test_find_myers_ukkonen_reads.cpp)

ADD_TEST(run_test_curve_smoothing test_curve_smoothing)
ADD_TEST(run_test_find_myers_ukkonen_reads test_find_myers_ukkonen_reads)

# EOF
