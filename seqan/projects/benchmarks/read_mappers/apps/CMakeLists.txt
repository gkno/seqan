CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


add_executable (wit_builder wit_builder.cpp)
add_dependencies (wit_builder ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (wit_builder rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (verify_sam verify_sam.cpp)
add_dependencies (verify_sam ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (verify_sam rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (do_search do_search.cpp)
add_dependencies (do_search ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (do_search rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (reweight_wit reweight_wit.cpp)
add_dependencies (reweight_wit ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (reweight_wit rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (compare_sam_wit compare_sam_wit.cpp)
add_dependencies (compare_sam_wit ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (compare_sam_wit rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (sam_info sam_info.cpp)
add_dependencies (sam_info ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (sam_info rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (fastq_info fastq_info.cpp)
add_dependencies (fastq_info ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (fastq_info rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

add_executable (sam_quality_scale sam_quality_scale.cpp)
add_dependencies (sam_quality_scale ${SEQAN_TARGET})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  target_link_libraries (sam_quality_scale rt)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

