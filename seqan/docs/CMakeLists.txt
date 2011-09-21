find_package (PythonInterp)

set(SEQAN_DOCS_PATH "${CMAKE_SOURCE_DIR}/docs")

add_custom_target(docs
                  WORKING_DIRECTORY ${SEQAN_DOCS_PATH}
                  COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_DOCS_PATH}/main.py ${CMAKE_SOURCE_DIR}/core/include ${CMAKE_SOURCE_DIR}/extras/include -d ${SEQAN_DOCS_PATH}/concepts -d ${SEQAN_DOCS_PATH}/pages -o ${CMAKE_BINARY_DIR}/docs/html)

install(DIRECTORY ${CMAKE_BINARY_DIR}/docs/html
        DESTINATION share/seqan/docs)
