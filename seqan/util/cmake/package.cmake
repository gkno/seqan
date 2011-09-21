INCLUDE(InstallRequiredSystemLibraries)

# NOTE that you have to run "make docs" before running cpack.  The reason
# is that we cannot add dependencies to the install target at the moment.
# See: http://public.kitware.com/Bug/view.php?id=8438

# ===========================================================================
# Archive Packages (.tar & .tar.bz2)
# ===========================================================================

SET(CPACK_GENERATOR "ZIP;TBZ2")
SET(CPACK_PACKAGE_NAME "seqan")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
SET(CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION "1.3.2")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "2")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

# Should be the last include.
INCLUDE(CPack)

# ===========================================================================
# Debian Packages
# ===========================================================================

SET(CPACK_GENERATOR "DEB")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn - The C++ library for sequence analysis.")

# Should be the last include.
INCLUDE(CPack)
