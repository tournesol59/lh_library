cmake_minimum_required (VERSION 2.6)
project (ident05_main)

# The version number.
set (ident05_main_VERSION_MAJOR 1)
set (ident05_main_VERSION_MINOR 0)

# should we use our own math functions
option(USE_MYMATH "Use tutorial provided math implementation" ON)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/TutorialConfig.h.in"
  "${PROJECT_BINARY_DIR}/TutorialConfig.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
set (HOME_DIR "/home/frede")
include_directories ("${PROJECT_BINARY_DIR}",
	"${HOME_DIR}/Ipopt-3.12.6/include/coin/")

# add the MathFunctions library?
if (USE_MYMATH)
  include_directories ("${PROJECT_SOURCE_DIR}/MathFunctions")
  add_subdirectory (MathFunctions)
  set (EXTRA_LIBS ${EXTRA_LIBS} MathFunctions)
endif ()

# add the executable
add_executable (ident05_main ident05_test.cpp)
target_link_libraries (ident05_main  ${EXTRA_LIBS})
