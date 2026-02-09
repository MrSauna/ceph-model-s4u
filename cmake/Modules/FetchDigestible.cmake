# Fetch digestible - A modern C++ header-only T-Digest implementation
# https://github.com/SpirentOrion/digestible
#
# Provides streaming percentile calculation with O(1) memory.
# Usage: #include "digestible/digestible.h"

include(FetchContent)

FetchContent_Declare(
  digestible
  GIT_REPOSITORY https://github.com/SpirentOrion/digestible.git
  GIT_TAG        master
)

FetchContent_MakeAvailable(digestible)

# Create interface library for header-only usage
# Headers are in the 'include' subdirectory
add_library(digestible_lib INTERFACE)
target_include_directories(digestible_lib INTERFACE ${digestible_SOURCE_DIR}/include)

