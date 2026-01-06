include(FetchContent)

FetchContent_Declare(
  dmclock
  URL https://github.com/ceph/ceph/archive/6e7fabd003ff21be9c077e00018505401053837f.zip
  SOURCE_SUBDIR src/dmclock
)

FetchContent_MakeAvailable(dmclock)
