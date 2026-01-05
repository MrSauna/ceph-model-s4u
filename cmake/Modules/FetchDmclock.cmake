include(FetchContent)

FetchContent_Declare(
  dmclock
  URL https://github.com/ceph/ceph/archive/refs/tags/v19.2.3.zip
  SOURCE_SUBDIR src/dmclock
)

FetchContent_MakeAvailable(dmclock)
