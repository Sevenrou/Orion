
find_package(Boost 1.37 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(Orion
  arbrecube.cpp
  getmeminfo.cpp
  main.cpp
  readfile.cpp
  utils.cpp
  os/mac_clock_gettime.cpp
)

if(NOT APPLE)
  target_link_libraries(Orion rt)
endif()
