
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(Orion
  arbrecube.cpp
  getmeminfo.cpp
  main.cpp
  readfile.cpp
  utils.cpp
)

if(APPLE)
  add_executable(Orion os/mac_clock_gettime.cpp)
endif(APPLE)

target_link_libraries(Orion rt)

