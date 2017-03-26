cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} "-std=c++0x")

find_package(SymEngine 0.1.0 CONFIG) 

include_directories(${SYMENGINE_INCLUDE_DIRS})
add_executable(main main.cpp)
target_link_libraries(main ${SYMENGINE_LIBRARIES})
