cmake_minimum_required(VERSION 3.0)

project(mckean-vlasov)
set (Tutorial_VERSION_MAJOR 0)
set (Tutorial_VERSION_MINOR 0)

set(CMAKE_BUILD_TYPE "DEBUG")
# set(CMAKE_BUILD_TYPE "RELEASE")

find_package(PythonLibs 3 REQUIRED)
find_package(Boost COMPONENTS python3 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS}) # ${PYTHON_INCLUDE_DIRS}
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_library(yay SHARED integrator)
target_link_libraries(yay ${Boost_LIBRARIES})
