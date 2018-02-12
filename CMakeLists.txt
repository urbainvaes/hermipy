cmake_minimum_required(VERSION 3.10)

project(Test)
set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 0)
set (CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fpermissive")

find_package(PythonLibs 3 REQUIRED)
find_package(Boost COMPONENTS python3 REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})

find_library (Hermite NAMES hermite.so PATHS libhermite)
include_directories("libhermite/src")

file(GLOB_RECURSE Sources "src/*.cpp")
include_directories("src/")
add_executable(test.exe ${Sources})

set_target_properties(test.exe PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    )
target_link_libraries(test.exe ${Hermite} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${CMAKE_DL_LIBS})
