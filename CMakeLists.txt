cmake_minimum_required(VERSION 3.12)


project(TISCC)


####################################################
# General configuration

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors -Wall -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -g3 -DNDEBUG -D_LIBCPP_DEBUG=1 ")
set(CMAKE_CXX_FLAGS_RELEASE_INIT  "-O3")
set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_BUILD_TYPE Debug)

if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    message("Setting stdlib=libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
endif ()

####################################################
# Dependencies

include_directories(external/include)

###################################################
# Library setup

add_library(
        TISCClib
        src/pipeline.cpp
        src/gridmanager.cpp
        src/logicalqubit.cpp
        src/plaquette.cpp
        src/hardwaremodel.cpp
        src/instruction.cpp
       )

set(TISCCLIB_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")

target_include_directories(
    TISCClib
    PUBLIC
    ${TISCCLIB_INCLUDE_DIRS}
)

set_property(TARGET TISCClib PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(
        TISCC
        src/main.cpp)

target_link_libraries(
        TISCC PUBLIC TISCClib
)

