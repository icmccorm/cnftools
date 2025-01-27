cmake_minimum_required(VERSION 3.10)

project(cnftools VERSION 1.0)
include_directories("/opt/homebrew/Cellar/libarchive/3.5.2/include")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(LibArchive_INCLUDE_DIR "/usr/local/opt/libarchive/lib")
set(LIBARCHIVE_INCLUDE_DIR "/usr/local/opt/libarchive/lib")

set(CADICAL_CONFIGURE "./configure")
set(CADICAL_CONFIGURE_OPTION "-fPIC")
set(CADICAL_BUILD make -j 8)
if (SKIP_CADICAL) 
    message("Warning: SKIP_CADICAL suppresses rebuild of CaDiCaL")
    set(CADICAL_CONFIGURE "")
    set(CADICAL_CONFIGURE_OPTION "")
    set(CADICAL_BUILD "")
endif()

include(ExternalProject)
ExternalProject_Add(Cadical
    GIT_REPOSITORY https://github.com/arminbiere/cadical
    GIT_TAG "origin/master"
    PREFIX "${CMAKE_BINARY_DIR}/cadical"
    CONFIGURE_COMMAND "${CADICAL_CONFIGURE}" "${CADICAL_CONFIGURE_OPTION}"
    BUILD_COMMAND "${CADICAL_BUILD}"
    INSTALL_COMMAND ""
    #BUILD_BYPRODUCTS "${CMAKE_BINARY_DIR}/cadical/src/Cadical/build/libcadical.a"
    BUILD_IN_SOURCE 1 
)

add_library(solver STATIC IMPORTED)
add_dependencies(solver Cadical)
set_target_properties(solver PROPERTIES IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/cadical/src/Cadical/build/libcadical.a")

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIRS})
set(LIBS ${LIBS} md5 ${LibArchive_LIBRARIES})

include_directories(cnftools PUBLIC "${PROJECT_SOURCE_DIR}/src")

add_subdirectory("lib/md5")
add_subdirectory("src")

add_executable(cnftools src/Main.cc)
add_dependencies(cnftools solver)
target_link_libraries(cnftools PUBLIC ${LIBS} solver $<TARGET_OBJECTS:gates> $<TARGET_OBJECTS:util> $<TARGET_OBJECTS:transform> $<TARGET_OBJECTS:features>)

target_include_directories(cnftools PUBLIC "${PROJECT_SOURCE_DIR}")
