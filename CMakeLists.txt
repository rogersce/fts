CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(FTS)

option(ENABLE_STATIC "Build static (.a) library" ON)

set(CMAKE_VERBOSE_MAKEFILE off)
include_directories($ENV{HOME}/include src)
link_directories($ENV{HOME}/lib)
file(GLOB SRCS src/*.cpp)
file(GLOB HEADERS src/*.h)
file(GLOB SCRIPTS scripts/*.py)

add_library(fts SHARED ${SRCS})
target_link_libraries(fts cnpy)
target_link_libraries(fts fftw3)
find_library(silo NAMES "silo" "siloh5" PATHS ENV LD_LIBRARY_PATH)
target_link_libraries(fts ${silo})
install(TARGETS "fts" LIBRARY DESTINATION lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ENABLE_STATIC)
    add_library(fts-static STATIC ${SRCS})
    set_target_properties(fts-static PROPERTIES OUTPUT_NAME "fts")
    install(TARGETS "fts-static" ARCHIVE DESTINATION lib)
endif(ENABLE_STATIC)

install(FILES ${HEADERS} DESTINATION include/fts)
install(FILES ${SCRIPTS} DESTINATION scripts/fts PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

#set up all executables in models directory
message(STATUS "Configuring executables...")
file(GLOB EXES RELATIVE ${FTS_SOURCE_DIR}/models ${FTS_SOURCE_DIR}/models/*.cpp)
foreach(EXE ${EXES})
   STRING(REGEX REPLACE "[.]cpp" "" SHORTNAME ${EXE})
   message(STATUS "Configuring " ${SHORTNAME})
      
   add_executable(${SHORTNAME} models/${EXE})
   target_link_libraries(${SHORTNAME} fts)
endforeach(EXE ${EXES})
