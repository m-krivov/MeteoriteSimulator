cmake_minimum_required(VERSION 3.22)
project(Meteorites LANGUAGES CXX)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)
set(CMAKE_SUPPRESS_REGENERATION ON)
unset(CMAKE_INSTALL_PREFIX CACHE)

set(ROOT_DIR      "${CMAKE_SOURCE_DIR}")
set(BIN_DIR       "${ROOT_DIR}/Bin")
set(3RDPARTY_DIR  "${ROOT_DIR}/3rdParty")

if(NOT EXISTS ${BIN_DIR})
  file(MAKE_DIRECTORY "${BIN_DIR}")
endif()

macro(find_cpp_sources dir headers cpp)
  file(GLOB_RECURSE ${headers}   ${dir}/*.h)
  file(GLOB_RECURSE ${cpp}       ${dir}/*.cpp)

  if(MSVC)
    foreach(_source IN ITEMS ${${headers}})
      get_filename_component(_source_path "${_source}" PATH)
      file(RELATIVE_PATH _source_path_rel "${dir}" "${_source_path}")
      string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
      source_group("Header Files\\${_group_path}" FILES "${_source}")
    endforeach()

    foreach(_source IN ITEMS ${${cpp}})
      get_filename_component(_source_path "${_source}" PATH)
      file(RELATIVE_PATH _source_path_rel "${dir}" "${_source_path}")
      string(REPLACE "/" "\\" _group_path "${_source_path_rel}")
      source_group("Source Files\\${_group_path}" FILES "${_source}")
    endforeach()
  endif()
endmacro()

macro(configure_cpp_proj target proj_name)
  target_include_directories(${target} PUBLIC ${ROOT_DIR})
  set_target_properties(${target} PROPERTIES
                                  PREFIX ""
                                  OUTPUT_NAME ${proj_name})
  target_compile_features(${target} PRIVATE cxx_std_17)
  set(args ${ARGN})
  list(LENGTH args arg_num)
  if(${arg_num} GREATER 0)
    math(EXPR arg_num ${arg_num}-1)
    foreach(i RANGE 0 ${arg_num})
      list(GET args ${i} arg)
      get_target_property(arg_type ${arg} TYPE)
      if(arg_type STREQUAL "STATIC_LIBRARY")
        target_link_libraries(${target} ${arg})
      endif()
    endforeach(i)
  endif()

  get_target_property(target_type ${target} TYPE)
  if(target_type STREQUAL "EXECUTABLE" OR
     target_type STREQUAL "SHARED_LIBRARY" OR
     target_type STREQUAL "MODULE_LIBRARY")
    add_custom_command(TARGET ${target} POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${target}> "${BIN_DIR}")
  endif()
  if(MSVC AND target_type STREQUAL "EXECUTABLE")
    set_property(TARGET ${target} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${BIN_DIR}")
  endif()
endmacro()


find_cpp_sources(${ROOT_DIR}/Meteorites.Core METEORITES_CORE_H METEORITES_CORE_CPP)
add_library(Meteorites.Core STATIC ${METEORITES_CORE_H} ${METEORITES_CORE_CPP})
configure_cpp_proj(Meteorites.Core "meteorites.core")

find_cpp_sources(${ROOT_DIR}/Meteorites.CpuSolvers METEORITES_CPU_H METEORITES_CPU_CPP)
add_library(Meteorites.CpuSolvers STATIC ${METEORITES_CPU_H} ${METEORITES_CPU_CPP})
configure_cpp_proj(Meteorites.CpuSolvers "meteorites.cpu")

find_cpp_sources(${ROOT_DIR}/Meteorites.KnowledgeBase METEORITES_KNOWLEDGE_H METEORITES_KNOWLEDGE_CPP)
add_library(Meteorites.KnowledgeBase STATIC ${METEORITES_KNOWLEDGE_H} ${METEORITES_KNOWLEDGE_CPP})
configure_cpp_proj(Meteorites.KnowledgeBase "meteorites.knowledge")

find_cpp_sources(${ROOT_DIR}/Meteorites.App METEORITES_APP_H METEORITES_APP_CPP)
add_executable(Meteorites.App ${METEORITES_APP_H} ${METEORITES_APP_CPP})
configure_cpp_proj(Meteorites.App "meteorites"
                   Meteorites.KnowledgeBase
                   Meteorites.CpuSolvers
                   Meteorites.Core)

find_cpp_sources(${ROOT_DIR}/Meteorites.Benchmark METEORITES_BENCHMARK_H METEORITES_BENCHMARK_CPP)
add_executable(Meteorites.Benchmark ${METEORITES_BENCHMARK_H} ${METEORITES_BENCHMARK_CPP})
configure_cpp_proj(Meteorites.Benchmark "meteorites.benchmark"
                   Meteorites.KnowledgeBase
                   Meteorites.CpuSolvers
                   Meteorites.Core)


include(3rdParty.cmake)
find_cpp_sources(${ROOT_DIR}/Meteorites.Tests METEORITES_TESTS_H METEORITES_TESTS_CPP)
add_executable(Meteorites.Tests ${METEORITES_TESTS_H} ${METEORITES_TESTS_CPP})
configure_cpp_proj(Meteorites.Tests "meteorites.tests"
                   Meteorites.KnowledgeBase
                   Meteorites.CpuSolvers
                   Meteorites.Core
                   googletest)
