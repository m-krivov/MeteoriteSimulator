
# A home-made analogue of the 'FetchContent()' routine
# Yes, I have a few reason to avoid using CMake's ExternalProject or git's submodules
macro(fetch_3rd_party
      repo_url commit_hash target_name cmake_options
      include_subpath release_lib_subpath debug_lib_subpath)
  set(_3rd_party_arch        "${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_PROCESSOR}")
  set(_3rd_party_dir         "${3RDPARTY_DIR}/${target_name}/${commit_hash}/${_3rd_party_arch}")
  set(_3rd_party_include     "${_3rd_party_dir}/${include_subpath}")
  set(_3rd_party_bin         "${_3rd_party_dir}/_Bin")
  set(_3rd_party_release_lib "${_3rd_party_bin}/${release_lib_subpath}")
  set(_3rd_party_debug_lib   "${_3rd_party_bin}/${debug_lib_subpath}")

  if(NOT EXISTS ${_3rd_party_include} OR
     NOT EXISTS ${_3rd_party_release_lib} OR
     NOT EXISTS ${_3rd_party_debug_lib})
    message(STATUS "Building '${target_name}'. Please, be patient ...")
    file(REMOVE_RECURSE ${_3rd_party_dir})
    file(MAKE_DIRECTORY ${_3rd_party_dir})

    # Clone repo
    execute_process(COMMAND git clone ${repo_url} ${_3rd_party_dir} --quiet
                    RESULTS_VARIABLE _3rd_party_err_code)
    if(${_3rd_party_err_code})
      message(FATAL_ERROR "Failed to clone 3rd-party library: ${target_name}")
    endif()

    # Bind it to a revision that works fine
    execute_process(COMMAND git checkout ${commit_hash} --quiet
                    WORKING_DIRECTORY ${_3rd_party_dir}
                    RESULTS_VARIABLE _3rd_party_err_code)
    if(${_3rd_party_err_code})
      message(FATAL_ERROR "Failed to checkout 3rd-party library: ${target_name}")
    endif()

    # Compile the project in Debug and Release modes, extract libraries
    if(MSVC)
      if(NOT ${CMAKE_GENERATOR_PLATFORM} STREQUAL "")
        execute_process(COMMAND ${CMAKE_COMMAND}
                                -B "${_3rd_party_bin}"
                                -G ${CMAKE_GENERATOR}
                                -A ${CMAKE_GENERATOR_PLATFORM}
                                ${cmake_options}
                                -DCMAKE_BUILD_TYPE=Debug;Release
                                -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
                                -DMSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
                                --no-warn-unused-cli
                        WORKING_DIRECTORY ${_3rd_party_dir}
                        RESULTS_VARIABLE _3rd_party_err_code
                        OUTPUT_QUIET)
      else()
        execute_process(COMMAND ${CMAKE_COMMAND}
                                -B "${_3rd_party_bin}"
                                -G ${CMAKE_GENERATOR}
                                ${cmake_options}
                                -DCMAKE_BUILD_TYPE=Debug;Release
                                -DMSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>DLL"
                                --no-warn-unused-cli
                        WORKING_DIRECTORY ${_3rd_party_dir}
                        RESULTS_VARIABLE _3rd_party_err_code
                        OUTPUT_QUIET)
      endif()
    else()
      execute_process(COMMAND ${CMAKE_COMMAND}
                                -B "${_3rd_party_bin}"
                                -G ${CMAKE_GENERATOR}
                                ${cmake_options}
                                --no-warn-unused-cli
                        WORKING_DIRECTORY ${_3rd_party_dir}
                        RESULTS_VARIABLE _3rd_party_err_code
                        OUTPUT_QUIET)
    endif()
    if(${_3rd_party_err_code})
      message(FATAL_ERROR "Failed to configure 3rd-party library: ${target_name}")
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND}
                            --build "${_3rd_party_bin}"
                            --config Release
                            --parallel 6
                    WORKING_DIRECTORY ${_3rd_party_dir}
                    RESULTS_VARIABLE _3rd_party_err_code
                    OUTPUT_QUIET
                    )
    if(${_3rd_party_err_code})
      message(FATAL_ERROR "Failed to compile 3rd-party library: ${target_name} (Release)")
    endif()

    if(MSVC)
      execute_process(COMMAND ${CMAKE_COMMAND}
                              --build "${_3rd_party_bin}"
                              --config Debug
                              --parallel 6
                      WORKING_DIRECTORY ${_3rd_party_dir}
                      RESULTS_VARIABLE _3rd_party_err_code
                      OUTPUT_QUIET)
      if(${_3rd_party_err_code})
        message(FATAL_ERROR "Failed to compile 3rd-party library: ${target_name} (Debug)")
      endif()
    endif()

    if(NOT EXISTS ${_3rd_party_include} OR
       NOT EXISTS ${_3rd_party_release_lib} OR
       NOT EXISTS ${_3rd_party_debug_lib})
      message(FATAL_ERROR "Building rule for '${target_name}' is broken")
    endif()
  endif()

  add_library(${target_name} STATIC IMPORTED)
  target_include_directories(${target_name} SYSTEM INTERFACE ${_3rd_party_include})
  if(MSVC)
  set_target_properties(${target_name} PROPERTIES
                        IMPORTED_LOCATION_RELEASE ${_3rd_party_release_lib}
                        IMPORTED_LOCATION_DEBUG ${_3rd_party_debug_lib})
  else()
    set_target_properties(${target_name} PROPERTIES
                          IMPORTED_LOCATION ${_3rd_party_release_lib})
  endif()
  
endmacro()

# Downloads some release of a simple header-only library and creates a target for it
macro(download_3rd_party name src_url header_subdir)
  set(_TARGET_DIR "${3RDPARTY_DIR}/${name}")
  get_filename_component("${src_url}" _TARGET_EXT NAME_WLE)
  if(NOT EXISTS "${_TARGET_DIR}")
    message(STATUS "Downloading a header-only library (${name}) ...")
    file(DOWNLOAD "${src_url}" "${_TARGET_DIR}/src${_TARGET_EXT}")
    file(ARCHIVE_EXTRACT
         INPUT "${_TARGET_DIR}/src${_TARGET_EXT}"
         DESTINATION "${_TARGET_DIR}")
    file(REMOVE "${_TARGET_DIR}/src${_TARGET_EXT}")
  endif()

  add_library(${name} INTERFACE)
  if("${header_subdir}" STREQUAL "")
    target_include_directories(${name} SYSTEM INTERFACE "${_TARGET_DIR}")
  else()
    target_include_directories(${name} SYSTEM INTERFACE "${_TARGET_DIR}/${header_subdir}")
  endif()
endmacro()

if(NOT EXISTS ${3RDPARTY_DIR})
  file(MAKE_DIRECTORY "${3RDPARTY_DIR}")
endif()


### GoogleTest ###

if(MSVC)
  fetch_3rd_party("https://github.com/google/googletest.git"
                  "2f3e2e39cc4c399b66711e6b720bf22373e841b5"
                  googletest
                  -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
                  "googletest/include"
                  "lib/Release/gtest.lib"
                  "lib/Debug/gtestd.lib")
else()
  fetch_3rd_party("https://github.com/google/googletest.git"
                  "2f3e2e39cc4c399b66711e6b720bf22373e841b5"
                  googletest
                  -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
                  "googletest/include"
                  "lib/libgtest.a"
                  "lib/libgtest.a")
endif()
set(GTEST_ROOT  "${3RDPARTY_DIR}/${googletest}")
set(GTest_FOUND true)
include(GoogleTest)

### Indicators ###

set(INDICATORS_VER "2.3")
download_3rd_party(indicators
                   "https://github.com/p-ranav/indicators/archive/refs/tags/v${INDICATORS_VER}.zip"
                   "indicators-${INDICATORS_VER}/include")

### Matplot++ ###

if(${METEORITES_GNUPLOT})
  if(MSVC)
    fetch_3rd_party("https://github.com/alandefreitas/matplotplusplus.git"
                    "b2fed97fca6e7a5380efe14ec76e227c2087b77e"
                    matplotpp
                    -DMATPLOTPP_BUILD_EXAMPLES=OFF
                    "source"
                    "source/matplot/Release/matplot.lib"
                    "source/matplot/Debug/matplot.lib")
  else()
    fetch_3rd_party("https://github.com/alandefreitas/matplotplusplus.git"
                    "b2fed97fca6e7a5380efe14ec76e227c2087b77e"
                    matplotpp
                    -DMATPLOTPP_BUILD_EXAMPLES=OFF
                    "source"
                    "source/matplot/libmatplot.a"
                    "source/matplot/libmatplot.a")
  endif()
  # TODO: replace by introducing global variables for home directories
  file(COPY_FILE
       "${3RDPARTY_DIR}/matplotpp/b2fed97fca6e7a5380efe14ec76e227c2087b77e/${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_PROCESSOR}/_Bin/source/matplot/matplot/detail/exports.h"
       "${3RDPARTY_DIR}/matplotpp/b2fed97fca6e7a5380efe14ec76e227c2087b77e/${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_PROCESSOR}/source/matplot/detail/exports.h")
endif()
