
# A home-made analogue of the 'FetchContent()' routine
# Yes, I have a few reason to avoid using CMake's ExternalProject or git's submodules
macro(fetch_3rd_party
      repo_url commit_hash target_name
      include_subpath release_lib_subpath debug_lib_subpath)
  set(_3rd_party_dir         "${3RDPARTY_DIR}/${target_name}")
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
                                -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
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
                                -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
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
                    WORKING_DIRECTORY ${_3rd_party_dir}
                    RESULTS_VARIABLE _3rd_party_err_code
                    OUTPUT_QUIET)
    if(${_3rd_party_err_code})
      message(FATAL_ERROR "Failed to compile 3rd-party library: ${target_name} (Release)")
    endif()

    if(MSVC)
      execute_process(COMMAND ${CMAKE_COMMAND}
                              --build "${_3rd_party_bin}"
                              --config Debug
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

if(NOT EXISTS ${3RDPARTY_DIR})
  file(MAKE_DIRECTORY "${3RDPARTY_DIR}")
endif()


if(MSVC)
  fetch_3rd_party("https://github.com/google/googletest.git"
                  "2f3e2e39cc4c399b66711e6b720bf22373e841b5"
                  googletest
                  "googletest/include"
                  "lib/Release/gtest.lib"
                  "lib/Debug/gtestd.lib")
else()
  fetch_3rd_party("https://github.com/google/googletest.git"
                    "2f3e2e39cc4c399b66711e6b720bf22373e841b5"
                    googletest
                    "googletest/include"
                    "lib/libgtest.a"
                    "lib/libgtest.a")
endif()
set(GTEST_ROOT  "${3RDPARTY_DIR}/${googletest}")
set(GTest_FOUND true)
include(GoogleTest)
