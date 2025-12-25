
#########################################
###         DepCache internals        ###
### No need to invoke them explicitly ###
#########################################


# Lazy initialization for DepCache
macro(_depcache_init)
  if(NOT DEFINED DEPCACHE_DIR)
    set(DEPCACHE_DIR "${CMAKE_SOURCE_DIR}/DepCache")
  endif()
  if(NOT EXISTS ${DEPCACHE_DIR})
    file(MAKE_DIRECTORY "${DEPCACHE_DIR}")
  endif()

  if(NOT DEFINED DEPCACHE_PARALLEL)
    if(DEFINED CMAKE_BUILD_PARALLEL_LEVEL AND
       NOT "${CMAKE_BUILD_PARALLEL_LEVEL}" STREQUAL "")
      set(DEPCACHE_PARALLEL ${CMAKE_BUILD_PARALLEL_LEVEL})
    else()
      cmake_host_system_information(RESULT logical_cores QUERY NUMBER_OF_LOGICAL_CORES)
      set(DEPCACHE_PARALLEL ${logical_cores})
    endif()
  endif()

  set(DEPCACHE_OUTPUT_FILE "${DEPCACHE_DIR}/output.log")
  set(DEPCACHE_ERROR_FILE  "${DEPCACHE_DIR}/error.log")
  set(DEPCACHE_ARCH        "${CMAKE_CXX_COMPILER_ID}_${CMAKE_SYSTEM_PROCESSOR}")

  file(WRITE ${DEPCACHE_OUTPUT_FILE} "")
  file(WRITE ${DEPCACHE_ERROR_FILE}  "")
endmacro()


# A wrapper around 'execute_process()' that redirects streams and checks the error code
# Use 'ERROR_MESSAGE <text>' to provide a custom error message
function(_depcache_execute_process)
  set(one_value_args ERROR_MESSAGE OUTPUT_FILE ERROR_FILE RESULT_VARIABLE)
  cmake_parse_arguments(arg "" "${one_value_args}" "" ${ARGV})
  if(DEFINED arg_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "DepCache: failed to understand arguments '${ARGV}'")
  endif()

  # Fulfill potentially omitted arguments with proper defaults
  if(NOT DEFINED arg_ERROR_MESSAGE)
    set(arg_ERROR_MESSAGE "DepCache: a command was not executed successfully")
  endif()
  if(NOT DEFINED arg_OUTPUT_FILE)
    set(arg_OUTPUT_FILE ${DEPCACHE_OUTPUT_FILE})
  endif()
  if(NOT DEFINED arg_ERROR_FILE)
    set(arg_ERROR_FILE ${DEPCACHE_ERROR_FILE})
  endif()

  execute_process(${arg_UNPARSED_ARGUMENTS}
                  OUTPUT_FILE "${arg_OUTPUT_FILE}"
                  ERROR_FILE "${arg_ERROR_FILE}"
                  RESULT_VARIABLE exec_result)
  if(${exec_result})
    message(FATAL_ERROR "${arg_ERROR_MESSAGE}")
  endif()
endfunction()


# Clones a repo and checks out the specific commit
# May overwrite local changes, can be called multiple times, works offline (if nothing to do)
# Supports optional 4th argument: was the state of repo updated?
function(_depcache_clone repo commit dir)
  if(NOT EXISTS "${dir}")
    message(STATUS "DepCache: clonning a library '${target}' ...")
    _depcache_execute_process(COMMAND git clone ${repo} ${dir})
    set(repo_was_changed TRUE)
  else()
    set(repo_was_changed FALSE)
  endif()

  # A hack to speed up the configure step
  # Also it allows us to work without access to Internet
  # Note: it does not work with branches (commits only)
  if(EXISTS "${dir}/.git/HEAD")
    file(STRINGS "${dir}/.git/HEAD" cur_rev LIMIT_COUNT 1)
  endif()
  if(NOT "${cur_rev}" STREQUAL "${commit}")
    _depcache_execute_process(COMMAND git reset --hard ${commit}
                              WORKING_DIRECTORY ${dir})

    _depcache_execute_process(COMMAND git clean -fdx
                              WORKING_DIRECTORY ${dir})

    _depcache_execute_process(COMMAND git fetch
                              WORKING_DIRECTORY ${dir})

    _depcache_execute_process(COMMAND git checkout ${commit}
                              WORKING_DIRECTORY ${dir})

    set(repo_was_changed true)
  endif()

  # Notify the caller about possible changes
  if(${ARGC} EQUAL 4)
    set(${ARGV3} ${repo_was_changed})
  endif()
endfunction()


# Clones a repo, compiles it and creates the requested target
# The following arguments are expected:
#   target
#   REPO git_url
#   COMMIT commit_hash
#   CMAKE_OPTIONS option1 option2 ...
#   REPO_DIR directory
#   BUILD_DIR directory
#   MSVC_LIB_SUBPATH_DEBUG subpath
#   MSVC_LIB_SUBPATH_RELEASE subpath
#   GNU_LIB_SUBPATH subpath
#   INCLUDE_SUBPATH subpath
function(_depcache_make_target target)
  # Extract arguments, ensure they are correct
  set(one_value_args REPO COMMIT
                     REPO_DIR BUILD_DIR
                     MSVC_LIB_SUBPATH_DEBUG
                     MSVC_LIB_SUBPATH_RELEASE
                     GNU_LIB_SUBPATH
                     INCLUDE_SUBPATH)
  set(multi_value_args CMAKE_OPTIONS)
  cmake_parse_arguments(arg "" "${one_value_args}" "${multi_value_args}" ${ARGN})
  if(DEFINED arg_UNPARSED_ARGUMENTS OR DEFINED arg_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "DepCache: failed to parse arguments for target '${target}'")
  endif()
  
  set(include_dir        "${arg_REPO_DIR}/${arg_INCLUDE_SUBPATH}")
  set(msvc_lib_debug     "${arg_BUILD_DIR}/${arg_MSVC_LIB_SUBPATH_DEBUG}")
  set(msvc_lib_release   "${arg_BUILD_DIR}/${arg_MSVC_LIB_SUBPATH_RELEASE}")
  set(gnu_lib            "${arg_BUILD_DIR}/${arg_GNU_LIB_SUBPATH}")

  # Prepare directory with sources, check its previous status
  _depcache_clone(${arg_REPO} ${arg_COMMIT} ${arg_REPO_DIR})
  set(rebuild ON)
  if(MSVC)
    if(EXISTS ${include_dir} AND
       EXISTS ${msvc_lib_debug} AND
       EXISTS ${msvc_lib_release})
      set(rebuild OFF)
    endif()
  else()
    if(EXISTS ${include_dir} AND
       EXISTS ${gnu_lib})
      set(rebuild OFF)
    endif()
  endif()
  
  # Ok, need to compile or recompile the target
  if(${rebuild})
    message(STATUS "DepCache: building '${target}'. Please, be patient ...")
    
    # Configure the project for Debug and Release modes
    set(configure_options -B "${arg_BUILD_DIR}"
                          -G ${CMAKE_GENERATOR}
                          -DCMAKE_BUILD_TYPE="Debug;Release"
                          -DDEPCACHE_DIR=${DEPCACHE_DIR}
                          -DDEPCACHE_PARALLEL=${DEPCACHE_PARALLEL}
                          ${arg_CMAKE_OPTIONS}
                          --no-warn-unused-cli)
    if(MSVC)
      list(APPEND configure_options -DCMAKE_POLICY_DEFAULT_CMP0091=NEW
                                    -DMSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
      if(NOT "${CMAKE_GENERATOR_PLATFORM}" STREQUAL "")
        list(APPEND configure_options -A ${CMAKE_GENERATOR_PLATFORM})
      endif()
    endif()

    _depcache_execute_process(COMMAND ${CMAKE_COMMAND} ${configure_options}
                              WORKING_DIRECTORY "${arg_REPO_DIR}")
    if(${exec_result})
      message(FATAL_ERROR "DepCache: failed to configure target '${target}'")
    endif()

    # Build the target (Release for gcc, Debug/Release for Visual C++)
    _depcache_execute_process(COMMAND ${CMAKE_COMMAND}
                                      --build "${arg_BUILD_DIR}"
                                      --config Release
                                      --parallel ${DEPCACHE_PARALLEL}
                              WORKING_DIRECTORY ${arg_REPO_DIR})
    
    if(MSVC)
      _depcache_execute_process(COMMAND ${CMAKE_COMMAND}
                                        --build "${arg_BUILD_DIR}"
                                        --config Debug
                                        --parallel ${DEPCACHE_PARALLEL}
                                WORKING_DIRECTORY ${arg_REPO_DIR})
    endif()
  endif()

  # Final checks, create the requested target
  if(NOT EXISTS ${include_dir} OR
     (MSVC AND NOT EXISTS ${msvc_lib_debug}) OR
     (MSVC AND NOT EXISTS ${msvc_lib_release}) OR
     (NOT MSVC AND NOT EXISTS ${gnu_lib}))
    message(FATAL_ERROR "DepCache: building rule for '${target}' is broken, more details in logs")
  endif()
  
  add_library(${target} STATIC IMPORTED)
  target_include_directories(${target} SYSTEM INTERFACE "${include_dir}")
  if(MSVC)
    set_target_properties(${target} PROPERTIES
                          IMPORTED_LOCATION_DEBUG "${msvc_lib_debug}"
                          IMPORTED_LOCATION_RELEASE "${msvc_lib_release}")
  else()
    set_target_properties(${target} PROPERTIES
                          IMPORTED_LOCATION "${gnu_lib}")
  endif()
endfunction()


# Finds and extracts commit hash from the arguments
function(_depcache_extract_commit hash)
  set(one_value_args COMMIT)
  cmake_parse_arguments(arg "" "${one_value_args}" "" ${ARGN})
  if(NOT DEFINED arg_COMMIT OR "${arg_COMMIT}" STREQUAL "")
    message(FATAL_ERROR "DepCache: failed to extract argument 'COMMIT' from the string '${ARGN}'")
  endif()
  set(${hash} ${arg_COMMIT} PARENT_SCOPE)
endfunction()


#################################################
###               DepCache API                ###
###  Use these macros in your CMakeLists.txt  ###
#################################################


macro(depcache_repo target)
  _depcache_init()
  _depcache_extract_commit(_DEPCACHE_COMMIT ${ARGN})
  set("${target}_SOURCE_DIR"  "${DEPCACHE_DIR}/${target}/${_DEPCACHE_COMMIT}/${DEPCACHE_ARCH}")
  set("${target}_BUILD_DIR"   "${DEPCACHE_DIR}/${target}/${_DEPCACHE_COMMIT}/${DEPCACHE_ARCH}_build")
  _depcache_make_target(${target}
                        REPO_DIR "${${target}_SOURCE_DIR}"
                        BUILD_DIR "${${target}_BUILD_DIR}"
                        ${ARGN})
endmacro()


# Downloads some release of a simple header-only library and creates a target for it
macro(depcache_header_only target src_url header_subdir)
  _depcache_init()
  set("${target}_SOURCE_DIR"  "${DEPCACHE_DIR}/${target}")
    if(NOT EXISTS "${${target}_SOURCE_DIR}")
    message(STATUS "DepCache: downloading a header-only library '${target}' ...")
    file(DOWNLOAD "${src_url}" "${${target}_SOURCE_DIR}/src${_DEPCACHE_EXT}")
    get_filename_component("${src_url}" _DEPCACHE_EXT NAME_WLE)
    file(ARCHIVE_EXTRACT
         INPUT       "${${target}_SOURCE_DIR}/src${_DEPCACHE_EXT}"
         DESTINATION "${${target}_SOURCE_DIR}")
    file(REMOVE "${${target}_SOURCE_DIR}/src${_DEPCACHE_EXT}")
  endif()

  add_library(${target} INTERFACE)
  if("${header_subdir}" STREQUAL "")
    target_include_directories(${target} SYSTEM INTERFACE "${${target}_SOURCE_DIR}")
  else()
    target_include_directories(${target} SYSTEM INTERFACE "${${target}_SOURCE_DIR}/${header_subdir}")
  endif()
endmacro()


###################################################
###    Several presets for popular libraries    ###
###  Feel free to edit them or to add new ones  ###
###################################################


### GoogleTest ###
macro(depcache_googletest target)
  depcache_repo(${target}
                REPO                     https://github.com/google/googletest.git
                COMMIT                   2f3e2e39cc4c399b66711e6b720bf22373e841b5
                INCLUDE_SUBPATH          "googletest/include"
                MSVC_LIB_SUBPATH_DEBUG   "lib/Debug/gtestd.lib"
                MSVC_LIB_SUBPATH_RELEASE "lib/Release/gtest.lib"
                GNU_LIB_SUBPATH          "lib/libgtest.a")
  set(GTEST_ROOT  "${target}_REPO_DIR")
  set(GTest_FOUND TRUE)
endmacro()


### Indicators ###
depcache_header_only(indicators
                    "https://github.com/p-ranav/indicators/archive/refs/tags/v2.3.zip"
                    "indicators-2.3/include")


### Matplot++ ###
macro(depcache_matplotpp target)
  depcache_repo(${target}
                REPO                     https://github.com/alandefreitas/matplotplusplus.git
                COMMIT                   b2fed97fca6e7a5380efe14ec76e227c2087b77e
                CMAKE_OPTIONS            -DMATPLOTPP_BUILD_EXAMPLES=OFF -DMATPLOTPP_BUILD_WITH_SANITIZERS=OFF
                INCLUDE_SUBPATH          "source"
                MSVC_LIB_SUBPATH_DEBUG   "source/matplot/Debug/matplot.lib"
                MSVC_LIB_SUBPATH_RELEASE "source/matplot/Release/matplot.lib"
                GNU_LIB_SUBPATH          "source/matplot/libmatplot.a")

  add_library("${target}_nodesoup" STATIC IMPORTED)
  if(MSVC)
    set_target_properties(matplotpp_nodesoup PROPERTIES
                          IMPORTED_LOCATION_DEBUG   "${${target}_BUILD_DIR}/source/3rd_party/Debug/nodesoup.lib"
                          IMPORTED_LOCATION_RELEASE "${${target}_BUILD_DIR}/source/3rd_party/Release/nodesoup.lib")
  else()
    set_target_properties(matplotpp_nodesoup PROPERTIES
                          IMPORTED_LOCATION         "${${target}_BUILD_DIR}/source/3rd_party/libnodesoup.a")
  endif()
  file(COPY_FILE
       "${${target}_BUILD_DIR}/source/matplot/matplot/detail/exports.h"
       "${${target}_SOURCE_DIR}/source/matplot/detail/exports.h")
  target_link_libraries(${target} PRIVATE INTERFACE "${target}_nodesoup")
endmacro()
