# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.20)
project(AutoRemesherHeadless LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(AUTOREMESHER_BUILD_CLI "Build the GPL headless OBJ driver" ON)
option(AUTOREMESHER_BUILD_TESTS "Build geometry and driver regression tests" ON)
find_package(TBB REQUIRED)

# Preserve the existing test build's correction of shallow redirect paths.
file(GLOB redirects CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/AutoRemesher/*")
foreach(header IN LISTS redirects)
    file(READ "${header}" contents)
    string(REGEX MATCH "src/AutoRemesher/([^\"]+)" match "${contents}")
    if(NOT match)
        message(FATAL_ERROR "Invalid AutoRemesher redirect header: ${header}")
    endif()
    get_filename_component(name "${header}" NAME)
    # Preserve timestamps when the redirect is unchanged on reconfiguration.
    file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/overlay/AutoRemesher/${name}"
        CONTENT "#include \"${CMAKE_CURRENT_SOURCE_DIR}/src/AutoRemesher/${CMAKE_MATCH_1}\"\n"
        @ONLY)
endforeach()

file(GLOB engine CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/AutoRemesher/*.cpp")
file(GLOB isotropic CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/isotropicremesher/*.cpp")
add_library(autoremesher_core STATIC ${engine} ${isotropic}
    thirdparty/meshoptimizer/src/simplifier.cpp
    thirdparty/meshoptimizer/src/indexgenerator.cpp)
add_library(AutoRemesher::core ALIAS autoremesher_core)
set_target_properties(autoremesher_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(autoremesher_core PUBLIC cxx_std_14)
target_include_directories(autoremesher_core PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}/overlay"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/isotropicremesher"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/meshoptimizer/src")
target_link_libraries(autoremesher_core PUBLIC TBB::tbb)
if(APPLE)
    # Retain the established Mac headless solver and compiler options.
    target_compile_definitions(autoremesher_core PUBLIC AUTO_REMESHER_USE_ACCELERATE=1)
    target_compile_options(autoremesher_core PRIVATE -O3 -flto -funroll-loops)
    target_link_options(autoremesher_core INTERFACE -flto)
    target_link_libraries(autoremesher_core PUBLIC
        TBB::tbbmalloc TBB::tbbmalloc_proxy "-framework Accelerate")
endif()

if(AUTOREMESHER_BUILD_CLI)
    set(AUTOREMESHER_SOURCE_REVISION "" CACHE STRING "Source revision for builds without Git metadata")
    if(NOT AUTOREMESHER_SOURCE_REVISION)
        find_package(Git QUIET)
        if(Git_FOUND)
            execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}" rev-parse HEAD
                OUTPUT_VARIABLE AUTOREMESHER_SOURCE_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE revision_status ERROR_QUIET)
            execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}" status --porcelain --untracked-files=normal -- .
                OUTPUT_VARIABLE dirty OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE dirty_status ERROR_QUIET)
            if(NOT revision_status EQUAL 0 OR NOT dirty_status EQUAL 0)
                set(AUTOREMESHER_SOURCE_REVISION "unknown")
            elseif(dirty)
                string(APPEND AUTOREMESHER_SOURCE_REVISION "+dirty")
            endif()
        else()
            set(AUTOREMESHER_SOURCE_REVISION "unknown")
        endif()
    endif()
    if(NOT AUTOREMESHER_SOURCE_REVISION MATCHES "^([0-9a-f]+(\\+dirty)?|unknown)$")
        message(FATAL_ERROR "AUTOREMESHER_SOURCE_REVISION must be a hexadecimal commit, optionally +dirty, or unknown")
    endif()
    add_executable(AutoRemesherEngine cli/main.cpp)
    target_include_directories(AutoRemesherEngine PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_definitions(AutoRemesherEngine PRIVATE
        ENGINE_REVISION="${AUTOREMESHER_SOURCE_REVISION}"
        ENGINE_PATCH="local-singularity-1")
    target_link_libraries(AutoRemesherEngine PRIVATE AutoRemesher::core)
    if(MINGW)
        # Select the wide CRT entry point for wmain; only the CLI needs this.
        target_link_options(AutoRemesherEngine PRIVATE -municode)
    endif()
endif()

if(AUTOREMESHER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
