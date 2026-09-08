cmake_minimum_required(VERSION 3.20)
project(autoremesher_reference_tests LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(UPSTREAM "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
find_package(TBB REQUIRED)

# Correct upstream's shallow redirect paths in the test build's include overlay.
file(GLOB redirects "${UPSTREAM}/include/AutoRemesher/*")
foreach(header IN LISTS redirects)
    file(READ "${header}" contents)
    string(REGEX MATCH "src/AutoRemesher/([^\"]+)" match "${contents}")
    get_filename_component(name "${header}" NAME)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/overlay/AutoRemesher/${name}"
        "#include \"${UPSTREAM}/src/AutoRemesher/${CMAKE_MATCH_1}\"\n")
endforeach()
file(GLOB engine "${UPSTREAM}/src/AutoRemesher/*.cpp")
file(GLOB isotropic "${UPSTREAM}/thirdparty/isotropicremesher/*.cpp")
add_executable(reference_surface_test reference_surface_test.cpp ${engine} ${isotropic}
    "${UPSTREAM}/thirdparty/meshoptimizer/src/simplifier.cpp"
    "${UPSTREAM}/thirdparty/meshoptimizer/src/indexgenerator.cpp")
target_include_directories(reference_surface_test PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/overlay" "${UPSTREAM}/thirdparty/eigen"
    "${UPSTREAM}/thirdparty/isotropicremesher" "${UPSTREAM}/thirdparty/meshoptimizer/src")
target_link_libraries(reference_surface_test PRIVATE TBB::tbb)
if(APPLE)
    target_compile_definitions(reference_surface_test PRIVATE AUTO_REMESHER_USE_ACCELERATE=1)
    target_compile_options(reference_surface_test PRIVATE -O3 -flto -funroll-loops)
    target_link_options(reference_surface_test PRIVATE -flto)
    target_link_libraries(reference_surface_test PRIVATE TBB::tbbmalloc TBB::tbbmalloc_proxy "-framework Accelerate")
endif()
enable_testing()
add_test(NAME preparation_reference COMMAND reference_surface_test)
