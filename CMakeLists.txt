cmake_minimum_required(VERSION 3.14)

project(hyperbolic)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/heads/master.zip
    )
FetchContent_Declare(
    glm
    URL https://github.com/g-truc/glm/archive/refs/heads/master.zip
    )

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FetchContent_MakeAvailable(glm)

enable_testing()

add_executable(
    tests
    test/test_glm_lorentz.cpp
    test/test_high_precision.cpp
    test/test_honeycombs.cpp
    test/test_hyperbolic.cpp
    test/test_liealgebra.cpp
    test/test_rigid_body.cpp
    )

target_include_directories(
    tests PRIVATE
    "include"
)
target_link_libraries(
    tests
    gtest_main
    glm::glm
    )

include(GoogleTest)
gtest_discover_tests(tests)