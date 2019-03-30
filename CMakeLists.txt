cmake_minimum_required(VERSION 3.0)
project(circle_fit)
add_compile_options(-std=c++17 -Wall -Wpedantic -Werror -fverbose-asm)
find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)
find_package(gflags REQUIRED)

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/src/generated/lm_derivatives.hpp
    COMMAND python3 src/python/generate_lm_derivatives.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    DEPENDS src/generated/lm_derivatives.hpp.template src/python/generate_lm_derivatives.py src/python/codegenerator.py 
    COMMENT "Generating code for lm_derivates.cpp."
)
add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/src/generated/Q3.hpp ${PROJECT_SOURCE_DIR}/src/generated/M_C.hpp
    COMMAND python3 src/python/generate_taubin_determinant.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    DEPENDS src/generated/Q3.hpp.template src/generated/M_C.hpp.template src/python/generate_taubin_determinant.py src/python/codegenerator.py 
    COMMENT "Generating code for Q3.hpp and M_C.hpp."
)

add_library(circle_fit STATIC
include/circle_fit/common.h
include/circle_fit/gwaf_taubin.h
src/gwaf_taubin.cpp
include/circle_fit/geo_lm.h
src/geo_lm.cpp
)

set_property(SOURCE src/geo_lm.cpp APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/src/generated/lm_derivatives.hpp)
set_property(SOURCE src/gwaf_taubin.cpp APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/src/generated/M_C.hpp ${PROJECT_SOURCE_DIR}/src/generated/Q3.hpp)

target_include_directories(circle_fit PUBLIC
include
${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(circle_fit
glog
)



set(WITH_TESTS "off" CACHE BOOL "Build tests (requires googletest and GTEST_DIR variable to be set)")
if (WITH_TESTS)
    # Tests
    enable_testing()
    # find_package(Boost REQUIRED COMPONENTS unit_test_framework)
    add_subdirectory(${GTEST_DIR} googletest EXCLUDE_FROM_ALL)

    # Common unit tests
    add_executable(test_circle_fit test/test_circle_fit.cpp)
    target_link_libraries(test_circle_fit
    circle_fit
    )
    target_include_directories(test_circle_fit PRIVATE ${BOOST_INCLUDE_DIRS})
    add_test(NAME unit_tests COMMAND test_circle_fit)

    # GWAF Taubin specific unit tests
    add_executable(test_gwaf_taubin
    test/test_gwaf_taubin.cpp
    )

    target_link_libraries(test_gwaf_taubin
    circle_fit
    )
    # target_include_directories(test_circle_fit PRIVATE ${BOOST_INCLUDE_DIRS})
    add_test(NAME gwaf_taubin_tests COMMAND test_gwaf_taubin)

    # Geometric LM specific unit tests
    add_executable(test_geo_lm
    test/test_geo_lm.cpp
    )

    target_link_libraries(test_geo_lm
    circle_fit
    gtest_main
    gflags
    )
    target_include_directories(test_circle_fit PRIVATE ${BOOST_INCLUDE_DIRS})
    add_test(NAME geo_lm_tests COMMAND test_geo_lm)
endif()

