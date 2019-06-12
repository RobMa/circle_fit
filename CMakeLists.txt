cmake_minimum_required(VERSION 3.3)
project(circle_fit)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Download and unpack external dependencies
set(EXTERNAL_PROJECTS_DIR "${PROJECT_SOURCE_DIR}/extern")
configure_file(CMakeLists.txt.in external-projects/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external-projects" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external-projects" )

include_directories(extern/eigen extern/spdlog/include)

# Code Generation
add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/src/generated/lm_derivatives.hpp
    COMMAND python src/python/generate_lm_derivatives.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    DEPENDS src/generated/lm_derivatives.hpp.template src/python/generate_lm_derivatives.py src/python/codegenerator.py 
    COMMENT "Generating code for lm_derivates.cpp."
)
add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/src/generated/Q3.hpp ${PROJECT_SOURCE_DIR}/src/generated/M_C.hpp
    COMMAND python src/python/generate_taubin_determinant.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    DEPENDS src/generated/Q3.hpp.template src/generated/M_C.hpp.template src/python/generate_taubin_determinant.py src/python/codegenerator.py 
    COMMENT "Generating code for Q3.hpp and M_C.hpp."
)

set(SPDLOG_ACTIVE_LEVEL "SPDLOG_LEVEL_OFF" CACHE STRING "Set the logging level (e.g. SPDLOG_LEVEL_DEBUG, or SPDLOG_LEVEL_OFF)")

add_compile_options(-DSPDLOG_ACTIVE_LEVEL=${SPDLOG_ACTIVE_LEVEL})

# Library

add_library(circle_fit STATIC
include/circle_fit/numeric_types.h
include/circle_fit/dataset.h
include/circle_fit/circle_params.h
include/circle_fit/gwaf_taubin.h
src/gwaf_taubin.cpp
include/circle_fit/geo_lm.h
src/geo_lm.cpp
include/circle_fit/circle_fit.h
src/c_interface.cpp
)
set_property(SOURCE src/geo_lm.cpp APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/src/generated/lm_derivatives.hpp)
set_property(SOURCE src/gwaf_taubin.cpp APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_SOURCE_DIR}/src/generated/M_C.hpp ${PROJECT_SOURCE_DIR}/src/generated/Q3.hpp)
set_target_properties(circle_fit PROPERTIES CMAKE_POSITION_INDEPENDENT_CODE ON)
target_link_libraries(circle_fit PRIVATE)
target_include_directories(circle_fit PUBLIC include)

# Also create a dynamic library (E.g. used by python), that contains the same functions.
add_library(circle_fit_dynamic SHARED
include/circle_fit/numeric_types.h
include/circle_fit/dataset.h
include/circle_fit/circle_params.h
include/circle_fit/gwaf_taubin.h
src/gwaf_taubin.cpp
include/circle_fit/geo_lm.h
src/geo_lm.cpp
include/circle_fit/circle_fit.h
src/c_interface.cpp
include/circle_fit/c_interface.h
)
target_include_directories(circle_fit_dynamic PUBLIC include)


set(WITH_TESTS "off" CACHE BOOL "Build tests (requires googletest and GTEST_DIR variable to be set)")
if (WITH_TESTS)
    # Tests
    enable_testing()

    # Common unit tests
    add_executable(test_circle_fit test/test_circle_fit.cpp)
    target_link_libraries(test_circle_fit circle_fit)
    target_include_directories(test_circle_fit PRIVATE extern/Catch2/single_include)
    add_test(NAME unit_tests COMMAND test_circle_fit)

    # GWAF Taubin specific unit tests
    add_executable(test_gwaf_taubin test/test_gwaf_taubin.cpp)
    target_link_libraries(test_gwaf_taubin circle_fit)
    target_include_directories(test_gwaf_taubin PRIVATE extern/Catch2/single_include)
    add_test(NAME gwaf_taubin_tests COMMAND test_gwaf_taubin)

    # Geometric LM specific unit tests
    add_executable(test_geo_lm test/test_geo_lm.cpp)
    target_link_libraries(test_geo_lm circle_fit)
    target_include_directories(test_geo_lm PRIVATE extern/Catch2/single_include)
    add_test(NAME geo_lm_tests COMMAND test_geo_lm)
endif()

set(BUILD_MEX "off" CACHE BOOL "Builds mex files for MATLAB")
if (BUILD_MEX)
    find_package(Matlab)
    # get_property(STATIC_LIB_NAME TARGET circle_fit PROPERTY OUTPUT_LOCATION)
    add_custom_target(circle_fit_mex ALL
        COMMAND ${Matlab_ROOT_DIR}/bin/mex -R2018a -I${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/matlab/mex/circle_fit.c circle_fit${CMAKE_STATIC_LIBRARY_SUFFIX} -output circle_fit_combined
        COMMAND ${Matlab_ROOT_DIR}/bin/mex -R2018a -I${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/matlab/mex/circle_fit_lm.c circle_fit${CMAKE_STATIC_LIBRARY_SUFFIX}
        COMMAND ${Matlab_ROOT_DIR}/bin/mex -R2018a -I${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/matlab/mex/circle_fit_taubin.c circle_fit${CMAKE_STATIC_LIBRARY_SUFFIX}
        DEPENDS circle_fit
    )
endif()
