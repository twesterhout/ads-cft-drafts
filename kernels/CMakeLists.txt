cmake_minimum_required(VERSION 3.16)
project(lesson_01)

# Set up language settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

# Find Halide
find_package(Halide REQUIRED)
find_package(OpenMP)

# Generator
add_executable(lesson_01 lesson_01.cpp)
target_link_libraries(lesson_01 PRIVATE Halide::Halide)

add_executable(lesson_15.generator lesson_15.cpp)
target_link_libraries(lesson_15.generator PRIVATE Halide::Generator)

add_executable(einstein_equations.generator einstein_equations_generator.cpp)
target_link_libraries(einstein_equations.generator PRIVATE Halide::Generator)

add_halide_library(g_dd_ref
    FROM einstein_equations.generator GENERATOR g_dd_ref_generator)
add_halide_library(g_UU_ref
    FROM einstein_equations.generator GENERATOR g_UU_ref_generator)
add_halide_library(Dg_dd_ref
    FROM einstein_equations.generator GENERATOR Dg_dd_ref_generator)
add_halide_library(Dg_UU_ref
    FROM einstein_equations.generator GENERATOR Dg_UU_ref_generator)
add_halide_library(DDg_dd_ref
    FROM einstein_equations.generator GENERATOR DDg_dd_ref_generator)
add_halide_library(Gamma_Udd_ref
    FROM einstein_equations.generator GENERATOR Gamma_Udd_ref_generator)
add_halide_library(DGamma_dUdd_ref
    FROM einstein_equations.generator GENERATOR DGamma_dUdd_ref_generator)

add_executable(test_generator test_generator.cpp)
target_link_libraries(test_generator
    PRIVATE g_dd_ref g_UU_ref Dg_dd_ref Dg_UU_ref DDg_dd_ref Gamma_Udd_ref DGamma_dUdd_ref)
target_include_directories(test_generator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../third_party)
