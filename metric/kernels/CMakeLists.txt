cmake_minimum_required(VERSION 3.16)
project(ads-cft-kernels)

# Find Halide
find_package(Halide REQUIRED)

_Halide_get_platform_details(
    is_crosscompiling
    object_suffix
    static_library_suffix)

add_executable(einstein.generator einstein_equations.cpp)
target_link_libraries(einstein.generator PRIVATE Halide::Generator)

# We build a custom runtime where
#   1) ${Halide_HOST_TARGET} is used rather than 'host'.
#   2) a normal name is used rather than 'something.runtime.a' (either Cabal or
#   Nix can't handle dots in static library names)
add_custom_command(
    OUTPUT runtime.o
    COMMAND einstein.generator -r "runtime" -e object -o .
                                target=${Halide_HOST_TARGET}
    DEPENDS einstein.generator
    VERBATIM)
add_library(halide_runtime runtime${object_suffix})
set_target_properties(halide_runtime PROPERTIES LINKER_LANGUAGE CXX)

# Generate the kernels for 'host' using our custom runtime
add_halide_library(kernels_metric
    FROM einstein.generator
    GENERATOR metric_generator
    FUNCTION_NAME ads_cft_halide_compute_metric
    USE_RUNTIME halide_runtime
    TARGETS host)

add_halide_library(kernels_christoffel
    FROM einstein.generator
    GENERATOR christoffel_generator
    FUNCTION_NAME ads_cft_halide_compute_christoffel
    USE_RUNTIME halide_runtime
    TARGETS host)

add_halide_library(kernels_maxwell
    FROM einstein.generator
    GENERATOR maxwell_generator
    FUNCTION_NAME ads_cft_halide_compute_maxwell
    USE_RUNTIME halide_runtime
    TARGETS host)

add_halide_library(kernels_deturck
    FROM einstein.generator
    GENERATOR deturck_generator
    FUNCTION_NAME ads_cft_halide_compute_deturck
    USE_RUNTIME halide_runtime
    TARGETS host)

add_halide_library(kernels_equations
    FROM einstein.generator
    GENERATOR equations_generator
    FUNCTION_NAME ads_cft_halide_evaluate_equations
    USE_RUNTIME halide_runtime
    TARGETS host)

add_halide_library(kernels_horizon
    FROM einstein.generator
    GENERATOR horizon_boundary_generator
    FUNCTION_NAME ads_cft_halide_horizon_boundary_conditions
    USE_RUNTIME halide_runtime
    TARGETS host)

# add_halide_library(kernels_foo
#     FROM einstein.generator
#     GENERATOR foo_generator
#     FUNCTION_NAME ads_cft_halide_foo
#     USE_RUNTIME halide_runtime
#     TARGETS host)

# Generate a combined library with both the kernel and runtime.
# If shared library is used, we can load it into ghci.
add_library(ads_cft_kernels kernels_metric${object_suffix}
                            kernels_maxwell${object_suffix}
                            kernels_christoffel${object_suffix}
                            kernels_deturck${object_suffix}
                            kernels_equations${object_suffix}
                            kernels_horizon${object_suffix}
                            runtime${object_suffix})
set_target_properties(ads_cft_kernels PROPERTIES LINKER_LANGUAGE CXX)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_metric.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_maxwell.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_christoffel.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_deturck.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_equations.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_horizon.h DESTINATION include)
install(TARGETS ads_cft_kernels DESTINATION lib)
install(TARGETS einstein.generator DESTINATION bin)
