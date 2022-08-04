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

# Generate a combined library with both the kernel and runtime.
# If shared library is used, we can load it into ghci.
add_library(ads_cft_kernels kernels_metric${object_suffix}
                            runtime${object_suffix})
set_target_properties(ads_cft_kernels PROPERTIES LINKER_LANGUAGE CXX)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/kernels_metric.h DESTINATION include)
install(TARGETS ads_cft_kernels DESTINATION lib)
install(TARGETS einstein.generator DESTINATION bin)
