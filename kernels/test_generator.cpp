#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "Dg_dd_ref.h"
#include "Gamma_Udd_ref.h"
#include "g_UU_ref.h"
#include "g_dd_ref.h"
#include <HalideBuffer.h>
#include <doctest.h>
#include <memory>

using namespace Halide::Runtime;

auto open_file(char const* filename) {
    struct close_file {
        auto operator()(std::FILE* p) const noexcept -> void { std::fclose(p); }
    };
    auto fp = std::unique_ptr<std::FILE, close_file>{std::fopen(filename, "r")};
    if (!fp) { throw std::runtime_error{"Failed to open file for reading"}; }
    return fp;
}

auto load_4x4_matrix(char const* filename) -> Halide::Runtime::Buffer<double> {
    auto fp = open_file(filename);
    Halide::Runtime::Buffer<double> buffer{4, 4};
    for (auto i = 0; i < 4; ++i) {
        for (auto j = 0; j < 4; ++j) {
            double element;
            if (std::fscanf(fp.get(), "%lf", &element) != 1) {
                throw std::runtime_error{"Parsing failed"};
            }
            buffer(i, j) = element;
        }
    }
    return buffer;
}

auto load_4x4x4(char const* filename) -> Halide::Runtime::Buffer<double> {
    auto fp = open_file(filename);
    Halide::Runtime::Buffer<double> buffer{4, 4, 4};
    for (auto i = 0; i < 4; ++i) {
        for (auto j = 0; j < 4; ++j) {
            for (auto k = 0; k < 4; ++k) {
                double element;
                if (std::fscanf(fp.get(), "%lf", &element) != 1) {
                    throw std::runtime_error{"Parsing failed"};
                }
                buffer(i, j, k) = element;
            }
        }
    }
    return buffer;
}

auto obtain_filenames(char const* basename, int const test_case) {
    char buffer[128];
    std::snprintf(buffer, std::size(buffer), "%s/input_%02d.dat", basename, test_case);
    std::string input{buffer};
    std::snprintf(buffer, std::size(buffer), "%s/output_%02d.dat", basename, test_case);
    std::string output{buffer};
    return std::make_tuple(std::move(input), std::move(output));
}

auto read_input_ref(char const* filename) {
    auto fp = open_file(filename);
    double L, mu, z;
    if (std::fscanf(fp.get(), "%lf\t%lf\t%lf", &L, &mu, &z) != 3) {
        throw std::runtime_error{"Parsing failed"};
    }
    return std::tuple{L, mu, z};
}

template <class Function> auto test_ref_function(char const* basename, Function fn) -> void {
    for (auto test_case = 1; test_case <= 3; ++test_case) {
        auto const& [input, output] = obtain_filenames(basename, test_case);
        auto const [L, mu, z] = read_input_ref(input.c_str());
        auto const expected = load_4x4_matrix(output.c_str());

        Buffer<double> predicted{4, 4};
        fn(z, L, mu, predicted);
        for (auto i = 0; i < 4; ++i) {
            for (auto j = 0; j < 4; ++j) {
                REQUIRE(predicted(i, j) == doctest::Approx(expected(i, j)));
            }
        }
    }
}

TEST_CASE("testing g_dd_ref function") { test_ref_function("test/g_dd_ref", g_dd_ref); }
TEST_CASE("testing g_UU_ref function") { test_ref_function("test/g_UU_ref", g_UU_ref); }
TEST_CASE("testing Dg_dd_ref function") { test_ref_function("test/Dg_dd_ref", Dg_dd_ref); }

TEST_CASE("testing Gamma_Udd_ref function") {
    for (auto test_case = 1; test_case <= 1; ++test_case) {
        auto const& [input, output] = obtain_filenames("test/Gamma_Udd_ref", test_case);
        auto const [L, mu, z] = read_input_ref(input.c_str());
        auto const expected = load_4x4x4(output.c_str());

        Buffer<double> predicted{4, 4, 4};
        Gamma_Udd_ref(z, L, mu, predicted);
        for (auto i = 0; i < 4; ++i) {
            for (auto j = 0; j < 4; ++j) {
                for (auto k = 0; k < 4; ++k) {
                    REQUIRE(predicted(i, j, k) == doctest::Approx(expected(i, j, k)));
                }
            }
        }
    }
}
