#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "DDg_dd_ref.h"
#include "DGamma_dUdd_ref.h"
#include "Dg_UU_ref.h"
#include "Dg_dd_ref.h"
#include "Gamma_Udd.h"
#include "Gamma_Udd_ref.h"
#include "Xi_U.h"
#include "g_UU.h"
#include "g_UU_ref.h"
#include "g_dd.h"
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

template <size_t... Ns> struct looper;
template <size_t N> struct looper<N> {
    template <class Function> auto operator()(Function fn) const -> void {
        for (auto i = 0; i < N; ++i) {
            fn(i);
        }
    }
};
template <size_t N, size_t... Ns> struct looper<N, Ns...> {
    template <class Function> auto operator()(Function fn) const -> void {
        looper<N>{}([&](auto const i) { looper<Ns...>{}([&](auto... is) { fn(i, is...); }); });
    }
};
template <size_t... Ns, class Function> auto loop(Function fn) -> void {
    static_assert(sizeof...(Ns) > 0);
    looper<Ns...>{}(fn);
}

auto read_double(std::FILE* fp) -> double {
    double element;
    if (std::fscanf(fp, "%lf", &element) != 1) { throw std::runtime_error{"Parsing failed"}; }
    return element;
}

template <size_t... Ns> auto load_tensor(std::FILE* fp) -> Halide::Runtime::Buffer<double> {
    Halide::Runtime::Buffer<double> buffer{Ns...};
    loop<Ns...>([&](auto... is) { buffer(is...) = read_double(fp); });
    return buffer;
}
template <size_t... Ns> auto load_tensor(char const* filename) -> Halide::Runtime::Buffer<double> {
    auto fp = open_file(filename);
    return load_tensor<Ns...>(fp.get());
}

auto obtain_filenames(char const* basename, int const test_case) {
    char buffer[128];
    std::snprintf(buffer, std::size(buffer), "%s/input_%02d.dat", basename, test_case);
    std::string input{buffer};
    std::snprintf(buffer, std::size(buffer), "%s/output_%02d.dat", basename, test_case);
    std::string output{buffer};
    return std::make_tuple(std::move(input), std::move(output));
}

auto read_parameters(char const* filename) {
    auto fp = open_file(filename);
    Halide::Runtime::Buffer<double> DDQ{2, 2, 5};
    Halide::Runtime::Buffer<double> DQ{2, 5};
    Halide::Runtime::Buffer<double> Q{5};
    for (auto i = 0; i < 5; ++i) {
        DDQ(0, 0, i) = read_double(fp.get());
        DDQ(0, 1, i) = read_double(fp.get());
        DDQ(1, 0, i) = DDQ(0, 1, i);
        DDQ(1, 1, i) = read_double(fp.get());
        DQ(0, i) = read_double(fp.get());
        DQ(1, i) = read_double(fp.get());
        Q(i) = read_double(fp.get());
    }
    auto z = read_double(fp.get());
    auto mu = read_double(fp.get());
    auto L = read_double(fp.get());
    return std::tuple{z, std::move(Q), std::move(DQ), std::move(DDQ), L, mu};
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
        auto const expected = load_tensor<4, 4>(output.c_str());

        Buffer<double> predicted{4, 4};
        fn(z, L, mu, predicted);
        for (auto i = 0; i < 4; ++i) {
            for (auto j = 0; j < 4; ++j) {
                REQUIRE(predicted(i, j) == doctest::Approx(expected(i, j)));
            }
        }
    }
}

TEST_CASE("testing g_dd function") {
    auto [z, Q, DQ, DDQ, L, mu] = read_parameters("test/input_01.dat");
    auto const expected = load_tensor<4, 4>("test/g_dd/output_01.dat");
    Buffer<double> predicted{4, 4};
    g_dd(z, Q, L, mu, predicted);
    loop<4, 4>([&](auto... is) { REQUIRE(predicted(is...) == doctest::Approx(expected(is...))); });
}
TEST_CASE("testing g_UU function") {
    auto [z, Q, DQ, DDQ, L, mu] = read_parameters("test/input_01.dat");
    auto const expected = load_tensor<4, 4>("test/g_UU/output_01.dat");
    Buffer<double> predicted{4, 4};
    g_UU(z, Q, L, mu, predicted);
    loop<4, 4>([&](auto... is) { REQUIRE(predicted(is...) == doctest::Approx(expected(is...))); });
}
TEST_CASE("testing Gamma_Udd function") {
    auto [z, Q, DQ, DDQ, L, mu] = read_parameters("test/input_01.dat");
    auto const expected = load_tensor<4, 4, 4>("test/Gamma_Udd/output_01.dat");
    Buffer<double> predicted{4, 4, 4};
    Gamma_Udd(z, Q, DQ, L, mu, predicted);
    loop<4, 4, 4>(
        [&](auto... is) { REQUIRE(predicted(is...) == doctest::Approx(expected(is...))); });
}
TEST_CASE("testing Xi_U function") {
    auto [z, Q, DQ, DDQ, L, mu] = read_parameters("test/input_01.dat");
    auto const expected = load_tensor<4>("test/Xi_U/output_01.dat");
    Buffer<double> predicted{4};
    Xi_U(z, Q, DQ, L, mu, predicted);
    loop<4>([&](auto... is) { REQUIRE(predicted(is...) == doctest::Approx(expected(is...))); });
}

TEST_CASE("testing g_dd_ref function") { test_ref_function("test/g_dd_ref", g_dd_ref); }
TEST_CASE("testing g_UU_ref function") { test_ref_function("test/g_UU_ref", g_UU_ref); }
TEST_CASE("testing Dg_dd_ref function") { test_ref_function("test/Dg_dd_ref", Dg_dd_ref); }
TEST_CASE("testing Dg_UU_ref function") { test_ref_function("test/Dg_UU_ref", Dg_UU_ref); }
TEST_CASE("testing DDg_dd_ref function") { test_ref_function("test/DDg_dd_ref", DDg_dd_ref); }
TEST_CASE("testing Gamma_Udd_ref function") {
    for (auto test_case = 1; test_case <= 1; ++test_case) {
        auto const& [input, output] = obtain_filenames("test/Gamma_Udd_ref", test_case);
        auto const [L, mu, z] = read_input_ref(input.c_str());
        auto const expected = load_tensor<4, 4, 4>(output.c_str());

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
TEST_CASE("testing DGamma_dUdd_ref function") {
    for (auto test_case = 1; test_case <= 3; ++test_case) {
        auto const& [input, output] = obtain_filenames("test/DGamma_dUdd_ref", test_case);
        auto const [L, mu, z] = read_input_ref(input.c_str());
        auto const expected = load_tensor<4, 4, 4>(output.c_str());

        Buffer<double> predicted{4, 4, 4, 4};
        DGamma_dUdd_ref(z, L, mu, predicted);
        for (auto i = 0; i < 4; ++i) {
            for (auto j = 0; j < 4; ++j) {
                for (auto k = 0; k < 4; ++k) {
                    for (auto l = 0; l < 4; ++l) {
                        if (i == 3) {
                            REQUIRE(predicted(i, j, k, l) == doctest::Approx(expected(j, k, l)));
                        } else {
                            REQUIRE(predicted(i, j, k, l) == doctest::Approx(0.0));
                        }
                    }
                }
            }
        }
    }
}
