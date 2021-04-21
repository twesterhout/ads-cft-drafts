#include "Halide.h"

#include <cstdio>

using namespace Halide;

template <class InputDouble>
auto make_g_dd_ref(Expr z, InputDouble const& L, InputDouble const& mu) -> Func {
    auto const z2 = z * z;
    auto const z3 = z2 * z;
    auto const z4 = z2 * z2;
    auto const pre = L * L / z2;
    auto const mu2 = mu * mu;

    Var _mu{"mu"}, _nu{"nu"};
    Func g_dd{"g_dd_ref"};
    g_dd(_mu, _nu) = cast<double>(0);
    g_dd(0, 0) = pre * (z - 1) * (1 + z + z2 - z3 * mu2 / 2);
    g_dd(1, 1) = pre;
    g_dd(2, 2) = pre;
    g_dd(3, 3) = pre * 2 / (2 + z4 * mu2 - z3 * (2 + mu2));
    return g_dd;
}

class g_dd_ref_generator : public Halide::Generator<g_dd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_g_dd_ref(z, length, chemical_potential);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4);
        // clang-format on
    }
};

template <class InputDouble>
auto make_g_UU_ref(Expr z, InputDouble const& L, InputDouble const& mu) -> Func {
    auto z2 = z * z;
    auto z3 = z2 * z;
    auto z4 = z2 * z2;
    auto pre = z2 / (L * L);
    auto mu2 = mu * mu;

    Var _mu{"mu"}, _nu{"nu"};
    Func g_UU{"g_UU_ref"};
    g_UU(_mu, _nu) = cast<double>(0);
    g_UU(0, 0) = -2 * pre / (2 + z4 * mu2 - z3 * (2 + mu2));
    g_UU(1, 1) = pre;
    g_UU(2, 2) = pre;
    g_UU(3, 3) = pre * (2 + z4 * mu2 - z3 * (2 + mu2)) / 2;
    return g_UU;
}

class g_UU_ref_generator : public Halide::Generator<g_UU_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_g_UU_ref(z, length, chemical_potential);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4);
        // clang-format on
    }
};

template <class InputDouble>
auto make_Dg_dd_ref(Expr z, InputDouble const& L, InputDouble const& mu) -> Func {
    auto z2 = z * z;
    auto z3 = z2 * z;
    auto z4 = z2 * z2;
    auto pre = L * L / z3;
    auto mu2 = mu * mu;

    Var _mu{"mu"}, _nu{"nu"};
    Func Dg_dd{"Dg_dd_ref"};
    Dg_dd(_mu, _nu) = cast<double>(0);
    Dg_dd(0, 0) = pre * (4 - 2 * z4 * mu2 + z3 * (2 + mu2)) / 2;
    Dg_dd(1, 1) = -2 * pre;
    Dg_dd(2, 2) = -2 * pre;
    Dg_dd(3, 3) =
        -2 * pre * (4 + 6 * z4 * mu2 - 5 * z3 * (2 + mu2)) / pow(2 + z4 * mu2 - z3 * (2 + mu2), 2);
    return Dg_dd;
}

class Dg_dd_ref_generator : public Halide::Generator<Dg_dd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_Dg_dd_ref(z, length, chemical_potential);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4);
        // clang-format on
    }
};

template <class InputDouble>
auto make_DDg_dd_ref(Expr z, InputDouble const& L, InputDouble const& mu) -> Func {
    auto z2 = z * z;
    auto z3 = z2 * z;
    auto z4 = z2 * z2;
    auto pre = L * L / z3;
    auto mu2 = mu * mu;

    Var _mu{"mu"}, _nu{"nu"};
    Func Dg_dd{"Dg_dd_ref"};
    Dg_dd(_mu, _nu) = cast<double>(0);
    Dg_dd(0, 0) = pre * (4 - 2 * z4 * mu2 + z3 * (2 + mu2)) / 2;
    Dg_dd(1, 1) = -2 * pre;
    Dg_dd(2, 2) = -2 * pre;
    Dg_dd(3, 3) =
        -2 * pre * (4 + 6 * z4 * mu2 - 5 * z3 * (2 + mu2)) / pow(2 + z4 * mu2 - z3 * (2 + mu2), 2);
    return Dg_dd;
}

class DDg_dd_ref_generator : public Halide::Generator<DDg_dd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_DDg_dd_ref(z, length, chemical_potential);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4);
        // clang-format on
    }
};

auto make_Gamma_Udd(Func g_UU, Func Dg_dd) -> Func {
    Var mu{"mu"}, nu{"nu"}, lambda{"lambda"};
    RDom rho{0, 4, "rho"};

    Func Gamma_Udd;
    Gamma_Udd(lambda, mu, nu) =
        0.5f *
        sum(g_UU(lambda, rho) * (Dg_dd(mu, nu, rho) + Dg_dd(nu, mu, rho) - Dg_dd(rho, mu, nu)));
    return Gamma_Udd;
}

class Gamma_Udd_ref_generator : public Halide::Generator<Gamma_Udd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 3};

    void generate() {
        Var lambda, mu, nu;

        auto g_UU = make_g_UU_ref(z, length, chemical_potential);
        Func Dg_dd;
        Dg_dd(lambda, mu, nu) = cast<double>(0);
        Dg_dd(3, mu, nu) = make_Dg_dd_ref(z, length, chemical_potential)(mu, nu);

        output = make_Gamma_Udd(g_UU, Dg_dd);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4)
              .dim(2).set_bounds(0, 4);
        // clang-format on
    }
};

class EinsteinEquationsGenerator : public Halide::Generator<EinsteinEquationsGenerator> {
  public:
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Input<Buffer<double>> grid{"grid", 1};
    Output<Buffer<double>> gdd{"gdd", 3};

    // Typically you declare your Vars at this scope as well, so that
    // they can be used in any helper methods you add later.
    Var x, y, zi, mu, nu;

    Func make_gdd_ref() {
        auto z = grid(zi);
        auto z2 = z * z;
        auto z3 = z2 * z;
        auto z4 = z2 * z2;
        auto pre = length * length / z2;
        auto chemical_potential2 = chemical_potential * chemical_potential;

        Func gdd_ref;
        gdd_ref(zi, mu, nu) = cast<double>(0);
        gdd_ref(zi, 0, 0) = pre * (z - 1) * (1 + z + z2 - z3 * chemical_potential2 / 2);
        gdd_ref(zi, 1, 1) = pre;
        gdd_ref(zi, 2, 2) = pre;
        gdd_ref(zi, 3, 3) =
            pre * 2 / (2 + z4 * chemical_potential2 - z3 * (2 + chemical_potential2));
        return gdd_ref;
    }

    Func make_Dgdd_ref() {
        auto z = grid(zi);
        auto z2 = z * z;
        auto z3 = z2 * z;
        auto z4 = z2 * z2;
        auto pre = length * length / z3;
        auto chemical_potential2 = chemical_potential * chemical_potential;

        Func Dgdd_ref;
        Dgdd_ref(zi, mu, nu) = cast<double>(0);
        Dgdd_ref(zi, 0, 0) =
            pre * (4 - 2 * z4 * chemical_potential2 + z3 * (2 + chemical_potential2)) / 2;
        Dgdd_ref(zi, 1, 1) = -2 * pre;
        Dgdd_ref(zi, 2, 2) = -2 * pre;
        Dgdd_ref(zi, 3, 3) =
            -2 * pre * (4 + 6 * z4 * chemical_potential2 - 5 * z3 * (2 + chemical_potential2)) /
            pow(2 + z4 * chemical_potential2 - z3 * (2 + chemical_potential2), 2);
        return Dgdd_ref;
    }

    Func make_Gamma_Udd_ref() {
        Var lambda;
        RDom rho{0, 4};
        Func Gamma_Udd;
        Func gUU;
        Gamma_Udd(lambda, mu, nu) = cast<double>(0);
        Gamma_Udd(lambda, mu, nu) += 0.5f * gUU(lambda, mu);

        return Gamma_Udd;
    }

    Func make_DDgdd_ref() {
        auto z = grid(zi);
        auto z2 = z * z;
        auto z3 = z2 * z;
        auto z4 = z2 * z2;
        auto z6 = z3 * z3;
        auto z7 = z4 * z3;
        auto z8 = z4 * z4;
        auto pre = length * length / z4;
        auto chemical_potential2 = chemical_potential * chemical_potential;
        auto chemical_potential4 = chemical_potential2 * chemical_potential2;

        Func DDgdd_ref;
        DDgdd_ref(zi, mu, nu) = cast<double>(0);
        DDgdd_ref(zi, 0, 0) = -pre * (6 + z4 * chemical_potential2);
        DDgdd_ref(zi, 1, 1) = 6 * pre;
        DDgdd_ref(zi, 2, 2) = 6 * pre;
        DDgdd_ref(zi, 3, 3) = 4 * pre *
                              (12 + 16 * z4 * chemical_potential2 + 21 * z8 * chemical_potential4 -
                               18 * z3 * (2 + chemical_potential2) -
                               35 * z7 * chemical_potential2 * (2 + chemical_potential2) +
                               15 * z6 * pow(2 + chemical_potential2, 2)) /
                              pow(2 + z4 * chemical_potential2 - z3 * (2 + chemical_potential2), 3);
        return DDgdd_ref;
    }

    void generate() {
        gdd = make_gdd_ref();
        gdd.dim(0).set_min(0).dim(1).set_bounds(0, 4).dim(2).set_bounds(0, 4);
    }
};

// We compile this file along with tools/GenGen.cpp. That file defines
// an "int main(...)" that provides the command-line interface to use
// your generator class. We need to tell that code about our
// generator. We do this like so:
HALIDE_REGISTER_GENERATOR(EinsteinEquationsGenerator, einstein_equations_generator)
HALIDE_REGISTER_GENERATOR(g_dd_ref_generator, g_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(g_UU_ref_generator, g_UU_ref_generator)
HALIDE_REGISTER_GENERATOR(Dg_dd_ref_generator, Dg_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(Gamma_Udd_ref_generator, Gamma_Udd_ref_generator)
