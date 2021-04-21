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

auto make_Dg_UU_ref(Expr z, Expr L, Expr mu) -> Func {
    auto z2 = z * z;
    auto z3 = z2 * z;
    auto z4 = z2 * z2;
    auto pre = z / (L * L);
    auto mu2 = mu * mu;

    Var _mu{"mu"}, _nu{"nu"};
    Func Dg_UU{"Dg_UU_ref"};
    Dg_UU(_mu, _nu) = cast<double>(0);
    Dg_UU(0, 0) = pre * (-8 + 4 * z4 * mu2 - 2 * z3 * (2 + mu2)) / pow(z - 1, 2) /
                  pow(2 + 2 * z + 2 * z2 - z3 * mu2, 2);
    Dg_UU(1, 1) = 2 * pre;
    Dg_UU(2, 2) = 2 * pre;
    Dg_UU(3, 3) = pre * (4 + 6 * z4 * mu2 - 5 * z3 * (2 + mu2)) / 2;
    return Dg_UU;
}

class Dg_UU_ref_generator : public Halide::Generator<Dg_UU_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_Dg_UU_ref(z, length, chemical_potential);
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
    auto z6 = z3 * z3;
    auto z7 = z4 * z3;
    auto z8 = z4 * z4;
    auto pre = L * L / z4;
    auto mu2 = mu * mu;
    auto mu4 = mu2 * mu2;

    Var _mu{"mu"}, _nu{"nu"};
    Func DDg_dd{"Dg_dd_ref"};
    DDg_dd(_mu, _nu) = cast<double>(0);
    DDg_dd(0, 0) = -pre * (6 + z4 * mu2);
    DDg_dd(1, 1) = 6 * pre;
    DDg_dd(2, 2) = 6 * pre;
    DDg_dd(3, 3) = 4 * pre *
                   (12 + 16 * z4 * mu2 + 21 * z8 * mu4 - 18 * z3 * (2 + mu2) -
                    35 * z7 * mu2 * (2 + mu2) + 15 * z6 * pow(2 + mu2, 2)) /
                   pow(2 + z4 * mu2 - z3 * (2 + mu2), 3);
    return DDg_dd;
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

auto make_DGamma_dUdd(Func g_UU, Func Dg_UU, Func Dg_dd, Func DDg_dd) -> Func {
    Var mu{"mu"}, nu{"nu"}, lambda{"lambda"}, kappa{"kappa"};
    RDom rho{0, 4, "rho"};

    Func DGamma_dUdd;
    DGamma_dUdd(kappa, lambda, mu, nu) =
        0.5f * sum(Dg_UU(kappa, lambda, rho) *
                       (Dg_dd(mu, nu, rho) + Dg_dd(nu, mu, rho) - Dg_dd(rho, mu, nu)) +
                   g_UU(lambda, rho) * (DDg_dd(kappa, mu, nu, rho) + DDg_dd(kappa, nu, mu, rho) -
                                        DDg_dd(kappa, rho, mu, nu)));
    return DGamma_dUdd;
}

class DGamma_dUdd_ref_generator : public Halide::Generator<DGamma_dUdd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 4};

    void generate() {
        Var kappa, lambda, mu, nu;

        auto g_UU = make_g_UU_ref(z, length, chemical_potential);
        Func Dg_UU;
        Dg_UU(lambda, mu, nu) = cast<double>(0);
        Dg_UU(3, mu, nu) = make_Dg_UU_ref(z, length, chemical_potential)(mu, nu);

        Func Dg_dd;
        Dg_dd(lambda, mu, nu) = cast<double>(0);
        Dg_dd(3, mu, nu) = make_Dg_dd_ref(z, length, chemical_potential)(mu, nu);
        Func DDg_dd;
        DDg_dd(kappa, lambda, mu, nu) = cast<double>(0);
        DDg_dd(3, 3, mu, nu) = make_DDg_dd_ref(z, length, chemical_potential)(mu, nu);

        output = make_DGamma_dUdd(g_UU, Dg_UU, Dg_dd, DDg_dd);
        // clang-format off
        output.dim(0).set_bounds(0, 4)
              .dim(1).set_bounds(0, 4)
              .dim(2).set_bounds(0, 4)
              .dim(3).set_bounds(0, 4);
        // clang-format on
    }
};

HALIDE_REGISTER_GENERATOR(g_dd_ref_generator, g_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(g_UU_ref_generator, g_UU_ref_generator)
HALIDE_REGISTER_GENERATOR(Dg_dd_ref_generator, Dg_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(Dg_UU_ref_generator, Dg_UU_ref_generator)
HALIDE_REGISTER_GENERATOR(DDg_dd_ref_generator, DDg_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(Gamma_Udd_ref_generator, Gamma_Udd_ref_generator)
HALIDE_REGISTER_GENERATOR(DGamma_dUdd_ref_generator, DGamma_dUdd_ref_generator)
