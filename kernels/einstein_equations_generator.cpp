#include "Halide.h"
#include "generated_expressions.hpp"

#include <cstdio>

using namespace Halide;

auto make_g_dd(Expr z, Func Q, Expr L, Expr mu) -> Func {
    Var _mu{"mu"}, _nu{"nu"};
    Func output{"g_dd"};
    output(_mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_g_dd(output);
    return output;
}

auto make_Dg_ddd(Expr z, Func Q, Func DQ, Expr L, Expr mu) -> Func {
    Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
    Func output{"Dg_ddd"};
    output(_lambda, _mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_Dg_ddd(output);
    return output;
}

auto make_DDg_dddd(Expr z, Func Q, Func DQ, Func DDQ, Expr L, Expr mu) -> Func {
    Var _kappa{"kappa"}, _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
    Func output{"DDg_dddd"};
    output(_kappa, _lambda, _mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_DDg_dddd(output);
    return output;
}

auto make_g_UU(Expr z, Func Q, Expr L, Expr mu) -> Func {
    Var _mu{"mu"}, _nu{"nu"};
    Func output{"g_UU"};
    output(_mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_g_UU(output);
    return output;
}

auto make_Dg_dUU(Expr z, Func Q, Func DQ, Expr L, Expr mu) -> Func {
    Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
    Func output{"Dg_dUU"};
    output(_lambda, _mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_Dg_dUU(output);
    return output;
}

auto make_Gamma_Udd(Func g_UU, Func Dg_ddd) -> Func {
    Var mu{"mu"}, nu{"nu"}, lambda{"lambda"};
    RDom rho{0, 4, "rho"};

    Func Gamma_Udd;
    Gamma_Udd(lambda, mu, nu) =
        0.5f *
        sum(g_UU(lambda, rho) * (Dg_ddd(mu, nu, rho) + Dg_ddd(nu, mu, rho) - Dg_ddd(rho, mu, nu)));
    return Gamma_Udd;
}

auto make_DGamma_dUdd(Func g_UU, Func Dg_dUU, Func Dg_ddd, Func DDg_dddd) -> Func {
    Var mu{"mu"}, nu{"nu"}, lambda{"lambda"}, kappa{"kappa"};
    RDom rho{0, 4, "rho"};

    Func DGamma_dUdd;
    DGamma_dUdd(kappa, lambda, mu, nu) =
        0.5f *
        sum(Dg_dUU(kappa, lambda, rho) *
                (Dg_ddd(mu, nu, rho) + Dg_ddd(nu, mu, rho) - Dg_ddd(rho, mu, nu)) +
            g_UU(lambda, rho) * (DDg_dddd(kappa, mu, nu, rho) + DDg_dddd(kappa, nu, mu, rho) -
                                 DDg_dddd(kappa, rho, mu, nu)));
    return DGamma_dUdd;
}

auto make_Xi_U(Func g_UU, Func Gamma_Udd, Func Gamma_Udd_ref) -> Func {
    Var nu{"nu"}, mu{"mu"};
    RDom lambda{0, 4, "lambda"};
    Func Xi_U{"Xi_U"};
    Func temp;
    temp(mu, nu) =
        sum(g_UU(nu, lambda) * (Gamma_Udd(mu, nu, lambda) - Gamma_Udd_ref(mu, nu, lambda)));
    Xi_U(mu) = sum(temp(mu, lambda));
    return Xi_U;
}

auto make_DXi_dU(Func g_UU, Func Dg_dUU, Func Gamma_Udd, Func DGamma_dUdd, Func Gamma_Udd_ref,
                 Func DGamma_dUdd_ref) -> Func {
    Var mu{"mu"}, rho{"rho"};
    RDom nu{0, 4, "nu"};
    RDom lambda{0, 4, "lambda"};
    Func output;
    output(rho, mu) =
        sum(g_UU(nu, lambda) *
                (DGamma_dUdd(rho, mu, nu, lambda) - DGamma_dUdd_ref(rho, mu, nu, lambda)) +
            Dg_dUU(rho, nu, lambda) * (Gamma_Udd(mu, nu, lambda) - Gamma_Udd_ref(mu, nu, lambda)));
    return output;
}

auto make_Xi_d(Func g_dd, Func Xi_U) -> Func {
    Var mu{"mu"};
    RDom nu{0, 4, "nu"};
    Func output;
    output(mu) = g_dd(mu, nu) * Xi_U(nu);
    return output;
}

auto make_DXi_dd(Func g_dd, Func Dg_ddd, Func Xi_U, Func DXi_dU) -> Func {
    Var lambda{"lambda"}, mu{"mu"};
    RDom nu{0, 4, "nu"};
    Func output;
    output(lambda, mu) = Dg_ddd(lambda, mu, nu) * Xi_U(nu) + g_dd(mu, nu) * DXi_dU(lambda, nu);
    return output;
}

auto covariant_derivative(Func v_d, Func Dv_dd, Func Gamma_Udd) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom lambda{0, 4, "lambda"};
    Func output;
    output(mu, nu) = Dv_dd(mu, nu) - Gamma_Udd(lambda, nu, mu) * v_d(lambda);
    return output;
}

auto make_Q_ref() -> Func {
    Var i;
    Func output;
    output(i) = cast<double>(1);
    output(4) = cast<double>(0);
    return output;
}

auto make_DQ_ref() -> Func {
    Var i, j;
    Func output;
    output(i, j) = cast<double>(0);
    return output;
}

auto make_DDQ_ref() -> Func {
    Var i, j, k;
    Func output;
    output(i, j, k) = cast<double>(0);
    return output;
}

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

class g_dd_generator : public Halide::Generator<g_dd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_g_dd(z, Q, length, chemical_potential);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4);
    }
};
class g_UU_generator : public Halide::Generator<g_UU_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        output = make_g_UU(z, Q, length, chemical_potential);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4);
    }
};
class Gamma_Udd_generator : public Halide::Generator<Gamma_Udd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 3};

    void generate() {
        auto g_UU = make_g_UU(z, Q, length, chemical_potential);
        auto Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        output = make_Gamma_Udd(g_UU, Dg_ddd);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4).dim(2).set_bounds(0, 4);
    }
};
class Xi_U_generator : public Halide::Generator<Xi_U_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 1};

    void generate() {
        auto g_UU = make_g_UU(z, Q, length, chemical_potential);
        auto Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        auto Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
        auto Q_ref = make_Q_ref();
        auto DQ_ref = make_DQ_ref();
        auto Gamma_Udd_ref =
            make_Gamma_Udd(make_g_UU(z, Q_ref, length, chemical_potential),
                           make_Dg_ddd(z, Q_ref, DQ_ref, length, chemical_potential));
        output = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
        output.dim(0).set_bounds(0, 4);
    }
};
class DXi_dU_generator : public Halide::Generator<DXi_dU_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DQ", 3};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 1};

    void generate() {
        auto g_UU = make_g_UU(z, Q, length, chemical_potential);
        auto Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        auto Dg_dUU = make_Dg_dUU(z, Q, DQ, length, chemical_potential);
        auto DDg_dddd = make_DDg_dddd(z, Q, DQ, DDQ, length, chemical_potential);
        auto Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
        auto DGamma_dUdd = make_DGamma_dUdd(g_UU, Dg_dUU, Dg_ddd, DDg_dddd);

        auto Q_ref = make_Q_ref();
        auto DQ_ref = make_DQ_ref();
        auto DDQ_ref = make_DDQ_ref();
        auto g_UU_ref = make_g_UU(z, Q_ref, length, chemical_potential);
        auto Dg_ddd_ref = make_Dg_ddd(z, Q_ref, DQ_ref, length, chemical_potential);
        auto Dg_dUU_ref = make_Dg_dUU(z, Q_ref, DQ_ref, length, chemical_potential);
        auto DDg_dddd_ref = make_DDg_dddd(z, Q_ref, DQ_ref, DDQ_ref, length, chemical_potential);
        auto Gamma_Udd_ref = make_Gamma_Udd(g_UU_ref, Dg_ddd_ref);
        auto DGamma_dUdd_ref = make_DGamma_dUdd(g_UU_ref, Dg_dUU_ref, Dg_ddd_ref, DDg_dddd_ref);

        output = make_DXi_dU(g_UU, Dg_dUU, Gamma_Udd, DGamma_dUdd, Gamma_Udd_ref, DGamma_dUdd_ref);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4);
    }
};

class g_dd_ref_generator : public Halide::Generator<g_dd_ref_generator> {
  public:
    Input<double> z{"z"};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        Var mu;
        Func Q;
        Q(mu) = cast<double>(1);
        Q(4) = cast<double>(0);
        output = make_g_dd(z, Q, length, chemical_potential);
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

HALIDE_REGISTER_GENERATOR(g_dd_generator, g_dd_generator)
HALIDE_REGISTER_GENERATOR(g_UU_generator, g_UU_generator)
HALIDE_REGISTER_GENERATOR(Gamma_Udd_generator, Gamma_Udd_generator)
HALIDE_REGISTER_GENERATOR(Xi_U_generator, Xi_U_generator)

HALIDE_REGISTER_GENERATOR(g_dd_ref_generator, g_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(g_UU_ref_generator, g_UU_ref_generator)
HALIDE_REGISTER_GENERATOR(Dg_dd_ref_generator, Dg_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(Dg_UU_ref_generator, Dg_UU_ref_generator)
HALIDE_REGISTER_GENERATOR(DDg_dd_ref_generator, DDg_dd_ref_generator)
HALIDE_REGISTER_GENERATOR(Gamma_Udd_ref_generator, Gamma_Udd_ref_generator)
HALIDE_REGISTER_GENERATOR(DGamma_dUdd_ref_generator, DGamma_dUdd_ref_generator)
