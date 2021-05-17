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
    Var mu{"mu"};
    RDom lambda{0, 4, "lambda"};
    RDom nu{0, 4, "nu"};
    Func Xi_U{"Xi_U"};
    auto expr = g_UU(nu, lambda) * (Gamma_Udd(mu, nu, lambda) - Gamma_Udd_ref(mu, nu, lambda));
    Xi_U(mu) = sum(nu, sum(lambda, expr));
    return Xi_U;
}

auto make_DXi_dU(Func g_UU, Func Dg_dUU, Func Gamma_Udd, Func DGamma_dUdd, Func Gamma_Udd_ref,
                 Func DGamma_dUdd_ref) -> Func {
    Var mu{"mu"}, rho{"rho"};
    RDom nu{0, 4, "nu"};
    RDom lambda{0, 4, "lambda"};
    Func output;
    auto expr =
        g_UU(nu, lambda) *
            (DGamma_dUdd(rho, mu, nu, lambda) - DGamma_dUdd_ref(rho, mu, nu, lambda)) +
        Dg_dUU(rho, nu, lambda) * (Gamma_Udd(mu, nu, lambda) - Gamma_Udd_ref(mu, nu, lambda));
    output(rho, mu) = sum(nu, sum(lambda, expr));
    return output;
}

auto make_Xi_d(Func g_dd, Func Xi_U) -> Func {
    Var mu{"mu"};
    RDom nu{0, 4, "nu"};
    Func output;
    output(mu) = sum(g_dd(mu, nu) * Xi_U(nu));
    return output;
}

auto make_DXi_dd(Func g_dd, Func Dg_ddd, Func Xi_U, Func DXi_dU) -> Func {
    Var lambda{"lambda"}, mu{"mu"};
    RDom nu{0, 4, "nu"};
    Func output;
    output(lambda, mu) = sum(Dg_ddd(lambda, mu, nu) * Xi_U(nu) + g_dd(mu, nu) * DXi_dU(lambda, nu));
    return output;
}

auto covariant_derivative(Func v_d, Func Dv_dd, Func Gamma_Udd) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom lambda{0, 4, "lambda"};
    Func output;
    output(mu, nu) = Dv_dd(mu, nu) - sum(Gamma_Udd(lambda, nu, mu) * v_d(lambda));
    return output;
}

auto make_DivXi_dd(Func Xi_d, Func DXi_dd, Func Gamma_Udd) -> Func {
    Var mu{"mu"}, nu{"nu"};
    Func output;
    auto CovDXi_dd = covariant_derivative(Xi_d, DXi_dd, Gamma_Udd);
    output(mu, nu) = 0.5f * (CovDXi_dd(mu, nu) + CovDXi_dd(nu, mu));
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

auto make_R_Uddd(Func Gamma_Udd, Func DGamma_dUdd) -> Func {
    Var mu{"mu"}, nu{"nu"}, rho{"rho"}, sigma{"sigma"};
    RDom lambda{0, 4, "lambda"};
    Func output;
    Func temp;
    temp(nu, rho, sigma, mu) = sum(Gamma_Udd(lambda, nu, rho) * Gamma_Udd(sigma, mu, lambda));
    output(sigma, rho, mu, nu) = DGamma_dUdd(mu, sigma, nu, rho) - DGamma_dUdd(nu, sigma, mu, rho) +
                                 temp(nu, rho, sigma, mu) - temp(mu, rho, sigma, nu);
    return output;
}

auto make_R_dd(Func R_Uddd) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom rho{0, 4, "rho"};
    Func output;
    output(mu, nu) = sum(R_Uddd(rho, mu, rho, nu));
    return output;
}

auto make_F_dd(Expr z, Expr Psi, Func DPsi) -> Func {
    Var _mu{"mu"}, _nu{"nu"};
    Func output{"F_dd"};
    output(_mu, _nu) = cast<double>(0);
    FROM_EXPRESSIONS_F_dd(output);
    return output;
}

auto make_DF_ddd(Expr z, Expr Psi, Func DPsi, Func DDPsi) -> Func {
    Var mu{"mu"}, nu{"nu"}, rho{"rho"};
    Func output{"DF_ddd"};
    output(mu, nu, rho) = cast<double>(0);
    FROM_EXPRESSIONS_DF_ddd(output);
    return output;
}

auto make_F_Ud(Func F_dd, Func g_UU) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom rho{0, 4, "rho"};
    Func output{"F_Ud"};
    output(mu, nu) = sum(g_UU(mu, rho) * F_dd(rho, nu));
    return output;
}

auto make_F_UU(Func F_Ud, Func g_UU) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom rho{0, 4, "rho"};
    Func output{"F_UU"};
    output(mu, nu) = sum(F_Ud(mu, rho) * g_UU(rho, nu));
    return output;
}

auto make_DF_dUd(Func F_dd, Func DF_ddd, Func g_UU, Func Dg_dUU) -> Func {
    Var mu{"mu"}, nu{"nu"}, rho{"rho"};
    RDom lambda{0, 4, "lambda"};
    Func output{"DF_dUd"};
    output(rho, mu, nu) = sum(Dg_dUU(rho, mu, lambda) * F_dd(lambda, nu) +
                              g_UU(mu, lambda) * DF_ddd(rho, lambda, nu));
    return output;
}

auto make_DivF_d(Func F_Ud, Func DF_dUd, Func Gamma_Udd) -> Func {
    Var nu{"nu"};
    RDom lambda{0, 4, "lambda"}, rho{0, 4, "rho"};
    Func output;
    output(nu) = sum(rho, DF_dUd(rho, rho, nu) +
                              sum(lambda, Gamma_Udd(rho, rho, lambda) * F_Ud(lambda, nu)) -
                              sum(lambda, Gamma_Udd(lambda, rho, nu) * F_Ud(rho, lambda)));
    return output;
}

auto make_wave_equation(Func Df_d, Func DDf_dd, Func g_UU, Func Dg_dUU, Func Gamma_Udd, Expr V)
    -> Expr {
    RDom lambda{0, 4, "lambda"}, mu{0, 4, "mu"}, rho{0, 4, "rho"};
    return sum(rho, sum(mu, Dg_dUU(rho, rho, mu) * Df_d(mu) + g_UU(rho, mu) * DDf_dd(rho, mu) +
                                sum(lambda,
                                    Gamma_Udd(rho, rho, lambda) * g_UU(lambda, mu) * Df_d(mu)))) -
           V;
}

auto make_V(Expr z, Expr Phi, Expr Chi, Expr L) -> Expr {
    Expr FROM_EXPRESSIONS_V(V);
    return V;
}

auto make_DPhi_d(Expr z, Expr Phi, Func DPhi) -> Func {
    Var mu{"mu"};
    Func output{"DPhi_d"};
    output(mu) = cast<double>(0);
    FROM_EXPRESSIONS_DPhi_d(output);
    return output;
}

auto make_DDPhi_dd(Expr z, Expr Phi, Func DPhi, Func DDPhi) -> Func {
    Var mu{"mu"}, nu{"nu"};
    Func output{"DDPhi_dd"};
    output(mu, nu) = cast<double>(0);
    FROM_EXPRESSIONS_DDPhi_dd(output);
    return output;
}

auto make_DChi_d(Expr z, Expr Chi, Func DChi) -> Func {
    Var mu{"mu"};
    Func output{"DChi_d"};
    output(mu) = cast<double>(0);
    FROM_EXPRESSIONS_DChi_d(output);
    return output;
}

auto make_DDChi_dd(Expr z, Expr Chi, Func DChi, Func DDChi) -> Func {
    Var mu{"mu"}, nu{"nu"};
    Func output{"DDChi_dd"};
    output(mu, nu) = cast<double>(0);
    FROM_EXPRESSIONS_DDChi_dd(output);
    return output;
}

auto make_G_dd(Func R_dd, Func DPhi_d, Func DChi_d, Func F_dd, Func F_Ud, Func F_UU, Func g_dd,
               Expr V, Expr L) -> Func {
    Var mu{"mu"}, nu{"nu"};
    RDom lambda{0, 4, "lambda"}, rho{0, 4, "rho"};
    Func output{"G_dd"};
    output(mu, nu) =
        R_dd(mu, nu) + (3 / (L * L) + V) * g_dd(mu, nu) -
        (DPhi_d(mu) * DPhi_d(nu) + DChi_d(mu) * DChi_d(nu)) -
        (sum(lambda, -F_dd(mu, lambda) * F_Ud(lambda, nu)) -
         g_dd(mu, nu) / 4 * sum(lambda, sum(rho, F_dd(lambda, rho) * F_UU(lambda, rho))));
    return output;
}

// auto covariant_derivative(Func v_d, Func Dv_dd, Func Gamma_Udd) -> Func {
//     Var mu{"mu"}, nu{"nu"};
//     RDom lambda{0, 4, "lambda"};
//     Func output;
//     output(mu, nu) = Dv_dd(mu, nu) - sum(Gamma_Udd(lambda, nu, mu) * v_d(lambda));
//     return output;
// }

#if 0
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
#endif

template <class T> auto set_bounds(T& buffer) {
    for (auto i = 0; i < buffer.dimensions(); ++i) {
        buffer.dim(i).set_bounds(0, 4);
    }
}

class compute_all_generator : public Halide::Generator<compute_all_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DDQ", 3};
    Input<double> Psi{"Psi"};
    Input<Buffer<double>> DPsi{"DPsi", 1};
    Input<Buffer<double>> DDPsi{"DDPsi", 2};
    Input<double> Phi{"Phi"};
    Input<Buffer<double>> DPhi{"DPhi", 1};
    Input<Buffer<double>> DDPhi{"DDPhi", 2};
    Input<double> Chi{"Chi"};
    Input<Buffer<double>> DChi{"DChi", 1};
    Input<Buffer<double>> DDChi{"DDChi", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};

    Output<Buffer<double>> g_dd{"g_dd", 2};
    Output<Buffer<double>> g_UU{"g_UU", 2};
    Output<Buffer<double>> Dg_ddd{"Dg_ddd", 3};
    Output<Buffer<double>> Dg_dUU{"Dg_dUU", 3};
    Output<Buffer<double>> DDg_dddd{"DDg_dddd", 4};
    Output<Buffer<double>> Gamma_Udd{"Gamma_Udd", 3};
    Output<Buffer<double>> DGamma_dUdd{"DGamma_dUdd", 4};
    Output<Buffer<double>> Xi_U{"Xi_U", 1};
    Output<Buffer<double>> Xi_d{"Xi_d", 1};
    Output<Buffer<double>> DXi_dU{"DXi_dU", 2};
    Output<Buffer<double>> DXi_dd{"DXi_dd", 2};
    Output<Buffer<double>> DivXi_dd{"DivXi_dd", 2};
    Output<Buffer<double>> R_Uddd{"R_Uddd", 4};
    Output<Buffer<double>> R_dd{"R_dd", 2};
    Output<Buffer<double>> F_dd{"F_dd", 2};
    Output<Buffer<double>> DF_ddd{"DF_ddd", 3};
    Output<Buffer<double>> F_Ud{"F_Ud", 2};
    Output<Buffer<double>> DF_dUd{"DF_dUd", 3};
    Output<Buffer<double>> DivF_d{"DivF_d", 1};
    Output<Buffer<double>> EOM_Phi{"EOM_Phi", 1};
    Output<Buffer<double>> EOM_Chi{"EOM_Chi", 1};
    Output<Buffer<double>> G_dd{"G_dd", 2};

    void generate() {
        g_dd = make_g_dd(z, Q, length, chemical_potential);
        g_UU = make_g_UU(z, Q, length, chemical_potential);
        Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        Dg_dUU = make_Dg_dUU(z, Q, DQ, length, chemical_potential);
        DDg_dddd = make_DDg_dddd(z, Q, DQ, DDQ, length, chemical_potential);
        Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
        DGamma_dUdd = make_DGamma_dUdd(g_UU, Dg_dUU, Dg_ddd, DDg_dddd);

        auto Q_ref = make_Q_ref();
        auto DQ_ref = make_DQ_ref();
        auto DDQ_ref = make_DDQ_ref();
        auto g_UU_ref = make_g_UU(z, Q_ref, length, chemical_potential);
        auto Dg_ddd_ref = make_Dg_ddd(z, Q_ref, DQ_ref, length, chemical_potential);
        auto Dg_dUU_ref = make_Dg_dUU(z, Q_ref, DQ_ref, length, chemical_potential);
        auto DDg_dddd_ref = make_DDg_dddd(z, Q_ref, DQ_ref, DDQ_ref, length, chemical_potential);
        auto Gamma_Udd_ref = make_Gamma_Udd(g_UU_ref, Dg_ddd_ref);
        auto DGamma_dUdd_ref = make_DGamma_dUdd(g_UU_ref, Dg_dUU_ref, Dg_ddd_ref, DDg_dddd_ref);

        Xi_U = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
        Xi_d = make_Xi_d(g_dd, Xi_U);
        DXi_dU = make_DXi_dU(g_UU, Dg_dUU, Gamma_Udd, DGamma_dUdd, Gamma_Udd_ref, DGamma_dUdd_ref);
        DXi_dd = make_DXi_dd(g_dd, Dg_ddd, Xi_U, DXi_dU);
        DivXi_dd = make_DivXi_dd(Xi_d, DXi_dd, Gamma_Udd);
        R_Uddd = make_R_Uddd(Gamma_Udd, DGamma_dUdd);
        R_dd = make_R_dd(R_Uddd);
        F_dd = make_F_dd(z, Psi, DPsi);
        DF_ddd = make_DF_ddd(z, Psi, DPsi, DDPsi);
        F_Ud = make_F_Ud(F_dd, g_UU);
        DF_dUd = make_DF_dUd(F_dd, DF_ddd, g_UU, Dg_dUU);
        DivF_d = make_DivF_d(F_Ud, DF_dUd, Gamma_Udd);
        auto V = make_V(z, Phi, Chi, length);
        auto DPhi_d = make_DPhi_d(z, Phi, DPhi);
        auto DDPhi_dd = make_DDPhi_dd(z, Phi, DPhi, DDPhi);
        auto DChi_d = make_DChi_d(z, Chi, DChi);
        auto DDChi_dd = make_DDChi_dd(z, Chi, DChi, DDChi);
        {
            EOM_Phi(Var{"temp"}) = cast<double>(0);
            EOM_Phi(0) = make_wave_equation(DPhi_d, DDPhi_dd, g_UU, Dg_dUU, Gamma_Udd, V);
        }
        {
            EOM_Chi(Var{"temp"}) = cast<double>(0);
            EOM_Chi(0) = make_wave_equation(DChi_d, DDChi_dd, g_UU, Dg_dUU, Gamma_Udd, V);
        }

        auto F_UU = make_F_UU(F_Ud, g_UU);
        G_dd = make_G_dd(R_dd, DPhi_d, DChi_d, F_dd, F_Ud, F_UU, g_dd, V, length);

        set_bounds(g_dd);
        set_bounds(g_UU);
        set_bounds(Dg_ddd);
        set_bounds(Dg_dUU);
        set_bounds(DDg_dddd);
        set_bounds(Gamma_Udd);
        set_bounds(DGamma_dUdd);
        set_bounds(Xi_U);
        set_bounds(Xi_d);
        set_bounds(DXi_dU);
        set_bounds(DXi_dd);
        set_bounds(DivXi_dd);
        set_bounds(R_Uddd);
        set_bounds(R_dd);
        set_bounds(F_dd);
        set_bounds(DF_ddd);
        set_bounds(F_Ud);
        set_bounds(DF_dUd);
        set_bounds(DivF_d);
        set_bounds(G_dd);
    }
};

#if 0
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
class Dg_ddd_generator : public Halide::Generator<Dg_ddd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 3};

    void generate() {
        output = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4).dim(2).set_bounds(0, 4);
    }
};
class Dg_dUU_generator : public Halide::Generator<Dg_dUU_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 3};

    void generate() {
        output = make_Dg_dUU(z, Q, DQ, length, chemical_potential);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4).dim(2).set_bounds(0, 4);
    }
};
class DDg_dddd_generator : public Halide::Generator<DDg_dddd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DDQ", 3};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 4};

    void generate() {
        output = make_DDg_dddd(z, Q, DQ, DDQ, length, chemical_potential);
        output.dim(0)
            .set_bounds(0, 4)
            .dim(1)
            .set_bounds(0, 4)
            .dim(2)
            .set_bounds(0, 4)
            .dim(3)
            .set_bounds(0, 4);
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
class DGamma_dUdd_generator : public Halide::Generator<DGamma_dUdd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DDQ", 3};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 4};

    void generate() {
        auto g_UU = make_g_UU(z, Q, length, chemical_potential);
        auto Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        auto Dg_dUU = make_Dg_dUU(z, Q, DQ, length, chemical_potential);
        auto DDg_dddd = make_DDg_dddd(z, Q, DQ, DDQ, length, chemical_potential);
        auto Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
        output = make_DGamma_dUdd(g_UU, Dg_dUU, Dg_ddd, DDg_dddd);
        output.dim(0)
            .set_bounds(0, 4)
            .dim(1)
            .set_bounds(0, 4)
            .dim(2)
            .set_bounds(0, 4)
            .dim(3)
            .set_bounds(0, 4);
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
        auto g_UU_ref = make_g_UU(z, Q_ref, length, chemical_potential);
        auto Dg_ddd_ref = make_Dg_ddd(z, Q_ref, DQ_ref, length, chemical_potential);
        auto Gamma_Udd_ref = make_Gamma_Udd(g_UU_ref, Dg_ddd_ref);

        output = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
        output.dim(0).set_bounds(0, 4);
    }
};
class Xi_d_generator : public Halide::Generator<Xi_d_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 1};

    void generate() {
        auto g_dd = make_g_dd(z, Q, length, chemical_potential);
        auto g_UU = make_g_UU(z, Q, length, chemical_potential);
        auto Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential);
        auto Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);

        auto Q_ref = make_Q_ref();
        auto DQ_ref = make_DQ_ref();
        auto g_UU_ref = make_g_UU(z, Q_ref, length, chemical_potential);
        auto Dg_ddd_ref = make_Dg_ddd(z, Q_ref, DQ_ref, length, chemical_potential);
        auto Gamma_Udd_ref = make_Gamma_Udd(g_UU_ref, Dg_ddd_ref);

        auto Xi_U = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
        output = make_Xi_d(g_dd, Xi_U);
        output.dim(0).set_bounds(0, 4);
    }
};
class DXi_dd_generator : public Halide::Generator<DXi_dd_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DDQ", 3};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

    void generate() {
        auto g_dd = make_g_dd(z, Q, length, chemical_potential);
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

        auto Xi_U = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
        auto DXi_dU =
            make_DXi_dU(g_UU, Dg_dUU, Gamma_Udd, DGamma_dUdd, Gamma_Udd_ref, DGamma_dUdd_ref);
        output = make_DXi_dd(g_dd, Dg_ddd, Xi_U, DXi_dU);
        output.dim(0).set_bounds(0, 4).dim(1).set_bounds(0, 4);
    }
};
class DXi_dU_generator : public Halide::Generator<DXi_dU_generator> {
  public:
    Input<double> z{"z"};
    Input<Buffer<double>> Q{"Q", 1};
    Input<Buffer<double>> DQ{"DQ", 2};
    Input<Buffer<double>> DDQ{"DDQ", 3};
    Input<double> length{"length"};
    Input<double> chemical_potential{"chemical_potential"};
    Output<Buffer<double>> output{"output", 2};

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
HALIDE_REGISTER_GENERATOR(Dg_ddd_generator, Dg_ddd_generator)
HALIDE_REGISTER_GENERATOR(Dg_dUU_generator, Dg_dUU_generator)
HALIDE_REGISTER_GENERATOR(DDg_dddd_generator, DDg_dddd_generator)
HALIDE_REGISTER_GENERATOR(Gamma_Udd_generator, Gamma_Udd_generator)
HALIDE_REGISTER_GENERATOR(DGamma_dUdd_generator, DGamma_dUdd_generator)
HALIDE_REGISTER_GENERATOR(Xi_U_generator, Xi_U_generator)
HALIDE_REGISTER_GENERATOR(Xi_d_generator, Xi_d_generator)
HALIDE_REGISTER_GENERATOR(DXi_dU_generator, DXi_dU_generator)
HALIDE_REGISTER_GENERATOR(DXi_dd_generator, DXi_dd_generator)
#endif
HALIDE_REGISTER_GENERATOR(compute_all_generator, compute_all_generator)

// HALIDE_REGISTER_GENERATOR(g_dd_ref_generator, g_dd_ref_generator)
// HALIDE_REGISTER_GENERATOR(g_UU_ref_generator, g_UU_ref_generator)
// HALIDE_REGISTER_GENERATOR(Dg_dd_ref_generator, Dg_dd_ref_generator)
// HALIDE_REGISTER_GENERATOR(Dg_UU_ref_generator, Dg_UU_ref_generator)
// HALIDE_REGISTER_GENERATOR(DDg_dd_ref_generator, DDg_dd_ref_generator)
// HALIDE_REGISTER_GENERATOR(Gamma_Udd_ref_generator, Gamma_Udd_ref_generator)
// HALIDE_REGISTER_GENERATOR(DGamma_dUdd_ref_generator, DGamma_dUdd_ref_generator)
