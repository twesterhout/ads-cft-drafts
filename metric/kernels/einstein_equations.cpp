#include "generated_expressions.hpp"
#include <Halide.h>
#include <cstdio>

using namespace Halide;

auto to_float_checked(double x) -> float {
  auto const float_x = static_cast<float>(x);
  if (static_cast<double>(float_x) != x) {
    throw std::runtime_error{"inexact conversion: " + std::to_string(x) +
                             " cannot be exactly represented by float"};
  }
  return float_x;
}

auto operator*(double a, Expr b) -> Expr { return to_float_checked(a) * b; }
auto operator*(Expr b, double a) -> Expr { return b * to_float_checked(a); }
auto operator/(double a, Expr b) -> Expr { return to_float_checked(a) / b; }
auto operator/(Expr b, double a) -> Expr { return b / to_float_checked(a); }
auto operator+(double a, Expr b) -> Expr { return to_float_checked(a) + b; }
auto operator-(double a, Expr b) -> Expr { return to_float_checked(a) - b; }

template <class... Batch>
auto make_g_dd(Func z, Func Q, Expr L, Expr mu, Batch... i) -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"g_dd"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_g_dd(output);
  return output;
}

template <class... Batch>
auto make_g_UU(Func z, Func Q, Expr L, Expr mu, Batch... i) -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"g_UU"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_g_UU(output);
  return output;
}

class metric_generator : public Halide::Generator<metric_generator> {
public:
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 4};
  Input<Buffer<double>> _DQ{"DQ", 5};
  Input<Buffer<double>> _DDQ{"DDQ", 6};
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};

  Output<Buffer<double>> _out_g_dd{"g_dd", 5};
  Output<Buffer<double>> _out_g_UU{"g_UU", 5};

  auto generate() -> void {
    Var i_x{"i_x"}, i_y{"i_y"}, i_z{"i_z"};
    Var mu{"mu"}, nu{"nu"};

    Func z;
    z(i_x, i_y, i_z) = _z(i_z);

    auto g_dd = make_g_dd(z, _Q, _length, _chemical_potential, i_x, i_y, i_z);
    auto g_UU = make_g_UU(z, _Q, _length, _chemical_potential, i_x, i_y, i_z);

    _out_g_dd(i_x, i_y, i_z, mu, nu) = g_dd(i_x, i_y, i_z, mu, nu);
    _out_g_UU(i_x, i_y, i_z, mu, nu) = g_UU(i_x, i_y, i_z, mu, nu);
  }
};

HALIDE_REGISTER_GENERATOR(metric_generator, metric_generator)

#if 0
template <class... Batch>
auto make_Dg_ddd(Func z, Func Q, Func DQ, Expr L, Expr mu, Batch... i) -> Func {
  Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Dg_ddd"};
  output(i..., _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Dg_ddd(output);
  return output;
}

template <class... Batch>
auto make_DDg_dddd(Func z, Func Q, Func DQ, Func DDQ, Expr L, Expr mu,
                   Batch... i) -> Func {
  Var _kappa{"kappa"}, _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"DDg_dddd"};
  output(i..., _kappa, _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDg_dddd(output);
  return output;
}

template <class... Batch>
auto make_Gamma_Udd_ref(Func z, Func Q, Func DQ, Expr L, Expr mu, Batch... i)
    -> Func {
  Var _rho{"rho"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Gamma_Udd_ref"};
  output(i..., _rho, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Gamma_Udd_ref(output);
  return output;
}

template <class... Batch>
auto make_DGamma_dUdd_ref(Func z, Func Q, Func DQ, Func DDQ, Expr L, Expr mu,
                          Batch... i) -> Func {
  Var _lambda{"lambda"}, _rho{"rho"}, _mu{"mu"}, _nu{"nu"};
  Func output{"DGamma_dUdd_ref"};
  output(i..., _lambda, _rho, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_DGamma_dUdd_ref(output);
  return output;
}

template <class... Batch>
auto make_Dg_dUU(Func z, Func Q, Func DQ, Expr L, Expr mu, Batch... i) -> Func {
  Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Dg_dUU"};
  output(i..., _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Dg_dUU(output);
  return output;
}

auto make_Gamma_Udd(Func g_UU, Func Dg_ddd) -> Func {
  Var mu{"mu"}, nu{"nu"}, lambda{"lambda"};
  RDom rho{0, 4, "rho"};

  Func Gamma_Udd;
  Gamma_Udd(_, lambda, mu, nu) =
      0.5f * sum(g_UU(_, lambda, rho) *
                 (Dg_ddd(_, mu, nu, rho) + Dg_ddd(_, nu, mu, rho) -
                  Dg_ddd(_, rho, mu, nu)));
  return Gamma_Udd;
}

auto make_DGamma_dUdd(Func g_UU, Func Dg_dUU, Func Dg_ddd, Func DDg_dddd)
    -> Func {
  Var mu{"mu"}, nu{"nu"}, lambda{"lambda"}, kappa{"kappa"};
  RDom rho{0, 4, "rho"};

  Func DGamma_dUdd;
  DGamma_dUdd(_, kappa, lambda, mu, nu) =
      0.5f * sum(Dg_dUU(_, kappa, lambda, rho) *
                     (Dg_ddd(_, mu, nu, rho) + Dg_ddd(_, nu, mu, rho) -
                      Dg_ddd(_, rho, mu, nu)) +
                 g_UU(_, lambda, rho) * (DDg_dddd(_, kappa, mu, nu, rho) +
                                         DDg_dddd(_, kappa, nu, mu, rho) -
                                         DDg_dddd(_, kappa, rho, mu, nu)));
  return DGamma_dUdd;
}

auto make_Xi_U(Func g_UU, Func Gamma_Udd, Func Gamma_Udd_ref) -> Func {
  Var mu{"mu"};
  RDom lambda{0, 4, "lambda"};
  RDom nu{0, 4, "nu"};
  Func Xi_U{"Xi_U"};
  auto expr = g_UU(_, nu, lambda) *
              (Gamma_Udd(_, mu, nu, lambda) - Gamma_Udd_ref(_, mu, nu, lambda));
  Xi_U(_, mu) = sum(nu, sum(lambda, expr));
  return Xi_U;
}

auto make_DXi_dU(Func g_UU, Func Dg_dUU, Func Gamma_Udd, Func DGamma_dUdd,
                 Func Gamma_Udd_ref, Func DGamma_dUdd_ref) -> Func {
  Var mu{"mu"}, rho{"rho"};
  RDom nu{0, 4, "nu"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  auto expr = g_UU(_, nu, lambda) * (DGamma_dUdd(_, rho, mu, nu, lambda) -
                                     DGamma_dUdd_ref(_, rho, mu, nu, lambda)) +
              Dg_dUU(_, rho, nu, lambda) * (Gamma_Udd(_, mu, nu, lambda) -
                                            Gamma_Udd_ref(_, mu, nu, lambda));
  output(_, rho, mu) = sum(nu, sum(lambda, expr));
  return output;
}

auto make_Xi_d(Func g_dd, Func Xi_U) -> Func {
  Var mu{"mu"};
  RDom nu{0, 4, "nu"};
  Func output;
  output(_, mu) = sum(g_dd(_, mu, nu) * Xi_U(_, nu));
  return output;
}

auto make_DXi_dd(Func g_dd, Func Dg_ddd, Func Xi_U, Func DXi_dU) -> Func {
  Var lambda{"lambda"}, mu{"mu"};
  RDom nu{0, 4, "nu"};
  Func output;
  output(_, lambda, mu) = sum(Dg_ddd(_, lambda, mu, nu) * Xi_U(_, nu) +
                              g_dd(_, mu, nu) * DXi_dU(_, lambda, nu));
  return output;
}

auto covariant_derivative(Func v_d, Func Dv_dd, Func Gamma_Udd) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  output(_, mu, nu) =
      Dv_dd(_, mu, nu) - sum(Gamma_Udd(_, lambda, nu, mu) * v_d(_, lambda));
  return output;
}

auto make_DivXi_dd(Func Xi_d, Func DXi_dd, Func Gamma_Udd) -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output;
  auto CovDXi_dd = covariant_derivative(Xi_d, DXi_dd, Gamma_Udd);
  output(_, mu, nu) = 0.5f * (CovDXi_dd(_, mu, nu) + CovDXi_dd(_, nu, mu));
  return output;
}

auto make_R_Uddd(Func Gamma_Udd, Func DGamma_dUdd) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"}, sigma{"sigma"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  Func temp;
  temp(_, nu, rho, sigma, mu) =
      sum(Gamma_Udd(_, lambda, nu, rho) * Gamma_Udd(_, sigma, mu, lambda));
  output(_, sigma, rho, mu, nu) =
      DGamma_dUdd(_, mu, sigma, nu, rho) - DGamma_dUdd(_, nu, sigma, mu, rho) +
      temp(_, nu, rho, sigma, mu) - temp(_, mu, rho, sigma, nu);
  return output;
}

auto make_R_dd(Func R_Uddd) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output;
  output(_, mu, nu) = sum(R_Uddd(_, rho, mu, rho, nu));
  return output;
}

template <class... Batch>
auto make_F_dd(Func z, Func Q, Func DQ, Batch... i) -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"F_dd"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_F_dd(output);
  return output;
}

template <class... Batch>
auto make_DF_ddd(Func z, Func Q, Func DQ, Func DDQ, Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"};
  Func output{"DF_ddd"};
  output(i..., mu, nu, rho) = cast<double>(0);
  FROM_EXPRESSIONS_DF_ddd(output);
  return output;
}

auto make_F_Ud(Func F_dd, Func g_UU) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output{"F_Ud"};
  output(_, mu, nu) = sum(g_UU(_, mu, rho) * F_dd(_, rho, nu));
  return output;
}

auto make_F_UU(Func F_Ud, Func g_UU) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output{"F_UU"};
  output(_, mu, nu) = sum(F_Ud(_, mu, rho) * g_UU(_, rho, nu));
  return output;
}

auto make_DF_dUd(Func F_dd, Func DF_ddd, Func g_UU, Func Dg_dUU) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"};
  RDom lambda{0, 4, "lambda"};
  Func output{"DF_dUd"};
  output(_, rho, mu, nu) =
      sum(Dg_dUU(_, rho, mu, lambda) * F_dd(_, lambda, nu) +
          g_UU(_, mu, lambda) * DF_ddd(_, rho, lambda, nu));
  return output;
}

auto make_DivF_d(Func F_Ud, Func DF_dUd, Func Gamma_Udd) -> Func {
  Var nu{"nu"};
  RDom lambda{0, 4, "lambda"}, mu{0, 4, "rho"};
  Func output;
  output(_, nu) = sum(
      mu, DF_dUd(_, mu, mu, nu) +
              sum(lambda, Gamma_Udd(_, mu, mu, lambda) * F_Ud(_, lambda, nu)) -
              sum(lambda, Gamma_Udd(_, lambda, mu, nu) * F_Ud(_, mu, lambda)));
  return output;
}

auto make_wave_equation(Func Df_d, Func DDf_dd, Func g_UU, Func Dg_dUU,
                        Func Gamma_Udd, Func DV_f) -> Func {
  RDom lambda{0, 4, "lambda"}, mu{0, 4, "mu"}, rho{0, 4, "rho"};
  Func output;
  output(_) =
      sum(rho,
          sum(mu, Dg_dUU(_, rho, rho, mu) * Df_d(_, mu) +
                      g_UU(_, rho, mu) * DDf_dd(_, rho, mu) +
                      sum(lambda, Gamma_Udd(_, rho, rho, lambda) *
                                      g_UU(_, lambda, mu) * Df_d(_, mu)))) -
      DV_f(_);
  return output;
}

template <class... Batch>
auto make_V(Func z, Func Q, Expr L, Batch... i) -> Func {
  Func output;
  FROM_EXPRESSIONS_V(output);
  return output;
}

template <class... Batch>
auto make_DV_Phi(Func z, Func Q, Expr L, Batch... i) -> Func {
  Func output;
  FROM_EXPRESSIONS_DV_Phi(output);
  return output;
}

template <class... Batch>
auto make_DV_Chi(Func z, Func Q, Expr L, Batch... i) -> Func {
  Func output;
  FROM_EXPRESSIONS_DV_Chi(output);
  return output;
}

template <class... Batch>
auto make_DPhi_d(Func z, Func Q, Func DQ, Batch... i) -> Func {
  Var mu{"mu"};
  Func output{"DPhi_d"};
  output(i..., mu) = cast<double>(0);
  FROM_EXPRESSIONS_DPhi_d(output);
  return output;
}

template <class... Batch>
auto make_DDPhi_dd(Func z, Func Q, Func DQ, Func DDQ, Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output{"DDPhi_dd"};
  output(i..., mu, nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDPhi_dd(output);
  return output;
}

template <class... Batch>
auto make_DChi_d(Func z, Func Q, Func DQ, Batch... i) -> Func {
  Var mu{"mu"};
  Func output{"DChi_d"};
  output(i..., mu) = cast<double>(0);
  FROM_EXPRESSIONS_DChi_d(output);
  return output;
}

template <class... Batch>
auto make_DDChi_dd(Func z, Func Q, Func DQ, Func DDQ, Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output{"DDChi_dd"};
  output(i..., mu, nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDChi_dd(output);
  return output;
}

auto make_G_dd(Func R_dd, Func DPhi_d, Func DChi_d, Func F_dd, Func F_Ud,
               Func F_UU, Func g_dd, Func V, Expr L) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom lambda{0, 4, "lambda"}, rho{0, 4, "rho"};
  Func output{"G_dd"};
  output(_, mu, nu) =
      R_dd(_, mu, nu) + (3 / (L * L) + V(_)) * g_dd(_, mu, nu) -
      (DPhi_d(_, mu) * DPhi_d(_, nu) + DChi_d(_, mu) * DChi_d(_, nu)) -
      (sum(lambda, -F_dd(_, mu, lambda) * F_Ud(_, lambda, nu)) -
       g_dd(_, mu, nu) / 4 *
           sum(lambda, sum(rho, F_dd(_, lambda, rho) * F_UU(_, lambda, rho))));
  return output;
}

template <class T> auto set_bounds(T &buffer) {
  // buffer.dim(0).set_bounds(0, 1);
  for (auto i = 1; i < buffer.dimensions(); ++i) {
    buffer.dim(i).set_bounds(0, 4);
  }
}

class compute_all_generator : public Halide::Generator<compute_all_generator> {
public:
  Input<Buffer<double>> z{"z", 1};
  Input<Buffer<double>> Q{"Q", 2};
  Input<Buffer<double>> DQ{"DQ", 3};
  Input<Buffer<double>> DDQ{"DDQ", 4};
  Input<double> length{"length"};
  Input<double> chemical_potential{"chemical_potential"};

  Output<Buffer<double>> g_dd{"g_dd", 3};
  Output<Buffer<double>> g_UU{"g_UU", 3};
  Output<Buffer<double>> Dg_ddd{"Dg_ddd", 4};
  Output<Buffer<double>> Dg_dUU{"Dg_dUU", 4};
  Output<Buffer<double>> DDg_dddd{"DDg_dddd", 5};
  Output<Buffer<double>> Gamma_Udd{"Gamma_Udd", 4};
  Output<Buffer<double>> DGamma_dUdd{"DGamma_dUdd", 5};
  Output<Buffer<double>> Xi_U{"Xi_U", 2};
  Output<Buffer<double>> Xi_d{"Xi_d", 2};
  Output<Buffer<double>> DXi_dU{"DXi_dU", 3};
  Output<Buffer<double>> DXi_dd{"DXi_dd", 3};
  Output<Buffer<double>> DivXi_dd{"DivXi_dd", 3};
  Output<Buffer<double>> R_Uddd{"R_Uddd", 5};
  Output<Buffer<double>> R_dd{"R_dd", 3};
  Output<Buffer<double>> F_dd{"F_dd", 3};
  Output<Buffer<double>> DF_ddd{"DF_ddd", 4};
  Output<Buffer<double>> F_Ud{"F_Ud", 3};
  Output<Buffer<double>> DF_dUd{"DF_dUd", 4};
  Output<Buffer<double>> DivF_d{"DivF_d", 2};
  Output<Buffer<double>> EOM_Phi{"EOM_Phi", 2};
  Output<Buffer<double>> EOM_Chi{"EOM_Chi", 2};
  Output<Buffer<double>> G_dd{"G_dd", 3};

  void generate() {
    Var i{"i"};

    g_dd = make_g_dd(z, Q, length, chemical_potential, i);
    g_UU = make_g_UU(z, Q, length, chemical_potential, i);
    Dg_ddd = make_Dg_ddd(z, Q, DQ, length, chemical_potential, i);
    Dg_dUU = make_Dg_dUU(z, Q, DQ, length, chemical_potential, i);
    DDg_dddd = make_DDg_dddd(z, Q, DQ, DDQ, length, chemical_potential, i);
    Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
    DGamma_dUdd = make_DGamma_dUdd(g_UU, Dg_dUU, Dg_ddd, DDg_dddd);
    auto Gamma_Udd_ref =
        make_Gamma_Udd_ref(z, Q, DQ, length, chemical_potential, i);
    auto DGamma_dUdd_ref =
        make_DGamma_dUdd_ref(z, Q, DQ, DDQ, length, chemical_potential, i);
    Xi_U = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
    Xi_d = make_Xi_d(g_dd, Xi_U);
    DXi_dU = make_DXi_dU(g_UU, Dg_dUU, Gamma_Udd, DGamma_dUdd, Gamma_Udd_ref,
                         DGamma_dUdd_ref);
    DXi_dd = make_DXi_dd(g_dd, Dg_ddd, Xi_U, DXi_dU);
    DivXi_dd = make_DivXi_dd(Xi_d, DXi_dd, Gamma_Udd);

    R_Uddd = make_R_Uddd(Gamma_Udd, DGamma_dUdd);
    R_dd = make_R_dd(R_Uddd);
    F_dd = make_F_dd(z, Q, DQ, i);
    DF_ddd = make_DF_ddd(z, Q, DQ, DDQ, i);
    F_Ud = make_F_Ud(F_dd, g_UU);
    DF_dUd = make_DF_dUd(F_dd, DF_ddd, g_UU, Dg_dUU);
    DivF_d = make_DivF_d(F_Ud, DF_dUd, Gamma_Udd);
    auto DPhi_d = make_DPhi_d(z, Q, DQ, i);
    auto DDPhi_dd = make_DDPhi_dd(z, Q, DQ, DDQ, i);
    auto DChi_d = make_DChi_d(z, Q, DQ, i);
    auto DDChi_dd = make_DDChi_dd(z, Q, DQ, DDQ, i);
    {
      EOM_Phi(i, Var{"temp"}) = cast<double>(0);
      EOM_Phi(i, 0) =
          make_wave_equation(DPhi_d, DDPhi_dd, g_UU, Dg_dUU, Gamma_Udd,
                             make_DV_Phi(z, Q, length, i))(i);
    }
    {
      EOM_Chi(i, Var{"temp"}) = cast<double>(0);
      EOM_Chi(i, 0) =
          make_wave_equation(DChi_d, DDChi_dd, g_UU, Dg_dUU, Gamma_Udd,
                             make_DV_Chi(z, Q, length, i))(i);
    }

    auto F_UU = make_F_UU(F_Ud, g_UU);
    auto V = make_V(z, Q, length, i);
    G_dd = make_G_dd(R_dd, DPhi_d, DChi_d, F_dd, F_Ud, F_UU, g_dd, V, length);

    g_dd.compute_root();
    g_UU.compute_root();
    Dg_ddd.compute_root();
    Dg_dUU.compute_root();
    DDg_dddd.compute_root();
    Gamma_Udd.compute_root();
    DGamma_dUdd.compute_root();
    Gamma_Udd_ref.compute_root();
    DGamma_dUdd_ref.compute_root();
    Xi_U.compute_root();
    Xi_d.compute_root();
    DXi_dU.compute_root();
    DXi_dd.compute_root();
    DivXi_dd.compute_root();
    R_Uddd.compute_root();
    R_dd.compute_root();
    F_dd.compute_root();
    DF_ddd.compute_root();
    F_Ud.compute_root();
    DF_dUd.compute_root();
    DivF_d.compute_root();
    DPhi_d.compute_root();
    DDPhi_dd.compute_root();
    DChi_d.compute_root();
    DDChi_dd.compute_root();

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

class einstein_equations_generator
    : public Halide::Generator<einstein_equations_generator> {
public:
  // Input<Buffer<double>> t{"t", 1};
  // Input<Buffer<double>> x{"x", 1};
  // Input<Buffer<double>> y{"y", 1};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 4};
  Input<Buffer<double>> _DQ{"DQ", 5};
  Input<Buffer<double>> _DDQ{"DDQ", 6};
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};
  Output<Buffer<double>> _output{"output", 4};

  void generate() {
    Var i_x{"i_x"}, i_y{"i_y"}, i_z{"i_z"}, j{"j"};

    Func z;
    z(i_x, i_y, i_z) = _z(i_z);

    auto g_dd = make_g_dd(z, _Q, _length, _chemical_potential, i_x, i_y, i_z);
    auto g_UU = make_g_UU(z, _Q, _length, _chemical_potential, i_x, i_y, i_z);
    auto Dg_ddd =
        make_Dg_ddd(z, _Q, _DQ, _length, _chemical_potential, i_x, i_y, i_z);
    auto Dg_dUU =
        make_Dg_dUU(z, _Q, _DQ, _length, _chemical_potential, i_x, i_y, i_z);
    auto DDg_dddd = make_DDg_dddd(z, _Q, _DQ, _DDQ, _length,
                                  _chemical_potential, i_x, i_y, i_z);
    auto Gamma_Udd_ref = make_Gamma_Udd_ref(z, _Q, _DQ, _length,
                                            _chemical_potential, i_x, i_y, i_z);
    auto DGamma_dUdd_ref = make_DGamma_dUdd_ref(
        z, _Q, _DQ, _DDQ, _length, _chemical_potential, i_x, i_y, i_z);
    auto Gamma_Udd = make_Gamma_Udd(g_UU, Dg_ddd);
    auto DGamma_dUdd = make_DGamma_dUdd(g_UU, Dg_dUU, Dg_ddd, DDg_dddd);
    auto Xi_U = make_Xi_U(g_UU, Gamma_Udd, Gamma_Udd_ref);
    auto Xi_d = make_Xi_d(g_dd, Xi_U);
    auto DXi_dU = make_DXi_dU(g_UU, Dg_dUU, Gamma_Udd, DGamma_dUdd,
                              Gamma_Udd_ref, DGamma_dUdd_ref);
    auto DXi_dd = make_DXi_dd(g_dd, Dg_ddd, Xi_U, DXi_dU);
    auto DivXi_dd = make_DivXi_dd(Xi_d, DXi_dd, Gamma_Udd);
    auto R_Uddd = make_R_Uddd(Gamma_Udd, DGamma_dUdd);
    auto R_dd = make_R_dd(R_Uddd);

    auto F_dd = make_F_dd(z, _Q, _DQ, i_x, i_y, i_z);
    auto DF_ddd = make_DF_ddd(z, _Q, _DQ, _DDQ, i_x, i_y, i_z);
    auto F_Ud = make_F_Ud(F_dd, g_UU);
    auto DF_dUd = make_DF_dUd(F_dd, DF_ddd, g_UU, Dg_dUU);
    auto DivF_d = make_DivF_d(F_Ud, DF_dUd, Gamma_Udd);
    auto DPhi_d = make_DPhi_d(z, _Q, _DQ, i_x, i_y, i_z);
    auto DDPhi_dd = make_DDPhi_dd(z, _Q, _DQ, _DDQ, i_x, i_y, i_z);
    auto DChi_d = make_DChi_d(z, _Q, _DQ, i_x, i_y, i_z);
    auto DDChi_dd = make_DDChi_dd(z, _Q, _DQ, _DDQ, i_x, i_y, i_z);
    auto F_UU = make_F_UU(F_Ud, g_UU);
    auto V = make_V(z, _Q, _length, i_x, i_y, i_z);
    auto G_dd =
        make_G_dd(R_dd, DPhi_d, DChi_d, F_dd, F_Ud, F_UU, g_dd, V, _length);

    _output(i_x, i_y, i_z, j) = cast<double>(0);
    _output(i_x, i_y, i_z, 0) = G_dd(i_x, i_y, i_z, 0, 0);
    _output(i_x, i_y, i_z, 1) = G_dd(i_x, i_y, i_z, 1, 1);
    _output(i_x, i_y, i_z, 2) = G_dd(i_x, i_y, i_z, 2, 2);
    _output(i_x, i_y, i_z, 3) = G_dd(i_x, i_y, i_z, 3, 3);
    _output(i_x, i_y, i_z, 4) = G_dd(i_x, i_y, i_z, 1, 3);
    _output(i_x, i_y, i_z, 5) = DivF_d(i_x, i_y, i_z, 0);
    _output(i_x, i_y, i_z, 6) = make_wave_equation(
        DPhi_d, DDPhi_dd, g_UU, Dg_dUU, Gamma_Udd,
        make_DV_Phi(z, _Q, _length, i_x, i_y, i_z))(i_x, i_y, i_z);
    _output(i_x, i_y, i_z, 7) = make_wave_equation(
        DChi_d, DDChi_dd, g_UU, Dg_dUU, Gamma_Udd,
        make_DV_Chi(z, _Q, _length, i_x, i_y, i_z))(i_x, i_y, i_z);
  }
};
#endif

#if 0
HALIDE_REGISTER_GENERATOR(compute_all_generator, compute_all_generator)
HALIDE_REGISTER_GENERATOR(einstein_equations_generator,
                          einstein_equations_generator)
#endif
