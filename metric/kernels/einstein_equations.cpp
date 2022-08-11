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
auto make_g_dd(bool auto_schedule, Func z, Func Q, Expr L, Expr mu, Batch... i)
    -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"g_dd"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_g_dd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_g_UU(bool auto_schedule, Func z, Func Q, Expr L, Expr mu, Batch... i)
    -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"g_UU"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_g_UU(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_Dg_ddd(bool auto_schedule, Func z, Func Q, Func DQ, Expr L, Expr mu,
                 Batch... i) -> Func {
  Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Dg_ddd"};
  output(i..., _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Dg_ddd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_Dg_dUU(bool auto_schedule, Func z, Func Q, Func DQ, Expr L, Expr mu,
                 Batch... i) -> Func {
  Var _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Dg_dUU"};
  output(i..., _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Dg_dUU(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DDg_dddd(bool auto_schedule, Func z, Func Q, Func DQ, Func DDQ,
                   Expr L, Expr mu, Batch... i) -> Func {
  Var _kappa{"kappa"}, _lambda{"lambda"}, _mu{"mu"}, _nu{"nu"};
  Func output{"DDg_dddd"};
  output(i..., _kappa, _lambda, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDg_dddd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class metric_generator : public Halide::Generator<metric_generator> {
public:
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};

  Output<Buffer<double>> _out_g_dd{"g_dd", 3};
  Output<Buffer<double>> _out_g_UU{"g_UU", 3};
  Output<Buffer<double>> _out_Dg_ddd{"Dg_ddd", 4};
  Output<Buffer<double>> _out_Dg_dUU{"Dg_dUU", 4};
  Output<Buffer<double>> _out_DDg_dddd{"DDg_dddd", 5};

  auto generate() -> void {
    Var i{"i"};

    _out_g_dd =
        make_g_dd(auto_schedule, _z, _Q, _length, _chemical_potential, i);
    _out_g_UU =
        make_g_UU(auto_schedule, _z, _Q, _length, _chemical_potential, i);
    _out_Dg_ddd = make_Dg_ddd(auto_schedule, _z, _Q, _DQ, _length,
                              _chemical_potential, i);
    _out_Dg_dUU = make_Dg_dUU(auto_schedule, _z, _Q, _DQ, _length,
                              _chemical_potential, i);
    _out_DDg_dddd = make_DDg_dddd(auto_schedule, _z, _Q, _DQ, _DDQ, _length,
                                  _chemical_potential, i);
  }
};

HALIDE_REGISTER_GENERATOR(metric_generator, metric_generator)

template <class... Batch>
auto make_F_dd(bool auto_schedule, Func z, Func Q, Func DQ, Batch... i)
    -> Func {
  Var _mu{"mu"}, _nu{"nu"};
  Func output{"F_dd"};
  output(i..., _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_F_dd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DF_ddd(bool auto_schedule, Func z, Func Q, Func DQ, Func DDQ,
                 Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"};
  Func output{"DF_ddd"};
  output(i..., mu, nu, rho) = cast<double>(0);
  FROM_EXPRESSIONS_DF_ddd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_F_Ud(bool auto_schedule, Func F_dd, Func g_UU) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output{"F_Ud"};
  output(_, mu, nu) = sum(g_UU(_, mu, rho) * F_dd(_, rho, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_F_UU(bool auto_schedule, Func F_Ud, Func g_UU) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output{"F_UU"};
  output(_, mu, nu) = sum(F_Ud(_, mu, rho) * g_UU(_, rho, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_DF_dUd(bool auto_schedule, Func F_dd, Func DF_ddd, Func g_UU,
                 Func Dg_dUU) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"};
  RDom lambda{0, 4, "lambda"};
  Func output{"DF_dUd"};
  output(_, rho, mu, nu) =
      sum(Dg_dUU(_, rho, mu, lambda) * F_dd(_, lambda, nu) +
          g_UU(_, mu, lambda) * DF_ddd(_, rho, lambda, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_DivF_d(bool auto_schedule, Func F_Ud, Func DF_dUd, Func Gamma_Udd)
    -> Func {
  Var nu{"nu"};
  RDom lambda{0, 4, "lambda"}, mu{0, 4, "rho"};
  Func output;
  output(_, nu) = sum(
      mu, DF_dUd(_, mu, mu, nu) +
              sum(lambda, Gamma_Udd(_, mu, mu, lambda) * F_Ud(_, lambda, nu)) -
              sum(lambda, Gamma_Udd(_, lambda, mu, nu) * F_Ud(_, mu, lambda)));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class maxwell_generator : public Halide::Generator<maxwell_generator> {
public:
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};
  Input<Buffer<double>> _g_UU{"g_UU", 3};
  Input<Buffer<double>> _Dg_ddd{"Dg_ddd", 4};
  Input<Buffer<double>> _Dg_dUU{"Dg_dUU", 4};
  Input<Buffer<double>> _Gamma_Udd{"Gamma_Udd", 4};

  Output<Buffer<double>> _out_F_dd{"F_dd", 3};
  Output<Buffer<double>> _out_F_Ud{"F_Ud", 3};
  Output<Buffer<double>> _out_F_UU{"F_UU", 3};
  Output<Buffer<double>> _out_DivF_d{"DivF_d", 2};

  auto generate() -> void {
    Var i{"i"};

    _out_F_dd = make_F_dd(auto_schedule, _z, _Q, _DQ, i);
    _out_F_Ud = make_F_Ud(auto_schedule, _out_F_dd, _g_UU);
    _out_F_UU = make_F_UU(auto_schedule, _out_F_Ud, _g_UU);
    auto DF_ddd = make_DF_ddd(auto_schedule, _z, _Q, _DQ, _DDQ, i);
    auto DF_dUd = make_DF_dUd(auto_schedule, _out_F_dd, DF_ddd, _g_UU, _Dg_dUU);
    _out_DivF_d = make_DivF_d(auto_schedule, _out_F_Ud, DF_dUd, _Gamma_Udd);
  }
};

HALIDE_REGISTER_GENERATOR(maxwell_generator, maxwell_generator)

auto make_Gamma_Udd(bool auto_schedule, Func g_UU, Func Dg_ddd) -> Func {
  Var mu{"mu"}, nu{"nu"}, lambda{"lambda"};
  RDom rho{0, 4, "rho"};

  Func Gamma_Udd;
  Gamma_Udd(_, lambda, mu, nu) =
      0.5f * sum(g_UU(_, lambda, rho) *
                 (Dg_ddd(_, mu, nu, rho) + Dg_ddd(_, nu, mu, rho) -
                  Dg_ddd(_, rho, mu, nu)));
  if (!auto_schedule) {
    Gamma_Udd.compute_root();
  }
  return Gamma_Udd;
}

auto make_DGamma_dUdd(bool auto_schedule, Func g_UU, Func Dg_dUU, Func Dg_ddd,
                      Func DDg_dddd) -> Func {
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
  if (!auto_schedule) {
    DGamma_dUdd.compute_root();
  }
  return DGamma_dUdd;
}

class christoffel_generator : public Halide::Generator<christoffel_generator> {
public:
  Input<Buffer<double>> _g_UU{"g_UU", 3};
  Input<Buffer<double>> _Dg_ddd{"Dg_ddd", 4};
  Input<Buffer<double>> _Dg_dUU{"Dg_dUU", 4};
  Input<Buffer<double>> _DDg_dddd{"DDg_dddd", 5};

  Output<Buffer<double>> _out_Gamma_Udd{"Gamma_Udd", 4};
  Output<Buffer<double>> _out_DGamma_dUdd{"DGamma_dUdd", 5};

  auto generate() -> void {
    _out_Gamma_Udd = make_Gamma_Udd(auto_schedule, _g_UU, _Dg_ddd);
    _out_DGamma_dUdd =
        make_DGamma_dUdd(auto_schedule, _g_UU, _Dg_dUU, _Dg_ddd, _DDg_dddd);
  }
};

HALIDE_REGISTER_GENERATOR(christoffel_generator, christoffel_generator)

auto make_Xi_U(bool auto_schedule, Func g_UU, Func Gamma_Udd,
               Func Gamma_Udd_ref) -> Func {
  Var mu{"mu"};
  RDom lambda{0, 4, "lambda"};
  RDom nu{0, 4, "nu"};
  Func Xi_U{"Xi_U"};
  auto expr = g_UU(_, nu, lambda) *
              (Gamma_Udd(_, mu, nu, lambda) - Gamma_Udd_ref(_, mu, nu, lambda));
  Xi_U(_, mu) = sum(nu, sum(lambda, expr));
  if (!auto_schedule) {
    Xi_U.compute_root();
  }
  return Xi_U;
}

auto make_Xi_d(bool auto_schedule, Func g_dd, Func Xi_U) -> Func {
  Var mu{"mu"};
  RDom nu{0, 4, "nu"};
  Func output;
  output(_, mu) = sum(g_dd(_, mu, nu) * Xi_U(_, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_DXi_dU(bool auto_schedule, Func g_UU, Func Dg_dUU, Func Gamma_Udd,
                 Func DGamma_dUdd, Func Gamma_Udd_ref, Func DGamma_dUdd_ref)
    -> Func {
  Var mu{"mu"}, rho{"rho"};
  RDom nu{0, 4, "nu"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  auto expr = g_UU(_, nu, lambda) * (DGamma_dUdd(_, rho, mu, nu, lambda) -
                                     DGamma_dUdd_ref(_, rho, mu, nu, lambda)) +
              Dg_dUU(_, rho, nu, lambda) * (Gamma_Udd(_, mu, nu, lambda) -
                                            Gamma_Udd_ref(_, mu, nu, lambda));
  output(_, rho, mu) = sum(nu, sum(lambda, expr));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_DXi_dd(bool auto_schedule, Func g_dd, Func Dg_ddd, Func Xi_U,
                 Func DXi_dU) -> Func {
  Var lambda{"lambda"}, mu{"mu"};
  RDom nu{0, 4, "nu"};
  Func output;
  output(_, lambda, mu) = sum(Dg_ddd(_, lambda, mu, nu) * Xi_U(_, nu) +
                              g_dd(_, mu, nu) * DXi_dU(_, lambda, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto covariant_derivative(bool auto_schedule, Func v_d, Func Dv_dd,
                          Func Gamma_Udd) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  output(_, mu, nu) =
      Dv_dd(_, mu, nu) - sum(Gamma_Udd(_, lambda, nu, mu) * v_d(_, lambda));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_DivXi_dd(bool auto_schedule, Func Xi_d, Func DXi_dd, Func Gamma_Udd)
    -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output;
  auto CovDXi_dd = covariant_derivative(auto_schedule, Xi_d, DXi_dd, Gamma_Udd);
  output(_, mu, nu) = 0.5f * (CovDXi_dd(_, mu, nu) + CovDXi_dd(_, nu, mu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_Gamma_Udd_ref(bool auto_schedule, Func z, Func Q, Func DQ, Expr L,
                        Expr mu, Batch... i) -> Func {
  Var _rho{"rho"}, _mu{"mu"}, _nu{"nu"};
  Func output{"Gamma_Udd_ref"};
  output(i..., _rho, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_Gamma_Udd_ref(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DGamma_dUdd_ref(bool auto_schedule, Func z, Func Q, Func DQ, Func DDQ,
                          Expr L, Expr mu, Batch... i) -> Func {
  Var _lambda{"lambda"}, _rho{"rho"}, _mu{"mu"}, _nu{"nu"};
  Func output{"DGamma_dUdd_ref"};
  output(i..., _lambda, _rho, _mu, _nu) = cast<double>(0);
  FROM_EXPRESSIONS_DGamma_dUdd_ref(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class deturck_generator : public Halide::Generator<deturck_generator> {
public:
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};
  Input<Buffer<double>> _g_dd{"g_dd", 3};
  Input<Buffer<double>> _g_UU{"g_UU", 3};
  Input<Buffer<double>> _Dg_ddd{"Dg_ddd", 4};
  Input<Buffer<double>> _Dg_dUU{"Dg_dUU", 4};
  Input<Buffer<double>> _Gamma_Udd{"Gamma_Udd", 4};
  Input<Buffer<double>> _DGamma_dUdd{"DGamma_dUdd", 5};

  Output<Buffer<double>> _out_divXi_dd{"divXi_dd", 3};

  void generate() {
    Var i{"i"};

    auto Gamma_Udd_ref = make_Gamma_Udd_ref(auto_schedule, _z, _Q, _DQ, _length,
                                            _chemical_potential, i);
    auto DGamma_dUdd_ref = make_DGamma_dUdd_ref(
        auto_schedule, _z, _Q, _DQ, _DDQ, _length, _chemical_potential, i);
    auto Xi_U = make_Xi_U(auto_schedule, _g_UU, _Gamma_Udd, Gamma_Udd_ref);
    auto Xi_d = make_Xi_d(auto_schedule, _g_dd, Xi_U);
    auto DXi_dU = make_DXi_dU(auto_schedule, _g_UU, _Dg_dUU, _Gamma_Udd,
                              _DGamma_dUdd, Gamma_Udd_ref, DGamma_dUdd_ref);
    auto DXi_dd = make_DXi_dd(auto_schedule, _g_dd, _Dg_ddd, Xi_U, DXi_dU);

    _out_divXi_dd = make_DivXi_dd(auto_schedule, Xi_d, DXi_dd, _Gamma_Udd);
  }
};

HALIDE_REGISTER_GENERATOR(deturck_generator, deturck_generator)

auto make_R_Uddd(bool auto_schedule, Func Gamma_Udd, Func DGamma_dUdd) -> Func {
  Var mu{"mu"}, nu{"nu"}, rho{"rho"}, sigma{"sigma"};
  RDom lambda{0, 4, "lambda"};
  Func output;
  Func temp;
  temp(_, nu, rho, sigma, mu) =
      sum(Gamma_Udd(_, lambda, nu, rho) * Gamma_Udd(_, sigma, mu, lambda));
  output(_, sigma, rho, mu, nu) =
      DGamma_dUdd(_, mu, sigma, nu, rho) - DGamma_dUdd(_, nu, sigma, mu, rho) +
      temp(_, nu, rho, sigma, mu) - temp(_, mu, rho, sigma, nu);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

auto make_R_dd(bool auto_schedule, Func R_Uddd) -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom rho{0, 4, "rho"};
  Func output;
  output(_, mu, nu) = sum(R_Uddd(_, rho, mu, rho, nu));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DPhi_d(bool auto_schedule, Func z, Func Q, Func DQ, Batch... i)
    -> Func {
  Var mu{"mu"};
  Func output{"DPhi_d"};
  output(i..., mu) = cast<double>(0);
  FROM_EXPRESSIONS_DPhi_d(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DDPhi_dd(bool auto_schedule, Func z, Func Q, Func DQ, Func DDQ,
                   Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output{"DDPhi_dd"};
  output(i..., mu, nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDPhi_dd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DChi_d(bool auto_schedule, Func z, Func Q, Func DQ, Batch... i)
    -> Func {
  Var mu{"mu"};
  Func output{"DChi_d"};
  output(i..., mu) = cast<double>(0);
  FROM_EXPRESSIONS_DChi_d(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DDChi_dd(bool auto_schedule, Func z, Func Q, Func DQ, Func DDQ,
                   Batch... i) -> Func {
  Var mu{"mu"}, nu{"nu"};
  Func output{"DDChi_dd"};
  output(i..., mu, nu) = cast<double>(0);
  FROM_EXPRESSIONS_DDChi_dd(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_V(bool auto_schedule, Func z, Func Q, Expr L, Batch... i) -> Func {
  Func output;
  FROM_EXPRESSIONS_V(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DV_Phi(bool auto_schedule, Func z, Func Q, Expr L, Batch... i)
    -> Func {
  Func output;
  FROM_EXPRESSIONS_DV_Phi(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

template <class... Batch>
auto make_DV_Chi(bool auto_schedule, Func z, Func Q, Expr L, Batch... i)
    -> Func {
  Func output;
  FROM_EXPRESSIONS_DV_Chi(output);
  if (!auto_schedule) {
    output.compute_root();
  }
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

auto make_G_dd(bool auto_schedule, Func R_dd, Func DPhi_d, Func DChi_d,
               Func F_dd, Func F_Ud, Func F_UU, Func g_dd, Func V, Expr L)
    -> Func {
  Var mu{"mu"}, nu{"nu"};
  RDom lambda{0, 4, "lambda"}, rho{0, 4, "rho"};
  Func output{"G_dd"};
  output(_, mu, nu) =
      R_dd(_, mu, nu) + (3 / (L * L) - V(_)) * g_dd(_, mu, nu) -
      2 * (DPhi_d(_, mu) * DPhi_d(_, nu) + DChi_d(_, mu) * DChi_d(_, nu)) -
      (sum(lambda, -F_dd(_, mu, lambda) * F_Ud(_, lambda, nu)) -
       g_dd(_, mu, nu) / 4 *
           sum(lambda, sum(rho, F_dd(_, lambda, rho) * F_UU(_, lambda, rho))));
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class equations_generator : public Halide::Generator<equations_generator> {
public:
  Input<double> _length{"length"};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _z{"z", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};
  Input<Buffer<double>> _g_dd{"g_dd", 3};
  Input<Buffer<double>> _g_UU{"g_UU", 3};
  Input<Buffer<double>> _Dg_dUU{"Dg_dUU", 4};
  Input<Buffer<double>> _Gamma_Udd{"Gamma_Udd", 4};
  Input<Buffer<double>> _DGamma_dUdd{"DGamma_dUdd", 5};
  Input<Buffer<double>> _F_dd{"F_dd", 3};
  Input<Buffer<double>> _F_Ud{"F_Ud", 3};
  Input<Buffer<double>> _F_UU{"F_UU", 3};
  Input<Buffer<double>> _divF_d{"divF_d", 2};
  Input<Buffer<double>> _divXi_dd{"divXi_dd", 3};

  Output<Buffer<double>> _out_equations{"equations", 2};

  void generate() {
    Var i{"i"}, j{"j"};

    auto R_Uddd = make_R_Uddd(auto_schedule, _Gamma_Udd, _DGamma_dUdd);
    auto R_dd = make_R_dd(auto_schedule, R_Uddd);

    auto DPhi_d = make_DPhi_d(auto_schedule, _z, _Q, _DQ, i);
    auto DDPhi_dd = make_DDPhi_dd(auto_schedule, _z, _Q, _DQ, _DDQ, i);
    auto DChi_d = make_DChi_d(auto_schedule, _z, _Q, _DQ, i);
    auto DDChi_dd = make_DDChi_dd(auto_schedule, _z, _Q, _DQ, _DDQ, i);
    auto V = make_V(auto_schedule, _z, _Q, _length, i);

    auto G_dd = make_G_dd(auto_schedule, R_dd, DPhi_d, DChi_d, _F_dd, _F_Ud,
                          _F_UU, _g_dd, V, _length);

    _out_equations(i, j) = cast<double>(0);
    _out_equations(i, 0) = G_dd(i, 0, 0) - _divXi_dd(i, 0, 0);
    _out_equations(i, 1) = G_dd(i, 1, 1) - _divXi_dd(i, 1, 1);
    _out_equations(i, 2) = G_dd(i, 2, 2) - _divXi_dd(i, 2, 2);
    _out_equations(i, 3) = G_dd(i, 3, 3) - _divXi_dd(i, 3, 3);
    _out_equations(i, 4) = G_dd(i, 1, 3) - _divXi_dd(i, 1, 3);
    _out_equations(i, 5) = _divF_d(i, 0);
    _out_equations(i, 6) =
        make_wave_equation(DPhi_d, DDPhi_dd, _g_UU, _Dg_dUU, _Gamma_Udd,
                           make_DV_Phi(auto_schedule, _z, _Q, _length, i))(i);
    _out_equations(i, 7) =
        make_wave_equation(DChi_d, DDChi_dd, _g_UU, _Dg_dUU, _Gamma_Udd,
                           make_DV_Chi(auto_schedule, _z, _Q, _length, i))(i);
  }
};

HALIDE_REGISTER_GENERATOR(equations_generator, equations_generator)

template <class... Batch>
auto make_horizon(bool auto_schedule, Func Q, Func DQ, Func DDQ, Expr L,
                  Expr mu, Batch... i) -> Func {
  Var j{"j"};
  Func output{"horizon"};
  output(i..., j) = cast<double>(0);
  FROM_EXPRESSIONS_Horizon(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class horizon_boundary_generator
    : public Halide::Generator<horizon_boundary_generator> {
public:
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};

  Output<Buffer<double>> _out_horizon{"horizon", 2};

  void generate() {
    Var i{"i"};

    _out_horizon = make_horizon(auto_schedule, _Q, _DQ, _DDQ, _length,
                                _chemical_potential, i);
  }
};

HALIDE_REGISTER_GENERATOR(horizon_boundary_generator,
                          horizon_boundary_generator)

template <class... Batch>
auto make_conformal(bool auto_schedule, Func x, Func y, Func Q, Func DQ,
                    Func DDQ, Expr L, Expr mu, Expr V0, Expr k0, Expr theta,
                    Batch... i) -> Func {
  Var j{"j"};
  Func output{"conformal"};
  output(i..., j) = cast<double>(0);
  FROM_EXPRESSIONS_Conformal(output);
  if (!auto_schedule) {
    output.compute_root();
  }
  return output;
}

class conformal_boundary_generator
    : public Halide::Generator<conformal_boundary_generator> {
public:
  Input<double> _length{"length"};
  Input<double> _chemical_potential{"chemical_potential"};
  Input<double> _v0{"V0"};
  Input<double> _k0{"k0"};
  Input<double> _theta{"theta"};
  Input<Buffer<double>> _x{"x", 1};
  Input<Buffer<double>> _y{"y", 1};
  Input<Buffer<double>> _Q{"Q", 2};
  Input<Buffer<double>> _DQ{"DQ", 3};
  Input<Buffer<double>> _DDQ{"DDQ", 4};

  Output<Buffer<double>> _out_conformal{"conformal", 2};

  void generate() {
    Var i{"i"};

    _out_conformal =
        make_conformal(auto_schedule, _x, _y, _Q, _DQ, _DDQ, _length,
                       _chemical_potential, _v0, _k0, _theta, i);
  }
};

HALIDE_REGISTER_GENERATOR(conformal_boundary_generator,
                          conformal_boundary_generator)

#if 0
class foo_generator : public Halide::Generator<foo_generator> {
public:
  Output<Buffer<double>> _out{"out", 2};

  auto generate() -> void {
    Var i{"i"}, j{"j"};

    _out(i, j) = cast<double>(i);
  }
};

HALIDE_REGISTER_GENERATOR(foo_generator, foo_generator)
#endif

#if 0

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
