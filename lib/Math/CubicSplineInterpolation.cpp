#include "CubicSplineInterpolation.h"

using namespace Math;

namespace
{
  void validate_input(const std::vector<double> &x, const std::vector<double> &y)
  {
    if (x.size() != y.size())
      THROW_INVALID_ARGUMENT(
        "knot_points_x = " + std::to_string(x.size()) +
        " knot_points_y.size() = " + std::to_string(y.size())
      );

    if (y.size() < 2)
      THROW_INVALID_ARGUMENT("knot_points_y_.size() = " + std::to_string(y.size()));

    if (!std::is_sorted(x.begin(), x.end()))
      THROW_INVALID_ARGUMENT("knot_points_x is not sorted");
  }

  Eigen::VectorXd get_derivatives(const std::vector<double> &x, const std::vector<double> &y)
  {
    validate_input(x, y);

    /// algorithm based on wikipedia article:
    /// https://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline

    const auto n = (int) x.size();
    Eigen::MatrixXd coeff(n, n);
    coeff.setZero();
    Eigen::VectorXd d(n);
    d.setZero();
    // left boundary condition
    {
      const auto step = (x[1] - x[0]);
      coeff(0, 0) = 2.0 / step;
      coeff(0, 1) = 1.0 / step;

      d(0) = 3.0 * (y[1] - y[0]) / step / step;
    }

    // middle part
    for (auto i = 1; i < n - 1; ++i) {
      const auto step_back = x[i] - x[i - 1];
      const auto step_forward = x[i + 1] - x[i];
      coeff(i, i - 1) = 1.0 / step_back;
      coeff(i, i) = 2.0 * (1.0 / step_back + step_forward);
      coeff(i, i + 1) = 1.0 / step_forward;

      d(i) = 3.0 * ((y[i] - y[i - 1]) / step_back / step_back + (y[i + 1] - y[i]) / step_forward / step_forward);
    }

    // right boundary condition
    {
      const auto step = (x[n - 1] - x[n - 2]);
      coeff(n - 1, n - 2) = 1.0 / step;
      coeff(n - 1, n - 1) = 2.0 / step;

      d(n - 1) = 3.0 * (y[n - 1] - y[n - 2]) / step / step;
    }

    // TODO replace to use tri-diagonal linear equation solver
    Eigen::VectorXd res = coeff.colPivHouseholderQr().solve(d).cast<double>();

    return res;
  }
}

CubicSplineInterpolation::CubicSplineInterpolation(
  const std::vector<double> &knot_points_x_, const std::vector<double> &knot_points_y_
) : knot_points_x(knot_points_x_), knot_points_y(knot_points_y_),
    k(get_derivatives(knot_points_x_, knot_points_y_)) {}

double CubicSplineInterpolation::interpolate(double x) const
{
  if (x < knot_points_x.front())
    THROW_INVALID_ARGUMENT(
      "x = " + std::to_string(x) + " < knot_points_x.front() = " + std::to_string(knot_points_x.front())
    );
  if (x > knot_points_x.back())
    THROW_INVALID_ARGUMENT(
      "x = " + std::to_string(x) + " > knot_points_x.back() = " + std::to_string(knot_points_x.back())
    );

  const auto lower = std::upper_bound(knot_points_x.begin(), knot_points_x.end(), x);
  const auto &x1 = *(lower - 1);
  const auto &x2 = *lower;
  const auto distance = std::distance(knot_points_x.begin(), lower);
  const auto &y1 = knot_points_y[distance - 1];
  const auto &y2 = knot_points_y[distance];

  const auto &k1 = k[distance - 1];
  const auto &k2 = k[distance];

  const auto t = (x - x1) / (x2 - x1);
  const auto a = k1 * (x2 - x1) - (y2 - y1);
  const auto b = -k2 * (x2 - x1) + (y2 - y1);

  const auto res = (1.0 - t) * y1 + t * y2 + t * (1.0 - t) * (a * (1.0 - t) + b * t);

  return res;
}
