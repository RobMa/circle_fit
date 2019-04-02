#pragma once

#include <Eigen/Core>

namespace circle_fit{

using Real = double;
using Vec2 = Eigen::Matrix<Real, 2, 1>;
using Vec3 = Eigen::Matrix<Real, 3, 1>;
using Vec4 = Eigen::Matrix<Real, 4, 1>;
using Mat3 = Eigen::Matrix<Real, 3, 3>;
using Mat4 = Eigen::Matrix<Real, 4, 4>;
using Mat6 = Eigen::Matrix<Real, 6, 6>;
using VecX = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

}




