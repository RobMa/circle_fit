#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <circle_fit/common.h>

namespace circle_fit{

namespace gwaf_taubin{
CircleParams estimate_circle(const Dataset& data);
// void estimate_circle_normalized(const VecX& x, const VecX& y, double &est_r, double& est_x, double &est_y);
}

}
