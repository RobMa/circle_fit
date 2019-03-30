#pragma once
#include <circle_fit/common.h>

namespace circle_fit{

namespace geometric_lm{
CircleParams estimate_circle(const Dataset& data, const CircleParams& init, bool *converged=nullptr, const double tolerance = 1e-13, const int max_iter=50, double mu = 3, const double sigma = 10);
}

}