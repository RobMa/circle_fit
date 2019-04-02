#pragma once

#include "circle_params.h"
#include "gwaf_taubin.h"
#include "geo_lm.h"

namespace circle_fit
{

// struct DatasetView{
//     const double x[];
//     const double y[];
//     const int length;
//     DatasetView(const double x[], const double y[], const int length) : x{x}, y{y}, length{length} {}
// };

/**
 * @brief Least squares estimate of a circle.
 * Uses the GWAF-Taubin method as initialization for the Levenberg-Marquardt solver.
 * @param x 
 * @param y 
 * @return CircleParams 
 */
inline CircleParams estimate_circle(const VecX& x, const VecX& y, bool *converged=nullptr)
{
    NormalizedDataset data(x,y);
    CircleParams init = gwaf_taubin::estimate_circle(data);
    CircleParams est = geometric_lm::estimate_circle(data, init, converged);
    circle_undo_normalization(data, est);
    return est;
}

} // circle_fit

