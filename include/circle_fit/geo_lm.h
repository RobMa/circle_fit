#pragma once
#include <circle_fit/dataset.h>

namespace circle_fit{

namespace geometric_lm{
// CircleParams estimate_circle(const Dataset& data, const CircleParams& init, bool *converged=nullptr, const double tolerance = 1e-15, const int max_iter=50, double mu = 0.0001, const double sigma = 10);
CircleParams estimate_circle(const Dataset& data, const CircleParams& init, const int max_iter=50, double mu=1e-4, const double sigma_inc=10, const double sigma_dec=0.04);

struct IterationData
{
    double mse;
    double gradient[3];
    double hessian[9];
    double rxy[3];
    double gradient_norm;
    double hessian_norm;

    IterationData() = default;
    IterationData(double mse, const Vec3& gradient, const Mat3 hessian, const CircleParams& c) :
        mse{mse}
    {
        for(int i=0; i<3; i++) this->gradient[i] = gradient(i);
        for (int j=0; j<3; j++) for(int i=0; i<3; i++) this->hessian[i] = hessian(i,j);
        rxy[0] = c.r;
        rxy[1] = c.x;
        rxy[2] = c.y;
        gradient_norm = gradient.norm();
        hessian_norm = hessian.norm();
    }
};

using LmTraceData = std::vector<IterationData>;

LmTraceData estimate_circle_trace(const Dataset& data, const CircleParams& init, const int max_iter=50, double mu=1e-4, const double sigma_inc=10, const double sigma_dec=0.04);

}

}