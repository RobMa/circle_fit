#pragma once
#include <cmath>
#include "circle_params.h"

// #include "blaze/math/"
namespace circle_fit{

struct Dataset
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VecX x;
    VecX y;
    VecX z;

    int n;

    double Mx;
    double My;
    double Mz;
    double Mxx;
    double Myy;
    double Mzz;
    double Mxy;
    double Mxz;
    double Myz;

    Dataset() = default;
    Dataset(const Dataset&) = default;
    Dataset(Dataset&&) = default;

    Dataset(const VecX& x, const VecX& y) : x{x}, y{y}, z{x.array().square() + y.array().square()}
    {
        update_moments();
    }

    void update_moments()
    {
        n = x.size();
        Mx = x.sum();
        My = y.sum();
        Mz = z.sum();
        Mxx = x.array().square().sum();
        Myy = y.array().square().sum();
        Mzz = z.array().square().sum();
        Mxy = x.dot(y);
        Mxz = x.dot(z);
        Myz = y.dot(z);
    }    
};

struct NormalizedDataset : public Dataset
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double mean_x;
    double mean_y;
    double scale;

    NormalizedDataset(const VecX& x_in, const VecX& y_in)
    {
        mean_x = x_in.mean();
        mean_y = y_in.mean();
        x = x_in.array() - mean_x;
        y = y_in.array() - mean_y;
        scale = std::sqrt(std::max(x.squaredNorm(), y.squaredNorm()));
        this->x /= scale;
        this->y /= scale;
        z = x.array().square() + y.array().square();
        update_moments();
    } 
};

inline void circle_undo_normalization(const NormalizedDataset& dataset, CircleParams& circle)
{
    circle.x = circle.x*dataset.scale + dataset.mean_x;
    circle.y = circle.y*dataset.scale + dataset.mean_y;
    circle.r *= dataset.scale;
}

}

