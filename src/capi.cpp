#include <spdlog/spdlog.h>
#include <circle_fit/gwaf_taubin.h>
#include <circle_fit/geo_lm.h>
#include <circle_fit/circle_fit.h>


using namespace circle_fit;

extern "C"
{
int estimate_circle_taubin(double x[], double y[], int length, double *out_x, double *out_y, double *out_r)
{
    if(length < 3 || x==nullptr || y == nullptr || out_x == nullptr || out_y == nullptr || out_r == nullptr) {
        SPDLOG_ERROR("Invalid arguments{0}");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);

    NormalizedDataset data(x_vec, y_vec);
    CircleParams est = gwaf_taubin::estimate_circle(data);
    circle_undo_normalization(data, est);
    *out_x = est.x;
    *out_y = est.y;
    *out_r = est.r;
    return true;
}

int estimate_circle_lm(double x[], double y[], int length, double *out_x, double *out_y, double *out_r)
{
    if(length < 3 || x==nullptr || y == nullptr || out_x == nullptr || out_y == nullptr || out_r == nullptr) {
        SPDLOG_ERROR("Invalid arguments{0}");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);
    
    Dataset data(x_vec, y_vec);
    CircleParams init; init.r = 1; init.x = 0; init.y = 0;
    bool converged = false;
    CircleParams est = circle_fit::geometric_lm::estimate_circle(data, init, &converged);
    *out_x = est.x;
    *out_y = est.y;
    *out_r = est.r;
    return converged;
}

int estimate_circle(double x[], double y[], int length, double *out_x, double *out_y, double *out_r)
{
    if(length < 3 || x==nullptr || y == nullptr || out_x == nullptr || out_y == nullptr || out_r == nullptr) {
        SPDLOG_ERROR("Invalid arguments{0}");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);
    bool converged;
    CircleParams est = circle_fit::estimate_circle(x_vec, y_vec, &converged);
    *out_x = est.x;
    *out_y = est.y;
    *out_r = est.r;
    return converged;
}

}

