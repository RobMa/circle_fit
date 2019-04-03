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
        SPDLOG_ERROR("Invalid arguments");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);

    Dataset data(x_vec, y_vec);
    CircleParams est = gwaf_taubin::estimate_circle(data);
    // circle_undo_normalization(data, est);
    *out_x = est.x;
    *out_y = est.y;
    *out_r = est.r;
    return true;
}

int estimate_circle_lm(double x[], double y[], int length, double rxy_init[3], double rxy_out[3])
{
    if(length < 3 || x==nullptr || y == nullptr || rxy_init == nullptr || rxy_out == nullptr) {
        SPDLOG_ERROR("Invalid arguments");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);
    
    Dataset data(x_vec, y_vec);
    // CircleParams init; init.r = ((x_vec.array()-x_vec.mean()).square() + (y_vec.array()-y_vec.mean()).square()).maxCoeff() ; init.x = x_vec.mean(); init.y = y_vec.mean();
    // CircleParams init; init.r = 1; init.x = 0; init.y = 0;
    CircleParams init; init.r = rxy_init[0]; init.x = rxy_init[1]; init.y = rxy_init[2];
    CircleParams est = circle_fit::geometric_lm::estimate_circle(data, init);
    rxy_out[0] = est.r;
    rxy_out[1] = est.x;
    rxy_out[2] = est.y;    
    return true;
}

int estimate_circle(double x[], double y[], int length, double *out_x, double *out_y, double *out_r)
{
    if(length < 3 || x==nullptr || y == nullptr || out_x == nullptr || out_y == nullptr || out_r == nullptr) {
        SPDLOG_ERROR("Invalid arguments");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);
    CircleParams est = circle_fit::estimate_circle(x_vec, y_vec);
    *out_x = est.x;
    *out_y = est.y;
    *out_r = est.r;
    return true;
}

int estimate_circle_lm_trace(double x[], double y[], int length, double rxy_init[3], int *out_numiters, double **out_x, double **out_y, double **out_r, double **out_grad_norm, double **out_mse)
{
    if(length < 3 || x==nullptr || y == nullptr || rxy_init == nullptr || out_numiters == nullptr || out_x == nullptr || out_y == nullptr || out_r == nullptr || out_grad_norm == nullptr || out_mse == nullptr) {
        SPDLOG_ERROR("Invalid arguments");
        return false;
    }
    Eigen::Map<VecX>x_vec(x, length);
    Eigen::Map<VecX>y_vec(y, length);
    Dataset data(x_vec, y_vec);
    CircleParams init; init.r = rxy_init[0]; init.x = rxy_init[1]; init.y = rxy_init[2];
    
    circle_fit::geometric_lm::LmTraceData trace = circle_fit::geometric_lm::estimate_circle_trace(data, init);
    *out_numiters = trace.size();
    *out_x = new double[trace.size()]; //malloc(sizeof(double)*trace.size());
    *out_y = new double[trace.size()];
    *out_r = new double[trace.size()];
    *out_grad_norm = new double[trace.size()];
    *out_mse = new double[trace.size()];
    if (*out_x == nullptr || *out_y == nullptr || *out_r == nullptr || *out_grad_norm == nullptr || *out_mse == nullptr) {
        SPDLOG_ERROR("Memory allocation failed");
        return false;
    }
    SPDLOG_DEBUG("Writing data to output buffer, trace.size()={}...", trace.size());
    for(int i=0; i<trace.size(); i++) {
        // SPDLOG_DEBUG("x");
        (*out_x)[i] = trace[i].rxy[1];
        // SPDLOG_DEBUG("y");
        (*out_y)[i] = trace[i].rxy[2];
        // SPDLOG_DEBUG("r");
        (*out_r)[i] = trace[i].rxy[0];
        // SPDLOG_DEBUG("grad");
        (*out_grad_norm)[i] = trace[i].gradient_norm;
        // SPDLOG_DEBUG("mse");
        (*out_mse)[i] = trace[i].mse;
    }
    SPDLOG_DEBUG("Done.");
    return true;
}

}

