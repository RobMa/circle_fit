#include <spdlog/spdlog.h>
#include <Eigen/Dense>
#include <circle_fit/geo_lm.h>
#include <Eigen/Cholesky>


namespace circle_fit
{
using JMat = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using Diag3 = Eigen::DiagonalMatrix<double,3>;

#include "generated/lm_derivatives.hpp"

void shift_origin(VecX& x, VecX& y, Vec2& origin, ADThetaParams& omega)
{
    // double offset = x.array().square().mean();
    double offset = 1;
    if(omega.A() != 0) {
        offset = 1.0 / sqrt(8.0) / std::fabs(omega.A());
    }else{
        SPDLOG_ERROR("Invalid omega: A == 0");
    }
    SPDLOG_DEBUG("Shifting origin by x0 -= {}.", offset);
    origin[0] -= offset;
    x = x.array() + offset;    
    CircleParams circle(omega);
    circle.x += offset;
    omega = ADThetaParams(circle);
}

inline void undo_origin_shift(const Vec2& origin, CircleParams& circle) {
    circle.x += origin[0];
    circle.y += origin[1];
}

inline CircleParams undo_origin_shift(const Vec2& origin, CircleParams&& circle) {
    CircleParams c = std::move(circle);
    undo_origin_shift(origin, c);
    return c;
}

enum class SingularityStatus{
    Singular,
    AlmostSingular,
    Regular
};

SingularityStatus is_singular(const ADThetaParams& omega, double threshold = 0.1)
{
    double val = omega.A() * omega.D() * 4.0 + 1.0;
    SingularityStatus status = val > threshold ? SingularityStatus::Regular : 
        val < 0 ? SingularityStatus::Singular : SingularityStatus::AlmostSingular;
    if (status != SingularityStatus::Regular) {
        SPDLOG_DEBUG("Omega is (almost) singular: 4AD+1={0}", val);
    }
    return status;
}

double mean_squared_error(const VecX& x, const VecX& y, const ADThetaParams& omega)
{
    CircleParams c(omega);
    double mse = (((x.array()-c.x).square() + (y.array()-c.y).square()).sqrt() - c.r).square().sum() / x.size();
    return mse;
}

void mean_squared_error_derivatives(const VecX& x, const VecX& y, const ADThetaParams& omega, double& mse_out, Vec3& gradient_out, Mat3& hessian_out)
{
    gradient_out.setZero();
    hessian_out.setZero();
    mse_out = 0;
    for(int i=0; i<x.size(); i++) {
        double d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta;
        lm_derivatives(x[i], y[i], omega.A(), omega.D(), omega.Theta(), d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta);
        gradient_out(0) += d * d_dA;
        gradient_out(1) += d * d_dD;
        gradient_out(2) += d * d_dTheta;
        hessian_out(0,0) += (d_dA     * d_dA     + d * d_dA_dA        );
        hessian_out(0,1) += (d_dD     * d_dA     + d * d_dA_dD        );
        hessian_out(0,2) += (d_dTheta * d_dA     + d * d_dA_dTheta    );
        hessian_out(1,1) += (d_dD     * d_dD     + d * d_dD_dD        );
        hessian_out(1,2) += (d_dTheta * d_dD     + d * d_dD_dTheta    );
        hessian_out(2,2) += (d_dTheta * d_dTheta + d * d_dTheta_dTheta);
        // hessian_out(0,0) += (d_dA     * d_dA     ); //+ d * d_dA_dA        );
        // hessian_out(0,1) += (d_dD     * d_dA     ); //+ d * d_dA_dD        );
        // hessian_out(0,2) += (d_dTheta * d_dA     ); //+ d * d_dA_dTheta    );
        // hessian_out(1,1) += (d_dD     * d_dD     ); //+ d * d_dD_dD        );
        // hessian_out(1,2) += (d_dTheta * d_dD     ); //+ d * d_dD_dTheta    );
        // hessian_out(2,2) += (d_dTheta * d_dTheta ); //+ d * d_dTheta_dTheta);
        mse_out += d * d;
    }
    const double normalization = 1.0 / x.size();
    mse_out *= normalization;
    gradient_out *= normalization * 2.0; // Factor 2 comes from square derivative
    hessian_out *= normalization * 2.0;
    hessian_out(1,0) = hessian_out(0,1);
    hessian_out(2,0) = hessian_out(0,2);
    hessian_out(2,1) = hessian_out(1,2);
}

CircleParams geometric_lm::estimate_circle(const Dataset& data, const CircleParams& init, const int max_iter, double mu, const double sigma_inc, const double sigma_dec)
{
    VecX x(data.x);
    VecX y(data.y);
    ADThetaParams omega(init);
    Vec2 origin(0, 0);
    switch(is_singular(omega)){
        case SingularityStatus::Regular:
            break;
        case SingularityStatus::Singular:
            SPDLOG_ERROR("Initial circle params are singular.");
            return init;
        case SingularityStatus::AlmostSingular:
            shift_origin(x, y, origin, omega);
    }
    
    const double mu_max = 1e20;
    const double mu_min = 1e-20;

    int iter = 1;
    bool converged = false;
    double mse = 0; 
    while(iter <= max_iter && !converged) {
        Mat3 hessian; Vec3 gradient;
        mean_squared_error_derivatives(x, y, omega, mse, gradient, hessian);
        // mse = mean_squared_error(x,y,omega);
        SPDLOG_DEBUG("Iteration {}: MSE={}, mu={}, A={}, D={}, Theta={}, det|H|={}, ||grad E||={}", 
                    iter-1, mse, mu, omega.A(), omega.D(), omega.Theta(), hessian.determinant(), gradient.norm());

        // Decrement mu
        mu = std::max(mu * sigma_dec, mu_min);

        while(!converged)
        {
            // Increment mu
            mu = std::min(mu * sigma_inc, mu_max);

            // Update weight
            // omega.vec.noalias() -= (H + Mat3::Identity()*mu).jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);
            // omega.vec.noalias() -= (H + Mat3::Identity()*mu).ch().solve(gradient);
            // omega.vec.noalias() -= Eigen::JacobiSVD<Mat3>(H + (mu*Diag3(1.0+H.diagonal()[0], 1.0+H.diagonal()[1], 1.0+H.diagonal()[2])).toDenseMatrix(), Eigen::ComputeFullV|Eigen::ComputeFullU).solve(gradient);
            ADThetaParams omega_new;
            Mat3 nash_offset = mu * Diag3(1.0+hessian(0,0), 1.0+hessian(1,1), 1.0+hessian(2,2));
            omega_new.vec = omega.vec - Eigen::LLT<Mat3>(hessian + nash_offset).solve(gradient);
            if(is_singular(omega_new) != SingularityStatus::Regular)
            {
                // We are moving into a singularity in parameter space.
                // Shift the origin and try again.
                shift_origin(x, y, origin, omega);
                mean_squared_error_derivatives(x, y, omega, mse, gradient, hessian);
                continue;
            }
            double mse_new = mean_squared_error(x, y, omega_new);
            SPDLOG_TRACE("Iteration {}, mu={}, MSE_new={}, MSE_old={}", iter, mu, mse_new, mse);
            if(mse_new < mse) {
                omega = omega_new;
                iter++;
                break;
            }else if (mu == mu_max){
                converged = true;
            }
        }
    }
    return undo_origin_shift(origin, CircleParams(omega));
}

using geometric_lm::IterationData;
using geometric_lm::LmTraceData;

LmTraceData geometric_lm::estimate_circle_trace(const Dataset& data, const CircleParams& init, const int max_iter, double mu, const double sigma_inc, const double sigma_dec)
{
    LmTraceData trace;
    VecX x(data.x);
    VecX y(data.y);
    ADThetaParams omega(init);
    Vec2 origin(0, 0);
    switch(is_singular(omega)){
        case SingularityStatus::Regular:
            break;
        case SingularityStatus::Singular:
            SPDLOG_ERROR("Initial circle params are singular.");
            return trace;
        case SingularityStatus::AlmostSingular:
            shift_origin(x, y, origin, omega);
    }
    
    const double mu_max = 1e20;
    const double mu_min = 1e-20;

    int iter = 1;
    bool converged = false;
    double mse = 0; 
    while(iter <= max_iter && !converged) {
        Mat3 hessian; Vec3 gradient;
        mean_squared_error_derivatives(x, y, omega, mse, gradient, hessian);
        mse = mean_squared_error(x,y,omega);
        SPDLOG_DEBUG("Iteration {}: MSE={}, mu={}, A={}, D={}, Theta={}, det|H|={}, ||grad E||={}", 
                    iter-1, mse, mu, omega.A(), omega.D(), omega.Theta(), hessian.determinant(), gradient.norm());

        // Decrement mu
        mu = std::max(mu * sigma_dec, mu_min);

        while(!converged)
        {
            // Increment mu
            mu = std::min(mu * sigma_inc, mu_max);

            // Update weight
            // omega.vec.noalias() -= (H + Mat3::Identity()*mu).jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);
            // omega.vec.noalias() -= (H + Mat3::Identity()*mu).ch().solve(gradient);
            // omega.vec.noalias() -= Eigen::JacobiSVD<Mat3>(H + (mu*Diag3(1.0+H.diagonal()[0], 1.0+H.diagonal()[1], 1.0+H.diagonal()[2])).toDenseMatrix(), Eigen::ComputeFullV|Eigen::ComputeFullU).solve(gradient);
            ADThetaParams omega_new;
            Mat3 nash_offset = mu * Diag3(1.0+hessian(0,0), 1.0+hessian(1,1), 1.0+hessian(2,2));
            omega_new.vec = omega.vec - Eigen::LLT<Mat3>(hessian + nash_offset).solve(gradient);
            if(is_singular(omega_new) != SingularityStatus::Regular)
            {
                // We are moving into a singularity in parameter space.
                // Shift the origin and try again.
                shift_origin(x, y, origin, omega);
                continue;
            }
            double mse_new = mean_squared_error(x, y, omega_new);
            SPDLOG_TRACE("Iteration {}, mu={}, MSE_new={}, MSE_old={}", iter, mu, mse_new, mse);
            if(mse_new < mse) {
                omega = omega_new;
                iter++;
                trace.push_back(IterationData(
                    mse, gradient, hessian, undo_origin_shift(origin, CircleParams(omega))
                ));
                break;
            }else if (mu == mu_max){
                converged = true;
            }
        }
    }
    return trace;
}

}
