#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#include <spdlog/spdlog.h>
#include <Eigen/Dense>
#include <circle_fit/geo_lm.h>

namespace circle_fit
{
using JMat = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using Diag3 = Eigen::DiagonalMatrix<double,3>;

#include "generated/lm_derivatives.hpp"

void shift_origin(VecX& x, VecX& y, Vec2& origin, ADThetaParams& omega)
{
    double offset = 1.0 / sqrt(8.0) / std::fabs(omega.A()); // After shift: 4AD+1 = 0.5
    origin[0] -= offset;
    x = x.array() + offset;
    
    CircleParams circle(omega);
    circle.x += offset;
    omega = ADThetaParams(circle);
}

bool update_state(VecX& x, VecX& y, ADThetaParams& omega, ADThetaParams& omega_pre, Vec3& gradient, Mat3& hessian, Vec2& origin, double& E) {
    const double singularity_threshold = 0.1;
    double singularity_val = omega.A() * omega.D() * 4.0 + 1.0;
    if (singularity_val < 0) {
        shift_origin(x, y, origin, omega_pre);
        SPDLOG_DEBUG("Omega is singular: AD+1={0}, now (previous omega) AD+1={1}, origin_x={2}", singularity_val, omega_pre.A() * omega_pre.D() * 4.0 + 1.0, origin[0]);
        return false;
    }else if (singularity_val <= singularity_threshold) {
        shift_origin(x, y, origin, omega);
        SPDLOG_DEBUG("Close to singularity: pre AD+1={0}, now AD+1={1}, origin_x={2}", singularity_val, omega_pre.A() * omega_pre.D() * 4.0 + 1.0, origin[0]);
    }
    gradient.setZero();
    hessian.setZero();
    E = 0;
    for(int i=0; i<x.size(); i++) {
        double d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta;
        lm_derivatives(x[i], y[i], omega.A(), omega.D(), omega.Theta(), d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta);
        gradient(0) += d * d_dA;
        gradient(1) += d * d_dD;
        gradient(2) += d * d_dTheta;
        hessian(0,0) += (d_dA     * d_dA     + d * d_dA_dA        );
        hessian(0,1) += (d_dD     * d_dA     + d * d_dA_dD        );
        hessian(0,2) += (d_dTheta * d_dA     + d * d_dA_dTheta    );
        hessian(1,1) += (d_dD     * d_dD     + d * d_dD_dD        );
        hessian(1,2) += (d_dTheta * d_dD     + d * d_dD_dTheta    );
        hessian(2,2) += (d_dTheta * d_dTheta + d * d_dTheta_dTheta);
        E += d * d;
    }
    const double normalization = 2 * x.size();
    E /= normalization;
    gradient /= normalization;
    hessian /= normalization;
    hessian(1,0) = hessian(0,1);
    hessian(2,0) = hessian(0,2);
    hessian(2,1) = hessian(1,2);
    return true;
}

void undo_origin_shift(const Vec2& origin, CircleParams& circle) {
    circle.x += origin[0];
    circle.y += origin[1];
}

CircleParams geometric_lm::estimate_circle(const Dataset& data, const CircleParams& init, bool *converged, const double tolerance, const int max_iter, double mu, const double sigma)
{
    if(converged != nullptr) *converged = false;
    ADThetaParams omega(init);
    VecX x(data.x);
    VecX y(data.y);
    Mat3 H;
    Vec2 origin(0, 0);
    Vec3 gradient;
    double E;
    if(!update_state(x, y, omega, omega, gradient, H, origin, E)) {
        SPDLOG_ERROR("Initial circle params are singular.");
        return init;
    }
    double E_previous = E;
    ADThetaParams omega_previous = omega;
    Mat3 H_previous = H;
    Vec3 gradient_previous = gradient;
    const double mu_max = mu * std::pow(sigma, 10); // Break after 10 failed iterations.
    SPDLOG_DEBUG("Initial values: E={0}, mu={1}, A={2}, D={3}, Theta={4}, det|H|={5}, ||grad E||={6}", E, mu, omega.A(), omega.D(), omega.Theta(), H.determinant(), gradient.norm());

    int iter = 1;
    while(iter <= max_iter) {
        // Update weight
        omega.vec.noalias() -= (H + Mat3::Identity()*mu).bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);

        // Analyze error
        bool update_ok = update_state(x, y, omega, omega_previous, gradient, H, origin, E);
        bool iteration_ok = E < E_previous && update_ok;
        SPDLOG_DEBUG("Iteration {0}: OK, E={1}, mu={2}, A={3}, D={4}, Theta={5}, det|H|={6}, ||grad E||={7}", iteration_ok?"OK":"FAIL", E, mu, omega.A(), omega.D(), omega.Theta(), H.determinant(), gradient.norm());
        if(iteration_ok) {
            if(E < tolerance) {
                SPDLOG_DEBUG("Converged after {0} iterations.", iter);
                if(converged != nullptr) *converged = true;
                break;
            }
            E_previous = E;
            omega_previous = omega;
            gradient_previous = gradient;
            H_previous = H;
            iter++;
            mu /= sigma;
        }else
        {
            omega = omega_previous;
            gradient = gradient_previous;
            H = H_previous;
            mu *= sigma;
            iter++;
            if(mu > mu_max) {
                SPDLOG_ERROR("Did not converge, because mu > mu_max.");
                // std::cerr << "Did not converge, because mu > mu_max." << std::endl;
                break;
            }
        }     
    }
    SPDLOG_INFO("Stopped after {0} iterations, E={1}", iter, E);
    CircleParams circle(omega);
    undo_origin_shift(origin, circle);
    return circle;
}

}
