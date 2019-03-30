#include <Eigen/Dense>
#include <circle_fit/geo_lm.h>
#include <iostream>
#include <glog/logging.h>

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
        VLOG(1) << "Omega is singular: AD+1=" << singularity_val << ", now (previous omega) AD+1=" << omega_pre.A() * omega_pre.D() * 4.0 + 1.0 << ", origin_x="<<origin[0]<<std::endl;
        return false;
    }else if (singularity_val <= singularity_threshold) {
        shift_origin(x, y, origin, omega);
        VLOG(1) << "Close to singularity: pre AD+1=" << singularity_val << ", now AD+1=" << omega.A() * omega.D() * 4.0 + 1.0 << ", origin_x="<<origin[0]<< std::endl;
    }
    gradient.setZero();
    hessian.setZero();
    E = 0;
    for(int i=0; i<x.size(); i++) {
        double d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta;
        lm_derivatives(x[i], y[i], omega.A(), omega.D(), omega.Theta(), d, d_dA, d_dD, d_dTheta, d_dA_dA, d_dA_dD, d_dA_dTheta, d_dD_dD, d_dD_dTheta, d_dTheta_dTheta);
        gradient(0) += 2.0 * d * d_dA;
        gradient(1) += 2.0 * d * d_dD;
        gradient(2) += 2.0 * d * d_dTheta;
        hessian(0,0) += 2.0 * (d_dA     * d_dA     + d * d_dA_dA        );
        hessian(0,1) += 2.0 * (d_dD     * d_dA     + d * d_dA_dD        );
        hessian(0,2) += 2.0 * (d_dTheta * d_dA     + d * d_dA_dTheta    );
        hessian(1,1) += 2.0 * (d_dD     * d_dD     + d * d_dD_dD        );
        hessian(1,2) += 2.0 * (d_dTheta * d_dD     + d * d_dD_dTheta    );
        hessian(2,2) += 2.0 * (d_dTheta * d_dTheta + d * d_dTheta_dTheta);
        E += d * d;
    }
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
        LOG(ERROR) << "Initial circle params are singular." << std::endl;
        return init;
    }
    double E_previous = E;
    ADThetaParams omega_previous = omega;
    Mat3 H_previous = H;
    Vec3 gradient_previous = gradient;
    const double mu_max = mu * std::pow(sigma, 10); // Break after 10 failed iterations.
    DVLOG(1) << "Initial params: " << "A=" << omega.A() << ", D=" << omega.D() << ", Theta=" << omega.Theta() << ", Error E=" << E_previous << ", det|H|=" << H.determinant()<< ", ||G||="<< gradient.norm() << std::endl;

    int iter = 1;
    while(iter <= max_iter) {
        // Update weight
        omega.vec.noalias() -= (H + Mat3::Identity()*mu).bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);

        // Analyze error
        double E;
        bool update_ok = update_state(x, y, omega, omega_previous, gradient, H, origin, E);

        if(E < E_previous && update_ok) {
            // Successful iteration
            DVLOG(2) << "Iteration " << iter << ": OK, E=" << E << ", mu=" << mu << ", A=" << omega.A() << ", D=" << omega.D() << ", Theta=" << omega.Theta()<< ", det|H|=" << H.determinant()<< ", ||G||="<< gradient.norm() << std::endl;
            if(E < tolerance) {
                DVLOG(1) << "Converged after " << iter << " iterations."<<std::endl;
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
            DVLOG(2) << "Iteration " << iter << ": FAIL, E=" << E << ", mu=" << mu << ", A=" << omega.A() << ", D=" << omega.D() << ", Theta=" << omega.Theta()<< ", det|H|=" << H.determinant()<< ", ||G||="<< gradient.norm() << std::endl;
            omega = omega_previous;
            gradient = gradient_previous;
            H = H_previous;
            mu *= sigma;
            iter++;
            if(mu > mu_max) {
                std::cerr << "Did not converge, because mu > mu_max." << std::endl;
                break;
            }
        }     
    }
    DVLOG(1) << "Stopped after "<< iter-1 <<" iterations, E = " << E << std::endl;
    CircleParams circle(omega);
    undo_origin_shift(origin, circle);
    return circle;
}

}
