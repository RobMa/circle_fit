#include <spdlog/spdlog.h>
#include <cmath>
#include <Eigen/SVD>
#include <circle_fit/gwaf_taubin.h>

namespace circle_fit{

template<typename Context>
inline double newton_solver(const Context& c, const double x0 = 0, bool *converged = nullptr, const double max_iter = 12, const double abs_tol = 1e-12, const double rel_tol = 1e-12)
{
    if(converged != nullptr) {
        *converged = false;
    }
    double x = x0;
    double f, f_dot;
    for(int iter=1; iter<=max_iter; iter++) {
        c.evaluate(x, f, f_dot);
        
        double delta = f/f_dot;
        x -= delta;
        bool abs_tol_ok = abs_tol>0 && fabs(delta) < abs_tol;
        const double delta_rel = fabs(delta/x);
        bool rel_tol_ok = rel_tol>0 && delta_rel < rel_tol;
        SPDLOG_DEBUG("Iteration {}, x={}, f={}, f_dot={}, delta={}, delta_rel={}", iter, x, f, f_dot, delta, delta_rel);
        // std::cout << "iter="<<iter <<", x="<<x << ", f="<<f << ", f_dot="<<f_dot << ", delta="<<delta << ", delta_rel="<<delta_rel <<std::endl;
        if(abs_tol_ok && rel_tol_ok){
            SPDLOG_INFO("Converged after {} iterations, delta_abs={}, delta_rel={}", iter, delta, delta_rel);
            // std::cout << "Converged after " << iter << " iterations." << std::endl;
            if(converged != nullptr) {
                *converged = true;
            }
            break;
        }
    }
    return x;
}

struct SolvePoly3{
    Vec4 coeffs;
    SolvePoly3(const Vec4& coeffs) : coeffs{coeffs} {}
    inline void evaluate(const double x, double& f, double& f_dot) const
    {
        f = 0;
        f_dot = 0;
        double x_pow = 1;
        for (int i=0; i<4; i++) {
            f += coeffs[i] * x_pow;
            if(i<3){
                f_dot += double(i+1) * coeffs[i+1] *  x_pow;
                x_pow *= x;
            }            
        }
    }
};

#include "generated/M_C.hpp"
#include "generated/Q3.hpp"

inline Vec4 compute_A(const Dataset& data)
{
    Mat4 M; Mat4 C;
    calculate_M_C(data.Mx, data.My, data.Mz, data.Mxx, data.Mxy, data.Mxz, data.Myy, data.Myz, data.Mzz, data.n, M, C);
    Vec4 Q3;
    calculate_Q3(data.Mx, data.My, data.Mz, data.Mxx, data.Mxy, data.Mxz, data.Myy, data.Myz, data.Mzz, data.n, Q3);

    SolvePoly3 solve_poly3(Q3);
    bool newton_converged;
    double eta = newton_solver<SolvePoly3>(solve_poly3, 0, &newton_converged);
    if(!newton_converged) {
        SPDLOG_ERROR("Warning: Newton solver did not converge.");
    }
    eta = std::max(0.0, eta);
    SPDLOG_DEBUG("eta={}", eta);
    Mat4 Mdet = M - C*eta;
    Eigen::JacobiSVD<Mat4> jacobi = Mdet.jacobiSvd(Eigen::ComputeFullV); // Use SVD, because FullPivLU is not robust enough to extract the nullspace!
    SPDLOG_DEBUG("Smallest singular value = {}", jacobi.singularValues().data()[3]);
    const auto& eigenvector = jacobi.matrixV().col(3);
    double normalization = std::sqrt(eigenvector.transpose() * C * eigenvector);
    Vec4 A = eigenvector / normalization;
    return A;
}

CircleParams gwaf_taubin::estimate_circle(const Dataset& dataset)
{
    ABCDParams abcd; abcd.vec = compute_A(dataset);
    return CircleParams(abcd);
}

}
