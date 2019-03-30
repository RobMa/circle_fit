#include <circle_fit/gwaf_taubin.h>

namespace circle_fit
{

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
        const double rel_err = fabs(delta/x);
        bool rel_tol_ok = rel_tol>0 && rel_err < rel_tol;
        std::cout << "iter="<<iter <<", x="<<x << ", f="<<f << ", f_dot="<<f_dot << ", delta="<<delta << ", rel_err="<<rel_err <<std::endl;
        if(abs_tol_ok && rel_tol_ok){
            std::cout << "Converged after " << iter << " iterations." << std::endl;
            if(converged != nullptr) {
                *converged = true;
            }
            break;
        }
    }
    return x;
}

struct SolvePoly3{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
    double n = data.n;
    const double Mxx = data.Mxx, 
                 Myy = data.Myy, 
                 Mzz = data.Mzz,
                 Mx  = data.Mx, 
                 My  = data.My,
                 Mz  = data.Mz, 
                 Mxz = data.Mxz,
                 Mxy = data.Mxy, 
                 Myz = data.Myz;
    
    Mat4 M; Mat4 C;
    calculate_M_C(Mx, My, Mz, Mxx, Mxy, Mxz, Myy, Myz, Mzz, n, M, C);
    Vec4 Q3;
    calculate_Q3(Mx, My, Mz, Mxx, Mxy, Mxz, Myy, Myz, Mzz, n, Q3);

    SolvePoly3 solve_poly3(Q3);
    bool newton_converged;
    double eta = newton_solver<SolvePoly3>(solve_poly3, 0, &newton_converged);
    if(!newton_converged) {
        std::cout << "Warning: Newton solver did not converge." << std::endl;
    }
    eta = std::max(0.0, eta);
    std::cout << "eta="<<eta<<std::endl;
    Mat4 Mdet = M - C*eta;
    
    Eigen::JacobiSVD<Mat4> jacobi = Mdet.jacobiSvd(Eigen::ComputeFullV); // Use SVD, because FullPivLU is not robust enough to extract the nullspace!

    const auto& eigenvector = jacobi.matrixV().col(3);
    double normalization = std::sqrt(eigenvector.transpose() * C * eigenvector);
    Vec4 A = eigenvector / normalization;
    return A;
}

// inline void A_to_circle(Vec4 A_vec, double& r, double& x, double& y)
// {
//     const double A = A_vec[0], B = A_vec[1], C = A_vec[2], D = A_vec[3];
//     x = -B/2.0/A;
//     y = -C/2.0/A;
//     r = std::sqrt(std::fabs(x*x + y*y - D/A)); // Formula in paper is wrong
// }

CircleParams gwaf_taubin::estimate_circle(const Dataset& dataset)
{
    ABCDParams abcd; abcd.vec = compute_A(dataset);
    return CircleParams(abcd);
}

// void gwaf_taubin::estimate_circle_normalized(const VecX& x, const VecX& y, double &est_r, double& est_x, double &est_y)
// {
//     double mean_x = x.mean();
//     double mean_y = y.mean();    
//     VecX x_normalized = x.array() - mean_x;
//     VecX y_normalized = y.array() - mean_y;
//     double scale = std::sqrt(std::max(x_normalized.array().square().mean(), y_normalized.array().square().mean()));
//     x_normalized /= scale;
//     y_normalized /= scale;
//     estimate_circle(x_normalized, y_normalized, est_r, est_x, est_y);
//     est_r = est_r*scale;
//     est_x = est_x*scale + mean_x;
//     est_y = est_y*scale + mean_y;
// }

}
