#pragma once

#include "numeric_types.h"

namespace circle_fit
{

struct CircleParams; struct ABCDParams; struct ADThetaParams;
inline void circle_params_to_abcd_params(const CircleParams& circle, ABCDParams& out);
inline void abcd_params_to_circle_params(const ABCDParams& abcd, CircleParams& out);
inline void abcd_params_to_adtheta_params(const ABCDParams& abcd, ADThetaParams& out);
inline void adtheta_params_to_abcd_params(const ADThetaParams& adtheta, ABCDParams& out);

struct CircleParams
{
    double x;
    double y;
    double r;
    CircleParams() = default;
    inline explicit CircleParams(const ABCDParams& abcd);
    inline explicit CircleParams(const ADThetaParams& adtheta);
};

struct ABCDParams
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Vec4 vec;
    ABCDParams() = default;
    inline explicit ABCDParams(const CircleParams& circle);
    inline explicit ABCDParams(const ADThetaParams& adtheta);
    inline double A() const {return vec[0];};
    inline double B() const {return vec[1];};
    inline double C() const {return vec[2];};
    inline double D() const {return vec[3];};

    inline double E() const {return std::sqrt(1.0+4.0*A()*D());}

    inline double& A() {return vec[0];};
    inline double& B() {return vec[1];};
    inline double& C() {return vec[2];};
    inline double& D() {return vec[3];};
};

struct ADThetaParams
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Vec3 vec;
    ADThetaParams() = default;
    inline explicit ADThetaParams(const CircleParams& circle) {
        ABCDParams abcd(circle);
        abcd_params_to_adtheta_params(abcd, *this);
    }
    inline explicit ADThetaParams(const ABCDParams& abcd) {
        abcd_params_to_adtheta_params(abcd, *this);
    }

    inline double A() const {return vec[0];};
    inline double D() const {return vec[1];};
    inline double Theta() const {return vec[2];};
    inline double E() const {return std::sqrt(1.0+4.0*A()*D());}

    inline double& A() {return vec[0];};
    inline double& D() {return vec[1];};
    inline double& Theta() {return vec[2];};
};

inline CircleParams::CircleParams(const ABCDParams& abcd){
    abcd_params_to_circle_params(abcd, *this);
}
inline CircleParams::CircleParams(const ADThetaParams& adtheta) {
    ABCDParams abcd(adtheta);
    abcd_params_to_circle_params(abcd, *this);
}
inline ABCDParams::ABCDParams(const CircleParams& circle) {
    circle_params_to_abcd_params(circle, *this);
}
inline ABCDParams::ABCDParams(const ADThetaParams& adtheta) {
    adtheta_params_to_abcd_params(adtheta, *this);
}

inline void circle_params_to_abcd_params(const CircleParams& circle, ABCDParams& out)
{
    out.A() = 1.0/2.0/circle.r;
    out.B() = -2.0 * out.A() * circle.x;
    out.C() = -2.0 * out.A() * circle.y;
    out.D() = (out.B()*out.B() + out.C()*out.C() - 1 ) / 4.0 / out.A();
} 

inline void abcd_params_to_circle_params(const ABCDParams& abcd, CircleParams& out)
{
    out.x = -abcd.B()/2.0/abcd.A();
    out.y = -abcd.C()/2.0/abcd.A();
    out.r = std::sqrt(std::fabs(out.x*out.x + out.y*out.y - abcd.D()/abcd.A())); // Formula in paper is only applicable for normalized ABCD
}

inline void abcd_params_to_adtheta_params(const ABCDParams& abcd, ADThetaParams& out)
{
    out.A() = abcd.A();
    out.D() = abcd.D();
    out.Theta() = abcd.E() != 0 ? atan2(abcd.C()/abcd.E(), abcd.B()/abcd.E()) : 0;
}

inline void adtheta_params_to_abcd_params(const ADThetaParams& adtheta, ABCDParams& out)
{
    out.A() = adtheta.A();
    out.D() = adtheta.D();
    out.B() = std::cos(adtheta.Theta()) * adtheta.E();
    out.C() = std::sin(adtheta.Theta()) * adtheta.E();
}

}

