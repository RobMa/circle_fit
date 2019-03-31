#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <circle_fit/gwaf_taubin.h>
using namespace circle_fit;

TEST_CASE( "simple_test" )
{
    const double gt_r = 1000, gt_x = 10000, gt_y = -.012;
    const int number_of_values = 3;
    VecX angles(number_of_values); angles.setRandom(); angles *= M_PI;
    VecX x = gt_r * angles.array().cos() + gt_x;
    VecX y = gt_r * angles.array().sin() + gt_y;
    NormalizedDataset data_norm(x, y);
    CircleParams est = gwaf_taubin::estimate_circle(data_norm);
    circle_undo_normalization(data_norm, est);
    REQUIRE(est.r == Approx(gt_r).epsilon(1e-6));
    REQUIRE(est.x == Approx(gt_x).epsilon(1e-6));
    REQUIRE(est.y == Approx(gt_y).epsilon(1e-6));
}


