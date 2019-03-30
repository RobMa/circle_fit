#define BOOST_TEST_MODULE Test GWAF Taubin
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

#include <cstdlib>
#include <circle_fit/gwaf_taubin.h>
using namespace circle_fit;

BOOST_AUTO_TEST_CASE( simple_test )
{
    const double gt_r = 1000, gt_x = 10000, gt_y = -.012;
    const int number_of_values = 3;
    VecX angles(number_of_values); angles.setRandom(); angles *= M_PI;
    VecX x = gt_r * angles.array().cos() + gt_x;
    VecX y = gt_r * angles.array().sin() + gt_y;
    NormalizedDataset data_norm(x, y);
    CircleParams est = gwaf_taubin::estimate_circle(data_norm);
    circle_undo_normalization(data_norm, est);
    BOOST_TEST(est.r == gt_r, tt::tolerance(1e-6));
    BOOST_TEST(est.x == gt_x, tt::tolerance(1e-6));
    BOOST_TEST(est.y == gt_y, tt::tolerance(1e-6));
}


