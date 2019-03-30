#include <gtest/gtest.h>
#include <circle_fit/geo_lm.h>

using namespace circle_fit;

TEST(simple_test, simple_test)
{
    CircleParams gt; gt.r = 1; gt.x = 0; gt.y = 0;
    CircleParams init; init.r = 2; init.x = 0; init.y = 0;
    Dataset data(Vec3(1, 0, -1),
                    Vec3(0, 1, 0));
    CircleParams est = geometric_lm::estimate_circle(data, init);
    ASSERT_NEAR(est.r, gt.r, 1e-6); //), tt::tolerance(1e-6));
    ASSERT_NEAR(est.x, gt.x, 1e-6); //), tt::tolerance(1e-6));
    ASSERT_NEAR(est.y, gt.y, 1e-6); //), tt::tolerance(1e-6));
}

#include <glog/logging.h>
int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    google::ParseCommandLineFlags(&argc, &argv, false);
    google::InitGoogleLogging(argv[0]);
    return RUN_ALL_TESTS();
}

// #define BOOST_TEST_MODULE test circle fit common
// #include <boost/test/included/unit_test.hpp>
// #include <boost/test/data/test_case.hpp>
// #include <boost/test/data/monomorphic.hpp>
// namespace bdata = boost::unit_test::data;
// namespace tt = boost::test_tools;

// #include <circle_fit/geo_lm.h>

// using namespace circle_fit;

// BOOST_AUTO_TEST_CASE( simple_test )
// {
//     CircleParams gt; gt.r = 1; gt.x = 0; gt.y = 0;
//     CircleParams init; init.r = 2; init.x = 0; init.y = 0;
//     Dataset data(Vec3(1, 0, -1),
//                  Vec3(0, 1, 0));
//     CircleParams est = geometric_lm::estimate_circle(data, init);
//     BOOST_TEST(est.r == gt.r, tt::tolerance(1e-6));
//     BOOST_TEST(est.x == gt.x, tt::tolerance(1e-6));
//     BOOST_TEST(est.y == gt.y, tt::tolerance(1e-6));
// }

// BOOST_DATA_TEST_CASE( test_random_circle,
//     bdata::xrange(4, 6)
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(0.01, 1000)) )
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) )
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) ),
//     number_of_points, gt_r, gt_x, gt_y ) __attribute__((optnone))
// {
//     if( number_of_points > 5)return;
//     CircleParams gt; gt.r = gt_r; gt.x = gt_x; gt.y = gt_y;
//     CircleParams init;
//     VecX random_angles(number_of_points);
//     random_angles.setRandom() * M_PI;
//     NormalizedDataset data{random_angles.array().cos() * gt_r + gt_x, random_angles.array().sin() * gt_r + gt_y};
//     init.x = 0; init.y = 0; init.r = 2;
//     CircleParams est = geometric_lm::estimate_circle(data, init);
//     circle_undo_normalization(data, est);
    
    
//     BOOST_TEST(gt.r == est.r, tt::tolerance(1e-5));
//     BOOST_TEST(gt.x == est.x, tt::tolerance(1e-5));
//     BOOST_TEST(gt.y == est.y, tt::tolerance(1e-5));
// }

// #include <circle_fit/gwaf_taubin.h>

// BOOST_DATA_TEST_CASE( test_random_circle_precondition,
//     bdata::xrange(4, 6)
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(0.01, 1000)) )
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) )
//     ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) ),
//     number_of_points, gt_r, gt_x, gt_y ) __attribute__((optnone))
// {
//     // if( number_of_points > 3)return;
//     CircleParams gt; gt.r = gt_r; gt.x = gt_x; gt.y = gt_y;

//     VecX random_angles(number_of_points);
//     random_angles.setRandom() * M_PI;
//     NormalizedDataset data{random_angles.array().cos() * gt_r + gt_x, random_angles.array().sin() * gt_r + gt_y};
//     CircleParams init = gwaf_taubin::estimate_circle(data);
//     CircleParams est = geometric_lm::estimate_circle(data, init);
//     circle_undo_normalization(data, est);
    
    
//     BOOST_TEST(gt.r == est.r, tt::tolerance(1e-5));
//     BOOST_TEST(gt.x == est.x, tt::tolerance(1e-5));
//     BOOST_TEST(gt.y == est.y, tt::tolerance(1e-5));
// }


