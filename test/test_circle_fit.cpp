// #define REQUIRE_MODULE test circle fit common
// #include <boost/test/included/unit_test.hpp>
// #include <boost/test/data/test_case.hpp>
// #include <boost/test/data/monomorphic.hpp>
// namespace bdata = boost::unit_test::data;
// namespace tt = boost::test_tools;



#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <circle_fit/common.h>
using namespace circle_fit;

TEST_CASE( "circle_to_abcd" )
{
    CircleParams circle; circle.r = 1; circle.x = 0; circle.y = 0;
    ABCDParams abcd(circle);
    REQUIRE(abcd.B() == 0);
    REQUIRE(abcd.C() == 0);
    REQUIRE(abcd.A() == 0.5);
    REQUIRE(abcd.D() == -0.5);
}

TEST_CASE( "abcd_to_circle" )
{
    ABCDParams abcd; abcd.A() = 0.5; abcd.D() = -0.5; abcd.B() = 0; abcd.C() = 0;
    CircleParams circle(abcd);
    REQUIRE(circle.r == 1);
    REQUIRE(circle.x == 0);
    REQUIRE(circle.y == 0);
}

TEST_CASE( "adtheta_to_abcd" )
{
    ADThetaParams adtheta;adtheta.A() = 0.5; adtheta.D() = -0.5; adtheta.Theta() = 0;
    ABCDParams abcd(adtheta);
    REQUIRE(abcd.A() == 0.5);
    REQUIRE(abcd.D() == -0.5);
    REQUIRE(abcd.C() == 0);
    REQUIRE(abcd.B() == 0);
}

TEST_CASE( "adtheta_to_circle" )
{
    ADThetaParams adtheta;adtheta.A() = 0.5; adtheta.D() = -0.5; adtheta.Theta() = 0;
    CircleParams circle(adtheta);
    REQUIRE(circle.r == 1);
    REQUIRE(circle.x == 0);
    REQUIRE(circle.y == 0);
}

TEST_CASE ( "circle_to_adtheta" )
{
    CircleParams circle; circle.r = 1; circle.x = 0; circle.y = 0;
    ADThetaParams adtheta(circle);
    REQUIRE(adtheta.A()     == 0.5);
    REQUIRE(adtheta.D()     == -0.5);
    REQUIRE(adtheta.Theta() == 0);
}

TEST_CASE( "abcd_test_e" )
{
    ABCDParams abcd; abcd.A() = 4; abcd.D() = 0.5; abcd.B() = 37; abcd.C() = 42;
    REQUIRE(abcd.E() == 3.0);
}

TEST_CASE( "adtheta_test_e" )
{
    ADThetaParams adtheta; adtheta.A() = 4; adtheta.D() = 0.5; adtheta.Theta() = M_PI;
    REQUIRE(adtheta.E() == 3.0);
}



// BOOST_DATA_TEST_CASE(
//   circle_to_abcd_to_adtheta_to_circle,
//   bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(0.01, 1000)) )
//       ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) )
//       ^ bdata::random( (bdata::distribution=std::uniform_real_distribution<float>(-1000, 1000)) )
//       ^ bdata::xrange(32),
//   r, x, y, index )
// {
//     CircleParams orig; orig.r = r; orig.x = x; orig.y = y;
//     ABCDParams abcd(orig);
//     ADThetaParams adtheta(abcd);
//     CircleParams converted(adtheta);
//     REQUIRE(orig.r == converted.r, tt::tolerance(1e-6));
//     REQUIRE(orig.x == converted.x, tt::tolerance(1e-6));
//     REQUIRE(orig.y == converted.y, tt::tolerance(1e-6));
// }


