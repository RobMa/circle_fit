#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <circle_fit/circle_params.h>
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

TEST_CASE("circle_to_abcd_to_adtheta_to_circle")
{
    using namespace Catch::Generators;
    const auto repeats = 1000UL;
    auto r = take(repeats, random(0.01, 1000.0));
    auto x = take(repeats, random(-1000.0, 1000.0));
    auto y = take(repeats, random(-1000.0, 1000.0));
    do{
        CircleParams orig; orig.r = r.get(); orig.x = x.get(); orig.y = y.get();
        ABCDParams abcd(orig);
        ADThetaParams adtheta(abcd);
        CircleParams converted(adtheta);
        REQUIRE(orig.r == Approx(converted.r).epsilon(1e-6).margin(1e-6));
        REQUIRE(orig.x == Approx(converted.x).epsilon(1e-6).margin(1e-6));
        REQUIRE(orig.y == Approx(converted.y).epsilon(1e-6).margin(1e-6));
    }while(r.next() && x.next() && y.next());
}

