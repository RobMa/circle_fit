#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <circle_fit/geo_lm.h>
using namespace circle_fit;

TEST_CASE("simple_test")
{
    CircleParams gt; gt.r = 1; gt.x = 0; gt.y = 0;
    CircleParams init; init.r = 2; init.x = 0; init.y = 0;
    Dataset data(Vec3(1, 0, -1),
                 Vec3(0, 1, 0));
    CircleParams est = geometric_lm::estimate_circle(data, init);
    REQUIRE(Approx(est.r).epsilon(1e-6).margin(1e-6) == gt.r);
    REQUIRE(Approx(est.x).epsilon(1e-6).margin(1e-6) == gt.x);
    REQUIRE(Approx(est.y).epsilon(1e-6).margin(1e-6) == gt.y);
}

TEST_CASE("simple_test2")
{
    CircleParams gt; gt.r = 1; gt.x = 0; gt.y = 0;
    CircleParams init; init.r = 2; init.x = 1; init.y = 1;
    Dataset data(Vec3(1, 0, -1),
                 Vec3(0, 1, 0));
    CircleParams est = geometric_lm::estimate_circle(data, init);
    REQUIRE(Approx(est.r).epsilon(1e-6).margin(1e-6) == gt.r);
    REQUIRE(Approx(est.x).epsilon(1e-6).margin(1e-6) == gt.x);
    REQUIRE(Approx(est.y).epsilon(1e-6).margin(1e-6) == gt.y);
}

