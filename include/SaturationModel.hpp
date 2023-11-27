#pragma once


#include "../external/Nucleus/include/HotspotNucleus.hpp"


namespace SaturationModel
{
    double dsigma_d2b(double x1, double x2, double y1, double y2);

    double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    double dsigma_d2b_sqr_reduced(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    namespace Sampled
    {
        double dsigma_d2b(double x1, double x2, double y1, double y2, const HotspotNucleus* nucleus);

        double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus);

        double dsigma_d2b_sqr_reduced(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus);
    }
}