#pragma once


#include "../include/utilities.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "../include/NRPhoton.hpp"


namespace SaturationModel
{
    double dsigma_d2b(double x1, double x2, double y1, double y2);

    // double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    double dsigma_d2b_sqr_reduced(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

    namespace Sampled
    {
        double dsigma_d2b(double x1, double x2, double y1, double y2, const HotspotNucleus* nucleus);

        // double dsigma_d2b_sqr(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus);

        double dsigma_d2b_sqr_reduced(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2, const HotspotNucleus* nucleus);
    }

    namespace InternalHotspotAvg
    {
        inline uint num = 0;
        inline double inv_num;
        inline double inv_num_sqr;
        inline std::vector<HotspotNucleus> hn;
        void init(uint A, uint H, uint num, uint start_seed = g_seed);
        void clear();

        double dsigma_d2b(double x1, double x2, double y1, double y2);
        double dsigma_d2b_sqr_reduced(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);

        double sch_sbch(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);
        double ssbch(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);
        double scsbch(double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2);
    }
}