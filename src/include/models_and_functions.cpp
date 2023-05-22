#pragma once


#include<math.h>


namespace GBWModel {
    double T_times_sigma0 (double b1, double b2) {
        return exp( -(b1*b1+b2*b2)/(2.0*BG) );
    }
};


namespace Photon {
    double epsilon (double Q, double z) {

    }

    namespace Trans {

    };

    namespace Longi {

    };
};