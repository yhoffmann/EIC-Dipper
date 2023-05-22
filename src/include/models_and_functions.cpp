#pragma once


#include<math.h>


namespace GBWModel {
    double T_times_sigma0 (double b1, double b2) {
        return exp( -(b1*b1+b2*b2)/(2.0*BG) );
    }
};


namespace NRPhoton {
    double epsilon (double Q, double z) {

    }

    double wave_function (double r1, double r1, double Q, double z, bool t_not_l) {
        if (t_not_l) {
            return 0;
        } else {
            return 0;
        }
    }
};