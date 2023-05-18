#ifndef MODELS_AND_FUNCTIONS_CPP
#define MODELS_AND_FUNCIONS_CPP


#include<math.h>
#include"constants.cpp"


namespace GBW {
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


#endif