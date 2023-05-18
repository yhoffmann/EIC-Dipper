#ifndef COHERENT_CPP
#define COHERENT_CPP

namespace Coherent {
    namespace Trans {
        double AIntegrandFunction (double b1, double b2, double r1, double r2, double Q, double z) {
            return 0;
        }

        double AIntegrandWrapper (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z)
        {
            return AIntegrandFunction(b1,b2,r1,r2,Q,z);
        }
    }

    namespace Longi {
        double AIntegrandFunction (double b1, double b2, double r1, double r2, double Q, double z) {
            return 0;
        }

        double AIntegrandWrapper (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z)
        {
            return AIntegrandFunction(b1,b2,r1,r2,Q,z);
        }
    }
}

#endif