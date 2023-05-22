#pragma once

namespace Coherent {
    namespace Trans {
        double AIntegrandFunction (double b1, double b2, double r1, double r2, double Q, double z) {
            std::cout << "getting called" << std::endl;
            return exp(-b1*b1-b2*b2);
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