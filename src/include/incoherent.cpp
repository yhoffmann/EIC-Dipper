#pragma once


namespace Incoherent {
    namespace Trans {
        double AIntegrandFunction (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z) {
            return 0;
        }

        double AIntegrandWrapper (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z)
        {
            return AIntegrandFunction(b1,b2,bb1,bb2,r1,r2,rb1,rb2,Q,z);
        }
    }

    namespace Longi {
        double AIntegrandFunction (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z)
        {
            return 0;
        }

        double AIntegrandWrapper (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double z)
        {
            return AIntegrandFunction(b1,b2,bb1,bb2,r1,r2,rb1,rb2,Q,z);
        }
    }
}