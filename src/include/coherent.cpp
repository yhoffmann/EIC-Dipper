#pragma once

namespace Coherent {
    double AIntegrandFunction (double b1, double b2, double r1, double r2, double Q, double z) {
        std::cout << "getting called" << std::endl;
        return exp(-b1*b1-b2*b2);
    }
}


int coherent_integrand () {

}


std::vector<double> calculate_dsigma_dt_coherent () {
    //TODO implement, loop over trans/longi
}