#pragma once
#ifndef EIC_COHERENT_HPP_
#define EIC_COHERENT_HPP_

#include <tuple>

#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "IntegrationRoutines.hpp"

namespace Coherent {
double A_integrand_function(double b1, double b2, double r1, double r2,
                            double bb1, double bb2, double rb1, double rb2,
                            double Delta1, double Delta2);

int integrand(unsigned ndim, const double* xx, void* userdata, unsigned fdim,
              double* ff);
double dsigmadt(double Delta, double phi);
double dsigmadt(CubatureConfig* cuba_config,
                IntegrationConfig* integration_config);

double dsigmadt_test(double Delta, double phi = 0.0);
}  // namespace Coherent

namespace Coherent {
namespace Demirci {
double Psi0(double Delta, double m_sqr);
double Psi(double r, double phi, double lambda, double Delta, double m_sqr);
double Z(double r, double phi, double lambda, double Delta, double m_sqr,
         double epsilon);

int integrand(unsigned ndim, const double* xx, void* userdata, unsigned fdim,
              double* ff);
double dsigmadt(double Delta);
double dsigmadt(CubatureConfig* c_config, IntegrationConfig* i_config);
}  // namespace Demirci
}  // namespace Coherent

namespace Coherent {
namespace Sampled {
int integrand(unsigned ndim, const double* xx, void* userdata, unsigned fdim,
              double* ff);
std::tuple<double, double> sqrt_dsigmadt_single_event(
    double Delta, double phi, const HotspotNucleus& nucleus);
std::tuple<double, double> sqrt_dsigmadt_single_event(
    CubatureConfig* c_config, IntegrationConfig* i_config);
}  // namespace Sampled
}  // namespace Coherent

namespace Coherent {
namespace InternalHotspotAvg {
double dsdt_sch_sbch(double Delta);
double dsdt_sch_sbch(CubatureConfig* c_config, IntegrationConfig* i_config);

double dsdt(double Delta);
}  // namespace InternalHotspotAvg
}  // namespace Coherent

#endif  // EIC_COHERENT_HPP_
