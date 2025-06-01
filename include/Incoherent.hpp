#pragma once
#ifndef EIC_INCOHERENT_CPP_
#define EIC_INCOHERENT_CPP_

#include "IntegrationRoutines.hpp"

namespace Incoherent {
double dsigmadt(double Delta);
void reset_integration_range(double min[], double max[]);
double dsigmadt(CubatureConfig* cuba_config,
                IntegrationConfig* integration_config);
}  // namespace Incoherent

namespace Incoherent {
namespace Demirci {
double one_disconnected(double Delta);
double two_disconnected(double Delta);

double color_fluctuations(double Delta);
double color_fluctuations(CubatureConfig* c_config,
                          IntegrationConfig* i_config);

double hotspot_fluctuations(double Delta);

double dsigmadt(double Delta);
}  // namespace Demirci
}  // namespace Incoherent

namespace Incoherent {
namespace Sampled {
double dsigmadt_single_event(double Delta, const HotspotNucleus& nucleus);
double dsigmadt_single_event(CubatureConfig* c_config,
                             IntegrationConfig* i_config);
}  // namespace Sampled
}  // namespace Incoherent

namespace Incoherent {
namespace InternalHotspotAvg {
double A_tilde(double b1, double b2, double r1, double r2, double bb1,
               double bb2, double rb1, double rb2, double Delta);

double dsdt_ssbch(double Delta);
double dsdt_ssbch(CubatureConfig* c_config, IntegrationConfig* i_config);

double dsdt_scsbch(double Delta);
double dsdt_scsbch(CubatureConfig* c_config, IntegrationConfig* i_config);

double color_fluc(double Delta);
double hotspot_fluc(double Delta);
double dsdt(double Delta);
}  // namespace InternalHotspotAvg
}  // namespace Incoherent

#endif  // EIC_INCOHERENT_CPP_
