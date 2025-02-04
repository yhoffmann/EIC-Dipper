#pragma once


#include "IntegrationRoutines.hpp"


namespace Incoherent {
    double A_integrand_function_factor(double Q);
    double dsigmadt(double Q, double Delta);
    void reset_integration_range(double min[], double max[]);
    double dsigmadt(CubatureConfig* cuba_config, IntegrationConfig* integration_config);
}


namespace Incoherent { namespace Demirci
{
    double one_disconnected (double Q, double Delta);
    double two_disconnected (double Q, double Delta);

    double color_fluctuations (double Q, double Delta);
    double color_fluctuations (CubatureConfig* c_config, IntegrationConfig* i_config);

    double hotspot_fluctuations (double Q, double Delta);

    double dsigmadt (double Q, double Delta);
} }


namespace Incoherent { namespace Sampled
{
    double dsigmadt_single_event (double Q, double Delta, const HotspotNucleus& nucleus);
    double dsigmadt_single_event (CubatureConfig* c_config, IntegrationConfig* i_config);
} }


namespace Incoherent { namespace InternalHotspotAvg
{
    double A_tilde(double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2, double Q, double Delta);

    double dsdt_ssbch(double Q, double Delta);
    double dsdt_ssbch(CubatureConfig* c_config, IntegrationConfig* i_config);
    
    double dsdt_scsbch(double Q, double Delta);
    double dsdt_scsbch(CubatureConfig* c_config, IntegrationConfig* i_config);

    double color_fluc(double Q, double Delta);
    double hotspot_fluc(double Q, double Delta);
    double dsdt(double Q, double Delta);
} }