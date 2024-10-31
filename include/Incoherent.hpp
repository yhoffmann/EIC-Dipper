#pragma once


#include "IntegrationRoutines.hpp"


namespace Incoherent {
    double dsigmadt(double Q, double Delta);
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