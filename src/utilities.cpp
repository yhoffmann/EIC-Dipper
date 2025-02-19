#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../include/IntegrationRoutines.hpp"
#include "../include/NRPhoton.hpp"
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>


void init(int argc, char* argv[])
{
    g_time_program_start = std::chrono::high_resolution_clock::now();
    set_parameters(argc, argv);
}


double get_r_max (double m_Q)
{
    return R_RANGE_FACTOR/m_Q;
}


void import_interp_data_by_params (const Interpolator3D::DataGenerationConfig* config)
{
    std::string filepath;

    std::string m_path = std::to_string( int(m*100) );
    std::string rH_sqr_path = (rH_sqr >= 1.0) ? std::to_string( int(std::round(rH_sqr*100)) ) : "0"+std::to_string( int(std::round(rH_sqr*100)) );

    filepath = "interpolator-data/";

    filepath += "G_rH2_"+rH_sqr_path+"_m_0"+m_path+".dat";

    std::ifstream file_check(filepath);

    if (!config)
        return;

    if (!file_check)
    {
        file_check.close();
#ifndef _QUIET
        std::string answer;

        std::cout << "No data for this combination of m(=" << m << ") and rH2(=" << rH_sqr << "). Do you want to generate that data set? (y/n)\n"
            "(The path to the new file would be " << filepath << ". File size would be " << ((config->nx+3)*(config->ny+3)*(config->nz+3)+config->nx+config->ny+config->nz+13)*8.0*1.0e-6 << "MB.)" << std::endl;
        bool input_accepted = false;
        while (!input_accepted)
        {
            std::cin >> answer;
            if (answer=="y")
            {
                input_accepted = true;

                std::cout << "Generating new data set..." << std::endl;
    #endif    
                GBWModel::G_ip.generate_data(GBWModel::G_wrapper, config, true);
                GBWModel::G_ip.export_data(filepath);
    #ifndef _QUIET
            }
            else if (answer=="n")
            {
                input_accepted = true;

                std::cout << "Aborting" << std::endl;
                exit(21);
            }
            else
            {
                std::cout << "Please enter either \"y\" or \"n\"." << std::endl;
            }
        }
#endif
    }
    else
    {
_TEST_LOG("Starting interpolator data import from file: " << filepath)
    GBWModel::G_ip.import_data(filepath);
_TEST_LOG("Finished interpolator data import")
    }
}


void set_parameters (int argc, char** argv)
{
    _TEST_LOG("Setting parameters")

    const std::string error_message = "Invalid use. Valid flags are\n"
                "\t[-s <seed>] (rng seed)\n"
                "\t[--add-to-seed <number>] (add <number> to seed, use after -s)\n"
                "\t[-t, --threads <number>] (run any multithreaded operation with <number> threads)\n"
                "\t[-p] (print progress and intermediate values to stdout)\n"
                "\t[-Q <photon virtuality>]\n"
                "\t[-Delta <Delta>] (output only for this single Delta value)\n"
                "\t[-m <gluon mass>]\n"
                "\t[--g2mu02-factor <number>] (multiplies base value of g2mu02 by <number>)\n"
                "\t[--charm, -c] (select charm quark)\n"
                "\t[--bottom, -b] (select bottom quark)\n"
                "\t[-A <atomic number>]\n"
                "\t[-H <number of hotspots per nucleon>]\n"
                "\t[-rH2 <hotspot radius square>]\n"
                "\t[-Rp2 <nucleon radius square>]\n"
                "\t[-o <output filepath>]";

    for (int i=1; i<argc; i+=2)
    {
        std::istringstream flag(argv[i]);

        if (flag.str()=="-p")
        {
            g_monitor_progress = true;
            --i;
            continue;
        }
        else if (flag.str()=="--bottom" || flag.str()=="-b")
        {
            m_Q = m_b;
            e_Q = e_b;

            R_MAX = get_r_max(m_b);
            NRPhoton::set_wave_function_factor_T(m_b, e_b);

            quark_config = Quark::Bottom;

            --i;
            continue;
        }
        else if (flag.str()=="--charm" || flag.str()=="-c")
        {
            m_Q = m_c;
            e_Q = e_c;

            R_MAX = get_r_max(m_c);
            NRPhoton::set_wave_function_factor_T(m_c, e_c);

            quark_config = Quark::Charm;

            --i;
            continue;
        }

        if (i+1==argc)
        {
            std::cerr << error_message << std::endl;
            exit(23);
        }

        std::istringstream arg(argv[i+1]);
        std::string arg_string;

        if (flag.str()=="-o" && (arg >> arg_string))
        {
            g_filepath = arg_string;
            continue;
        }

        double arg_number;

        if ( !(arg >> arg_number) || !arg.eof() )
        {
            std::cerr << "You entered the argument \"" << arg.str() << "\" for the option " << "\"" << flag.str() << "\". Either the option or the argument are invalid." << std::endl;
            std::cerr << error_message << std::endl;
            exit(23);
        }

        if (flag.str()=="-m")
            m = arg_number;
        
        else if (flag.str()=="-H")
        {
            H = uint(std::round(arg_number));
            NH = arg_number;
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            // g_g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-rH2")
        {
            rH_sqr = arg_number;
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            // g_g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-Rp2")
        {
            R_sqr = arg_number;
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            // g_g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-A")
            A = uint(std::round(arg_number));

        else if (flag.str()=="-s")
            g_seed = uint(std::round(arg_number));

        else if (flag.str()=="--add-to-seed")
            g_seed += uint(std::round(arg_number));

        else if (flag.str()=="--g2mu02-factor")
        {
            g_g2mu02 *= arg_number;
            g_g2mu02_config_factor = arg_number;
    #ifdef _G2MU02
            std::cerr << "ERROR: --g2mu02-factor flag was given. This is probably not what you want as this has been compiled with G2MU02=1." << std::endl;
            exit(25);
    #endif
        }
        else if (flag.str()=="-t" || flag.str()=="--threads")
        {
            g_num_threads = uint(std::round(arg_number));
    #ifndef _PC2
            GBWModel::G_ip.set_num_threads(g_num_threads);
    #endif
        }
        else if (flag.str()=="-Q")
            Q = arg_number;

        else if (flag.str() == "-Delta")
        {
            g_Delta_single = arg_number;
            g_Delta_single_set = true;
        }
        else
        {
            std::cerr << error_message << std::endl;
            exit(23);
        }
    }
    import_interp_data_by_params();

_TEST_LOG("Parameters set")
}


void print_infos (std::ofstream& out, uint seed)
{
    out << "#m=" << m << std::endl;
    out << "#A=" << A << std::endl;
    out << "#H=" << H << std::endl;
    out << "#rH2=" << rH_sqr << std::endl;
    out << "#Rp2=" << R_sqr << std::endl;
    out << "#RC2=" << RC_sqr << std::endl;
    out << "#Seed=" << seed << std::endl;
}


void print_infos (std::ofstream& out, uint seed, const HotspotNucleus& nucleus)
{
    print_infos(out, seed);

    for (uint n=0, n_max=nucleus.get_atomic_num(); n<n_max; n++)
    {
        out << "#Nucleon " << n << std::endl;
        for (uint i=0, i_max=nucleus.get_num_hotspots_per_nucleon(); i<i_max; i++)
        {
            HotspotPos hotspot_pos = *nucleus.get_hotspot_pos(n, i);

            out << "#" << hotspot_pos.x << " " << hotspot_pos.y << std::endl;
        }
    }
}


uint get_unique_seed()
{
    std::thread::id thread_id_temp = std::this_thread::get_id();
    std::hash<std::thread::id> hash;

    uint seed = hash(thread_id_temp) * GET_PROCESS_ID() / time(0);
    if (seed==0)
        seed = hash(thread_id_temp)*GET_PROCESS_ID() - time(0);

    return seed;
}


std::string get_default_filepath_from_parameters()
{
    std::string filepath = "data/samples/";
#ifdef _G2MU02
    filepath += "g2mu02/";
#endif
    filepath += (char)quark_config;

#ifndef _G2MU02
    filepath += (g_g2mu02_config_factor >= 1.0) ? std::to_string( int(std::round(g_g2mu02_config_factor*10.0)) ) : "0"+std::to_string( int(std::round(g_g2mu02_config_factor*10.0)) );
#endif
    filepath += 
#ifndef _DILUTE
        "/de/";
#else
        "/di/";
#endif

    filepath += std::to_string(g_seed);

    return filepath;
}


double sin_zeros (uint n)
{
    return double(n)*PI;
}


double cos_zeros (uint n)
{
    return (-0.5+double(n))*PI;
}