#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../include/IntegrationRoutines.hpp"
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


double sqr (double x)
{
    return x*x;
}


double bessel_K_safe (int n, double x)
{
    if (x==0)
        x = 1.0e-20;
    return gsl_sf_bessel_Kn(n,x);
}


void set_import_filepath_by_parameters (std::string& filepath, const DataGenerationConfig* config)
{
    std::string m_path = std::to_string( int(m*100) );
    std::string rH_sqr_path = (rH_sqr >= 1.0) ?  std::to_string( int(std::round(rH_sqr*100)) ) : "0"+std::to_string( int(std::round(rH_sqr*100)) );

    filepath = "InterpolatorData/";

    filepath += "G_rH2_"+rH_sqr_path+"_m_0"+m_path+".dat";
#ifndef _QUIET
    std::ifstream file_check (filepath);

    if (!config)
        return;

    if (!file_check)
    {
        file_check.close();
        std::string answer;

        std::cout << "No data for this combination of m(=" << m << ") and rH2(=" << rH_sqr << "). Do you want to generate that data set? (y/n)\n(The path to the new file would be " << filepath << ". File size would be (at most) " << (config->nx+config->ny+config->nz+2+config->nx*config->ny*config->nz)*17.0/1.0e6 << "MB.)" << std::endl;
        bool input_accepted = false;
        while (!input_accepted)
        {
            std::cin >> answer;
            if (answer=="y")
            {
                input_accepted = true;

                std::cout << "Generating new data set..." << std::endl;
                GBWModel::G_ip.generate_data(GBWModel::G_wrapper, config, true);
                GBWModel::G_ip.export_data(filepath);
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
    }
    else 
    {
        return;
        file_check.close();
        std::cout << "File with name " << filepath << " for m=" << m << " found. Do you want to make that the active filepath? (y/n)" << std::endl;

        std::string answer;

        std::cin >> answer;

        if (answer=="y")
            return;
        else
            exit(22);
    }
#endif
}


void set_parameters (int argc, char** argv)
{
    const std::string error_message = "Invalid option. Valid flags are\n"
                "[-p] (progress monitor)\n"
                "[-m <gluon mass>]\n"
                "[-A <atomic number>]\n"
                "[-H <number of hotspots per nucleon>]\n"
                "[-rH2 <hotspot radius square>]\n"
                "[-Rp2 <nucleon radius square>]\n"
                "[-o <output filepath>]";

    for (int i=1; i<argc; i+=2)
    {
        std::istringstream flag(argv[i]);

        if (flag.str()=="-p")
        {
            progress_monitor_global = true;
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
            filepath_global = arg_string;
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
            //g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-rH2")
        {
            rH_sqr = arg_number;
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            //g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-Rp2")
        {
            R_sqr = arg_number;
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            //g2mu02 = g2mu02_factor*RC_sqr/NH;
        }
        else if (flag.str()=="-A")
            A = uint(std::round(arg_number));
        
        else
        {
            std::cerr << error_message << std::endl;
            exit(23);
        }
    }
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

    uint seed = hash(thread_id_temp)*getpid()/time(0);
    if (seed==0)
        seed = hash(thread_id_temp)*getpid()-time(0);

    return seed;
}