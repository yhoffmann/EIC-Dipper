#include "../include/utilities.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../include/IntegrationRoutines.hpp"
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>


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

        std::cout << "No data for this combination of m(=" << m << ") and rH2(=" << rH_sqr << "). Do you want to generate that data set? (y/n)\n(The path to the new file would be " << filepath << ". File size would be (at most) " << (config->n_x+config->n_y+config->n_z+2+config->n_x*config->n_y*config->n_z)*17.0/1.0e6 << "MB.)" << std::endl;
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
    std::string error_message = "Invalid option. Valid flags are\n"
                "[-pm] (progress monitor)\n"
                "[-m <gluon mass>]\n"
                "[-A <atomic number>]\n"
                "[-H <number of hotspots per nucleon>]\n"
                "[-rH2 <hotspot radius square>]\n"
                "[-Rp2 <nucleon radius square>]\n"
                "[-o <output filepath>]";

    for (int i=1; i<argc; i+=2)
    {
        if (!strcmp(argv[i], "-pm"))
        {
            progress_monitor_global = true;
            --i;
            continue;
        }
        else if (i+1==argc)
        {
            std::cerr << error_message << std::endl;
            exit(23);
        }
        else if (!strcmp(argv[i], "-o"))
        {
            filepath_global = std::string(argv[i+1]);
            continue;
        }

        if (std::atof(argv[i+1])==0.0 && std::atoi(argv[i+1])==0)
        {
            std::cerr << "Please enter a valid number other than 0" << std::endl;
            std::cerr << error_message << std::endl;
            exit(23);
        }

        if (!strcmp(argv[i], "-m"))
            m = std::atof(argv[i+1]);
        else if (!strcmp(argv[i], "-H"))
        {
            H = std::atoi(argv[i+1]);
            NH = double(H);
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            g2mu02 = g2mu02_factor*rH_sqr/NH;
        }
        else if (!strcmp(argv[i], "-rH2"))
        {
            rH_sqr = std::atof(argv[i+1]);
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
            g2mu02 = g2mu02_factor*rH_sqr/NH;
        }
        else if (!strcmp(argv[i], "-Rp2"))
        {
            R_sqr = std::atof(argv[i+1]);
            RC_sqr = rH_sqr + (NH-1.0)/NH*R_sqr;
        }
        else if (!strcmp(argv[i], "-A"))
            A = std::atoi(argv[i+1]);
        else
        {
            std::cerr << error_message << std::endl;
            exit(23);
        }
    }
}


void print_infos (std::ofstream& out)
{
    out << "#m=" << m << std::endl;
    out << "#A=" << A << std::endl;
    out << "#H=" << H << std::endl;
    out << "#NH=" << NH << std::endl;
    out << "#rH2=" << rH_sqr << std::endl;
    out << "#Rp2=" << R_sqr << std::endl;
    out << "#RC2=" << RC_sqr << std::endl;
}