#include "../include/utilities.h"

#include <gsl/gsl_sf.h>
#include <iostream>
#include "../include/constants.h"
#include <fstream>
#include "../include/GBWModel.h"
#include "../include/IntegrationRoutines.h"

double bessel_K_safe (int n, double x) {
    if (x==0) {x = 1.0e-20;}
    return gsl_sf_bessel_Kn(n,x);
}


void set_import_filepath_by_m (std::string& filepath, DataGenerationConfig* config)
{
    std::string m_path;
    m_path = std::to_string(int(100*m));

    filepath = (m_path.size()<=2) ? filepath+"_m_0"+m_path+".dat" : filepath+"_m_"+m_path+".dat";

    std::ifstream file_check (filepath);

    if (!file_check)
    {
        file_check.close();
        std::string answer;

        std::cout << "No data for m=" << m << ". Do you want to generate that data set? (y/n)\n(The path to the new file would be " << filepath << ". File would have size ~" << config->n_x*config->n_y*config->n_z*61.0/1.0e6 << "MB.)" << std::endl;
        bool input_accepted = false;
        while (!input_accepted)
        {
            std::cin >> answer;
            if (answer=="y")
            {
                input_accepted = true;

                std::cout << "Generating new data set..." << std::endl;
                GBWModel::G_ip.generate_data(GBWModel::G_wrapper,config,true);
                GBWModel::G_ip.export_data(filepath);
            }
            else if (answer=="n")
            {
                input_accepted = true;

                std::cout << "Aborting" << std::endl;
                exit(0);
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
        bool accepted = false;
        while (!accepted)
        {
            std::cout << "File with name " << filepath << " for m=" << m << " found. Do you want to read data from that file? (y/n)" << std::endl;

            std::string answer;

            std::cin >> answer;

            if (answer=="y") accepted = true;
            else exit(0);
        }
    }
}