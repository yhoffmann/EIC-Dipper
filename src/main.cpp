#include <string>
#include <random>
#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../include/include.h"
#include "../Interpolation3D/include/Interpolator3D.h" // TODO change include path to the file in this directory


double rng01() {
    return drand48();
}


int main (int argc, char** argv)
{
    double x1 = std::atof(argv[1]);
    double x2 = std::atof(argv[2]);
    double y1 = std::atof(argv[3]);
    double y2 = std::atof(argv[4]);

    std::string filepath = "InterpolatorData/G";
    set_import_filepath_by_m(filepath);


    //dc.x_grid_spacing = "log";
    //dc.y_grid_spacing = "log";
    //dc.z_grid_spacing = "linear";

    GBWModel::G_ip.import_data(filepath);
/*
    double counter = 0.0;
    while (true) {
        double bx = (5+counter)*2*(rng01()-0.5);
        double by = (5+counter)*2*(rng01()-0.5);
        double bbx = (5+counter)*2*(rng01()-0.5);
        double bby = (5+counter)*2*(rng01()-0.5);
        double rx = (5+counter)*2*(rng01()-0.5);
        double ry = (5+counter)*2*(rng01()-0.5);
        double rbx = (5+counter)*2*(rng01()-0.5);
        double rby = (5+counter)*2*(rng01()-0.5);

        double printthis = SaturationModel::DD(bx,by,bbx,bby,rx,ry,rbx,rby);

        //std::cout << bvalue << "\t" << bbvalue << "\t" << _a << "\t" << _b << "\t" << _c << "\t" << _d << "\t" << GBWModel::dsigma_dip_d2b(bx,by,rx,ry) << std::endl;
        if (printthis>1.0 || printthis<0.0)
        {
            std::cout << "val " << printthis << "\n" << bx << ", " << by << ", " << bbx << ", " << bby << "," << rx << ", " << ry << ", " << rbx << ", " << rby << std::endl;
        }
        
        //std::this_thread::sleep_for(std::chrono::milliseconds(150));
        //counter += 0.1;
    }
*/
    
/*
    const int imax = 1000;
    double results[imax][3];
    #pragma omp parallel for schedule(dynamic,1)
    for (int i=0; i<imax; i++)
    {
        double x = -10.0+20.0*double(i)/double(imax-1);
        results[i][0] = x;
        results[i][1] = GBWModel::G(x,x2,y1,y2);
        results[i][2] = GBWModel::G_by_integration(x,x2,y1,y2);
        std::cout << i << std::endl;
    }

    std::ofstream out;
    out.open("Data/G_interpolation_vs_integration.dat");
    if (!out.is_open()) { std::cout << "couldnt open" << std::endl; exit(0); }

    for (int i=0; i<imax; i++)
    {
        out << std::setprecision(10) << results[i][0] << " " << results[i][1] << " " << results[i][2] << std::endl;
    }

    out.close();
*/



    return 0;
}