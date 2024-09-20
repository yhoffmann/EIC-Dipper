#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <unistd.h>
#include "../include/IntegrationRoutines.hpp"
#include "../include/Output.hpp"
#include "../include/Coherent.hpp"
#include "../include/Incoherent.hpp"
#include "../include/constants.hpp"
#include "../include/GBWModel.hpp"
#include "../external/Nucleus/include/HotspotNucleus.hpp"
#include "../external/thread-pool/include/ThreadPool.hpp"


namespace Output
{
    double get_default_Q()
    {
        return std::sqrt(0.1);
    }

    std::vector<double> get_default_Delta_vec()
    {
        return std::vector<double>{0.001, 0.01, 0.04, 0.08, 0.12, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.3, 3.0, 4.0}; //{0.001, 0.002, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.08, 0.09, 0.12, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.7, 3.0, 3.3, 3.7, 4.0}; //
    }

    std::vector<double> get_default_phi_vec()
    {
        std::vector<double> default_phi_vec;
        uint num_angles = 48;
        for (uint i=0; i<num_angles; i++)
            default_phi_vec.push_back(PI*double(i)/double(num_angles)); // no num_angles-1 because we do not want to reach the end // NOTE if the calculation of the coherent cross section ever changes to not use sin and cos anymore but something more complicated, this needs to be adjusted to 2pi again and also the writing of data to file needs to be changed (currently I am using that cos and sin are symmetric/antisymmetric to print 2 results for one calculation/angle)

        return default_phi_vec;
    }

    void dsigmadt (bool do_coherent, bool do_incoherent, std::string output_file = "")
    {
        double default_Q = get_default_Q();
        std::vector<double> default_Delta_vec = get_default_Delta_vec();
        std::vector<double> default_phi_vec = std::vector<double>{0.0};//get_default_phi_vec();//

        dsigmadt(do_coherent, do_incoherent, default_Q, default_Delta_vec, default_phi_vec, output_file);
    }

    void dsigmadt (bool do_coherent, bool do_incoherent, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec, std::string filepath)
    {
        if (phi_vec.size()==0)
            phi_vec = std::vector<double>{0.0};

        double coherent_results[Delta_vec.size()][phi_vec.size()];
        double coherent_avg[Delta_vec.size()];

        double incoherent[Delta_vec.size()];

        #pragma omp parallel for schedule(static, 1)
        for (uint i=0; i<Delta_vec.size(); i++)
        {
            coherent_avg[i] = 0.0;
            
            if (do_incoherent)
                incoherent[i] = Incoherent::dsigmadt(Q, Delta_vec[i]);

            for (uint j=0, jmax=phi_vec.size(); j<jmax; j++)
            {
                if (progress_monitor_global)
                    std::cout << Q << " " << Delta_vec[i] << " " << phi_vec[j] << std::endl;

                if (do_coherent)
                {
                    coherent_results[i][j] = Coherent::dsigmadt_test(Q, Delta_vec[i], phi_vec[j]);
                    coherent_avg[i] += coherent_results[i][j]/double(jmax);
                }
            }
        }

        std::ofstream out;

        if (filepath == "")
            filepath = "Data/dsigma_dt.dat";

        out.open(filepath);
        if (!out.is_open())
            exit(20);

        out << std::setprecision(10);
        out << "#Q, Delta, Coher, Incoher; " << std::endl;
        out << "#1, 2,           3,4" << std::endl;
        out.flush();

        for (uint i=0; i<Delta_vec.size(); i++)
            out << Delta_vec[i] << " " << Q << " " << coherent_avg[i] << " " << incoherent[i] << std::endl;
        out << std::endl;

        out.close();
        out.open(filepath + ".all");
        if (!out.is_open())
            exit(20);

        for (uint i=0; i<Delta_vec.size(); i++)
        {
            out << Delta_vec[i] << " " << Q << "   ";
            for (uint j=0; j<phi_vec.size(); j++)
                out << phi_vec[j] << " " << coherent_results[i][j] << " " << incoherent[i] << "   ";
            out << std::endl;
        }
        out << std::endl;

        out.close();
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed)
    {
        double default_Q = get_default_Q();
        std::vector<double> default_Delta_vec = get_default_Delta_vec();
        std::vector<double> default_phi_vec = get_default_phi_vec();

        dsigmadt_nucleus(atomic_num, num_hotspots, seed, default_Q, default_Delta_vec, default_phi_vec);
    }


    void dsigmadt_nucleus (uint atomic_num, uint num_hotspots, uint seed, double Q, std::vector<double> Delta_vec, std::vector<double> phi_vec)
    {
        import_interp_data_by_params(interpolator_filepath);

        double coherent_results_real[Delta_vec.size()*phi_vec.size()];
        double coherent_results_imag[Delta_vec.size()*phi_vec.size()];
        double incoherent_results_real[Delta_vec.size()];
        double incoherent_results_imag[Delta_vec.size()];

        if (seed==0)
            seed = get_unique_seed();

        std::mt19937 rng(seed);

        HotspotNucleus nucleus(atomic_num, num_hotspots, rng);

        ThreadPool pool(num_threads);

        for (uint Delta_index=0, Delta_size=Delta_vec.size(); Delta_index<Delta_size; ++Delta_index)
        {
            double Delta = Delta_vec[Delta_index];
            pool.enq_job([Delta_index, Q, Delta, &nucleus, &incoherent_results_real, &incoherent_results_imag] {
                incoherent_results_real[Delta_index] = Incoherent::Sampled::dsigmadt_single_event(Q, Delta, nucleus);
                incoherent_results_imag[Delta_index] = 0.0;
            });
        }

        for (uint Delta_index=0, Delta_size=Delta_vec.size(); Delta_index<Delta_size; ++Delta_index)
        {
            for (uint phi_index=0, phi_size=phi_vec.size(); phi_index<phi_size; ++phi_index)
            {
                double Delta = Delta_vec[Delta_index];
                double phi = phi_vec[phi_index];
                pool.enq_job([Delta_index, phi_index, Q, Delta, phi_size, phi, &nucleus, &coherent_results_real, &coherent_results_imag] {
                    auto [coh_real, coh_imag] = Coherent::Sampled::sqrt_dsigmadt_single_event(Q, Delta, phi, nucleus);
                    
                    coherent_results_real[Delta_index*phi_size + phi_index] = coh_real;
                    coherent_results_imag[Delta_index*phi_size + phi_index] = coh_imag;
                });
            }
        }

        pool.await();
        pool.stop();

        std::ofstream out(get_default_filepath_from_parameters()+"_Amplitude.dat");

        if (!out.is_open())
            exit(20);

        out << std::setprecision(10);
        print_infos(out, seed, nucleus);

        out << "##Delta,   Q,        A Co real,Co imag,  A2 Inco\n";

        for (uint Delta_index=0, Delta_size=Delta_vec.size(); Delta_index<Delta_size; ++Delta_index)
        {
            out << Delta_vec[Delta_index] << " " << Q << "   ";
            for (uint phi_index=0, phi_size=phi_vec.size(); phi_index<phi_size; ++phi_index)
                out << phi_vec[phi_index] << " " << coherent_results_real[Delta_index*phi_size + phi_index] << " " << coherent_results_imag[Delta_index*phi_size + phi_index] << " " << incoherent_results_real[Delta_index] << " " << incoherent_results_imag[Delta_index] << "   ";
            for (uint phi_index=0, phi_size=phi_vec.size(); phi_index<phi_size; ++phi_index)
                out << phi_vec[phi_index]+PI << " " << -coherent_results_real[Delta_index*phi_size + phi_index] << " " << coherent_results_imag[Delta_index*phi_size + phi_index] << " " << incoherent_results_real[Delta_index] << " " << incoherent_results_imag[Delta_index] << "   ";
            out << std::endl;
        }
        out << std::endl;

        out.close();
    }


    void dsigmadt_demirci(std::string filepath) 
    {
        dsigmadt_demirci(get_default_Q(), filepath);
    }


    void dsigmadt_demirci (double Q, std::string filepath)
    {
        const uint DELTA_VEC_SIZE = 50;
        std::vector<double> Delta_vec;
        for (uint i=0; i<DELTA_VEC_SIZE; ++i)
            Delta_vec.push_back(std::sqrt(25.0)*double(i)/double(DELTA_VEC_SIZE-1) + 0.001);

        double coherent[DELTA_VEC_SIZE];
        double color_fluc[DELTA_VEC_SIZE];
        double hotspot_fluc[DELTA_VEC_SIZE];

        #pragma omp parallel for
        for (uint i=0; i<DELTA_VEC_SIZE; ++i)
        {
            coherent[i] = Coherent::Demirci::dsigmadt(Q, Delta_vec[i]);
            color_fluc[i] = Incoherent::Demirci::color_fluctuations(Q, Delta_vec[i]);
            hotspot_fluc[i] = Incoherent::Demirci::hotspot_fluctuations(Q, Delta_vec[i]);

            if(progress_monitor_global)
                std::cout << i << std::endl;
        }

        std::ofstream out(filepath);
        if (!out.is_open())
            exit(20);

        out << "t Co Inco Color Hotspot" << std::endl;
        for (uint i=0; i<DELTA_VEC_SIZE; ++i)
            out << sqr(Delta_vec[i]) << " " << coherent[i] << " " << color_fluc[i]+hotspot_fluc[i] << " " << color_fluc[i] << " " << hotspot_fluc[i] << std::endl;

        out.close();
    }


    void G (uint num_points, std::string filepath)
    {
        uint num_angles = 8;
        double angles[num_angles];
        for (uint i=0; i<num_angles; ++i)
            angles[i] = double(i)/double(num_angles)*PI;

        double r = 1.0;

        double results[num_points][1+2*num_angles];

        #pragma omp parallel for schedule(dynamic,1)
        for (unsigned int i=0; i<num_points; ++i)
        {
            double x1, x2, y1, y2;
            double b = 100.0*double(i)/double(num_points-1);

            results[i][0] = b;
            for (uint j=0; j<num_angles; ++j)
            {
                double b1 = 0.0;
                double b2 = b;

                double r1 = r*cos(angles[j]);
                double r2 = r*sin(angles[j]);

                x1 = b1 + 0.5*r1;
                x2 = b2 + 0.5*r2;
                y1 = b1 - 0.5*r1;
                y2 = b2 - 0.5*r2;

                results[i][1+2*j] = GBWModel::G(x1, x2, y1, y2);
                results[i][2+2*j] = GBWModel::G_by_integration(x1, x2, y1, y2);
            }
            std::cout << i << std::endl;
        }

        std::ofstream out;
        out.open(filepath);
        if (!out.is_open())
        {
#ifndef _QUIET
            std::cout << "couldnt open" << std::endl;
#endif
            exit(20);
        }

        out << std::setprecision(10);

        out << "x ";
        for (uint i=0; i<num_angles; ++i)
        {
            out << "angle" << i << "interp angle" << i << "integration ";
        }
        out << std::endl;

        for (unsigned int i=0; i<num_points; i++)
        {
            out << results[i][0] << " ";
            for (uint j=0; j<num_angles; ++j)
            {
                out << results[i][1+2*j] << " " << results[i][2+2*j] << " ";
            }
            out << std::endl;
        }

        out.close();
    }


    void hotspot_nucleus_thickness_1d (uint atomic_num, uint num_hotspots_per_nucleon, uint num_samples, uint num_points, uint seed, std::string filepath)
    {
        if (seed==0)
            seed = get_unique_seed();

        std::mt19937 rng(seed);

        HotspotNucleus hn(atomic_num, num_hotspots_per_nucleon, rng);

        double* thickness = new(std::nothrow) double [num_points];
        double* x = new(std::nothrow) double [num_points];
        if (thickness==nullptr || x==nullptr)
            exit(24);

        double x_max = hn.get_mean_bulk_radius()+10.0*hn.get_mean_surface_diffusiveness(), x_min = -x_max;
        double inverse_x_divisor = 1.0/double(num_points-1);
        for (uint i=0; i<num_points; i++)
        {
            thickness[i] = 0.0;
            x[i] = x_min+(x_max-x_min)*double(i)*inverse_x_divisor;
        }
        
        for (uint i=0; i<num_samples; i++)
        {
            for (uint j=0; j<num_points; j++)
            {
                thickness[j] += hn.get_hotspot_thickness(x[j], 0.0);
            }
            hn.sample_nucleon_pos();
        }

        double inverse_thickness_divisor = 1.0/double(num_samples);
        for (uint i=0; i<num_points; i++)
            thickness[i] *= inverse_thickness_divisor;

        if (filepath=="")
            filepath = "Data/hotspot_nucleus_thickness_1d.dat";

        std::ofstream out(filepath);
        if (!out.is_open())
            exit(20);

        for (uint i=0; i<num_points; i++)
            out << x[i] << " " << thickness[i] << std::endl;

        out.close();

        delete[] thickness;
        delete[] x;
    }


    void hotspot_nucleus_thickness_avg (uint atomic_num, uint num_hotspots_per_nucleon, uint start_seed, uint num_events, std::string filepath)
    {
        const uint size_x = 2e2;
        const uint size_y = 2e2;

        const double xmin = -5.0;
        const double xmax = 5.0;
        const double ymin = -5.0;
        const double ymax = 5.0;

        const double inverse_num = 1.0/double(num_events);

        double* thickness = new double [size_x*size_y];
        double* thickness_stddev = new double [size_x*size_y]; 
        for (uint j=0; j<size_y; ++j)
            for (uint k=0; k<size_x; ++k)
            {
                thickness[j*size_x + k] = 0.0;
                thickness_stddev[j*size_x + k] = 0.0;
            }
        
        double* x = new double [size_x];
        double* y = new double [size_y];

        HotspotPos* pos = new HotspotPos [atomic_num*num_hotspots_per_nucleon*num_events];

        for (uint i=0; i<num_events; ++i)
        {
            uint seed = start_seed + i;

            std::mt19937 rng(seed);
            HotspotNucleus hn(atomic_num, num_hotspots_per_nucleon, rng);
            hn.set_hotspot_size(std::sqrt(rH_sqr));
            hn.set_nucleon_size(std::sqrt(R_sqr));
            for (uint n=0; n<atomic_num; ++n)
                for (uint h=0; h<num_hotspots_per_nucleon; ++h)
                    pos[i*atomic_num*num_hotspots_per_nucleon + n*num_hotspots_per_nucleon + h] = *hn.get_hotspot_pos(n, h);

            for (uint j=0; j<size_y; ++j)
            {
                y[j] = ymin+(ymax-ymin)*double(j)/double(size_y-1);
                for (uint k=0; k<size_x; ++k)
                {
                    x[k] = xmin+(xmax-xmin)*double(k)/double(size_x-1);

                    thickness[j*size_x + k] += hn.get_hotspot_thickness(x[k], y[j])*inverse_num;
                }
            }
        }


        for (uint i=0; i<num_events; ++i)
        {
            uint seed = start_seed + i;
            std::mt19937 rng(seed);

            HotspotNucleus hn(atomic_num, num_hotspots_per_nucleon, rng);
            hn.set_hotspot_size(std::sqrt(rH_sqr));
            hn.set_nucleon_size(std::sqrt(R_sqr));
            for (uint j=0; j<size_y; ++j)
            {
                y[j] = ymin+(ymax-ymin)*double(j)/double(size_y-1);
                for (uint k=0; k<size_x; ++k)
                {
                    x[k] = xmin+(xmax-xmin)*double(k)/double(size_x-1);

                    thickness_stddev[j*size_x + k] += sqr(hn.get_hotspot_thickness(x[k], y[j]) - thickness[j*size_x + k])*inverse_num;
                }
            }
        }

        for (uint j=0; j<size_y; ++j)
        {
            y[j] = ymin+(ymax-ymin)*double(j)/double(size_y-1);
            for (uint k=0; k<size_x; ++k)
            {
                x[k] = xmin+(xmax-xmin)*double(k)/double(size_x-1);

                thickness_stddev[j*size_x + k] = std::sqrt(thickness_stddev[j*size_x + k]);
            }
        }

        std::ofstream out(filepath);
        if (!out.is_open())
            exit(20);
        // out << "x y thickness" << std::endl;
        for (uint j=0; j<size_y; ++j)
            for (uint k=0; k<size_x; ++k)
                out << x[k] << " " << y[j] << " " << thickness[j*size_x + k] << std::endl;
        out.close();

        out.open(filepath + ".stddev");
        if (!out.is_open())
            exit(20);
        for (uint j=0; j<size_y; ++j)
            for (uint k=0; k<size_x; ++k)
                out << x[k] << " " << y[j] << " " << thickness_stddev[j*size_x + k] << std::endl;
        out.close();

        out.open(filepath + ".pos");
        if (!out.is_open())
            exit(20);
        out << "x y" << std::endl;
        for (uint i=0; i<num_events; ++i)
            for (uint n=0; n<atomic_num; ++n)
                for (uint h=0; h<num_hotspots_per_nucleon; ++h)
                    out << pos[i*atomic_num*num_hotspots_per_nucleon + n*num_hotspots_per_nucleon + h].x << " " << pos[i*atomic_num*num_hotspots_per_nucleon + n*num_hotspots_per_nucleon + h].y << std::endl;
        
        delete[] thickness;
        delete[] thickness_stddev;
        delete[] x;
        delete[] y;
        delete[] pos;
    }
}
