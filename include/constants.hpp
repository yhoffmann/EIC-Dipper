#pragma once


#include <math.h>


// DEFINING GLOBAL PARAMETERS //

const double PI = M_PI;
const double alpha_em = 1.0/137.036; //unit 1
const double e = std::sqrt(4.0*PI*alpha_em); //unit 1
const double Nc = 3.0;
const double CF = (Nc*Nc-1.0)/(2.0*Nc);
const double Qs0 = 1.0; //GeV
const double A_c = 0.211; //in GeV3/2

const double e_c = 2.0/3.0;
const double e_b = -1.0/3.0;
const double e_t = 2.0/3.0;
const double m_c = 1.275; // GeV
const double m_b = 4.18; // GeV


// Unit conversion factors
const double hbarc = 0.1973; //GeV fm

const double fm_to_GeVm1 = 1.0/hbarc;
const double GeV_to_fmm1 = 1.0/hbarc;
const double GeVm1_to_fm = hbarc;
const double fmm1_to_GeV = hbarc;
const double fm2_to_nb = 1.0e7;
const double nb_to_fm2 = 1.0e-7;


// not exactly constant constants
inline double Q = std::sqrt(0.1);
inline double m_Q = m_c;
inline double e_Q = e_c;

inline double A_Q = A_c; // in GeV3/2

inline double NH = 1.0;
inline double m = 0.22;

inline double rH_sqr = 0.7;
inline double R_sqr = 3.3;

inline double RC_sqr = rH_sqr + (NH-1)/NH*R_sqr;

inline double sigma0 = 2.0*PI*RC_sqr;

inline uint A = 1;
inline uint H = 3;


// Integration ranges
const double B_RANGE_FACTOR = 15.0;
const double R_RANGE_FACTOR = 20.0;

inline double B_MAX = B_RANGE_FACTOR*std::sqrt(2.0*RC_sqr);
inline double R_MAX = R_RANGE_FACTOR/m_Q;

const double g2mu02_demirci = std::sqrt(43.22);
const double g2mu02_factor = g2mu02_demirci/(2.9)*3.0;
inline double g_g2mu02_config_factor = 1.0;
inline double g_g2mu02 = g2mu02_demirci;


// useful consts for speed
const double Ncsqrm1 = Nc*Nc-1.0;
const double Ncsqrm1_inverse = 1.0/Ncsqrm1;
const double twoNc = 2.0*Nc;
const double twoNc_inverse = 1.0/twoNc;
const double twoCF_inverse = 1.0/(2.0*CF);