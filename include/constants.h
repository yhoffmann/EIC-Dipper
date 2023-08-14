#pragma once


#include <math.h>
#include "../include/utilities.h"

// DEFINING GLOBAL PARAMETERS //
const double PI = M_PI;
const double alpha_em = 1.0/137.036; //unit 1
const double e = std::sqrt(4.0*PI*alpha_em); //unit 1
const double e_c = 2.0/3.0;
const double e_t = 2.0/3.0;
const double e_b = -1.0/3.0;
const double Nc = 3.0;
const double CF = (Nc*Nc-1.0)/(2.0*Nc);
const double Qs0 = 1.0; //GeV
const double g4 = 1.0;
const double mu0 = 1.0;
const double m = 1.5;

const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV3/2

const double e_Q = e_c;
const double A_Q = A_c; // in GeV3/2
const double epsilon = 1.0; // in GeV  // not used, see epsilonFunc

const double sqrt_2m_c_Nc = std::sqrt(2.0*m_Q_c*Nc); // in GeV

// Unit conversion factors
const double hbarc = 0.1973; //GeV fm

const double fmToGeVm1 = 1.0/hbarc;
const double GeVTofmm1 = 1.0/hbarc;
const double GeVm1Tofm = hbarc;
const double fmm1ToGeV = hbarc;
const double fm2TonB = 1.0e7;
const double nBTofm2 = 1.0e-7;

const double BG = 4.0; //GeVm2 https://physics.nist.gov/cgi-bin/cuu/Value?rp
//const double BGinfm2 = BG*GeVm1Tofm*GeVm1Tofm;
const double sigma0 = 2.0*PI*BG; //GeVm2

// Integration ranges
const double B_RANGE_FACTOR = 15.0;
const double R_RANGE_FACTOR = 20.0;

const double BMAX = std::sqrt(B_RANGE_FACTOR*2.0*BG);