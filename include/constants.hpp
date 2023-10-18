#pragma once


#include <math.h>
#include "../include/utilities.hpp"

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

const double m_Q_c = 1.275; //in GeV
const double A_c = 0.211; //in GeV3/2

const double e_Q = e_c;
const double A_Q = A_c; // in GeV3/2
//const double epsilon = 1.0; // in GeV  // not used, see epsilonFunc

const double sqrt_2m_c_Nc = std::sqrt(2.0*m_Q_c*Nc); // in GeV

// Unit conversion factors
const double hbarc = 0.1973; //GeV fm

static const double fmToGeVm1 = 1.0/hbarc;
static const double GeVTofmm1 = 1.0/hbarc;
static const double GeVm1Tofm = hbarc;
static const double fmm1ToGeV = hbarc;
static const double fm2TonB = 1.0e7;
static const double nBTofm2 = 1.0e-7;

// not exactly constant constants
inline double Nq = 1.0;
inline double m = 0.22;

inline double rH_sqr = 0.7;
inline double R_sqr = 3.3;

inline double RC_sqr = rH_sqr + (Nq-1)/Nq*R_sqr;

inline double BG = 4.0; //GeVm2 https://physics.nist.gov/cgi-bin/cuu/Value?rp
inline double sigma0 = 2.0*PI*BG; //GeVm2

inline uint A = 0;

// Integration ranges
const double B_RANGE_FACTOR = 15.0;
const double R_RANGE_FACTOR = 20.0;

const double B_MAX = std::sqrt(B_RANGE_FACTOR*2.0*BG);
const double R_MAX = R_RANGE_FACTOR/m_Q_c;

const double g2mu02_demirci = std::sqrt(43.22);
const double g2mu02 = sigma0;