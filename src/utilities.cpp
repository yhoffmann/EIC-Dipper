#include "../include/utilities.h"

#include <gsl/gsl_sf.h>


double bessel_K_safe (int n, double x) {
    if (x==0) {x = 1.0e-20;}
    return gsl_sf_bessel_Kn(n,x);
}