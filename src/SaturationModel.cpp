#include "../include/SaturationModel.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/GBWModel.h"

//#include "../../../libs/EigenLib/Eigen/Eigen"
//#include "../../../libs/EigenLib/unsupported/Eigen/MatrixFunctions"


namespace SaturationModel {
    double D (double G) {
        return exp(G);
    }

    double DD (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
        double a = DDCorrelationMatrixElements::a(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        double b = DDCorrelationMatrixElements::b(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        double c = DDCorrelationMatrixElements::c(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        double d = DDCorrelationMatrixElements::d(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

        double Sqrt = std::sqrt(4.0*b*c+sqr(a-d));
        double factor = (a-d+2.0/Nc*b)/Sqrt;

        double ret = 0.5* ( (1.0+factor)*exp((a+d+Sqrt)/2.0) + (1.0-factor)*exp((a+d-Sqrt)/2.0) );

        //double ret = exp(aplusd/2.0) * ( cosh(0.5*Sqrt) + factor*sinh(0.5*Sqrt) );

        if (gsl_isnan(ret) && 1==0) {
            std::cerr << ret << " " << x1 << " " << x2 << " " << y1 << " " << y2 << " " << xb1 << " " << xb2 << " " << yb1 << " " << yb2 << std::endl;
            std::cout << 4.0*b*c+sqr(a-d) << " " << a+d << " " << factor;
            std::cerr << "DipoleDipole is nan. Aborting" << std::endl;
            exit(0);
        }

        return ret;
    }

    
    double dsigma_d2b (double x1, double x2, double y1, double y2) {
        return 2.0*(1.0 - exp(GBWModel::G(x1,x2,y1,y2)));
    }

    double dsigma_d2b_sqr (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
        using namespace DDCorrelationMatrixElements;

        double G_xy = GBWModel::G(x1,x2,y1,y2);
        double G_xbyb = GBWModel::G(xb1,xb2,yb1,yb2);
        
        return 4.0 * (1.0 -D(G_xy) -D(G_xbyb) +DD(x1,x2,y1,y2,xb1,xb2,yb1,yb2));
    }

    /*namespace DDCorrelationMatrixElements {
        using namespace GBWModel;

        double G_xy (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(x(b1,r1),x(b2,r2),y(b1,r1),y(b2,r2));
        }

        double G_xyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(x(b1,r1),x(b2,r2),y(bb1,rb1),y(bb2,rb2));
        }

        double G_xby (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(x(bb1,rb1),x(bb2,rb2),y(b1,r1),y(b2,r2));
        }

        double G_xbyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(x(bb1,rb1),x(bb2,rb2),y(bb1,rb1),y(bb2,rb2));
        }

        double G_xxb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(x(b1,r1),x(b2,r2),x(bb1,rb1),x(bb2,rb2));
        }

        double G_yyb (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G(y(b1,r1),y(b2,r2),y(bb1,rb1),y(bb2,rb2));
        }

        double T1 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G_xxb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_yyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xby(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
        }

        double T2 (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G_xy(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xbyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - G_xby(b1,b2,r1,r2,bb1,bb2,rb1,rb2);
        }


        double a (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G_xy(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_xbyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - (T1(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(Nc*Nc-1);
        }

        double b (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return (T1(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(2*CF);
        }

        double c (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {       
            return (T2(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(2*CF);
        }

        double d (double b1, double b2, double r1, double r2, double bb1, double bb2, double rb1, double rb2) {
            return G_xxb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) + G_yyb(b1,b2,r1,r2,bb1,bb2,rb1,rb2) - (T2(b1,b2,r1,r2,bb1,bb2,rb1,rb2))/(Nc*Nc-1);
        }
    }*/

    namespace DDCorrelationMatrixElements {
        using namespace GBWModel;

        double T_xy_ybxb (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double G_xxb = G(x1,x2,xb1,xb2);
            double G_yyb = G(y1,y2,yb1,yb2);
            double G_xyb = G(x1,x2,yb1,yb2);
            double G_ybx = G(y1,y2,xb1,xb2);

            return G_xxb+G_yyb-G_xyb-G_ybx;
        }

        double T_xxb_yby (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double G_xy = G(x1,x2,y1,y2);
            double G_xbyb = G(xb1,xb2,yb1,yb2);
            double G_xyb = G(x1,x2,yb1,yb2);
            double G_ybx = G(y1,y2,xb1,xb2);

            return G_xy+G_xbyb-G_xyb-G_ybx;
        }

        double a (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double G_xy = G(x1,x2,y1,y2);
            double G_ybxb = G(xb1,xb2,yb1,yb2);
            double T = T_xy_ybxb(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

            return G_xy+G_ybxb-T/(sqr(Nc)-1.0);
        }

        double b (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double T = T_xy_ybxb(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

            return T/(2.0*CF);
        }

        double c (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double T = T_xxb_yby(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

            return T/(2.0*CF);
        }

        double d (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
            double G_xxb = G(x1,x2,xb1,xb2);
            double G_yyb = G(y1,y2,yb1,yb2);
            double T = T_xxb_yby(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

            return G_xxb+G_yyb-T/(sqr(Nc)-1.0);
        }
    }
/*
    double DDEigen (double x1, double x2, double y1, double y2, double xb1, double xb2, double yb1, double yb2) {
        Eigen::MatrixXd M(2,2);

        M(0,0) = DDCorrelationMatrixElements::a(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        M(0,1) = DDCorrelationMatrixElements::b(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        M(1,0) = DDCorrelationMatrixElements::c(x1,x2,y1,y2,xb1,xb2,yb1,yb2);
        M(1,1) = DDCorrelationMatrixElements::d(x1,x2,y1,y2,xb1,xb2,yb1,yb2);

        Eigen::VectorXd lVec(2); lVec(0)=1.0; 	lVec(1)=0.0;
        Eigen::VectorXd rVec(2); rVec(0)=Nc*Nc; rVec(1)=Nc;

        return lVec.dot(M.exp()*rVec);
    }
*/
}