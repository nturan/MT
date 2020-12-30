#pragma once

/******************************************************************************/
/* This file includes some helper functions, structures, declarations for     */
/* fortran routines etc.                                                      */
/******************************************************************************/


#include <cmath>     /* some common mathematical functions                    */
#include <iostream>  /* input-output stream                                   */
#include <fstream>   /* file stream                                           */
#include <string>    
#include <gsl/gsl_sf_dilog.h>  
#include <vector>
#include <complex>


using namespace std;

/******************************************************************************/
/* This structure was implemented for easy readability.                       */
struct divTerm {
    double eps_2;   /* eps^-2 pole                                              */
    double eps_1;   /* eps^-1 pole                                              */
    double eps0;    /* finite term                                              */
    double eps1;    /* terms with factor eps^1                                  */
    double eps2;    /* terms with factor eps^2                                  */

  /* hier we define several operator overloadings                               */

    divTerm operator+ (const divTerm& rhs) const {
        /* divTerm + divTerm                                                          */
        divTerm res = { eps_2 + rhs.eps_2,
                       eps_1 + rhs.eps_1,
                       eps0 + rhs.eps0,
                       eps1 + rhs.eps1,
                       eps2 + rhs.eps2 };
        return res;
    }

    divTerm operator- (const divTerm& rhs) const {
        /* divTerm - divTerm                                                          */
        divTerm res = { eps_2 - rhs.eps_2,
                       eps_1 - rhs.eps_1,
                        eps0 - rhs.eps0,
                        eps1 - rhs.eps1,
                        eps2 - rhs.eps2 };
        return res;
    }

    divTerm operator* (const divTerm& rhs) const {
        /* divTerm * divTerm                                                          */
        divTerm res = {
    eps_2 * rhs.eps0 + eps_1 * rhs.eps_1 + eps0 * rhs.eps_2,
    eps_2 * rhs.eps1 + eps_1 * rhs.eps0 + eps0 * rhs.eps_1 + eps1 * rhs.eps_2,
    eps_2 * rhs.eps2 + eps_1 * rhs.eps1 + eps0 * rhs.eps0 + eps1 * rhs.eps_1 + eps2 * rhs.eps_2,
    eps_1 * rhs.eps2 + eps0 * rhs.eps1 + eps1 * rhs.eps0 + eps2 * rhs.eps_1,
     eps0 * rhs.eps2 + eps1 * rhs.eps1 + eps2 * rhs.eps0 };
        return res;
    }

    divTerm operator* (const double& rhs) const {
        /* divTerm * double                                                           */
        divTerm res = { eps_2 * rhs,
                       eps_1 * rhs,
                        eps0 * rhs,
                        eps1 * rhs,
                        eps2 * rhs };
        return res;
    }

    divTerm operator+ (const double& rhs) const {
        /* divTerm * double                                                           */
        divTerm res = { eps_2,
                       eps_1,
                        eps0 + rhs,
                        eps1,
                        eps2 };
        return res;
    }

    divTerm operator- (const double& rhs) const {
        /* divTerm * double                                                           */
        divTerm res = { eps_2,
                       eps_1,
                        eps0 - rhs,
                        eps1,
                        eps2 };
        return res;
    }

    divTerm operator/ (const double& rhs) const {
        /* divTerm / double                                                           */
        divTerm res = { eps_2 / rhs,
                       eps_1 / rhs,
                        eps0 / rhs,
                        eps1 / rhs,
                        eps2 / rhs };
        return res;
    }

    divTerm operator+= (const divTerm rhs) const {
        return *this + rhs;
    }


    divTerm operator-= (const divTerm rhs) const {
        return *this - rhs;
    }

    friend divTerm operator* (const double& lhs, const divTerm rhs) {
        /* double * divTerm                                                           */
        divTerm res = { rhs.eps_2 * lhs,
                       rhs.eps_1 * lhs,
                        rhs.eps0 * lhs,
                        rhs.eps1 * lhs,
                        rhs.eps2 * lhs };
        return res;
    }

    friend divTerm operator+ (const double& lhs, const divTerm rhs) {
        divTerm res = { rhs.eps_2,
                        rhs.eps_1,
                         rhs.eps0 + lhs,
                         rhs.eps1,
                         rhs.eps2 };
        return res;
    }

    friend divTerm operator- (const double& lhs, const divTerm rhs) {
        divTerm res = { -rhs.eps_2,
                           -rhs.eps_1,
                       lhs - rhs.eps0,
                           -rhs.eps1,
                           -rhs.eps2 };
        return res;
    }

    divTerm xEps() {
        /* This function multiplies divTerm with eps.  Practivally it moves           */
        /* the terms to the right.                                                    */
        return (divTerm) { 0.0, eps_2, eps_1, eps0, eps1 };
    }


    divTerm qEps() {
        /* This function divides divTerm with eps.  Practivally it moves              */
        /* the terms to the left.                                                     */
        return (divTerm) { eps_1, eps0, eps1, eps2, 0.0 };
    }
};


/******************************************************************************/
/* The following functions, methods are implemented for convenience, easy     */
/* readability.                                                               */
double ln(const double x);            /* natural logarithm                  */
double ln2(const double x);           /* square of the natural logarithm    */
double Li2(const double x);           /* dilogarithm                        */
/* functions to calculate phase space parameters                              */
double calcInvariantMass(const double vectorP[4]);
double calcTransversalMomentum(const double vectorP[4]);
double calcRapidity(const double vectorP[4]);
double calcPseudoRapidity(const double vectorP[4]);
double calcAzimuth(const double vectorP[4]);
double calcEnergy(const double vectorP[4]);
double dotp(const double p1[4], const double p2[4]); /* scalar 4 product    */
double dot3p(const double p1[4], const double p2[4]);/* scalar 3 product    */
double magnitude3(const double p[4]); /* magnitude of the vector components */
/* cosine between vector components of two 4 vectors                          */
double cos3p(const double p1[4], const double p2[4]);
/* lorentz boost along z axis                                                 */
void boostZ(const double pOld[4], double beta, double pNew[4]);
/* print 4 vector                                                             */
string printP(const double p[4]);
