#include "helper.h"


using namespace std;



double calcInvariantMass(const double vectorP[4]) {
    return sqrt(dotp(vectorP, vectorP));
}

double calcTransversalMomentum(const double vectorP[4]) {
    double  pT = sqrt(pow(vectorP[1], 2) + pow(vectorP[2], 2));

    if (isnan(pT)) {
        return 0.0;
    }
    else {
        return pT;
    }
}

double calcRapidity(const double vectorP[4]) {
    double pZ = vectorP[3];
    double E = vectorP[0];

    double rapidity = 0.5 * ln((E + pZ) / (E - pZ));
    if (isnan(rapidity)) {
        cout << rapidity << endl;
        cout << (E + pZ) / (E - pZ) << endl;
        exit(1);
        return 100.0;
    }
    else {
        return rapidity;
    }
}

double calcPseudoRapidity(const double vectorP[4]) {
    double theta = acos(vectorP[3] / magnitude3(vectorP));


    double eta = -ln(tan(theta / 2.0));

    if (isnan(eta)) {
        return 0.0;
    }
    else {
        return eta;
    }

}

double calcAzimuth(const double vectorP[4]) {
    double phi = atan2(vectorP[2], vectorP[1]);

    if (isnan(phi)) {
        return 0.0;
    }
    else {
        return phi;
    }
}


double calcEnergy(const double vectorP[4]) {
    return vectorP[0];
}

void setMomenta(double p[4], double E,
    double px, double py, double pz) {
    p[0] = E;
    p[1] = px;
    p[2] = py;
    p[3] = pz;
}

void boostZ(const double pOld[4], double beta, double pNew[4]) {


    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    pNew[0] = gamma * pOld[0] - gamma * beta * pOld[3];
    pNew[1] = pOld[1];
    pNew[2] = pOld[2];
    pNew[3] = gamma * pOld[3] - gamma * beta * pOld[0];
}


double ln(const double x) { return log(x); }

double ln2(const double x) { return log(x) * log(x); }

double Li2(const double x) { return gsl_sf_dilog(x); }

double dotp(const double p1[4], const double p2[4]) {
    /**
       calculates scalar product of lorentz 4-vectors.
    */
    return p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2] - p1[3] * p2[3];
}


double dot3p(const double p1[4], const double p2[4]) {
    /**
       calculates scalar product of vector components of
       lorentz 4 vectors.
    */
    return p1[1] * p2[1] + p1[2] * p2[2] + p1[3] * p2[3];

}
double magnitude(const double p[4]) {
    /**
       returns magnitude of 4-vector.
    */
    return sqrt(dotp(p, p));
}

double magnitude3(const double p[4]) {
    /**
       returns magnitude of vector part of the 4-vector.
    */
    return sqrt(dot3p(p, p));
}

double cos3p(const double p1[4], const double p2[4]) {
    return dot3p(p1, p2) / magnitude3(p1) / magnitude3(p2);
}

void BoostInGeneratlDirection(const double p_old[4], double beta[4], double p_new[4]) {
    double gamma = beta[0];
    double beta_p_old = dot3p(p_old, beta);
    double beta_squared = dot3p(beta, beta);
    p_new[0] =                      gamma * p_old[0] - gamma * beta_p_old;
    p_new[1] = p_old[1] - gamma * beta[1] * p_old[0] + (gamma - 1) * beta[1] * beta_p_old / beta_squared;
    p_new[2] = p_old[2] - gamma * beta[2] * p_old[0] + (gamma - 1) * beta[2] * beta_p_old / beta_squared;
    p_new[3] = p_old[3] - gamma * beta[3] * p_old[0] + (gamma - 1) * beta[3] * beta_p_old / beta_squared;
}

string printP(const double p[4]) {
    return "(" + to_string(p[0]) + ", " + to_string(p[1]) + ", "
        + to_string(p[2]) + ", " + to_string(p[3]) + ")";
}

string logP(const double p[4]) {
    return to_string(p[0]) + "," + to_string(p[1]) + ","
        + to_string(p[2]) + "," + to_string(p[3]);
}
