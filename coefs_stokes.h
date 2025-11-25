#ifndef COEFS_STOKES_H
#define COEFS_STOKES_H

#include <vector>

double avg(double phiP, double phi);
double convx(double rho, double dx, double un, double ue, double us, double uw, double vn, double vs, double q);
double convy(double rho, double dx, double vn, double ve, double vs, double vw, double ue, double uw, double q);
double diffx(double uN, double uE, double uS, double uW, double uP, double mu, double q);
double diffy(double vN, double vE, double vS, double vW, double vP, double mu, double q);
double Ppcoef(int N, int i, int j, double q, int M);
double Tpcoef(int N, int i, int j, double q, int M);
double Ppcoefreborn(int N, int i, int j, double q, int M);

#endif
