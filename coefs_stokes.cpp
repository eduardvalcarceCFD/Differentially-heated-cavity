#include "coefs_stokes.h"
#include <cmath>
#include <iostream>

double avg(double phiP, double phi){
    return phiP*0.5 + phi*0.5;
}

double convx(double rho, double dx, double un, double ue, double us, double uw, double vn, double vs, double q){
    return rho*dx*(vn*un + ue*ue*q - vs*us - uw*uw*q);
}

double convy(double rho, double dx, double vn, double ve, double vs, double vw, double ue, double uw, double q){
    return rho*dx*(vn*vn + ue*ve*q - vs*vs - vw*uw*q);
}

double diffx(double uN, double uE, double uS, double uW, double uP, double mu, double q){
    return mu*((uN+uS)/q + (uE+uW)*q-2*(q+1/q)*uP);
}

double diffy(double vN, double vE, double vS, double vW, double vP, double mu, double q){
    return mu*((vN+vS)/q + (vE+vW)*q-2*(q+1/q)*vP);
}

double Ppcoef(int N, int i, int j, double q, int M){
    return 2*(q+1/q)-static_cast<double>(i==1||i==M)/q-static_cast<double>(j==1||j==N)*q; 
}

double Tpcoef(int N, int i, int j, double q, int M){
    return 2*(q+1/q)+static_cast<double>(i==1||i==M)/q+static_cast<double>(j==1||j==N)*q; 
}


double Ppcoefreborn(int N, int i, int j, double q, int M){
    return 2*(q+1/q)+(i==1||i==M)/q+(j==1||j==N)*q;
}
