#ifndef MAINFX_H
#define MAINFX_H

#include <vector>

void PredUfill(double M, double N, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &RU, double rho, double dx, double dy, double dt, double mu);
void PredVfill(double M, double N, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &PredV
    , std::vector<std::vector<double>> &RV, std::vector<std::vector<double>> &Temp, std::vector<std::vector<double>> &nuTemp 
    , double dt, double beta, double g, double rho, double dx, double dy, double mu);
void ConjugateGradient(std::vector<double> &PressureVec, std::vector<double> &R, std::vector<double> &ct, double M, double N, double epsilon, double maxiter
    , std::vector<double> &P, std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP,int k);
void ErrorPrint(std::vector<std::vector<double>> Pressure, double M, double N, double rho, double dx, double dt, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &PredV, double &Error, double q, int k, double &StationaryError);
void TempChange(std::vector<std::vector<double>> &Temp, double N, double M, std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &MatU
    , std::vector<std::vector<double>> &nuTemp, double dt, double lambda, double dy, double dx, double rho, double cp
    , double &StationaryError);
void NewVelocity(std::vector<std::vector<double>> &MatV, std::vector<std::vector<double>> &MatU, std::vector<std::vector<double>> &PredU
    , std::vector<std::vector<double>> &PredV, double &Vmax, double &vmax, double &Umax, double &umax, int k, double dt, int M
    , int N, double rho, std::vector<std::vector<double>> &Pressure, double L, int TSteps, double &StationaryError
    , std::vector<std::vector<double>>& MatPreU, std::vector<std::vector<double>>& MatPreV);
void GaussSeidel(int M, int N, double dy, double dx, double rho, double dt 
    , std::vector<std::vector<double>> &Pressure, std::vector<std::vector<double>>& PredU
    , std::vector<std::vector<double>>& PredV, int maxiter, double Omega, double epsilon
    , double Error, bool B,std::vector<double>& CoefP, std::vector<double>& CoefN
    , std::vector<double>& CoefE, std::vector<double>& CoefS, std::vector<double>& CoefW, int k);

inline int idx(int i, int j, int N){
    return i*N + j;
}

std::vector<double> matVecMult2DGrid(const std::vector<double> &vec, std::vector<double> &CoefN
    , std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP,int M, int N);
double dotP(const std::vector<double>& vec1, const std::vector<double>& vec2);
double maxAbsValue(const std::vector<double>& vec);
std::vector<double> CombLin(const std::vector<double>& vec1, const std::vector<double>& vec2, double a, double b);
void getCoefs(std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP, int M, int N, double q);


double avg(double phiP, double phi);
double convx(double rho, double dx, double un, double ue, double us, double uw, double vn, double vs, double q);
double convy(double rho, double dx, double vn, double ve, double vs, double vw, double ue, double uw, double q);
double diffx(double uN, double uE, double uS, double uW, double uP, double mu, double q);
double diffy(double vN, double vE, double vS, double vW, double vP, double mu, double q);
double Ppcoef(int N, int i, int j, double q, int M);
double Tpcoef(int N, int i, int j, double q, int M);
double getSafe(const std::vector<std::vector<double>>& matrix, int i, int j, int M, int N);
double Ppcoefreborn(int N, int i, int j, double q, int M);

//proper functions actually used in main are at the top
#endif