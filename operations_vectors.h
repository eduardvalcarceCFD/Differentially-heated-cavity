#ifndef OPERATIONS_VECTORS_H
#define OPERATIONS_VECTORS_H

#include <vector>

inline int idx(int i, int j, int N);
double getSafe(const std::vector<std::vector<double>>& matrix, int i, int j, int M, int N);
std::vector<double> matVecMult2DGrid(const std::vector<double> &vec, std::vector<double> &CoefN
    , std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP,int M, int N);
double dotP(const std::vector<double>& vec1, const std::vector<double>& vec2);
double maxAbsValue(const std::vector<double>& vec);
std::vector<double> CombLin(const std::vector<double>& vec1, const std::vector<double>& vec2, double a, double b);
void getCoefs(std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP, int M, int N, double q);

#endif
