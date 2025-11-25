#include "operations_vectors.h"
#include <cmath>
#include <iostream>

inline int idx(int i, int j, int N) {
    return i * N + j; //incidentalment al programa ho fare servir amb N, no N+2
}                     //perque tot i que realment la matriu es M+2 x N+2, el vector es N*M  

double getSafeMat(const std::vector<std::vector<double>>& matrix, int i, int j, int M, int N) {
    if (i < 0 || i > M ) {
        return 0.0;  // Treat as zero if out-of-bounds
    }
    if(j < 0 || j > N){
        return 0.0;  // Treat as zero if out-of-bounds
    }
    return matrix[i][j];
}

double getSafeVec(const std::vector<double>& vec, int i, int L) {
    if (i < 0 || i > L-1 ) {
        return 0.0;  // Treat as zero if out-of-bounds
    }
    return vec[i];
}

 
std::vector<double> matVecMult2DGrid(const std::vector<double> &vec, std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS
    , std::vector<double> &CoefW, std::vector<double> &CoefP,int M, int N){

    double size = vec.size();
    if(size!=CoefN.size() || size!=CoefE.size() || size!=CoefS.size() || size!=CoefW.size() || size!=CoefP.size()){
        std::cout << "Error, to do the matrix products with the 5 coefficient vectors they must all be the same size" << std::endl;
        return {};
    }

    std::vector<double> output(size,0);

    for (int i = 1; i < M+1; i++) {
        for (int j = 1; j < N+1; j++) {
            int k = idx(i-1, j-1, N); //la matriu de pressio realment es M+2 x N+2 
            // Diagonal 
            output[k] += CoefP[k] * vec[k];
            // Neighbours
            if(i > 1) output[k] += CoefN[k]*vec[idx(i-2, j-1, N)];  // N, la q=dy/dx hauria d'estar ja al vector dels coeficients
            if(i < M) output[k] += CoefS[k]*vec[idx(i, j-1, N)];  // S
            if(j > 1) output[k] += CoefW[k]*vec[idx(i-1, j-2, N)];  // W
            if(j < N) output[k] += CoefE[k]*vec[idx(i-1, j, N)];  // E
        }
    }

    return output;
}

double dotP(const std::vector<double>& vec1, const std::vector<double>& vec2){

    if(vec1.size()!=vec2.size()){
        std::cout << "Error, dot products need same size vectors." << std::endl;
        return 1;
    }
    double sum = 0.0;
    for (int i = 0; i < vec1.size(); i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}        

double maxAbsValue(const std::vector<double>& vec) {
    double maxVal = 0.0;  // Initialize to 0 (or use std::numeric_limits<double>::lowest() for a more general approach)
    
    for (double val : vec) {
        maxVal = std::max(maxVal, std::abs(val));  // Compare current value with max absolute value so far
    }
    
    return maxVal;
}


std::vector<double> CombLin(const std::vector<double>& vec1, const std::vector<double>& vec2, double a, double b) { //a*vec1 + b*vec2
    // Ensure both vectors are of the same size
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: Vectors must be of the same size!" << std::endl;
        return {};  // Return empty vector on error
    }

    std::vector<double> result(vec1.size());  // Initialize result vector with the same size

    // Perform element-wise operation (either addition or subtraction)
    for(int i = 0; i < vec1.size();i++){
        result[i] = a*vec1[i] + b*vec2[i];
    }  
    return result;
}

void getCoefs(std::vector<double> &CoefN, std::vector<double> &CoefE, std::vector<double> &CoefS, std::vector<double> &CoefW, std::vector<double> &CoefP
    , int M, int N, double q){

    double size = CoefN.size();
    if(size!=CoefE.size() || size!=CoefS.size() || size!=CoefW.size() || size!=CoefP.size()){
        std::cout << "Error, to get the 5 coefficient vectors they must all be the same size" << std::endl;
        return;
    }

    for(int i = 1; i < M+1; i++){
        for(int j = 1; j < N+1; j++){
            if(i > 1){ CoefN[idx(i-1,j-1,N)] = 1/q; }
            if(j < N){ CoefE[idx(i-1,j-1,N)] = q; }
            if(i < M){ CoefS[idx(i-1,j-1,N)] = 1/q; }
            if(j > 1){ CoefW[idx(i-1,j-1,N)] = q; }
            CoefP[idx(i-1,j-1,N)] = 0 - getSafeVec(CoefN,idx(i-1,j-1,N),N*M) 
                                  - getSafeVec(CoefE,idx(i-1,j-1,N),N*M) 
                                  - getSafeVec(CoefS,idx(i-1,j-1,N),N*M) 
                                  - getSafeVec(CoefW,idx(i-1,j-1,N),M*N);
            //This relies on the CoefN,CoefE,..., vectors being initialized at 0 previously
        }
    }
} 