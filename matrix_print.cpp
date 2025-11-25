#include "matrix_print.h"
#include <cmath>
#include <fstream>
#include <iostream>




void writeInnerMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    int rows = matrix.size();
    if (rows < 3) return; // No inner region
    int cols = matrix[0].size();
    if (cols < 3) return;

    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            file << matrix[i][j];
            if (j != cols - 2) file << " "; // space except after last column
        }
        file << "\n";
    }

    file.close();
}



void write_full_matrix_to_file(const std::vector<std::vector<double>>& matrix) {
    const std::string filename = "UV.txt";
    std::ofstream file(filename);

    int rows = matrix.size();
    if (rows == 0) return;
    int cols = matrix[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << matrix[i][j] << " ";
        }
        file << "\n";
    }

    file.close();
}

void write_full_matrix_to_fileTemp(const std::vector<std::vector<double>>& matrix) {
    const std::string filename = "Temp.txt";
    std::ofstream file(filename);

    int rows = matrix.size();
    if (rows == 0) return;
    int cols = matrix[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << matrix[i][j] << " ";
        }
        file << "\n";
    }

    file.close();
}

void write_centralcolumn(const std::vector<std::vector<double>>& matrix, int N, double dy) {
    //to print the central column to check the values as asked.
    // by our choice, (N+1)/2 should always be integer.
    const std::string filename = "outputcol.txt";
    std::ofstream file(filename);

    int rows = matrix.size();
    if (rows == 0) return;

    for (int i = 1; i < rows-1; ++i) {
        file <<  1-(2*i-1)*dy/2  << " " << matrix[i][floor((N+2)/2)] << " ";
        file << "\n";
    }

    file.close();
}

void write_centralfile(const std::vector<std::vector<double>>& matrix, int M, double dx) {
    //to print the central column to check the values as asked.
    // by our choice, (N+1)/2 should always be integer.
    const std::string filename = "outputrow.txt";
    std::ofstream file(filename);

    int columns = matrix[0].size();
    if (columns == 0) return;

    for (int i = 1; i < columns-1; ++i) {
        file <<  (2*i-1)*dx/2  << " " << matrix[floor((M+2)/2)][i] << " ";
        file << "\n";
    }

    file.close();
}