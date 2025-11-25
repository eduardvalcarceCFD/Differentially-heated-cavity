#ifndef MATRIX_PRINT_H
#define MATRIX_PRINT_H

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

void writeInnerMatrixToFile(const std::vector<std::vector<double>>& matrix, const std::string& filename);
void write_full_matrix_to_file(const std::vector<std::vector<double>>& matrix);
void write_full_matrix_to_fileTemp(const std::vector<std::vector<double>>& matrix);
void write_centralcolumn(const std::vector<std::vector<double>>& matrix, int N, double dy);
void write_centralfile(const std::vector<std::vector<double>>& matrix, int M, double dx);

#endif
