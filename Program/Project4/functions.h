#ifndef FUNCTIONS
#define FUNCTIONS
#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
using namespace std;

void initializeSameValueMatrix(int **matrix, int row, int column, double value, double &M, double &E, double J);
void initializeRamdomSpinMatrix(int **matrix, int row, int column);
void metropolis(int **spinMatrix, int matrixDimension, double temperature, double J, double &E, double &M,
                mt19937 mt, uniform_int_distribution<int> intDist, uniform_real_distribution<double> realDist);
void output(double temperature, double E, ofstream &outFile);
#endif // FUNCTIONS

