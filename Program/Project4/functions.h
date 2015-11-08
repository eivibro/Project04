#ifndef FUNCTIONS
#define FUNCTIONS
#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

void initializeSameValueMatrix(int **matrix, int row, int column, double value, double &M, double &E, double J);
void initializeRamdomSpinMatrix(int **matrix, int row, int column, double &M, double &E, double J, mt19937 &mt, uniform_int_distribution<int> &myDist);
void metropolis(int **spinMatrix, int matrixDimension, double J, double &E, double &M,
                mt19937 &mt, uniform_int_distribution<int> &intDist, uniform_real_distribution<double> &realDist, double exponentials[]);
void output(double temperature, int mcCycles, int L, double average[], ofstream &outFile);
#endif // FUNCTIONS

