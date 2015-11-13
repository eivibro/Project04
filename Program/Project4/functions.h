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
void initializeRamdomSpinMatrix(int **matrix, int row, int column, double &M, int spin_indices[], double &E, double J, mt19937 &mt, uniform_int_distribution<int> &myDist);
void metropolis(int **spinMatrix, int matrixDimension, double J, int spin_indices[], int &accepted_configurations, double &E, double &M,
                mt19937 &mt, uniform_int_distribution<int> &intDist, uniform_real_distribution<double> &realDist, double exponentials[]);
void output(double temperature, int mcCycles, int L, double average[], ofstream &outFile);
void outputc(double temperature, int mcCycles, int L, double average[], int accepted_configurations, ofstream &outFile);
//void spin_indices
#endif // FUNCTIONS

