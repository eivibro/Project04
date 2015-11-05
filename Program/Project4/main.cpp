#include <iostream>
#include "lib.h"
#include "functions.h"
#include <unittest++/UnitTest++.h>
#include <random>

using namespace std;

int main()
{
    int L = 2;
    int mcCycles = 100;
    //Random generator engine and distributions
    mt19937 mt(4724);
    uniform_int_distribution<int> intDist(0,1);
    uniform_real_distribution<double> realDist(0,1);
    double E,M;
    E = M = 0;
    double startTemperature = 1.0;
    double endTemperature = 20;
    double tempStep = 0.1;
    double J = 1;
    int **spins;
    spins = (int **)matrix(L, L, sizeof(double));
    initializeSameValueMatrix(spins, L, L, 1, M, E, J);
    cout << "Initial energy:            " << E << endl;
    cout << "Initial magnetic momentum: " << M << endl;
    ofstream outfile;
    outfile.open("output.txt");
    double averageE = 0;
    for(double t = startTemperature; t <= endTemperature; t+=tempStep){
        averageE = 0;
        for(int i = 0; i < mcCycles; i++){
            metropolis(spins, L, t, J, E, M, mt, intDist, realDist);
            averageE += E;

        }
        cout << "Temperature = " << t << endl;
        output(t,averageE/mcCycles, outfile);

    }
    outfile.close();
    free_matrix((void **) spins);
    return 0;
    //return UnitTest::RunAllTests();
}




