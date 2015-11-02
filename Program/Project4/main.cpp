#include <iostream>
#include "lib.h"
#include "functions.h"
#include <unittest++/UnitTest++.h>
#include <random>

using namespace std;

int main()
{
    int L = 2;
    //Random generator engine and distributions
    mt19937 mt(4724);
    uniform_int_distribution<int> intDist(0,1);
    uniform_real_distribution<double> realDist(0,1);
    double E,M;
    E = M = 0;
    double startTemperature = 0.1;
    double endTemperature = 30;
    double tempStep = 1;
    double J = 1;
    int **spins;
    spins = (int **)matrix(L, L, sizeof(double));
    initializeSameValueMatrix(spins, L, L, 1, M, E, J);
    cout << "Initial energy:            " << E << endl;
    cout << "Initial magnetic momentum: " << M << endl;
    ofstream outfile;
    outfile.open("output.txt");
    for(double t = startTemperature; t <= endTemperature; t+=tempStep){
        for(int i = 0; i < 20; i++){
            metropolis(spins, L, t, J, E, M, mt, intDist, realDist);
        }
        cout << "Temperature = " << t << endl;
        cout << spins[0][0] << " " << spins[0][1] << endl;
        cout << spins[1][0] << " " << spins[1][1] << endl;
        output(t,E, outfile);

    }
    free_matrix((void **) spins);
    return UnitTest::RunAllTests();
}




