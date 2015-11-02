#include "functions.h"
//Initialize all elements of the matrix with the same value
void initializeSameValueMatrix(int **spins, int row, int column, double
                               value, double &M, double &E, double J)
{
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++ ){
            spins[i][j] = value;
        }
    }
    M = row*column*value;
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++){
            E -= (double)spins[i][j]*(spins[i][(j+column-1)%column
                    ] + spins[(i+row-1)%row][j]);
        }
    }
    E *= J;
}


//Initialize all elements of the matrix with random 1's or -1's
void initializeRamdomSpinMatrix(int **matrix, int row, int column)
{
    mt19937 mt(1393);
    uniform_int_distribution<int> myDist(0,1);

    for(int i = 0; i < row; i++){
        for(int j = 0; j < row; j++){
            matrix[i][j] = 2*myDist(mt)-1;
        }
    }
}

void metropolis(int **spins, int L, double temperature, double J,
                double &E, double &M, mt19937 mt, uniform_int_distribution<int> intDist,
                uniform_real_distribution<double> realDist)
{
    //double beta = 1./(1.381e-23*temperature);
    double exponentials[3];
    double deltaE;
    double energies[3] = {0, 4*J, 8*J};
    int i,j;
    for(int i = 0; i < 3; i++){
        exponentials[i] = exp(-energies[i]/temperature);
    }
    for(int k = 0; k < L*L; k++){
        i = intDist(mt);
        j = intDist(mt);
        deltaE = 2*J*spins[i][j]*(spins[i][(j+L-1)%L]+spins[i][(j+L+1)%L]
                +spins[(i+L-1)%L][j]+spins[(i+L+1)%L][j]);
        if (deltaE < 0){
            spins[i][j] *= -1;
            E += (double) deltaE;
            M += (double) 2*spins[i][j];
        }else{
            if(realDist(mt) <= exponentials[(int)((2*deltaE+8*J)/(8*J) - 1)]){
                spins[i][j] *= -1;
                E += (double) deltaE;
                M += (double) 2*spins[i][j];
            }
        }
    }
}

void output(double temperature, double E, ofstream &outFile)
{
    outFile << temperature << "\t" << E << endl;
}
