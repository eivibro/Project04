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
void initializeRamdomSpinMatrix(int **spins, int row, int column, double &M, int spin_indices[],
                                double &E, double J, mt19937 &mt, uniform_int_distribution<int> &myDist)
{
    for(int i = 0; i < row; i++){
        for(int j = 0; j < row; j++){
            spins[i][j] = 2*myDist(mt)-1;
        }
    }
    for(int i = 0; i < row; i++){
        for(int j = 0; j < column; j++){
            E -= (double)spins[i][j]*(spins[i][spin_indices[j+2]] + spins[spin_indices[i+2]][j]);
            M += spins[i][j];
        }
    }
    E *= J;
}

//Function running one Monte Carlo cycle
void metropolis(int **spins, int L, double J, int spin_indices[], int &accepted_configurations,
                double &E, double &M, mt19937 &mt, uniform_int_distribution<int> &intDist,
                uniform_real_distribution<double> &realDist, double exponentials[])
{
    accepted_configurations = 0;
    int i,j;
    double deltaE;
    for(int k = 0; k < L*L; k++){
        i = intDist(mt);
        j = intDist(mt);
        deltaE = 2*J*spins[i][j]*(spins[i][spin_indices[j]]+spins[i][spin_indices[j+2]]
                +spins[spin_indices[i]][j]+spins[spin_indices[i+2]][j]);
        if (deltaE < 0 || realDist(mt) <= exponentials[(int)deltaE]){
            accepted_configurations++;
            spins[i][j] *= -1;
            E += (double) deltaE;
            M += (double) 2*spins[i][j];
        }
    }
}

//Methods to write to file
void output(double temperature, int mcCycles, int L, double average[], ofstream &outFile)
{
    double inverse_number_spins = 1./(L*L);
    double ExpEnergy = average[0]/mcCycles;
    double ExpEnergy2 = average[1]/mcCycles;
    double ExpAbsMagnetization = average[2]/mcCycles;
    double ExpMagnetization2 = average[3]/mcCycles;
    double ExpMagnetization = average[4]/mcCycles;

    double heatCapacity = (ExpEnergy2-ExpEnergy*ExpEnergy)/(temperature*temperature);
    double suceptibility = (ExpMagnetization2-ExpAbsMagnetization*ExpAbsMagnetization)/temperature;
    outFile << setw(15) << setprecision(8) << temperature;
    outFile << setw(15) << setprecision(8) << ExpEnergy*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << ExpAbsMagnetization*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << heatCapacity*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << suceptibility*inverse_number_spins << endl;
 }

void outputc(double temperature, int mcCycles, int L,
             double average[], int accepted_configurations, ofstream &outFile)
{
    double inverse_number_spins = 1./(L*L);
    double ExpEnergy = average[0]/mcCycles;
    double ExpEnergy2 = average[1]/mcCycles;
    double ExpAbsMagnetization = average[2]/mcCycles;
    double ExpMagnetization2 = average[3]/mcCycles;
    double ExpMagnetization = average[4]/mcCycles;

    double heatCapacity = (ExpEnergy2-ExpEnergy*ExpEnergy)/(temperature*temperature);
    double suceptibility = (ExpMagnetization2-ExpAbsMagnetization*ExpAbsMagnetization)/temperature;
    double suceptibilityAbs = (ExpMagnetization2-ExpAbsMagnetization*ExpAbsMagnetization)/temperature;
    outFile << setw(15) << setprecision(8) << mcCycles;
    outFile << setw(15) << setprecision(8) << ExpEnergy*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << ExpAbsMagnetization*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << heatCapacity*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << suceptibility*inverse_number_spins;
    outFile << setw(15) << setprecision(8) << accepted_configurations << endl;
 }
