#include <iostream>
#include "lib.h"
#include "functions.h"
#include <unittest++/UnitTest++.h>
#include <random>
#include <mpi.h>
#include <QElapsedTimer>


using namespace std;

int main(int nargs, char* args[])
{

    //Random generator engine and distributions
//    mt19937 mt(4724);
//    uniform_int_distribution<int> intDist(0,L-1);
//    uniform_real_distribution<double> realDist(0,1);

    //Initializing
//    int mcCycles = 1000000;
//    double E,M;
//    E = M = 0;
//    double startTemperature = 1.0;
//    double endTemperature = 3;
//    double tempStep = 0.1;
//    double J = 1;
//    int **spins;
//    spins = (int **)matrix(L, L, sizeof(double));
//    initializeSameValueMatrix(spins, L, L, 1, M, E, J);
//    cout << "Initial energy:            " << E << endl;
//    cout << "Initial magnetic momentum: " << M << endl;
//    ofstream outfile;
//    outfile.open("output.txt");
//    double *average;
//    average = new double [5];
//    double exponentials[17];
//    for(double t = startTemperature; t <= endTemperature; t+=tempStep){
//        for(int i = 0; i < 17; i+=8){
//            exponentials[i] = exp(-(double)i*J/t);
//        }
//        for(int i = 0; i < 5; i++){
//            average[i] = 0;
//        }
//        for(int i = 0; i < mcCycles; i++){
//            metropolis(spins, L, J, E, M, mt, intDist, realDist, exponentials);
//            average[0] += E;
//            average[1] += E*E;
//            average[2] += fabs(M);
//            average[3] += M*M;
//            average[4] += M;
//        }
//        cout << "Temperature = " << t << endl;
//        output(t, mcCycles, L, average, outfile);
//    }
//    outfile.close();
//    free_matrix((void **) spins);

    //Paralellizing

    //Initializing
    int L = 2;
    uniform_int_distribution<int> intDist(0,L-1);
    uniform_real_distribution<double> realDist(0,1);
    string fileJoining, fileRemoving;
    int mcCycles = 1000000;
    double E,M;
    E = M = 0;
    double startTemperature = 1.0;
    double endTemperature = 3;
    double tempStep = 0.1;
    double J = 1;
    int **spins, my_rank, num_of_threads;
    spins = (int **)matrix(L, L, sizeof(double));
    initializeSameValueMatrix(spins, L, L, 1, M, E, J);
    cout << "Initial energy:            " << E << endl;
    cout << "Initial magnetic momentum: " << M << endl;
    double *average;
    average = new double [5];
    double exponentials[17];

    MPI_Init (&nargs, &args);;
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double temperature_per_thread = (endTemperature-startTemperature)/(double)num_of_threads;
    cout << temperature_per_thread << endl;
    double start_temp = startTemperature+temperature_per_thread*my_rank;
    double end_temp = start_temp+temperature_per_thread;
    if(my_rank == num_of_threads-1){
        end_temp = endTemperature+tempStep;
    }
    ofstream outfile;
    outfile.open("output"+to_string(my_rank)+".txt");
    cout << "My rank: " << my_rank << ", and start_temp: " << start_temp << ", and end_temp: " << end_temp << endl;

    QElapsedTimer timer;
    timer.start();
    mt19937 mt(343+my_rank);
    for(double t = start_temp; t <= end_temp; t+=tempStep){
        for(int i = 0; i < 17; i+=8){
            exponentials[i] = exp(-(double)i*J/t);
        }
        for(int i = 0; i < 5; i++){
            average[i] = 0;
        }
        for(int i = 0; i < mcCycles; i++){
            metropolis(spins, L, J, E, M, mt, intDist, realDist, exponentials);
            average[0] += E;
            average[1] += E*E;
            average[2] += fabs(M);
            average[3] += M*M;
            average[4] += M;
        }
        cout << "Temperature = " << t << endl;
        output(t, mcCycles, L, average, outfile);
    }
    outfile.close();
    free_matrix((void **) spins);

    //File editing
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        system("echo >> joinedResults.txt");
        system("rm joinedResults.txt");
        for(int i = 0; i < num_of_threads; i++){
            fileJoining = "cat output"+to_string(i)+".txt >> joinedResults.txt";
            fileRemoving = "rm output"+to_string(i)+".txt";
            cout << fileJoining << endl;
            system(fileJoining.c_str());
            system(fileRemoving.c_str());
        }
    }
    MPI_Finalize ();
    cout << "Elapsed time: " << timer.elapsed()/1000. << endl;
    return 0;
}




