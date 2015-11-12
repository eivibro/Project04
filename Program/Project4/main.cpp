#include <iostream>
#include "lib.h"
#include "functions.h"
#include <random>
#include <mpi.h>
#include <QElapsedTimer>

using namespace std;

int main(int nargs, char* args[])
{
    //Initializing
    int L = 100;
    uniform_int_distribution<int> intDist(0,L-1);
    uniform_int_distribution<int> intDistSingle(0,1);
    uniform_real_distribution<double> realDist(0,1);
    int **spins, spin_indices[L+2],mcCycles = 20000;
    int skip = 10000;
    double E,M;
    E = M = 0;
    double *average;
    average = new double [5];
    double exponentials[9];
    double J = 1;
    int accepted_configurations = 0;
    spins = (int **)matrix(L, L, sizeof(double));
    //Array to deal with periodic boundary conditions
    for(int i = 1; i<L+1; i++){
        spin_indices[i] = i-1;
    }
    spin_indices[0] = L-1;
    spin_indices[L+1] = 0;
    for(int i = 0; i < L+2; i++){
        cout << spin_indices[i] << ", ";
    }
    cout << endl;
    initializeSameValueMatrix(spins, L, L, 1, M, E, J);
    //mt19937 mt(372);
    //initializeRamdomSpinMatrix(spins, L, L, M, spin_indices, E, J, mt, intDistSingle);


//    for(int i = 0; i < L; i++){
//        for(int j = 0; j < L; j++){
//            cout << spins[i][j] <<  " ";
//        }
//        cout << endl;
//    }
//        cout << "Initial magnetization: " << M << endl;
//        cout << "Initial Energy: " << E << endl;

//    ofstream outc;
//    outc.open("taskCOutR.txt");
//    double temperature = 2.4;
//    for(int i = 0; i < 17; i+=8){
//        exponentials[i] = exp(-(double)i*J/temperature);
//    }
//    for(int i = 0; i < 5; i++){
//        average[i] = 0;
//    }
//    for(int i = 1; i <= mcCycles; i++){
//        metropolis(spins, L, J, spin_indices, accepted_configurations,
//                   E, M, mt, intDist, realDist, exponentials);
//        average[0] += E;
//        average[1] += E*E;
//        average[2] += fabs(M);
//        average[3] += M*M;
//        average[4] += M;
//        outputc(temperature, i, L, average, accepted_configurations, outc);
//    }
//    outc.close();


    //Parallel version
    string fileJoining, fileRemoving;
    int my_rank, num_of_cpus;


    QElapsedTimer timer;
    timer.start();

    MPI_Init (&nargs, &args);;
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_cpus);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double startTemperature = 2.0;
    double endTemperature = 2.6;
    int t_step = 16;
    double deltaT = (endTemperature - startTemperature) / (t_step-1);
    vector<double> temperatures(t_step, 0);
    for(int i=0; i<t_step; i++) {
        temperatures[i] = startTemperature + i*deltaT;
        cout << temperatures[i] << endl;
    }
    int numTemperaturesPerProcessor = t_step / num_of_cpus;
    int TindexStart = my_rank*numTemperaturesPerProcessor;
    int TindexStop = (my_rank+1)*numTemperaturesPerProcessor;
    if(my_rank == num_of_cpus-1) TindexStop = t_step-1;
    // 0,1,2
    // 3,4,5
    // 6,7,8
    // 9,10,11,12

//    double temperatures_per_thread = (endTemperature-startTemperature)/t_step;
//    double start_temp = startTemperature+int(temperatures_per_thread)*t_step*my_rank/num_of_threads;
//    double end_temp = start_temp+int(temperatures_per_thread)*t_step/num_of_threads;
//    if(my_rank == num_of_threads-1){
//        end_temp = endTemperature;
//    }


    ofstream outfile;
    outfile.open("output"+to_string(my_rank)+".txt");
    cout << "My rank: " << my_rank << ", and start_temp: " << temperatures[TindexStart]
            << ", and end_temp: " << temperatures[TindexStop] << endl;

    mt19937 mt(343+my_rank);
    for(int i = 0; i < 9; i++){
        exponentials[i] = 0;
    }
    for(double t = temperatures[TindexStart]; t < temperatures[TindexStop]; t+=deltaT){
        for(int i = 0; i < 9; i+=4){
            exponentials[i] = exp(-(double)i*J/t);
        }
        for(int i = 0; i < 5; i++){
            average[i] = 0;
        }
        //Running Monte Carlo cycles before starting to record data
        for(int i = 0; i < skip; i++){
            metropolis(spins, L, J, spin_indices, accepted_configurations, E, M, mt, intDist, realDist, exponentials);
        }
        //Recording data
        for(int i = 0; i < mcCycles-skip; i++){
            metropolis(spins, L, J, spin_indices, accepted_configurations, E, M, mt, intDist, realDist, exponentials);
            average[0] += E;
            average[1] += E*E;
            average[2] += fabs(M);
            average[3] += M*M;
            average[4] += M;
        }
        cout << "Temperature = " << t << endl;
        output(t, mcCycles-skip, L, average, outfile);
    }
    outfile.close();
    free_matrix((void **) spins);

    //File editing
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        system("echo >> joinedResults.txt");
        system("rm joinedResults.txt");
        for(int i = 0; i < num_of_cpus; i++){
            fileJoining = "cat output"+to_string(i)+".txt >> joinedResults.txt";
            fileRemoving = "rm output"+to_string(i)+".txt";
            cout << fileJoining << endl;
            system(fileJoining.c_str());
            system(fileRemoving.c_str());
        }
    }
    MPI_Finalize();
    cout << "Elapsed time: " << timer.elapsed()/1000. << endl;
    return 0;
}




