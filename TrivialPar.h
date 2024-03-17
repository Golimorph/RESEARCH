//
//  MCsimulation.h
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-08-30.
//
//
//#pragma once

#ifndef TrivialPar_h
#define TrivialPar_h



#include <chrono>
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include "Crystal.h"
#include "Parameters.h"
#include "Ewald.h"
#include <iomanip>
#include <ctime>
#include "rand_gen.h"
#include "mpi.h"





class TrivialPar {
public:
    TrivialPar(Parameters *parameters, Crystal *crystal, Ewald *ewald , Eigen::VectorXd *spinOrientation);
    Eigen::MatrixXd run();
    Eigen::MatrixXd run(int myrank, int size);
    void test();
    void saveToFile(std::string path);
    void saveStateToFile(std::string path, int i);
    
    
    Eigen::MatrixXd getSumCorrSiSj();
    Eigen::MatrixXd getData();
    
    
    

private:
    friend class JobStarter;
    
    Parameters *parameters;
    rand_class rng_s;
    Ewald *ewald;
    //----------------------------
    void MCsimRun(unsigned long long int stepsPerBin);
    void findElist();
    inline void AddToCorrelationSiSj();
    bool findLoop();
    double Energy();
    double singleSweepList();
    double loopSweepList();
    //----------------------------
    Eigen::RowVector3d Magnetization();
    std::vector<int> shortSpinLoop;
    //----------------------------
    int datapoints;
    int Nspins;
    int bins;
    int stepsAcceptedInBin;
    int loopStepsAcceptedInBin;
    double Jnnxy;
    //double Jnnz;
    double Jnnz1;
    double Jnnz2;
    double JnnzStart;
    double JnnzEnd;
    double Jnnn;
    double D;
    double T;
    double B;
    double muNorm;
    double E;
    signed long long int MCstepsPerDatapoint;
    signed long long int stepsPerBin;
    //----------------------------
    Eigen::MatrixXi NNlist;
    Eigen::MatrixXi NNNlist;
    Eigen::MatrixXi XYNNlist;
    Eigen::MatrixXi ZNNlist;
    Eigen::MatrixXi Z1NNlist;
    Eigen::MatrixXi Z2NNlist;
    Eigen::MatrixXi spinTetraList;
    Eigen::MatrixXi tetraSpinList;
    Eigen::MatrixXd lattice;
    Eigen::MatrixXd spinDirection;
    Eigen::MatrixXd diamondLattice;
    //Eigen::MatrixXd ewaldMatrix;
    Eigen::MatrixXd data;
    Eigen::VectorXd *spinOrientation;
    Eigen::VectorXd Elist;
    Eigen::Vector3d fieldDirection;
    Eigen::RowVector3d M;
    Eigen::MatrixXd sumCorrSiSj;
    
 
};


#endif /* TrivialPar_h */










