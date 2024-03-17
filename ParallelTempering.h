//
//  MCsimulation.h
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-08-30.
//
//
//#pragma once

#ifndef ParallelTempering_h
#define ParallelTempering_h



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





class ParallelTempering {
public:
    
    ParallelTempering(Parameters parameters, Crystal crystal, Ewald ewald , Eigen::VectorXd *spinOrientation, int myrank, int nproc);
    Eigen::MatrixXd run();
    void saveStateToFile(std::string path, int i);
    
    
private:
    
    
    
    
    
    Parameters parameters;
    Crystal crystal;
    Ewald ewald;
    rand_class rng_s;
    
    
    //----------------------------
    void MCsimRun();
    void findElist();
    bool isEqual(double a, double b, double tol);
    bool findLoop();
    double Energy();
    double singleSweepList();
    double loopSweepList();
    void findTemperatureList();
    void replicaExchange(int oddEven);
     //----------------------------
    Eigen::RowVector3d Magnetization();
    std::vector<int> shortSpinLoop;
    //----------------------------
    int myrank;
    int nproc;
    int datapoints;
    int Nspins;
    int bins;
    int stepsAcceptedInBin;
    int loopStepsAcceptedInBin;
    double Jnnxy;
    double Jnnz1;
    double Jnnz2;
    double Jnnn;
    double D;
    double T;
    double B;
    double muNorm;
    double E;
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
    Eigen::MatrixXd ewaldMatrix;
    Eigen::MatrixXd data;
    Eigen::VectorXd *spinOrientation;
    Eigen::VectorXd Elist;
    Eigen::Vector3d fieldDirection;
    Eigen::VectorXi temperatureList;
    Eigen::VectorXi rankList;
    Eigen::RowVector3d M;
 
    
    
    
    
};


#endif /* ParallelTempering_h */










