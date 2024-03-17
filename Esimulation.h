//
//  Esimulation.h
//  exactPyrochloreDipole
//
//  Created by Richard Edberg on 2017-11-23.
//  Copyright © 2017 Richard Edberg. All rights reserved.


#ifndef Esimulation_h
#define Esimulation_h
#include "Parameters.h"
#include "Crystal.h"
#include "Ewald.h"

#include <iomanip>
#include <stdio.h>
#include <iostream>
#include "Eigen/Dense"
#include <fstream>
#include <limits.h>


class Esimulation{
public:
    Esimulation(Parameters parameters, Crystal crystal, Ewald ewald);
    
   
    Eigen::MatrixXd getDataMatrix();
    Eigen::MatrixXd run();
   //void saveToFile();
    
    
private:
    
    Eigen::RowVectorXi increase(Eigen::RowVectorXi spinOrientation, int i);
    double Energy();
    
    
    Eigen::MatrixXd data;
    
    
    Parameters parameters;
    Crystal crystal;
    Ewald ewald;
    
    double T;
    double B;
    double JnnzStart;
    double JnnzEnd;
    
    
    
    double Jnnxy;
    double Jnnz;
    double Jnnz1;
    double Jnnz2;
    double Jnnn;
    double muNorm;
    Eigen::Vector3d fieldDirection;
    
    
    Eigen::MatrixXd lattice;
    Eigen::MatrixXd spinDirection;
    Eigen::RowVectorXi spinOrientation;
    Eigen::MatrixXi NNlist;
    Eigen::MatrixXi NNNlist;
    Eigen::MatrixXi XYNNlist;
    Eigen::MatrixXi ZNNlist;
    Eigen::MatrixXi Z1NNlist;
    Eigen::MatrixXi Z2NNlist;
    Eigen::MatrixXd ewaldMatrix;
    

    
    
    
    
  
    
};










#endif /* Esimulation_h */
