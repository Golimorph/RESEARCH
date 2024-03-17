//
//  Ewald.h
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-10-23.
//

#ifndef Ewald_h
#define Ewald_h

#include "mpi.h"
#include "Parameters.h"
#include "Crystal.h"
#include "Eigen/Dense"
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iomanip>


class Ewald{
public:
    Ewald(Parameters parameters, Crystal crystal);
    
    
   
   
    Eigen::MatrixXd getEwaldMatrix();
    

    
private:
    Crystal crystal;
    Parameters parameters;
    friend class TrivialPar;
    friend class Microcanonical;
    
    void findEwaldMatrix();
    void save();
    void load();
    
    
    
    int myrank;
    int nproc;
    
    double A(double r);
    double B(double r);
    double alpha;
    double directSum(int i,int j, int numberOfMirrors);
    double reciprocalSum(int i,int j, int numberOfMirrors);
    double selfSum(int i,int j);
    void addEwaldDemagBoundaryTermToEwaldMatrix();
    
    
    
   
   
    
    Eigen::MatrixXd ewaldMatrix;
    bool ewaldMatrixFound;
   
};





























#endif /* Ewald_h */
