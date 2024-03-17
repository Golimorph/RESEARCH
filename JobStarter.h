//
//  JobStarter.h
//
//
//  Created by Richard Edberg on 2018-03-02.
//


#ifndef JobStarter_h
#define JobStarter_h


#include "mpi.h"
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include <stdlib.h>
#include <ctime>
#include "Parameters.h"
#include "Crystal.h"
#include "Ewald.h"
#include "TrivialPar.h"
#include "Microcanonical.h"
#include "ParallelTempering.h"
#include "FileWriter.h"
#include <unistd.h>

class JobStarter{
public:
    explicit JobStarter (int argc, char *argv[], std::string pathToDataFolder);
    
    void start();
    void startParallelTempering();
    void startTrivialMPI();
    void startMicrocanonical();
    
private:
    
    int myrank;
    int nproc;
    std::string pathToDataFolder;
    std::string argFilename;    
    
    
    
    
    
    
};


#endif /* JobStarter_h */

















