//
//  Parameters.h
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-10-05.
//



#ifndef Parameters_h
#define Parameters_h
#include <stdio.h>
#include <math.h>

#include <complex>  
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include "Eigen/Dense"


class Parameters
{
public:
    Parameters(std::string path, std::string parametersFilename);
    
    void load();
    std::string print() const; 
    
    std::string getPathToDataFolder() const;
    std::string getOutputFilename() const;
    std::string getPIndex() const;
    
    bool getAddEwaldDemagTerm() const;
    bool getLoadCrystal() const;
    bool getLoadEwald() const;
    bool getAddDateStamp() const;
    bool getWriteCorrelationMatrix() const;
    bool getReadConfiguration() const;
    
    bool getLogSteps() const;
    int getCellType() const;
    int getRunTime() const;
    int getMaxEwaldMirrors() const;
    int getCellX() const;
    int getCellY() const;
    int getCellZ() const;
    int getDataPoints() const;
    int getBins() const;
    int getTimingSteps() const;
    int getReplicaExchanges() const;
    signed long long int getMCsteps() const;
   
    
    double getEqPercentage() const;
    double getTstart() const;
    double getTend() const;
    double getJnnxy() const;
    double getJnnz1() const;
    double getJnnz2() const;
    double getJnnn() const;
    double getD() const;
    double getkb() const;
    double getBstart() const;
    double getBend() const;
    double getMu() const;
    double getMuNorm() const;
    double getUnitCellSideLength() const;
    
    
    
    Eigen::Vector3d getPDirection() const;
    Eigen::Vector3d getFieldDirection() const;
    Eigen::Matrix3d getFmatrix() const;
    Eigen::Matrix3d getDemagMatrix() const;
    
private:
    void findRuntime();
    void CheckInputandGetParameters();
    void findFieldDirection();
    void findPDirection();
    void createFmatrixFromOneMkappa();
    void findCellDimensions();
    void checkInput();
    
    std::string path;
    std::string parametersFilename;
    std::string outputFilename;
    std::string endTimeStr;
    std::string fieldIndex;
    std::string pIndex;
   
    
    
    
    
    
    bool addEwaldDemagTerm;
    bool loadCrystal;
    bool loadEwald;
    bool addDateStamp;
    bool logSteps;
    bool writeCorrelationMatrix;
    bool readConfiguration;
   
    
    int runTime;
    int cellSide;
    int cellType;
    int cellX;
    int cellY;
    int cellZ;
    int MaxEwaldMirrors;
    int dataPoints;
    int replicaExchanges;
    int bins;
    int timingSteps;
    signed long long int MCsteps;
  
  
    double oneMkappa;
    double fieldPhi;
    double fieldTheta;
    double unitCellSideLength;
    double Tstart;
    double Tend;
    double Jnnxy;
    double Jnnz1;
    double Jnnz2;
    double JnnzStart;
    double JnnzEnd;
    double Jnnn;
    double D;
    double kb;
    double Bstart;
    double Bend;
    double mu;
    
    double mubDivkb;
    double eqPercentage;
    
    Eigen::Matrix3d demagMatrix;
    Eigen::Matrix3d F;
    Eigen::Vector3d fieldDirection;
    Eigen::Vector3d pDirection;
  
};

#endif /* Parameters_h */


