//
//  Crystal.h
//
//
//  Created by Richard Edberg on 2017-07-20.
//


#ifndef Crystal_h
#define Crystal_h

#include "Parameters.h"
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>


class Crystal {
public:
    explicit Crystal (Parameters parameters);
    

    void findLattice();
    void findUpDiamondLattice();
    void findDownDiamondLattice();
    void findSpinDirection();
    void findLists();
    void findNN();
    void findNNN();
    void findTetrahedronLists();
    void saveLists();
    void loadLists();
    
    int getNspins() const;
    int getUspins() const;
    
    double getNNdist() const;
    double getCellSideLengthX() const;
    double getCellSideLengthY() const;
    double getCellSideLengthZ() const;
    double getCellVolume() const;
   
    Eigen::MatrixXi getNNlist() const;
    Eigen::MatrixXi getNNNlist() const;
    Eigen::MatrixXi getXYNNlist() const;
    Eigen::MatrixXi getZNNlist() const;
    Eigen::MatrixXi getZ1NNlist() const;
    Eigen::MatrixXi getZ2NNlist() const;
    Eigen::MatrixXi getSpinTetraList() const;
    Eigen::MatrixXi getTetraSpinList() const;
    Eigen::MatrixXd getLattice() const;
    Eigen::MatrixXd getUpDiamondLattice() const;
    Eigen::MatrixXd getDownDiamondLattice() const;
    Eigen::MatrixXd getSpinDirection() const;
    Eigen::MatrixXd getCrystalSideTranslation() const;
    
    Parameters parameters;

private:
    
   
    
    
    
    Eigen::MatrixXd unitCell() const;
    bool isEqual(double a, double b, double tol);
    
    
    

    bool centering;
    bool deform;
    
    int Uspins;
    int Nspins;
  
    
    double NNdist;
    double NNNdist;
    double unitCellXSideLength;
    double unitCellYSideLength;
    double unitCellZSideLength;
    
    
    
    double tetraCornerCenterDist;
    
    Eigen::MatrixXi NNlist;
    Eigen::MatrixXi NNNlist;
    Eigen::MatrixXi XYNNlist;
    Eigen::MatrixXi ZNNlist;
    Eigen::MatrixXi Z1NNlist;
    Eigen::MatrixXi Z2NNlist;
    Eigen::MatrixXi spinTetraList;
    Eigen::MatrixXi tetraSpinList;
    
    Eigen::Matrix3d crystalSideTranslation;
    Eigen::Matrix3d deformedCrystalSideTranslation;
    
    Eigen::MatrixXd lattice;
    Eigen::MatrixXd deformedLattice;
    Eigen::MatrixXd upDiamondLattice;
    Eigen::MatrixXd deformedUpDiamondLattice;
    Eigen::MatrixXd downDiamondLattice;
    Eigen::MatrixXd deformedDownDiamondLattice;
    Eigen::MatrixXd spinDirection;
    Eigen::MatrixXd deformedSpinDirection;
 
};


#endif /* Crystal_h */
