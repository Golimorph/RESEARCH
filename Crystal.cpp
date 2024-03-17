//
//  Crystal.cpp
//  
//
//  Created by Richard Edberg on 2017-07-20.
//
//This class is intended to be used for creating an nxnxn lattice of the pyrochlore type. It can also load and save to files in order to save time.


#include "Crystal.h"


Crystal::Crystal(Parameters parameters)
:
parameters(parameters)
{
    
    
    
    Uspins=parameters.getCellType();
    Nspins=Uspins*parameters.getCellX()*parameters.getCellY()*parameters.getCellZ();
    
    if(Uspins==16)
    {
        unitCellXSideLength=1;
        unitCellYSideLength=1;
        unitCellZSideLength=1;
    }else if(Uspins==24)
    {
        unitCellXSideLength=1./sqrt(2);
        unitCellYSideLength=sqrt(3./2);
        unitCellZSideLength=sqrt(3);
    }else{
        std::cerr<<"MCsimulation: Something wrong with number of Uspins";
    }
    
    
    NNdist=1/sqrt(8);
    NNNdist=sqrt(3)*NNdist;
    tetraCornerCenterDist=sqrt(3./8.)*NNdist;
    centering=true;//puts the lattice with the origin in its center (rather than the corner)
    deform=true;//deforms the crystal and spins using the F matrix, keep true
    
   

    
    //pyrochlore-16cell-specific matrix size initialization
    lattice = Eigen::MatrixXd(Nspins,3);
    upDiamondLattice = Eigen::MatrixXd(Nspins/4,3);
    downDiamondLattice = Eigen::MatrixXd(Nspins/4,3);
    spinDirection = Eigen::MatrixXd(Nspins,3);
    NNlist= -Eigen::MatrixXi::Ones(Nspins,6);
    NNNlist= -Eigen::MatrixXi::Ones(Nspins,12);
    XYNNlist =  -Eigen::MatrixXi::Ones(Nspins,7);
    ZNNlist = -Eigen::MatrixXi::Ones(Nspins,7);
    Z1NNlist = -Eigen::MatrixXi::Ones(Nspins,7);
    Z2NNlist = -Eigen::MatrixXi::Ones(Nspins,7);
    
    
    spinTetraList= -Eigen::MatrixXi::Ones(Nspins,2);
    tetraSpinList= -Eigen::MatrixXi::Ones(Nspins/2,4);
    
    
    
    
    
    
    //the deformed lattice and spins
    deformedLattice=lattice;
    deformedSpinDirection=spinDirection;
    deformedUpDiamondLattice=upDiamondLattice;
    deformedDownDiamondLattice=downDiamondLattice;
    //translation vectors
    crystalSideTranslation.row(0)<<parameters.getCellX()*unitCellXSideLength,0.,0.;
    crystalSideTranslation.row(1)<<0.,parameters.getCellY()*unitCellYSideLength,0.;
    crystalSideTranslation.row(2)<<0.,0.,parameters.getCellZ()*unitCellZSideLength;
    
    for(int i=0; i<3; ++i)
    {
        deformedCrystalSideTranslation.row(i)= (parameters.getFmatrix()*(crystalSideTranslation.row(i)).transpose()).transpose();
    }
 
    
    findLattice();
    findUpDiamondLattice();
    findDownDiamondLattice();
    findSpinDirection();
    std::cout<<"Crystal initialized"<<"\n";
    

     
}

Eigen::MatrixXd Crystal::getCrystalSideTranslation() const
{
    return deformedCrystalSideTranslation;
}
double Crystal::getCellSideLengthX() const
{
    return (deformedCrystalSideTranslation.row(0)).norm();
}

double Crystal::getCellSideLengthY() const
{
    return (deformedCrystalSideTranslation.row(1)).norm();
}

double Crystal::getCellSideLengthZ() const
{
    return (deformedCrystalSideTranslation.row(2)).norm();
}

double Crystal::getCellVolume() const
{
    return parameters.getCellX()*parameters.getCellY()*parameters.getCellZ()*unitCellXSideLength*unitCellYSideLength*unitCellZSideLength*(parameters.getFmatrix()).determinant();
}
Eigen::MatrixXd Crystal::getLattice() const
{
    if(deform)
    {
        return deformedLattice;
    }
    return lattice;
}
Eigen::MatrixXd Crystal::getUpDiamondLattice() const
{
    if(deform)
    {
        return deformedUpDiamondLattice;
    }
    return upDiamondLattice;
}





Eigen::MatrixXd Crystal::getDownDiamondLattice() const
{
    if(deform)
    {
        return deformedDownDiamondLattice;
    }

    return downDiamondLattice;
}



Eigen::MatrixXi Crystal::getSpinTetraList() const
{
    return spinTetraList;
}







Eigen::MatrixXi Crystal::getTetraSpinList() const
{
    return tetraSpinList;
}






Eigen::MatrixXd Crystal::getSpinDirection() const
{
    if(deform)
    {
        return deformedSpinDirection;
    }
    return spinDirection;
}








Eigen::MatrixXi Crystal::getNNlist() const
{
    return NNlist;
}

Eigen::MatrixXi Crystal::getNNNlist() const
{
    return NNNlist;
}
Eigen::MatrixXi Crystal::getXYNNlist() const
{
    return XYNNlist;
}
Eigen::MatrixXi Crystal::getZNNlist() const
{
    return ZNNlist;
}
Eigen::MatrixXi Crystal::getZ1NNlist() const
{
    return Z1NNlist;
}
Eigen::MatrixXi Crystal::getZ2NNlist() const
{
    return Z2NNlist;
}


int Crystal::getNspins() const
{
    return Nspins;
}
int Crystal::getUspins() const
{
    return Uspins;
}
double Crystal::getNNdist() const{
    return NNdist;
}









void Crystal::findLists()
{//This void finds lists that take a long time to compute
    findNN();
    //findNNN();
    findTetrahedronLists();
   
    std::cout<<"Crystal lists found\n";

}













Eigen::MatrixXd  Crystal::unitCell() const{
    //This function returns a 3xunitSpins matrix with the points of the Crystal unit cell with NN distance = 1.
    Eigen::MatrixXd unitCell(Uspins,3);
    //normalized positions with respect to one fourth of the cell cellSide
    //First tetrahedron
    
    if(Uspins==16)
    {
        unitCell.row(0) << 0,0,0;
        unitCell.row(1) << 1,1,0;
        unitCell.row(2) << 0,1,1;
        unitCell.row(3) << 1,0,1;
        //Second tetrahedron
        unitCell.row(4) << 2,2,0;
        unitCell.row(5) << 3,3,0;
        unitCell.row(6) << 2,3,1;
        unitCell.row(7) << 3,2,1;
        //Third tetrahedron
        unitCell.row(8) << 0,2,2;
        unitCell.row(9) << 1,3,2;
        unitCell.row(10)<< 0,3,3;
        unitCell.row(11)<< 1,2,3;
        //Fourth tetrahedron
        unitCell.row(12)<< 2,0,2;
        unitCell.row(13)<< 3,1,2;
        unitCell.row(14)<< 2,1,3;
        unitCell.row(15)<< 3,0,3;
        unitCell=unitCell/4.;
    }
        
    if(Uspins==24){
        unitCell(0,0) =0.  ; unitCell(0,1) =0.    ; unitCell(0,2) =0.  ;
        unitCell(1,0) =1./4; unitCell(1,1) =1. /4 ; unitCell(1,2) =0.  ;
        unitCell(2,0) =3./4; unitCell(2,1) =1. /4 ; unitCell(2,2) =0.  ;
        unitCell(3,0) =1./2; unitCell(3,1) =1. /2 ; unitCell(3,2) =0.  ;
     
        unitCell(4,0) =1./4; unitCell(4,1) =3. /4 ; unitCell(4,2) =0.  ;
        unitCell(5,0) =3./4; unitCell(5,1) =3. /4 ; unitCell(5,2) =0.  ;
        unitCell(6,0) =0.  ; unitCell(6,1) =1. /6 ; unitCell(6,2) =1./6;
        unitCell(7,0) =1./2; unitCell(7,1) =2. /3 ; unitCell(7,2) =1./6;
     
        unitCell(8,0) =1./4; unitCell(8,1) =1. /12; unitCell(8,2) =1./3;
        unitCell(9,0) =0.  ; unitCell(9,1) =1. /3 ; unitCell(9,2) =1./3;
        unitCell(10,0)=3./4; unitCell(10,1)=1. /12; unitCell(10,2)=1./3;
        unitCell(11,0)=1./4; unitCell(11,1)=7. /12; unitCell(11,2)=1./3;
     
        unitCell(12,0)=3./4; unitCell(12,1)=7. /12; unitCell(12,2)=1./3;
        unitCell(13,0)=1./2; unitCell(13,1)=10./12; unitCell(13,2)=1./3;
        unitCell(14,0)=0.  ; unitCell(14,1)=1. /2 ; unitCell(14,2)=1./2;
        unitCell(15,0)=1./2; unitCell(15,1)=0.    ; unitCell(15,2)=1./2;
     
        unitCell(16,0)=1./2; unitCell(16,1)=1. /6 ; unitCell(16,2)=2./3;
        unitCell(17,0)=1./4; unitCell(17,1)=5. /12; unitCell(17,2)=2./3;
        unitCell(18,0)=3./4; unitCell(18,1)=5. /12; unitCell(18,2)=2./3;
        unitCell(19,0)=0.  ; unitCell(19,1)=2. /3 ; unitCell(19,2)=2./3;
     
        unitCell(20,0)=1./4; unitCell(20,1)=11./12; unitCell(20,2)=2./3;
        unitCell(21,0)=3./4; unitCell(21,1)=11./12; unitCell(21,2)=2./3;
        unitCell(22,0)=1./2; unitCell(22,1)=1. /3 ; unitCell(22,2)=5./6;
        unitCell(23,0)=0.  ; unitCell(23,1)=5. /6 ; unitCell(23,2)=5./6;
        
        for(int i=0; i<24; ++i)
        {
            unitCell(i,0)=unitCell(i,0)*1./sqrt(2);
            unitCell(i,1)=unitCell(i,1)*sqrt(3./2);
            unitCell(i,2)=unitCell(i,2)*sqrt(3);
        }
    }
  
    return unitCell;
}



void Crystal::findLattice()
{   //This void calls unitCell and generates a cellSide*cellSide*cellSide cubic crystal with NNdist=1/sqrt(8)
    
    Eigen::MatrixXd unitCell = Crystal::unitCell();
    
    for(int z=0;z<parameters.getCellZ();z++){
        for(int y=0;y<parameters.getCellY();y++){
            for(int x=0;x<parameters.getCellX();x++){
                Eigen::MatrixXd cell=unitCell;
                Eigen::RowVector3d T(unitCellXSideLength*x, unitCellYSideLength*y, unitCellZSideLength*z);
                cell.rowwise() += T;
                
                
                
                for(int i=0; i<cell.rows(); ++i)
                {
                    lattice.row(x*unitCell.rows()+y*parameters.getCellX()*unitCell.rows()+z*parameters.getCellX()*parameters.getCellY()*unitCell.rows()+i)=cell.row(i);
                }
                
                
            }
        }
    }
    //Fixes so that the lattice is centered in the origin
    if(centering) lattice.rowwise() -= lattice.colwise().mean();
    
    for(int i=0; i<Nspins; ++i)
    {
        deformedLattice.row(i)=(parameters.getFmatrix()*((lattice.row(i)).transpose())).transpose();
    }
}












void Crystal::findUpDiamondLattice()
{//this void uses the diamond cell to create an upp diamond lattice
   
    
    
    findLattice();
    
    if(Uspins==16){
        for(int i=0; i<lattice.rows()/4; ++i)
        {
            //Each tetrahedron
            upDiamondLattice.row(i) = (lattice.block<4,3>(4*i,0)).colwise().mean();
        }
    }
    if(Uspins==24){
        for(int uCellNo=0; uCellNo<Nspins/Uspins; ++uCellNo)
        {
            //Each tetrahedron
            upDiamondLattice.row(Uspins*uCellNo/4)   = (lattice.row(0+Uspins*uCellNo)+lattice.row(1+Uspins*uCellNo)+(lattice.row(2+Uspins*uCellNo)-crystalSideTranslation.row(0)/parameters.getCellX())+lattice.row(6+Uspins*uCellNo))/4;
            upDiamondLattice.row(Uspins*uCellNo/4+1) = (lattice.row(3+Uspins*uCellNo)+lattice.row(4+Uspins*uCellNo)+lattice.row(5+Uspins*uCellNo)+lattice.row(7+Uspins*uCellNo))/4;
            upDiamondLattice.row(Uspins*uCellNo/4+2)   = (lattice.row(9+Uspins*uCellNo)+lattice.row(11+Uspins*uCellNo)+(lattice.row(12+Uspins*uCellNo)-crystalSideTranslation.row(0)/parameters.getCellX())+lattice.row(14+Uspins*uCellNo))/4;
            upDiamondLattice.row(Uspins*uCellNo/4+3)   = (lattice.row(8+Uspins*uCellNo)+lattice.row(10+Uspins*uCellNo)+(lattice.row(13+Uspins*uCellNo)-crystalSideTranslation.row(1)/parameters.getCellY())+lattice.row(15+Uspins*uCellNo))/4;
            upDiamondLattice.row(Uspins*uCellNo/4+4) = (lattice.row(16+Uspins*uCellNo)+lattice.row(17+Uspins*uCellNo)+lattice.row(18+Uspins*uCellNo)+lattice.row(22+Uspins*uCellNo))/4;
            upDiamondLattice.row(Uspins*uCellNo/4+5)   = (lattice.row(19+Uspins*uCellNo)+lattice.row(20+Uspins*uCellNo)+(lattice.row(21+Uspins*uCellNo)-crystalSideTranslation.row(0)/parameters.getCellX())+lattice.row(23+Uspins*uCellNo))/4;
            
            
        }
    }
    
    
    
    
    for(int i=0; i<upDiamondLattice.rows(); ++i)
    {
        deformedUpDiamondLattice.row(i)=(parameters.getFmatrix()*((upDiamondLattice.row(i)).transpose())).transpose();
    }
    
}


void Crystal::findDownDiamondLattice()
{//this void uses the upp diamond lattice to create a down diamond lattice
    findUpDiamondLattice();
    
    downDiamondLattice= upDiamondLattice;
    
    if(Uspins==16)
    {
        for(int i=0; i<downDiamondLattice.rows() ; ++i )
        {
            downDiamondLattice.row(i)=downDiamondLattice.row(i)+2*(lattice.row(3)-(lattice.row(0)+lattice.row(1)+lattice.row(2)+lattice.row(3))/4);
            deformedDownDiamondLattice.row(i)= (parameters.getFmatrix()*((downDiamondLattice.row(i)).transpose()));
        }
    }
    if(Uspins==24)
    {
        for(int i=0; i<downDiamondLattice.rows() ; ++i )
        {
            downDiamondLattice(i,0)=upDiamondLattice(i,0);
            downDiamondLattice(i,1)=upDiamondLattice(i,1);
            downDiamondLattice(i,2)=upDiamondLattice(i,2)+2*tetraCornerCenterDist;
            deformedDownDiamondLattice.row(i)= (parameters.getFmatrix()*((downDiamondLattice.row(i)).transpose()));
        }
    }
}


void Crystal::findTetrahedronLists()
{//This void calls functions generates tetrahedron to spin and spin to tetrahedron lists
    
    
    if(Uspins==16)
    {
        Eigen::MatrixXi upList(lattice.rows()/4,4);
        Eigen::MatrixXi downList(lattice.rows()/4,4);
        spinTetraList = Eigen::MatrixXi(lattice.rows(),2);
        tetraSpinList = Eigen::MatrixXi(lattice.rows()/2,4);
        for(int i=0; i<lattice.rows(); ++i)
        {
            upList(i/4,i%4)=i;
            spinTetraList(i,0)=i/4;
        }
        double distanceFactor=1.1;
        double spinNumber=0;
        for(int i=0; i<downDiamondLattice.rows(); ++i)
        {
            spinNumber=0;
            for(int j=0; j<lattice.rows(); ++j)
            {
                for (int a=-1; a<2; a++)
                {
                    for (int b=-1; b<2; b++)
                    {
                        for (int c=-1; c<2; c++)
                        {
                            
                            if((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-downDiamondLattice.row(i)).norm()<tetraCornerCenterDist*distanceFactor)
                            {
                                downList(i,spinNumber++)=j;
                                spinTetraList(j,1)= (int) upList.rows()+i;
                            }
                        }
                    }
                }
            }
        }
        tetraSpinList << upList,
                         downList;
    }
    if(Uspins==24)
    {
        
        Eigen::MatrixXi upList(lattice.rows()/4,4);
        Eigen::MatrixXi downList(lattice.rows()/4,4);
        double distanceFactor=1.1;
        double spinNumber=0;
        for(int i=0; i<downDiamondLattice.rows(); ++i)
        {
            spinNumber=0;
            for(int j=0; j<lattice.rows(); ++j)
            {
                for (int a=-1; a<2; a++)
                {
                    for (int b=-1; b<2; b++)
                    {
                        for (int c=-1; c<2; c++)
                        {
                            if((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-upDiamondLattice.row(i)).norm()<tetraCornerCenterDist*distanceFactor)
                            {
                                upList(i,spinNumber++)=j;
                                spinTetraList(j,0)= (int) i;
                            }
                        }
                    }
                }
            }
            
            
            spinNumber=0;
            for(int j=0; j<lattice.rows(); ++j)
            {
                for (int a=-1; a<2; a++)
                {
                    for (int b=-1; b<2; b++)
                    {
                        for (int c=-1; c<2; c++)
                        {
                            if((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-downDiamondLattice.row(i)).norm()<tetraCornerCenterDist*distanceFactor)
                            {
                                downList(i,spinNumber++)=j;
                                spinTetraList(j,1)= (int) upList.rows()+i;
                            }
                        }
                    }
                }
            }
            
            
            
        }
        tetraSpinList << upList,
                         downList;
    }
  
}







void Crystal::findSpinDirection()
{
    //This void returns a matrix with the direction of each spin generated from the difference of the dual diamond lattice and the prochlore lattice
    spinDirection = Eigen::MatrixXd(Nspins,3);
    int j=-1;
    if(Uspins==16)
    {
        for (int i=0;i<Nspins;++i){
            if(i%4==0) j++;
            spinDirection.row(i) = upDiamondLattice.row(j) - lattice.row(i);
            spinDirection.row(i) =  (spinDirection.row(i))/((spinDirection.row(i)).norm());
        }
    }else if(Uspins==24)
    {
       
        for(int uCellNo=0; uCellNo<Nspins/Uspins; ++uCellNo)
        {
            spinDirection(   Uspins*uCellNo,0)=  0.          ;spinDirection(    Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection(    Uspins*uCellNo,2)=  1./3;
            spinDirection( 1+ Uspins*uCellNo,0)=  -sqrt(2./3.) ;spinDirection( 1+ Uspins*uCellNo,1)=   - sqrt(2.)/3. ;spinDirection( 1+ Uspins*uCellNo,2)= 1./3;
            spinDirection( 2+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection( 2+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection( 2+ Uspins*uCellNo,2)=  1./3;
            spinDirection( 3+ Uspins*uCellNo,0)=  0.          ;spinDirection( 3+ Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection( 3+ Uspins*uCellNo,2)=  1./3;
            
            spinDirection( 4+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection( 4+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection( 4+ Uspins*uCellNo,2)=  1./3;
            spinDirection( 5+ Uspins*uCellNo,0)= - sqrt(2./3.) ;spinDirection( 5+ Uspins*uCellNo,1)=  -  sqrt(2.)/3. ;spinDirection( 5+ Uspins*uCellNo,2)= 1./3;
            spinDirection( 6+ Uspins*uCellNo,0)=  0.          ;spinDirection( 6+ Uspins*uCellNo,1)=    0.          ;spinDirection( 6+ Uspins*uCellNo,2)= - 1.;
            spinDirection( 7+ Uspins*uCellNo,0)=  0.          ;spinDirection( 7+ Uspins*uCellNo,1)=    0.          ;spinDirection( 7+ Uspins*uCellNo,2)=  -1.;
            
            spinDirection( 8+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection( 8+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection( 8+ Uspins*uCellNo,2)=  1./3;
            spinDirection( 9+ Uspins*uCellNo,0)=  0.          ;spinDirection( 9+ Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection( 9+ Uspins*uCellNo,2)=  1./3;
            spinDirection(10+ Uspins*uCellNo,0)= - sqrt(2./3.) ;spinDirection(10+ Uspins*uCellNo,1)=    -sqrt(2.)/3. ;spinDirection(10+ Uspins*uCellNo,2)= 1./3;
            spinDirection(11+ Uspins*uCellNo,0)=  -sqrt(2./3.) ;spinDirection(11+ Uspins*uCellNo,1)=    -sqrt(2.)/3. ;spinDirection(11+ Uspins*uCellNo,2)= 1./3;
            
            spinDirection(12+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection(12+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection(12+ Uspins*uCellNo,2)=  1./3;
            spinDirection(13+ Uspins*uCellNo,0)=  0.          ;spinDirection(13+ Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection(13+ Uspins*uCellNo,2)=  1./3;
            spinDirection(14+ Uspins*uCellNo,0)=  0.          ;spinDirection(14+ Uspins*uCellNo,1)=    0.          ;spinDirection(14+ Uspins*uCellNo,2)= -1.;
            spinDirection(15+ Uspins*uCellNo,0)=  0.          ;spinDirection(15+ Uspins*uCellNo,1)=    0.          ;spinDirection(15+ Uspins*uCellNo,2)=  -1.;
            
            spinDirection(16+ Uspins*uCellNo,0)=  0.          ;spinDirection(16+ Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection(16+ Uspins*uCellNo,2)=  1./3;
            spinDirection(17+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection(17+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection(17+ Uspins*uCellNo,2)=  1./3;
            spinDirection(18+ Uspins*uCellNo,0)=  -sqrt(2./3.) ;spinDirection(18+ Uspins*uCellNo,1)=   - sqrt(2.)/3. ;spinDirection(18+ Uspins*uCellNo,2)= 1./3;
            spinDirection(19+ Uspins*uCellNo,0)=  0.          ;spinDirection(19+ Uspins*uCellNo,1)= 2.*sqrt(2.)/3. ;spinDirection(19+ Uspins*uCellNo,2)=  1./3;
            
            spinDirection(20+ Uspins*uCellNo,0)= - sqrt(2./3.) ;spinDirection(20+ Uspins*uCellNo,1)=   - sqrt(2.)/3. ;spinDirection(20+ Uspins*uCellNo,2)= 1./3;
            spinDirection(21+ Uspins*uCellNo,0)=  sqrt(2./3.) ;spinDirection(21+ Uspins*uCellNo,1)=   -sqrt(2.)/3. ;spinDirection(21+ Uspins*uCellNo,2)=  1./3;
            spinDirection(22+ Uspins*uCellNo,0)=  0.          ;spinDirection(22+ Uspins*uCellNo,1)=    0.          ;spinDirection(22+ Uspins*uCellNo,2)=  -1.;
            spinDirection(23+ Uspins*uCellNo,0)=  0.          ;spinDirection(23+ Uspins*uCellNo,1)=    0.          ;spinDirection(23+ Uspins*uCellNo,2)=  -1.;
        }
    }
   
    for (int i=0;i<Nspins;i++)
    {
        deformedSpinDirection.row(i)= (parameters.getFmatrix()*((spinDirection.row(i)).transpose()));
        deformedSpinDirection.row(i)= (deformedSpinDirection.row(i))/((deformedSpinDirection.row(i)).norm());
    }
    
   
    
}





void Crystal::findNN()
{   //This function returns an integer matrix NNlist, such that NNlist(n,m) is the m:th neighbor of the n:th particle
    NNlist= -Eigen::MatrixXi::Ones(lattice.rows(),6);//if no neighbor, the default value is -1
    double distanceFactor=1.1;
    Eigen::Vector3d pDirection=parameters.getPDirection();
    
    int NNnumber=0;
    int XYnumber=0;
    int Znumber=0;
    int Z1number=0;
    int Z2number=0;
    
   
    for(int i=0; i<lattice.rows() ;++i)
    {
        
        NNnumber=0;
        XYnumber=0;
        Znumber=0;
        Z1number=0;
        Z2number=0;
        
    
        //neighbors in other cells
        
        for(int j=0; j<lattice.rows();j++)
        {
            for (int a=-1; a<2; a++)
            {
                for (int b=-1; b<2; b++)
                {
                    for (int c=-1; c<2; c++)
                    {
                        
                        if(((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).norm()>0) &&
                           ((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).norm()<NNdist*distanceFactor))
                        {
                            //std::cout<<NNnumber<<"\n";
                            NNlist(i,NNnumber)=j;
                            NNnumber++;
                         
                            if(isEqual((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).dot(pDirection),0,0.001))
                            {//check if Z component differs
                                XYNNlist(i,XYnumber++)=j;
                            }else{
                                ZNNlist(i,Znumber++)=j;
                                
                                
                                if(  isEqual( std::abs((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).dot(pDirection)),
                                              (lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).norm(),0.001) ){
                                    
                                    Z1NNlist(i,Z1number++)=j;
                                }else{
                                    Z2NNlist(i,Z2number++)=j;
                                }
                                 
                                 
                            }
                        }
                    }
                }
            }
        }
    }
    /*
    for(int i=0; i<NNlist.rows(); i++)
    {
        std::cout<<i<<": "<<NNlist.row(i)<<"\n";
    }
     */
}














void Crystal::findNNN()
{   //This function returns an integer matrix NNNlist, such that NNNlist(n,m) is the m:th next nearest neighbor of the n:th particle
    findLattice();
    NNNlist= -Eigen::MatrixXi::Ones(lattice.rows(),12);//if no neighbor, the default value is -1
    double distanceFactor=1.1;
    //translation vectors of the crystal in order to find the neigbours in cells closeby
    //six different neighbors
    int NNNnumber=0;
    for(int i=0; i<lattice.rows() ;++i)
    {
        NNNnumber=0;
        //neighbors in other cells
        
        for(int j=0; j<lattice.rows();j++)
        {
            for (int a=-1; a<2; a++)
            {
                for (int b=-1; b<2; b++)
                {
                    for (int c=-1; c<2; c++)
                    {
                        if(((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).norm()>NNdist*distanceFactor) &&
                           ((lattice.row(j)+a*crystalSideTranslation.row(0)+b*crystalSideTranslation.row(1)+c*crystalSideTranslation.row(2)-lattice.row(i)).norm()<NNNdist*distanceFactor)  )
                        {
                            NNNlist(i,NNNnumber)=j;
                            NNNnumber++;
                        }
                    }
                }
            }
        }
    }
    /*
    for(int i=0; i<NNNlist.rows(); i++)
    {
        std::cout<<i<<": "<<NNNlist.row(i)<<"\n";
    }
     */
}














































void Crystal::saveLists()
{
    std::string path= parameters.getPathToDataFolder();
    //This void saves matrix data to files in the "path", it saves it with the name according to its type.
    std::string dataType[]={"NNlist", "NNNlist", "tetraSpinList","spinTetraList", "XYNNlist", "ZNNlist", "Z1NNlist", "Z2NNlist"};
    for(int type=0; type<8; type++)
    {
        
        std::string fileName=dataType[type]+"_"+std::to_string(Uspins)+"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())+".txt";
        
        
    
        if(dataType[type]=="XYNNlist" || dataType[type]=="ZNNlist" || dataType[type]=="Z1NNlist" || dataType[type]=="Z2NNlist")
        {
            fileName=dataType[type]+"_"+parameters.getPIndex()+"_"+std::to_string(Uspins)+"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())+".txt";
        }
       
        
        
        
        
        std::ofstream file(path+"crystalData/"+fileName);
        //exit if unable to create file
        if(!file){
            std::cerr<<"Crystal::saveToFile() file could not be opened"<<"\n";
            exit(EXIT_FAILURE);
        }
        
        if(dataType[type]=="NNlist")         file << NNlist;
        if(dataType[type]=="NNNlist")        file << NNNlist;
        if(dataType[type]=="spinTetraList")  file << spinTetraList;
        if(dataType[type]=="tetraSpinList")  file << tetraSpinList;
        if(dataType[type]=="XYNNlist")       file << XYNNlist;
        if(dataType[type]=="ZNNlist")        file << ZNNlist;
        if(dataType[type]=="Z1NNlist")        file << Z1NNlist;
        if(dataType[type]=="Z2NNlist")        file << Z2NNlist;
        
        
    }
    std::cout<<"Crystal Saved"<<"\n";
}






void Crystal::loadLists()
{
    //This void loads matrix data from files in the "path", into the class variables according to its type.
    
    
    std::string path= parameters.getPathToDataFolder();
    std::string dataType[]={"NNlist", "NNNlist", "tetraSpinList","spinTetraList", "XYNNlist", "ZNNlist", "Z1NNlist" , "Z2NNlist"};
    int rowsInSpecificType []={Uspins, Uspins,    Uspins/2,       Uspins,         Uspins,      Uspins, Uspins, Uspins};
    std::string line;
    for(int type=0; type<8; type++)
    {
        std::string filename=dataType[type]+"_"+std::to_string(Uspins)+"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())+".txt";
        
      
        if(dataType[type]=="XYNNlist" || dataType[type]=="ZNNlist" || dataType[type]=="Z1NNlist" || dataType[type]=="Z2NNlist")
        {
            filename=dataType[type]+"_"+parameters.getPIndex()+"_"+std::to_string(Uspins)+"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())+".txt";
        }
        
        
        
        std::ifstream file(path+"crystalData/"+filename);
        if(!file){
            std::cerr<<"Crystal::loadFromFile() Crystal file could not be opened"<<"\n";
            exit(EXIT_FAILURE);
        }
        for(int i=0; i<rowsInSpecificType[type]*parameters.getCellX()*parameters.getCellY()*parameters.getCellZ(); i++)
        {
            getline(file,line);
            std::istringstream is(line);
            double n;
            int j=0;
            while( is >> n )
            {
                if(dataType[type]=="NNlist")   NNlist(i,j)=n;
                if(dataType[type]=="NNNlist") NNNlist(i,j)=n;
                if(dataType[type]=="spinTetraList") spinTetraList(i,j)=n;
                if(dataType[type]=="tetraSpinList") tetraSpinList(i,j)=n;
                if(dataType[type]=="XYNNlist") XYNNlist(i,j)=n;
                if(dataType[type]=="ZNNlist")  ZNNlist(i,j)=n;
                if(dataType[type]=="ZNNlist")  Z1NNlist(i,j)=n;
                if(dataType[type]=="ZNNlist")  Z2NNlist(i,j)=n;
                j++;
            }
        }
    }
    std::cout<<"Crystal loaded"<<"\n";
}





bool Crystal::isEqual(double a, double b, double tol)
{
    return std::abs(a-b)<tol;
}















