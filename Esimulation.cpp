//
//  Esimulation.cpp
//  exactPyrochloreDipole
//
//  Created by Richard Edberg on 2017-11-23.
//  Copyright © 2017 Richard Edberg. All rights reserved.
//


#include "Esimulation.h"



Esimulation::Esimulation(Parameters parameters, Crystal crystal, Ewald ewald)
:parameters(parameters),
crystal(crystal),
ewald(ewald)
{
    Jnnxy=parameters.getJnnxy();
    
    Jnnz1=parameters.getJnnz1();
    Jnnz2=parameters.getJnnz2();
    
    
    Jnnn=parameters.getJnnn();
    muNorm=parameters.getMuNorm();
    fieldDirection= parameters.getFieldDirection();
    
    
    lattice=crystal.getLattice();
    spinDirection=crystal.getSpinDirection();
    NNlist=crystal.getNNlist();
    NNNlist=crystal.getNNNlist();
    XYNNlist=crystal.getXYNNlist();
    ZNNlist=crystal.getZNNlist();
    Z1NNlist=crystal.getZ1NNlist();
    Z2NNlist=crystal.getZ2NNlist();
    ewaldMatrix=ewald.getEwaldMatrix();
    
    int datapoints=parameters.getDataPoints();
    data= Eigen::MatrixXd(datapoints,18);

  
}



double Esimulation::Energy()
{//This void calculates the energy of the system
    int Nspins=16;
    double E=0;
    Eigen::RowVector3d a;
    Eigen::Vector3d b;
    for(int i=0; i<Nspins; ++i)
    {
        //NNXY------------------------------------------------------------------------------
        for(int j=0; XYNNlist(i,j)>-1; j++)
        {
            a=spinDirection.row(i)*spinOrientation(i);
            b=spinDirection.row(XYNNlist(i,j)).transpose()*spinOrientation(XYNNlist(i,j));
            E += 0.5*Jnnxy*a*b;//0.5 to remove double counting
        }
        
        
        //NNZ1--------------------------------------------------------------------------------
        for(int j=0; Z1NNlist(i,j)>-1 ; j++)
        {
            a=spinDirection.row(i)*spinOrientation(i);
            b=spinDirection.row(Z1NNlist(i,j)).transpose()*spinOrientation(Z1NNlist(i,j));
            E += 0.5*Jnnz1*a*b;//0.5 to remove double counting
        }
        //NNZ2--------------------------------------------------------------------------------
        for(int j=0; Z2NNlist(i,j)>-1 ; j++)
        {
            a=spinDirection.row(i)*spinOrientation(i);
            b=spinDirection.row(Z2NNlist(i,j)).transpose()*spinOrientation(Z2NNlist(i,j));
            E += 0.5*Jnnz2*a*b;//0.5 to remove double counting
        }
         
        
        //Dipole------------------------------------------------------------------------------
        for(int j=0; j<lattice.rows(); ++j)
        {
            E+= 0.5*ewaldMatrix(i,j)*spinOrientation(i)*spinOrientation(j);//0.5 to remove double counting
        }
        //Field-------------------------------------------------------------------------------
        E += -spinOrientation(i)*muNorm*B*(spinDirection.row(i)).dot(fieldDirection);
    }
    
    
    return E;
}













Eigen::MatrixXd Esimulation::run()
{
    
    
    //looping parameters
   
    double Tstart= parameters.getTstart();
    double Tend=   parameters.getTend();
    double Bstart= parameters.getBstart();
    double Bend=   parameters.getBend();
   
    
    
    T=Tstart;
    B=Bstart;
    Jnnz=JnnzStart;
  
    
    
     int datapoints=parameters.getDataPoints();
    //initializing crystal quantities
    spinOrientation = - Eigen::RowVectorXi::Ones(lattice.rows());
   
    
    //partition function quantities
    double E=0;
    Eigen::RowVector3d M;
    M<<0,0,0;
    double boltzmann=0;
    double Z=0;
    double ZE=0;
    double ZEsqr=0;
    double Cv=0;
    double ZM=0;
    double ZMx=0;
    double ZMy=0;
    double ZMz=0;
    double ZMxsqr=0;
    double ZMysqr=0;
    double ZMzsqr=0;
    double ZMsqr=0;
    double chi=0;
    double ZMalongField=0;
    double ZMalongFieldSqr=0;
    
    
   
    for(int k=0; k<datapoints; k++)
    {
        boltzmann=0;
        Z=0;
        ZE=0;
        ZEsqr=0;
        Cv=0;
        ZM=0;
        ZMx=0;
        ZMy=0;
        ZMz=0;
        ZMxsqr=0;
        ZMysqr=0;
        ZMzsqr=0;
        ZMsqr=0;
        chi=0;
        ZMalongField=0;
        ZMalongFieldSqr=0;
       
        

        B=(Bend-Bstart)*k/(double)datapoints+Bstart;
        T=(Tend-Tstart)*k/(double)datapoints+Tstart;
        Jnnz=(JnnzEnd-JnnzStart)*k/(double)datapoints+JnnzStart;
        
        
        std::cout<<"Progress: "<<std::setprecision(2) << 100*k/datapoints <<"%\n";
        
       
        for(int i=0; i<pow(2,lattice.rows()); i++)
        {
            
            E = Energy();
            
            
            M << 0,0,0;
            for(int n=0; n<lattice.rows(); n++)
            {
                M += spinDirection.row(n)*spinOrientation(n);
            }
            
     
          
            boltzmann=exp(-E/T); //kb=1
            Z     +=                   boltzmann;
            ZE    +=                 E*boltzmann;
            ZEsqr +=          pow(E,2)*boltzmann;
            ZM    +=          M.norm()*boltzmann;
            ZMsqr +=   M.squaredNorm()*boltzmann;
            ZMx +=                  M(0)*boltzmann;// zero average checked
            ZMy +=                  M(1)*boltzmann;
            ZMz +=                  M(2)*boltzmann;
            ZMxsqr +=        pow(M(0),2)*boltzmann;// sum of these become ZMsqr checked
            ZMysqr +=        pow(M(1),2)*boltzmann;
            ZMzsqr +=        pow(M(2),2)*boltzmann;
            
            ZMalongField += M.dot(fieldDirection)*boltzmann;
            ZMalongFieldSqr += pow(M.dot(fieldDirection),2)*boltzmann;
            
            spinOrientation = increase(spinOrientation, 0);
            
        }
        
        Cv =(ZEsqr/Z-pow(ZE/Z,2))/pow(T,2);
        chi=(ZMsqr/Z)/T; //This is equal to chix+chiy+chiz, because <Mx^2>+<My^2>+<Mz^2>=<Mx^2+My^2+Mz^2>
        
        
        
        data(k,0)=T;
        data(k,1)=B;
        data(k,2)=(ZE/Z);
        data(k,3)=Cv;
        data(k,4)=(ZM/Z);
        data(k,5)=chi;
        data(k,6)=(ZMx/Z);
        data(k,7)=(ZMy/Z);
        data(k,8)=(ZMz/Z);
        data(k,9)=std::numeric_limits<double>::quiet_NaN();
        data(k,10)=std::numeric_limits<double>::quiet_NaN();
        data(k,11)=Jnnz;
        data(k,12)=std::numeric_limits<double>::quiet_NaN();;
        data(k,13)=std::numeric_limits<double>::quiet_NaN();;
        data(k,14)=std::numeric_limits<double>::quiet_NaN();;
        data(k,15)=ZMalongField/Z;
        data(k,16)=ZMalongFieldSqr/Z;
        data(k,17)= log(Z)+(ZE/Z)/T;
    }
    return data;
}










Eigen::RowVectorXi Esimulation::increase(Eigen::RowVectorXi spinOrientation, int i)
{
    //This recursive function changes the state of a the matrix as if it was a binary number
    
    if(i>=spinOrientation.size()) return spinOrientation;
    if (spinOrientation(i)<0)
    {
        spinOrientation(i)=1;
        return spinOrientation;
    }
    
    spinOrientation(i)=-1;
    return  increase(spinOrientation,i+1);
   
}




/* old, do not use.
void Esimulation::saveToFile()
{
    //This void saves the matrix data to the file data.txt
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::string timeString=std::to_string(now->tm_year-100)+std::to_string(now->tm_mon+1)+std::to_string(now->tm_mday)+std::to_string(now->tm_hour)+std::to_string(now->tm_min);
    
    //std::ofstream file(path+"exactData"+timeString+".txt"); //open in constructor, ios::out standard, also ios::app
    std::ofstream file(parameters.getPathToDataFolder()+parameters.getOutputFilename()+".txt");
    //exit if unable to create file
    if(!file){
        std::cerr<<"File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file<<parameters.print()<<"\n";
    file<<"SimulationResult:\n";
    file<<"data labels not added here";
    file <<std::setprecision(16)<< data;
    std::cout<<"File Saved"<<"\n";
    
}
*/


Eigen::MatrixXd Esimulation::getDataMatrix()
{
    
    return data;
}



























