//
//  Ewald.cpp
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-10-23.
//


//The interaction of particle i and particle j has a contribution to the total potential energy. This contribution can be written in terms of a matrix. This matrix is called the "Ewald" matrix,
//and the contribution is V=1/2*∑V_ij 





#include "Ewald.h"



Ewald::Ewald(Parameters parameters, Crystal crystal)
:
crystal(crystal),
parameters(parameters),
ewaldMatrixFound(false)
{
    alpha=(double) sqrt((16*M_PI)/((double)crystal.getLattice().rows()));
    ewaldMatrix=Eigen::MatrixXd(crystal.getLattice().rows(),crystal.getLattice().rows()); //initialize with the number of particles
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    
    if(myrank==0)
    {
        if(parameters.getLoadEwald())
        {
            load();
        }else
        {
            findEwaldMatrix();
            save();
        }
    }
    ewaldMatrix=ewaldMatrix*parameters.getD();
    
    
    
    MPI_Bcast(ewaldMatrix.data(), (int) ewaldMatrix.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    ewaldMatrixFound=true;
    
    
}






Eigen::MatrixXd Ewald::getEwaldMatrix()
{
    if(!ewaldMatrixFound)
    {
        std::cerr<<"Error: Ewald matrix not initialized\n";
        exit(EXIT_FAILURE);
    }
    //prefactor using input D and NNdist already in constructor
    return ewaldMatrix;
    
}







double Ewald::A(double r)
{
    return erfc(sqrt(alpha)*r)*pow(r,-3)+2*sqrt(alpha/M_PI)*exp(-alpha*pow(r,2))*pow(r,-2);
}

double Ewald::B(double r)
{
    return 3*erfc(sqrt(alpha)*r)*pow(r,-5)+2*sqrt(alpha/M_PI)*(2*alpha+3*pow(r,-2))*exp(-alpha*pow(r,2))*pow(r,-2);
}





void Ewald::findEwaldMatrix()
{
    ewaldMatrix=Eigen::MatrixXd(crystal.getLattice().rows(),crystal.getLattice().rows());
    double V = crystal.getCellVolume();
    double diagonalElement=0;
    double diagonalElementOld=0;
    //double exactDiagonalElement=-pow(parameters.getNNdist(),3)*(4*M_PI)/(3*V);
    double diagonalDifference=1;
    double tolerance=pow(10,-14);
    //find number of mirrors
    int mirrors=0;
    bool mirrorsExceeded=false;
    while(abs(diagonalDifference)>tolerance && !mirrorsExceeded)
    {
        diagonalElementOld=diagonalElement;
        diagonalElement=pow(crystal.getNNdist(),3)*((4*M_PI/V)*reciprocalSum(0,0, mirrors) - selfSum(0, 0) +directSum(0, 0, mirrors));
        diagonalDifference=diagonalElement-diagonalElementOld;
        mirrors++;
        if(mirrors>parameters.getMaxEwaldMirrors())
        {
            std::cerr<<"Error Ewald: number of mirrors too high, stop at "<< parameters.getMaxEwaldMirrors()<<"\n";
            mirrorsExceeded=true;
            mirrors--;
        }
    }
    /*
    std::cout<<"Ewald: Number of mirrors: "<<mirrors<<", diagonalDifference: "<<std::setprecision(16)<<diagonalDifference<<"\n";
    std::cout<<"Ewald: My diagonal element       = "<<std::setprecision(16)<<diagonalElement<<"\n";
    std::cout<<"Ewald: The exact diagonal element= "<<std::setprecision(16)<< exactDiagonalElement<<"\n";
    */
    //Find elements
    for(int i=0; i<crystal.getLattice().rows(); i++)
    {
        for(int j=i; j<crystal.getLattice().rows(); j++)
        {
            if(i!=j)
            {
                ewaldMatrix(i,j)= pow(crystal.getNNdist(),3)*((4*M_PI/V)*reciprocalSum(i,j, mirrors) - selfSum(i, j) +directSum(i, j, mirrors));
                ewaldMatrix(j,i)= ewaldMatrix(i,j);
            }else
            {
                ewaldMatrix(i,i)=0;
            }
        }
    }
    
    

    if(parameters.getAddEwaldDemagTerm())
    {//add Demag extra energy term to Ewald matrix.
        for(int i=0; i<crystal.getLattice().rows(); i++)
        {
            for(int j=i; j<crystal.getLattice().rows(); j++)
            {
                Eigen::Vector3d mui =     crystal.getSpinDirection().row(i);
                Eigen::Vector3d muj =     crystal.getSpinDirection().row(j);
                
                if(i!=j)
                {
                    ewaldMatrix(i,j) += pow(crystal.getNNdist(),3)*((4*M_PI/V))*mui.dot(parameters.getDemagMatrix()*muj);
                    ewaldMatrix(j,i) += pow(crystal.getNNdist(),3)*((4*M_PI/V))*mui.dot(parameters.getDemagMatrix()*muj);
                }else
                {
                    ewaldMatrix(i,i)=0;
                }
            }
        }
    }
    
    ewaldMatrixFound=true;



}



















double Ewald::directSum(int i, int j, int numberOfMirrors)
{
    //This function should make the direct Ewald summation with numberOfMirrors mirrors, that is, a cube of side length 1+2*numberOfMirrors (thats how we keep the middle in the middle)
    Eigen::Vector3d mui =     crystal.getSpinDirection().row(i);
    Eigen::Vector3d muj =     crystal.getSpinDirection().row(j);
    Eigen::MatrixXd lattice = crystal.getLattice();
    Eigen::Vector3d rijn = lattice.row(j)-lattice.row(i);
    Eigen::MatrixXd crystalSideTranslation = crystal.getCrystalSideTranslation();
    
    
    
    double sum = 0;
    double muimuj= (mui.transpose())*muj;
    double muirijn=(mui.transpose())*rijn;
    double mujrijn=(muj.transpose())*rijn;
    for(int nz=-numberOfMirrors; nz<=numberOfMirrors; nz++)
    {
        for(int ny=-numberOfMirrors; ny<=numberOfMirrors; ny++)
        {
            for(int nx=-numberOfMirrors; nx<=numberOfMirrors; nx++)
            {
                if(!(i==j && nx==0 && ny==0 && nz==0))
                {
                    rijn=lattice.row(j)-lattice.row(i)-nx*crystalSideTranslation.row(0)-ny*crystalSideTranslation.row(1)-nz*crystalSideTranslation.row(2);
                    muimuj= (mui.transpose())*muj;
                    muirijn=(mui.transpose())*rijn;
                    mujrijn=(muj.transpose())*rijn;
                    sum += muimuj*A(rijn.norm())-muirijn*mujrijn*B(rijn.norm());
                }
            }
        }
    }
    
    return sum;
}







double Ewald::reciprocalSum(int i, int j, int numberOfMirrors)
{
    Eigen::Vector3d mui =     crystal.getSpinDirection().row(i);
    Eigen::Vector3d muj =     crystal.getSpinDirection().row(j);
    Eigen::MatrixXd lattice = crystal.getLattice();
    Eigen::RowVector3d rij = lattice.row(j)-lattice.row(i);
    Eigen::Vector3d k;
    k<<1,1,1; //initialize
    
    double muik = mui.transpose()*k;
    double mujk = muj.transpose()*k;
    double rijk = rij*k;
    double ksqr= pow(k.norm(),2);
    double sum=0;
    

    
    
    
    Eigen::Vector3d L1=(crystal.getCrystalSideTranslation().row(0)).transpose();
    Eigen::Vector3d L2=(crystal.getCrystalSideTranslation().row(1)).transpose();
    Eigen::Vector3d L3=(crystal.getCrystalSideTranslation().row(2)).transpose();
    Eigen::Vector3d kx= (2*M_PI*L2.cross(L3))/(L1.dot(L2.cross(L3)));
    Eigen::Vector3d ky= (2*M_PI*L3.cross(L1))/(L1.dot(L2.cross(L3)));
    Eigen::Vector3d kz= (2*M_PI*L1.cross(L2))/(L1.dot(L2.cross(L3)));
    
    
    for(int nz=-numberOfMirrors; nz<=numberOfMirrors; nz++)
    {
        for(int ny=-numberOfMirrors; ny<=numberOfMirrors; ny++)
        {
            for(int nx=-numberOfMirrors; nx<=numberOfMirrors; nx++)
            {
                if(!(nx==0 && ny==0 && nz==0))
                {
                    
                    //old k
                    //k<< nx*(2*M_PI/crystal.getCellSideLengthX()), ny*(2*M_PI/crystal.getCellSideLengthY()), nz*(2*M_PI/crystal.getCellSideLengthZ());
                    
                    
                    k= nx*kx+ny*ky+nz*kz;
                    
                    
                    muik = mui.transpose()*k;
                    mujk = muj.transpose()*k;
                    rijk = rij*k;
                    ksqr = pow(k.norm(),2);
                    sum += ((muik*mujk)/(double)ksqr)*exp(-ksqr/(4*alpha))*cos(rijk);
                }
            }
        }
    }
    return sum;
}



double Ewald::selfSum(int i, int j)
{
    if (i==j)
    {
        return ((4*M_PI)/(double)3)*pow(alpha/M_PI,3/(double)2);
    }
    return 0;
}






void Ewald::save()
{
    //This void saves the matrix data to the file Ewald_x_x.txt
    
    
    std::string path = parameters.getPathToDataFolder()+"ewaldData/";
    std::string filePathAndName=path+"Ewald"
            +"_"+std::to_string(parameters.getCellType())
            +"_"+parameters.getPIndex()
            +"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())
            +"_"+std::to_string(parameters.getFmatrix()(0,0))+"_"+std::to_string(parameters.getFmatrix()(0,1))+"_"+std::to_string(parameters.getFmatrix()(0,2))
            +"_"+std::to_string(parameters.getFmatrix()(1,0))+"_"+std::to_string(parameters.getFmatrix()(1,1))+"_"+std::to_string(parameters.getFmatrix()(1,2))
            +"_"+std::to_string(parameters.getFmatrix()(2,0))+"_"+std::to_string(parameters.getFmatrix()(2,1))+"_"+std::to_string(parameters.getFmatrix()(2,2));
    
    if(parameters.getAddEwaldDemagTerm())
    {
        filePathAndName=filePathAndName+"_DemagAdded_"
        +"_"+std::to_string(parameters.getDemagMatrix()(0,0))+"_"+std::to_string(parameters.getDemagMatrix()(0,1))+"_"+std::to_string(parameters.getDemagMatrix()(0,2))
        +"_"+std::to_string(parameters.getDemagMatrix()(1,0))+"_"+std::to_string(parameters.getDemagMatrix()(1,1))+"_"+std::to_string(parameters.getDemagMatrix()(1,2))
        +"_"+std::to_string(parameters.getDemagMatrix()(2,0))+"_"+std::to_string(parameters.getDemagMatrix()(2,1))+"_"+std::to_string(parameters.getDemagMatrix()(2,2));
    }
    
    filePathAndName=filePathAndName+".txt";
    std::ofstream file(filePathAndName);
   
    //exit if unable to create file
    if(!file){
        std::cerr<<"Error: Ewald file could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file << std::setprecision(16)<<  ewaldMatrix;
    std::cout<<"Ewald Matrix Saved"<<"\n";
}



void Ewald::load()
{
    //This void the Ewald Matrix from files in the "path", into the class variable Ewald.
    
    
    std::string path= parameters.getPathToDataFolder()+"ewaldData/";
    std::string filename="Ewald_"
        +std::to_string(parameters.getCellType())
        +"_"+parameters.getPIndex()
        +"_"+std::to_string(parameters.getCellX())+"_"+std::to_string(parameters.getCellY())+"_"+std::to_string(parameters.getCellZ())
        +"_"+std::to_string(parameters.getFmatrix()(0,0))+"_"+std::to_string(parameters.getFmatrix()(0,1))+"_"+std::to_string(parameters.getFmatrix()(0,2))
        +"_"+std::to_string(parameters.getFmatrix()(1,0))+"_"+std::to_string(parameters.getFmatrix()(1,1))+"_"+std::to_string(parameters.getFmatrix()(1,2))
        +"_"+std::to_string(parameters.getFmatrix()(2,0))+"_"+std::to_string(parameters.getFmatrix()(2,1))+"_"+std::to_string(parameters.getFmatrix()(2,2));
    
    if(parameters.getAddEwaldDemagTerm())
    {
        filename=filename+"_DemagAdded_"
        +"_"+std::to_string(parameters.getDemagMatrix()(0,0))+"_"+std::to_string(parameters.getDemagMatrix()(0,1))+"_"+std::to_string(parameters.getDemagMatrix()(0,2))
        +"_"+std::to_string(parameters.getDemagMatrix()(1,0))+"_"+std::to_string(parameters.getDemagMatrix()(1,1))+"_"+std::to_string(parameters.getDemagMatrix()(1,2))
        +"_"+std::to_string(parameters.getDemagMatrix()(2,0))+"_"+std::to_string(parameters.getDemagMatrix()(2,1))+"_"+std::to_string(parameters.getDemagMatrix()(2,2));
    }
    
    filename=filename+".txt";
    
    std::ifstream file(path+filename);
    if(!file){
        std::cerr<<"Ewald matrix file could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line;
    
    
    for(int i=0; i<crystal.getLattice().rows() ; i++)
    {
            getline(file,line);
            std::istringstream is(line);
            double n;
            int j=0;
        while( is >> std::setprecision(16) >> n )
            {
                ewaldMatrix(i,j)=n;
                j++;
            }
    }
    ewaldMatrixFound=true;
    std::cout<<"Ewald matrix loaded"<<"\n";
}



