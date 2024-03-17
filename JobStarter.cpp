//
//  main.cpp
//  NNMCPyrochlore
//
//  Created by Richard Edberg on 2018-03-02.
//  Copyright © 2017 Richard Edberg. All rights reserved.
//

#include "JobStarter.h"
#include <string>

//========================================================================

JobStarter::JobStarter(int argc, char *argv[], std::string pathToDataFolder)
:pathToDataFolder(pathToDataFolder)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    std::cout<<"number of arguments= "<<argc<<"\n";
    if(argc==1)
    {
        argFilename="parametersMPI.txt";
    }else
    {
        argFilename=argv[1];
    }
    
}

void JobStarter::start()
{
    //startParallelTempering();
    //startTrivialMPI();
    startMicrocanonical();

}










void JobStarter::startTrivialMPI()
{
    //load
    /******************************************************************************************************************/
   
    std::string parametersFilename= argFilename;
    std::cout<<parametersFilename<<"\n";
// std::string parametersFilename="parametersMPI.txt";
    Parameters parameters(pathToDataFolder, parametersFilename);
    if(myrank==0) std::cout<<parameters.print()<<"\n";
    Crystal crystal(parameters);
    if(parameters.getLoadCrystal())
    {
        crystal.loadLists();
    }else
    {
        crystal.findLists();
        crystal.saveLists();
    }
    Ewald* ewald= new Ewald(parameters, crystal);
    
   
    
    
    //Trivial MPI
    /******************************************************************************************************************/
    Eigen::VectorXd spinOrientation = Eigen::VectorXd::Ones(crystal.getLattice().rows());
    
    
    
    
    //read configuration (optional)-------------------------------------------------------------------------------------------------------------------------------------a
    if(parameters.getReadConfiguration())
    {
        std::ifstream inputFile;
        std::string s_reader, t_reader;
        inputFile.open((parameters.getPathToDataFolder()+parameters.getOutputFilename()+"_spinConfiguration_"+std::to_string(myrank)+".txt").c_str());
        if (!inputFile.is_open()){
            std::cerr<<("input spin configuration not found, exiting\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int i=0;
        while (getline(inputFile,s_reader))
        {
            std::istringstream r_pick(s_reader);
            if(i>=crystal.getNspins())
            {
                std::cerr<<"JobStarter::Error with experimental read spinConfiguration dimensions.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }else{
                std::string tempString;
                r_pick>>tempString;
                spinOrientation(i)=std::stod(tempString);
                }
            i++;
        }
        inputFile.close();
        std::cout<<"Configuration "<< myrank <<" loaded. "<< std::to_string(i) <<"spins.\n";
    }
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------a'
    
    
    
    
    
    
    TrivialPar mcSimulation(&parameters, &crystal, ewald, &spinOrientation);
    FileWriter fileWriter(&parameters);
    int numberOfdataCols = (int) mcSimulation.getData().cols();
    Eigen::MatrixXd data = Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins(),numberOfdataCols);
    
    
  
    if (myrank == 0)
    {
        Eigen::MatrixXd totalData= Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins()*nproc,numberOfdataCols);
        Eigen::MatrixXd totalSumCorrSiSj = Eigen::MatrixXd::Zero(crystal.getNspins(),crystal.getNspins());
        
       
        clock_t begin = std::clock();
        mcSimulation.run(myrank,nproc);
        delete ewald;
        data=mcSimulation.getData();
        totalSumCorrSiSj+=mcSimulation.sumCorrSiSj;
        
        
        totalData.block(0,0,data.rows(),data.cols())=data;
        
        
        for(int ranknumber=1; ranknumber<nproc; ++ranknumber)
        {
            MPI_Recv(data.data(), (int)data.size(), MPI_DOUBLE, ranknumber, 4711, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(mcSimulation.sumCorrSiSj.data(), (int) mcSimulation.sumCorrSiSj.size(), MPI_INT, ranknumber, 4712, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalSumCorrSiSj+=mcSimulation.sumCorrSiSj;
            totalData.block(ranknumber*data.rows(),0,data.rows(),data.cols())=data;
        }
        
        Eigen::MatrixXd sortedTotalData = Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins()*nproc,numberOfdataCols);
        
        
        int Z= parameters.getDataPoints()*parameters.getBins();
        int dataSortedRow=0;
        for(int I=0; I<Z; ++I)
        {
            for(int R=0;  R< nproc; ++R)
            {
                sortedTotalData.row(dataSortedRow++)=totalData.row(I+Z*R);
            }
        }
    
        
        
        
        
        totalSumCorrSiSj= totalSumCorrSiSj/(nproc*((parameters.getMCsteps()/parameters.getDataPoints())/parameters.getBins())*parameters.getBins()*parameters.getDataPoints());
        
        
        clock_t end = std::clock();
        long long int elapsed_secs = (long long int) double(end - begin)/CLOCKS_PER_SEC;
        fileWriter.writeDataMatrix(&sortedTotalData,elapsed_secs, nproc);
        if(parameters.getWriteCorrelationMatrix())
        {
            fileWriter.append("corrSiSj:", &totalSumCorrSiSj);
        }
        std::cout<<"Elapsed time: "<<elapsed_secs/3600<<"h:"<<(elapsed_secs%3600)/60<<"m:"<<elapsed_secs%60<<"s\n";
    }
    else
    {
        
        mcSimulation.run(myrank,nproc);
        delete ewald;
        
        MPI_Send(mcSimulation.getData().data(), (int)mcSimulation.getData().size(), MPI_DOUBLE, 0, 4711, MPI_COMM_WORLD);
        
        
        
        MPI_Send(mcSimulation.sumCorrSiSj.data(), (int) mcSimulation.sumCorrSiSj.size(), MPI_INT, 0, 4712, MPI_COMM_WORLD);
    }
    
    
    //-----write configuration to file.
    std::ofstream file((parameters.getPathToDataFolder()+parameters.getOutputFilename()+"_spinConfiguration_"+std::to_string(myrank)+".txt").c_str());
    //exit if unable to create file
    if(!file){
        std::cerr<<"FileWriter:: File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file<<spinOrientation;
    //-----end write configuration to file

    MPI_Finalize();
  

}




void JobStarter::startMicrocanonical()
{
    //load
    /******************************************************************************************************************/
   
    std::string parametersFilename= argFilename;
    std::cout<<parametersFilename<<"\n";
// std::string parametersFilename="parametersMPI.txt";
    Parameters parameters(pathToDataFolder, parametersFilename);
    if(myrank==0) std::cout<<parameters.print()<<"\n";
    Crystal crystal(parameters);
    if(parameters.getLoadCrystal())
    {
        crystal.loadLists();
    }else
    {
        crystal.findLists();
        crystal.saveLists();
    }
    Ewald* ewald= new Ewald(parameters, crystal);
    
   
    
    
    //Trivial MPI
    /******************************************************************************************************************/
    Eigen::VectorXd spinOrientation = Eigen::VectorXd::Ones(crystal.getLattice().rows());
    
    
    
    
    //read configuration (optional)-------------------------------------------------------------------------------------------------------------------------------------a
    if(parameters.getReadConfiguration())
    {
        std::ifstream inputFile;
        std::string s_reader, t_reader;
        inputFile.open((parameters.getPathToDataFolder()+parameters.getOutputFilename()+"_spinConfiguration_"+std::to_string(myrank)+".txt").c_str());
        if (!inputFile.is_open()){
            std::cerr<<("input spin configuration not found, exiting\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int i=0;
        while (getline(inputFile,s_reader))
        {
            std::istringstream r_pick(s_reader);
            if(i>=crystal.getNspins())
            {
                std::cerr<<"JobStarter::Error with experimental read spinConfiguration dimensions.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }else{
                std::string tempString;
                r_pick>>tempString;
                spinOrientation(i)=std::stod(tempString);
                }
            i++;
        }
        inputFile.close();
        std::cout<<"Configuration "<< myrank <<" loaded. "<< std::to_string(i) <<"spins.\n";
    }
    
    
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------a'
    
    
    
    
    
    
    Microcanonical mcSimulation(&parameters, &crystal, ewald, &spinOrientation);
    FileWriter fileWriter(&parameters);
    int numberOfdataCols = (int) mcSimulation.getData().cols();
    Eigen::MatrixXd data = Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins(),numberOfdataCols);
    
    
  
    if (myrank == 0)
    {
        Eigen::MatrixXd totalData= Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins()*nproc,numberOfdataCols);
        Eigen::MatrixXd totalSumCorrSiSj = Eigen::MatrixXd::Zero(crystal.getNspins(),crystal.getNspins());
        
       
        clock_t begin = std::clock();
        mcSimulation.run(myrank,nproc);
        delete ewald;
        data=mcSimulation.getData();
        totalSumCorrSiSj+=mcSimulation.sumCorrSiSj;
        
        
        totalData.block(0,0,data.rows(),data.cols())=data;
        
        
        for(int ranknumber=1; ranknumber<nproc; ++ranknumber)
        {
            MPI_Recv(data.data(), (int)data.size(), MPI_DOUBLE, ranknumber, 4711, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(mcSimulation.sumCorrSiSj.data(), (int) mcSimulation.sumCorrSiSj.size(), MPI_INT, ranknumber, 4712, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            totalSumCorrSiSj+=mcSimulation.sumCorrSiSj;
            totalData.block(ranknumber*data.rows(),0,data.rows(),data.cols())=data;
        }
        
        Eigen::MatrixXd sortedTotalData = Eigen::MatrixXd::Zero(parameters.getDataPoints()*parameters.getBins()*nproc,numberOfdataCols);
        
        
        int Z= parameters.getDataPoints()*parameters.getBins();
        int dataSortedRow=0;
        for(int I=0; I<Z; ++I)
        {
            for(int R=0;  R< nproc; ++R)
            {
                sortedTotalData.row(dataSortedRow++)=totalData.row(I+Z*R);
            }
        }
    
        
        
        
        
        totalSumCorrSiSj= totalSumCorrSiSj/(nproc*((parameters.getMCsteps()/parameters.getDataPoints())/parameters.getBins())*parameters.getBins()*parameters.getDataPoints());
        
        
        clock_t end = std::clock();
        long long int elapsed_secs = (long long int) double(end - begin)/CLOCKS_PER_SEC;
        fileWriter.writeDataMatrix(&sortedTotalData,elapsed_secs, nproc);
        if(parameters.getWriteCorrelationMatrix())
        {
            fileWriter.append("corrSiSj:", &totalSumCorrSiSj);
        }
        std::cout<<"Elapsed time: "<<elapsed_secs/3600<<"h:"<<(elapsed_secs%3600)/60<<"m:"<<elapsed_secs%60<<"s\n";
    }
    else
    {
        
        mcSimulation.run(myrank,nproc);
        delete ewald;
        
        MPI_Send(mcSimulation.getData().data(), (int)mcSimulation.getData().size(), MPI_DOUBLE, 0, 4711, MPI_COMM_WORLD);
        
        
        
        MPI_Send(mcSimulation.sumCorrSiSj.data(), (int) mcSimulation.sumCorrSiSj.size(), MPI_INT, 0, 4712, MPI_COMM_WORLD);
    }
    
    
    //-----write configuration to file.
    std::ofstream file((parameters.getPathToDataFolder()+parameters.getOutputFilename()+"_spinConfiguration_"+std::to_string(myrank)+".txt").c_str());
    //exit if unable to create file
    if(!file){
        std::cerr<<"FileWriter:: File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file<<spinOrientation;
    //-----end write configuration to file

    MPI_Finalize();
  

}





















void JobStarter::startParallelTempering()
{
    
    
    /*
    std::string parametersFilename="parametersMPI.txt";
    Parameters parameters(pathToDataFolder, parametersFilename);
    if(myrank==0) parameters.print();
    
    
    
    //-----------------------------------------------------------------------------------------------
    Crystal crystal(parameters);
    if(parameters.getLoadCrystal())
    {
        crystal.loadLists();
    }else
    {
        crystal.findLists();
        crystal.saveLists();
    }
    //---------------------------------------------------------------------------------------------------
    Ewald ewald(parameters, crystal);
    
    //-----------------------------------------------------------------------------------------------------------
    
    
    
    Eigen::VectorXd spinOrientation = Eigen::VectorXd::Ones(crystal.getLattice().rows());
    ParallelTempering parallelTempering(parameters, crystal, ewald, &spinOrientation, myrank, nproc);
    FileWriter fileWriter(parameters);
   
    //MPI---------------------------------------------------------------------------------------------------
    
   
    
    
    parallelTempering.run();
    MPI_Finalize();
     
     */
    
}

















