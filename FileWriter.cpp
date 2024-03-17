//
//  FileWriter.cpp
//  PyrochloreG
//
//  Created by Richard Edberg on 2018-02-02.
//

#include "FileWriter.h"



FileWriter::FileWriter(Parameters *parameters)
:
parameters(parameters)
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    timeString=std::to_string(now->tm_year-100)+"."+std::to_string(now->tm_mon+1)+"."+std::to_string(now->tm_mday)+"."+std::to_string(now->tm_hour)+"."+std::to_string(now->tm_min)+"."+std::to_string(now->tm_sec);
}

void FileWriter::append(std::string delimeterLine, Eigen::MatrixXd *data)
{
    std::ofstream outfile;
    
    std::string filePathAndName=(*parameters).getPathToDataFolder()+(*parameters).getOutputFilename()+".txt";
    if((*parameters).getAddDateStamp()) filePathAndName=(*parameters).getPathToDataFolder()+(*parameters).getOutputFilename()+timeString+".txt";
    outfile.open(filePathAndName, std::ios_base::app);
    outfile << "\n"+delimeterLine+"\n";
    outfile << std::setprecision(16)<<(*data);
    
}




void FileWriter::writeDataMatrix(Eigen::MatrixXd *data, long long int elapsed_secs, int nproc)
{//This void saves the matrix data to a file
    
    std::string filePathAndName=(*parameters).getPathToDataFolder()+(*parameters).getOutputFilename()+".txt";
    if((*parameters).getAddDateStamp()) filePathAndName=(*parameters).getPathToDataFolder()+(*parameters).getOutputFilename()+timeString+".txt";
    std::ofstream file(filePathAndName);
    
    //exit if unable to create file
    if(!file){
        std::cerr<<"FileWriter:: File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file<<"Elapsed time: "<<elapsed_secs/3600<<"h:"<<(elapsed_secs%3600)/60<<"m:"<<elapsed_secs%60<<"s\n";
    file<<"MCsteps= "<<(*parameters).getMCsteps()*nproc<<"\n";
    file<<"bins= "<<(*parameters).getBins()*nproc<<"\n";
    file<<"nproc= "<<nproc<<"\n";
    file<<(*parameters).print();
    file<<"SimulationResult:\n";
    file<<"                 T                      B                   Eavg                EsqrAvg                   Mavg                MsqrAvg                     Mxavg                     Myavg                     Mzavg          singleAccRate            loopAccRate                  Jnnz2               MxSqrAvg               MySqrAvg               MzSqrAvg       MalongFieldavg        MalongFieldSqravg      S\n"
    ;
    file <<std::setprecision(16)<< (*data) << "\n";
    file <<"endofSimulationResult\n";
    std::cout<<"File Saved"<<"\n";
}












