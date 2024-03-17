//
//  Parameters.cpp
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-10-05.
//


#include "Parameters.h"

Parameters::Parameters(std::string path, std::string parametersFilename)
:path(path),
parametersFilename(parametersFilename)
{
    writeCorrelationMatrix=0;
    MCsteps=1;
    fieldTheta=std::numeric_limits<double>::quiet_NaN();
    fieldPhi=std::numeric_limits<double>::quiet_NaN();
    oneMkappa=std::numeric_limits<double>::quiet_NaN();
    
    demagMatrix=std::numeric_limits<double>::quiet_NaN()*Eigen::Matrix3d::Zero(3, 3);
    
    
    
    addEwaldDemagTerm= false;
    readConfiguration=false;
    load();
   
}


void Parameters::load()
{//This void loads the simulation parameters from a file in path
    
    std::ifstream inputFile;
    std::string s_reader, t_reader;
    //std::cout<<path+parametersFilename;
    inputFile.open((path+parametersFilename).c_str());
    if (!inputFile.is_open()){
        std::cerr<<("Parameters file not found, exiting\n");
        
        exit(EXIT_FAILURE);
    }
    
    
    bool FmatrixFound=false;
    bool DemagMatrixFound=false;
    while (getline(inputFile,s_reader))
    {
        std::istringstream r_pick(s_reader);
        r_pick >> t_reader;
        if(t_reader=="outputFilename=") r_pick >> outputFilename;
        if(t_reader=="writeCorrelationMatrix=") r_pick >> writeCorrelationMatrix;
        if(t_reader=="cellType=") r_pick >> cellType;
        if(t_reader=="cellSide=") r_pick >> cellSide;
        if(t_reader=="cellX=") r_pick >> cellX;
        if(t_reader=="cellY=") r_pick >> cellY;
        if(t_reader=="cellZ=") r_pick >> cellZ;
        if(t_reader=="fieldPhi=") r_pick >> fieldPhi;
        if(t_reader=="fieldTheta=") r_pick >> fieldTheta;
        if(t_reader=="dataPoints=") r_pick >> dataPoints;
        if(t_reader=="Tstart=") r_pick >> Tstart;
        if(t_reader=="Tend=") r_pick >> Tend;
        if(t_reader=="MCsteps=") r_pick >> MCsteps;
        if(t_reader=="FinishTime=") r_pick>> endTimeStr;
        if(t_reader=="runTime=") r_pick>> runTime;
        if(t_reader=="bins=") r_pick >> bins;
        if(t_reader=="Jnnxy=") r_pick >> Jnnxy;
        if(t_reader=="Jnnz1=") r_pick >> Jnnz1;
        if(t_reader=="Jnnz2=") r_pick >> Jnnz2;
        if(t_reader=="Jnnn=") r_pick >> Jnnn;
        if(t_reader=="D=") r_pick >> D;
        if(t_reader=="kb=") r_pick >> kb;
        if(t_reader=="Bstart=") r_pick >> Bstart;
        if(t_reader=="Bend=") r_pick >> Bend;
        if(t_reader=="unitCellSideLength=") r_pick >> unitCellSideLength;
        if(t_reader=="eqPercentage=") r_pick >> eqPercentage;
        if(t_reader=="addDateStamp=") r_pick >> addDateStamp;
        if(t_reader=="loadCrystal=") r_pick >> loadCrystal;
        if(t_reader=="loadEwald=") r_pick >> loadEwald;
        if(t_reader=="fieldIndex=") r_pick >> fieldIndex;
        if(t_reader=="pIndex=") r_pick >> pIndex;
        if(t_reader=="timingSteps=") r_pick >> timingSteps;
        if(t_reader=="mu=") r_pick >> mu;
        if(t_reader=="addEwaldDemagTerm=") r_pick >> addEwaldDemagTerm;
        if(t_reader=="replicaExchanges=") r_pick >> replicaExchanges;
        if(t_reader=="logSteps=") r_pick >> logSteps;
        if(t_reader=="mubDivkb=") r_pick >> mubDivkb;
        if(t_reader=="MaxEwaldMirrors=") r_pick >> MaxEwaldMirrors;
        if(t_reader=="readConfiguration=") r_pick >> readConfiguration;
        if(t_reader=="oneMkappa=") r_pick >> oneMkappa;
        if(t_reader=="DeformationMatrixF:" && !FmatrixFound){
            for(int i=0; i<3; ++i)
            {
                getline(inputFile,s_reader);
                std::istringstream r_pick(s_reader);
                for(int j=0; j<3; ++j)
                {
                    r_pick >> F(i,j);
                }
            }
            FmatrixFound=true;
        }
        if(t_reader=="DemagMatrix:" && !DemagMatrixFound){
            for(int i=0; i<3; ++i)
            {
                getline(inputFile,s_reader);
                std::istringstream r_pick(s_reader);
                for(int j=0; j<3; ++j)
                {
                    r_pick >> demagMatrix(i,j);
                }
            }
            DemagMatrixFound=true;
        }
    }
    inputFile.close();
    CheckInputandGetParameters();//check so that input was correct and set quantities;
    std::cout<<"Parameters loaded\n";
}








void Parameters::CheckInputandGetParameters()
{   //This void checks so that input was correct and does some arithmetic in order to find parameters from parameters, such as finding RunTime and NNdist
    findCellDimensions();
    findFieldDirection();
    findPDirection();
    findRuntime();
    checkInput();
    createFmatrixFromOneMkappa();
    
}


void Parameters::createFmatrixFromOneMkappa()
{
    
    
    if(oneMkappa==oneMkappa)//check if nan
    {
        
        Eigen::Matrix3d D;
        Eigen::Matrix3d P;
        Eigen::Matrix3d invP;
        
        const std::string s001="001";
        const std::string s110="110";
        const std::string s111="111";
        
        
        if(pIndex.compare(s001)==0)
        {
            D<< 1,0,0,
            0,1,0,
            0,0,oneMkappa;
            
            F=D;
        }else if(pIndex.compare(s110)==0)
        {
            P<<1/sqrt(2),1/sqrt(2),0,
              -1/sqrt(2),1/sqrt(2),0,
               0,0,1;
            invP<<1/sqrt(2),-1/sqrt(2),0
                 ,1/sqrt(2),1/sqrt(2), 0,
                  0,0,1;
            D<< 1,0,0,
                0,oneMkappa,0,
                0,0,1;
            F= P*D*invP;
        }else if(pIndex.compare(s111)==0)
        {
        
            P<<1/sqrt(2), 1/sqrt(6), 1/sqrt(3),
              -1/sqrt(2),  1/sqrt(6), 1/sqrt(3),
               0,        -2/sqrt(6), 1/sqrt(3);
            invP<<1/sqrt(2),-1/sqrt(2), 0,
                  1/sqrt(6), 1/sqrt(6), -2/sqrt(6),
                  1/sqrt(3), 1/sqrt(3), 1/sqrt(3);
            D<< 1,0,0,
            0,1,0,
            0,0,oneMkappa;
            F= P*D*invP;
        }else
        {
            std::cerr<<"Parameters.cpp:: error with pIndex input.\n";
            exit(EXIT_FAILURE);
        }
    }else
    {
        std::cout<<"oneKappa input not detected, preceding with manual F matrix input.\n";
    }
    
    
}





void Parameters::checkInput()
{
    if(eqPercentage<0)
    {
        std::cerr<<"eqPercentage incorrect format\n";
        exit(EXIT_FAILURE);
    }
    if(runTime<1 && MCsteps <0)
    {
        std::cerr<<"Error: Invalid runTime or MCsteps combination, check finish time"<<"\n";
        exit(EXIT_FAILURE);
    }
}







void Parameters::findCellDimensions()
{
    if(cellSide>0)
    {
        cellX=cellSide;
        cellY=cellSide;
        cellZ=cellSide;
    }
}








void Parameters::findPDirection()
{
    double a=std::stoi(pIndex.substr(0,1));
    double b=std::stoi(pIndex.substr(1,1));
    double c=std::stoi(pIndex.substr(2,1));
    pDirection<<a,b,c;
    pDirection=pDirection/pDirection.norm();
}




void Parameters::findFieldDirection()
{
    
    double b1=std::stoi(fieldIndex.substr(1,1));
    double b2=std::stoi(fieldIndex.substr(2,1));
    double b3=std::stoi(fieldIndex.substr(3,1));
    fieldDirection<<b1,b2,b3;
    fieldDirection=fieldDirection/fieldDirection.norm();
    
    //if missalignment, alhpa>0
    if(fieldPhi==fieldPhi && fieldTheta==fieldTheta)
    {
        b1= std::sin(M_PI/180*fieldTheta)*std::cos(M_PI/180*fieldPhi);
        b2= std::sin(M_PI/180*fieldTheta)*std::sin(M_PI/180*fieldPhi);
        b3= std::cos(M_PI/180*fieldTheta);
        fieldDirection<<b1,b2,b3;
        fieldDirection=fieldDirection/fieldDirection.norm();
    }else
    {
        fieldPhi=std::numeric_limits<double>::quiet_NaN();
        fieldTheta=std::numeric_limits<double>::quiet_NaN();
    }
}





void Parameters::findRuntime()
{
    if(MCsteps==-1 && runTime==-1)
    {
        time_t t = time(0);   // get time now
        struct tm * startingTime = localtime( & t );
        
        std::string timeString=std::to_string(startingTime->tm_year+1900)+"-"+std::to_string(startingTime->tm_mon+1)+"-"+std::to_string(startingTime->tm_mday)+"-"+std::to_string(startingTime->tm_hour)+":"+std::to_string(startingTime->tm_min)+":"+std::to_string(startingTime->tm_sec);
        
        int startYear= startingTime->tm_year+1900;
        int startMonth= startingTime->tm_mon+1;
        int startDay= startingTime->tm_mday;
        int startHour= startingTime->tm_hour;
        int startMinute= startingTime->tm_min;
        int startSecond= startingTime->tm_sec;
        
        
        
        int endYear= std::stoi(endTimeStr.substr(0,4));
        int endMonth= std::stoi(endTimeStr.substr(5,2));
        int endDay= std::stoi(endTimeStr.substr(8,2));
        int endHour= std::stoi(endTimeStr.substr(11,2));
        int endMinute= std::stoi(endTimeStr.substr(14,2));
        int endSecond= std::stoi(endTimeStr.substr(17,2));
        
        /*
         std::cout<<"endYear="<<endYear<<"\n";
         std::cout<<"endMonth="<<endMonth<<"\n";
         std::cout<<"endDay="<<endDay<<"\n";
         std::cout<<"endHour="<<endHour<<"\n";
         std::cout<<"endMinute="<<endMinute<<"\n";
         std::cout<<"endSecond="<<endSecond<<"\n";
         */
        
        //std::cout<<"diffyears: "<<endYear-startYear<<"\n";
        //std::cout<<"diffmonth: "<<endMonth-startMonth<<"\n";
        //std::cout<<"diffday: "<<endDay-startDay<<"\n";
        //std::cout<<"diffhour: "<<endHour-startHour<<"\n";
        //std::cout<<"diffmin: "<<endMinute-startMinute<<"\n";
        //std::cout<<"diffsec: "<<endSecond-startSecond<<"\n";
        
        
        runTime=(endYear-startYear)*12*30*24*60*60+(endMonth-startMonth)*30*24*60*60+(endDay-startDay)*24*60*60+(endHour-startHour)*60*60+(endMinute-startMinute)*60+(endSecond-startSecond);
    }
}





Eigen::Matrix3d Parameters::getDemagMatrix() const
{
    return demagMatrix;
}




bool Parameters::getAddEwaldDemagTerm() const
{
    return addEwaldDemagTerm;
}





int Parameters::getReplicaExchanges() const
{
    return replicaExchanges;
}



Eigen::Vector3d Parameters::getPDirection() const
{
    return pDirection;
}

std::string Parameters::getPIndex() const
{
    return pIndex;
}

bool Parameters::getLogSteps() const
{
    return logSteps;
}



double Parameters::getUnitCellSideLength() const
{
    return unitCellSideLength;
}
int Parameters::getMaxEwaldMirrors() const
{
    return MaxEwaldMirrors;
}
double Parameters::getMu() const
{
    return mu;
}
Eigen::Matrix3d Parameters::getFmatrix() const
{
    return F;
}
std::string Parameters::getPathToDataFolder() const
{
    return path+"data/";
}
int Parameters::getRunTime() const
{
    return runTime;
}
int Parameters::getTimingSteps() const
{
    return timingSteps;
}
bool Parameters::getLoadCrystal() const
{
    return loadCrystal;
}
bool Parameters::getLoadEwald() const
{
    return loadEwald;
}
bool Parameters::getAddDateStamp() const
{
    return addDateStamp;
}
double Parameters::getEqPercentage() const
{
    return eqPercentage;
}
std::string Parameters::getOutputFilename() const
{
    return outputFilename;
}
int Parameters::getCellType() const
{
    if(cellType==16 || cellType==24)
    {
        return cellType;
    }
    else
    {
        std::cerr<<"Input error cellType\n";
        return 0;
    }
}
int Parameters::getCellX() const
{
    return cellX;
}
int Parameters::getCellY() const
{
    return cellY;
}
int Parameters::getCellZ() const
{
    return cellZ;
}
double Parameters::getTstart() const
{
    return Tstart;
}
double Parameters::getTend() const
{
    return Tend;
}
double Parameters::getJnnxy() const
{
    return Jnnxy;
}
double Parameters::getJnnz1() const
{
    return Jnnz1;
}
double Parameters::getJnnz2() const
{
    return Jnnz2;
}


double Parameters::getJnnn() const
{
    return Jnnn;
}
double Parameters::getD() const
{
    return D;
}
double Parameters::getMuNorm() const
{
    return mu*mubDivkb;
}
double Parameters::getBstart() const
{
    return Bstart;
}

double Parameters::getBend() const
{
    return Bend;
}
Eigen::Vector3d Parameters::getFieldDirection() const
{
    return fieldDirection;
}
int Parameters::getDataPoints() const
{
    return dataPoints;
}


signed long long int Parameters::getMCsteps() const
{
    return MCsteps;
}




bool Parameters::getReadConfiguration() const
{
    return readConfiguration;
}



int Parameters::getBins() const
{
    return bins;
}

bool Parameters::getWriteCorrelationMatrix() const
{
    return writeCorrelationMatrix;
}




std::string Parameters::print() const
{   //This void prints all parameters
    
    std::string output="";
    
    
    
    output+="File Parameters-----------------------------\n";
    output+="\tloadCrystal= "+std::to_string(getLoadCrystal())+"\n";
    output+="\tloadEwald= "+std::to_string(getLoadEwald())+"\n";
    output+="\tpath= "+path+"\n";
    output+="\toutputFileName= "+outputFilename+"\n";
    output+="writeCorrelationMatrix= "+std::to_string(getWriteCorrelationMatrix())+"\n";
    output+="Cell Parameters-----------------------------\n";
    output+="\tcellType= "+std::to_string(getCellType())+"\n";
    output+="\tcellX= "+std::to_string(getCellX())+"\n";
    output+="\tcellY= "+std::to_string(getCellY())+"\n";
    output+="\tcellZ= "+std::to_string(getCellZ())+"\n";
    output+="\tNspins= "+std::to_string(cellType*cellX*cellY*cellZ)+"\n";
    output+="\tDeformationMatrixF:\n";
    output+="\t\t"+std::to_string(getFmatrix()(0,0))+"\t"+std::to_string(getFmatrix()(0,1))+"\t"+std::to_string(getFmatrix()(0,2))+"\n\t\t"+std::to_string(getFmatrix()(1,0))+"\t"+std::to_string(getFmatrix()(1,1))+"\t"+std::to_string(getFmatrix()(1,2))+"\n\t\t"+std::to_string(getFmatrix()(2,0))+"\t"+std::to_string(getFmatrix()(2,1))+"\t"+std::to_string(getFmatrix()(2,2))+"\n";
    output+="\toneMkappa="+std::to_string(oneMkappa)+"\n";
    output+="MC Parameters-------------------------------\n";
    output+="\tDuration:\n";
    output+="\t\tMCstepsPerCore= "+std::to_string(getMCsteps())+"\n";
    output+="\tStatistics:\n";
    output+="\t\treadConfiguration="+std::to_string(getReadConfiguration())+"\n";
    output+="\t\tdataPoints= "+std::to_string(getDataPoints())+"\n";
    output+="\t\treplicaExchanges= "+std::to_string(getReplicaExchanges())+"\n";
    output+="\t\teqPercentage= "+std::to_string(getEqPercentage())+"\n";
    output+="\t\tbinsPerCore= "+std::to_string(getBins())+"\n";
    output+="\t\tlogSteps= "+std::to_string(getLogSteps())+"\n";
    output+="\tEwald:\n";
    output+="\t\tMaxEwaldMirrors= "+std::to_string(getMaxEwaldMirrors())+"\n";
    output+="Physical Parameters-------------------------\n";
    output+="\taddEwaldDemagTerm= "+std::to_string(getAddEwaldDemagTerm())+"\n";
    output+="\tDemagMatrix:\n";
    output+="\t\t"+std::to_string(getDemagMatrix()(0,0))+"\t"+std::to_string(getDemagMatrix()(0,1))+"\t"+std::to_string(getDemagMatrix()(0,2))+"\n\t\t"+std::to_string(getDemagMatrix()(1,0))+"\t"+std::to_string(getDemagMatrix()(1,1))+"\t"+std::to_string(getDemagMatrix()(1,2))+"\n\t\t"+std::to_string(getDemagMatrix()(2,0))+"\t"+std::to_string(getDemagMatrix()(2,1))+"\t"+std::to_string(getDemagMatrix()(2,2))+"\n";
    
    
    output+="\tTstart= "+std::to_string(getTstart())+"\n";
    output+="\tTend= "+std::to_string(getTend())+"\n";
    output+="\tBstart= "+std::to_string(getBstart())+"\n";
    output+="\tBend= "+std::to_string(getBend())+"\n";
    output+="\tJnnxy= "+std::to_string(getJnnxy())+"\n";
    output+="\tJnnz1= "+std::to_string(getJnnz1())+"\n";
    output+="\tJnnz2= "+std::to_string(getJnnz2())+"\n";
    output+="\tD= "+std::to_string(getD())+"\n";
    output+="\tJnnn= "+std::to_string(getJnnn())+"\n";
    output+="\tmu= "+std::to_string(mu)+"\n";
    output+="\tmubDivkb= "+std::to_string(mubDivkb)+"\n";
    output+="\tpDirection:\t"+std::to_string(getPDirection()(0))+"\t"+std::to_string(getPDirection()(1))+"\t"+std::to_string(getPDirection()(2))+"\n";
    output+="\tfieldDirection:\t"+std::to_string(getFieldDirection()(0))+"\t"+std::to_string(getFieldDirection()(1))+"\t"+std::to_string(getFieldDirection()(2))+"\n";
    output+="\tfieldTheta= "+std::to_string(fieldTheta)+"\n";
    output+="\tfieldPhi= "+std::to_string(fieldPhi)+"\n";
    
    
    
    
    return output;
    
    
    
}



