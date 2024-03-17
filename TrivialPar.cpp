//
//  MCsimulation.cpp
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-08-30.
//
//

#include "TrivialPar.h"










TrivialPar::TrivialPar(Parameters *parameters, Crystal *crystal, Ewald *ewald , Eigen::VectorXd *spinOrientation)
:
parameters(parameters),
ewald(ewald),
spinOrientation(spinOrientation)
{
    NNlist=(*crystal).getNNlist();
    lattice=(*crystal).getLattice();
    spinDirection=(*crystal).getSpinDirection();
    
    
    diamondLattice=(*crystal).getUpDiamondLattice();
    datapoints=(*parameters).getDataPoints();

    data= Eigen::MatrixXd((*parameters).getDataPoints()*(*parameters).getBins(),18);
    Jnnxy=(*parameters).getJnnxy();
    //Jnnz=(*parameters).getJnnzStart();
    Jnnz1=(*parameters).getJnnz1();
    Jnnz2=(*parameters).getJnnz2();
    Jnnn=(*parameters).getJnnn();
    D=(*parameters).getD();
    B=(*parameters).getBstart();
    T=(*parameters).getTstart();
    Nspins = (*crystal).getNspins();
    Elist = Eigen::VectorXd::Zero(Nspins);
    
    
    fieldDirection=(*parameters).getFieldDirection();
    muNorm = (*parameters).getMuNorm();
    NNlist=(*crystal).getNNlist();
    XYNNlist=(*crystal).getXYNNlist();
    ZNNlist=(*crystal).getZNNlist();
    Z1NNlist=(*crystal).getZ1NNlist();
    Z2NNlist=(*crystal).getZ2NNlist();
    NNNlist=(*crystal).getNNNlist();
    
    
   
    spinTetraList=(*crystal).getSpinTetraList();
    tetraSpinList=(*crystal).getTetraSpinList();
    
    //ewaldMatrix=ewald.getEwaldMatrix();
    
    rng_s.read_seed_local((*parameters).getPathToDataFolder()+"randomSeed/"+"inputSeed_MC.txt");
    
    
    
    sumCorrSiSj = Eigen::MatrixXd::Zero(Nspins, Nspins);
    

}




void TrivialPar::findElist(){
    //this void creates Elist given the current configuration so that it contains the energy for each particle. The full energy is the sum of the half-field energy plus Elist divided by 2(double counting)
    Elist = Eigen::VectorXd::Zero(Nspins);
    Eigen::RowVector3d a;
    Eigen::Vector3d b;
    for(int i=0; i<Nspins; ++i)
    {
        //NNXY------------------------------------------------------------------------------
        for(int j=0; XYNNlist(i,j)>-1; j++)
        {
            a=spinDirection.row(i)*(*spinOrientation)(i);
            b=spinDirection.row(XYNNlist(i,j)).transpose()*(*spinOrientation)(XYNNlist(i,j));
            Elist(i) += Jnnxy*a*b;
        }
        //NNZ1--------------------------------------------------------------------------------
        for(int j=0; Z1NNlist(i,j)>-1 ; j++)
        {
            a=spinDirection.row(i)*(*spinOrientation)(i);
            b=spinDirection.row(Z1NNlist(i,j)).transpose()*(*spinOrientation)(Z1NNlist(i,j));
            Elist(i) += Jnnz1*a*b;
        }
        //NNZ2--------------------------------------------------------------------------------
        for(int j=0; Z2NNlist(i,j)>-1 ; j++)
        {
            a=spinDirection.row(i)*(*spinOrientation)(i);
            b=spinDirection.row(Z2NNlist(i,j)).transpose()*(*spinOrientation)(Z2NNlist(i,j));
            Elist(i) += Jnnz2*a*b;
        }
        //Dipole------------------------------------------------------------------------------
        for(int j=0; j<lattice.rows(); ++j)
        {
            Elist(i)+= (*ewald).ewaldMatrix(i,j)*(*spinOrientation)(i)*(*spinOrientation)(j);
        }
        //Field-------------------------------------------------------------------------------
        Elist(i) += -(*spinOrientation)(i)*muNorm*B*(spinDirection.row(i)).dot(fieldDirection);
    }
}


















Eigen::MatrixXd TrivialPar::run()
{//This void essentially calls the MCsimRun simulation. It first calls MCSimRun with a small nubmer of steps in order to make an estimate of the time needed to complete a simulation.(dont use)
    
    bool timingUsed = ((*parameters).getMCsteps()<0);
    
    if(!timingUsed) MCsimRun(((*parameters).getMCsteps()/(*parameters).getDataPoints())/(*parameters).getBins());
    if(timingUsed)
    {
        //Timing
        double timingSteps= (*parameters).getTimingSteps();
        long long microseconds=0;
        auto start=std::chrono::high_resolution_clock::now();
        std::cout<<"Timing the simulation:\n";
        MCsimRun(timingSteps);
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        stepsPerBin=(unsigned long long int) timingSteps*(double)((pow(10,6)*(*parameters).getRunTime()/(double)microseconds));
        std::cout<<"Simulating: "<<"\n";
        MCsimRun(stepsPerBin);
    }
    return data;
}




Eigen::MatrixXd TrivialPar::run(int myrank, int nproc)
{//Same as ::run() but with myrank input for use in MPI applications.
    
    
    if (myrank == 0){
        //seed with different random numbers for each rank.
        rng_s.read_seed_local((*parameters).getPathToDataFolder()+"randomSeed/"+"inputSeed_MC.txt");
        for (int i = 1; i < nproc; ++i){
            std::vector<unsigned int> seed_to_send = rng_s.get_seed_vector();
            MPI_Send(seed_to_send.data(), (int)seed_to_send.size(), MPI_UNSIGNED, i, 10, MPI_COMM_WORLD);
        }
    }else if(myrank != 0){
        //seed with different random numbers for each rank.
        std::vector<unsigned int> seed_to_receive(4,0);
        MPI_Status status;
        MPI_Recv(seed_to_receive.data(), (int)seed_to_receive.size(), MPI_UNSIGNED, 0, 10, MPI_COMM_WORLD, &status);
        rng_s.set_seed(seed_to_receive[0],seed_to_receive[1],seed_to_receive[2],seed_to_receive[3]);
    }
    MCsimRun(((*parameters).getMCsteps()/(*parameters).getDataPoints())/(*parameters).getBins());
   
    return data;
}

















Eigen::RowVector3d TrivialPar::Magnetization()
{
    Eigen::RowVector3d Mag(0,0,0);
    for(int n=0; n<lattice.rows(); n++)
    {
        Mag += spinDirection.row(n)*(*spinOrientation)(n);
    }
    return Mag;
}







double TrivialPar::Energy()
{//This function calculates the total energy of the system, which is not simply the sum of Elist, due to the fact that field is not double counted.
    double E=0;
    findElist();
    E+=Elist.sum()/2;
    //adding an extra half of the field energy since this should not be divided by 2 in the summation.
    for(int i=0; i<Nspins; ++i)
    {
        E += -0.5*(*spinOrientation)(i)*muNorm*B*(spinDirection.row(i)).dot(fieldDirection);
    }
    return E;
}











void TrivialPar::MCsimRun(unsigned long long int stepsPerBin)
{//This void carries out the metropolis algorithm on the lattice. it makes a sweep from "Tstart" to "Tend" and saves "TdataPoints" in a matrix "data"
    double eqPercentage= (*parameters).getEqPercentage();
    double Tstart=       (*parameters).getTstart();
    double Tend=         (*parameters).getTend();
    double Bstart=       (*parameters).getBstart();
    double Bend=         (*parameters).getBend();
    //double JnnzStart=    (*parameters).getJnnzStart();
    //double JnnzEnd=      (*parameters).getJnnzEnd();
 
    
    
    
    
    
    //Initializing datapoint variables----------------------------
    bins=(*parameters).getBins();
    MCstepsPerDatapoint=stepsPerBin*bins;
    stepsAcceptedInBin=0;
    loopStepsAcceptedInBin=0;
    T=Tstart;
    B=Bstart;
    //Jnnz=JnnzStart;
    
    //Initializing "directly" observable quantites-----------------
    M=Magnetization();
    findElist();
    E=Energy();
    //Initializing the local sums used in each bin-------------------
    double Esum=0;
    double EsqrSum=0;
    Eigen::Vector3d MvecSum(0,0,0);
    Eigen::Vector3d MSqrVecSum(0,0,0);
    double MalongFieldSum=0;
    double MalongFieldSqrSum=0;
    double Msum=0;
    double MsqrSum=0;
    // double Zapprox=0; not working dont use
    
    //There are 3 for-loops, the datapoint loop, bin loop and MCstep loop
    for(int k=0; k < datapoints; ++k)
    {
        std::cout<<"Progress: "<<100*(k/(double)datapoints)<<"%"<<"\n";
        
        
        
        
        //Jnnz=(JnnzEnd-JnnzStart)*k/(double)datapoints+JnnzStart;
       
        
        if((*parameters).getLogSteps())
        {
            B=(Bend-Bstart)*log(k+1)/(double)log(datapoints+1)+Bstart;
            T=(Tend-Tstart)*log(k+1)/(double)log(datapoints+1)+Tstart;
        }else
        {
            B=(Bend-Bstart)*k/(double)datapoints+Bstart;
            T=(Tend-Tstart)*k/(double)datapoints+Tstart;
        }
       
        
        findElist();
        //E=Elist.sum()/2;
        E=Energy();
        
        for(int step=0; step<bins*stepsPerBin*(eqPercentage/100); step++)
        {//Equilibrate the system small for loop.
            E+=singleSweepList();
         //   E+=loopSweepList();
        }
        for(int bin=0; bin<bins; bin++)
        {
            //std::cout<<"E="<<E<<"\n";
            MvecSum<<0,0,0;
            MSqrVecSum<<0,0,0;
            MalongFieldSum=0;
            MalongFieldSqrSum=0;
            Msum=0;//resetting all local sums for each bin
            Esum=0;
            MsqrSum=0;
            EsqrSum=0;
            stepsAcceptedInBin=0;
            loopStepsAcceptedInBin=0;
            // Zapprox=0; not working dont use
            for(unsigned int step=0; step<stepsPerBin; step++)
            {
                E+=singleSweepList(); //spin flip
            //    E+=loopSweepList();   //loop flip
                Esum+=E;//summing of local quantities in order to create an average
                EsqrSum+=E*E;
                MvecSum+=M;
                MSqrVecSum.array()+= M.array()*M.array();
                MalongFieldSum+=M.dot(fieldDirection);
                MalongFieldSqrSum+=pow(M.dot(fieldDirection),2);
                Msum+=M.norm();
                MsqrSum+=pow(M.norm(),2);
                AddToCorrelationSiSj();
                //Zapprox=Zapprox+exp(-E/T); not working dont use
            }
            //Collecting all averages over the bin
            data(k*bins+bin,0)=T;
            data(k*bins+bin,1)=B;
            data(k*bins+bin,2)=Esum/stepsPerBin;
            data(k*bins+bin,3)=EsqrSum/stepsPerBin;
            data(k*bins+bin,4)=Msum/stepsPerBin;
            data(k*bins+bin,5)=MsqrSum/stepsPerBin;
            data(k*bins+bin,6)=MvecSum(0)/stepsPerBin;
            data(k*bins+bin,7)=MvecSum(1)/stepsPerBin;
            data(k*bins+bin,8)=MvecSum(2)/stepsPerBin;
            data(k*bins+bin,9)=stepsAcceptedInBin/((double)stepsPerBin*Nspins);
            data(k*bins+bin,10)=loopStepsAcceptedInBin/((double)stepsPerBin*Nspins);
            data(k*bins+bin,11)=Jnnz2;
            data(k*bins+bin,12)=MSqrVecSum(0)/stepsPerBin;
            data(k*bins+bin,13)=MSqrVecSum(1)/stepsPerBin;
            data(k*bins+bin,14)=MSqrVecSum(2)/stepsPerBin;
            data(k*bins+bin,15)=MalongFieldSum/stepsPerBin;
            data(k*bins+bin,16)=MalongFieldSqrSum/stepsPerBin;
            //data(k*bins+bin,17)= Nspins*log(2*Zapprox/stepsPerBin)+(Esum/stepsPerBin)/T;//not working dont use
        }
    }
}















double TrivialPar::singleSweepList()
{
    //This function performs Nspins single spin flips on random positions and returns the change in energy of the system if some spins are flipped
    
    //generic energy and vectors for vector products
    double dE=0;
    double de=0;
    double r=1;
    int i=-1;
    int j=0;
    int k=0;
    Eigen::RowVector3d a;
    Eigen::Vector3d b;
    
    
    for(int N=0; N<Nspins; ++N)
    {
        i= rng_s.r_int(Nspins);
        r= rng_s.r_double1();
        de= -2*Elist(i);
        if(r<exp(-de/T)){//Acceptance condition
            dE+=de;
            Elist(i)= -Elist(i);//updating field and interaction energies for i.
   
            //Updating the energy both for its neighbors (j and k)
            //NNXY--------------------------------------------------------------------------------------------------------------
            for(j=0; XYNNlist(i,j)>-1; ++j)
            {
                Elist(XYNNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(XYNNlist(i,j)))*Jnnxy*(*spinOrientation)(i)*(*spinOrientation)(XYNNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
            }
            //NNZ1---------------------------------------------------------------------------------------------------------------
            for(j=0; Z1NNlist(i,j)>-1; ++j)
            {
                Elist(Z1NNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(Z1NNlist(i,j)))*Jnnz1*(*spinOrientation)(i)*(*spinOrientation)(Z1NNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
            }
            //NNZ2---------------------------------------------------------------------------------------------------------------
            for(j=0; Z2NNlist(i,j)>-1; ++j)
            {
                Elist(Z2NNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(Z2NNlist(i,j)))*Jnnz2*(*spinOrientation)(i)*(*spinOrientation)(Z2NNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
            }
            //Dipolar------------------------------------------------------------------------------------------------------------
            for(k=0; k<Nspins; ++k)
            {// add/subtract energy by one column from the Ewald Matrix
                if(i!=k) Elist(k)+= -2*(*ewald).ewaldMatrix(i,k)*(*spinOrientation)(i)*(*spinOrientation)(k);
            }
            M+= -2*spinDirection.row(i)*(*spinOrientation)(i);
            (*spinOrientation)(i)=-(*spinOrientation)(i);
            stepsAcceptedInBin++;
        }
    }
    return dE;
}














double TrivialPar::loopSweepList()
{//This void tries making Nspins loop flips, if a loop fails to be created due to too many excitations in the crystal, the specific attempt is aborted
    
    double de=0;
    double dE=0;
    double r;
    double loopEnergy=0;
    int i=-1;
    int j=-1;
    int k=-1;
    
    
    
    for(int N=0; N<Nspins; ++N)
    {
        de=0;
        if(findLoop())
        {
            loopEnergy=0;
            for(unsigned int MemberChainNumber=0; MemberChainNumber<shortSpinLoop.size(); ++MemberChainNumber)
            {
                i=shortSpinLoop.at(MemberChainNumber);
                loopEnergy+= Elist(i);
                for(unsigned int MemberChainNumber2=0; MemberChainNumber2<shortSpinLoop.size(); ++MemberChainNumber2)
                {
                    if(MemberChainNumber!=MemberChainNumber2)
                    {//Compensating loopEnergy so that we do not count the energy against other spins in the loop
                        j=shortSpinLoop.at(MemberChainNumber2);
                        if((XYNNlist.row(i).array() == j).any()) loopEnergy+=  -(spinDirection.row(i)).dot(spinDirection.row(j))*Jnnxy*(*spinOrientation)(i)*(*spinOrientation)(j);
                        if(( Z1NNlist.row(i).array() == j).any()) loopEnergy+=  -(spinDirection.row(i)).dot(spinDirection.row(j))*Jnnz1*(*spinOrientation)(i)*(*spinOrientation)(j);
                        if(( Z2NNlist.row(i).array() == j).any()) loopEnergy+=  -(spinDirection.row(i)).dot(spinDirection.row(j))*Jnnz2*(*spinOrientation)(i)*(*spinOrientation)(j);
                        loopEnergy+=  -(*ewald).ewaldMatrix(i,j)*(*spinOrientation)(i)*(*spinOrientation)(j);
                    }
                }
            }
            de=-2*loopEnergy;
            r=rng_s.r_double();
            if(r<exp(-de/T)){
                dE+=de;
                loopStepsAcceptedInBin++;
                for(unsigned int MemberChainNumber=0; MemberChainNumber<shortSpinLoop.size(); ++MemberChainNumber)
                {
                    i=shortSpinLoop.at(MemberChainNumber);
                    Elist(i)= -Elist(i); //Updating for i
                    //Updating for its neighbors (j and k)
                    //NNXY--------------------------------------------------------------------------------------------------------------
                    for(j=0; XYNNlist(i,j)>-1; ++j)
                    {
                        Elist(XYNNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(XYNNlist(i,j)))*Jnnxy*(*spinOrientation)(i)*(*spinOrientation)(XYNNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
                    }
                    //NNZ1---------------------------------------------------------------------------------------------------------------
                    for(j=0; Z1NNlist(i,j)>-1; ++j)
                    {
                        Elist(Z1NNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(Z1NNlist(i,j)))*Jnnz1*(*spinOrientation)(i)*(*spinOrientation)(Z1NNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
                    }
                    //NNZ2---------------------------------------------------------------------------------------------------------------
                    for(j=0; Z2NNlist(i,j)>-1; ++j)
                    {
                        Elist(Z2NNlist(i,j))+=  -2*(spinDirection.row(i)).dot(spinDirection.row(Z2NNlist(i,j)))*Jnnz2*(*spinOrientation)(i)*(*spinOrientation)(Z2NNlist(i,j)); //we flip the sign of this energy contribution, and thus the factor 2
                    }
                    //Dipolar------------------------------------------------------------------------------------------------------------
                    for(k=0; k<Nspins; ++k)
                    {// add/subtract energy by one column from the Ewald Matrix
                        if(i!=k) Elist(k)+= -2*(*ewald).ewaldMatrix(i,k)*(*spinOrientation)(i)*(*spinOrientation)(k);
                    }
                    M+= -2*spinDirection.row(i)*(*spinOrientation)(i);
                    (*spinOrientation)(i)=-(*spinOrientation)(i);
                }
            }
        }
    }
    return dE;
}









bool TrivialPar::findLoop()
{
    std::vector<int> spinLoop;
    std::vector<int> tetraLoop;
    
    
    //random starting point
    int spin=rng_s.r_int(Nspins);
    int upDownTetra=rng_s.r_2(); //up or down lattice when walking away from the first spin
    int tetrahedron= spinTetraList(spin,upDownTetra);
   
    bool loopFound=false;
    spinLoop.push_back(spin);
    tetraLoop.push_back(tetrahedron);
    int intersectionNumber=0;
    int interiorChoicesOfnextSpin[4] = {-1,-1,-1,-1};
    int interiorSpinCounter=0;
    
    do
    {
        interiorSpinCounter=0;
        for(int s=0; s<4; ++s)
        {
            if((*spinOrientation)(spinLoop.back())!=(*spinOrientation)(tetraSpinList(tetrahedron,s)))
            {
                interiorChoicesOfnextSpin[interiorSpinCounter++]=s;
            }
        }
        if(interiorSpinCounter!=2)
        {
            //std::cout<<"Entered bad tetrahedron, abort and return false\n";
            return loopFound;
        }
        
        
        
        spin=tetraSpinList(tetrahedron,interiorChoicesOfnextSpin[rng_s.r_2()]);
        spinLoop.push_back(spin);
        upDownTetra=(upDownTetra+1)%2;
        tetrahedron=spinTetraList(spin,upDownTetra);
        tetraLoop.push_back(tetrahedron);
        
        
        for(unsigned int i=0; i<tetraLoop.size()-1; ++i)
        {
            if(tetraLoop.at(i)==tetrahedron)
            {
                loopFound=true;
                intersectionNumber=i;
            }
        }
        
        
    }while(!loopFound);
    shortSpinLoop = std::vector<int>(spinLoop.begin()+intersectionNumber+1, spinLoop.end());
    return loopFound;
}















// This void adds the current value of SiSj to the correlation matrix
/*************************************************************************************/
inline void TrivialPar::AddToCorrelationSiSj()
{
    sumCorrSiSj.matrix().noalias() += (*spinOrientation).matrix()*(*spinOrientation).matrix().transpose();
}

































void TrivialPar::saveToFile(std::string path)
{
    
    std::cerr<<"WARNING TRIVIALAR::SAVETOFILE CALLED, THIS VOID IS OLD AND NOT USED. PROGRAM WILL EXIT.";
    exit(EXIT_FAILURE);
    
    /*
    //This void saves the matrix data to the file data.txt
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::string timeString=std::to_string(now->tm_year-100)+std::to_string(now->tm_mon+1)+std::to_string(now->tm_mday)+std::to_string(now->tm_hour)+std::to_string(now->tm_min);
    
    
    
    std::string filePathAndName=path+(*parameters).getOutputFilename()+timeString+".txt";
    if(!(*parameters).getAddDateStamp()) filePathAndName=path+(*parameters).getOutputFilename()+".txt";
    std::ofstream file(filePathAndName);
  
    //exit if unable to create file
    if(!file){
        std::cerr<<"File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file<<"cellX= "+std::to_string((*parameters).getCellX())<<"\n";
    file<<"cellY= "+std::to_string((*parameters).getCellY())<<"\n";
    file<<"cellZ= "+std::to_string((*parameters).getCellZ())<<"\n";
    file<<"Nspins= "+ std::to_string(lattice.rows())<<"\n";
    file<<"MCsteps= "+ std::to_string(MCstepsPerDatapoint*datapoints)<<"\n";
    //file<<"elapsed time: "<<elapsed_secs/3600<<"h:"<<(elapsed_secs%3600)/60<<"m:"<<elapsed_secs%60<<"s\n";
    file<<"eqPercentage= " +std::to_string((int) (*parameters).getEqPercentage())<<"%\n";
    file<<"Bins= "+ std::to_string((*parameters).getBins())<<"\n";
    file<<"MaxEwaldMirrors= "<<(*parameters).getMaxEwaldMirrors()<<"\n";
    file<<"deformation matrix F:\n";
    file<<(*parameters).getFmatrix()<<"\n";
    file<<"dataPoints= "+ std::to_string((*parameters).getDataPoints())<<"\n";
    file<<"Tstart= "+ std::to_string((*parameters).getTstart())<<"\n";
    file<<"Tend= "+ std::to_string((*parameters).getTend())<<"\n";
    file<<"File= \""+(*parameters).getOutputFilename() +"\"\n";
    file<<"Model: NN + Dipole\n";
    file<<"Jnnxy= "+ std::to_string(Jnnxy)<<" K\n";
    file<<"JnnzStart= "+ std::to_string(JnnzStart)<<" K\n";
    file<<"JnnzEnd= "+ std::to_string(JnnzEnd)<<" K\n";
    file<<"Jnnn= "+ std::to_string(Jnnn)<<" K \n";
    file<<"D= "+ std::to_string(D)<<" K \n";
    file<<"Bstart= "+ std::to_string((*parameters).getBstart())<<"\n";
    file<<"Bend= "+std::to_string((*parameters).getBend())<<"\n";
    file<<"unitCellSideLength= " + std::to_string((*parameters).getUnitCellSideLength())<<" Å\n";
    file<<"SimulationResult:\n";
    file<<"                     T             B              Eavg                EsqrAvg                   Mavg                MsqrAvg                     Mx                     My                     Mz          singleAccRate            loopAccRate\n"
;
    file <<std::setprecision(16)<< data;
    std::cout<<"File Saved"<<"\n";
     */
}






void TrivialPar::saveStateToFile(std::string path, int i)
{
    //This void the current spin state to a file state_i
    
    
    std::string filePathAndName=path+(*parameters).getOutputFilename()+"_state_"+std::to_string(i)+".txt";
    std::ofstream file(filePathAndName);
    
    //exit if unable to create file
    if(!file){
        std::cerr<<"File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file <<std::setprecision(16)<< (*spinOrientation);
    std::cout<<"File Saved"<<"\n";
}













void TrivialPar::test()
{
    
    while(!findLoop())
    {
        std::cout<<"finding a loop";
    }
    std::cout<<"found a loop";
    
    while(true)
    {
        for(unsigned int i=0; i<shortSpinLoop.size(); ++i)
        {
            (*spinOrientation)(shortSpinLoop.at(i)) = -(*spinOrientation)(shortSpinLoop.at(i));
        }
    }
     
}





//This void returns the sumCorrSiSj matrix
/*****************************************************************************************************/
Eigen::MatrixXd TrivialPar::getSumCorrSiSj()
{
    return sumCorrSiSj;
}


//This void returns the data matrix
/*****************************************************************************************************/
Eigen::MatrixXd TrivialPar::getData()
{
    return data;
}
















