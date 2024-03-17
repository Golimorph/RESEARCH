//
//  MCsimulation.cpp
//  PyrochloreG
//
//  Created by Richard Edberg on 2017-08-30.
//
//

#include "ParallelTempering.h"



ParallelTempering::ParallelTempering(Parameters parameters, Crystal crystal, Ewald ewald , Eigen::VectorXd *spinOrientation, int myrankArg, int nprocArg)
:
parameters(parameters),
crystal(crystal),
ewald(ewald),
spinOrientation(spinOrientation)
{

    NNlist=crystal.getNNlist();
    lattice=crystal.getLattice();
    spinDirection=crystal.getSpinDirection();
    diamondLattice=crystal.getUpDiamondLattice();
    datapoints=parameters.getDataPoints();
    
    myrank = myrankArg;
    nproc = nprocArg;
    
    
    if(nproc%2!=0)
    {
        if(myrank==0)
        {
            std::cout<<"Error: ParalelTempering::number of processes not a multiple of 2";
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    
    data= Eigen::MatrixXd(parameters.getDataPoints()*parameters.getBins(),12);
    temperatureList = Eigen::VectorXi(nproc);
    rankList = Eigen::VectorXi(nproc);
   
    Jnnxy=parameters.getJnnxy();
    
    Jnnz1=parameters.getJnnz1();
    Jnnz2=parameters.getJnnz2();
    Jnnn=parameters.getJnnn();
    D=parameters.getD();
    B=parameters.getBstart();
    T=parameters.getTstart();
    Nspins = crystal.getNspins();
    Elist = Eigen::VectorXd::Zero(Nspins);
    
    
    fieldDirection=parameters.getFieldDirection();
    muNorm = parameters.getMuNorm();
    NNlist=crystal.getNNlist();
    XYNNlist=crystal.getXYNNlist();
    ZNNlist=crystal.getZNNlist();
    Z1NNlist=crystal.getZ1NNlist();
    Z2NNlist=crystal.getZ2NNlist();
    NNNlist=crystal.getNNNlist();
    
    
    
    spinTetraList=crystal.getSpinTetraList();
    tetraSpinList=crystal.getTetraSpinList();
    
    ewaldMatrix=ewald.getEwaldMatrix();
    
    rng_s.read_seed_local(parameters.getPathToDataFolder()+"randomSeed/"+"inputSeed_MC.txt");
    
    
    
    
    
    
}




void ParallelTempering::findElist(){
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
            Elist(i)+= ewaldMatrix(i,j)*(*spinOrientation)(i)*(*spinOrientation)(j);
        }
        //Field-------------------------------------------------------------------------------
        Elist(i) += -(*spinOrientation)(i)*muNorm*B*(spinDirection.row(i)).dot(fieldDirection);
    }
}


















Eigen::MatrixXd ParallelTempering::run()
{//Same as ::run() but with myrank input for use in MPI applications.
    
    
    if (myrank == 0){
        //seed with different random numbers for each rank.
        rng_s.read_seed_local(parameters.getPathToDataFolder()+"randomSeed/"+"inputSeed_MC.txt");
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
    
   
    MCsimRun();
    
    return data;
}















void ParallelTempering::MCsimRun()
{//This void carries out the metropolis algorithm on the lattice. It makes nproc data points and uses parallel tempering between them
    
    
    double MPIstartTime = MPI_Wtime();
    double relativeTime = MPI_Wtime()-MPIstartTime;
    double runTime = parameters.getRunTime();
    double eqPercentage= parameters.getEqPercentage();
    double Tstart=       parameters.getTstart();
    double Tend=         parameters.getTend();
    bins=parameters.getBins();
    stepsAcceptedInBin=0;
    loopStepsAcceptedInBin=0;
    
    
    
    
    
    if(parameters.getLogSteps())
    {
         T=(Tstart-Tend)*log(myrank+1)/(double)log(nproc+1)+Tend;
        
    }else
    {
         T=(Tstart-Tend)*myrank/(double)nproc+Tend;
    }
    
    
    
    temperatureList(myrank)=myrank;
    
    
    findTemperatureList();
    
    //Initializing "directly" observable quantites-----------------
    M=Magnetization();
    findElist();
    E=Energy();
    //Initializing the local sums used in each bin-------------------
    double Esum=0;
    double EsqrSum=0;
    Eigen::Vector3d MvecSum(0,0,0);
    double Msum=0;
    double MsqrSum=0;
    //Equilibration----------------------------------------------------------------
    MPIstartTime = MPI_Wtime();
    relativeTime=0;
    while(relativeTime<runTime*eqPercentage/100)
    {
        
        relativeTime= MPI_Wtime()-MPIstartTime;
        for(int i=0; i<100;++i)
        {
            E+=singleSweepList(); //spin flip
            E+=loopSweepList();   //loop flip
        }
    }
    //-----------------------------------------------------------------------------
    
    
    
    
    double binStartTime= MPI_Wtime();
    double relativeTimeBin=0;
    int bin=0;
    long int stepsPerBin=0;
    int binaryNumber=0;
    while(relativeTime<runTime)
    {
        relativeTime= MPI_Wtime()-MPIstartTime;
        MvecSum<<0,0,0;
        Msum=0;//resetting all local sums for each bin
        Esum=0;
        MsqrSum=0;
        EsqrSum=0;
        stepsAcceptedInBin=0;
        loopStepsAcceptedInBin=0;
        stepsPerBin=0;
        std::cout<<"relativeTime= "<<relativeTime<<"\n";
        
        
        binStartTime= MPI_Wtime();
        relativeTimeBin=0;
        
        replicaExchange(++binaryNumber%2);
        while(relativeTimeBin<runTime/(double)bins)
        {
            relativeTime= MPI_Wtime()-MPIstartTime;
            std::cout<<"relativeTime= "<<relativeTime<<"\n";
            std::cout<<"\t\t\t relativeTimeBin= "<<relativeTimeBin<<"\n";
            std::cout<<"\t\t\t\t\t\t"<<myrank<<"ttemperatureList= "<<temperatureList.transpose()<<"\n";
            if(myrank==0) std::cout<<"\t\t\t\t\t\t\t\t\t"<<"T_rank0= "<<T<<"\n";
            relativeTimeBin= MPI_Wtime()-binStartTime;
            for(long int i=0; i<10000;++i)
            {
                E+=singleSweepList(); //spin flip
                E+=loopSweepList();   //loop flip
                Esum+=E;//summing of local quantities in order to create an average
                EsqrSum+=pow(E,2);
                MvecSum+=M;
                Msum+=M.norm();
                MsqrSum+=pow(M.norm(),2);
                stepsPerBin++;
            }
        }
        data(bin,0)=T;
        data(bin,1)=B;
        data(bin,2)=Esum/stepsPerBin;
        data(bin,3)=EsqrSum/stepsPerBin;
        data(bin,4)=Msum/stepsPerBin;
        data(bin,5)=MsqrSum/stepsPerBin;
        data(bin,6)=MvecSum(0)/stepsPerBin;
        data(bin,7)=MvecSum(1)/stepsPerBin;
        data(bin,8)=MvecSum(2)/stepsPerBin;
        data(bin,9)=stepsAcceptedInBin/((double)stepsPerBin*Nspins);
        data(bin,10)=loopStepsAcceptedInBin/((double)stepsPerBin*Nspins);
        data(bin,11)=Jnnz1;
        bin++;
        
      
        if(myrank==0){
            std::cout<<myrank<<", changeNow------------------------------------------------------!!\n";
            std::cout<<relativeTime<<"\n";
            std::cout<<temperatureList<<"\n";
            //std::cout<<"T0= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
        }
//        if(myrank==1){
//            std::cout<<"\t\t\tT1= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==2){
//            std::cout<<"\t\t\t\t\t\tT2= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==3){
//            std::cout<<"\t\t\t\t\t\t\t\t\tT3= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==4){
//            std::cout<<"\t\t\t\t\t\t\t\t\t\t\t\tT4= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==5){
//            std::cout<<"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT5= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==6){
//            std::cout<<"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT6= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==7){
//            std::cout<<"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT7= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
//        if(myrank==8){
//            std::cout<<"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT8= "<<T<<", E= "<<E<<", time="<<MPI_Wtime()-MPIstartTime<<"\n";
//        }
   
        
   
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    
}













Eigen::RowVector3d ParallelTempering::Magnetization()
{
    Eigen::RowVector3d Mag(0,0,0);
    for(int n=0; n<lattice.rows(); n++)
    {
        Mag += spinDirection.row(n)*(*spinOrientation)(n);
    }
    return Mag;
}
double ParallelTempering::Energy()
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



double ParallelTempering::singleSweepList()
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
                if(i!=k) Elist(k)+= -2*ewaldMatrix(i,k)*(*spinOrientation)(i)*(*spinOrientation)(k);
            }
            M+= -2*spinDirection.row(i)*(*spinOrientation)(i);
            (*spinOrientation)(i)=-(*spinOrientation)(i);
            stepsAcceptedInBin++;
        }
    }
    return dE;
}














double ParallelTempering::loopSweepList()
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
                        loopEnergy+=  -ewaldMatrix(i,j)*(*spinOrientation)(i)*(*spinOrientation)(j);
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
                        if(i!=k) Elist(k)+= -2*ewaldMatrix(i,k)*(*spinOrientation)(i)*(*spinOrientation)(k);
                    }
                    M+= -2*spinDirection.row(i)*(*spinOrientation)(i);
                    (*spinOrientation)(i)=-(*spinOrientation)(i);
                }
            }
        }
    }
    return dE;
}









bool ParallelTempering::findLoop()
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














void ParallelTempering::replicaExchange(int oddEven){
    //This void carries out the replica exchange update for two adjacent temperatures. If exp(abs(dE)/T)<r is met, it makes a swap of the temperatures of the two processes. it compares pairwise odd and even temperatures
    int n1rank=rankList((temperatureList(myrank)-1+nproc)%nproc);
    int n2rank=rankList((temperatureList(myrank)+1)%nproc);
    double n2E=NAN;
    double n2T=NAN;
    Eigen::MatrixXd n2data(data.rows(), data.cols());
    std::cout<<"replicaExchange called!--------------------------------------------------------\n";
    int accepted=0;
    if((temperatureList(myrank)+1)%2==oddEven)
    {
        MPI_Send(&E, 1, MPI_DOUBLE, n1rank, 11, MPI_COMM_WORLD);
        MPI_Send(&T, 1, MPI_DOUBLE, n1rank, 12, MPI_COMM_WORLD);
        MPI_Send(data.data(), (int) data.size(), MPI_DOUBLE, n1rank, 13, MPI_COMM_WORLD);
        MPI_Recv(&E, 1, MPI_DOUBLE, n1rank, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&T, 1, MPI_DOUBLE, n1rank, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(data.data(), (int) data.size(), MPI_DOUBLE, n1rank, 16, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }else
    {
        MPI_Recv(&n2E, 1, MPI_DOUBLE, n2rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&n2T, 1, MPI_DOUBLE, n2rank, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n2data.data(), (int) n2data.size(), MPI_DOUBLE, n2rank, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(true)//replica exchange criteria.
        {
            MPI_Send(&E, 1, MPI_DOUBLE, n2rank, 14, MPI_COMM_WORLD);
            MPI_Send(&T, 1, MPI_DOUBLE, n2rank, 15, MPI_COMM_WORLD);
            MPI_Send(data.data(), (int) data.size(), MPI_DOUBLE, n2rank, 16, MPI_COMM_WORLD);
            T=n2T;
            data=n2data;
            accepted=1;
        }else
        {
            MPI_Send(&n2E, 1, MPI_DOUBLE, n2rank, 14, MPI_COMM_WORLD);
            MPI_Send(&n2T, 1, MPI_DOUBLE, n2rank, 15, MPI_COMM_WORLD);
            MPI_Send(n2data.data(), (int) n2data.size(), MPI_DOUBLE, n2rank, 16, MPI_COMM_WORLD);
        }
        
        if(accepted==1&&myrank==0)
        {
            temperatureList(0)=temperatureList(n2rank);
            temperatureList(n2rank)=temperatureList(0);
            rankList(temperatureList(0))=0;
            rankList(temperatureList(n2rank))=n2rank;
        }
        if(myrank!=0) MPI_Send(&accepted, 1, MPI_INT, 0, 17, MPI_COMM_WORLD);
    }
    if(myrank==0)
    {
        for(int Ti=oddEven; Ti<nproc; Ti=Ti+2)
        {
            if(rankList(Ti)!=0) MPI_Recv(&accepted, 1, MPI_INT, rankList(Ti), 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(accepted==1)
            {
                int temporary=rankList(Ti);
                rankList(Ti)=rankList((Ti+1)%nproc);
                rankList((Ti+1)%nproc)=temporary;
            }
        }
        
        for(int i=0; i<nproc; ++i)
        {
            temperatureList(rankList(i))=i;
        }
    }
    MPI_Bcast(temperatureList.data(), (int)temperatureList.size(), MPI_INT, 0, MPI_COMM_WORLD);
    
    if(myrank>0){
        for(int i=0; i<nproc; ++i)
        {
            rankList(temperatureList(i))=i;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}






























void ParallelTempering::findTemperatureList()
{
    
    
    if(myrank==0)
    {
        for(int rankNumber=1; rankNumber<nproc; ++rankNumber)
        {
            MPI_Recv(&temperatureList(rankNumber), 1, MPI_INT, rankNumber, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }else
    {
        MPI_Send(&temperatureList(myrank), 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
    }
    
    
    
    MPI_Bcast(temperatureList.data(), (int) temperatureList.size(), MPI_INT, 0, MPI_COMM_WORLD);
   
    
    
    
    for(int rankNumber=0; rankNumber<nproc; ++rankNumber)
    {
        rankList(temperatureList(rankNumber))=rankNumber;
    }
    
    
   
  
   
    for(int rankNumber=0; rankNumber<nproc; ++rankNumber)
    {
        rankList(temperatureList(rankNumber))=rankNumber;
    }
    

}






















bool ParallelTempering::isEqual(double a, double b, double tol)
{
    return std::abs(a-b)<tol;
}



void ParallelTempering::saveStateToFile(std::string path, int i)
{
    //This void the current spin state to a file state_i
    
    
    std::string filePathAndName=path+parameters.getOutputFilename()+"_state_"+std::to_string(i)+".txt";
    std::ofstream file(filePathAndName);
    
    //exit if unable to create file
    if(!file){
        std::cerr<<"File could not be opened"<<"\n";
        exit(EXIT_FAILURE);
    }
    file <<std::setprecision(16)<< (*spinOrientation);
    std::cout<<"File Saved"<<"\n";
}














