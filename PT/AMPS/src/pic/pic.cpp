//====================================================
//$Id$
//====================================================
//the general functions for the pic solver

#include "pic.h"


//sampling variables

long int PIC::LastSampleLength=0,PIC::CollectingSampleCounter=0,PIC::RequiredSampleLength=100,PIC::DataOutputFileNumber=0;
int PIC::SamplingMode=_RESTART_SAMPLING_MODE_;

//====================================================
//perform one time step
void PIC::TimeStep() {
   double ParticleMovingTime,InjectionBoundaryTime,ParticleExchangeTime,IterationExecutionTime,SamplingTime,StartTime=MPI_Wtime();

   //Collect and exchange the run's statictic information
   static const int nRunStatisticExchangeIterationsMin=5,nRunStatisticExchangeIterationsMax=500,nRunStatisticExchangeTime=120;
   static long int nTotalIterations=0,nInteractionsAfterRunStatisticExchange=0;
   static int nExchangeStatisticsIterationNumberSteps=10;

   nTotalIterations++;
   nInteractionsAfterRunStatisticExchange++;
   PIC::Parallel::IterationNumberAfterRebalancing++;

  //sampling of the particle data
  PIC::Sampling::Sampling();
  SamplingTime=MPI_Wtime()-StartTime;

  //injection boundary conditions
  InjectionBoundaryTime=MPI_Wtime();
  PIC::BC::InjectionBoundaryConditions();
  InjectionBoundaryTime=MPI_Wtime()-InjectionBoundaryTime;

  //move existing particles
  ParticleMovingTime=MPI_Wtime();
  PIC::Mover::MoveParticles();
  ParticleMovingTime=MPI_Wtime()-ParticleMovingTime;

  //check the consistence of the particles lists
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  PIC::ParticleBuffer::CheckParticleList();
#endif

  IterationExecutionTime=MPI_Wtime()-StartTime;

  //syncrinie processors and exchnge particle data
  ParticleExchangeTime=MPI_Wtime();
  PIC::Parallel::ExchangeParticleData();
  ParticleExchangeTime=MPI_Wtime()-ParticleExchangeTime;


  struct cExchangeStatisticData {
    double TotalInterationRunTime;
    double IterationExecutionTime;
    long int TotalParticlesNumber;
    double ParticleExchangeTime;
    double SamplingTime;
    double InjectionBoundaryTime;
    double ParticleMovingTime;
    double Latency;
    long int recvParticleCounter,sendParticleCounter;
    long int nInjectedParticles;
  };

  struct cStatExchangeControlParameters {
    int nExchangeStatisticsIterationNumberSteps;
    bool WallTimeExeedsLimit;
  };


  if (nInteractionsAfterRunStatisticExchange==nExchangeStatisticsIterationNumberSteps) { //collect and exchenge the statistical data of the run
    cExchangeStatisticData localRunStatisticData;
    int thread;
    cExchangeStatisticData *ExchangeBuffer=NULL;

    localRunStatisticData.TotalInterationRunTime=MPI_Wtime()-StartTime;
    localRunStatisticData.IterationExecutionTime=IterationExecutionTime;
    localRunStatisticData.TotalParticlesNumber=PIC::ParticleBuffer::NAllPart;
    localRunStatisticData.ParticleExchangeTime=ParticleExchangeTime;
    localRunStatisticData.SamplingTime=SamplingTime;
    localRunStatisticData.ParticleMovingTime=ParticleMovingTime;
    localRunStatisticData.InjectionBoundaryTime=InjectionBoundaryTime;
    localRunStatisticData.Latency=PIC::Parallel::Latency;
    localRunStatisticData.recvParticleCounter=PIC::Parallel::recvParticleCounter;
    localRunStatisticData.sendParticleCounter=PIC::Parallel::sendParticleCounter;
    localRunStatisticData.nInjectedParticles=PIC::BC::nInjectedParticles;


    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf("\nPIC: (%i/%i %i:%i:%i), Iteration: %ld  (%ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,nTotalIterations,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);



      ExchangeBuffer=new cExchangeStatisticData[PIC::Mesh::mesh.nTotalThreads];
      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_COMM_WORLD);


      //output the data
      long int nTotalModelParticles=0,nTotalInjectedParticels=0;
      double MinExecutionTime=localRunStatisticData.IterationExecutionTime,MaxExecutionTime=localRunStatisticData.IterationExecutionTime,MaxLatency=0.0,MeanLatency=0.0;

      printf("Description:\n");
      printf("1:\t Thread\n");
      printf("2:\t Total Particle's number\n");
      printf("3:\t Total Interation Time\n");
      printf("4:\t Iteration Execution Time\n");
      printf("5:\t Sampling Time\n");
      printf("6:\t Injection Boundary Time\n");
      printf("7:\t Particle Moving Time\n");
      printf("8:\t Particle Exchange Time\n");
      printf("9:\t Latency\n");
      printf("10:\t Send Particles\n");
      printf("11:\t Recv Particles\n");
      printf("12:\t nInjected Particls\n");

      printf("1\t 2\t 3\t\t 4\t\t 5\t\t 6\t\t 7\t\t 8\t\t 9\t\t 10\t 11\t 12\n");

//      printf("Thread, Total Particle's number, Total Interation Time, Iteration Execution Time, Sampling Time, Injection Boundary Time, Particle Moving Time, Particle Exchange Time, Latency, Send Particles, Recv Particles, nInjected Particls\n");

      for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
        printf("%i\t %ld\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %ld\t %ld\t %ld\n",thread,ExchangeBuffer[thread].TotalParticlesNumber,ExchangeBuffer[thread].TotalInterationRunTime,
            ExchangeBuffer[thread].IterationExecutionTime,ExchangeBuffer[thread].SamplingTime,ExchangeBuffer[thread].InjectionBoundaryTime,ExchangeBuffer[thread].ParticleMovingTime,
            ExchangeBuffer[thread].ParticleExchangeTime,ExchangeBuffer[thread].Latency,ExchangeBuffer[thread].sendParticleCounter,
            ExchangeBuffer[thread].recvParticleCounter,ExchangeBuffer[thread].nInjectedParticles);

        nTotalModelParticles+=ExchangeBuffer[thread].TotalParticlesNumber;
        nTotalInjectedParticels+=ExchangeBuffer[thread].nInjectedParticles;
        MeanLatency+=ExchangeBuffer[thread].Latency;

        if (MinExecutionTime>ExchangeBuffer[thread].IterationExecutionTime) MinExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxExecutionTime<ExchangeBuffer[thread].IterationExecutionTime) MaxExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxLatency<ExchangeBuffer[thread].Latency) MaxLatency=ExchangeBuffer[thread].Latency;
      }

      MeanLatency/=PIC::Mesh::mesh.nTotalThreads;
      PIC::Parallel::CumulativeLatency+=MeanLatency*nInteractionsAfterRunStatisticExchange;

      printf("Total number of particles: %ld\n",nTotalModelParticles);
      printf("Total number of injected particles: %ld\n",nTotalInjectedParticels);
      printf("Iteration Execution Time: min=%e, max=%e\n",MinExecutionTime,MaxExecutionTime);
      printf("Latency: max=%e,mean=%e;CumulativeLatency=%e,Iterations after rebalabcing=%ld,Rebalancing Time=%e\n",MaxLatency,MeanLatency,PIC::Parallel::CumulativeLatency,PIC::Parallel::IterationNumberAfterRebalancing,PIC::Parallel::RebalancingTime);

      //check the elapsed walltime for the alarm
      if (PIC::Alarm::AlarmInitialized==true) {
        if (MPI_Wtime()-PIC::Alarm::StartTime>PIC::Alarm::RequestedExecutionWallTime) PIC::Alarm::WallTimeExeedsLimit=true;
        printf("Execution walltime: %e sec\n",MPI_Wtime()-PIC::Alarm::StartTime);
      }


      //determine the new value for the interation number between exchanging of the run stat data
      if (localRunStatisticData.TotalInterationRunTime>0.0) {
        nExchangeStatisticsIterationNumberSteps=(int)(nRunStatisticExchangeTime/localRunStatisticData.TotalInterationRunTime);
        if (nExchangeStatisticsIterationNumberSteps<nRunStatisticExchangeIterationsMin) nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMin;
        if (nExchangeStatisticsIterationNumberSteps>nRunStatisticExchangeIterationsMax) nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMax;
      }
      else nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMin;

      cStatExchangeControlParameters StatExchangeControlParameters;
      StatExchangeControlParameters.nExchangeStatisticsIterationNumberSteps=nExchangeStatisticsIterationNumberSteps;
      StatExchangeControlParameters.WallTimeExeedsLimit=PIC::Alarm::WallTimeExeedsLimit;

      MPI_Bcast(&StatExchangeControlParameters,sizeof(cStatExchangeControlParameters),MPI_CHAR,0,MPI_COMM_WORLD);


      delete [] ExchangeBuffer;
    }
    else {
      cStatExchangeControlParameters StatExchangeControlParameters;

      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_COMM_WORLD);
      MPI_Bcast(&StatExchangeControlParameters,sizeof(cStatExchangeControlParameters),MPI_CHAR,0,MPI_COMM_WORLD);

      nExchangeStatisticsIterationNumberSteps=StatExchangeControlParameters.nExchangeStatisticsIterationNumberSteps;
      PIC::Alarm::WallTimeExeedsLimit=StatExchangeControlParameters.WallTimeExeedsLimit;
    }

    //collect the run time execution statistics from individual models
    if ((PIC::ExchangeExecutionStatisticsFunctions.size()!=0)&&(nInteractionsAfterRunStatisticExchange!=0)) {
      CMPI_channel pipe(10000);
      vector<PIC::fExchangeExecutionStatistics>::iterator fptr;

      if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
      else pipe.openSend(0);

      for (fptr=PIC::ExchangeExecutionStatisticsFunctions.begin();fptr!=PIC::ExchangeExecutionStatisticsFunctions.end();fptr++) (*fptr)(&pipe,nInteractionsAfterRunStatisticExchange);

      if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
      else pipe.closeSend();
    }

    //redistribute the processor load and check the mesh afterward
#if _PIC_EMERGENCY_LOAD_REBALANCING_MODE_ == _PIC_MODE_ON_
    int EmergencyLoadRebalancingFlag=false;

    if (PIC::Mesh::mesh.ThisThread==0) if (PIC::Parallel::CumulativeLatency>PIC::Parallel::EmergencyLoadRebalancingFactor*PIC::Parallel::RebalancingTime) EmergencyLoadRebalancingFlag=true;

    MPI_Bcast(&EmergencyLoadRebalancingFlag,1,MPI_INT,0,MPI_COMM_WORLD);

    if (EmergencyLoadRebalancingFlag==true) {
      if (PIC::Mesh::mesh.ThisThread==0) printf("Load Rebalancing.....  begins\n");


      MPI_Barrier(MPI_COMM_WORLD);
      PIC::Parallel::RebalancingTime=MPI_Wtime();

      PIC::Mesh::mesh.CreateNewParallelDistributionLists(_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
      PIC::Parallel::IterationNumberAfterRebalancing=0,PIC::Parallel::CumulativeLatency=0.0;

      MPI_Barrier(MPI_COMM_WORLD);
      PIC::Parallel::RebalancingTime=MPI_Wtime()-PIC::Parallel::RebalancingTime;
      if (PIC::Mesh::mesh.ThisThread==0) printf("Load Rebalancing.....  done\n");
    }
#elif _PIC_EMERGENCY_LOAD_REBALANCING_MODE_ == _PIC_MODE_OFF_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    nInteractionsAfterRunStatisticExchange=0;
  }

}
//====================================================
//the general sampling procedure
void PIC::Sampling::Sampling() {
  int s,i,j,k,idim;
  long int LocalCellNumber,ptr,ptrNext;


  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  char *SamplingData;
  double *v,LocalParticleWeight,Speed2,v2;

  double *sampledVelocityOffset,*sampledVelocity2Offset;


  //temporary buffer for sampled data
  char  tempSamplingBuffer[PIC::Mesh::sampleSetDataLength]; //[320];

//  double *tempSamplingBuffer=new double [PIC::Mesh::sampleSetDataLength];  //// (double*)alloca(sizeof(double)*PIC::Mesh::sampleSetDataLength);
//  if (tempSamplingBuffer==NULL) exit(__LINE__,__FILE__,"Error: cannot allocate the buffer in teh stack");

  //temporaty buffer to store the copy of the particle
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];   /////[73]; //[PIC::ParticleBuffer::ParticleDataLength];


  //global offsets

  const int sampledParticleWeghtRelativeOffset=PIC::Mesh::sampledParticleWeghtRelativeOffset;
  const int sampledParticleNumberRelativeOffset=PIC::Mesh::sampledParticleNumberRelativeOffset;
  const int sampledParticleNumberDensityRelativeOffset=PIC::Mesh::sampledParticleNumberDensityRelativeOffset;
  const int sampledParticleVelocityRelativeOffset=PIC::Mesh::sampledParticleVelocityRelativeOffset;
  const int sampledParticleVelocity2RelativeOffset=PIC::Mesh::sampledParticleVelocity2RelativeOffset;
  const int sampleSetDataLength=PIC::Mesh::sampleSetDataLength;
  const int sampledParticleSpeedRelativeOffset=PIC::Mesh::sampledParticleSpeedRelativeOffset;
  const int collectingCellSampleDataPointerOffset=PIC::Mesh::collectingCellSampleDataPointerOffset;
  const int ParticleDataLength=PIC::ParticleBuffer::ParticleDataLength;



  //parallel efficientcy measure
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
  long int TreeNodeTotalParticleNumber;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double TreeNodeProcessingTime;
#endif


  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndex[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  int centerNodeIndexCounter;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          centerNodeIndex[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        }

        if (DIM==1) break;
      }

      if ((DIM==1)||(DIM==2)) break;
    }
  }

#if _PIC_SAMPLE_PARTICLE_DATA_MODE_ == _PIC_SAMPLE_PARTICLE_DATA_MODE__BETWEEN_ITERATIONS_
  //go through the 'local nodes'
  while (node!=NULL) {
    block=node->block;






#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    TreeNodeTotalParticleNumber=0;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    TreeNodeProcessingTime=MPI_Wtime();
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_OFF_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the otion is not recognized");
#endif

    //sample the distrivution function
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
    if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::SampleDistributionFnction(node);
#endif



    {
      {
        for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
            LocalCellNumber=centerNodeIndex[centerNodeIndexCounter];
            cell=block->GetCenterNode(LocalCellNumber);
            ptr=cell->FirstCellParticle;

            //===================    DEBUG ==============================

//            if (cell->Temp_ID==56119) {
//              cout << __LINE__ << __FILE__ << endl;
//            }

            //====================   END DEBUG ==========================


            if (ptr!=-1) {

            SamplingData=cell->GetAssociatedDataBufferPointer()+/*PIC::Mesh::*/collectingCellSampleDataPointerOffset;
            memcpy((void*)tempSamplingBuffer,(void*)SamplingData,/*PIC::Mesh::*/sampleSetDataLength);

            ptr=cell->FirstCellParticle;
            ptrNext=ptr;

            //===================    DEBUG ==============================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
            if (cell->Measure<=0.0) {
              cout << __FILE__<< __LINE__ << endl;
              exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
            }
#endif
            //===================   END DEBUG ==============================

            while (ptrNext!=-1) {
              ptr=ptrNext;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
              TreeNodeTotalParticleNumber++;
#endif




              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

              memcpy((void*)tempParticleData,(void*)ParticleData,/*PIC::ParticleBuffer::*/ParticleDataLength);


              Speed2=0.0;

              s=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)tempParticleData); ///ParticleData);
              v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)tempParticleData); ///ParticleData);

              ptrNext=PIC::ParticleBuffer::GetNext((PIC::ParticleBuffer::byte*)tempParticleData);  //ParticleData);


              LocalParticleWeight=block->GetLocalParticleWeight(s);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)tempParticleData); //ParticleData);


//              move sample memory into cash


              /*

              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberRelativeOffset))+=1;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;


              sampledVelocityOffset=3*s+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocityRelativeOffset);
              sampledVelocity2Offset=3*s+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocity2RelativeOffset);

              #pragma unroll(3)
              for (idim=0;idim<3;idim++) {
                v2=v[idim]*v[idim];
                Speed2+=v2;

//                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocityRelativeOffset))+=v[idim]*LocalParticleWeight;
//                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocity2RelativeOffset))+=v2*LocalParticleWeight;

                *(idim+sampledVelocityOffset)+=v[idim]*LocalParticleWeight;
                *(idim+sampledVelocity2Offset)+=v2*LocalParticleWeight;

              }

              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleSpeedRelativeOffset))+=sqrt(Speed2)*LocalParticleWeight;

           */



              *(s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
              *(s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleNumberRelativeOffset))+=1;
              *(s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;


              sampledVelocityOffset=3*s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleVelocityRelativeOffset);
              sampledVelocity2Offset=3*s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleVelocity2RelativeOffset);

              for (idim=0;idim<3;idim++) {
                v2=v[idim]*v[idim];
                Speed2+=v2;

//                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocityRelativeOffset))+=v[idim]*LocalParticleWeight;
//                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocity2RelativeOffset))+=v2*LocalParticleWeight;

                *(idim+sampledVelocityOffset)+=v[idim]*LocalParticleWeight;
                *(idim+sampledVelocity2Offset)+=v2*LocalParticleWeight;

              }

              *(s+(double*)(tempSamplingBuffer+/*PIC::Mesh::*/sampledParticleSpeedRelativeOffset))+=sqrt(Speed2)*LocalParticleWeight;


            }

            memcpy((void*)SamplingData,(void*)tempSamplingBuffer,/*PIC::Mesh::*/sampleSetDataLength);

            }
          }
       }
    }


    //Sample the parallel load: the total particle number
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    node->ParallelLoadMeasure+=TreeNodeTotalParticleNumber;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    node->ParallelLoadMeasure+=MPI_Wtime()-TreeNodeProcessingTime;
#endif


    node=node->nextNodeThisThread;
  }
#elif _PIC_SAMPLE_PARTICLE_DATA_MODE_ == _PIC_SAMPLE_PARTICLE_DATA_MODE__DURING_PARTICLE_MOTION_
  //do nothing
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

  //sample user defined data
  if (PIC::IndividualModelSampling::SamplingProcedure.size()!=0) {
    int nfunc,nfuncTotal=PIC::IndividualModelSampling::SamplingProcedure.size();

    for (nfunc=0;nfunc<nfuncTotal;nfunc++) PIC::IndividualModelSampling::SamplingProcedure[nfunc]();
  }

  //Increment the sample length
  CollectingSampleCounter++;

  //Output the data flow file if needed
  if ((CollectingSampleCounter==RequiredSampleLength)||(PIC::Alarm::WallTimeExeedsLimit==true)) {
    //exchnge the sampling data
    PIC::Mesh::mesh.ParallelBlockDataExchange();

    //check different sampling modes
    if (SamplingMode==_RESTART_SAMPLING_MODE_) {
      PIC::Mesh::switchSamplingBuffers();
      LastSampleLength=CollectingSampleCounter;
      CollectingSampleCounter=0;

      for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        block=node->block;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
              PIC::Mesh::flushCollectingSamplingBuffer(block->GetCenterNode(LocalCellNumber));
            }

            if (DIM==1) break;
          }

          if ((DIM==1)||(DIM==2)) break;
        }
      }

      //flush sampling buffers in the internal surfaces installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
  #elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();
//      PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere* Sphere;

      cInternalSphericalData *Sphere;

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
//        Sphere=(PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere*)PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->GetSurfaceDataPointer();

        Sphere=PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface);

        PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
      }


  #else
      exit(__LINE__,__FILE__,"Error: unknown option");
  #endif


    }
    else if (SamplingMode==_ACCUMULATE_SAMPLING_MODE_) {
      exit(__LINE__,__FILE__,"Error: the sampling mode '_ACCUMULATE_SAMPLING_MODE_' is not implemented");
    }
    else exit(__LINE__,__FILE__,"Error: the sampling mode is not defined");

    //print output file
    char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld.dat",ChemSymbol,s,DataOutputFileNumber);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("done.\n");
        fflush(stdout);
      }

      //print the sampled distribution function into a file
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
      if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) {
        sprintf(fname,"pic.distribution.%s.s=%i.out=%ld",ChemSymbol,s,DataOutputFileNumber);
        PIC::DistributionFunctionSample::printDistributionFunction(fname,s);
      }
#endif

      //print into the file the sampled data of the internal surfaces installed into the mesh
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
#elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
        sprintf(fname,"pic.Sphere=%ld.%s.s=%i.out=%ld.dat",iSphericalSurface,ChemSymbol,s,DataOutputFileNumber);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("printing output file: %s.........",fname);
          fflush(stdout);
        }

        PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->PrintSurfaceData(fname,s);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("done.\n");
          fflush(stdout);
        }
      }


#else
      exit(__LINE__,__FILE__,"Error: unknown option");
#endif
    }

    //print the statistic information for the run
    CMPI_channel pipe(1000000);
    long int nTotalSimulatedParticles;

    nTotalSimulatedParticles=PIC::ParticleBuffer::NAllPart;

    if (PIC::Mesh::mesh.ThisThread!=0) {
      pipe.openSend(0);

      pipe.send(nTotalSimulatedParticles);

      pipe.closeSend();
    }
    else {
      pipe.openRecvAll();

      for (int thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
        nTotalSimulatedParticles+=pipe.recv<long int>(thread);
      }

      printf("Model run statistics:\n");
      printf("The total number of model particles: %ld\n",nTotalSimulatedParticles);

      pipe.closeRecvAll();
    }

    //Finish the code execution if the walltime exeeds the limit
    if (PIC::Alarm::WallTimeExeedsLimit==true) PIC::Alarm::FinishExecution();

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
    if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::flushSamplingBuffers();
#endif

    //increment the output file number
    DataOutputFileNumber++;

    //redistribute the processor load and check the mesh afterward
#if _PIC_SAMPLING_BREAK_LOAD_REBALANCING_MODE_ == _PIC_MODE_ON_
    if (PIC::Parallel::IterationNumberAfterRebalancing!=0) {
      MPI_Barrier(MPI_COMM_WORLD);
      PIC::Parallel::RebalancingTime=MPI_Wtime();

      PIC::Mesh::mesh.CreateNewParallelDistributionLists(_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
      PIC::Parallel::IterationNumberAfterRebalancing=0,PIC::Parallel::CumulativeLatency=0.0;

      MPI_Barrier(MPI_COMM_WORLD);
      PIC::Parallel::RebalancingTime=MPI_Wtime()-PIC::Parallel::RebalancingTime;
    }
#elif _PIC_SAMPLING_BREAK_LOAD_REBALANCING_MODE_ == _PIC_MODE_OFF_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif


//=====================  DEBUG   ========================
    /*
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    //output the data into a file to check the new cells' distribution
    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld-redistributed-load-CompleteSample.dat",ChemSymbol,s,DataOutputFileNumber-1);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("done.\n");
        fflush(stdout);
      }
    }

    PIC::Mesh::switchSamplingBuffers();

    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld-redistributed-load-TempSample.dat",ChemSymbol,s,DataOutputFileNumber-1);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("done.\n");
        fflush(stdout);
      }
    }

    PIC::Mesh::switchSamplingBuffers();
#endif
*/
//=====================  END DEBUG ===========================


  }





}


//====================================================
//the run time signal and exeptions handler
void PIC::SignalHandler(int sig) {
  cout << "Signal is intersepted: thread=" << PIC::Mesh::mesh.ThisThread << endl;

  switch(sig) {
  case SIGFPE :
    cout << "Signal=SIGFPE" << endl;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown signal");
  }

  exit(__LINE__,__FILE__,"Error: exit in the signal handler");
}
//====================================================
//init the particle solver
void PIC::Init() {

  //init MPI variables
  MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
  MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);

  //init sections of the particle solver

  //set up the signal handler
  signal(SIGFPE,SignalHandler);


  //init ICES
#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
#define _PIC_INIT_ICES_
#endif

#ifdef _PIC_INIT_ICES_
  PIC::ICES::Init();
#endif
}



//====================================================
//set up the particle weight and time step
void PIC::Mesh::cDataBlockAMR::SetLocalParticleWeight(double weight, int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  *(spec+(double *)(associatedDataPointer+LocalParticleWeightOffset))=weight;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  if (PIC::ParticleWeightTimeStep::GlobalParticleWeight==NULL) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]>weight)) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]=weight;
  }

  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalParticleWeight(int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalParticleWeightOffset);
  return *res;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  return PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

//set up the particle time step
void PIC::Mesh::cDataBlockAMR::SetLocalTimeStep(double dt, int spec) {
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  *(spec+(double *)(associatedDataPointer+LocalTimeStepOffset))=dt;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_

  if (PIC::ParticleWeightTimeStep::GlobalTimeStep==NULL) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]>dt)) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]=dt;
  }


  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalTimeStep(int spec) {

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalTimeStepOffset);
  return *res;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  return  PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif


}




