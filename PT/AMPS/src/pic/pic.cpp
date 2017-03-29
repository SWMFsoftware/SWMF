//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the general functions for the pic solver

#include "pic.h"
#include "PhotolyticReactions.h"


//sampling variables

long int PIC::LastSampleLength=0,PIC::CollectingSampleCounter=0,PIC::DataOutputFileNumber=0;
int PIC::SamplingMode=_RESTART_SAMPLING_MODE_;
long int PIC::RequiredSampleLength=100;

//external sampling procedure
const int PIC::Sampling::ExternalSamplingLocalVariables::nMaxSamplingRoutines=128;
int PIC::Sampling::ExternalSamplingLocalVariables::SamplingRoutinesRegistrationCounter=0;
PIC::Sampling::ExternalSamplingLocalVariables::fSamplingProcessor *PIC::Sampling::ExternalSamplingLocalVariables::SamplingProcessor=NULL;
PIC::Sampling::ExternalSamplingLocalVariables::fPrintOutputFile *PIC::Sampling::ExternalSamplingLocalVariables::PrintOutputFile=NULL;


int PIC::Sampling::minIterationNumberForDataOutput=0;
bool *PIC::Sampling::SaveOutputDataFile=NULL;


//====================================================
//perform one time step
int PIC::TimeStep() {
   double UserDefinedMPI_RoutineExecutionTime=0.0,ParticleMovingTime,PhotoChemistryTime=0.0,InjectionBoundaryTime,ParticleExchangeTime,IterationExecutionTime,SamplingTime,StartTime=MPI_Wtime();
   double ParticleCollisionTime=0.0,BackgroundAtmosphereCollisionTime=0.0;
   double UserDefinedParticleProcessingTime=0.0;
   static double summIterationExecutionTime=0.0;

   //Set the exit error code
   ExitErrorCode=_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_;

   //update the local block list
   ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__UpdateBlockTable_;
   DomainBlockDecomposition::UpdateBlockTable();

   //recover the sampling data from the sampling data restart file, print the TECPLOT files and quit
   if (_PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__MODE_==_PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__MODE_ON_) {
     static bool RestartFileReadFlag=false;

     if (RestartFileReadFlag==false) {
       RestartFileReadFlag=true;

       //prepare the list of the files/speces which data fille be revobered and saved
       list<pair<string,list<int> > > SampledDataRecoveryTable;

       if (Restart::SamplingData::DataRecoveryManager==NULL) {
         //no user-defined fucntion to create the 'SampledDataRecoveringTable' is provided
         pair<string,list<int> > Table;

         for (int s=0;s<PIC::nTotalSpecies;s++) Table.second.push_back(s);

         Table.first=Restart::SamplingData::RestartFileName;
         SampledDataRecoveryTable.push_back(Table);
       }
       else Restart::SamplingData::DataRecoveryManager(SampledDataRecoveryTable,Restart::SamplingData::minReadFileNumber,Restart::SamplingData::maxReadFileNumber);

       //iterate through SampledDataRecoveryTable
       list<pair<string,list<int> > >::iterator RecoveryEntry;

       for (RecoveryEntry=SampledDataRecoveryTable.begin();RecoveryEntry!=SampledDataRecoveryTable.end();RecoveryEntry++) {
         Restart::SamplingData::Read(RecoveryEntry->first.c_str());
         MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

         for (list<int>::iterator s=RecoveryEntry->second.begin();s!=RecoveryEntry->second.end();s++) {
           char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

           PIC::MolecularData::GetChemSymbol(ChemSymbol,*s);
           sprintf(fname,"RECOVERED.%s.%s.s=%i.dat",RecoveryEntry->first.c_str(),ChemSymbol,*s);
           PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,*s);

           //preplot the recovered file if needed
           if (Restart::SamplingData::PreplotRecoveredData==true) {
             char cmd[_MAX_STRING_LENGTH_PIC_];

             sprintf(cmd,"preplot %s",fname);
             system(cmd);
           }
         }
       }

       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

       if (_PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__EXECUTION_MODE_ == _PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__EXECUTION_MODE__STOP_) {
         if (PIC::ThisThread==0) {
           printf("The sampled data from restart file/files\n");
           for (RecoveryEntry=SampledDataRecoveryTable.begin();RecoveryEntry!=SampledDataRecoveryTable.end();RecoveryEntry++) printf("%s\n",RecoveryEntry->first.c_str());
           printf("is successfully recoved. Execution is complete. See you later :-)\n");
         }

         MPI_Finalize();
         exit(EXIT_SUCCESS);
       }
     }
   }

   //recover the particle data restart file
   if (_PIC_READ_PARTICLE_DATA_RESTART_FILE__MODE_ == _PIC_READ_PARTICLE_DATA_RESTART_FILE__MODE_ON_) {
     static bool RestartFileReadFlag=false;

     if (RestartFileReadFlag==false) {
       RestartFileReadFlag=true;

       Restart::ReadParticleData(Restart::recoverParticleDataRestartFileName);
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
     }
   }

   //save the particle data restart file
   if (_PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ == _PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ON_) {
     static long int IterationCounter=0;
     static int saveRestartFileCounter=0;

     IterationCounter++;

     if (IterationCounter%Restart::ParticleRestartAutosaveIterationInterval==0) {
       char fname[_MAX_STRING_LENGTH_PIC_];

       if (Restart::ParticleDataRestartFileOverwriteMode==true) sprintf(fname,"%s",Restart::saveParticleDataRestartFileName);
       else sprintf(fname,"%s.restart=%i",Restart::saveParticleDataRestartFileName,saveRestartFileCounter);

       Restart::SaveParticleData(fname);
       saveRestartFileCounter++;

       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
     }
   }


   //Collect and exchange the run's statictic information
//   static const int nRunStatisticExchangeIterationsMin=5,nRunStatisticExchangeIterationsMax=500,nRunStatisticExchangeTime=120;
   static long int nTotalIterations=0,nInteractionsAfterRunStatisticExchange=0;
   static int nExchangeStatisticsIterationNumberSteps=10;

   nTotalIterations++;
   nInteractionsAfterRunStatisticExchange++;
   PIC::Parallel::IterationNumberAfterRebalancing++;

  //sampling of the particle data
   ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__Sampling_;
   SamplingTime=MPI_Wtime();
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   Sampling::CatchOutLimitSampledValue();
#endif
#endif

  PIC::Sampling::Sampling();

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   Sampling::CatchOutLimitSampledValue();
#endif
#endif

  SamplingTime=MPI_Wtime()-SamplingTime;
//#endif

  //injection boundary conditions
  InjectionBoundaryTime=MPI_Wtime();

  //inject particle through the domain's boundaries
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__InjectionBoundaryConditions_;
  PIC::BC::InjectionBoundaryConditions();

  //inject particles into the volume of the domain
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  if (PIC::VolumeParticleInjection::nRegistratedInjectionProcesses!=0) PIC::BC::nTotalInjectedParticles+=PIC::VolumeParticleInjection::InjectParticle();
#endif

  //call a user-defined injection function
  if (PIC::BC::UserDefinedParticleInjectionFunction!=NULL) {
    ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__UserDefinedParticleInjectionFunction_;
    PIC::BC::nTotalInjectedParticles+=PIC::BC::UserDefinedParticleInjectionFunction();
  }

  //the extra injection process by the exosphere model (src/models/exosphere)
  if (BC::ExosphereModelExtraInjectionFunction!=NULL) {
    ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__ExosphereModelExtraInjectionFunction_;
    PIC::BC::nTotalInjectedParticles+=BC::ExosphereModelExtraInjectionFunction();
  }

  InjectionBoundaryTime=MPI_Wtime()-InjectionBoundaryTime;


  //simulate particle collisions
#if _PIC__PARTICLE_COLLISION_MODEL__MODE_ == _PIC_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__ParticleCollisionModel_;
  ParticleCollisionTime=MPI_Wtime();

  #if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__NTC_
  PIC::MolecularCollisions::ParticleCollisionModel::ntc();
  #elif _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__MF_
  PIC::MolecularCollisions::ParticleCollisionModel::mf();
  #elif _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__USER_DEFINED_
  PIC::MolecularCollisions::ParticleCollisionModel::ntc();
  #else
  exit(__LINE__,__FILE__,"Error: the option is not implemented");
  #endif

  ParticleCollisionTime=MPI_Wtime()-ParticleCollisionTime;
#endif

  //simulate collisions with the background atmosphere
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__BackgroundAtmosphere_;
  BackgroundAtmosphereCollisionTime=MPI_Wtime();

  #if _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_ == _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL__PARTICLE_COLLISIONS_
  MolecularCollisions::BackgroundAtmosphere::CollisionProcessor();
  #elif _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_ == _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL__STOPPING_POWER_
  MolecularCollisions::BackgroundAtmosphere::StoppingPowerProcessor();
  #else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
  #endif //_PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_

  BackgroundAtmosphereCollisionTime=MPI_Wtime()-BackgroundAtmosphereCollisionTime;
#endif //_PIC_BACKGROUND_ATMOSPHERE_MODE_

  //particle photochemistry model
  #if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__PhotolyticReactions_;
  PhotoChemistryTime=MPI_Wtime();
  ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  PhotoChemistryTime=MPI_Wtime()-PhotoChemistryTime;
  #endif  //_PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_

  //perform user-define processing of the model particles
#if _PIC_USER_PARTICLE_PROCESSING__MODE_ == _PIC_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__UserParticleProcessing_;
  UserDefinedParticleProcessingTime=MPI_Wtime();
  PIC::UserParticleProcessing::Processing();
  UserDefinedParticleProcessingTime=MPI_Wtime()-UserDefinedParticleProcessingTime;
#endif

  //move existing particles
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__MoveParticles_;
  ParticleMovingTime=MPI_Wtime();
  PIC::Mover::MoveParticles();
  ParticleMovingTime=MPI_Wtime()-ParticleMovingTime;

  //check the consistence of the particles lists
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__CheckParticleList_;
  PIC::ParticleBuffer::CheckParticleList();
#endif

  IterationExecutionTime=MPI_Wtime()-StartTime;
  summIterationExecutionTime+=IterationExecutionTime;

  //syncrinie processors and exchnge particle data
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__ExchangeParticleData_;
  ParticleExchangeTime=MPI_Wtime();
  PIC::Parallel::ExchangeParticleData();
  ParticleExchangeTime=MPI_Wtime()-ParticleExchangeTime;

  //in case the dust model is turned on: re-distribute the dust particles in the velocity groups
  #if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if (_PIC_MODEL__DUST__ADJUST_VELOCITY_GROUP__MODE_ == _PIC_MODE_ON_) {
      ElectricallyChargedDust::GrainVelocityGroup::AdjustParticleVelocityGroup();
    }
  #endif //_PIC_MODEL__DUST__MODE_

  //incase OpenMP is used: rebalance the list of the available particles between OpenMP threads
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__RebalanceParticleList_;
  PIC::ParticleBuffer::Thread::RebalanceParticleList();
  #endif


  //update the total number of the sampled trajecotries
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__UpdateTrajectoryCounter_;
  PIC::ParticleTracker::UpdateTrajectoryCounter();
#endif

  //call user defined MPI procedure
#if _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_MODE_ == _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_MODE__ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__UserDefinedMPI_RoutineExecutionTime_;
  UserDefinedMPI_RoutineExecutionTime=MPI_Wtime();
  _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_();
  UserDefinedMPI_RoutineExecutionTime=MPI_Wtime()-UserDefinedMPI_RoutineExecutionTime;
#endif

  //update the glabal time counter if needed
#if _PIC_GLOBAL_TIME_COUNTER_MODE_ == _PIC_MODE_ON_
  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__CPLR_DATAFILE_MULTIFILE_;
  PIC::SimulationTime::Update();
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
  //update data
  if(PIC::CPLR::DATAFILE::MULTIFILE::IsTimeToUpdate()){
    //check if reached the last file
    if(PIC::CPLR::DATAFILE::MULTIFILE::ReachedLastFile){
      //end simulation if needed
      if(PIC::CPLR::DATAFILE::MULTIFILE::BreakAtLastFile)
	return _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_;
    }
    else //load data file
      PIC::CPLR::DATAFILE::MULTIFILE::UpdateDataFile();
  }
#if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
  // update field lines
  PIC::FieldLine::Update();
#endif//_PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
  
#endif//_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
#endif//_PIC_GLOBAL_TIME_COUNTER_MODE_ == _PIC_MODE_ON_

  ExitErrorCode=_PIC__EXIT_CODE__LAST_BLOCK__ExchangeStatisticData_;

  struct cExchangeStatisticData {
    double TotalInterationRunTime;
    double IterationExecutionTime;
    long int TotalParticlesNumber;
    double ParticleExchangeTime;
    double SamplingTime;
    double ParticleCollisionTime;
    double BackgroundAtmosphereCollisionTime;
    double InjectionBoundaryTime;
    double ParticleMovingTime;
    double PhotoChemistryTime;
    double Latency;
    long int recvParticleCounter,sendParticleCounter;
    long int nInjectedParticles;
    double UserDefinedMPI_RoutineExecutionTime;
    double UserDefinedParticleProcessingTime;
  };

  struct cStatExchangeControlParameters {
    int nExchangeStatisticsIterationNumberSteps;
    bool WallTimeExeedsLimit;
  };


  if (nInteractionsAfterRunStatisticExchange==nExchangeStatisticsIterationNumberSteps) { //collect and exchenge the statistical data of the run
    cExchangeStatisticData localRunStatisticData;
    int thread;
    cExchangeStatisticData *ExchangeBuffer=NULL;
    long int nInjectedParticleExchangeBuffer[PIC::nTotalThreads];
    double ParticleProductionRateExchangeBuffer[PIC::nTotalThreads];
    double ParticleMassProductionRateExchangeBuffer[PIC::nTotalThreads];

    localRunStatisticData.TotalInterationRunTime=MPI_Wtime()-StartTime;
    localRunStatisticData.IterationExecutionTime=IterationExecutionTime;
    localRunStatisticData.TotalParticlesNumber=PIC::ParticleBuffer::GetAllPartNum();
    localRunStatisticData.ParticleExchangeTime=ParticleExchangeTime;
    localRunStatisticData.SamplingTime=SamplingTime;
    localRunStatisticData.ParticleCollisionTime=ParticleCollisionTime;
    localRunStatisticData.BackgroundAtmosphereCollisionTime=BackgroundAtmosphereCollisionTime;
    localRunStatisticData.ParticleMovingTime=ParticleMovingTime;
    localRunStatisticData.PhotoChemistryTime=PhotoChemistryTime;
    localRunStatisticData.InjectionBoundaryTime=InjectionBoundaryTime;
    localRunStatisticData.Latency=PIC::Parallel::Latency;
    localRunStatisticData.recvParticleCounter=PIC::Parallel::recvParticleCounter;
    localRunStatisticData.sendParticleCounter=PIC::Parallel::sendParticleCounter;
    localRunStatisticData.nInjectedParticles=PIC::BC::nTotalInjectedParticles;
    localRunStatisticData.UserDefinedMPI_RoutineExecutionTime=UserDefinedMPI_RoutineExecutionTime;
    localRunStatisticData.UserDefinedParticleProcessingTime=UserDefinedParticleProcessingTime;

    PIC::BC::nTotalInjectedParticles=0;


    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      fprintf(PIC::DiagnospticMessageStream,"\n$PREFIX: (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,nTotalIterations,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);



      ExchangeBuffer=new cExchangeStatisticData[PIC::Mesh::mesh.nTotalThreads];
      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);


      //output the data
      long int nTotalModelParticles=0;
      double nTotalInjectedParticels=0.0;
      double MinExecutionTime=localRunStatisticData.IterationExecutionTime,MaxExecutionTime=localRunStatisticData.IterationExecutionTime,MaxLatency=0.0,MeanLatency=0.0;

      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Description:\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:1:\t Thread\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:2:\t Total Particle's number\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:3:\t Total Interation Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:4:\t Iteration Execution Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:5:\t Sampling Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:6:\t Injection Boundary Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:7:\t Particle Moving Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:8:\t Photo Chemistry Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:9:\t Particle Exchange Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:10:\t Latency\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:11:\t Send Particles\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:12:\t Recv Particles\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:13:\t nInjected Particls\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:14:\t User Defined MPI Routine - Execution Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:15:\t Particle Collision Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:16:\t Background Atmosphere Collision Time\n");
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:17:\t User Defined Particle Processing Time\n");


//      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:1\t 2\t 3\t\t 4\t\t 5\t\t 6\t\t 7\t\t 8\t\t 9\t\t 10\t 11\t 12\t 13\t\t 14\t\t 15\t 16\t 17\n");


      fprintf(PIC::DiagnospticMessageStream,"$PREFIX: ");
      for (int i=1;i<=17;i++) fprintf(PIC::DiagnospticMessageStream,"%12d ",i);
      fprintf(PIC::DiagnospticMessageStream,"\n");


      for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
        fprintf(PIC::DiagnospticMessageStream,"$PREFIX: %12d %12d %10e %10e %10e %10e %10e %10e %10e %10e %12d %12d %10e %10e %10e %10e %10e\n",
            thread,
            (int)ExchangeBuffer[thread].TotalParticlesNumber,
            ExchangeBuffer[thread].TotalInterationRunTime,
            ExchangeBuffer[thread].IterationExecutionTime,
            ExchangeBuffer[thread].SamplingTime,
            ExchangeBuffer[thread].InjectionBoundaryTime,
            ExchangeBuffer[thread].ParticleMovingTime,
            ExchangeBuffer[thread].PhotoChemistryTime,
            ExchangeBuffer[thread].ParticleExchangeTime,
            ExchangeBuffer[thread].Latency,

            (int)ExchangeBuffer[thread].sendParticleCounter,
            (int)ExchangeBuffer[thread].recvParticleCounter,
            ExchangeBuffer[thread].nInjectedParticles/double(((nExchangeStatisticsIterationNumberSteps!=0) ? nExchangeStatisticsIterationNumberSteps : 1)),
            ExchangeBuffer[thread].UserDefinedMPI_RoutineExecutionTime,
            ExchangeBuffer[thread].ParticleCollisionTime,
            ExchangeBuffer[thread].BackgroundAtmosphereCollisionTime,
            ExchangeBuffer[thread].UserDefinedParticleProcessingTime);

        nTotalModelParticles+=ExchangeBuffer[thread].TotalParticlesNumber;
        nTotalInjectedParticels+=ExchangeBuffer[thread].nInjectedParticles;
        MeanLatency+=ExchangeBuffer[thread].Latency;

        if (MinExecutionTime>ExchangeBuffer[thread].IterationExecutionTime) MinExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxExecutionTime<ExchangeBuffer[thread].IterationExecutionTime) MaxExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxLatency<ExchangeBuffer[thread].Latency) MaxLatency=ExchangeBuffer[thread].Latency;
      }

      MeanLatency/=PIC::Mesh::mesh.nTotalThreads;
      PIC::Parallel::CumulativeLatency+=MeanLatency*nInteractionsAfterRunStatisticExchange;
      if (nExchangeStatisticsIterationNumberSteps!=0) nTotalInjectedParticels/=nExchangeStatisticsIterationNumberSteps;

      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Total number of particles: %ld\n",nTotalModelParticles);

      //exchange statistics of the particle production
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Species dependent particle production:\nSpecie\tInjected Particles\tProductionRate\tMassProductionRate\n");

      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        double c=0.0;

        MPI_Gather(PIC::BC::nInjectedParticles+spec,1,MPI_LONG,nInjectedParticleExchangeBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Gather(PIC::BC::ParticleProductionRate+spec,1,MPI_DOUBLE,ParticleProductionRateExchangeBuffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Gather(PIC::BC::ParticleMassProductionRate+spec,1,MPI_DOUBLE,ParticleMassProductionRateExchangeBuffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

        for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
          nInjectedParticleExchangeBuffer[0]+=nInjectedParticleExchangeBuffer[thread];
          ParticleProductionRateExchangeBuffer[0]+=ParticleProductionRateExchangeBuffer[thread];
          ParticleMassProductionRateExchangeBuffer[0]+=ParticleMassProductionRateExchangeBuffer[thread];
        }

        if (nExchangeStatisticsIterationNumberSteps!=0) {
          c=double(nInjectedParticleExchangeBuffer[0])/nExchangeStatisticsIterationNumberSteps;
          ParticleProductionRateExchangeBuffer[0]/=nExchangeStatisticsIterationNumberSteps;
          ParticleMassProductionRateExchangeBuffer[0]/=nExchangeStatisticsIterationNumberSteps;
        }

        fprintf(PIC::DiagnospticMessageStream,"$PREFIX:%i\t%e\t%e\t%e\t(%s)\n",spec,c,ParticleProductionRateExchangeBuffer[0],ParticleMassProductionRateExchangeBuffer[0],PIC::MolecularData::GetChemSymbol(spec));

        PIC::BC::nInjectedParticles[spec]=0,PIC::BC::ParticleProductionRate[spec]=0.0,PIC::BC::ParticleMassProductionRate[spec]=0.0;
      }

      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Total number of injected particles: %e\n",nTotalInjectedParticels);
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Iteration Execution Time: min=%e, max=%e\n",MinExecutionTime,MaxExecutionTime);
      fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Latency: max=%e,mean=%e;CumulativeLatency=%e,Iterations after rebalabcing=%ld,Rebalancing Time=%e\n",MaxLatency,MeanLatency,PIC::Parallel::CumulativeLatency,PIC::Parallel::IterationNumberAfterRebalancing,PIC::Parallel::RebalancingTime);

      //flush the stream
      fflush(PIC::DiagnospticMessageStream);

      //check the elapsed walltime for the alarm
      if (PIC::Alarm::AlarmInitialized==true) {
        if (MPI_Wtime()-PIC::Alarm::StartTime>PIC::Alarm::RequestedExecutionWallTime) PIC::Alarm::WallTimeExeedsLimit=true;
        fprintf(PIC::DiagnospticMessageStream,"$PREFIX:Execution walltime: %e sec\n",MPI_Wtime()-PIC::Alarm::StartTime);
      }


      //determine the new value for the interation number between exchanging of the run stat data
      if (localRunStatisticData.TotalInterationRunTime>0.0) {
        nExchangeStatisticsIterationNumberSteps=(int)(_PIC_RUNTIME_STAT_OUTPUT__TIME_INTERVAL_/localRunStatisticData.TotalInterationRunTime);
        if (nExchangeStatisticsIterationNumberSteps<_PIC_RUNTIME_STAT_OUTPUT__MIN_ITERATION_NUMBER_) nExchangeStatisticsIterationNumberSteps=_PIC_RUNTIME_STAT_OUTPUT__MIN_ITERATION_NUMBER_;
        if (nExchangeStatisticsIterationNumberSteps>_PIC_RUNTIME_STAT_OUTPUT__MAX_ITERATION_NUMBER_) nExchangeStatisticsIterationNumberSteps=_PIC_RUNTIME_STAT_OUTPUT__MAX_ITERATION_NUMBER_;
      }
      else nExchangeStatisticsIterationNumberSteps=_PIC_RUNTIME_STAT_OUTPUT__MIN_ITERATION_NUMBER_;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      nExchangeStatisticsIterationNumberSteps=_PIC_RUNTIME_STAT_OUTPUT__MAX_ITERATION_NUMBER_;
#endif

      cStatExchangeControlParameters StatExchangeControlParameters;
      StatExchangeControlParameters.nExchangeStatisticsIterationNumberSteps=nExchangeStatisticsIterationNumberSteps;
      StatExchangeControlParameters.WallTimeExeedsLimit=PIC::Alarm::WallTimeExeedsLimit;

      MPI_Bcast(&StatExchangeControlParameters,sizeof(cStatExchangeControlParameters),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);

      //flush the diagnostric stream
      fflush(PIC::DiagnospticMessageStream);

      delete [] ExchangeBuffer;
    }
    else {
      cStatExchangeControlParameters StatExchangeControlParameters;

      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);

      //exchange statistics of the particle production
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        MPI_Gather(PIC::BC::nInjectedParticles+spec,1,MPI_LONG,nInjectedParticleExchangeBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Gather(PIC::BC::ParticleProductionRate+spec,1,MPI_DOUBLE,ParticleProductionRateExchangeBuffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Gather(PIC::BC::ParticleMassProductionRate+spec,1,MPI_DOUBLE,ParticleMassProductionRateExchangeBuffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

        PIC::BC::nInjectedParticles[spec]=0,PIC::BC::ParticleProductionRate[spec]=0.0,PIC::BC::ParticleMassProductionRate[spec]=0.0;
      }


      MPI_Bcast(&StatExchangeControlParameters,sizeof(cStatExchangeControlParameters),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);

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

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    EmergencyLoadRebalancingFlag=true;
#endif
#endif

    //check if the number of iterations between the load rebalancing exeeds the minimum number '_PIC_DYNAMIC_LOAD_BALANCING__MIN_ITERATION_BETWEEN_LOAD_REBALANCING_'
#ifdef _PIC_DYNAMIC_LOAD_BALANCING__MIN_ITERATION_BETWEEN_LOAD_REBALANCING_
    if (PIC::Parallel::IterationNumberAfterRebalancing<_PIC_DYNAMIC_LOAD_BALANCING__MIN_ITERATION_BETWEEN_LOAD_REBALANCING_) EmergencyLoadRebalancingFlag=false;
#endif

    MPI_Bcast(&EmergencyLoadRebalancingFlag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    if (EmergencyLoadRebalancingFlag==true) {
      if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"Load Rebalancing.....  begins\n");

      //correct the node's load balancing measure
      #if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
      if (summIterationExecutionTime>0.0)  {
        //normalize the load to the summ of the execution time
        double c,norm=0.0;
        int nLocalNode;

        for (nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) norm+=DomainBlockDecomposition::BlockTable[nLocalNode]->ParallelLoadMeasure;

        for (nLocalNode=0,c=summIterationExecutionTime/norm;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
          DomainBlockDecomposition::BlockTable[nLocalNode]->ParallelLoadMeasure*=c;
        }
      }

      summIterationExecutionTime=0.0;
      #endif //_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

      //start the rebalancing procedure
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Parallel::RebalancingTime=MPI_Wtime();

      PIC::Mesh::mesh.CreateNewParallelDistributionLists(_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
      PIC::Parallel::IterationNumberAfterRebalancing=0,PIC::Parallel::CumulativeLatency=0.0;
      PIC::DomainBlockDecomposition::UpdateBlockTable();

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Parallel::RebalancingTime=MPI_Wtime()-PIC::Parallel::RebalancingTime;
      if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"Load Rebalancing.....  done\n");
    }
#elif _PIC_EMERGENCY_LOAD_REBALANCING_MODE_ == _PIC_MODE_OFF_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    nInteractionsAfterRunStatisticExchange=0;

    //fflush the diagnostic output stream
    if (strcmp(PIC::DiagnospticMessageStreamName,"stdout")!=0) {
      fflush(PIC::DiagnospticMessageStream);
    }
  }

  return _PIC_TIMESTEP_RETURN_CODE__SUCCESS_;

}
//====================================================
//the general sampling procedure
void PIC::Sampling::Sampling() {
  int s,i,j,k,idim;
  long int LocalCellNumber,ptr,ptrNext;
  
  #if _PIC_SAMPLING_MODE_ == _PIC_MODE_ON_ //<-- begining of the particle sample section

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData,*ParticleDataNext;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  char *SamplingData;
  double *v,LocalParticleWeight,Speed2,v2;
  
  double *sampledVelocityOffset,*sampledVelocity2Offset;
  
  //the total number of the sampled particles to compare with the number of the partticles in the buffer
  long int nTotalSampledParticles=0;
  
  //temporary buffer for sampled data
  char  tempSamplingBuffer[PIC::Mesh::sampleSetDataLength]; 
  
  //temporaty buffer to store the copy of the particle
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];


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
//#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
  long int TreeNodeTotalParticleNumber;
//#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double TreeNodeProcessingTime;
//#endif

  //the table of cells' particles
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  

#if _PIC_SAMPLE_PARTICLE_DATA_MODE_ == _PIC_SAMPLE_PARTICLE_DATA_MODE__BETWEEN_ITERATIONS_
  //go through the 'local nodes'
//  while (node!=NULL) {

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) reduction(+:nTotalSampledParticles) default(none) \
  private (s,i,j,k,idim,LocalCellNumber,ptr,ptrNext,ParticleData, ParticleDataNext,cell,block,SamplingData,v,LocalParticleWeight,Speed2,v2,sampledVelocityOffset,node, \
    sampledVelocity2Offset,tempSamplingBuffer,tempParticleData,FirstCellParticleTable,TreeNodeProcessingTime,TreeNodeTotalParticleNumber)  \
  shared (PIC::Mesh::mesh,sampledParticleWeghtRelativeOffset,sampledParticleNumberRelativeOffset,sampledParticleNumberDensityRelativeOffset,sampledParticleVelocityRelativeOffset, \
     DomainBlockDecomposition::nLocalBlocks,DomainBlockDecomposition::BlockTable,sampledParticleVelocity2RelativeOffset,sampleSetDataLength,sampledParticleSpeedRelativeOffset,collectingCellSampleDataPointerOffset,ParticleDataLength, \
     PIC::IDF::_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_, PIC::IDF::_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_,PIC::IDF::_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_,  \
     PIC::Mesh::DatumParticleParallelVelocity2,PIC::Mesh::DatumParticleParallelVelocity,PIC::Mesh::DatumParticleSpeed,PIC::Mesh::DatumParticleVelocity2,PIC::Mesh::DatumParticleVelocity, \
     PIC::Mesh::DatumNumberDensity,PIC::Mesh::DatumParticleNumber,PIC::Mesh::DatumParticleWeight, PIC::Sampling::constNormalDirection__SampleParallelTangentialTemperature, \
     PIC::IDF::nTotalVibtationalModes, cout)

#endif //_COMPILATION_MODE_
  for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=DomainBlockDecomposition::BlockTable[nLocalNode];

//  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;
    
    memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
    
    
    
    
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    TreeNodeTotalParticleNumber=0;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    TreeNodeProcessingTime=MPI_Wtime();
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_OFF_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the otion is not recognized");
#endif


    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
	  
	  ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

	  if (ptr!=-1) {
	    
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            cell=block->GetCenterNode(LocalCellNumber);
	    
            SamplingData=
	      cell->GetAssociatedDataBufferPointer() + 
	      collectingCellSampleDataPointerOffset;

            memcpy((void*)tempSamplingBuffer,
		   (void*)SamplingData,sampleSetDataLength);
	    
	    
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
            for (s=0;s<PIC::nTotalSpecies;s++) PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingData+sampledParticleWeghtRelativeOffset)),1,__LINE__,__FILE__);
#endif
#endif


            ptrNext=ptr;
            ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

            //===================    DEBUG ==============================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
            if (cell->Measure<=0.0) {
              cout << "$PREFIX:" << __FILE__<< __LINE__ << std::endl;
              exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
            }

#endif
            //===================   END DEBUG ==============================

            while (ptrNext!=-1) {
              ptr=ptrNext;
              ParticleData=ParticleDataNext;
              nTotalSampledParticles++;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
              TreeNodeTotalParticleNumber++;
#endif


              memcpy((void*)tempParticleData,
		     (void*)ParticleData,ParticleDataLength);

              ptrNext=PIC::ParticleBuffer::GetNext((PIC::ParticleBuffer::byte*)tempParticleData);

              //================ Prefetch particle data
              if (ptrNext!=-1) {
                ParticleDataNext=PIC::ParticleBuffer::GetParticleDataPointer(ptrNext);

#if _PIC_MEMORY_PREFETCH_MODE_ == _PIC_MEMORY_PREFETCH_MODE__ON_
                int iPrefetch,iPrefetchMax=1+(int)(ParticleDataLength/_PIC_MEMORY_PREFETCH__CHACHE_LINE_);

                for (iPrefetch=0;iPrefetch<iPrefetchMax;iPrefetch++) {
                  _mm_prefetch(iPrefetch*_PIC_MEMORY_PREFETCH__CHACHE_LINE_+(char*)ptrNext,_MM_HINT_T1);
                }
#endif
              }


              //================ End prefetch particle data

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING__PARTICLE_COORDINATES_ == _PIC_DEBUGGER_MODE_ON_
              //check position of the particle:
              //the particle must be within the computational domain and outside of the internal boundarues
              list<cInternalBoundaryConditionsDescriptor>::iterator InternalBoundaryListIterator;

              for (InternalBoundaryListIterator=PIC::Mesh::mesh.InternalBoundaryList.begin();InternalBoundaryListIterator!=PIC::Mesh::mesh.InternalBoundaryList.end();InternalBoundaryListIterator++) {
                cInternalBoundaryConditionsDescriptor InternalBoundaryDescriptor;
                double x[3];
                double *x0Sphere,radiusSphere;

                PIC::ParticleBuffer::GetX(x,ptr);
                InternalBoundaryDescriptor=*InternalBoundaryListIterator;

                switch (InternalBoundaryDescriptor.BondaryType) {
                case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
                  #if DIM == 3
                  ((cInternalSphericalData*)(InternalBoundaryDescriptor.BoundaryElement))->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

                  if (pow(x[0]-x0Sphere[0],2)+pow(x[1]-x0Sphere[1],2)+pow(x[2]-x0Sphere[2],2)<pow(radiusSphere-PIC::Mesh::mesh.EPS,2)) {
                    cout << "$PREFIX:" << __FILE__ << "@" << __LINE__ << "Sphere: x0=" << x0Sphere[0] << ", " << x0Sphere[1] << ", " << x0Sphere[2] << ", R=" << radiusSphere << endl;
                    cout << "$PREFIX:" << __FILE__ << "@" << __LINE__ << "Particle Position: x=" << x[0] << ", " << x[1] << ", " << x[2] << endl;
                    cout << "$PREFIX:" << __FILE__ << "@" << __LINE__ << "Particle Distance from the center of the sphere: " << sqrt(pow(x[0]-x0Sphere[0],2)+pow(x[1]-x0Sphere[1],2)+pow(x[2]-x0Sphere[2],2)) << endl;

                    exit(__LINE__,__FILE__,"Error: particle inside spherical body");
                  }
                  #endif

		  
                  break;
                default:
                  exit(__LINE__,__FILE__,"Error: not implemented");
                }

                for (idim=0;idim<DIM;idim++) if ((x[idim]<node->xmin[idim])||(x[idim]>node->xmax[idim])) exit(__LINE__,__FILE__,"Error: particle is outside of the block");
              }
#endif
#endif

	    
              Speed2=0.0;

              s=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)tempParticleData);
              v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)tempParticleData);

              LocalParticleWeight=block->GetLocalParticleWeight(s);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)tempParticleData);

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
              PIC::Debugger::CatchOutLimitValue(v,DIM,__LINE__,__FILE__);
              PIC::Debugger::CatchOutLimitValue(LocalParticleWeight,__LINE__,__FILE__);
#endif
#endif


	      // sample data
	      cell->SampleDatum(PIC::Mesh::DatumParticleWeight,
				LocalParticleWeight, s);
	      cell->SampleDatum(PIC::Mesh::DatumParticleNumber, 1.0, s);
	      cell->SampleDatum(PIC::Mesh::DatumNumberDensity, 
				LocalParticleWeight/cell->Measure, s);

	      double miscv2[3];
	      for (idim=0;idim<3;idim++) {
		v2=v[idim]*v[idim];
		Speed2+=v2;
		miscv2[idim]=v2;
	      }
	      cell->SampleDatum(PIC::Mesh::DatumParticleVelocity,
				v, s, LocalParticleWeight);
	      cell->SampleDatum(PIC::Mesh::DatumParticleVelocity2,
				miscv2, s, LocalParticleWeight);
	      cell->SampleDatum(PIC::Mesh::DatumParticleSpeed,
				sqrt(Speed2), s, LocalParticleWeight);

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
              PIC::Debugger::CatchOutLimitValue((s+(double*)(tempSamplingBuffer+sampledParticleWeghtRelativeOffset)),1,__LINE__,__FILE__); 
              PIC::Debugger::CatchOutLimitValue(sampledVelocityOffset,DIM,__LINE__,__FILE__);
              PIC::Debugger::CatchOutLimitValue(sampledVelocity2Offset,DIM,__LINE__,__FILE__);
              PIC::Debugger::CatchOutLimitValue((s+(double*)(tempSamplingBuffer+sampledParticleSpeedRelativeOffset)),1,__LINE__,__FILE__);
#endif
#endif


              //sample the data for calculation the normal and tangential kinetic temepratures
            #if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
              //do nothing
            #else
              //calcualte the direction of the normal
              double l[3];

              #if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__FUNSTION_CALCULATED_NORMAL_DIRECTION_
              exit(__LINE__,__FILE__,"error: not developed: need a function for calculation of the local direction for the normal for calcualtion of the parallel and tengential temperatures");
              #elif _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__CONSTANT_DIRECTION_ORIGIN_
              double c,*x;

              x=PIC::ParticleBuffer::GetX((PIC::ParticleBuffer::byte*)tempParticleData);

              l[0]=x[0]-constNormalDirection__SampleParallelTangentialTemperature[0];
              l[1]=x[1]-constNormalDirection__SampleParallelTangentialTemperature[1];
              l[2]=x[2]-constNormalDirection__SampleParallelTangentialTemperature[2];

              c=sqrt((l[0]*l[0])+(l[1]*l[1])+(l[2]*l[2]));

              l[0]/=c,l[1]/=c,l[2]/=c;
              #else
              exit(__LINE__,__FILE__,"Error: the option is not defined");
              #endif

              double vParallel=(l[0]*v[0])+(l[1]*v[1])+(l[2]*v[2]);
	      
	      cell->SampleDatum(PIC::Mesh::DatumParticleParallelVelocity,
				vParallel, s, LocalParticleWeight);
	      cell->SampleDatum(PIC::Mesh::DatumParticleParallelVelocity2,
				vParallel*vParallel, s, LocalParticleWeight);
#endif

              //sample data for the internal degrees of freedom model

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_OFF_
              //do mothing
#elif _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
              //save the total particle weight used for further interpolation
              *(s+(double*)(SamplingData+IDF::_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_))+=LocalParticleWeight;

              //save rotational energy
              *(s+(double*)(SamplingData+IDF::_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_))+=IDF::GetRotE((PIC::ParticleBuffer::byte*)tempParticleData)*LocalParticleWeight;

              //save vibrational evergy
              for (int nmode=0;nmode<IDF::nTotalVibtationalModes[s];nmode++) {
                *(s+(double*)(SamplingData+IDF::_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s]))+=IDF::GetVibE(nmode,(PIC::ParticleBuffer::byte*)tempParticleData)*LocalParticleWeight;
              }

              //save the population of the first two vibrational levels for qLB
              #if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL_ == _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODEL__QLB_
              {
                int VibrationLevel=*((int*)(tempParticleData+IDF::qLB::nVibLevel));

                switch (VibrationLevel) {
                case 0:
                  *(s+(double*)(SamplingData+IDF::qLB::_VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_))+=LocalParticleWeight;
                  break;
                case 1:
                  *(s+(double*)(SamplingData+IDF::qLB::_VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_))+=LocalParticleWeight;
                  break;
                }
              }
              #endif


#else
              exit(__LINE__,__FILE__,"the option is not defined");
#endif

              //call sampling procedures of indivudual models
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
              ElectricallyChargedDust::Sampling::SampleParticleData(tempParticleData,LocalParticleWeight, SamplingData, s);
#endif

	      //call sampling procedures of indivudual models

#if _PIC_MOVER_INTEGRATOR_MODE_ == _PIC_MOVER_INTEGRATOR_MODE__GUIDING_CENTER_
	      PIC::Mover::GuidingCenter::Sampling::SampleParticleData(tempParticleData,LocalParticleWeight, SamplingData, s);//tempSamplingBuffer, s);
#endif //_PIC_MOVER_INTEGRATOR_MODE_ 

#if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
	      PIC::FieldLine::Sampling(ptr,LocalParticleWeight, SamplingData);
#endif//_PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_ 

              //call user defined particle sampling procedure
#ifdef _PIC_USER_DEFING_PARTICLE_SAMPLING_
              _PIC_USER_DEFING_PARTICLE_SAMPLING_(tempParticleData,LocalParticleWeight,SamplingData,s);
#endif

#ifdef _PIC_USER_DEFING_PARTICLE_SAMPLING__NODE_ //call a user-defind particle sampling procedure with passing the node information
              _PIC_USER_DEFING_PARTICLE_SAMPLING__NODE_(tempParticleData,LocalParticleWeight,SamplingData,s,node);
#endif

            }

	    //            memcpy((void*)SamplingData,(void*)tempSamplingBuffer,/*PIC::Mesh::*/sampleSetDataLength);

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


//    node=node->nextNodeThisThread;
  }
#elif _PIC_SAMPLE_PARTICLE_DATA_MODE_ == _PIC_SAMPLE_PARTICLE_DATA_MODE__DURING_PARTICLE_MOTION_
  //do nothing
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   CatchOutLimitSampledValue();
#endif
#endif

   //check if the number of sampled particles coinsides with the number of particles in the buffer
   if (nTotalSampledParticles!=ParticleBuffer::GetAllPartNum()) exit(__LINE__,__FILE__,"The number of the sampled particles is different from that in the particle buffer");

  //sample user defined data
  if (PIC::IndividualModelSampling::SamplingProcedure.size()!=0) {
    int nfunc,nfuncTotal=PIC::IndividualModelSampling::SamplingProcedure.size();

    for (nfunc=0;nfunc<nfuncTotal;nfunc++) PIC::IndividualModelSampling::SamplingProcedure[nfunc]();

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   CatchOutLimitSampledValue();
#endif
#endif
  }

  //sample local data sets of the user defined functions
  for (int nfunc=0;nfunc<PIC::Sampling::ExternalSamplingLocalVariables::SamplingRoutinesRegistrationCounter;nfunc++) PIC::Sampling::ExternalSamplingLocalVariables::SamplingProcessor[nfunc]();

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   CatchOutLimitSampledValue();
#endif
#endif


  //sample the distribution functions
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
  if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::SampleDistributionFnction();
  if (PIC::EnergyDistributionSampleRelativistic::SamplingInitializedFlag==true) PIC::EnergyDistributionSampleRelativistic::SampleDistributionFnction();
  if (PIC::ParticleFluxDistributionSample::SamplingInitializedFlag==true) PIC::ParticleFluxDistributionSample::SampleDistributionFnction();
  if (_PIC_COUPLER_MODE_!=_PIC_COUPLER_MODE__OFF_) if (PIC::PitchAngleDistributionSample::SamplingInitializedFlag==true) PIC::PitchAngleDistributionSample::SampleDistributionFnction();
#endif

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
   CatchOutLimitSampledValue();
#endif
#endif

  //Sample size distribution parameters of dust grains
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SampleDistributionFnction();
  ElectricallyChargedDust::Sampling::FluxMap::Sampling();
#endif

  //END OF THE PARTICLE SAMPLING SECTION
  #endif //<-- end of the particle sampling section

  //Increment the sample length
  CollectingSampleCounter++;

  //Output the data flow file if needed
  if ((CollectingSampleCounter==RequiredSampleLength)||(PIC::Alarm::WallTimeExeedsLimit==true)) {
    //exchnge the sampling data
    PIC::Mesh::mesh.ParallelBlockDataExchange();

    //check different sampling modes
    #if _PIC_SAMPLING_MODE_ == _PIC_MODE_ON_
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
/*      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();
//      PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere* Sphere;

      cInternalSphericalData *Sphere;

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
//        Sphere=(PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere*)PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->GetSurfaceDataPointer();

        Sphere=PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface);

//        PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
      }*/


  #else
      exit(__LINE__,__FILE__,"Error: unknown option");
  #endif


    }
    else if (SamplingMode==_ACCUMULATE_SAMPLING_MODE_) {
      exit(__LINE__,__FILE__,"Error: the sampling mode '_ACCUMULATE_SAMPLING_MODE_' is not implemented");
    }
    else exit(__LINE__,__FILE__,"Error: the sampling mode is not defined");

    #else // <-- #if _PIC_SAMPLING_MODE_ == _PIC_MODE_ON_
    LastSampleLength=CollectingSampleCounter;
    CollectingSampleCounter=0;
    #endif // <-- #if _PIC_SAMPLING_MODE_ == _PIC_MODE_ON_

    //print the error messages accumulated by the Mover
    if (_PIC_MOVER__UNKNOWN_ERROR_IN_PARTICLE_MOTION__STOP_EXECUTION_ == _PIC_MODE_OFF_) {
      Mover::Sampling::Errors::PrintData();
    }

    //print output file
    char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

    if (LastSampleLength>=minIterationNumberForDataOutput) {

/*----------------------------------  BEGIN OUTPUT OF THE DATA FILES SECTION --------------------------------*/

      #if _CUT_CELL__TRIANGULAR_FACE__USER_DATA__MODE_ == _ON_AMR_MESH_
      //output the data sampled on the triangulated cut-cells
      sprintf(fname,"%s/amps.cut-cell.surface-data.out=%ld.dat",OutputDataFileDirectory,DataOutputFileNumber);
      CutCell::PrintSurfaceData(fname);
      #endif

      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      //print sampled particle trajectories
      sprintf(fname,"%s/amps.TrajectoryTracking.out=%ld",OutputDataFileDirectory,DataOutputFileNumber);
      PIC::ParticleTracker::OutputTrajectory(fname);
      #endif

      #if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
      //print sampled data along field lines
      sprintf(fname,
              "%s/amps.FieldLines.out=%ld.dat",
              OutputDataFileDirectory,DataOutputFileNumber);
      PIC::FieldLine::Output(fname, false);
      #endif

      //print the macroscopic parameters of the flow
      for (s=0;s<PIC::nTotalSpecies;s++) if (SaveOutputDataFile[s]==true) {
        PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
        sprintf(fname,"%s/pic.%s.s=%i.out=%ld.dat",OutputDataFileDirectory,ChemSymbol,s,DataOutputFileNumber);

        if (PIC::Mesh::mesh.ThisThread==0) {
          fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........",fname);
          fflush(stdout);
        }

        if (DataOutputFileNumber>=FirstPrintedOutputFile) {
          if (_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE_==_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE__TECPLOT_ASCII_) {
            PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);
          }
        }

        if (PIC::Mesh::mesh.ThisThread==0) {
          fprintf(PIC::DiagnospticMessageStream,"done.\n");
          fflush(stdout);
        }

      //print the sampled distribution function into a file
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
        if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) {
          sprintf(fname,"%s/pic.distribution.%s.s=%i.out=%ld",OutputDataFileDirectory,ChemSymbol,s,DataOutputFileNumber);
          PIC::DistributionFunctionSample::printDistributionFunction(fname,s);
        }

        if (PIC::EnergyDistributionSampleRelativistic::SamplingInitializedFlag==true) {
          sprintf(fname,"%s/pic.energy-distribution.%s.s=%i.out=%ld",OutputDataFileDirectory,ChemSymbol,s,DataOutputFileNumber);
          PIC::EnergyDistributionSampleRelativistic::printDistributionFunction(fname,s);
        }

        if (PIC::ParticleFluxDistributionSample::SamplingInitializedFlag==true) {
          sprintf(fname,"%s/pic.flux.%s.s=%i.out=%ld.dat",OutputDataFileDirectory,ChemSymbol,s,DataOutputFileNumber);
          PIC::ParticleFluxDistributionSample::printMacroscopicParameters(fname,s);
        }

	if (PIC::PitchAngleDistributionSample::SamplingInitializedFlag==true) {
          sprintf(fname,"%s/pic.pitch_angle.%s.s=%i.out=%ld.dat",OutputDataFileDirectory,ChemSymbol,s,DataOutputFileNumber);
          PIC::PitchAngleDistributionSample::printDistributionFunction(fname,s);
        }

#endif
      }

      //save the sampling data restart file in case when the macroscopic data are downloaded from remote host for post-processing
      if (_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE_==_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE__SAMPLING_DATA_RESTART_FILE_) {
        sprintf(fname,"%s/pic.SamplingDataRestart.out=%ld.dat",OutputDataFileDirectory,DataOutputFileNumber);
        PIC::Restart::SamplingData::Save(fname);
      }

      //print the sampled local data sets of the user defined functions
      for (int nfunc=0;nfunc<PIC::Sampling::ExternalSamplingLocalVariables::SamplingRoutinesRegistrationCounter;nfunc++) PIC::Sampling::ExternalSamplingLocalVariables::PrintOutputFile[nfunc](DataOutputFileNumber);

      //print the sampled total production rate due to volume injection
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
      if (PIC::VolumeParticleInjection::SourceRate!=NULL) if (LastSampleLength>=minIterationNumberForDataOutput) {
         double buffer[PIC::nTotalSpecies*PIC::nTotalThreads];

         MPI_Gather(PIC::VolumeParticleInjection::SourceRate,PIC::nTotalSpecies,MPI_DOUBLE,buffer,PIC::nTotalSpecies,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

         if (PIC::ThisThread==0) {
           for (int thread=1;thread<PIC::nTotalThreads;thread++) for (s=0;s<PIC::nTotalSpecies;s++) buffer[s]+=buffer[thread*PIC::nTotalSpecies+s];
           fprintf(PIC::DiagnospticMessageStream,"Total sources rate by volume production injection models: \n Species \t Source rate [s^{-1}]\n");

           for (s=0;s<PIC::nTotalSpecies;s++) {
             PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
             fprintf(PIC::DiagnospticMessageStream,"%s (s=%i):\t %e\n",ChemSymbol,s,buffer[s]/LastSampleLength);
           }

           fprintf(PIC::DiagnospticMessageStream,"\n");
         }

         if (SamplingMode==_RESTART_SAMPLING_MODE_) for (s=0;s<PIC::nTotalSpecies;s++) PIC::VolumeParticleInjection::SourceRate[s]=0.0;
      }
#endif

      //print into the file the sampled data of the internal surfaces installed into the mesh
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
#elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();

      for (s=0;s<PIC::nTotalSpecies;s++) for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
        PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
        sprintf(fname,"%s/pic.Sphere=%ld.%s.s=%i.out=%ld.dat",OutputDataFileDirectory,iSphericalSurface,ChemSymbol,s,DataOutputFileNumber);

        if (PIC::Mesh::mesh.ThisThread==0) {
          fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........",fname);
          fflush(stdout);
        }

        PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->PrintSurfaceData(fname,s);

        if (PIC::Mesh::mesh.ThisThread==0) {
          fprintf(PIC::DiagnospticMessageStream,"done.\n");
          fflush(stdout);
        }
      }

      //flush sampled surface data
      if (SamplingMode==_RESTART_SAMPLING_MODE_) {
        for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
          cInternalSphericalData *Sphere=PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface);
          PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
        }

      }
      else exit(__LINE__,__FILE__,"Error: the surface sampling is implemented only for the case of SamplingMode==_RESTART_SAMPLING_MODE_");


#else
      exit(__LINE__,__FILE__,"Error: unknown option");
#endif

      //clean user defined sampling buffers if needed
#if _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA__MODE_ == _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA__MODE__ON_
     _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA_();
#endif
/*----------------------------------  END OUTPUT OF THE DATA FILES SECTION --------------------------------*/
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

      fprintf(PIC::DiagnospticMessageStream,"Model run statistics:\n");
      fprintf(PIC::DiagnospticMessageStream,"The total number of model particles: %ld\n",nTotalSimulatedParticles);

      pipe.closeRecvAll();
    }



#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
    if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::flushSamplingBuffers();
    if (PIC::EnergyDistributionSampleRelativistic::SamplingInitializedFlag==true) PIC::EnergyDistributionSampleRelativistic::flushSamplingBuffers();
    if (PIC::ParticleFluxDistributionSample::SamplingInitializedFlag==true) PIC::ParticleFluxDistributionSample::flushSamplingBuffers();
#endif

    //Sample size distribution parameters of dust grains
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::printDistributionFunction(DataOutputFileNumber);
    ElectricallyChargedDust::Sampling::FluxMap::PrintSurfaceData(DataOutputFileNumber);
#endif


    //Finish the code execution if the walltime exeeds the limit
    if (PIC::Alarm::WallTimeExeedsLimit==true) {
      //save the restart file
      if (_PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ == _PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ON_) {
        char fname[_MAX_STRING_LENGTH_PIC_];

        if (Restart::ParticleDataRestartFileOverwriteMode==true) sprintf(fname,"%s",Restart::saveParticleDataRestartFileName);
        else sprintf(fname,"%s.Final",Restart::saveParticleDataRestartFileName);

        Restart::SaveParticleData(fname);
        MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      }

      //finish exdcution
      PIC::Alarm::FinishExecution();
    }


    //increment the output file number
    DataOutputFileNumber++;

    //redistribute the processor load and check the mesh afterward
#if _PIC_SAMPLING_BREAK_LOAD_REBALANCING_MODE_ == _PIC_MODE_ON_
    if (PIC::Parallel::IterationNumberAfterRebalancing!=0) {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Parallel::RebalancingTime=MPI_Wtime();

      PIC::Mesh::mesh.CreateNewParallelDistributionLists(_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
      PIC::Parallel::IterationNumberAfterRebalancing=0,PIC::Parallel::CumulativeLatency=0.0;
      PIC::DomainBlockDecomposition::UpdateBlockTable();

      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
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
        fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        fprintf(PIC::DiagnospticMessageStream,"done.\n");
        fflush(stdout);
      }
    }

    PIC::Mesh::switchSamplingBuffers();

    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld-redistributed-load-TempSample.dat",ChemSymbol,s,DataOutputFileNumber-1);

      if (PIC::Mesh::mesh.ThisThread==0) {
        fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        fprintf(PIC::DiagnospticMessageStream,"done.\n");
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
  cout << "$PREFIX:Signal is intersepted: thread=" << PIC::Mesh::mesh.ThisThread << endl;

  switch(sig) {
  case SIGFPE :
    cout << "$PREFIX:Signal=SIGFPE" << endl;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown signal");
  }

  exit(__LINE__,__FILE__,"Error: exit in the signal handler");
}

//====================================================
void PIC::InitMPI() {

  //check is MPI is initialized
  int initialized;

  MPI_Initialized(&initialized);

  if (!initialized) {
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
    MPI_Init(NULL, NULL);

    #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    int provided;
    MPI_Init_thread(NULL,NULL,MPI_THREAD_FUNNELED,&provided);

    #else
    #error Unknown option
    #endif //_COMPILATION_MODE_

    MPI_GLOBAL_COMMUNICATOR=MPI_COMM_WORLD;
  }

  //init MPI variables
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

  ::ThisThread=ThisThread;
  ::TotalThreadsNumber=nTotalThreads;

  //determine the total number of the OpenMP threads
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int nThreadsOpenMP,i;

  #pragma omp parallel shared(nThreadsOpenMP)
  {
    #pragma omp single
    {
      nTotalThreadsOpenMP=omp_get_num_threads();
    }
  }
#endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

  if (ThisThread==0) {
    printf("$PREFIX: The total number of the MPI processes=%i\n",nTotalThreads);
    printf("$PREFIX: The total number of the OpenMP threads per each MPI process=%i\n",nTotalThreadsOpenMP);
  }
}

//init the particle solver
void PIC::Init_BeforeParser() {

  //initiate MPI
  InitMPI();

  //intialize the interpolation module
  InterpolationRoutines::Init();

/*
  //init the particle buffer
  //the default value of 'RequestedParticleBufferLength' is -1, which can be chganged by ampsConfig.pl
  //if 'RequestedParticleBufferLength' is -1 than the particle buffer is not initialized
  const int RequestedParticleBufferLength=-1;

  if (RequestedParticleBufferLength!=-1) PIC::ParticleBuffer::Init(RequestedParticleBufferLength);
*/


  //set up the DiagnospticMessageStream
  if (strcmp(PIC::DiagnospticMessageStreamName,"stdout")!=0) {
    char cmd[_MAX_STRING_LENGTH_PIC_];

    if (PIC::ThisThread==0) {

      /*
      struct stat s;
      int err;

      err=stat(PIC::DiagnospticMessageStreamName,&s);

      if (err==-1) {
        //the directory does not exists -> need to create it
        std::string path=std::string(PIC::DiagnospticMessageStreamName);
        std::string t,CreatedDirectories;
        size_t pos=0;

        while ((pos=path.find("/"))!=std::string::npos) {
          t=path.substr(0,pos);
          path.erase(0,1+pos);

          CreatedDirectories+=t;
          CreatedDirectories+="/";

          err=stat(CreatedDirectories.c_str(),&s);

          if (err==-1) {
            if (mkdir(CreatedDirectories.c_str(),0777)==-1) {
              printf("$PREFIX:Cannot create the output directory\n");
              exit(__LINE__,__FILE__);
            }
          }
        }

        if (mkdir(PIC::DiagnospticMessageStreamName,0777)==-1) {
          printf("$PREFIX:Cannot create the output directory\n");
          exit(__LINE__,__FILE__);
        }
      }*/

      //remove the content of the output directory
      sprintf(cmd,"mkdir -p %s",PIC::DiagnospticMessageStreamName);
      system(cmd);

      sprintf(cmd,"rm -rf %s/*",PIC::DiagnospticMessageStreamName);
      system(cmd);
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    sprintf(PIC::DiagnospticMessageStreamName,"%s/thread=%i.log",PIC::DiagnospticMessageStreamName,PIC::ThisThread);
    PIC::DiagnospticMessageStream=fopen(PIC::DiagnospticMessageStreamName,"w");
    PIC::Mesh::mesh.DiagnospticMessageStream=PIC::DiagnospticMessageStream;
  }

  //create the output directory if needed
  if (strcmp(PIC::OutputDataFileDirectory,".")!=0) {
    char cmd[_MAX_STRING_LENGTH_PIC_];

    if (PIC::ThisThread==0) {
/*
      struct stat s;
      int err;

      err=stat(PIC::OutputDataFileDirectory,&s);

      if (err==-1) {
        //the directory does not exists -> need to create it
        std::string path=std::string(PIC::OutputDataFileDirectory);
        std::string t,CreatedDirectories;
        size_t pos=0;

        while ((pos=path.find("/"))!=std::string::npos) {
          t=path.substr(0,pos);
          path.erase(0,1+pos);

          CreatedDirectories+=t;
          CreatedDirectories+="/";

          err=stat(CreatedDirectories.c_str(),&s);

          if (err==-1) {
            if (mkdir(CreatedDirectories.c_str(),0777)==-1) {
              printf("$PREFIX:Cannot create the output directory\n");
              exit(__LINE__,__FILE__);
            }
          }
        }


        if (mkdir(PIC::OutputDataFileDirectory,0777)==-1) {
          printf("$PREFIX:Cannot create the output directory\n");
          exit(__LINE__,__FILE__);
        }
      }
      */

      //remove the content of the output directory
      sprintf(cmd,"mkdir -p %s",PIC::OutputDataFileDirectory);
      system(cmd);

      sprintf(cmd,"rm -rf %s/*",PIC::OutputDataFileDirectory);
      system(cmd);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if (ThisThread==0) {
    time_t TimeValue=time(NULL);
    tm *ct=localtime(&TimeValue);

    fprintf(PIC::DiagnospticMessageStream,"\n$PREFIX: (%i/%i %i:%i:%i), Initialization of the code\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
    fprintf(PIC::DiagnospticMessageStream,"$PREFIX: Simulation Target - %s \n",_MACRO_STR_VALUE_(_TARGET_));
  }




  /*-----------------------------------  The assembly code is not compiled with gcc 4.6: BEGIN  ---------------------------------

  //check the type of the CPU and the size of the cash line
  int CacheLineSize = -1;

  asm ( "mov $5, %%eax\n\t"   // EAX=80000005h: L1 Cache and TLB Identifiers
        "cpuid\n\t"
        "mov %%eax, %0"       // eax into CacheLineSize
        : "=r"(CacheLineSize)   // output
        :                     // no input
        : "%eax"              // clobbered register
       );

  if (ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"Cache line size if %i bytes\n",CacheLineSize);



  char VendorSign[13];
//  char *VendorSign_dword1=VendorSign+4,*VendorSign_dword2=VendorSign+8;
  unsigned int w0,w1,w2;

  __asm {
    xor eax,eax
    cpuid

    mov [w0],ebx
    mov [w1],edx
    mov [w2],ecx
  }

  memcpy(VendorSign,&w0,4);
  memcpy(VendorSign+4,&w1,4);
  memcpy(VendorSign+8,&w2,4);

  VendorSign[12]='\0';

  if (strcmp(VendorSign,"GenuineIntel")==0) {
    //intell processor
    //check the size of the cache line
    //function 80000006h return the size of the cache line in register ecx [bits 7:0]; ecx [bits 31:16] L2 cache size

    __asm {
      mov eax,0x80000006
      cpuid

      mov eax,ecx
      and eax,0xff
      mov [w0],eax

      mov eax,ecx
      and eax,0xff00
      mov [w1],eax
    }

    if (ThisThread==0) {
      fprintf(PIC::DiagnospticMessageStream,"CPU Manifacturer: INTEL\nCache line size is %i bytes, L2 cache size is %i KB\n",(int)w0,(int)w1);
    }

  }
  else if (strcmp(VendorSign,"AuthenticAMD")==0) {
    //AMD processor
    exit(__LINE__,__FILE__,"Error: not implemented");
  }
  else {
    cout << "Unknown type of CPU: the vendor string is \"" << VendorSign <<"\"" << endl;
    exit(__LINE__,__FILE__,"Error: unknown processor");
  }

  -----------------------------------  The assembly code is not compiled with gcc 4.6: END  ---------------------------------*/

  //init sections of the particle solver

  //set up the signal handler
  signal(SIGFPE,SignalHandler);


  //init coupler 
  if ((_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__DATAFILE_)||(_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__T96_)) {
    PIC::CPLR::DATAFILE::Init();
  }


  //init the background atmosphere model
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
  PIC::MolecularCollisions::BackgroundAtmosphere::Init_BeforeParser();
#endif

  //init the particle collision procedure
#if _PIC__PARTICLE_COLLISION_MODEL__MODE_ == _PIC_MODE_ON_
  PIC::MolecularCollisions::ParticleCollisionModel::Init();
#endif

  //init the model of internal degrees of freedom
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::Init();
#endif

  //Init the photolytic reaction model
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  ::PhotolyticReactions::Init();
#endif

  //Interpolation routines
  if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    if (_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_) {
      //in case AMPS' second order interpolation routine is used:
      //the number of cells in a block MUST be even
      if ((_BLOCK_CELLS_X_%2!=0)||(_BLOCK_CELLS_Y_%2!=0)||(_BLOCK_CELLS_Z_%2!=0)) exit(__LINE__,__FILE__,"Error: in case the second order interpolation routine is used the number of cells in a block MUST be even");
    }
  }

}

void PIC::Init_AfterParser() {
  int i,j,k;
  long int LocalCellNumber;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;

  //Interpolation routines
  if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    //set the interpolation retine for constructing of the stencil when output the model data file
    if (_PIC_OUTPUT__CELL_CORNER_INTERPOLATION_STENCIL_MODE_ == _PIC_OUTPUT__CELL_CORNER_INTERPOLATION_STENCIL_MODE__LINEAR_) {
      PIC::Mesh::mesh.GetCenterNodesInterpolationCoefficients=PIC::Mesh::GetCenterNodesInterpolationCoefficients;
    }
  }

  //flush the sampling buffers
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::flushCollectingSamplingBuffer(block->GetCenterNode(LocalCellNumber));
          PIC::Mesh::flushCompletedSamplingBuffer(block->GetCenterNode(LocalCellNumber));
        }

        if (DIM==1) break;
      }

      if ((DIM==1)||(DIM==2)) break;
    }
  }

  //init the vector of flag that determins output of the data files
  PIC::Sampling::SaveOutputDataFile=new bool[PIC::nTotalSpecies];
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) PIC::Sampling::SaveOutputDataFile[spec]=true;

  //init the counter of the injected particles and the injection rates
  PIC::BC::nInjectedParticles=new long int[PIC::nTotalSpecies];
  PIC::BC::ParticleProductionRate=new double [PIC::nTotalSpecies];
  PIC::BC::ParticleMassProductionRate=new double [PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    PIC::BC::nInjectedParticles[spec]=0,PIC::BC::ParticleProductionRate[spec]=0.0;
    PIC::BC::ParticleMassProductionRate[spec]=0.0;
  }

  //init the model of volume particle injections
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  PIC::VolumeParticleInjection::Init();
#endif

  //init the background atmosphere model
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
  PIC::MolecularCollisions::BackgroundAtmosphere::Init_AfterParser();
#endif

  //init particle trajectory sampling
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::Init();
#endif

  //init the velocity distribution sample procedure
#if _PIC_VELOCITY_DISTRIBUTION_SAMPLING_MODE_ == _PIC_MODE_ON_
  PIC::DistributionFunctionSample::Init();
#endif

  //init the energy distribution sample procedure
#if _PIC_ENERGY_DISTRIBUTION_SAMPLING_RELATIVISTIC_MODE_ == _PIC_MODE_ON_
  PIC::EnergyDistributionSampleRelativistic::Init();
#endif

  //init the pitch angle distribution sample procedure
#if _PIC_PITCH_ANGLE_DISTRIBUTION_SAMPLING_MODE_ == _PIC_MODE_ON_
  PIC::PitchAngleDistributionSample::Init();
#endif

  //init the particle buffer
  //the default value of 'RequestedParticleBufferLength' is -1, which can be chganged by ampsConfig.pl
  //if 'RequestedParticleBufferLength' is -1 than the particle buffer is not initialized
  const int RequestedParticleBufferLength=-1;

  if (RequestedParticleBufferLength!=-1) PIC::ParticleBuffer::Init(RequestedParticleBufferLength);

  //copy the cut face information between the neighouring nodes 
  if (Mesh::IrregularSurface::nCutFaceInformationCopyAttempts!=0) {
    for (int i=0;i<Mesh::IrregularSurface::nCutFaceInformationCopyAttempts;i++) Mesh::IrregularSurface::CopyCutFaceInformation(); 
  }

  //init the ray tracking module if needed
  PIC::RayTracing::Init();

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
  return 0.0;
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




