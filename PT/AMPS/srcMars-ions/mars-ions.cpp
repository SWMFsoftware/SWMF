//$Id$

#include "mars-ions.h"
#include "pic.h"



// user defined global time step



double MarsIon::UserGlobalTimeStep = -1.0;

char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MARS";

//  injection boundary condition
double MarsIon::InjectionVelocity[3] = {26.3E1, 0.0, -2.3E1};
double MarsIon::InjectionNDensity    = 0.18E-6;
double MarsIon::InjectionTemperature = .6519;

// computational domain size
double MarsIon::DomainXMin[3] = {-2.25E14,-2.25E14,-2.25E14};
double MarsIon::DomainXMax[3] = { 2.25E14, 2.25E14, 2.25E14};
double MarsIon::DomainDXMin   = 1.8E13;
double MarsIon::DomainDXMax   = 1.8E13;


// OUTPUT ---------------------------------------------------------------------
int MarsIon::Output::TotalDataLength = 0;

int MarsIon::Output::OplusSource::RateOffset=-1;
int MarsIon::Output::OplusSource::BulkVelocityOffset=-1;
int MarsIon::Output::OplusSource::TemperatureOffset=-1;


//prameters of the flux sampling sphres
bool MarsIon::Sampling::SphericalShells::SamplingMode=false;
double MarsIon::Sampling::SphericalShells::SampleSphereRadii[MarsIon::Sampling::SphericalShells::nSphericalShells];
double MarsIon::Sampling::SphericalShells::minSampledEnergy=0.0;
double MarsIon::Sampling::SphericalShells::maxSampledEnergy=0.0;
double MarsIon::Sampling::SphericalShells::dLogEnergy=0.0;
int MarsIon::Sampling::SphericalShells::nSampledLevels=0;
cInternalSphericalData MarsIon::Sampling::SphericalShells::SamplingSphericlaShell[MarsIon::Sampling::SphericalShells::nSphericalShells];
int Earth_Sampling_Fluency_nSampledLevels=0;

void MarsIon::Output::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"OplusSource: Rate\",\"OplusSource: Bulk Velocity[0]\",\"OplusSource: Bulk Velocity[1]\",\"OplusSource: Bulk Velocity[2]\",\"OplusSource: Temperature\"");
}

void MarsIon::Output::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode){

  double S1=0.0, S2[3]={0.0,0.0,0.0}, S3=0.0;
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {

    S1+=(*((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+OplusSource::RateOffset)))*InterpolationCoeficients[i];

    for(idim=0 ; idim<3; idim++)
      S2[idim]+=(*(idim+(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+OplusSource::BulkVelocityOffset)))*InterpolationCoeficients[i];

    S3+=(*((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+OplusSource::TemperatureOffset)))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::RateOffset,&S1,sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::BulkVelocityOffset,S2,3*sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::TemperatureOffset,&S3,sizeof(double));

}

void MarsIon::Output::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode){
  double t;

  //Source Rate
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::RateOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //Source Bulk Velocity
  for(int idim=0; idim < 3; idim++){
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(idim+(double*)(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::BulkVelocityOffset));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }

  //Source Temperature
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+OplusSource::TemperatureOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

}

int MarsIon::Output::RequestDataBuffer(int offset){
  OplusSource::RateOffset=offset;
  TotalDataLength = 1;
  offset+=sizeof(double);

  OplusSource::BulkVelocityOffset=offset;
  TotalDataLength+=3;
  offset+=3*sizeof(double);

  OplusSource::TemperatureOffset=offset;
  TotalDataLength++;
  offset+=sizeof(double);

  return TotalDataLength*sizeof(double);
}


void MarsIon::Output::Init() {
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
}






//initialization of the model
void MarsIon::Init_BeforeParser() {

  //init the output module of the model
  MarsIon::Output::Init();

  //set the injection function
  PIC::BC::UserDefinedParticleInjectionFunction=SourceProcesses::InjectParticles;

  //set output of the particle flux sampling
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SamplingManager,Sampling::PrintManager);
}

void MarsIon::Init_AfterParser() {
  //set sampling of the partilce flux
  if (Sampling::SphericalShells::SamplingMode==true) {
    //set parameters of the fluency sampling
    Earth_Sampling_Fluency_nSampledLevels=Sampling::SphericalShells::nSampledLevels;
    Sampling::SphericalShells::dLogEnergy=log10(Sampling::SphericalShells::maxSampledEnergy/Sampling::SphericalShells::minSampledEnergy)/Sampling::SphericalShells::nSampledLevels;

    //init the shells
     for (int iShell=0;iShell<Sampling::SphericalShells::nSphericalShells;iShell++) {
       double sx0[3]={0.0,0.0,0.0};

       Sampling::SphericalShells::SamplingSphericlaShell[iShell].PrintDataStateVector=Sampling::SphericalShells::Output::PrintDataStateVector;
       Sampling::SphericalShells::SamplingSphericlaShell[iShell].PrintTitle=Sampling::SphericalShells::Output::PrintTitle;
       Sampling::SphericalShells::SamplingSphericlaShell[iShell].PrintVariableList=Sampling::SphericalShells::Output::PrintVariableList;

       Sampling::SphericalShells::SamplingSphericlaShell[iShell].SetSphereGeometricalParameters(sx0,Sampling::SphericalShells::SampleSphereRadii[iShell]);
       Sampling::SphericalShells::SamplingSphericlaShell[iShell].Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,1,Sampling::SphericalShells::SamplingSphericlaShell+iShell);
     }
   }

}




//init the background data from that loaded from TECPLOT
void MarsIon::InitBackgroundData() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;
  int i,j,k,LocalCellNumber,idim,thread;
  char *data;
  double *xCell;

  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    node=(thread==PIC::ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    while (node!=NULL) {
      block=node->block;

      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
        LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        cell=block->GetCenterNode(LocalCellNumber);

        if (cell!=NULL) {
          data=cell->GetAssociatedDataBufferPointer();
          xCell=cell->GetX();

          *(double*)(data+Output::OplusSource::RateOffset)=SourceProcesses::GetCellInjectionRate(_O_PLUS_SPEC_,xCell);
          for (idim=0;idim<3;idim++) *(double*)(data+Output::OplusSource::BulkVelocityOffset+idim*sizeof(double))=0.0;
          *(double*)(data+Output::OplusSource::TemperatureOffset)=100.0;
        }
      }

      node=node->nextNodeThisThread;
    }
  }
}


//particle mover: call relativistic Boris, ans sample particle flux
int MarsIon::ParticleMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double xInit[3],xFinal[3];
  int res,iShell;

  PIC::ParticleBuffer::GetX(xInit,ptr);

  switch (PIC::ParticleBuffer::GetI(ptr)) {
  case _ELECTRON_SPEC_:
    res=PIC::Mover::GuidingCenter::Mover_SecondOrder(ptr,dtTotal,startNode);
    break;
  default:
    res=PIC::Mover::Relativistic::Boris(ptr,dtTotal,startNode);
  }

  if ((Sampling::SphericalShells::SamplingMode==true)&&(res==_PARTICLE_MOTION_FINISHED_)) {
    //trajectory integration was successfully completed
    double r0,r1;

    PIC::ParticleBuffer::GetX(xFinal,ptr);
    r0=sqrt(xInit[0]*xInit[0]+xInit[1]*xInit[1]+xInit[2]*xInit[2]);
    r1=sqrt(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2]);

    //search for a shell that was intersected by the trajectory segment
    for (iShell=0;iShell<Sampling::SphericalShells::nSphericalShells;iShell++) if ((r0-Sampling::SphericalShells::SampleSphereRadii[iShell])*(r1-Sampling::SphericalShells::SampleSphereRadii[iShell])<=0.0) {
      //the trajectory segment has intersected a sampling shell
      //get the point of intersection
      double l[3],t=-1.0,a=0.0,b=0.0,c=0.0,d4,t1,t2;
      int idim;

      for (idim=0;idim<3;idim++) {
        a+=pow(xFinal[idim]-xInit[idim],2);
        b+=xInit[idim]*(xFinal[idim]-xInit[idim]);
        c+=pow(xInit[idim],2);
      }

      c-=pow(Sampling::SphericalShells::SampleSphereRadii[iShell],2);
      d4=sqrt(b*b-a*c);

      if ((isfinite(d4)==true)&&(a!=0.0)) {
        //the determinent is real
        t1=(-b+d4)/a;
        t2=(-b-d4)/a;

        if ((0.0<t1)&&(t1<1.0)) t=t1;
        if ((0.0<t2)&&(t2<1.0)) if ((t<0.0)||(t2<t)) t=t2;

        if (t>0.0) {
          //the point of intersection is found. Sample flux and rigidity at that location
          double xIntersection[3];
          long int iSurfaceElementNumber,iZenithElement,iAzimuthalElement;

          for (idim=0;idim<3;idim++) xIntersection[idim]=xInit[idim]+t*(xFinal[idim]-xInit[idim]);

          Sampling::SphericalShells::SamplingSphericlaShell[iShell].GetSurfaceElementProjectionIndex(xIntersection,iZenithElement,iAzimuthalElement);
          iSurfaceElementNumber=Sampling::SphericalShells::SamplingSphericlaShell[iShell].GetLocalSurfaceElementNumber(iZenithElement,iAzimuthalElement);

          //increment sampling counters
          double ParticleWeight,Rigidity;
          int spec;


          spec=PIC::ParticleBuffer::GetI(ptr);
          ParticleWeight=startNode->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);


          #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
          double *t;
          t=Sampling::SamplingSphericlaShell[iShell].Flux[spec]+iSurfaceElementNumber;

          #pragma omp atomic
          *t+=ParticleWeight;
          #else
          Sampling::SphericalShells::SamplingSphericlaShell[iShell].Flux[spec][iSurfaceElementNumber]+=ParticleWeight;
          #endif

          //Fluence Sampling Energy level
          int iEnergyLevel=-1;
          double Energy,Speed,v[3];

          PIC::ParticleBuffer::GetV(v,ptr);
          Speed=Vector3D::Length(v);
          Energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

          if ((Sampling::SphericalShells::minSampledEnergy<Energy)&&(Energy<Sampling::SphericalShells::maxSampledEnergy)) {
            iEnergyLevel=(int)(log10(Energy/Sampling::SphericalShells::minSampledEnergy)/Sampling::SphericalShells::dLogEnergy);

            if ((iEnergyLevel<0)||(iEnergyLevel>=Sampling::SphericalShells::nSampledLevels)) iEnergyLevel=-1;
          }
          else iEnergyLevel=-1;

          if (r0>Sampling::SphericalShells::SampleSphereRadii[iShell]) {
            //the particle moves toward hte Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            t=Sampling::SamplingSphericlaShell[iShell].ParticleFluxDown[spec]+iSurfaceElementNumber;

            #pragma omp atomic
            *t+=ParticleWeight;

            if (iEnergyLevel!=-1) {
              t=Sampling::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber]+iEnergyLevel;

              #pragma omp atomic
              *t+=ParticleWeight;

              //Sampling::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            }


            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].ParticleFluxDown[spec][iSurfaceElementNumber]+=ParticleWeight;
            if (iEnergyLevel!=-1) Sampling::SphericalShells::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_

          }
          else {
            //the particle moves outward the Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            t=Sampling::SamplingSphericlaShell[iShell].ParticleFluxUp[spec]+iSurfaceElementNumber;

            #pragma omp atomic
            *t+=ParticleWeight;

            if (iEnergyLevel!=-1) {
              t=Sampling::SamplingSphericlaShell[iShell].ParticleFluencyUp[spec][iSurfaceElementNumber]+iEnergyLevel;

              #pragma omp atomic
              *t+=ParticleWeight;

              //Sampling::SamplingSphericlaShell[i].ParticleFluencyUp[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            }

            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].ParticleFluxUp[spec][iSurfaceElementNumber]+=ParticleWeight;
            if (iEnergyLevel!=-1) Sampling::SphericalShells::SamplingSphericlaShell[iShell].ParticleFluencyUp[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_
          }

          //sampling of the speed of the particle
          if (r0>Sampling::SphericalShells::SampleSphereRadii[iShell]) {
            //the particle moves toward the Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            double *t,c;

            t=&Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleDown[spec][iSurfaceElementNumber].SampledParticleWeight;

            #pragma omp atomic
            *t+=ParticleWeight;

            t=&Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleDown[spec][iSurfaceElementNumber].SampledParticleSpeed;
            c=Speed*ParticleWeight;

            #pragma omp atomic
            *t+=c;

            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleDown[spec][iSurfaceElementNumber].SampledParticleWeight+=ParticleWeight;
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleDown[spec][iSurfaceElementNumber].SampledParticleSpeed+=Speed*ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_
          }
          else {
            //the particle moves outward of the Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            double *t,c;

            t=&Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleUp[spec][iSurfaceElementNumber].SampledParticleWeight;

            #pragma omp atomic
            *t+=ParticleWeight;

            t=&Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleUp[spec][iSurfaceElementNumber].SampledParticleSpeed;
            c=Speed*ParticleWeight;

            #pragma omp atomic
            *t+=c;

            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleUp[spec][iSurfaceElementNumber].SampledParticleWeight+=ParticleWeight;
            Sampling::SphericalShells::SamplingSphericlaShell[iShell].VelocitySampleUp[spec][iSurfaceElementNumber].SampledParticleSpeed+=Speed*ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_
          }

          //calculate particle rigidity
          double ElectricCharge;

          ElectricCharge=fabs(PIC::MolecularData::GetElectricCharge(spec));

          if (ElectricCharge>0.0) {
            Rigidity=Relativistic::Speed2Momentum(Speed,PIC::MolecularData::GetMass(spec))/PIC::MolecularData::GetElectricCharge(spec);

            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            #pragma omp critical
            #endif
            {
              if ((Sampling::SphericalShells::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber]<0.0) || (Rigidity<Sampling::SphericalShells::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber])) {
                Sampling::SphericalShells::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber]=Rigidity;
              }
            }

          }
        }
      }

    }
  }
}












/*
// Loss -------------------------------------------------------------------------------------

double OH::Loss::LifeTime(double *x, int spec, long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node){

  double PlasmaNumberDensity, PlasmaPressure, PlasmaTemperature;
  double PlasmaBulkVelocity[3];

  double lifetime=0.0;

  PIC::CPLR::InitInterpolationStencil(x,node);

  PlasmaNumberDensity = PIC::CPLR::GetBackgroundPlasmaNumberDensity();
  PlasmaPressure      = PIC::CPLR::GetBackgroundPlasmaPressure();
  PlasmaTemperature   = PlasmaPressure / (Kbol * PlasmaNumberDensity);
  PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);

  PhotolyticReactionAllowedFlag=true;

  // the model has hydrogen only
  if (spec!=_H_SPEC_) {
    PhotolyticReactionAllowedFlag=false;
    return -1.0;
  }

  // velocity of a particle
  double v[3];
  PIC::ParticleBuffer::GetV(v,ptr);
  //-------------------
  switch (spec) {
  case _H_SPEC_:
    lifetime= ChargeExchange::LifeTime(_H_SPEC_, v, PlasmaBulkVelocity, PlasmaTemperature, PlasmaNumberDensity);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown specie");
  }

  return lifetime;

}

int OH::Loss::ReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
  //for one lost particle one new particle is generated
  //----------------------------------------------------------------------
  //inject the products of the reaction
  double ParentTimeStep,ParentParticleWeight;

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  ParentParticleWeight=0.0;
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
  ParentTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif


  //account for the parent particle correction factor
  ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

  //the particle buffer used to set-up the new particle data
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  //copy the state of the initial parent particle into the new-daugher particle
  PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*)tempParticleData,ParticleData);

  // injection time of the particle (after the beginning of AMPS' time step)
  double ModelParticleInjectionRate,TimeCounter=0.0,TimeIncrement,ProductWeightCorrection=1.0;
  ModelParticleInjectionRate=1.0/ParentTimeStep;
  TimeIncrement=-log(rnd())/ModelParticleInjectionRate * rnd();
  TimeCounter += TimeIncrement;

  //generate a particle
  // new particle comes from solar wind and has velocity ~ plasma bulk velocity
  double PlasmaBulkVelocity[3];
  {
    PIC::CPLR::InitInterpolationStencil(xFinal,node);
    PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);

    // charge exchange process transfers momentum and energy to plasma
    PIC::Mesh::cDataCenterNode *CenterNode;
    char *offset;

    CenterNode=PIC::Mesh::Search::FindCell(xFinal); ///node->block->GetCenterNode(nd);
    offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
    double v2 = 0.0, plasmav2 = 0.0;
    double c = ParentParticleWeight
                   / PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]
                   / CenterNode->Measure;
    *((double*)(offset+OH::Output::ohSourceDensityOffset)) += 0.0;
    for(int idim=0; idim<3; idim++){
      *(idim + (double*)(offset+OH::Output::ohSourceMomentumOffset)) +=
	c*_MASS_(_H_)*(vFinal[idim]-PlasmaBulkVelocity[idim]);
      v2      +=vFinal[idim]*vFinal[idim];
      plasmav2+=PlasmaBulkVelocity[idim]*PlasmaBulkVelocity[idim];
    }
    *((double*)(offset+OH::Output::ohSourceEnergyOffset)) +=
      c*0.5*_MASS_(_H_)*(v2-plasmav2);
  }

  PIC::ParticleBuffer::SetX(xFinal,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetV(PlasmaBulkVelocity,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
#endif

  //get and injection into the system the new model particle
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,ParentTimeStep-TimeCounter,node,true);

  return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;

}


void OH::Init_BeforeParser(){
  Exosphere::Init_BeforeParser();
  OH::Output::Init();

  //set the coupling procedure
  PIC::CPLR::SWMF::SendCenterPointData.push_back(Coupling::Send);
}


//-----------------------------------------------------------------------------
//substitutes for Exosphere functions
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_],Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_],Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_];

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=1;
  }

  return nVariables;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;


  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    //    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  //do nothing
}

double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {
  ReemissionParticleFraction=0.0;

  return 1.0;
}


double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {


  return 100.0;
}
*/
