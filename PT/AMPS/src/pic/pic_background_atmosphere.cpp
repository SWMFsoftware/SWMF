//====================================================
//$Id$
//====================================================
//the functions that model colecular collisions with the background atmosphere

#include "pic.h"

#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_

int PIC::MolecularCollisions::BackgroundAtmosphere::LocalTotalCollisionFreqSamplingOffset=-1;

//init the model
void PIC::MolecularCollisions::BackgroundAtmosphere::Init_BeforeParser() {
  //set up the output of the model
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
  PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);
}

//Request Sampling Buffers
int PIC::MolecularCollisions::BackgroundAtmosphere::RequestSamplingData(int offset) {
  LocalTotalCollisionFreqSamplingOffset=offset;
  return PIC::nTotalSpecies*sizeof(double);
}

void PIC::MolecularCollisions::BackgroundAtmosphere::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Total Background Atmosphere Collison Freq [per a real particle,per sec]\"");
}

void PIC::MolecularCollisions::BackgroundAtmosphere::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

  struct cDataExchengeBuffer {
    double TotalCollisionFreq;
  };

  cDataExchengeBuffer buffer;

  if (pipe->ThisThread==CenterNodeThread) {
    buffer.TotalCollisionFreq=*(DataSetNumber+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv((char*)&buffer,sizeof(cDataExchengeBuffer),CenterNodeThread);
    }

    if (PIC::LastSampleLength!=0) buffer.TotalCollisionFreq/=PIC::LastSampleLength;

    fprintf(fout,"%e  ",buffer.TotalCollisionFreq);
  }
  else {
    pipe->send((char*)&buffer,sizeof(cDataExchengeBuffer));
  }
}

void PIC::MolecularCollisions::BackgroundAtmosphere::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i,spec;
  double c;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //interpolate the local productions rate
    double InterpoaltedTotalCollisionFreq=0.0;

    //interpolate the sampled data
    for (i=0;i<nInterpolationCoeficients;i++) {
      c=InterpolationCoeficients[i];

      InterpoaltedTotalCollisionFreq+=c*(*(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
    }

    //stored the interpolated data in the associated data buffer
    *(spec+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpoaltedTotalCollisionFreq;
  }
}


/*
void PIC::MolecularCollisions::BackgroundAtmosphere::CollisionProcessor() {
  int i,j,k;

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  int centerNodeIndexCounter,LocalCellNumber;

  const int ParticleBufferLength=10000;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataCenterNode *cell;

  //the buffer of particles that occuping the local cell
  int s,ptr,nParticles[PIC::nTotalSpecies],spec,idim;
  PIC::ParticleBuffer::byte *ParticleData[PIC::nTotalSpecies][ParticleBufferLength],*pdata;
  long int ParticleList[PIC::nTotalSpecies][ParticleBufferLength];
  bool ParticlesStoredInBuffer[PIC::nTotalSpecies];
  double WeightCorrectionMin[PIC::nTotalSpecies],WeightCorrection,WeightCorrectionMax[PIC::nTotalSpecies],cr2,cellCumulativeParticleWeight[PIC::nTotalSpecies];
  double vModelParticle[3],vBackgroundParticle[3];
  PIC::Mesh::cDataBlockAMR *block;

  //the temporary particle that represents the background atmosphere
  long int tempBackgroundAtmosphereParticle;
  PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData;

  tempBackgroundAtmosphereParticle=PIC::ParticleBuffer::GetNewParticle();
  BackgroundAtmosphereParticleData=PIC::ParticleBuffer::GetParticleDataPointer(tempBackgroundAtmosphereParticle);

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  //loop through all nodes and cells on the currect processor
  while (node!=NULL) {
    block=node->block;

    for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
      //sort the particle from the cell
      for (s=0;s<PIC::nTotalSpecies;s++) nParticles[s]=0,ParticlesStoredInBuffer[s]=true,WeightCorrectionMin[s]=-1.0,WeightCorrectionMax[s]=-1.0,cellCumulativeParticleWeight[s]=0.0;

      LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
      cell=block->GetCenterNode(LocalCellNumber);
      ptr=cell->FirstCellParticle;

      while (ptr!=-1) {
        pdata=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        s=PIC::ParticleBuffer::GetI(pdata);

        WeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pdata);
        if ((WeightCorrectionMin[s]<0.0)||(WeightCorrection<WeightCorrectionMin[s])) WeightCorrectionMin[s]=WeightCorrection;
        if (WeightCorrectionMax[s]<WeightCorrection) WeightCorrectionMax[s]=WeightCorrection;

        if (nParticles[s]<ParticleBufferLength-1) {
          ParticleData[s][nParticles[s]]=pdata;
          ParticleList[s][nParticles[s]]=ptr;
        }
        else {
          ParticlesStoredInBuffer[s]=false;
        }

        ++nParticles[s];
        cellCumulativeParticleWeight[s]+=WeightCorrection;
        ptr=PIC::ParticleBuffer::GetNext(ptr);
      }

      //simulate collisions with the background atmosphere
      for (s=0;s<PIC::nTotalSpecies;s++) if (nParticles[s]!=0) for (int BackgroundSpecieNumber=0;BackgroundSpecieNumber<GetTotalNumberBackgroundSpecies();BackgroundSpecieNumber++) {
        const int nTotalTest=5*nParticles[s];

        double SigmaCr,SigmaCrMax=0.0,WeightCorrection=0.0;
        int nTest,nptr;

        //calcualte the majorant collision frequentcy
        for (nTest=0;nTest<nTotalTest;nTest++) {
          //select randomly a model particle from the cell
          nptr=(int)(rnd()*nParticles[s]);

          if (ParticlesStoredInBuffer[s]==true) pdata=ParticleData[s][nptr]; //take the particle from the buffer
          else { //get the particle directly from the list
            ptr=cell->FirstCellParticle;

            while (ptr!=-1) {
              pdata=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              spec=PIC::ParticleBuffer::GetI(pdata);

              if (s==spec) {
                if (nptr==0) break;
                else --nptr;
              }

                ptr=PIC::ParticleBuffer::GetNext(pdata);
            }
          }


          //generate particle that represents the background atmosphere and calcualte the cross section
          GenerateBackgoundAtmosphereParticle(BackgroundAtmosphereParticleData,BackgroundSpecieNumber,cell,node);

          PIC::ParticleBuffer::GetV(vModelParticle,pdata);
          PIC::ParticleBuffer::GetV(vBackgroundParticle,BackgroundAtmosphereParticleData);
          for (idim=0,cr2=0.0;idim<3;idim++) cr2+=pow(vModelParticle[idim]-vBackgroundParticle[idim],2);

          SigmaCr=GetCollisionCrossSectionBackgoundAtmosphereParticle(pdata,s,cell,node,BackgroundAtmosphereParticleData,BackgroundSpecieNumber)*sqrt(cr2);
          if (SigmaCr>SigmaCrMax) SigmaCrMax=SigmaCr;
        }

        //calcualte the collision cross section
        double MajorantCollisionFreq;

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__ISOTROPIC_
        MajorantCollisionFreq=2*GetCellMeanBackgroundNumberDensity(BackgroundSpecieNumber,cell,node)*cellCumulativeParticleWeight[s]*SigmaCrMax/WeightCorrectionMin[s];
#elif _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_
        MajorantCollisionFreq=2*GetCellMaximumBackgroundNumberDensity(BackgroundSpecieNumber,cell,node)*cellCumulativeParticleWeight[s]*SigmaCrMax/WeightCorrectionMin[s];
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

        //model collisions with the background atmosphere
        double timeCounter,localTimeStep;
        double massModelParticle,massBackgroundParticle,Vrel[3]={0.0,0.0,0.0},Vcm[3]={0.0,0.0,0.0},am,mr;

        timeCounter=0.0;
        localTimeStep=node->block->GetLocalTimeStep(s);

        massModelParticle=PIC::MolecularData::GetMass(s);
        massBackgroundParticle=GetBackgroundMolecularMass(BackgroundSpecieNumber);

        mr=massModelParticle*massBackgroundParticle/(massModelParticle+massBackgroundParticle);
        am=massModelParticle+massBackgroundParticle;

        if (MajorantCollisionFreq>0.0) while ((timeCounter-=log(rnd())/MajorantCollisionFreq)<localTimeStep) {
          do {
            //select randomly a model particle from the cell
            nptr=(int)(rnd()*nParticles[s]);

            if (ParticlesStoredInBuffer[s]==true) {
              pdata=ParticleData[s][nptr]; //take the particle from the buffer
              ptr=ParticleList[s][nptr];
            }
            else { //get the particle directly from the list
              ptr=cell->FirstCellParticle;

              while (ptr!=-1) {
                pdata=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
                spec=PIC::ParticleBuffer::GetI(pdata);

                if (s==spec) {
                  if (nptr==0) break;
                  else --nptr;
                }

                ptr=PIC::ParticleBuffer::GetNext(ptr);
              }
            }

            WeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pdata);
          }
          while (WeightCorrection/WeightCorrectionMax[s]<rnd());

          //generate particle that represents the background atmosphere and calcualte the cross section
          GenerateBackgoundAtmosphereParticle(BackgroundAtmosphereParticleData,BackgroundSpecieNumber,cell,node);

          PIC::ParticleBuffer::GetV(vModelParticle,pdata);
          PIC::ParticleBuffer::GetV(vBackgroundParticle,BackgroundAtmosphereParticleData);

          for (idim=0,cr2=0.0;idim<3;idim++) {
            Vrel[idim]=vModelParticle[idim]-vBackgroundParticle[idim];
            Vcm[idim]=(massModelParticle*vModelParticle[idim]+massBackgroundParticle*vBackgroundParticle[idim])/am;

            cr2+=pow(Vrel[idim],2);
          }

          SigmaCr=GetCollisionCrossSectionBackgoundAtmosphereParticle(pdata,s,cell,node,BackgroundAtmosphereParticleData,BackgroundSpecieNumber)*sqrt(cr2);

          //check if the collision is possible
          if ((rnd()*SigmaCrMax>SigmaCr)||(WeightCorrectionMin[s]/WeightCorrection<rnd())) continue;

          //account for the variation of the collision frequentcy within the cell
#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_
          double xParticle[3];

          PIC::ParticleBuffer::GetX(xParticle,pdata);
          if (GetCellLocalBackgroundNumberDensity(xParticle,BackgroundSpecieNumber,cell,node)/GetCellMaximumBackgroundNumberDensity(BackgroundSpecieNumber,cell,node)<rnd()) continue;
#endif

          //redistribute the relative velocity of the collided particles
          double Vrc,V[3];
          double CosKsi,SinKsi,CosEps,SinEps,D,c;

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_
          CosKsi=2.0*rnd()-1.0;
          SinKsi=sqrt(1.0-CosKsi*CosKsi);
#elif _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__USER_DEFINED_
          double VelocityScatteringAngle;

          VelocityScatteringAngle=UserDefinedVelocityScatteringAngle();
          CosKsi=cos(VelocityScatteringAngle);
          SinKsi=sin(VelocityScatteringAngle);
#else
          exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

          c=2*Pi*rnd();
          SinEps=sin(c);
          CosEps=cos(c);

          D=sqrt(Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

          if (D>1.0E-6) {
            Vrc=sqrt(Vrel[0]*Vrel[0]+Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

            V[0]=CosKsi*Vrel[0]+SinKsi*SinEps*D;
            V[1]=CosKsi*Vrel[1]+SinKsi*(Vrc*Vrel[2]*CosEps-Vrel[0]*Vrel[1]*SinEps)/D;
            V[2]=CosKsi*Vrel[2]-SinKsi*(Vrc*Vrel[1]*CosEps+Vrel[0]*Vrel[2]*SinEps)/D;
          }
          else {
            V[0]=CosKsi*Vrel[0];
            V[1]=SinKsi*CosEps*Vrel[0];
            V[2]=SinKsi*SinEps*Vrel[0];
          }

          Vrel[0]=V[0];
          Vrel[1]=V[1];
          Vrel[2]=V[2];

          //the collision between the model particle and the particle from the background atmosphere has occured
          for (idim=0;idim<3;idim++) {
            vModelParticle[idim]=Vcm[idim]+massBackgroundParticle/am*Vrel[idim];
            vBackgroundParticle[idim]=Vcm[idim]-massModelParticle/am*Vrel[idim];
          }

          //check if the 'background' particle should be kept in the system
          PIC::ParticleBuffer::SetV(vBackgroundParticle,BackgroundAtmosphereParticleData);

          if (KeepBackgroundAtmosphereParticle(BackgroundAtmosphereParticleData)==true) {
            long int p;

            p=PIC::ParticleBuffer::GetNewParticle(cell->FirstCellParticle);
            PIC::ParticleBuffer::CloneParticle(p,tempBackgroundAtmosphereParticle);

            //set particle weight
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
            double wcorrection;
            int bspec;

            bspec=PIC::ParticleBuffer::GetI(p);
            wcorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pdata)*block->GetLocalParticleWeight(s)/block->GetLocalParticleWeight(bspec);
            PIC::ParticleBuffer::SetIndividualStatWeightCorrection(wcorrection,p);

            MajorantCollisionFreq*=1.0+WeightCorrection/cellCumulativeParticleWeight[s];
            cellCumulativeParticleWeight[s]+=WeightCorrection;
            ++nParticles[s];

            if ((ParticlesStoredInBuffer[s]==true)&&(nParticles[s]<=ParticleBufferLength-1)) {
              ParticleData[s][nParticles[s]-1]=PIC::ParticleBuffer::GetParticleDataPointer(p);
              ParticleList[s][nParticles[s]-1]=p;
            }
            else ParticlesStoredInBuffer[s]=false;

#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      exit(__LINE__,__FILE__,"Error: not implementd");
#else
            exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif
          }

          //set new velocity vector to the model particle
          PIC::ParticleBuffer::SetV(vModelParticle,pdata);

          //check if the model particle should be removed from the system
          if (KeepBackgroundAtmosphereParticle(pdata)==false) {
            //the particle should be removed
            long int next,prev;

            next=PIC::ParticleBuffer::GetNext(pdata);
            prev=PIC::ParticleBuffer::GetPrev(pdata);

            //remove the particle from the particle buffer 'ParticleData[][]' and adjust the colision frequentcy
            if ((ParticlesStoredInBuffer[s]==true)&&(nParticles[s]!=1)) {
              ParticleData[s][nptr]=ParticleData[s][nParticles[s]-1];
              ParticleList[s][nptr]=ParticleList[s][nParticles[s]-1];
            }

            WeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pdata);
            MajorantCollisionFreq*=1.0-WeightCorrection/cellCumulativeParticleWeight[s];
            cellCumulativeParticleWeight[s]-=WeightCorrection;
            --nParticles[s];

            //reconnect particles from the list
            if (prev==-1) cell->FirstCellParticle=next;
            else PIC::ParticleBuffer::SetNext(next,prev);

            if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

            PIC::ParticleBuffer::DeleteParticle(ptr);

            //if no particles of the given specie has left than break the loop that models the collisions
            if ((nParticles[s]==0)||(MajorantCollisionFreq<=0.0)) break;
          }

        }
      }

    }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

    node=node->nextNodeThisThread;
  }

  //delete the temporary particle representing the background atmosphere
  PIC::ParticleBuffer::DeleteParticle(tempBackgroundAtmosphereParticle);
}
*/


void PIC::MolecularCollisions::BackgroundAtmosphere::CollisionProcessor() {
  int i,j,k;

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  //sample collisino frequentcy
  double TotalProspectiveCollisionParticleWeight[PIC::nTotalSpecies],TotalOccurringCollisionParticleWeight[PIC::nTotalSpecies];

  int centerNodeIndexCounter,LocalCellNumber;

  const int ParticleBufferLength=50000;

  struct cCollidingParticleList {
    long int Particle;
    double CollisionTimeFraction;
  };

  cCollidingParticleList CollidingParticleList[ParticleBufferLength];
  long int nCollidingParticles;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataCenterNode *cell;

  //the buffer of particles that occuping the local cell
  long int modelParticle;
  int BackgroundSpecieNumber,spec,idim;
  PIC::ParticleBuffer::byte *modelParticleData;
  double vModelParticle[3],xModelParticle[3],vBackgroundParticle[3],particleCollisionTime,cr2;
  PIC::Mesh::cDataBlockAMR *block;

  //the temporary particle that represents the background atmosphere
  long int tempBackgroundAtmosphereParticle;
  PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData;

  tempBackgroundAtmosphereParticle=PIC::ParticleBuffer::GetNewParticle();
  BackgroundAtmosphereParticleData=PIC::ParticleBuffer::GetParticleDataPointer(tempBackgroundAtmosphereParticle);

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  //loop through all nodes and cells on the currect processor
  double MajorantCollisionFreq,SigmaCr,SigmaCrMax=0.0;
  double timeCounter,localTimeStep,TranslationalEnergy;
  double massModelParticle,massBackgroundParticle,Vrel[3]={0.0,0.0,0.0},Vcm[3]={0.0,0.0,0.0},am;


//--------------------------   DEBUG --------------------------
/*
static long int nCall=0;
*/
//--------------------------   END DEBUg -----------------------

  while (node!=NULL) {
    block=node->block;

//--------------------------   DEBUG --------------------------
/*
nCall++;

if ((nCall==1452339-1)&&(PIC::ThisThread==1)) {
  cout << __FILE__ << "$" << __LINE__ << endl;
}
*/

//--------------------------   END DEBUg -----------------------



    for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
      LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
      cell=block->GetCenterNode(LocalCellNumber);
      modelParticle=cell->FirstCellParticle;

      //init the list of colliding particles
      nCollidingParticles=0;
      modelParticle=cell->FirstCellParticle;

      for (spec=0;spec<PIC::nTotalSpecies;spec++) TotalProspectiveCollisionParticleWeight[spec]=0.0,TotalOccurringCollisionParticleWeight[spec]=0.0;


      while (modelParticle!=-1) {
        if (nCollidingParticles==ParticleBufferLength) exit(__LINE__,__FILE__,"Error: the value of 'ParticleBufferLength' is exeeded - too many particles in a cells. Increase the value of 'ParticleBufferLength'");

        CollidingParticleList[nCollidingParticles].Particle=modelParticle;
        CollidingParticleList[nCollidingParticles].CollisionTimeFraction=1.0;

        ++nCollidingParticles;
        modelParticle=PIC::ParticleBuffer::GetNext(modelParticle);
      }

_StartParticleCollisionLoop_:
      while (nCollidingParticles!=0) {
        modelParticle=CollidingParticleList[0].Particle;
        modelParticleData=PIC::ParticleBuffer::GetParticleDataPointer(modelParticle);

        spec=PIC::ParticleBuffer::GetI(modelParticleData);
        PIC::ParticleBuffer::GetV(vModelParticle,modelParticleData);
        PIC::ParticleBuffer::GetX(xModelParticle,modelParticleData);

        massModelParticle=PIC::MolecularData::GetMass(spec);
        localTimeStep=node->block->GetLocalTimeStep(spec);
        particleCollisionTime=CollidingParticleList[0].CollisionTimeFraction*localTimeStep;

        TotalProspectiveCollisionParticleWeight[spec]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

        //exclude the particle from the list of the model particle that still needs to be collided
        if (nCollidingParticles!=1) CollidingParticleList[0]=CollidingParticleList[nCollidingParticles-1];
        --nCollidingParticles;

        //simulate collisions with the background atmosphere
        for (BackgroundSpecieNumber=0;BackgroundSpecieNumber<GetTotalNumberBackgroundSpecies();BackgroundSpecieNumber++) {
          SigmaCrMax=GetSigmaCrMax(spec,BackgroundSpecieNumber,modelParticleData);
          MajorantCollisionFreq=GetBackgroundNumberDensity(BackgroundSpecieNumber,xModelParticle)*SigmaCrMax;


//TEST!!!!!!!!!!!!!!!
//MajorantCollisionFreq=(sqrt(xModelParticle[0]*xModelParticle[0]+xModelParticle[1]*xModelParticle[1]+xModelParticle[2]*xModelParticle[2])<700.0E3+_RADIUS_(_TARGET_)) ? 4E5*1.0E6*SigmaCrMax : 0.0;



#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_
          //if the majorant frequentcy large, reivalute it with the mean value of the background density calcualted along the trajectory of the particle with the sped xdMeanDensity;
          if (false) { //(MajorantCollisionFreq*particleCollisionTime>50.0) {
            double xProbe[3],meanValueOld,meanValueNew;
            int iDensityAveragingLevel;

            const int maxDensityAveragingLevels=6;


            for (idim=0;idim<DIM;idim++) xProbe[idim]=xModelParticle[idim]-vModelParticle[idim]*particleCollisionTime;
            meanValueOld=0.5*(GetBackgroundNumberDensity(BackgroundSpecieNumber,xModelParticle)+GetBackgroundNumberDensity(BackgroundSpecieNumber,xProbe));


            for (iDensityAveragingLevel=1;iDensityAveragingLevel<=maxDensityAveragingLevels;iDensityAveragingLevel++) {
              int idim,i,nTestPoints;

              nTestPoints=(1<<iDensityAveragingLevel);
              meanValueNew=meanValueOld*(1+(1<<(iDensityAveragingLevel-1)));

              for (i=1;i<=nTestPoints;i+=2) {
                for (idim=0;idim<DIM;idim++) xProbe[idim]=xModelParticle[idim]-vModelParticle[idim]*double(i)*particleCollisionTime/nTestPoints;     //'-' because the avaraging clong the 'past' trajectory

                meanValueNew+=GetBackgroundNumberDensity(BackgroundSpecieNumber,xProbe);
              }

              meanValueNew/=(1+nTestPoints);
              if (fabs(1.0-meanValueOld/meanValueNew)<0.2) break;
            }

            MajorantCollisionFreq=meanValueNew*SigmaCrMax;
          }
#endif


          timeCounter=0.0;

          massBackgroundParticle=GetBackgroundMolecularMass(BackgroundSpecieNumber);
          am=massModelParticle+massBackgroundParticle;

          if (MajorantCollisionFreq>0.0) while ((timeCounter-=log(rnd())/MajorantCollisionFreq)<particleCollisionTime) {
            //generate particle that represents the background atmosphere and calcualte the cross section
            GenerateBackgoundAtmosphereParticle(BackgroundAtmosphereParticleData,BackgroundSpecieNumber,cell,node);
            PIC::ParticleBuffer::GetV(vBackgroundParticle,BackgroundAtmosphereParticleData);

//TEST!!!!!!!!!!!!!!!
for (idim=0;idim<3;idim++) vBackgroundParticle[idim]=0.0;

            for (idim=0,cr2=0.0;idim<3;idim++) {
              Vrel[idim]=vModelParticle[idim]-vBackgroundParticle[idim];
              Vcm[idim]=(massModelParticle*vModelParticle[idim]+massBackgroundParticle*vBackgroundParticle[idim])/am;

              cr2+=pow(Vrel[idim],2);
            }

            TranslationalEnergy=0.5*massModelParticle*massBackgroundParticle/(massModelParticle+massBackgroundParticle)*cr2;
            SigmaCr=GetCollisionCrossSectionBackgoundAtmosphereParticle(spec,BackgroundSpecieNumber,modelParticleData,BackgroundAtmosphereParticleData,TranslationalEnergy,cr2)*sqrt(cr2);

            //check if the collision is possible
            if (rnd()*SigmaCrMax>SigmaCr) continue;

            //sample collision frecuentcy
            TotalOccurringCollisionParticleWeight[spec]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

            //redistribute the relative velocity of the collided particles
            double Vrc,V[3];
            double CosKsi,SinKsi,CosEps,SinEps,D,c;

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_
            CosKsi=2.0*rnd()-1.0;
            SinKsi=sqrt(1.0-CosKsi*CosKsi);
#elif _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__USER_DEFINED_
            double VelocityScatteringAngle;

            VelocityScatteringAngle=UserDefinedVelocityScatteringAngle(Vrel,TranslationalEnergy,spec,BackgroundSpecieNumber);
            CosKsi=cos(VelocityScatteringAngle);
            SinKsi=sin(VelocityScatteringAngle);
#else
            exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

            c=2*Pi*rnd();
            SinEps=sin(c);
            CosEps=cos(c);

            D=sqrt(Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

            if (D>1.0E-6) {
              Vrc=sqrt(Vrel[0]*Vrel[0]+Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

              V[0]=CosKsi*Vrel[0]+SinKsi*SinEps*D;
              V[1]=CosKsi*Vrel[1]+SinKsi*(Vrc*Vrel[2]*CosEps-Vrel[0]*Vrel[1]*SinEps)/D;
              V[2]=CosKsi*Vrel[2]-SinKsi*(Vrc*Vrel[1]*CosEps+Vrel[0]*Vrel[2]*SinEps)/D;
            }
            else {
              V[0]=CosKsi*Vrel[0];
              V[1]=SinKsi*CosEps*Vrel[0];
              V[2]=SinKsi*SinEps*Vrel[0];
            }

            Vrel[0]=V[0];
            Vrel[1]=V[1];
            Vrel[2]=V[2];

            //the collision between the model particle and the particle from the background atmosphere has occured
            for (idim=0;idim<3;idim++) {
              vModelParticle[idim]=Vcm[idim]+massBackgroundParticle/am*Vrel[idim];
              vBackgroundParticle[idim]=Vcm[idim]-massModelParticle/am*Vrel[idim];
            }

            //update velocities of the particles
            PIC::ParticleBuffer::SetV(vBackgroundParticle,BackgroundAtmosphereParticleData);
            PIC::ParticleBuffer::SetV(vModelParticle,modelParticle);

            //check if the 'background' particle should be kept in the system
            if (false) { ///(KeepBackgroundAtmosphereParticle(BackgroundAtmosphereParticleData)==true) {
              long int newParticle;
              PIC::ParticleBuffer::byte *newParticleData;

              newParticle=PIC::ParticleBuffer::GetNewParticle(cell->FirstCellParticle);
              PIC::ParticleBuffer::CloneParticle(newParticle,tempBackgroundAtmosphereParticle);
              newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

              //set hte position of the new particle to be the position of the original model particle
              PIC::ParticleBuffer::SetX(xModelParticle,newParticleData);

              //add the new particle to teh list of the particles that needes to be collided
              if (nCollidingParticles>ParticleBufferLength-1) exit(__LINE__,__FILE__,"Error: the particle buffer is overflown");

              CollidingParticleList[nCollidingParticles].Particle=newParticle;
              CollidingParticleList[nCollidingParticles].CollisionTimeFraction=1.0-timeCounter/localTimeStep;
              ++nCollidingParticles;

              //set particle weight
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
              double wcorrection;
              int bspec;

              bspec=PIC::ParticleBuffer::GetI(newParticleData);
              wcorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)*block->GetLocalParticleWeight(spec)/block->GetLocalParticleWeight(bspec);
              PIC::ParticleBuffer::SetIndividualStatWeightCorrection(wcorrection,newParticleData);
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      exit(__LINE__,__FILE__,"Error: not implementd");
#else
            exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif
           }

            //check if the model particle should be removed from the system
            if (KeepBackgroundAtmosphereParticle(modelParticleData)==false) {
              //the particle should be removed
              long int next,prev;

              next=PIC::ParticleBuffer::GetNext(modelParticleData);
              prev=PIC::ParticleBuffer::GetPrev(modelParticleData);

              //reconnect particles from the list
              if (prev==-1) cell->FirstCellParticle=next;
              else PIC::ParticleBuffer::SetNext(next,prev);

              if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

              PIC::ParticleBuffer::DeleteParticle(modelParticle);
              goto _StartParticleCollisionLoop_;
            }
          }
        }

      }

      //sample total collision frequentcy
      for (spec=0;spec<PIC::nTotalSpecies;spec++) if (TotalProspectiveCollisionParticleWeight[spec]>0.0) {
        *(spec+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=
          TotalOccurringCollisionParticleWeight[spec]/TotalProspectiveCollisionParticleWeight[spec]/node->block->GetLocalTimeStep(spec);
      }

    }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

    node=node->nextNodeThisThread;
  }

  //delete the temporary particle representing the background atmosphere
  PIC::ParticleBuffer::DeleteParticle(tempBackgroundAtmosphereParticle);
}


void PIC::MolecularCollisions::BackgroundAtmosphere::RemoveThermalBackgroundParticles() {

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;
  long int next,prev,ptr,centerNodeIndexCounter,LocalCellNumber;
  PIC::ParticleBuffer::byte *pdata;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;
    int i,j,k;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];


  //filter particles
  while (node!=NULL) {
    for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
      //sort the particle from the cell
      LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
      ptr=node->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;

      while (ptr!=-1) {
        pdata=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

        next=PIC::ParticleBuffer::GetNext(pdata);
        prev=PIC::ParticleBuffer::GetPrev(pdata);

        if (KeepBackgroundAtmosphereParticle(pdata)==false) {
          //the particle should be removed
          //reconnect particles from the list
          if (prev==-1) node->block->GetCenterNode(LocalCellNumber)->FirstCellParticle=next;
          else PIC::ParticleBuffer::SetNext(next,prev);

          if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

          PIC::ParticleBuffer::DeleteParticle(ptr);
        }

        ptr=next;
      }
    }

    node=node->nextNodeThisThread;
  }
}


#endif
