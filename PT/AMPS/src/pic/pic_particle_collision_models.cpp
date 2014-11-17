//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//modeling of collision between the model particles

#include "pic.h"

//sampling offset of the collision frequentcy
int PIC::MolecularCollisions::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset=-1;


//request the sampling data
int PIC::MolecularCollisions::ParticleCollisionModel::RequestSamplingData(int offset) {
  int SamplingLength=0;

#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
  CollsionFrequentcySampling::SamplingBufferOffset=offset+SamplingLength;
  SamplingLength+=sizeof(double)*(PIC::nTotalSpecies*PIC::nTotalSpecies);
#endif

  return SamplingLength;
}

void PIC::MolecularCollisions::ParticleCollisionModel::PrintVariableList(FILE* fout,int DataSetNumber) {
  int s;

  for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,", \"Coll Freq(%s) [m^{-3} s^{-1}]\"",PIC::MolecularData::GetChemSymbol(s));
}

void PIC::MolecularCollisions::ParticleCollisionModel::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode)  {
  double CollisionFrequentcy[PIC::nTotalSpecies*PIC::nTotalSpecies];
  int i,s;
  double *SamplingBuffer;

  for (i=0;i<PIC::nTotalSpecies*PIC::nTotalSpecies;i++) CollisionFrequentcy[i]=0.0;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);

    for (s=0;s<PIC::nTotalSpecies*PIC::nTotalSpecies;s++) CollisionFrequentcy[s]+=SamplingBuffer[s]*InterpolationCoeficients[i];
  }

  SamplingBuffer=(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);
  for (s=0;s<PIC::nTotalSpecies*PIC::nTotalSpecies;s++) SamplingBuffer[s]=CollisionFrequentcy[s]/((PIC::LastSampleLength!=0) ? PIC::LastSampleLength : 1);
}

void PIC::MolecularCollisions::ParticleCollisionModel::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double t;
  double *SamplingBuffer=(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);
  int s;

  for (s=0;s<PIC::nTotalSpecies;s++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(SamplingBuffer+CollsionFrequentcySampling::Offset(DataSetNumber,s));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}

void PIC::MolecularCollisions::ParticleCollisionModel::Init() {
#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
#endif
}


//Non-time Counter collision procedure
void PIC::MolecularCollisions::ParticleCollisionModel::ntc() {
  int s,s0,s1,i,j,k;

  const int SigmaCrMax_nTest=100;
  const double SigmaCrMax_SafetyMargin=1.3;

  struct cParticleDataList {
    double vel[3];
    bool ValueChangedFlag;
    PIC::ParticleBuffer::byte *ParticleData;
  };


  int nParticleNumber[PIC::nTotalSpecies],nMaxSpecParticleNumber,cnt;
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],FirstCellParticle,ptr;
  double cellMeasure;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataBlockAMR *block;

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  //particle lists
  int ParticleDataListLength=100;
  cParticleDataList *s0ParticleDataList=new cParticleDataList[ParticleDataListLength];
  cParticleDataList *s1ParticleDataList=new cParticleDataList[ParticleDataListLength];  //s0ParticleDataList,s1ParticleDataList are used to store the actual data
  cParticleDataList *s0List=NULL,*s1List=NULL; //s0List and s1List are used to access s0ParticleDataList and s1ParticleDataList

  //sampling data
  PIC::Mesh::cDataCenterNode *cell;

#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
  char *SamplingData;
#endif


  //simulate particle's collisions
  while (node!=NULL) {
    block=node->block;
    memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            FirstCellParticle=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
            cellMeasure=cell->Measure;

#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
            SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
#endif

            if (FirstCellParticle!=-1) {
              //simulate collision in the cell
              //1. coalcualte the number of the model particles in the cell and re-allocate the particle lists if needed
              for (s=0;s<PIC::nTotalSpecies;s++) nParticleNumber[s]=0;
              for (ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) nParticleNumber[PIC::ParticleBuffer::GetI(ptr)]++;
              for (nMaxSpecParticleNumber=0,s=0;s<PIC::nTotalSpecies;s++) if (nMaxSpecParticleNumber<nParticleNumber[s]) nMaxSpecParticleNumber=nParticleNumber[s];

              //re-allocate the particle list if needed
              if (nMaxSpecParticleNumber>ParticleDataListLength) {
                delete [] s1ParticleDataList;
                delete [] s0ParticleDataList;

                ParticleDataListLength=std::max(1.2*nMaxSpecParticleNumber,2.0*ParticleDataListLength);
                s0ParticleDataList=new cParticleDataList[ParticleDataListLength];
                s1ParticleDataList=new cParticleDataList[ParticleDataListLength];
              }


              //loop through the first species
              double m0,m1,am;
              double LocalParticleWeight_s0,LocalParticleWeight_s1;
              double LocalTimeStep_s0,LocalTimeStep_s1;
              double minParticleWeightCorrection_s0,minParticleWeightCorrection_s1,sumWeightCorrection_s0,sumWeightCorrection_s1;
              PIC::ParticleBuffer::byte *ParticleData;

              for (s0=0;s0<PIC::nTotalSpecies;s0++) if (nParticleNumber[s0]!=0) {
                m0=PIC::MolecularData::GetMass(s0);
                LocalParticleWeight_s0=block->GetLocalParticleWeight(s0);
                LocalTimeStep_s0=PIC::ParticleWeightTimeStep::LocalTimeStep(s0,node);
                sumWeightCorrection_s0=0.0;

                //populate the particle list
                s0List=s0ParticleDataList;

                for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s0) {
                  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

                  PIC::ParticleBuffer::GetV(s0ParticleDataList[cnt].vel,ParticleData);
                  s0ParticleDataList[cnt].ParticleData=ParticleData;
                  s0ParticleDataList[cnt].ValueChangedFlag=false;

                  #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
                  double wc=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

                  if (cnt==0) minParticleWeightCorrection_s0=wc;
                  else if (wc<minParticleWeightCorrection_s0) minParticleWeightCorrection_s0=wc;

                  sumWeightCorrection_s0+=wc;
                  #else
                  sumWeightCorrection_s0+=1.0;
                  #endif

                  cnt++;
                }

                //loop through the second speices
                for (s1=s0;s1<PIC::nTotalSpecies;s1++) if (nParticleNumber[s1]!=0) {
                  m1=PIC::MolecularData::GetMass(s1);
                  am=m0+m1;
                  LocalParticleWeight_s1=block->GetLocalParticleWeight(s1);
                  LocalTimeStep_s1=PIC::ParticleWeightTimeStep::LocalTimeStep(s1,node);
                  sumWeightCorrection_s1=0.0;

                  //populate the list
                  if (s0==s1) s1List=s0List,minParticleWeightCorrection_s1=minParticleWeightCorrection_s0,sumWeightCorrection_s1=sumWeightCorrection_s0;
                  else {
                    s1List=s1ParticleDataList;

                    for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s1) {
                      ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

                      PIC::ParticleBuffer::GetV(s1ParticleDataList[cnt].vel,ParticleData);
                      s1ParticleDataList[cnt].ParticleData=ParticleData;
                      s1ParticleDataList[cnt].ValueChangedFlag=false;

                      #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
                      double wc=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

                      if (cnt==0) minParticleWeightCorrection_s1=wc;
                      else if (wc<minParticleWeightCorrection_s1) minParticleWeightCorrection_s1=wc;

                      sumWeightCorrection_s1+=wc;
                      #else
                      sumWeightCorrection_s1+=1.0;
                      #endif

                      cnt++;
                    }
                  }

                  //simulate collsions between the pair of species
                  double SigmaCrMax=0.0,SigmaCr,ancoll;
                  long int ncoll;
                  double v0[3],v1[3],cr;

                  //1.Evaluate the maximum value of SigmaCr
                  for (int ntest=0;ntest<SigmaCrMax_nTest;ntest++) {
                    memcpy(v0,s0List[(int)(rnd()*nParticleNumber[s0])].vel,3*sizeof(double));
                    memcpy(v1,s1List[(int)(rnd()*nParticleNumber[s1])].vel,3*sizeof(double));
                    cr=sqrt(pow(v1[0]-v0[0],2)+pow(v1[1]-v0[1],2)+pow(v1[2]-v0[2],2));

                    #if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__HS_
                    SigmaCr=cr*PIC::MolecularData::MolecularModels::HS::GetTotalCrossSection(s0,s1);
                    #elif _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__USER_DEFINED_
                    SigmaCr=cr*PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::GetTotalCrossSection(v0,s0,v1,s1,block,cell);
                    #else
                    exit(__LINE__,__FILE__,"Error: not implemented");
                    #endif

                    if (SigmaCr>SigmaCrMax) SigmaCrMax=SigmaCr;
                  }

                  SigmaCrMax*=SigmaCrMax_SafetyMargin;

                  //2.Evaluate the prospective number of collisions
                  double maxLocalTimeStep,minLocalParticleWeight;

                  maxLocalTimeStep=max(LocalTimeStep_s0,LocalTimeStep_s1);
                  minLocalParticleWeight=min(LocalParticleWeight_s0*minParticleWeightCorrection_s0,LocalParticleWeight_s1*minParticleWeightCorrection_s1);

                  if (s0==s1) ancoll=0.5*sumWeightCorrection_s0*(sumWeightCorrection_s0-1.0)*LocalParticleWeight_s0*SigmaCrMax*LocalTimeStep_s0/cellMeasure;
                  else ancoll=(sumWeightCorrection_s0*LocalParticleWeight_s0)*(sumWeightCorrection_s1*LocalParticleWeight_s1)*
                      SigmaCrMax*maxLocalTimeStep/minLocalParticleWeight/cellMeasure;

                  ncoll=(long int)ancoll;
                  ancoll-=ncoll;
                  if (rnd()<ancoll) ncoll++;


                  //3. Collision Limiting


                  //4. simulate collisions
                  int s0ptr,s1ptr,idim;
                  double vrel[3],vcm[3];

                  while (ncoll-->0) {
                    s0ptr=(int)((int)(rnd()*nParticleNumber[s0]));

                    if (s0!=s1) s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
                    else {
                      do {
                        s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
                      }
                      while (s0ptr==s1ptr);
                    }

                    memcpy(v0,s0List[s0ptr].vel,3*sizeof(double));
                    memcpy(v1,s1List[s1ptr].vel,3*sizeof(double));

                    for (cr=0.0,idim=0;idim<3;idim++) {
                      vrel[idim]=v1[idim]-v0[idim];
                      vcm[idim]=(m1*v1[idim]+m0*v0[idim])/(m1+m0);
                      cr+=vrel[idim]*vrel[idim];
                    }

                    cr=sqrt(cr);

                    #if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__HS_
                    SigmaCr=cr*PIC::MolecularData::MolecularModels::HS::GetTotalCrossSection(s0,s1);
                    #elif _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__USER_DEFINED_
                    SigmaCr=cr*PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::GetTotalCrossSection(v0,s0,v1,s1,block,cell);
                    #else
                    exit(__LINE__,__FILE__,"Error: not implemented");
                    #endif

                    if (rnd()*SigmaCrMax>=SigmaCr) continue;

                    //determine weather the properties of particles were updated
                    double pUpdate_s0,pUpdate_s1;
                    bool UpdateFlag[2];

                    pUpdate_s0=minLocalParticleWeight/LocalParticleWeight_s0 * LocalTimeStep_s0/maxLocalTimeStep;
                    pUpdate_s1=minLocalParticleWeight/LocalParticleWeight_s1 * LocalTimeStep_s1/maxLocalTimeStep;

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
                    pUpdate_s0/=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData);
                    pUpdate_s1/=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData);
#endif

                    UpdateFlag[0]=(rnd()<pUpdate_s0) ? true :false;
                    UpdateFlag[1]=(rnd()<pUpdate_s1) ? true :false;

                    //model the internal energy exchange
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
                    double crInit=cr;

                    PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData,cr,UpdateFlag,cell);
                    for (int idim=0;idim<3;idim++) vrel[idim]=(cr>1.0E-10) ? vrel[idim]*cr/crInit : 0.0;
#endif



                    //the collision is considered to be true
                    PIC::MolecularCollisions::VelocityScattering::HS::VelocityAfterCollision(vrel,s0,s1);


                    //calcualte the new value of the particle's velocities
                    for (int idim=0;idim<3;idim++) {
                      v1[idim]=vcm[idim]+m0/am*vrel[idim];
                      v0[idim]=vcm[idim]-m1/am*vrel[idim];
                    }


                    //update the velocities in the lists
                    if (UpdateFlag[0]==true) {
                      s0List[s0ptr].ValueChangedFlag=true;
                      memcpy(s0List[s0ptr].vel,v0,3*sizeof(double));

#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
                      int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s0,s1);

                      *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData)*LocalParticleWeight_s0/LocalTimeStep_s0/cellMeasure;
#endif
                    }

                    if (UpdateFlag[1]==true) {
                      s1List[s1ptr].ValueChangedFlag=true;
                      memcpy(s1List[s1ptr].vel,v1,3*sizeof(double));

#if _PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_
                      int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s1,s0);

                      *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData)*LocalParticleWeight_s1/LocalTimeStep_s1/cellMeasure;
#endif
                    }
                  }

                  //update the velocities of the species 's1'
                  if (s0!=s1) {
                    for (cnt=0;cnt<nParticleNumber[s1];cnt++) if (s1List[cnt].ValueChangedFlag==true) {
                      PIC::ParticleBuffer::SetV(s1List[cnt].vel,s1List[cnt].ParticleData);
                    }
                  }
                }

                //update velocities of species 's0'
                for (cnt=0;cnt<nParticleNumber[s0];cnt++) if (s0List[cnt].ValueChangedFlag==true) {
                  PIC::ParticleBuffer::SetV(s0List[cnt].vel,s0List[cnt].ParticleData);
                }
              }

            } //end simulation of collisons in a cell


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

  delete [] s1ParticleDataList;
  delete [] s0ParticleDataList;
}
