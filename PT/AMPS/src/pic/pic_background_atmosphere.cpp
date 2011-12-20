//====================================================
//$Id$
//====================================================
//the functions that model colecular collisions with the background atmosphere

#include "pic.h"

#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_

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

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__UNIFORM_
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


//==================== DEBUG ===================
/*
        double tmp,tmp1,tmp2,xtmp[3];

        cell->GetX(xtmp);

        tmp=GetCellMaximumBackgroundNumberDensity(BackgroundSpecieNumber,cell,node);
        tmp1=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_CO2(xtmp);
        tmp2=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_O(xtmp);
*/
//==================== END DEBUG ===============




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

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__UNIFORM_
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
