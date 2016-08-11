//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that controls execution of the boundary conditions

#include "pic.h"


//the list of blocks where the injection BCs are applied

PIC::BC::fBlockInjectionIndicator PIC::BC::BlockInjectionBCindicatior=NULL;
PIC::BC::fBlockInjectionBC PIC::BC::userDefinedBoundingBlockInjectionFunction=NULL;
list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > PIC::BC::boundingBoxInjectionBlocksList;
long int* PIC::BC::nInjectedParticles=NULL;
double *PIC::BC::ParticleProductionRate=NULL;
double *PIC::BC::ParticleMassProductionRate=NULL;
long int PIC::BC::nTotalInjectedParticles=0;
PIC::BC::fUserDefinedParticleInjectionFunction PIC::BC::UserDefinedParticleInjectionFunction=NULL;

//the extra injection process by the exosphere model (src/models/exosphere)
PIC::BC::fExosphereModelExtraInjectionFunction PIC::BC::ExosphereModelExtraInjectionFunction=NULL;

//speces that are injected into the computational domain with PIC::BC::ExternalBoundary::OpenFlow::Inject()
bool PIC::BC::ExternalBoundary::OpenFlow::BoundaryInjectionFlag[PIC::nTotalSpecies];
list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > PIC::BC::ExternalBoundary::OpenFlow::InjectionBlocksList;
bool PIC::BC::ExternalBoundary::OpenFlow::InjectionBlocksListInitFlag=false;

//====================================================
//create the list of blocks where the injection BCs are applied
void PIC::BC::InitBoundingBoxInjectionBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  if (startNode==PIC::Mesh::mesh.rootTree) {
    if (boundingBoxInjectionBlocksList.size()!=0) exit(__LINE__,__FILE__,"Error: reinitialization of the 'boundingBoxInjectionBlocksList' list");
    if (BlockInjectionBCindicatior==NULL) exit(__LINE__,__FILE__,"Error: 'BlockInjectionBCindicatior' is not defined");
  }


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (BlockInjectionBCindicatior(startNode)==true) boundingBoxInjectionBlocksList.push_back(startNode);
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (int i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) InitBoundingBoxInjectionBlockList(downNode);
  }
}

//====================================================
//the function controls the overall execution of the injection boundary conditions
void PIC::BC::InjectionBoundaryConditions() {


  //inject particle from the boundary of the domain with the "Free flow" injection model
  if (_PIC_BC__OPEN_FLOW_INJECTION__MODE_ == _PIC_BC__OPEN_FLOW_INJECTION__MODE_ON_) {
    PIC::BC::ExternalBoundary::OpenFlow::Inject();
  }

  //model the particle injection through the face of the bounding box
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  if (userDefinedBoundingBlockInjectionFunction!=NULL) for (nodeptr=boundingBoxInjectionBlocksList.begin(),end=boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (node->Thread==PIC::Mesh::mesh.ThisThread) {
      nTotalInjectedParticles+=userDefinedBoundingBlockInjectionFunction(node);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
      EndTime=MPI_Wtime();
      node->ParallelLoadMeasure+=EndTime-StartTime;
      StartTime=EndTime;
#endif
    }
  }


  //particle injection from the internal surface
#ifndef _INTERNAL_BOUNDARY_MODE_
  exit(__LINE__,__FILE__,"Error: macroscopic variable '_INTERNAL_BOUNDARY_MODE_' is not defined. It souvle be defined in 'meshAMRdef.h'");
#endif

#if  _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  //reset the processing flag in the internal surface boundaries
  list<cInternalBoundaryConditionsDescriptor>::iterator descriptor;

  for (descriptor=PIC::Mesh::mesh.InternalBoundaryList.begin();descriptor!=PIC::Mesh::mesh.InternalBoundaryList.end();descriptor++) {
    switch (descriptor->BondaryType) {
    case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
      ((cInternalSphericalData*)(descriptor->BoundaryElement))->ProcessedBCflag=false;
      break;
    case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
      ((cInternalCircleData*)(descriptor->BoundaryElement))->ProcessedBCflag=false;
      break;
    case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
      ((cInternalSphere1DData*)(descriptor->BoundaryElement))->ProcessedBCflag=false;
      break;
    case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
      ((cInternalRotationBodyData*)(descriptor->BoundaryElement))->ProcessedBCflag=false;
      break;
    case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
      ((cInternalNastranSurfaceData*)(descriptor->BoundaryElement))->ProcessedBCflag=false;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
    }
  }

  //go through the blocks on the currect processor and execute the boundary condition procesure
  cInternalBoundaryConditionsDescriptor *bc;
  cInternalSphericalData *Sphere;
  cInternalCircleData *Circle;
  cInternalSphere1DData *Sphere1D;
  cInternalRotationBodyData *RotationBody;
  cInternalNastranSurfaceData *NastranSurface;

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) if ((bc=node->InternalBoundaryDescriptorList)!=NULL) {
    for (;bc!=NULL;bc=bc->nextInternalBCelement) {
      switch (bc->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
        Sphere=(cInternalSphericalData*)(bc->BoundaryElement);

        if (Sphere->ProcessedBCflag==false) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          StartTime=MPI_Wtime();
#endif

          Sphere->ProcessedBCflag=true;
          if (Sphere->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Sphere->InjectionBoundaryCondition(bc->BondaryType,(void*)Sphere);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
#endif
        }

        break;
      case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
        Circle=(cInternalCircleData*)(bc->BoundaryElement);

        if (Circle->ProcessedBCflag==false) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          StartTime=MPI_Wtime();
#endif

          Circle->ProcessedBCflag=true;
          if (Circle->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Circle->InjectionBoundaryCondition(bc->BondaryType,(void*)Circle);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
#endif
        }

        break;
      case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
        Sphere1D=(cInternalSphere1DData*)(bc->BoundaryElement);

        if (Sphere1D->ProcessedBCflag==false) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          StartTime=MPI_Wtime();
#endif

          Sphere1D->ProcessedBCflag=true;
          if (Sphere1D->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Sphere1D->InjectionBoundaryCondition(bc->BondaryType,(void*)Sphere1D);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
#endif
        }
        break;
      case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
        RotationBody=(cInternalRotationBodyData*)(bc->BoundaryElement);
	
        if (RotationBody->ProcessedBCflag==false) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          StartTime=MPI_Wtime();
#endif
	  
          RotationBody->ProcessedBCflag=true;
          if (RotationBody->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=RotationBody->InjectionBoundaryCondition(bc->BondaryType,(void*)RotationBody);
	  
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
#endif
	     }
	     break;

      case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
        NastranSurface=(cInternalNastranSurfaceData*)(bc->BoundaryElement);

        if (NastranSurface->ProcessedBCflag==false) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          StartTime=MPI_Wtime();
#endif

          NastranSurface->ProcessedBCflag=true;
          if (NastranSurface->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=NastranSurface->InjectionBoundaryCondition(bc->BondaryType,(void*)NastranSurface);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
          node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
#endif
       }
       break;

      default:
        exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
      }
    }
  }

#endif
}

//====================================================
//calculate the rate of injection of particles with the Maxwellian distribution
double PIC::BC::CalculateInjectionRate_MaxwellianDistribution(double NumberDesnity,double Temp,const double *BulkVelocity,double *ExternalNormal,int spec) {
  double res=0.0,vNorm,beta,vv,sc,cc,aa;
  int idim;

#if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if (spec<0) exit(__LINE__,__FILE__,"Error: negative value of the 'spec' variable");
#endif

  for (vv=0.0,vNorm=0.0,idim=0;idim<DIM;idim++) vv+=BulkVelocity[idim]*BulkVelocity[idim],vNorm+=BulkVelocity[idim]*ExternalNormal[idim];

  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*Temp));
  vv=sqrt(vv);
  sc=vv*beta;

  if ((sc>10.0)&&(vNorm<0.0)) res=fabs(vNorm)*NumberDesnity;
  else if ((sc>10.0)&&(vNorm>0.0)) res=0.0;
  else {
    if (vv>1.0E-6) {
      cc=-vNorm/vv;
      aa=(exp(-sc*sc*cc*cc)+sqrtPi*sc*cc*(1.0+erf(sc*cc)))/(2.0*sqrtPi);
    }
    else aa=1.0/(2.0*sqrtPi);

    res=fabs(aa)*NumberDesnity/beta;
  }

  return res;
}

//====================================================
//Free flow injection boundary conditions
void PIC::BC::ExternalBoundary::OpenFlow::InitBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  bool ExternalFaces[6];
  int nface;

  if ((startNode==PIC::Mesh::mesh.rootTree)&&(InjectionBlocksListInitFlag==true)) return;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
      for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
        InjectionBlocksList.push_back(startNode);
        break;
      }
    }
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (int i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) InitBlockList(downNode);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) InjectionBlocksListInitFlag=true;
}


int PIC::BC::ExternalBoundary::OpenFlow::InjectBlock(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,int nInjectionFace) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;
  double v[3];
  int nInjectionFaceBegin,nInjectionFaceEnd;

  //the total number of cells along each direction
  static const int iFaceIndex[6]={1,1, 0,0, 0,0};
  static const int jFaceIndex[6]={2,2, 2,2, 1,1};
  static const int iFaceMax[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  double ModelParticlesInjectionRate,SurfaceArea;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);

    if (nInjectionFace==-1) nInjectionFaceBegin=0,nInjectionFaceEnd=2*DIM;
    else nInjectionFaceBegin=nInjectionFace,nInjectionFaceEnd=nInjectionFace+1;

    for (nface=nInjectionFaceBegin;nface<nInjectionFaceEnd;nface++) if (ExternalFaces[nface]==true) for (int i=0;i<iFaceMax[iFaceIndex[nface]];i++) for (int j=0;j<iFaceMax[jFaceIndex[nface]];j++) {
      int iCell,jCell,kCell,nd;
      PIC::Mesh::cDataCenterNode *CenterNode;

      SurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      switch (nface) {
      case 0:
        iCell=0,jCell=i,kCell=j;
        SurfaceArea/=_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
        break;
      case 1:
        iCell=_BLOCK_CELLS_X_-1,jCell=i,kCell=j;
        SurfaceArea/=_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
        break;
      case 2:
        iCell=i,jCell=0,kCell=j;
        SurfaceArea/=_BLOCK_CELLS_X_*_BLOCK_CELLS_Z_;
        break;
      case 3:
        iCell=i,jCell=_BLOCK_CELLS_Y_-1,kCell=j;
        SurfaceArea/=_BLOCK_CELLS_X_*_BLOCK_CELLS_Z_;
        break;
      case 4:
        iCell=i,jCell=j,kCell=0;
        SurfaceArea/=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_;
        break;
      case 5:
        iCell=i,jCell=j,kCell=_BLOCK_CELLS_Z_-1;
        SurfaceArea/=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_;
        break;
      default:
        exit(__LINE__,__FILE__,"Oops.... the face counted is out of range");
      }

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) exit(__LINE__,__FILE__,"Error: the cell is not alocated");

      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      double BulkVelocity[3],Temp,NumberDensity;

      CenterNode->GetBulkVelocity(BulkVelocity,spec);
      CenterNode->GetTranslationalTemperature(&Temp,spec);
      NumberDensity=CenterNode->GetNumberDensity(spec);

      ModelParticlesInjectionRate=((NumberDensity!=0.0)&&(Temp!=0.0)) ? PIC::BC::CalculateInjectionRate_MaxwellianDistribution(NumberDensity,Temp,BulkVelocity,ExternalNormal,spec) : 0.0;

      if (ModelParticlesInjectionRate>0.0) {
        ModelParticlesInjectionRate*=SurfaceArea/ParticleWeight;

        PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

        //adjust the coordinate frame
        switch (nface) {
        case 0: case 1:
          for (idim=0;idim<3;idim++) {
            e0[idim]/=_BLOCK_CELLS_Y_,e1[idim]/=_BLOCK_CELLS_Z_;
            x0[idim]+=e0[idim]*i+e1[idim]*j;
          }

          break;
        case 2: case 3:
          for (idim=0;idim<3;idim++) {
            e0[idim]/=_BLOCK_CELLS_X_,e1[idim]/=_BLOCK_CELLS_Z_;
            x0[idim]+=e0[idim]*i+e1[idim]*j;
          }

          break;
        case 4: case 5:
          for (idim=0;idim<3;idim++) {
            e0[idim]/=_BLOCK_CELLS_X_,e1[idim]/=_BLOCK_CELLS_Y_;
            x0[idim]+=e0[idim]*i+e1[idim]*j;
          }

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the face counter is out of range");
        }


        //inject new particles into the domain
        while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
          //generate the new particle position on the face
          for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

          //generate a particle
          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
          nInjectedParticles++;

          //sample the source rate
          PIC::BC::nInjectedParticles[spec]++;
          PIC::BC::ParticleProductionRate[spec]+=ParticleWeight/LocalTimeStep;
          PIC::BC::ParticleMassProductionRate[spec]+=ParticleWeight*PIC::MolecularData::GetMass(spec)/LocalTimeStep;

          //generate particles' velocity
          PIC::Distribution::InjectMaxwellianDistribution(v,BulkVelocity,Temp,ExternalNormal,spec,-1);

          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

          //inject the particle into the system
          _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
        }
      }
    }
  }

  return nInjectedParticles;
}

void PIC::BC::ExternalBoundary::OpenFlow::Inject() {
  int spec,nInjectedParticles=0;

  if (InjectionBlocksListInitFlag==false) InitBlockList();

  //model the particle injection through the face of the bounding box
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

  #if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
  #endif

  for (nodeptr=InjectionBlocksList.begin(),end=InjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (node->Thread==PIC::Mesh::mesh.ThisThread) for (spec=0;spec<PIC::nTotalSpecies;spec++) if (BoundaryInjectionFlag[spec]==true) {
      nInjectedParticles+=InjectBlock(spec,node);

      //account for the source rate
      double LocalTimeStep=node->block->GetLocalTimeStep(spec);

      PIC::BC::nTotalInjectedParticles+=nInjectedParticles;
      PIC::BC::nInjectedParticles[spec]+=nInjectedParticles;

      PIC::BC::ParticleProductionRate[spec]+=nInjectedParticles/LocalTimeStep;
      PIC::BC::ParticleMassProductionRate[spec]+=nInjectedParticles/LocalTimeStep*PIC::MolecularData::GetMass(spec);

      #if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
      EndTime=MPI_Wtime();
      node->ParallelLoadMeasure+=EndTime-StartTime;
      StartTime=EndTime;
      #endif
    }
  }
}

//=========================================================================================================
//get the direction of the gas flow at the boundary of the computational domain
int PIC::BC::ExternalBoundary::ExternalBoundaryFlowDirection(int spec, int nface, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,nd,idim;
  double AveragedVelocity[3]={0.0,0.0,0.0};
  double BulkVelocity[3],ExternalNormal[3],c;

  //the total number of cells along each direction
  static const int iFaceIndex[6]={1,1, 0,0, 0,0};
  static const int jFaceIndex[6]={2,2, 2,2, 1,1};
  static const int iFaceMax[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  for (int i=0;i<iFaceMax[iFaceIndex[nface]];i++) for (int j=0;j<iFaceMax[jFaceIndex[nface]];j++) {
    int iCell,jCell,kCell,nd;
    PIC::Mesh::cDataCenterNode *CenterNode;

    switch (nface) {
    case 0:
      iCell=0,jCell=i,kCell=j;
      break;
    case 1:
      iCell=_BLOCK_CELLS_X_-1,jCell=i,kCell=j;
      break;
    case 2:
      iCell=i,jCell=0,kCell=j;
      break;
    case 3:
      iCell=i,jCell=_BLOCK_CELLS_Y_-1,kCell=j;
      break;
    case 4:
      iCell=i,jCell=j,kCell=0;
      break;
    case 5:
      iCell=i,jCell=j,kCell=_BLOCK_CELLS_Z_-1;
      break;
    default:
      exit(__LINE__,__FILE__,"Oops....");
    }

    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
    if ((CenterNode=node->block->GetCenterNode(nd))==NULL) exit(__LINE__,__FILE__,"Error: the cell is not alocated");

    CenterNode->GetBulkVelocity(BulkVelocity,spec);
    for (idim=0;idim<3;idim++) AveragedVelocity[idim]+=BulkVelocity[idim];
  }

  node->GetExternalNormal(ExternalNormal,nface);
  for (idim=0,c=0.0;idim<3;idim++) c+=AveragedVelocity[idim]*ExternalNormal[idim];

  int res=_PIC__EXTERNAL_BOUNDARY_FLOW_DIRECTION__UNDEFINED_;

  if (c>0.0) res=_PIC__EXTERNAL_BOUNDARY_FLOW_DIRECTION__OUTWARD_;
  if (c<0.0) res=_PIC__EXTERNAL_BOUNDARY_FLOW_DIRECTION__INWARD_;

  return res;
}












