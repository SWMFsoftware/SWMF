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
long int PIC::BC::nTotalInjectedParticles=0;

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

//  nInjectedParticles=0;

  //model the particle injection through the face of the bounding box
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  for (nodeptr=boundingBoxInjectionBlocksList.begin(),end=boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
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
    default:
      exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
    }
  }

  //go through the blocks on the currect processor and execute the boundary condition procesure
  cInternalBoundaryConditionsDescriptor *bc;
  cInternalSphericalData *Sphere;
  cInternalCircleData *Circle;
  cInternalSphere1DData *Sphere1D;

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
          if (Sphere->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Sphere->InjectionBoundaryCondition((void*)Sphere);

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
          if (Circle->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Circle->InjectionBoundaryCondition((void*)Circle);

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
          if (Sphere1D->InjectionBoundaryCondition!=NULL) nTotalInjectedParticles+=Sphere1D->InjectionBoundaryCondition((void*)Sphere1D);

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


