//====================================================
//$Id$
//====================================================
//the functions that controls the particle motion

#include "pic.h"


PIC::Mover::fSpeciesDependentParticleMover *PIC::Mover::MoveParticleTimeStep=NULL;

//====================================================
//init the particle mover
void PIC::Mover::Init() {
  int s;

  //allocate the mover functions array
  if ((MoveParticleTimeStep!=NULL)||(PIC::nTotalSpecies==0)) exit(__LINE__,__FILE__,"Error: the initialization of PIC::Mover is failed");

  MoveParticleTimeStep=new fSpeciesDependentParticleMover[PIC::nTotalSpecies];
  for (s=0;s<PIC::nTotalSpecies;s++) MoveParticleTimeStep[s]=NULL;
}

//====================================================
//move all existing particles
void PIC::Mover::MoveParticles() {
  int s,i,j,k;
  long int LocalCellNumber,ParticleList,ptr;
  double LocalTimeStep;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  //move existing particles
  while (node!=NULL) {

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            ParticleList=node->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;

            while (ParticleList!=-1) {
              ptr=ParticleList;
              ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);
              s=PIC::ParticleBuffer::GetI(ptr);
              LocalTimeStep=node->block->GetLocalTimeStep(s);

              MoveParticleTimeStep[s](ptr,LocalTimeStep,node);
            }

          }
       }
    }

    node=node->nextNodeThisThread;
  }


  //update the particle lists (The local processor)
  PIC::Mesh::cDataCenterNode *cell;

  /*
  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  while (node!=NULL) {
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            cell=node->block->GetCenterNode(LocalCellNumber);

            cell->FirstCellParticle=cell->tempParticleMovingList;
            cell->tempParticleMovingList=-1;
          }
       }
    }

    node=node->nextNodeThisThread;
  }
  */

  //update the particle lists (The connecting and local processors)
  for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
    node=(thread==PIC::Mesh::mesh.ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    if (node==NULL) continue;

    while (node!=NULL) {
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
              cell=node->block->GetCenterNode(LocalCellNumber);

              cell->FirstCellParticle=cell->tempParticleMovingList;
              cell->tempParticleMovingList=-1;
            }
         }
      }

      node=node->nextNodeThisThread;
    }
  }

}


//====================================================
//not forces, constant time step, constant particle weight
int PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  double dtMin=-1.0,dtTemp;
  PIC::ParticleBuffer::byte *ParticleData;
  double *v,*x,*xminBlock,*xmaxBlock;
  int spec,iNeibNode[3],neibNodeDirection,idim,nface_dtMin=-1;




//########## DEBUG ############

static long int nCallCounter=0;

nCallCounter++;

//########## END DEBUG ##########



  //the descriptors of the internal surfaces
  cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*InternalBoundaryDescriptor_dtMin=NULL,*lastInternalBoundaryDescriptor=NULL;

  //spherical internal surface
  cInternalSphericalData *Sphere;
  double radiusSphere,*x0Sphere;
  double a,b,c,d,dx,dy,dz,sqrt_d,dt1;

#define _UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_        0
#define _BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_       1
#define _INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_  2

  int ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;


  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  v=PIC::ParticleBuffer::GetV(ParticleData);
  x=PIC::ParticleBuffer::GetX(ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);



//#######  DEBUG #########
double xpinit[3];
for (idim=0;idim<3;idim++) xpinit[idim]=x[idim];
//####### END DEBUG #########

  while (dtTotal>0.0) {
    xminBlock=startNode->xmin;
    xmaxBlock=startNode->xmax;

    dtMin=dtTotal,nface_dtMin=-1;
    InternalBoundaryDescriptor_dtMin=NULL;
    ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;


    //Calculate the time of flight to the nearest block's face
#if DIM == 1
    exit(__LINE__,__FILE__,"not implemented");
#elif DIM == 2
    exit(__LINE__,__FILE__,"not implemented");
#elif DIM == 3

    for (idim=0;idim<DIM;idim++) if (v[idim]!=0.0) {
      //nface=0,2,4
      dtTemp=(xminBlock[idim]-x[idim])/v[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=2*idim,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

      //nface=1,3,5
      dtTemp=(xmaxBlock[idim]-x[idim])/v[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1+2*idim,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

    }

#else
    exit(__LINE__,__FILE__,"Error: unknown value of DIM");
#endif


    //Calculate the time of flight to the nearest internal surface
    for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {
      if (InternalBoundaryDescriptor==lastInternalBoundaryDescriptor) continue;

      switch (InternalBoundaryDescriptor->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
        Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

        dx=x[0]-x0Sphere[0],dy=x[1]-x0Sphere[1],dz=x[2]-x0Sphere[2];
        a=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
        b=2.0*(v[0]*dx+v[1]*dy+v[2]*dz);
        c=dx*dx+dy*dy+dz*dz-radiusSphere*radiusSphere;

        d=b*b-4.0*a*c;

        if (d<0.0) {
          if (4.0*a*pow(PIC::Mesh::mesh.EPS,2)>-d) d=0.0; //equvalent to |EPS/particle speed| > sqrt(|d|)/(2a) -> the particle is within the distance of EPS from the surface's boundary
          else continue;
        }

        if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

        sqrt_d=sqrt(d);
        dtTemp=-(b+sqrt_d)/(2.0*a);
        dt1=-2.0*c/(b+sqrt_d);
        if ((dtTemp<0.0)||((dt1>0.0)&&(dt1<dtTemp))) dtTemp=dt1;

        if ((0.0<dtTemp)&&(dtTemp<dtMin)) {
          dtMin=dtTemp,InternalBoundaryDescriptor_dtMin=InternalBoundaryDescriptor;
          ParticleIntersectionCode=_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_;
        }

        break;
      default:
        exit(__LINE__,__FILE__,"Error: undetermined internal boundary type");
      }
    }


    //advance the particle's position
    for (idim=0;idim<DIM;idim++) x[idim]+=dtMin*v[idim];

    //adjust the particle moving time
    dtTotal-=dtMin;

    //interaction with the faces of the block and internal surfaces
    if (ParticleIntersectionCode==_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      lastInternalBoundaryDescriptor=InternalBoundaryDescriptor_dtMin;
      PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction(spec,ptr,x,v,dtTotal,startNode,InternalBoundaryDescriptor_dtMin);
    }
    else if (ParticleIntersectionCode==_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      iNeibNode[0]=0,iNeibNode[1]=0,iNeibNode[2]=0;
      neibNodeDirection=(int)(nface_dtMin/2);
      iNeibNode[neibNodeDirection]=(nface_dtMin-2*neibNodeDirection==0) ? -1 : 1;


      double xProbe[3]={x[0],x[1],x[2]};
      xProbe[neibNodeDirection]+=(startNode->xmax[neibNodeDirection]-startNode->xmin[neibNodeDirection])*iNeibNode[neibNodeDirection]/100.0;
      newNode=PIC::Mesh::mesh.findTreeNode(xProbe,startNode);


//      newNode=PIC::Mesh::mesh.getNeibNode(iNeibNode[0],iNeibNode[1],iNeibNode[2],startNode);





      if (newNode==NULL) {
        //the particle left the computational domain
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

      xminBlock=newNode->xmin;
      xmaxBlock=newNode->xmax;

      #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
      //check if the new particle coordiname is within the new block
      for (idim=0;idim<DIM;idim++) if ((x[idim]<xminBlock[idim]-PIC::Mesh::mesh.EPS)||(x[idim]>xmaxBlock[idim]+PIC::Mesh::mesh.EPS)) {
        exit(__LINE__,__FILE__,"Error: the new particles' coordinates are outside of the block");
      }
      #endif

      //move the particle's position exactly to the block's face
      for (idim=0;idim<DIM;idim++) {
        if (x[idim]<=xminBlock[idim]) x[idim]=xminBlock[idim]*(1.0+1.0E-8);
        if (x[idim]>=xmaxBlock[idim]) x[idim]=xmaxBlock[idim]*(1.0-1.0E-8);
      }

      x[neibNodeDirection]=(iNeibNode[neibNodeDirection]==-1) ? xmaxBlock[neibNodeDirection] : xminBlock[neibNodeDirection];

      //reserve the place for particle's cloning

      //adjust the value of 'startNode'
      startNode=newNode;
    }



  }




  //the particle is still within the computational domain:
  //place it to the local list of particles related to the new block and cell
  long int LocalCellNumber;
  int i,j,k;

  PIC::Mesh::cDataCenterNode *cell;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=startNode->block->GetCenterNode(LocalCellNumber);


  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;

  return _PARTICLE_MOTION_FINISHED_;
}

int PIC::Mover::UniformWeight_UniformTimeStep_noForce(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  PIC::ParticleBuffer::byte *ParticleData;
  double v[3],x[3];
  int idim;

  #if  _SIMULATION_TIME_STEP_MODE_ ==  _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #endif



 //################  DEBUG ###################
//  return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt,startNode);

//################# END DEBUG ################



  //Check if the startNode has cut cells
  #if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
  if (startNode->InternalBoundaryDescriptorList!=NULL) {
     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt,startNode);
  }
  #endif



//######## DEBUG ########

/*
if (ptr==1212) {
  std::cout << __LINE__ << endl;
}
*/

//######## END DEBUG ####



  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(v,ParticleData);
  PIC::ParticleBuffer::GetX(x,ParticleData);

  //advance the particle positions
  #if DIM == 3
  for (idim=0;idim<DIM;idim++) x[idim]+=dt*v[idim];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //determine the new particle location
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  newNode=PIC::Mesh::mesh.findTreeNode(x,startNode);


  if (newNode==NULL) {
    //the particle left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  //Check if the newNode has cut cells
  #if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
  if (newNode->InternalBoundaryDescriptorList!=NULL) {
     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt,startNode);
  }
  #endif

  //the particle is still within the computational domain:
  //place it to the local list of particles related to the new block and cell
  long int LocalCellNumber;
  int i,j,k;

  PIC::Mesh::cDataCenterNode *cell;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,newNode))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=newNode->block->GetCenterNode(LocalCellNumber);

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetX(x,ParticleData);

  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;

  return _PARTICLE_MOTION_FINISHED_;
}
