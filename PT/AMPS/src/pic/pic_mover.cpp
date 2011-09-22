//====================================================
//$Id$
//====================================================
//the functions that controls the particle motion

#include "pic.h"


PIC::Mover::fSpeciesDependentParticleMover *PIC::Mover::MoveParticleTimeStep=NULL;
PIC::Mover::fTotalParticleAcceleration PIC::Mover::TotalParticleAcceleration=PIC::Mover::TotalParticleAcceleration_default;
PIC::Mover::fSpeciesDependentParticleMover_BoundaryInjection *PIC::Mover::MoveParticleBoundaryInjection=NULL;

//====================================================
//init the particle mover
void PIC::Mover::Init() {
  int s;

  //allocate the mover functions array
  if ((MoveParticleTimeStep!=NULL)||(PIC::nTotalSpecies==0)) exit(__LINE__,__FILE__,"Error: the initialization of PIC::Mover is failed");

  MoveParticleTimeStep=new fSpeciesDependentParticleMover[PIC::nTotalSpecies];
  MoveParticleBoundaryInjection=new fSpeciesDependentParticleMover_BoundaryInjection[PIC::nTotalSpecies];
  for (s=0;s<PIC::nTotalSpecies;s++) MoveParticleTimeStep[s]=NULL,MoveParticleBoundaryInjection[s]=NULL;
}

//====================================================
//the default function for the particle acceleration
void PIC::Mover::TotalParticleAcceleration_default(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  for (int idim=0;idim<3;idim++) accl[idim]=0.0;
}
//====================================================
//move all existing particles
void PIC::Mover::MoveParticles() {
  int s,i,j,k;
  long int LocalCellNumber,ParticleList,ptr;
  double LocalTimeStep;


  /*
  //the table of increments for accessing the cells in the block

  static int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  static bool initTableFlag=false;
  static int centerNodeIndexIncrement[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nStartCenterNodeIndex=-1,nTotalCenterNodes=-1;
  int centerNodeIndexCounter;

  if (initTableFlag==false) {
    int previousIndex=-1,newIndex;

    initTableFlag=true,nTotalCenterNodes=0,nStartCenterNodeIndex=PIC::Mesh::mesh.getCenterNodeLocalNumber(0,0,0);

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {

      centerNodeIndexTable[nTotalCenterNodes]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

      if (previousIndex==-1) {
        previousIndex=nStartCenterNodeIndex;
      }
      else {
        newIndex=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        centerNodeIndexIncrement[nTotalCenterNodes-1]=newIndex-previousIndex;

        previousIndex=newIndex;
      }

      nTotalCenterNodes++;
    }
  }
  */



  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  int centerNodeIndexCounter;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
  }



  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  //move existing particles
  while (node!=NULL) {

    /*
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      */
    {
      {
//        for (LocalCellNumber=nStartCenterNodeIndex,centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;LocalCellNumber+=centerNodeIndexIncrement[centerNodeIndexCounter++]) {

        for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {

            LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
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

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

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

      /*
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

              */
      {
        {
//          for (LocalCellNumber=nStartCenterNodeIndex,centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;LocalCellNumber+=centerNodeIndexIncrement[centerNodeIndexCounter++]) {

          for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {

              LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
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
int PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  double dtMin=-1.0,dtTemp;
  PIC::ParticleBuffer::byte *ParticleData;
  double *v,*x,*xminBlock,*xmaxBlock;
  int iNeibNode[3],neibNodeDirection,idim,nface_dtMin=-1;

  long int LocalCellNumber;
  int i,j,k;

  PIC::Mesh::cDataCenterNode *cell;
  bool MovingTimeFinished=false;


#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  double xInit[3];
#endif
  //=====================  DEBUG =========================




//########## DEBUG ############

static long int nCallCounter=0;

nCallCounter++;

if ((nCallCounter==1057317)&&(PIC::Mesh::mesh.ThisThread==5)) {
  cout << __FILE__ << "@" << __LINE__ << endl;
}

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
//  spec=PIC::ParticleBuffer::GetI(ParticleData);



//#######  DEBUG #########
//double xpinit[3];
//for (idim=0;idim<3;idim++) xpinit[idim]=x[idim];

//=====================  DEBUG =========================


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=startNode->block->GetCenterNode(LocalCellNumber);


  if (cell->Measure<=0.0) {
    cout << __FILE__<< __LINE__ << endl;
    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif


/*
if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1737.0E3) {
  cout << __FILE__ << __LINE__ << endl;
}

*/

  if ((nCallCounter==1057317)&&(PIC::Mesh::mesh.ThisThread==5)) {
    cout << __FILE__ << "@" << __LINE__ << endl;
  }
//===================== END DEBUG ==================



//####### END DEBUG #########

//  while (dtTotal>0.0) {
  while (MovingTimeFinished==false) {
MovingLoop:

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    //check the consistency of the particle mover
    int iTemp,jTemp,kTemp;

    if (PIC::Mesh::mesh.fingCellIndex(x,iTemp,jTemp,kTemp,startNode,false)==-1) {
      exit(__LINE__,__FILE__,"Error: the cell is not found");
    }

    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif



    xminBlock=startNode->xmin;
    xmaxBlock=startNode->xmax;

    MovingTimeFinished=true;
    dtMin=dtTotal,nface_dtMin=-1;
    InternalBoundaryDescriptor_dtMin=NULL;
    ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

#if _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ == _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ON_
    double accl[3],vv;
    int spec;

    spec=PIC::ParticleBuffer::GetI(ParticleData);
    TotalParticleAcceleration(accl,spec,ptr,x,v,startNode);

    a=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    vv=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (a*dtMin>0.2*vv) dtMin=0.2*vv/a;
#endif


    //Calculate the time of flight to the nearest block's face
#if DIM == 1
    exit(__LINE__,__FILE__,"not implemented");
#elif DIM == 2
    exit(__LINE__,__FILE__,"not implemented");
#elif DIM == 3

    for (idim=0;idim<DIM;idim++) if (fabs(v[idim])>0.0) {
      //nface=0,2,4
      dtTemp=(xminBlock[idim]-x[idim])/v[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=2*idim,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      //nface=1,3,5
      dtTemp=(xmaxBlock[idim]-x[idim])/v[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1+2*idim,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

    }

#else
    exit(__LINE__,__FILE__,"Error: unknown value of DIM");
#endif


    //Calculate the time of flight to the nearest internal surface
    if (FirstBoundaryFlag==false) for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {
      if (InternalBoundaryDescriptor==lastInternalBoundaryDescriptor) continue;

      switch (InternalBoundaryDescriptor->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
        Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

        dx=x[0]-x0Sphere[0],dy=x[1]-x0Sphere[1],dz=x[2]-x0Sphere[2];
        a=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
        b=2.0*(v[0]*dx+v[1]*dy+v[2]*dz);
        c=dx*dx+dy*dy+dz*dz-radiusSphere*radiusSphere;

        if (c<0.0) {
          /*
          double r=sqrt(dx*dx+dy*dy+dz*dz);

          if (radiusSphere-r>PIC::Mesh::mesh.EPS) {
            cout << "r=" << r << ",  Rsphere=" << radiusSphere << endl;
            exit(__LINE__,__FILE__,"Error: the particle inside the sphere");
          }
          else c=0.0;
          */

          //the particle is inside the sphese
          //1. project the particle on the surface of the spehre
          //2. apply boundary conditions

          double l=0.0;
          int code;

          l=sqrt(dx*dx+dy*dy+dz*dz);
          l=(radiusSphere+PIC::Mesh::mesh.EPS)/l;

          x[0]=x0Sphere[0]+l*dx;
          x[1]=x0Sphere[1]+l*dy;
          x[2]=x0Sphere[2]+l*dz;
          startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);

          FirstBoundaryFlag=true;

          code=Sphere->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
          if (code==_PARTICLE_DELETED_ON_THE_FACE_) return _PARTICLE_LEFT_THE_DOMAIN_;

          goto MovingLoop;
        }

        d=b*b-4.0*a*c;

        if (d<0.0) {
          if (4.0*a*pow(PIC::Mesh::mesh.EPS,2)>-d) d=0.0; //equvalent to |EPS/particle speed| > sqrt(|d|)/(2a) -> the particle is within the distance of EPS from the surface's boundary
          else continue;
        }

        if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

        sqrt_d=sqrt(d);
        dtTemp=-(b+sqrt_d)/(2.0*a);
//        if ((dtTemp>0.0)&&(dtTemp*dtTemp*a<PIC::Mesh::mesh.EPS*PIC::Mesh::mesh.EPS)) dtTemp=-1.0;

        dt1=-2.0*c/(b+sqrt_d);
        if ((dtTemp<0.0)||((dt1>0.0)&&(dt1<dtTemp))) {
          dtTemp=dt1;
//          if ((dtTemp>0.0)&&(dtTemp*dtTemp*a<PIC::Mesh::mesh.EPS*PIC::Mesh::mesh.EPS)) dtTemp=-1.0;
        }

        if ((0.0<dtTemp)&&(dtTemp<dtMin)) {
          dtMin=dtTemp,InternalBoundaryDescriptor_dtMin=InternalBoundaryDescriptor;
          ParticleIntersectionCode=_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
        }

        break;
      default:
        exit(__LINE__,__FILE__,"Error: undetermined internal boundary type");
      }
    }


#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //check if a photolytic reaction is possible and get the time interval before the transformation occures
    int PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(x,ptr,spec,dtMin);
#endif

    //advance the particle's position
    for (idim=0;idim<DIM;idim++) {

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
      xInit[idim]=x[idim];
#endif

      x[idim]+=dtMin*v[idim];
    }

    FirstBoundaryFlag=false;

#if _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ == _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ON_
    v[0]+=dtMin*accl[0],v[1]+=dtMin*accl[1],v[2]+=dtMin*accl[2];
#endif

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //model the photolytic transformation
    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) {
      int specInit=spec;

      PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable[spec](xInit,x,ptr,spec,ParticleData);

      //adjust the value of the dtLeft to match the time step for the species 'spec'
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      dtTotal*=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]/PIC::ParticleWeightTimeStep::GlobalTimeStep[specInit];
#else
      exit(__LINE__,__FILE__,"Error: not implemeted");
#endif


      if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

      MovingTimeFinished=false;
      ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;
    }
#endif

    //adjust the particle moving time
    dtTotal-=dtMin;

    //interaction with the faces of the block and internal surfaces
    if (ParticleIntersectionCode==_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      int code;

      FirstBoundaryFlag=true;

      lastInternalBoundaryDescriptor=InternalBoundaryDescriptor_dtMin;
      code=((cInternalSphericalData*)(InternalBoundaryDescriptor_dtMin->BoundaryElement))->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor_dtMin->BoundaryElement);


//      code=PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction(spec,ptr,x,v,dtTotal,startNode,InternalBoundaryDescriptor_dtMin);


      if (code==_PARTICLE_DELETED_ON_THE_FACE_) return _PARTICLE_LEFT_THE_DOMAIN_;
      startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
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

//      for (idim=0;idim<DIM;idim++) if ((x[idim]<xminBlock[idim]-PIC::Mesh::mesh.EPS)||(x[idim]>xmaxBlock[idim]+PIC::Mesh::mesh.EPS)) {

      if ((x[0]<xminBlock[0]-PIC::Mesh::mesh.EPS)||(x[0]>xmaxBlock[0]+PIC::Mesh::mesh.EPS)
#if DIM > 1
          || (x[1]<xminBlock[1]-PIC::Mesh::mesh.EPS)||(x[1]>xmaxBlock[1]+PIC::Mesh::mesh.EPS)
#endif
#if DIM > 2
          || (x[2]<xminBlock[2]-PIC::Mesh::mesh.EPS)||(x[2]>xmaxBlock[2]+PIC::Mesh::mesh.EPS)
#endif
      ) {
        exit(__LINE__,__FILE__,"Error: the new particles' coordinates are outside of the block");
      }
      #endif

      //move the particle's position exactly to the block's face
      for (idim=0;idim<DIM;idim++) {
        if (x[idim]<=xminBlock[idim]) x[idim]=xminBlock[idim]+PIC::Mesh::mesh.EPS;
        if (x[idim]>=xmaxBlock[idim]) x[idim]=xmaxBlock[idim]-PIC::Mesh::mesh.EPS;
      }

//      x[neibNodeDirection]=(iNeibNode[neibNodeDirection]==-1) ? xmaxBlock[neibNodeDirection] : xminBlock[neibNodeDirection];

      //reserve the place for particle's cloning

      //adjust the value of 'startNode'
      startNode=newNode;
    }
    else startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);



  }




  //the particle is still within the computational domain:
  //place it to the local list of particles related to the new block and cell
//  long int LocalCellNumber;
//  int i,j,k;

//  PIC::Mesh::cDataCenterNode *cell;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=startNode->block->GetCenterNode(LocalCellNumber);


  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;



  //=====================  DEBUG =========================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if ((cell->Measure<=0.0)&&(startNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << __FILE__<< __LINE__ << endl;


    cout << "Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

     double r,rprobe[3]={0.0,0.0,0.0};
     int di,dj,dk;


     cout << "x particle=";
     for (r=0.0,idim=0;idim<DIM;idim++) {
       r+=pow(x[idim],2);
       cout << x[idim] << " ";
     }

     cout << ", |x|= " << sqrt(r) << endl;

     for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
       startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

       for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
       cout << "Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
     }


    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif
  //===================   END DEBUG ==============================


  /*

  if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1737.0E3) {
    cout << __FILE__ << __LINE__ << endl;
  }
*/
  //===================== END DEBUG ==================




  return _PARTICLE_MOTION_FINISHED_;
}

int PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode)  {
  return UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(ptr,dtTotal,startNode,false);
}

int PIC::Mover::UniformWeight_UniformTimeStep_noForce(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  PIC::ParticleBuffer::byte *ParticleData;
  double v[3],x[3],vinit[3],xinit[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  long int LocalCellNumber;
  int i,j,k;
  PIC::Mesh::cDataCenterNode *cell;


  //########## DEBUG ############

  static long int nCallCounter=0;

  nCallCounter++;


  if (nCallCounter==5053963548) {
    cout << __LINE__ << __FILE__<< endl;
  }
  //########## END DEBUG ##########




  #if  _SIMULATION_TIME_STEP_MODE_ ==  _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #endif


  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(v,ParticleData);
  PIC::ParticleBuffer::GetX(x,ParticleData);

  //=====================  DEBUG =========================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(startNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

    double r,rprobe[3]={0.0,0.0,0.0};
    int di,dj,dk,idim;


    cout << "x particle=";
    for (r=0.0,idim=0;idim<DIM;idim++) {
      r+=pow(x[idim],2);
      cout << x[idim] << " ";
    }

    cout << ", |x|= " << sqrt(r) << endl;

    for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
      startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

      for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
      cout << "Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
    }

    PIC::Mesh::mesh.InitCellMeasure(startNode);

    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif


  //===================== END DEBUG ==================


#if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
  if (startNode->InternalBoundaryDescriptorList!=NULL) {
   return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt,startNode);
  }
#endif


  double dtLeft=0.0;
  double accl[3];
  int spec;

  spec=PIC::ParticleBuffer::GetI(ParticleData);
  dtLeft=dt;

  while (dtLeft>0.0) {
    double aa,vv;

    dt=dtLeft;
    dtLeft=0.0;

    TotalParticleAcceleration(accl,spec,ptr,x,v,startNode);

    aa=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    vv=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    if (aa*dt>0.2*vv) {
      dtLeft=dt;
      dt=0.2*vv/aa;
      dtLeft-=dt;
    }

    memcpy(vinit,v,3*sizeof(double));
    memcpy(xinit,x,3*sizeof(double));

    //Check if the startNode has cut cells
    /*
#if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
    if (startNode->InternalBoundaryDescriptorList!=NULL) {
      PIC::ParticleBuffer::SetV(v,ParticleData);
      PIC::ParticleBuffer::SetX(x,ParticleData);

     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt+dtLeft,startNode);
    }
#endif
*/


    //check the occurence of photolytic reactions
    //1. check if a reaction is possible
    //2. move the particle for the time interval before the reaction has occured
    //3. model the particle transformation

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //check if a photolytic reaction is possible and get the time interval before the transformation occures
    double dtInit;
    int PhotolyticReactionsReturnCode;

    dtInit=dt;
    PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(x,ptr,spec,dt);

    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) dtLeft+=dtInit-dt;
#endif

    //advance the particle positions
    #if DIM == 3
    x[0]+=dt*v[0],x[1]+=dt*v[1],x[2]+=dt*v[2];

    //first order trajectory integration
    #if _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ == _PARTICLE_TRAJECTORY_FORCE_INTEGRTAION_MODE_ON_
    v[0]+=dt*accl[0],v[1]+=dt*accl[1],v[2]+=dt*accl[2];
    #endif

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //model the photolytic transformation
    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) {
      int specInit=spec;

      PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable[spec](xinit,x,ptr,spec,ParticleData);

      //adjust the value of the dtLeft to match the time step for the species 'spec'
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      dtLeft*=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]/PIC::ParticleWeightTimeStep::GlobalTimeStep[specInit];
#else
      exit(__LINE__,__FILE__,"Error: not implemeted");
#endif

      if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }
    }
#endif


    #else
    exit(__LINE__,__FILE__,"not implemented");
    #endif

    //determine the new particle location
    newNode=PIC::Mesh::mesh.findTreeNode(x,startNode);


    if (newNode==NULL) {
      //the particle left the computational domain
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //Check if the newNode has cut cells
    #if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
    if (newNode->InternalBoundaryDescriptorList!=NULL) {
      PIC::ParticleBuffer::SetV(vinit,ParticleData);
      PIC::ParticleBuffer::SetX(xinit,ParticleData);

     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt+dtLeft,startNode);
    }
    #endif

    startNode=newNode;
  }

  //the particle is still within the computational domain:
  //place it to the local list of particles related to the new block and cell
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,newNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=newNode->block->GetCenterNode(LocalCellNumber);

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetX(x,ParticleData);

  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;


  //=====================  DEBUG =========================


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if ((cell->Measure<=0.0)&&(newNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << __FILE__<< __LINE__ << endl;
    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif



  /*
  if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1737.0E3) {
    cout << __FILE__ << __LINE__ << endl;
  }
*/
  //===================== END DEBUG ==================



  return _PARTICLE_MOTION_FINISHED_;
}
