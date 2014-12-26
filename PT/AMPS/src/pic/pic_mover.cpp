//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that controls the particle motion

#include "pic.h"


//PIC::Mover::fSpeciesDependentParticleMover *PIC::Mover::MoveParticleTimeStep=NULL;
//PIC::Mover::fTotalParticleAcceleration PIC::Mover::TotalParticleAcceleration=PIC::Mover::TotalParticleAcceleration_default;
//PIC::Mover::fSpeciesDependentParticleMover_BoundaryInjection *PIC::Mover::MoveParticleBoundaryInjection=NULL;
PIC::Mover::fProcessOutsideDomainParticles PIC::Mover::ProcessOutsideDomainParticles=NULL;
PIC::Mover::fProcessTriangleCutFaceIntersection PIC::Mover::ProcessTriangleCutFaceIntersection=NULL;

//====================================================
//init the particle mover
void PIC::Mover::Init() {
/*
  int s;


  //allocate the mover functions array
  if ((MoveParticleTimeStep!=NULL)||(PIC::nTotalSpecies==0)) exit(__LINE__,__FILE__,"Error: the initialization of PIC::Mover is failed");

  MoveParticleTimeStep=new fSpeciesDependentParticleMover[PIC::nTotalSpecies];
  MoveParticleBoundaryInjection=new fSpeciesDependentParticleMover_BoundaryInjection[PIC::nTotalSpecies];
  for (s=0;s<PIC::nTotalSpecies;s++) MoveParticleTimeStep[s]=NULL,MoveParticleBoundaryInjection[s]=NULL;
  */
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
  long int /*LocalCellNumber,*/ ParticleList,ptr;
  double LocalTimeStep;



//  return;


  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

//  int centerNodeIndexCounter;

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
  PIC::Mesh::cDataBlockAMR *block;

  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  //move existing particles
  while (node!=NULL) {

    block=node->block;
    memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            /*LocalCellNumber=*/  PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            ParticleList=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            /*
    {
      {
//        for (LocalCellNumber=nStartCenterNodeIndex,centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;LocalCellNumber+=centerNodeIndexIncrement[centerNodeIndexCounter++]) {

        for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {

            LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
            ParticleList=node->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;
*/
            while (ParticleList!=-1) {
              ptr=ParticleList;
              ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);
              s=PIC::ParticleBuffer::GetI(ptr);
              LocalTimeStep=block->GetLocalTimeStep(s);

//              MoveParticleTimeStep[s](ptr,LocalTimeStep,node);
              _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,LocalTimeStep,node);
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
//  PIC::Mesh::cDataCenterNode *cell;

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

  for (k=0;k<_BLOCK_CELLS_Z_;k++) {
     for (j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {

          //LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

          FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
        }
     }
  }


  for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
    node=(thread==PIC::Mesh::mesh.ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    if (node==NULL) continue;

    while (node!=NULL) {

/*
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {

              //LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

              FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
            }
         }
      }
      */

      block=node->block;
      memcpy(block->FirstCellParticleTable,block->tempParticleMovingListTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
      memcpy(block->tempParticleMovingListTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));




      /*
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


      */

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
  double v[3],x[3]={0.0,0.0,0.0},*xminBlock,*xmaxBlock;
  int iNeibNode[3],neibNodeDirection,idim,nface_dtMin=-1;

  long int LocalCellNumber;
  int i,j,k;

  PIC::Mesh::cDataCenterNode *cell;
  bool MovingTimeFinished=false;


#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  double xInit[3];
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
  double xInit[3];
#endif
  //=====================  DEBUG =========================




//########## DEBUG ############

static long int nCallCounter=0;

nCallCounter++;

/*
if ((nCallCounter==1057317)&&(PIC::Mesh::mesh.ThisThread==5)) {
  cout << __FILE__ << "@" << __LINE__ << endl;
}
*/

//########## END DEBUG ##########



  //the descriptors of the internal surfaces
  cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*InternalBoundaryDescriptor_dtMin=NULL,*lastInternalBoundaryDescriptor=NULL;

  //spherical internal surface
#if DIM == 3
  cInternalSphericalData *Sphere;
  cInternalRotationBodyData *Nucleus;
#elif DIM == 2
  exit(__LINE__,__FILE__,"not yet");
#else
  cInternalSphere1DData *Sphere1D;
#endif

  double radiusSphere,*x0Sphere;
  double a,b,c,d,dx,dy,dz,sqrt_d,dt1;

  double *x0Nucleus,*lNucleus;
  double xmin,xmax,rSurface;

#define _UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_        0
#define _BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_       1
#define _INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_  2

  int ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;


  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(v,ParticleData);
  PIC::ParticleBuffer::GetX(x,ParticleData);
//  spec=PIC::ParticleBuffer::GetI(ParticleData);



//#######  DEBUG #########
//double xpinit[3];
//for (idim=0;idim<3;idim++) xpinit[idim]=x[idim];

//=====================  DEBUG =========================


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

  {
    double r[3]={0.0,0.0,0.0};

    r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(r,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  }

#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if (cell->Measure<=0.0) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;
    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif


/*
if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1737.0E3) {
  cout << __FILE__ << __LINE__ << endl;
}

*/

  /*
  if ((nCallCounter==1057317)&&(PIC::Mesh::mesh.ThisThread==5)) {
    cout << __FILE__ << "@" << __LINE__ << endl;
  }
  */
//===================== END DEBUG ==================



//####### END DEBUG #########

//  while (dtTotal>0.0) {
  while (MovingTimeFinished==false) {
MovingLoop:

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    //check the consistency of the particle mover
int iTemp,jTemp,kTemp;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_


    if (PIC::Mesh::mesh.fingCellIndex(x,iTemp,jTemp,kTemp,startNode,false)==-1) {
      exit(__LINE__,__FILE__,"Error: the cell is not found");
    }
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

    {
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      if (PIC::Mesh::mesh.fingCellIndex(r,iTemp,jTemp,kTemp,startNode,false)==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }

#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif

    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif



    xminBlock=startNode->xmin;
    xmaxBlock=startNode->xmax;

    MovingTimeFinished=true;
    dtMin=dtTotal,nface_dtMin=-1;
    InternalBoundaryDescriptor_dtMin=NULL;
    ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    double accl[3],vv;
    int spec;

    spec=PIC::ParticleBuffer::GetI(ParticleData);
    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(accl,spec,ptr,x,v,startNode);

    a=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    vv=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (a*dtMin>0.2*vv) dtMin=0.2*vv/a;
#endif


    //Calculate the time of flight to the nearest block's face
#if DIM == 1

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    if (fabs(v[0])>0.0) {
      dtTemp=(xminBlock[0]-x[0])/v[0];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=(xmaxBlock[0]-x[0])/v[0];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

    //nface=0 -> rmin
    a=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    b=2.0*v[0]*x[0];
    c=x[0]*x[0]-xminBlock[0]*xminBlock[0];
    d=b*b-4.0*a*c;

    if (d>0.0) {
      if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

      sqrt_d=sqrt(d);
      dtTemp=-(b+sqrt_d)/(2.0*a);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=-2.0*c/(b+sqrt_d);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }

    //nface=1 -> rmax
    c=x[0]*x[0]-xmaxBlock[0]*xmaxBlock[0];
    d=b*b-4.0*a*c;

    if (d>0.0) {
      if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

      sqrt_d=sqrt(d);
      dtTemp=-(b+sqrt_d)/(2.0*a);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=-2.0*c/(b+sqrt_d);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }

#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

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
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_: case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:

#if DIM == 3
        Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#elif DIM == 2
        exit(__LINE__,__FILE__,"not implemented");
#else
        Sphere1D=(cInternalSphere1DData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere1D->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#endif

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
#if DIM == 3
          code=Sphere->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#elif DIM == 2
        exit(__LINE__,__FILE__,"not implemented");
#else
        code=Sphere1D->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#endif

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
      case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
        Nucleus=(cInternalRotationBodyData*)(InternalBoundaryDescriptor->BoundaryElement);
	Nucleus->GetSphereGeometricalParameters(x0Nucleus,lNucleus,xmin,xmax);

        Nucleus->SurfaceCurve(rSurface,x[0]);
        if(xmin<=x[0] && xmax>=x[0] && rSurface>sqrt(x[1]*x[1]+x[2]*x[2])) {
	  //the particle is inside the nucleus                                                                                                                  
	  PIC::ParticleBuffer::DeleteParticle(ptr);
          return _PARTICLE_LEFT_THE_DOMAIN_;
        }

        break;
      default:
        exit(__LINE__,__FILE__,"Error: undetermined internal boundary type");
      }
    }


#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //check if a photolytic reaction is possible and get the time interval before the transformation occures
    int PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(x,ptr,spec,dtMin,startNode);
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
    bool TransformationTimeStepLimitFlag=false;

    int GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_(x,v,spec,ptr,ParticleData,dtMin,TransformationTimeStepLimitFlag,startNode);
#endif

    //advance the particle's position
    for (idim=0;idim<3;idim++) {

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
      xInit[idim]=x[idim];
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
      xInit[idim]=x[idim];
#endif

      x[idim]+=dtMin*v[idim];
    }

    FirstBoundaryFlag=false;

#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    v[0]+=dtMin*accl[0],v[1]+=dtMin*accl[1],v[2]+=dtMin*accl[2];
#endif

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //model the photolytic transformation
    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) {
      int specInit=spec;

      //PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable[specInit](xInit,x,ptr,spec,ParticleData);
      PhotolyticReactionsReturnCode=_PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(xInit,x,ptr,spec,ParticleData);

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
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
   //model the generic particle transformation
   if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_) {
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
     int specInit=spec;
#endif

     GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xInit,x,v,spec,ptr,ParticleData,dtMin,startNode);   //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

     //adjust the value of the dtLeft to match the time step for the species 'spec'
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
     dtTotal*=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]/PIC::ParticleWeightTimeStep::GlobalTimeStep[specInit];
#else
     exit(__LINE__,__FILE__,"Error: not implemeted");
#endif


     if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_) {
       PIC::ParticleBuffer::DeleteParticle(ptr);
       return _PARTICLE_LEFT_THE_DOMAIN_;
     }

//     MovingTimeFinished=false;
//     ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;
   }
#endif

    //adjust the particle moving time
    dtTotal-=dtMin;

    //interaction with the faces of the block and internal surfaces
    if (ParticleIntersectionCode==_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      int code;

      FirstBoundaryFlag=true;

      lastInternalBoundaryDescriptor=InternalBoundaryDescriptor_dtMin;

#if DIM == 3
      code=((cInternalSphericalData*)(InternalBoundaryDescriptor_dtMin->BoundaryElement))->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor_dtMin->BoundaryElement);
#elif DIM == 2
      exit(__LINE__,__FILE__,"not implemented");
#else
      code=((cInternalSphere1DData*)(InternalBoundaryDescriptor_dtMin->BoundaryElement))->ParticleSphereInteraction(spec,ptr,x,v,dtTotal,(void*)startNode,InternalBoundaryDescriptor_dtMin->BoundaryElement);
#endif

//      code=PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction(spec,ptr,x,v,dtTotal,startNode,InternalBoundaryDescriptor_dtMin);


      if (code==_PARTICLE_DELETED_ON_THE_FACE_) return _PARTICLE_LEFT_THE_DOMAIN_;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      startNode=PIC::Mesh::mesh.findTreeNode(r,startNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif
    }
    else if (ParticleIntersectionCode==_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      iNeibNode[0]=0,iNeibNode[1]=0,iNeibNode[2]=0;
      neibNodeDirection=(int)(nface_dtMin/2);
      iNeibNode[neibNodeDirection]=(nface_dtMin-2*neibNodeDirection==0) ? -1 : 1;


      double xProbe[3]={x[0],x[1],x[2]};
      xProbe[neibNodeDirection]+=(startNode->xmax[neibNodeDirection]-startNode->xmin[neibNodeDirection])*iNeibNode[neibNodeDirection]/100.0;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xProbe,startNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double rProbe[3]={0.0,0.0,0.0};

      rProbe[0]=sqrt(xProbe[0]*xProbe[0]+xProbe[1]*xProbe[1]+xProbe[2]*xProbe[2]);
      newNode=PIC::Mesh::mesh.findTreeNode(rProbe,startNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


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

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
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
      #endif

      //move the particle's position exactly to the block's face
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      for (idim=0;idim<DIM;idim++) {
        if (x[idim]<=xminBlock[idim]) x[idim]=xminBlock[idim]+PIC::Mesh::mesh.EPS;
        if (x[idim]>=xmaxBlock[idim]) x[idim]=xmaxBlock[idim]-PIC::Mesh::mesh.EPS;
      }
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

      if (r2<=xminBlock[0]*xminBlock[0]) {
        double c=(xminBlock[0]+PIC::Mesh::mesh.EPS)/sqrt(r2);
        x[0]*=c,x[1]*=c,x[2]*=c;
      }

      if (r2>=xmaxBlock[0]*xmaxBlock[0]){
        double c=(xmaxBlock[0]-PIC::Mesh::mesh.EPS)/sqrt(r2);
        x[0]*=c,x[1]*=c,x[2]*=c;
      }
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif
//      x[neibNodeDirection]=(iNeibNode[neibNodeDirection]==-1) ? xmaxBlock[neibNodeDirection] : xminBlock[neibNodeDirection];

      //reserve the place for particle's cloning

      //adjust the value of 'startNode'
      startNode=newNode;
    }
    else {
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      startNode=PIC::Mesh::mesh.findTreeNode(r,startNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif
    }



  }




  //the particle is still within the computational domain:
  //place it to the local list of particles related to the new block and cell
//  long int LocalCellNumber;
//  int i,j,k;

//  PIC::Mesh::cDataCenterNode *cell;


  //Rotate particle position and velocity when symmetry is accounted
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
  //do nothing
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
  double r,l,v1[3],cosTz,sinTz,cosTy,sinTy;

  r=sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));

  if (r>1.0E-20) {
    double xfinal[3];

    xfinal[0]=x[0]/r,xfinal[1]=x[1]/r,xfinal[2]=x[1]/r;
    l=sqrt(pow(xfinal[0],2)+pow(xfinal[1],2));

    if (l>1.0E-20) {
      cosTz=xfinal[0]/l,sinTz=xfinal[1]/l;
      cosTy=l,sinTy=xfinal[2];
    }
    else cosTz=1.0,sinTz=0.0,sinTy=xfinal[2],cosTy=0.0;

    v1[0]=cosTy*cosTz*v[0]+cosTy*sinTz*v[1]+sinTy*v[2];
    v1[1]=-sinTz*v[0]+cosTz*v[1];
    v1[2]=-sinTy*cosTz*v[0]-sinTy*sinTz*v[1]+cosTy*v[2];

    v[0]=v1[0],v[1]=v1[1],v[2]=v1[2];
    x[0]=r,x[1]=0.0,x[2]=0.0;
  }
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif


  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
//  cell=startNode->block->GetCenterNode(LocalCellNumber);

  PIC::Mesh::cDataBlockAMR *block=startNode->block;
  long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetX(x,ParticleData);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;


  /*
  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;
*/


  //=====================  DEBUG =========================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(startNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;


    cout << "$PREFIX:Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

     double r,rprobe[3]={0.0,0.0,0.0};
     int di,dj,dk;


     cout << "$PREFIX:x particle=";
     for (r=0.0,idim=0;idim<DIM;idim++) {
       r+=pow(x[idim],2);
       cout << x[idim] << " ";
     }

     cout << ", |x|= " << sqrt(r) << endl;

     for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
       startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

       for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
       cout << "$PREFIX:Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
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
  double v[3],x[3]={0.0,0.0,0.0},vinit[3],xinit[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode=NULL;
  long int LocalCellNumber;
  int i,j,k;
  PIC::Mesh::cDataCenterNode *cell;


  //########## DEBUG ############

  static long int nCallCounter=0;

  nCallCounter++;


  /*
  if (nCallCounter==5053963) {
    cout << __LINE__ << __FILE__<< endl;
  }
  */
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
    cout << "$PREFIX:Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

    double r,rprobe[3]={0.0,0.0,0.0};
    int di,dj,dk,idim;


    cout << "$PREFIX:x particle=";
    for (r=0.0,idim=0;idim<DIM;idim++) {
      r+=pow(x[idim],2);
      cout << x[idim] << " ";
    }

    cout << ", |x|= " << sqrt(r) << endl;

    for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
      startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

      for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
      cout << "$PREFIX:Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
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

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(accl,spec,ptr,x,v,startNode);

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
    PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(x,ptr,spec,dt,startNode);

    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) dtLeft+=dtInit-dt;
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
    int GenericParticleTransformationReturnCode;
    bool TransformationTimeStepLimitFlag=false;

    dtLeft+=dt;
    GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_(x,v,spec,ptr,ParticleData,dt,TransformationTimeStepLimitFlag,startNode);
    dtLeft-=dt;
#endif

    //advance the particle positions

    x[0]+=dt*v[0],x[1]+=dt*v[1],x[2]+=dt*v[2];

    //first order trajectory integration
    #if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    v[0]+=dt*accl[0],v[1]+=dt*accl[1],v[2]+=dt*accl[2];
    #endif

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //model the photolytic transformation
    if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTION_OCCURES_) {
      int specInit=spec;

//      PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable[specInit](xinit,x,ptr,spec,ParticleData);
      PhotolyticReactionsReturnCode=_PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(xinit,x,ptr,spec,ParticleData);

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
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
    //model the generic particle transformation
    if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_) {
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      int specInit=spec;
#endif

      GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xinit,x,v,spec,ptr,ParticleData,dt,startNode);  //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

      //adjust the value of the dtLeft to match the time step for the species 'spec'
 #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      dtLeft*=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]/PIC::ParticleWeightTimeStep::GlobalTimeStep[specInit];
 #else
      exit(__LINE__,__FILE__,"Error: not implemeted");
 #endif


      if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }
    }


#endif




    //determine the new particle location
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    newNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
    double r[3]={0.0,0.0,0.0};

    r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    newNode=PIC::Mesh::mesh.findTreeNode(r,startNode);
#else
    exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


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
  //Rotate particle position and velocity when symmetry is accounted
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
  //do nothing
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
  double r,l,v1[3],cosTz,sinTz,cosTy,sinTy;

  r=sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));

  if (r>1.0E-20) {
    double xfinal[3];

    xfinal[0]=x[0]/r,xfinal[1]=x[1]/r,xfinal[2]=x[1]/r;
    l=sqrt(pow(xfinal[0],2)+pow(xfinal[1],2));

    if (l>1.0E-20) {
      cosTz=xfinal[0]/l,sinTz=xfinal[1]/l;
      cosTy=l,sinTy=xfinal[2];
    }
    else cosTz=1.0,sinTz=0.0,sinTy=xfinal[2],cosTy=0.0;

    v1[0]=cosTy*cosTz*v[0]+cosTy*sinTz*v[1]+sinTy*v[2];
    v1[1]=-sinTz*v[0]+cosTz*v[1];
    v1[2]=-sinTy*cosTz*v[0]-sinTy*sinTz*v[1]+cosTy*v[2];

    v[0]=v1[0],v[1]=v1[1],v[2]=v1[2];
    x[0]=r,x[1]=0.0,x[2]=0.0;
  }
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif


  //place it to the local list of particles related to the new block and cell
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,newNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  //cell=newNode->block->GetCenterNode(LocalCellNumber);

  PIC::Mesh::cDataBlockAMR *block=startNode->block;
  long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetX(x,ParticleData);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;

  /*
  PIC::ParticleBuffer::SetNext(cell->tempParticleMovingList,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (cell->tempParticleMovingList!=-1) PIC::ParticleBuffer::SetPrev(ptr,cell->tempParticleMovingList);
  cell->tempParticleMovingList=ptr;
  */


  //=====================  DEBUG =========================


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  cell=newNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(newNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;
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



//======================================================================================
int PIC::Mover::UniformWeight_UniformTimeStep_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0},vMiddle[3],xMiddle[3],vFinal[3],xFinal[3],dt2;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode=NULL,*middleNode=NULL;
  long int LocalCellNumber;
  int i,j,k;
  PIC::Mesh::cDataCenterNode *cell;
  bool IntegrationInterrupted=true;

  #if  _SIMULATION_TIME_STEP_MODE_ ==  _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  exit(__LINE__,__FILE__,"Error: the function cannot be applied for such configuration");
  #endif


  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(startNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "$PREFIX:Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

    double r,rprobe[3]={0.0,0.0,0.0};
    int di,dj,dk,idim;


    cout << "$PREFIX:x particle=";
    for (r=0.0,idim=0;idim<DIM;idim++) {
      r+=pow(xInit[idim],2);
      cout << xInit[idim] << " ";
    }

    cout << ", |x|= " << sqrt(r) << endl;

    for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
      startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

      for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
      cout << "$PREFIX:Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
    }

    PIC::Mesh::mesh.InitCellMeasure(startNode);

    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif




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

  while (IntegrationInterrupted==true) {
    double aa,vv;

    dt=dtLeft;
    dtLeft=0.0;
    IntegrationInterrupted=false;

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(accl,spec,ptr,xInit,vInit,startNode);

    aa=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    vv=sqrt(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2]);

    if (aa*dt>2.0*vv) {
      dtLeft=dt;
      dt=2.0*vv/aa;
      dtLeft-=dt;
      IntegrationInterrupted=true;
    }


    //predictor
    dt2=dt/2.0;
    xMiddle[0]=xInit[0]+dt2*vInit[0];
    xMiddle[1]=xInit[1]+dt2*vInit[1];
    xMiddle[2]=xInit[2]+dt2*vInit[2];

    vMiddle[0]=vInit[0]+dt2*accl[0];
    vMiddle[1]=vInit[1]+dt2*accl[1];
    vMiddle[2]=vInit[2]+dt2*accl[2];

    //corrector

    //in the case when a symmetry has to be considered, transfer the particle position and velcoty accordinly
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    middleNode=PIC::Mesh::mesh.findTreeNode(xMiddle,startNode);
    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(accl,spec,ptr,xMiddle,vMiddle,newNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
    exit(__LINE__,__FILE__,"not implemented");
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

    if (middleNode==NULL) {
      //the particle left the computational domain
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //Check if the newNode has cut cells
    #if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
    if (middleNode->InternalBoundaryDescriptorList!=NULL) {
      PIC::ParticleBuffer::SetV(vInit,ParticleData);
      PIC::ParticleBuffer::SetX(xInit,ParticleData);

     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt+dtLeft,startNode);
    }
    #endif


    xFinal[0]=xInit[0]+dt*vMiddle[0];
    xFinal[1]=xInit[1]+dt*vMiddle[1];
    xFinal[2]=xInit[2]+dt*vMiddle[2];

    vFinal[0]=vInit[0]+dt*accl[0];
    vFinal[1]=vInit[1]+dt*accl[1];
    vFinal[2]=vInit[2]+dt*accl[2];

    //in the case when a symmetry has to be considered, transfer the particle position and velcoty accordinly
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
    exit(__LINE__,__FILE__,"not implemented");
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif


    if (newNode==NULL) {
      //the particle left the computational domain
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //Check if the newNode has cut cells
    #if _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
    if (middleNode->InternalBoundaryDescriptorList!=NULL) {
      PIC::ParticleBuffer::SetV(vInit,ParticleData);
      PIC::ParticleBuffer::SetX(xInit,ParticleData);

     return UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,dt+dtLeft,startNode);
    }
    #endif

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
exit(__LINE__,__FILE__,"not implemented");
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
exit(__LINE__,__FILE__,"not implemented");
#endif

    //prepare for the next pass of the loop
    startNode=newNode;
    memcpy(vInit,vFinal,3*sizeof(double));
    memcpy(xInit,xFinal,3*sizeof(double));
  }

  //place it to the local list of particles related to the new block and cell
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,newNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");


  PIC::Mesh::cDataBlockAMR *block=startNode->block;
  long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;


  //=====================  DEBUG =========================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  cell=newNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(newNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;
    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif
  //====================  END DEBUG ======================

  return _PARTICLE_MOTION_FINISHED_;
}

int PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL,*middleNode=NULL;
  double dtMin=-1.0,dtTemp,dtMinInit2,dtMinInit;
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0},vMiddle[3],xMiddle[3],vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
  int idim;
//  int iNeibNode[3],neibNodeDirection;

  long int LocalCellNumber;
  int i,j,k,spec;

  PIC::Mesh::cDataCenterNode *cell;
  bool MovingTimeFinished=false;

  //the descriptors of the internal surfaces
  cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor,*InternalBoundaryDescriptor_dtMin=NULL,*lastInternalBoundaryDescriptor=NULL;
  CutCell::cTriangleFace *IntersectionFace=NULL;
  CutCell::cTriangleFace *lastIntersectedTriangleFace=NULL;

  //the description of the boundaries of the block faces
  struct cExternalBoundaryFace {
    double norm[3];
    int nX0[3];
    double e0[3],e1[3],x0[3];
    double lE0,lE1;
  };

  static bool initExternalBoundaryFaceTable=false;

  static cExternalBoundaryFace ExternalBoundaryFaceTable[6]={
      {{-1.0,0.0,0.0}, {0,0,0}, {0,1,0},{0,0,1},{0,0,0}, 0.0,0.0}, {{1.0,0.0,0.0}, {1,0,0}, {0,1,0},{0,0,1},{0,0,0}, 0.0,0.0},
      {{0.0,-1.0,0.0}, {0,0,0}, {1,0,0},{0,0,1},{0,0,0}, 0.0,0.0}, {{0.0,1.0,0.0}, {0,1,0}, {1,0,0},{0,0,1},{0,0,0}, 0.0,0.0},
      {{0.0,0.0,-1.0}, {0,0,0}, {1,0,0},{0,1,0},{0,0,0}, 0.0,0.0}, {{0.0,0.0,1.0}, {0,0,1}, {1,0,0},{0,1,0},{0,0,0}, 0.0,0.0}
  };

  if (initExternalBoundaryFaceTable==false) {
    initExternalBoundaryFaceTable=true;

    for (int nface=0;nface<6;nface++) {
      double cE0=0.0,cE1=0.0;

      for (int idim=0;idim<3;idim++) {
        ExternalBoundaryFaceTable[nface].x0[idim]=(ExternalBoundaryFaceTable[nface].nX0[idim]==0) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim];

        cE0+=pow(((ExternalBoundaryFaceTable[nface].e0[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
        cE1+=pow(((ExternalBoundaryFaceTable[nface].e1[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
      }

      ExternalBoundaryFaceTable[nface].lE0=sqrt(cE0);
      ExternalBoundaryFaceTable[nface].lE1=sqrt(cE1);
    }
  }

  //spherical internal surface
#if DIM == 3
  cInternalSphericalData *Sphere;
  cInternalRotationBodyData *Nucleus;
#elif DIM == 2
  exit(__LINE__,__FILE__,"not yet");
#else
  cInternalSphere1DData *Sphere1D;
#endif

  double radiusSphere,*x0Sphere;
  double a,b,c,d,dx,dy,dz,sqrt_d,dt1;

  double *x0Nucleus,*lNucleus;
  double xmin,xmax,rSurface;

#define _UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_        0
#define _BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_       1
#define _INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_  2
#define _BOUNDARY_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_    3

  int ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);


//=====================  DEBUG ==================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_

  if (PIC::ParticleBuffer::IsParticleAllocated(ParticleData)==false) {
	  exit(__LINE__,__FILE__,"Error: an unallocated particle is intercepted");
  }

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

  {
    double r[3]={0.0,0.0,0.0};

    r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(r,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
  }

#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if (cell->Measure<=0.0) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;

    double  vol=-1,xmin[3],xmax[3],xmiddle[3];

    xmin[0]=startNode->xmin[0]+i*(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_;
    xmin[1]=startNode->xmin[1]+j*(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_;
    xmin[2]=startNode->xmin[2]+k*(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_;

    xmax[0]=startNode->xmin[0]+(i+1)*(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_;
    xmax[1]=startNode->xmin[1]+(j+1)*(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_;
    xmax[2]=startNode->xmin[2]+(k+1)*(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_;

    vol=CutCell::GetRemainedBlockVolume(xmin,xmax,PIC::Mesh::mesh.EPS,1.0E-2,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,startNode->FirstTriangleCutFace);
    cout << "$PREFIX: recalculated volume: vol=" << vol << "(" << __FILE__ << "@" << __LINE__ << ")" << endl;
    cout << " xInit=" << xInit[0] << ", " << xInit[1] << ", " << xInit[2] << endl;

    if (CutCell::CheckPointInsideDomain(xInit,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
      cout << "$PREFIX: xInit is outside of the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
    }
    else {
      cout << "$PREFIX: xInit is inside of the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
    }

    xmiddle[0]=startNode->xmin[0]+(i+0.5)*(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_;
    xmiddle[1]=startNode->xmin[1]+(j+0.5)*(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_;
    xmiddle[2]=startNode->xmin[2]+(k+0.5)*(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_;

    cout << " xMiddle=" << xmiddle[0] << ", " << xmiddle[1] << ", " << xmiddle[2] << endl;

    if (CutCell::CheckPointInsideDomain(xmiddle,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
      cout << "$PREFIX: middle point of the cell is outside of the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
    }
    else {
      cout << "$PREFIX: middle point of the cell is inside of the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
    }

    for (int ii=0;ii<2;ii++) for (int jj=0;jj<2;jj++) for (int kk=0;kk<2;kk++) {
      double t[3];

      t[0]=(ii==0) ? xmin[0] : xmax[0];
      t[1]=(jj==0) ? xmin[1] : xmax[1];
      t[2]=(kk==0) ? xmin[2] : xmax[2];

      cout << "node (" << ii << "," << jj << "," << kk << ") x=" << t[0] << ", " << t[1] << ", " << t[2] << ": ";

      if (CutCell::CheckPointInsideDomain(t,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
        cout << "Outside of the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
      }
      else {
        cout << "Inside the domain (" << __FILE__ << "@" << __LINE__ << ")" << endl;
      }

    }

    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif




  static long int nCall=0;


  nCall++;


/*  if ((nCall==117771926)||(ptr==72941)) {
    cout << __FILE__ << "@" << __LINE__ << endl;
  }*/



//===================== END DEBUG ==================



  while (MovingTimeFinished==false) {
MovingLoop:



//    DEBUG
/*{
  double R,R1;

  R=sqrt(pow(xInit[1],2)+pow(xInit[2],2));
  R1=R;
}*/


//===================== DEBUG ==================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
   //check the consistency of the particle mover
   int iTemp,jTemp,kTemp;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    if (PIC::Mesh::mesh.fingCellIndex(xInit,iTemp,jTemp,kTemp,startNode,false)==-1) {
      exit(__LINE__,__FILE__,"Error: the cell is not found");
    }
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

    {
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      if (PIC::Mesh::mesh.fingCellIndex(r,iTemp,jTemp,kTemp,startNode,false)==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }

#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif

    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif
    //===================== END DEBUG ==================


    memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
    memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));

    MovingTimeFinished=true;
    dtMin=dtTotal;
    InternalBoundaryDescriptor_dtMin=NULL;
    ParticleIntersectionCode=_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_;

#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
    double acclInit[3],acclMiddle[3],vv;

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclInit,spec,ptr,xInit,vInit,startNode);

    a=sqrt(acclInit[0]*acclInit[0]+acclInit[1]*acclInit[1]+acclInit[2]*acclInit[2]);
    vv=sqrt(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2]);
    if (a*dtMin>2.0*vv) dtMin=2.0*vv/a;
#endif


    //Integrate the equations of motion
    //predictor
    dtMinInit=dtMin;
    dtMinInit2=dtMinInit/2.0;
    xMiddle[0]=xInit[0]+dtMinInit2*vInit[0];
    xMiddle[1]=xInit[1]+dtMinInit2*vInit[1];
    xMiddle[2]=xInit[2]+dtMinInit2*vInit[2];

    vMiddle[0]=vInit[0]+dtMinInit2*acclInit[0];
    vMiddle[1]=vInit[1]+dtMinInit2*acclInit[1];
    vMiddle[2]=vInit[2]+dtMinInit2*acclInit[2];

    //corrector
    //in the case when a symmetry has to be considered, transfer the particle position and velcoty accordinly
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    middleNode=PIC::Mesh::mesh.findTreeNode(xMiddle,startNode);

    if (middleNode==NULL) {
      int idim,nface,nIntersectionFace=-1;
      double cx,cv,r0[3],dt,dtIntersection=-1.0;

      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      if (ProcessOutsideDomainParticles!=NULL) {
        //determine through which face the particle left the domain

        for (nface=0;nface<6;nface++) {
          for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
            r0[idim]=xInit[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
            cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
            cv+=vInit[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
          }

          if (cv>0.0) {
            dt=-cx/cv;

            if ((dtIntersection<0.0)||(dt<dtIntersection)) {
              double cE0=0.0,cE1=0.0;

              for (idim=0;idim<3;idim++) {
                c=r0[idim]+dt*vInit[idim];

                cE0+=c*ExternalBoundaryFaceTable[nface].e0[idim],cE1+=c*ExternalBoundaryFaceTable[nface].e1[idim];
              }

              if ((cE0<-PIC::Mesh::mesh.EPS)||(cE0>ExternalBoundaryFaceTable[nface].lE0+PIC::Mesh::mesh.EPS) || (cE1<-PIC::Mesh::mesh.EPS)||(cE1>ExternalBoundaryFaceTable[nface].lE1+PIC::Mesh::mesh.EPS)) continue;

              nIntersectionFace=nface,dtIntersection=dt;
            }
          }
        }

        if (nIntersectionFace==-1) exit(__LINE__,__FILE__,"Error: cannot find the face of the intersection");

        for (idim=0;idim<3;idim++) {
          xInit[idim]+=dtIntersection*vInit[idim]-ExternalBoundaryFaceTable[nface].norm[idim]*PIC::Mesh::mesh.EPS;
          vInit[idim]+=dtIntersection*acclInit[idim];
        } 

        startNode=PIC::Mesh::mesh.findTreeNode(xInit,startNode);
        if (startNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the node");

        code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,startNode);
        memcpy(vFinal,vInit,3*sizeof(double));
        memcpy(xFinal,xInit,3*sizeof(double));
      }

      if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }
      else if (code==_PARTICLE_REJECTED_ON_THE_FACE_) {
        dtMin=dtIntersection;
        dtTotal-=dtMin;
        newNode=startNode;
        goto ProcessPhotoChemistry;
      }
      else  exit(__LINE__,__FILE__,"Error: not implemented");
    }

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
    exit(__LINE__,__FILE__,"not implemented");
#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif




    //Calculate the time of flight to the nearest block's face
#if _PIC__PARTICLE_MOVER__CHECK_BLOCK_FACE_INTERSECTION__MODE_ == _PIC_MODE_ON_

#if DIM == 1

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_

    !!!!! NOT ADJESTED FOR THE NEW PROCEDURE !!!!!

    if (fabs(v[0])>0.0) {
      dtTemp=(xminBlock[0]-x[0])/v[0];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=(xmaxBlock[0]-x[0])/v[0];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_

    !!!!! NOT ADJESTED FOR THE NEW PROCEDURE !!!!!

    //nface=0 -> rmin
    a=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    b=2.0*v[0]*x[0];
    c=x[0]*x[0]-xminBlock[0]*xminBlock[0];
    d=b*b-4.0*a*c;

    if (d>0.0) {
      if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

      sqrt_d=sqrt(d);
      dtTemp=-(b+sqrt_d)/(2.0*a);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=-2.0*c/(b+sqrt_d);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=0,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }

    //nface=1 -> rmax
    c=x[0]*x[0]-xmaxBlock[0]*xmaxBlock[0];
    d=b*b-4.0*a*c;

    if (d>0.0) {
      if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

      sqrt_d=sqrt(d);
      dtTemp=-(b+sqrt_d)/(2.0*a);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      dtTemp=-2.0*c/(b+sqrt_d);
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,nface_dtMin=1,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }

#else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

#elif DIM == 2
    exit(__LINE__,__FILE__,"not implemented");
#elif DIM == 3

    for (idim=0;idim<DIM;idim++) if (fabs(vMiddle[idim])>0.0) {
      //nface=0,2,4
      dtTemp=(xminBlock[idim]-xInit[idim])/vMiddle[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;

      //nface=1,3,5
      dtTemp=(xmaxBlock[idim]-xInit[idim])/vMiddle[idim];
      if ((0.0<dtTemp)&&(dtTemp<dtMin)) dtMin=dtTemp,ParticleIntersectionCode=_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
    }

#else
    exit(__LINE__,__FILE__,"Error: unknown value of DIM");
#endif
#endif


    //Calculate the time of flight to the nearest internal surface
    if (FirstBoundaryFlag==false) for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {
      if (InternalBoundaryDescriptor==lastInternalBoundaryDescriptor) continue;

      switch (InternalBoundaryDescriptor->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_: case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:

#if DIM == 3
        Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#elif DIM == 2
        exit(__LINE__,__FILE__,"not implemented");
#else
        Sphere1D=(cInternalSphere1DData*)(InternalBoundaryDescriptor->BoundaryElement);
        Sphere1D->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#endif

        dx=xInit[0]-x0Sphere[0],dy=xInit[1]-x0Sphere[1],dz=xInit[2]-x0Sphere[2];
        a=vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2];
        b=2.0*(vInit[0]*dx+vInit[1]*dy+vInit[2]*dz);
        c=dx*dx+dy*dy+dz*dz-radiusSphere*radiusSphere;

        if (c<0.0) {
          //the particle is inside the sphese
          //1. project the particle on the surface of the spehre
          //2. apply boundary conditions

          double l=0.0;
          int code;

          l=sqrt(dx*dx+dy*dy+dz*dz);
          l=(radiusSphere+PIC::Mesh::mesh.EPS)/l;

          xInit[0]=x0Sphere[0]+l*dx;
          xInit[1]=x0Sphere[1]+l*dy;
          xInit[2]=x0Sphere[2]+l*dz;
          startNode=PIC::Mesh::mesh.findTreeNode(xInit,startNode);

          FirstBoundaryFlag=true;
#if DIM == 3
          code=Sphere->ParticleSphereInteraction(spec,ptr,xInit,vInit,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#elif DIM == 2
        exit(__LINE__,__FILE__,"not implemented");
#else
        code=Sphere1D->ParticleSphereInteraction(spec,ptr,xInit,vInit,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#endif

          if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
            PIC::ParticleBuffer::DeleteParticle(ptr);
            return _PARTICLE_LEFT_THE_DOMAIN_;
          }

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

        dt1=-2.0*c/(b+sqrt_d);
        if ((dtTemp<0.0)||((dt1>0.0)&&(dt1<dtTemp))) {
          dtTemp=dt1;
        }

        if ((0.0<dtTemp)&&(dtTemp<dtMin)) {
          dtMin=dtTemp,InternalBoundaryDescriptor_dtMin=InternalBoundaryDescriptor;
          ParticleIntersectionCode=_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
        }

        break;
      case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
        Nucleus=(cInternalRotationBodyData*)(InternalBoundaryDescriptor->BoundaryElement);
	Nucleus->GetSphereGeometricalParameters(x0Nucleus,lNucleus,xmin,xmax);


	exit(__LINE__,__FILE__,"calculation of the time of flight is not implemented");


        Nucleus->SurfaceCurve(rSurface,xInit[0]);
	if(xmin<=xInit[0] && xmax>=xInit[0] && rSurface>sqrt(xInit[1]*xInit[1]+xInit[2]*xInit[2])) {
	  //the particle is inside the nucleus                                                                                                                  
	  PIC::ParticleBuffer::DeleteParticle(ptr);
          return _PARTICLE_LEFT_THE_DOMAIN_;
        }

        break;
      case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
    	  //do nothing
    	break;
      default:
        exit(__LINE__,__FILE__,"Error: undetermined internal boundary type");
      }
    }

    //check intersection of the particle trajectory with the cut-faces
//    cTriangleFace *IntersectionFace=NULL;

    if (startNode->FirstTriangleCutFace!=NULL) {
      CutCell::cTriangleFaceDescriptor *t;
      double TimeOfFlight;

      for (t=startNode->FirstTriangleCutFace;t!=NULL;t=t->next) if (t->TriangleFace!=lastIntersectedTriangleFace) {
        if (t->TriangleFace->RayIntersection(xInit,vInit,TimeOfFlight,PIC::Mesh::mesh.EPS)==true) {
          if (TimeOfFlight<dtMin) {
            dtMin=TimeOfFlight;
            IntersectionFace=t->TriangleFace;
            ParticleIntersectionCode=_BOUNDARY_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_,MovingTimeFinished=false;
          }
        }
      }
    }

    //adjust the particle moving time
    dtTotal-=dtMin;

    if (ParticleIntersectionCode!=_BOUNDARY_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) lastIntersectedTriangleFace=NULL;

    //advance the particle's position and velocity
    //interaction with the faces of the block and internal surfaces
    if (ParticleIntersectionCode==_BOUNDARY_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
/*      xFinal[0]=xInit[0]+dtMin*vInit[0]+PIC::Mesh::mesh.EPS*IntersectionFace->ExternalNormal[0];
      xFinal[1]=xInit[1]+dtMin*vInit[1]+PIC::Mesh::mesh.EPS*IntersectionFace->ExternalNormal[1];
      xFinal[2]=xInit[2]+dtMin*vInit[2]+PIC::Mesh::mesh.EPS*IntersectionFace->ExternalNormal[2];*/

//      if (dtTotal==0.0) dtTotal=1.0E-7*dtMin; //make sure that the particle doe not stay on the boundary face at the end of the motion

      //adjust 'dtMin' such that the
      double x0Face[3],FaceNorm[3];
      bool ExitFlag=true;
      int code;

      memcpy(x0Face,IntersectionFace->x0Face,3*sizeof(double));
      memcpy(FaceNorm,IntersectionFace->ExternalNormal,3*sizeof(double));

      double r0=(xInit[0]-x0Face[0])*FaceNorm[0] + (xInit[1]-x0Face[1])*FaceNorm[1] + (xInit[2]-x0Face[2])*FaceNorm[2];

      if (r0>PIC::Mesh::mesh.EPS) {
        do {
          ExitFlag=true;

          xFinal[0]=xInit[0]+dtMin*vInit[0];
          xFinal[1]=xInit[1]+dtMin*vInit[1];
          xFinal[2]=xInit[2]+dtMin*vInit[2];

          if ( ((xFinal[0]-x0Face[0])*FaceNorm[0] + (xFinal[1]-x0Face[1])*FaceNorm[1] + (xFinal[2]-x0Face[2])*FaceNorm[2]) < PIC::Mesh::mesh.EPS) {
            dtMin*=(1.0-1.0E-3);
            ExitFlag=false;
          }
        }
        while (ExitFlag==false);
      }
      else if (r0<0.0) {
        if (r0>-PIC::Mesh::mesh.EPS) {
          for (int idim=0;idim<DIM;idim++) xInit[idim]+=(0.0*PIC::Mesh::mesh.EPS-r0)*FaceNorm[idim];

          memcpy(xFinal,xInit,3*sizeof(double));
          dtMin=0.0;
        }
        else {
          double dtRecalculated,xIntersection[3],xIntersectionLocal[2];
          bool IntersectionCode;

          IntersectionCode=IntersectionFace->RayIntersection(xInit,vInit,dtRecalculated,PIC::Mesh::mesh.EPS);

          if (IntersectionCode==true) {
            for (int idim=0;idim<3;idim++) xIntersection[idim]=xInit[idim]+dtRecalculated*vInit[idim];

            IntersectionFace->GetProjectedLocalCoordinates(xIntersectionLocal,xIntersection);

            printf("$PREFIX: ptr=%ld, r0=%e, dtMin=%e, dtRecalculated=%e (%i,%s)\n",ptr,r0,dtMin,dtRecalculated,__LINE__,__FILE__);
            printf("$PREFIX: The \"global\" coordinates of the intersection: xIntersection[]=(%e,%e,%e)  (%i,%s)\n",xIntersection[0],xIntersection[1],xIntersection[2],__LINE__,__FILE__);
            printf("$PREFIX: The \"local\" coordinates of the intersection: xIntersectionLocal[]=(%e,%e)  (%i,%s)\n",xIntersectionLocal[0],xIntersectionLocal[1],__LINE__,__FILE__);
            printf("$PREFIX: e0Length*xLocal[0]=%e, e1Length*xLocal[1]=%e, EPS=%e (%i,%s)\n",IntersectionFace->e0Length*xIntersectionLocal[0],IntersectionFace->e1Length*xIntersectionLocal[1],PIC::Mesh::mesh.EPS,__LINE__,__FILE__);
            printf("$PREFIX: e0Length=%e, e1Length=%e (%i,%s)\n",IntersectionFace->e0Length,IntersectionFace->e1Length,__LINE__,__FILE__);
          }
          else {
            printf("$PREFIX: the intersection is not found (%i,%s)\n",__LINE__,__FILE__);
          }

          if (CutCell::CheckPointInsideDomain(xInit,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
             exit(__LINE__,__FILE__,"The point is outside of the domain");
          }

          exit(__LINE__,__FILE__,"error: the point is inside the body");
        }
      }
      else {
        memcpy(xFinal,xInit,3*sizeof(double));
        dtMin=0.0;
      }

/*      xFinal[0]=xInit[0]+dtMin*vInit[0];
      xFinal[1]=xInit[1]+dtMin*vInit[1];
      xFinal[2]=xInit[2]+dtMin*vInit[2];*/


      vFinal[0]=vInit[0]+dtMin*acclInit[0];
      vFinal[1]=vInit[1]+dtMin*acclInit[1];
      vFinal[2]=vInit[2]+dtMin*acclInit[2];

      lastIntersectedTriangleFace=IntersectionFace;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


      code=(ProcessTriangleCutFaceIntersection!=NULL) ? ProcessTriangleCutFaceIntersection(ptr,xFinal,vFinal,IntersectionFace) : _PARTICLE_DELETED_ON_THE_FACE_;


/*      double c=vFinal[0]*IntersectionFace->ExternalNormal[0]+vFinal[1]*IntersectionFace->ExternalNormal[1]+vFinal[2]*IntersectionFace->ExternalNormal[2];
      vFinal[0]-=2.0*c*IntersectionFace->ExternalNormal[0];
      vFinal[1]-=2.0*c*IntersectionFace->ExternalNormal[1];
      vFinal[2]-=2.0*c*IntersectionFace->ExternalNormal[2];*/

      if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

    }
    else if (startNode->FirstTriangleCutFace!=NULL) {
    	//use the first order integration if 'startNode' contains cut-faces, but no intersection with them is determened for the 1st order scheme

        xFinal[0]=xInit[0]+dtMin*vInit[0];
        xFinal[1]=xInit[1]+dtMin*vInit[1];
        xFinal[2]=xInit[2]+dtMin*vInit[2];

        vFinal[0]=vInit[0]+dtMin*acclInit[0];
        vFinal[1]=vInit[1]+dtMin*acclInit[1];
        vFinal[2]=vInit[2]+dtMin*acclInit[2];

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif

    }



    else if (ParticleIntersectionCode==_INTERNAL_SPHERE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      int code;

      xFinal[0]=xInit[0]+dtMin*vInit[0];
      xFinal[1]=xInit[1]+dtMin*vInit[1];
      xFinal[2]=xInit[2]+dtMin*vInit[2];

      vFinal[0]=vInit[0]+dtMin*acclInit[0];
      vFinal[1]=vInit[1]+dtMin*acclInit[1];
      vFinal[2]=vInit[2]+dtMin*acclInit[2];

      FirstBoundaryFlag=true;

      lastInternalBoundaryDescriptor=InternalBoundaryDescriptor_dtMin;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2]);
      newNode=PIC::Mesh::mesh.findTreeNode(r,middleNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif

#if DIM == 3
      code=((cInternalSphericalData*)(InternalBoundaryDescriptor_dtMin->BoundaryElement))->ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,InternalBoundaryDescriptor_dtMin->BoundaryElement);
#elif DIM == 2
      exit(__LINE__,__FILE__,"not implemented");
#else
      code=((cInternalSphere1DData*)(InternalBoundaryDescriptor_dtMin->BoundaryElement))->ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,InternalBoundaryDescriptor_dtMin->BoundaryElement);
#endif


      if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

/*
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      double r[3]={0.0,0.0,0.0};

      r[0]=sqrt(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2]);
      newNode=PIC::Mesh::mesh.findTreeNode(r,middleNode);
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif
*/
    }
    else if (ParticleIntersectionCode==_BLOCK_FACE_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {

      xFinal[0]=xInit[0]+dtMin*vMiddle[0];
      xFinal[1]=xInit[1]+dtMin*vMiddle[1];
      xFinal[2]=xInit[2]+dtMin*vMiddle[2];

      vFinal[0]=vInit[0]+dtMin*acclMiddle[0];
      vFinal[1]=vInit[1]+dtMin*acclMiddle[1];
      vFinal[2]=vInit[2]+dtMin*acclMiddle[2];

      FirstBoundaryFlag=false;

      //check if the particle is outside ofthe block - force the particle outside of the block if its needed
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      double EPS=PIC::Mesh::mesh.EPS;

      for (idim=0;idim<DIM;idim++) {
        if (xFinal[idim]<xminBlock[idim]+EPS) xFinal[idim]=xminBlock[idim]-EPS;
        if (xFinal[idim]>xmaxBlock[idim]-EPS) xFinal[idim]=xmaxBlock[idim]+EPS;

      }

      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
exit(__LINE__,__FILE__,"Error: not implemented");
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif


      //reserve the place for particle's cloning:
      //if a particle crossed a face, the time step and particle weight are changed
      //adjust the value of the dtLeft to match the time step for the species 'spec'
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    //do nothing
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
      double WeightCorrectionFactor,TimeStepRatio;

      if (newNode!=NULL) {
        TimeStepRatio=newNode->block->GetLocalTimeStep(spec)/startNode->block->GetLocalTimeStep(spec);
        WeightCorrectionFactor=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)*startNode->block->GetLocalParticleWeight(spec)*TimeStepRatio/newNode->block->GetLocalParticleWeight(spec);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,ParticleData);

        dtTotal*=TimeStepRatio;
      }

      #else
      exit(__LINE__,__FILE__,"Error: option is not recognized");
#endif

    }
    else if (ParticleIntersectionCode==_UNDEFINED_MIN_DT_INTERSECTION_CODE_UTSNFTT_) {
      xFinal[0]=xInit[0]+dtMin*vMiddle[0];
      xFinal[1]=xInit[1]+dtMin*vMiddle[1];
      xFinal[2]=xInit[2]+dtMin*vMiddle[2];

      vFinal[0]=vInit[0]+dtMin*acclMiddle[0];
      vFinal[1]=vInit[1]+dtMin*acclMiddle[1];
      vFinal[2]=vInit[2]+dtMin*acclMiddle[2];

      FirstBoundaryFlag=false;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
exit(__LINE__,__FILE__,"not implemented");
#else
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
#endif

    }
    else {
      exit(__LINE__,__FILE__,"Error: the option is nor defined");
    }

    if (newNode==NULL) {
       int idim,nface,nIntersectionFace=-1;
       double cx,cv,r0[3],dt,dtIntersection=-1.0;

       //the particle left the computational domain
       int code=_PARTICLE_DELETED_ON_THE_FACE_;

       //call the function that process particles that leaved the coputational domain
       if (ProcessOutsideDomainParticles!=NULL) {
         //determine through which face the particle left the domain

         for (nface=0;nface<6;nface++) {
           for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
             r0[idim]=xInit[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
             cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
             cv+=vMiddle[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
           }

           if (cv>0.0) {
             dt=-cx/cv;

             if ((dtIntersection<0.0)||(dt<dtIntersection)) {
               double cE0=0.0,cE1=0.0;

               for (idim=0;idim<3;idim++) {
                 c=r0[idim]+dt*vMiddle[idim];

                 cE0+=c*ExternalBoundaryFaceTable[nface].e0[idim],cE1+=c*ExternalBoundaryFaceTable[nface].e1[idim];
               }

               if ((cE0<-PIC::Mesh::mesh.EPS)||(cE0>ExternalBoundaryFaceTable[nface].lE0+PIC::Mesh::mesh.EPS) || (cE1<-PIC::Mesh::mesh.EPS)||(cE1>ExternalBoundaryFaceTable[nface].lE1+PIC::Mesh::mesh.EPS)) continue;

               nIntersectionFace=nface,dtIntersection=dt;
             }
           }
         }

         if (nIntersectionFace==-1) exit(__LINE__,__FILE__,"Error: cannot find the face of the intersection");

         for (idim=0;idim<3;idim++) {
           xInit[idim]+=dtIntersection*vMiddle[idim]-ExternalBoundaryFaceTable[nface].norm[idim]*PIC::Mesh::mesh.EPS;
           vInit[idim]+=dtIntersection*acclMiddle[idim];
         } 

         newNode=PIC::Mesh::mesh.findTreeNode(xInit,middleNode);
         if (newNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the node");

         code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,newNode);
         memcpy(vFinal,vInit,3*sizeof(double));
         memcpy(xFinal,xInit,3*sizeof(double));
       }

       if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
         PIC::ParticleBuffer::DeleteParticle(ptr);
         return _PARTICLE_LEFT_THE_DOMAIN_;
       }
       else if (code==_PARTICLE_REJECTED_ON_THE_FACE_) {
         dtTotal+=dtMin-dtIntersection;
         dtMin=dtIntersection;
       }
       else  exit(__LINE__,__FILE__,"Error: not implemented");
     }


    //check the possible photolytic reactions
ProcessPhotoChemistry:
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    //model the photolytic transformation
    if (PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(xFinal,ptr,spec,dtMin,newNode)==_PHOTOLYTIC_REACTION_OCCURES_) {
      int PhotolyticReactionsReturnCode,specInit=spec;

      //PhotolyticReactionsReturnCode=PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable[specInit](xInit,xFinal,ptr,spec,ParticleData);
      PhotolyticReactionsReturnCode=_PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(xInit,xFinal,ptr,spec,ParticleData);

      //adjust the value of the dtLeft to match the time step for the species 'spec'
      if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_) {
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }
      else if (PhotolyticReactionsReturnCode==_PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_) {
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        dtTotal*=newNode->block->GetLocalTimeStep(spec)/newNode->block->GetLocalTimeStep(specInit);

        //check the probability for the new-species particle tostay in the system
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        double Correction,Rate;

        Rate=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)*newNode->block->GetLocalParticleWeight(specInit)/newNode->block->GetLocalTimeStep(specInit);
        Correction=Rate*newNode->block->GetLocalTimeStep(spec)/newNode->block->GetLocalParticleWeight(spec);

        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(Correction,ParticleData);
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
        exit(__LINE__,__FILE__,"Accounting for the possible difference in the time steps and weights have not been inplemented when the individual particle weight corraction factor is off");

#else
     exit(__LINE__,__FILE__,"The option is unknown");
#endif

        #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
        #if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__CHEMISTRY_ == _PIC_MODE_ON_
        PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData);
        #endif
        #endif
      }
      else (__LINE__,__FILE__,"Error: the option is unknown");
    }
#elif _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
    //model the generic particle transformation
    int GenericParticleTransformationReturnCode,specInit=spec;

    GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode);   //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

    if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //adjust the value of the dtLeft to match the time step for the species 'spec'
    if (spec!=specInit) {
      dtTotal*=newNode->block->GetLocalTimeStep(spec)/newNode->block->GetLocalTimeStep(specInit);

      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      #if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__CHEMISTRY_ == _PIC_MODE_ON_
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData);
      #endif
      #endif
    }
#endif





    //adjust the value of 'startNode'
    startNode=newNode;
    memcpy(vInit,vFinal,3*sizeof(double));
    memcpy(xInit,xFinal,3*sizeof(double));

    //save the trajectory point
    #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::RecordTrajectoryPoint(xInit,vInit,spec,ParticleData);
    #endif


  }

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  #if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData);
  #endif
  #endif


  //check if the particle is outside of the internal surfaces. In a case if the particle is inside an internal surface -> correct its position and exdcute the boundary condition procedure
  for (InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;InternalBoundaryDescriptor!=NULL;InternalBoundaryDescriptor=InternalBoundaryDescriptor->nextInternalBCelement) {

    switch (InternalBoundaryDescriptor->BondaryType) {
    case _INTERNAL_BOUNDARY_TYPE_SPHERE_: case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:

#if DIM == 3
      Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);
      Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#elif DIM == 2
      exit(__LINE__,__FILE__,"not implemented");
#else
      Sphere1D=(cInternalSphere1DData*)(InternalBoundaryDescriptor->BoundaryElement);
      Sphere1D->GetSphereGeometricalParameters(x0Sphere,radiusSphere);
#endif

      dx=xFinal[0]-x0Sphere[0],dy=xFinal[1]-x0Sphere[1],dz=xFinal[2]-x0Sphere[2];
      c=dx*dx+dy*dy+dz*dz-radiusSphere*radiusSphere;

      if (c<0.0) {
        //the particle is inside the sphese
        //1. project the particle on the surface of the spehre
        //2. apply boundary conditions

        double l=0.0;
        int code;

        l=sqrt(dx*dx+dy*dy+dz*dz);
        l=(radiusSphere+PIC::Mesh::mesh.EPS)/l;

        xFinal[0]=x0Sphere[0]+l*dx;
        xFinal[1]=x0Sphere[1]+l*dy;
        xFinal[2]=x0Sphere[2]+l*dz;
        startNode=PIC::Mesh::mesh.findTreeNode(xFinal,startNode);

        FirstBoundaryFlag=true;
#if DIM == 3
        code=Sphere->ParticleSphereInteraction(spec,ptr,xInit,vInit,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#elif DIM == 2
        exit(__LINE__,__FILE__,"not implemented");
#else
        code=Sphere1D->ParticleSphereInteraction(spec,ptr,xInit,vInit,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);
#endif

        if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
          PIC::ParticleBuffer::DeleteParticle(ptr);
          return _PARTICLE_LEFT_THE_DOMAIN_;
        }
      }
      break;
    case _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_:
      Nucleus=(cInternalRotationBodyData*)(InternalBoundaryDescriptor->BoundaryElement);
      Nucleus->GetSphereGeometricalParameters(x0Nucleus,lNucleus,xmin,xmax);

      Nucleus->SurfaceCurve(rSurface,xFinal[0]);
      if(xmin<=xFinal[0] && xmax>=xFinal[0] && rSurface>sqrt(xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])) {
	//the particle is inside the nucleus                                                                                                                     
	PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

      break;
    case _INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_:
#if  _PIC_CONTROL_PARTICLE_INSIDE_NASTRAN_SURFACE_ == _PIC_MODE_ON_
        if (CutCell::CheckPointInsideDomain(xFinal,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
            exit(__LINE__,__FILE__,"The point is outside of the domain");
        }
#endif
    	break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not recognized");
    }
  }



  //Rotate particle position and velocity when symmetry is accounted
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
  //do nothing
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
  double r,l,v1[3],cosTz,sinTz,cosTy,sinTy;

  r=sqrt(pow(xFinal[0],2)+pow(xFinal[1],2)+pow(xFinal[2],2));

  if (r>1.0E-20) {
    double xfinal[3];

    xfinal[0]=xFinal[0]/r,xfinal[1]=xFinal[1]/r,xfinal[2]=xFinal[1]/r;
    l=sqrt(pow(xfinal[0],2)+pow(xfinal[1],2));

    if (l>1.0E-20) {
      cosTz=xfinal[0]/l,sinTz=xfinal[1]/l;
      cosTy=l,sinTy=xfinal[2];
    }
    else cosTz=1.0,sinTz=0.0,sinTy=xfinal[2],cosTy=0.0;

    v1[0]=cosTy*cosTz*vFinal[0]+cosTy*sinTz*vFinal[1]+sinTy*vFinal[2];
    v1[1]=-sinTz*vFinal[0]+cosTz*vFinal[1];
    v1[2]=-sinTy*cosTz*vFinal[0]-sinTy*sinTz*vFinal[1]+cosTy*vFinal[2];

    vFinal[0]=v1[0],vFinal[1]=v1[1],vFinal[2]=v1[2];
    xFinal[0]=r,xFinal[1]=0.0,xFinal[2]=0.0;
  }
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif


  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");


  //check if the block is allocated
  if (startNode->block==NULL) {
    cout << "$PREFIX: Most probably the time step is too large. Error at " << __FILE__ << "@" << __LINE__  << endl;
    cout << "$PREFIX: startNode->block==NULL" << endl;
    cout << "$PREFIX: ThisThread=" << PIC::ThisThread << ", startNode->Thread=" << startNode->Thread << endl;
    cout << "$PREFIX: x=" << xFinal[0] << ", " << xFinal[1] << ", " << xFinal[2] << endl;
    cout << "$PREFIX: v=" << vFinal[0] << ", " << vFinal[1] << ", " << vFinal[2] << endl;
    cout << "$PREFIX: startNode->xmin=" << startNode->xmin[0] << ", " << startNode->xmin[1] << ", " << startNode->xmin[2] << endl;
    cout << "$PREFIX: startNode->xmax=" << startNode->xmax[0] << ", " << startNode->xmax[1] << ", " << startNode->xmax[2] << endl;
  }


  PIC::Mesh::cDataBlockAMR *block=startNode->block;
  long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;



  //=====================  DEBUG =========================
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_

  if (PIC::ParticleBuffer::IsParticleAllocated(ParticleData)==false) {
	  exit(__LINE__,__FILE__,"Error: an unallocated particle is intercepted");
  }



  cell=startNode->block->GetCenterNode(LocalCellNumber);

  if ((cell->Measure<=0.0)&&(startNode->Thread==PIC::Mesh::mesh.ThisThread)) {
    cout << "$PREFIX:" << __FILE__<< __LINE__ << endl;


    cout << "$PREFIX:Error: cell has zero volume (" <<__FILE__<< "@" << __LINE__ << ")" << endl;

     double r,rprobe[3]={0.0,0.0,0.0};
     int di,dj,dk;


     cout << "$PREFIX:x particle=";
     for (r=0.0,idim=0;idim<DIM;idim++) {
       r+=pow(xFinal[idim],2);
       cout << xFinal[idim] << " ";
     }

     cout << ", |x|= " << sqrt(r) << endl;

     for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
       startNode->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

       for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
       cout << "$PREFIX:Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << sqrt(r) << endl;
     }


     double  vol=-1,xmin[3],xmax[3];

     if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");

     xmin[0]=startNode->xmin[0]+i*(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_;
     xmin[1]=startNode->xmin[1]+j*(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_;
     xmin[2]=startNode->xmin[2]+k*(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_;

     xmax[0]=startNode->xmin[0]+(i+1)*(startNode->xmax[0]-startNode->xmin[0])/_BLOCK_CELLS_X_;
     xmax[1]=startNode->xmin[1]+(j+1)*(startNode->xmax[1]-startNode->xmin[1])/_BLOCK_CELLS_Y_;
     xmax[2]=startNode->xmin[2]+(k+1)*(startNode->xmax[2]-startNode->xmin[2])/_BLOCK_CELLS_Z_;

     vol=CutCell::GetRemainedBlockVolume(xmin,xmax,PIC::Mesh::mesh.EPS,1.0E-2,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,startNode->FirstTriangleCutFace);

     if (CutCell::CheckPointInsideDomain(xFinal,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0*PIC::Mesh::mesh.EPS)==false) {
         exit(__LINE__,__FILE__,"The point is outside of the domain");
     }

    exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
  }
#endif
  //===================   END DEBUG ==============================



  return _PARTICLE_MOTION_FINISHED_;
}

int PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode)  {
  return UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder(ptr,dtTotal,startNode,false);
}

