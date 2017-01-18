//$Id$
//relativistic particle trajectory integration

/*
 * pic_mover_relativistic_boris.cpp
 *
 *  Created on: Sep 25, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"

//advance particle location
int PIC::Mover::Relativistic::Boris(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  double dtTemp;
  PIC::ParticleBuffer::byte *ParticleData;
  double gamma,vInit[3],xInit[3]={0.0,0.0,0.0},vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
  double mass,uMinus[3],pFinal[3],ElectricCharge,c,E[3],B[3],QdT_over_twoM;
  int idim;

  long int LocalCellNumber;
  int i,j,k,spec;

  PIC::Mesh::cDataCenterNode *cell;
  bool MovingTimeFinished=false;

  double xmin,xmax;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

  ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
  mass=PIC::MolecularData::GetMass(spec);
  gamma=1.0/sqrt(1.0-(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2])/(SpeedOfLight*SpeedOfLight));


  static long int nCall=0;
  nCall++;

  memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
  memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));

  //calculate fields acting upon the particle
  PIC::CPLR::InitInterpolationStencil(xInit,startNode);
  PIC::CPLR::GetBackgroundElectricField(E);
  PIC::CPLR::GetBackgroundMagneticField(B);

  //convert velocity into momentum and advance particle half time step
  QdT_over_twoM=ElectricCharge*dtTotal/(2.0*mass);
  for (idim=0;idim<3;idim++) uMinus[idim]=gamma*vInit[idim]+QdT_over_twoM*E[idim];

  //first rotation
  double t[3],s[3],uPrime[3],uPlus[3],l=0.0;

  gamma=sqrt(1.0+(uMinus[0]*uMinus[0]+uMinus[1]*uMinus[1]+uMinus[2]*uMinus[2])/(SpeedOfLight*SpeedOfLight));

  for (idim=0;idim<3;idim++) {
    t[idim]=QdT_over_twoM/gamma*B[idim];
    l+=pow(t[idim],2);
  }

  Vector3D::CrossProduct(uPrime,uMinus,t);
  for (idim=0;idim<3;idim++) uPrime[idim]+=uMinus[idim];

  //second rotation
  for (idim=0;idim<3;idim++) s[idim]=2.0*t[idim]/(1.0+l);

  Vector3D::CrossProduct(uPlus,uPrime,s);
  for (idim=0;idim<3;idim++) uPlus[idim]+=uMinus[idim];


  //second half-time electric field acceleration
  double uFinal[3];

  for (idim=0;idim<3;idim++) uFinal[idim]=uPlus[idim]+QdT_over_twoM*E[idim];
  gamma=sqrt(1.0+(uFinal[0]*uFinal[0]+uFinal[1]*uFinal[1]+uFinal[2]*uFinal[2])/(SpeedOfLight*SpeedOfLight));

  for (idim=0;idim<3;idim++) {
    vFinal[idim]=uFinal[idim]/gamma;
    xFinal[idim]=xInit[idim]+vFinal[idim]*dtTotal;
  }

  //interaction with the faces of the block and internal surfaces
  //check whether the particle trajectory is intersected the spherical body
#if _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_
  double rFinal2;
  static double rSphere=max(_RADIUS_(_TARGET_),Exosphere::Planet->Radius);

  //if the particle is inside the sphere -> apply the boundary condition procedure
  if ((rFinal2=xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])<rSphere*rSphere) {
    double r=sqrt(rFinal2);
    int code;

    static cInternalSphericalData_UserDefined::fParticleSphereInteraction ParticleSphereInteraction=
        ((cInternalSphericalData*)(PIC::Mesh::mesh.InternalBoundaryList.front().BoundaryElement))->ParticleSphereInteraction;
    static void* BoundaryElement=PIC::Mesh::mesh.InternalBoundaryList.front().BoundaryElement;

    //move the particle location at the surface of the sphere
    for (int idim=0;idim<DIM;idim++) xFinal[idim]*=rSphere/r;

    //determine the block of the particle location
    newNode=PIC::Mesh::mesh.findTreeNode(xFinal,startNode);

    //apply the boundary condition
    code=ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,BoundaryElement);

    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
  }
  else {
    newNode=PIC::Mesh::mesh.findTreeNode(xFinal,startNode);
  }
#else
  newNode=PIC::Mesh::mesh.findTreeNode(xFinal,startNode);
#endif //_TARGET_ == _TARGET_NONE_


  if (newNode==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

    //call the function that process particles that leaved the coputational domain
    switch(code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }

  //save the trajectory point
 #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
   PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)newNode);

   #if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_
   PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)newNode);
   #endif
 #endif


#if _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
    //model the generic particle transformation
    int GenericParticleTransformationReturnCode,specInit=spec;

    GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xInit,xFinal,vFinal,spec,ptr,ParticleData,dtTotal,startNode);   //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

    if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //adjust the value of the dtLeft to match the time step for the species 'spec'
    if (spec!=specInit) {
      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      #if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__CHEMISTRY_ == _PIC_MODE_ON_
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)startNode);
      #endif
      #endif
    }
#endif //_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_




  //finish the trajectory integration procedure
  PIC::Mesh::cDataBlockAMR *block;
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
  PIC::Debugger::CatchOutLimitValue(vFinal,DIM,__LINE__,__FILE__);
  PIC::Debugger::CatchOutLimitValue(xFinal,DIM,__LINE__,__FILE__);
#endif
#endif

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,newNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

  if ((block=newNode->block)==NULL) {
    exit(__LINE__,__FILE__,"Error: the block is empty. Most probably hte tiime step is too long");
  }

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  tempFirstCellParticlePtr=block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  tempFirstCellParticlePtr=block->GetTempParticleMovingListTableThread(omp_get_thread_num(),i,j,k);
#else
#error The option is unknown
#endif

  tempFirstCellParticle=(*tempFirstCellParticlePtr);

  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  *tempFirstCellParticlePtr=ptr;

  return _PARTICLE_MOTION_FINISHED_;
}





