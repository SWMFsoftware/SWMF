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
int PIC::Mover::Relativistic::Boris(long int ptr,double dtTotalIn,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  PIC::ParticleBuffer::byte *ParticleData;
  double gamma,vInit[3],xInit[3]={0.0,0.0,0.0},vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
  double mass,uMinus[3],ElectricCharge,E[3],B[3],QdT_over_twoM;
  int idim,i,j,k,spec;

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

      for (idim=0;idim<3;idim++) {
        ExternalBoundaryFaceTable[nface].x0[idim]=(ExternalBoundaryFaceTable[nface].nX0[idim]==0) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim];

        cE0+=pow(((ExternalBoundaryFaceTable[nface].e0[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
        cE1+=pow(((ExternalBoundaryFaceTable[nface].e1[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh.rootTree->xmin[idim] : PIC::Mesh::mesh.rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
      }

      ExternalBoundaryFaceTable[nface].lE0=sqrt(cE0);
      ExternalBoundaryFaceTable[nface].lE1=sqrt(cE1);
    }
  }

  static long int nCall=0;
  nCall++;

  #if _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_
  double dtTotalInit=dtTotalIn;
  #endif

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

  ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
  mass=PIC::MolecularData::GetMass(spec);

  if (dtTotalIn==0.0) {
    memcpy(xFinal,xInit,3*sizeof(double));
    memcpy(vFinal,vInit,3*sizeof(double));
  }
  else while (dtTotalIn>0.0) {
    gamma=1.0/sqrt(1.0-(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2])/(SpeedOfLight*SpeedOfLight));

    memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
    memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));

    //calculate fields acting upon the particle
    PIC::CPLR::InitInterpolationStencil(xInit,startNode);
    PIC::CPLR::GetBackgroundElectricField(E);
    PIC::CPLR::GetBackgroundMagneticField(B);

    //limit the time integration period with the gyroqrecuency of the particle
    double dt,GyroFreq,dtMax;

    if (Vector3D::Length(B)>1.0E-25) {
      GyroFreq=::Relativistic::GetGyroFrequency(vInit,mass,ElectricCharge,B);
      dtMax=1.0/GyroFreq;
      dt=(dtMax<dtTotalIn) ? dtMax : dtTotalIn;
    }
    else dt=dtTotalIn;

    dtTotalIn-=dt;

    #if _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_ == _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__ENABLED_
    if  (BackwardTimeIntegrationMode==_PIC_MODE_ON_) {
      for (idim=0;idim<3;idim++) B[idim]=-B[idim],E[idim]=-E[idim];
    }
    #endif

    //convert velocity into momentum and advance particle half time step
    QdT_over_twoM=ElectricCharge*dt/(2.0*mass);
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

    #if _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_ == _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__ENABLED_
    switch (BackwardTimeIntegrationMode) {
    case _PIC_MODE_OFF_:
      for (idim=0;idim<3;idim++) {
        vFinal[idim]=uFinal[idim]/gamma;
        xFinal[idim]=xInit[idim]+vFinal[idim]*dt;
      }
      break;
    case _PIC_MODE_ON_:
      for (idim=0;idim<3;idim++) {
        vFinal[idim]=uFinal[idim]/gamma;
        xFinal[idim]=xInit[idim]-vFinal[idim]*dt;
      }
    }
    #else  //_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_
    for (idim=0;idim<3;idim++) {
      vFinal[idim]=uFinal[idim]/gamma;
      xFinal[idim]=xInit[idim]+vFinal[idim]*dt;
    }
    #endif //_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_

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
      code=ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dt,(void*)newNode,BoundaryElement);

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
      #if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_
      //do nothing -> the particle deleting code is already set
      #else
       //call the function that process particles that leaved the coputational domain
//       if (ProcessOutsideDomainParticles!=NULL) {
         //determine through which face the particle left the domain

      int nface,nIntersectionFace;
      double tVelocityIncrement,cx,cv,r0[3],dtEffective,vEffective[3]={0.5*(vInit[0]+vFinal[0]),0.5*(vInit[1]+vFinal[1]),0.5*(vInit[2]+vFinal[2])},c,dtIntersection=-1.0;

      switch (BackwardTimeIntegrationMode) {
      case _PIC_MODE_ON_:
        for (idim=0;idim<3;idim++) vEffective[idim]=xInit[idim]-xFinal[idim];
        break;
      default:
        for (idim=0;idim<3;idim++) vEffective[idim]=xFinal[idim]-xInit[idim];
      }

      for (nface=0;nface<6;nface++) {

       #if _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_ == _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__ENABLED_
        switch (BackwardTimeIntegrationMode) {
        case _PIC_MODE_ON_:
          for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
            r0[idim]=xFinal[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
            cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
            cv+=vEffective[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
          }

          dtEffective=(cv<0.0) ? -cx/cv : -1.0;

          break;
        default:
          for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
            r0[idim]=xInit[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
            cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
            cv+=vEffective[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
          }

          dtEffective=(cv>0.0) ? -cx/cv : -1.0;
        }
        #else  //_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_
        for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
          r0[idim]=xInit[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
          cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
          cv+=vEffective[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
        }

        dtEffective=(cv>0.0) ? -cx/cv : -1.0;
        #endif  //_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_

        if (dtEffective>0.0) {
          if ((dtIntersection<0.0)||(dtEffective<dtIntersection)&&(dtEffective>0.0)) {
            double cE0=0.0,cE1=0.0;

            for (idim=0;idim<3;idim++) {
              c=r0[idim]+dtEffective*vEffective[idim];

              cE0+=c*ExternalBoundaryFaceTable[nface].e0[idim],cE1+=c*ExternalBoundaryFaceTable[nface].e1[idim];
            }

            if ((cE0<-PIC::Mesh::mesh.EPS)||(cE0>ExternalBoundaryFaceTable[nface].lE0+PIC::Mesh::mesh.EPS) || (cE1<-PIC::Mesh::mesh.EPS)||(cE1>ExternalBoundaryFaceTable[nface].lE1+PIC::Mesh::mesh.EPS)) continue;

            nIntersectionFace=nface,dtIntersection=dtEffective;
          }
        }
      }

      if (nIntersectionFace==-1) exit(__LINE__,__FILE__,"Error: cannot find the face of the intersection");


      switch (BackwardTimeIntegrationMode) {
      case _PIC_MODE_ON_:
        for (idim=0;idim<3;idim++) {
          xInit[idim]=xFinal[idim]+dtIntersection*(xInit[idim]-xFinal[idim])-ExternalBoundaryFaceTable[nface].norm[idim]*PIC::Mesh::mesh.EPS;
          vInit[idim]=vFinal[idim]+dtIntersection*(vInit[idim]-vFinal[idim]);
        }
        break;
      default:
        for (idim=0;idim<3;idim++) {
          xInit[idim]+=dtIntersection*(xFinal[idim]-xInit[idim])-ExternalBoundaryFaceTable[nface].norm[idim]*PIC::Mesh::mesh.EPS;
          vInit[idim]+=dtIntersection*(vFinal[idim]-vInit[idim]);
        }
      }

      newNode=PIC::Mesh::mesh.findTreeNode(xInit,startNode);

      if (newNode==NULL) {
        //the partcle is outside of the domain -> correct particle location and determine the newNode;
        double xmin[3],xmax[3];
        int ii;

        memcpy(xmin,PIC::Mesh::mesh.xGlobalMin,3*sizeof(double));
        memcpy(xmax,PIC::Mesh::mesh.xGlobalMax,3*sizeof(double));

        for (ii=0;ii<3;ii++) {
          if (xmin[ii]>=xInit[ii]) xInit[ii]=xmin[ii]+PIC::Mesh::mesh.EPS;
          if (xmax[ii]<=xInit[ii]) xInit[ii]=xmax[ii]-PIC::Mesh::mesh.EPS;
        }

        newNode=PIC::Mesh::mesh.findTreeNode(xInit,startNode);

        if (newNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the node");
      }

      switch (_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_:
        code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,newNode);
        break;

      case  _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_:
        exit(__LINE__,__FILE__,"Error: not implemented");
        break;

      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_:
        //reflect the particle back into the domain
        {
          double c=0.0;
          for (int idim=0;idim<3;idim++) c+=ExternalBoundaryFaceTable[nIntersectionFace].norm[idim]*vInit[idim];
          for (int idim=0;idim<3;idim++) vInit[idim]-=2.0*c*ExternalBoundaryFaceTable[nIntersectionFace].norm[idim];
        }

        code=_PARTICLE_REJECTED_ON_THE_FACE_;
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }



      memcpy(vFinal,vInit,3*sizeof(double));
      memcpy(xFinal,xInit,3*sizeof(double));
#endif //_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE


      switch(code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }



    startNode=newNode;
    memcpy(xInit,xFinal,3*sizeof(double));
    memcpy(vInit,vFinal,3*sizeof(double));
  }  //end of the particle Trajectory Integration Loop
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

    GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xInit,xFinal,vFinal,spec,ptr,ParticleData,dtTotalInit,startNode);   //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

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

  if (PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

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





