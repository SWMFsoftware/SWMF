//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the module that performs column integration

#include "pic.h"

PIC::ColumnIntegration::cBoundingBoxFace PIC::ColumnIntegration::BoundingBoxFace[6]={ {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0}, {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0},
    {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0}, {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0},
    {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0}, {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0}};

bool PIC::ColumnIntegration::InitializedFlag=false;

//====================================================
//init the module
void PIC::ColumnIntegration::Init() {
  double *xmin,*xmax,e0Length,e1Length;
  int nface,i;

  if (PIC::Mesh::mesh.rootTree==NULL) exit(__LINE__,__FILE__,"Error: the mesh should be initialized before this point");
  if (DIM!=3) exit(__LINE__,__FILE__,"Error: the model is implemented only for DIM=3");

  InitializedFlag=true;

  xmin=PIC::Mesh::mesh.xGlobalMin;
  xmax=PIC::Mesh::mesh.xGlobalMax;

  //set up the definitions of the bounding faces
  for (nface=0;nface<6;nface++) {
    e0Length=0.0,e1Length=0.0;

    for (i=0;i<3;i++) {
      BoundingBoxFace[nface].x[i]=(x0PlaneNodeIndex[nface][i]==0) ? xmin[i] : xmax[i];

      BoundingBoxFace[nface].e0[i]=((x1PlaneNodeIndex[nface][i]==0) ? xmin[i] : xmax[i]) - BoundingBoxFace[nface].x[i];
      BoundingBoxFace[nface].e1[i]=((x2PlaneNodeIndex[nface][i]==0) ? xmin[i] : xmax[i]) - BoundingBoxFace[nface].x[i];

      BoundingBoxFace[nface].Normal[i]=PlaneNormal[nface][i];
      e0Length+=pow(BoundingBoxFace[nface].e0[i],2);
      e1Length+=pow(BoundingBoxFace[nface].e1[i],2);
    }

    BoundingBoxFace[nface].e0Length=sqrt(e0Length);
    BoundingBoxFace[nface].e1Length=sqrt(e1Length);
  }
}


//====================================================
//find initial and final points of the column integration
bool PIC::ColumnIntegration::FindIntegrationLimits(double *x0,double *l,double& IntegrationPathLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xStartNode,double *xFinish,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xFinishNode) {
  double t,lPerp,c,c0,c1;
  int nface,idim;
  vector<double> IntersectionTime;

  xStartNode=NULL,xFinishNode=NULL;

  //initialize the model if needed
  if (InitializedFlag==false) Init();

  //determine the intersection time for all faces
  for (nface=0;nface<6;nface++) {
     for (idim=0,lPerp=0.0;idim<3;idim++) lPerp+=l[idim]*BoundingBoxFace[nface].Normal[idim];

     if (fabs(lPerp)>0.0) {
       for (idim=0,c=0.0;idim<3;idim++) c+=(x0[idim]-BoundingBoxFace[nface].x[idim])*BoundingBoxFace[nface].Normal[idim];
       t=-c/lPerp;

       if (t>0.0) {
         //check if the intersection point within the face
         for (idim=0,c0=0.0,c1=0.0;idim<3;idim++) {
           c0+=(x1PlaneNodeIndex[nface][idim]-x0PlaneNodeIndex[nface][idim])*(x0[idim]+t*l[idim]-BoundingBoxFace[nface].x[idim]);
           c1+=(x2PlaneNodeIndex[nface][idim]-x0PlaneNodeIndex[nface][idim])*(x0[idim]+t*l[idim]-BoundingBoxFace[nface].x[idim]);
         }

         if ((c0>-PIC::Mesh::mesh.EPS)&&(c0<BoundingBoxFace[nface].e0Length+PIC::Mesh::mesh.EPS) && (c1>-PIC::Mesh::mesh.EPS)&&(c1<BoundingBoxFace[nface].e1Length+PIC::Mesh::mesh.EPS)) {
           IntersectionTime.push_back(t);
         }

       }

     }
  }

  //determine intersections with the internals od the computational domain
  list<cInternalBoundaryConditionsDescriptor>::iterator InternalBoundaryDescriptor;
  double radiusSphere,*x0Sphere;
  double a,b,d,dx,dy,dz,sqrt_d;

  //spherical internal surface
#if DIM == 3
  cInternalSphericalData *Sphere;
#elif DIM == 2
  exit(__LINE__,__FILE__,"not yet");
#else
  cInternalSphere1DData *Sphere1D;
#endif



  for (InternalBoundaryDescriptor=PIC::Mesh::mesh.InternalBoundaryList.begin();InternalBoundaryDescriptor!=PIC::Mesh::mesh.InternalBoundaryList.end();InternalBoundaryDescriptor++) {
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

      dx=x0[0]-x0Sphere[0],dy=x0[1]-x0Sphere[1],dz=x0[2]-x0Sphere[2];
      a=l[0]*l[0]+l[1]*l[1]+l[2]*l[2];
      b=2.0*(l[0]*dx+l[1]*dy+l[2]*dz);
      c=dx*dx+dy*dy+dz*dz-radiusSphere*radiusSphere;

      if (c<0.0) {
        exit(__LINE__,__FILE__,"Error: the starting point of the column integration is within a internals spherical body");
      }

      d=b*b-4.0*a*c;

      if (d<0.0) continue;
      if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

      sqrt_d=sqrt(d);

      t=-(b+sqrt_d)/(2.0*a);
      if (t>0.0) IntersectionTime.push_back(t);

      t=-2.0*c/(b+sqrt_d);
      if (t>0.0) IntersectionTime.push_back(t);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
    }

  }

  //check if any intersections with the boundary of the domain have found
  if (IntersectionTime.size()==0) return false;

  //sort the intersection time
  sort(IntersectionTime.begin(),IntersectionTime.end());


  //determine weather the point 'x0' is within the computational domain
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* x0Node;

  x0Node=PIC::Mesh::mesh.findTreeNode(x0);

  if (x0Node==NULL) {
    //the point 'x0' is outside of the domain
    for (idim=0;idim<3;idim++) xStart[idim]=x0[idim]+IntersectionTime[0]*l[idim],xFinish[idim]=x0[idim]+IntersectionTime[1]*l[idim];

    //check if the 'start' and 'finish' nodes are within the domain
    while ((xStartNode=PIC::Mesh::mesh.findTreeNode(xStart))==NULL) {
      for (idim=0;idim<3;idim++) xStart[idim]+=(xFinish[idim]-xStart[idim])/100000;
    }

    while ((xFinishNode=PIC::Mesh::mesh.findTreeNode(xFinish))==NULL) {
      for (idim=0;idim<3;idim++) xFinish[idim]-=(xFinish[idim]-xStart[idim])/100000;
    }

    IntegrationPathLength=(IntersectionTime[1]-IntersectionTime[0])*sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
  }
  else {
    //the point 'x0' in within the domain
    xStartNode=x0Node;

    for (idim=0;idim<3;idim++) xStart[idim]=x0[idim],xFinish[idim]=x0[idim]+IntersectionTime[0]*l[idim];

    while ((xFinishNode=PIC::Mesh::mesh.findTreeNode(xFinish))==NULL) {
      for (idim=0;idim<3;idim++) xFinish[idim]-=(xFinish[idim]-xStart[idim])/100000;
    }

    IntegrationPathLength=IntersectionTime[0]*sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
  }


  return true;
}

//====================================================
void PIC::ColumnIntegration::GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xStartNode,double *l,double IntegrationPathLength,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)) {
  double lNormalized[3],c,x[3],dl=0.0,IntegratedPath=0.0,a0[ResultVectorLength],a1[ResultVectorLength];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xNode=xStartNode;
  int idim,i;
//  MPI_Status status;
  int IntegrationFinished=false;

  //the ratio between the step of the integration procedure and the local cell size
  static const double IntegrationStep2CellSizeRatio=0.3;

  for (i=0;i<ResultVectorLength;i++) ResultVector[i]=0.0,a0[i]=0.0,a1[i]=0.0;

  //normalize the pointing vector
  for (c=0.0,idim=0;idim<3;idim++) x[idim]=xStart[idim],c+=pow(l[idim],2);
  for (c=sqrt(c),idim=0;idim<3;idim++) lNormalized[idim]=l[idim]/c;

  double ExchangeBuffer[4];
  memcpy(ExchangeBuffer,x,3*sizeof(double));
  memcpy(ExchangeBuffer+3,&IntegrationPathLength,sizeof(double));

  MPI_Bcast(ExchangeBuffer,4,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  memcpy(x,ExchangeBuffer,3*sizeof(double));
  memcpy(&IntegrationPathLength,ExchangeBuffer+3,sizeof(double));


//  MPI_Bcast(&IntegrationPathLength,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
//  MPI_Bcast(x,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  xStartNode=PIC::Mesh::mesh.findTreeNode(x,xStartNode);

  //get the first value of the integrand function
  if (xNode->Thread==PIC::ThisThread) {
    Integrand(a0,ResultVectorLength,x,xNode);
//    if (PIC::ThisThread!=0) MPI_Send(a0,ResultVectorLength,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);
  }
//  else if (PIC::ThisThread==0) MPI_Recv(a0,ResultVectorLength,MPI_DOUBLE,xNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);


  while ((IntegratedPath<IntegrationPathLength)&&(IntegrationFinished==false)) {
    dl=IntegrationStep2CellSizeRatio*xNode->GetCharacteristicCellSize();

    if (dl>IntegrationPathLength-IntegratedPath) {
      dl=IntegrationPathLength-IntegratedPath;
      IntegrationFinished=true;
    }

//    if (PIC::ThisThread==0) for (idim=0;idim<3;idim++) x[idim]+=dl*lNormalized[idim];
//    MPI_Bcast(x,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
//    MPI_Bcast(&IntegrationFinished,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    for (idim=0;idim<3;idim++) x[idim]+=dl*lNormalized[idim];
    xNode=PIC::Mesh::mesh.findTreeNode(x,xNode);

    if (xNode==NULL) break;

    if (xNode->Thread==PIC::ThisThread) {
      Integrand(a1,ResultVectorLength,x,xNode);
//      if (PIC::ThisThread!=0) MPI_Send(a1,ResultVectorLength,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);
    }
//    else if (PIC::ThisThread==0) MPI_Recv(a1,ResultVectorLength,MPI_DOUBLE,xNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
    else for (i=0;i<ResultVectorLength;i++) a1[i]=0.0;

    for (i=0;i<ResultVectorLength;i++) {
      ResultVector[i]+=0.5*(a0[i]+a1[i])*dl;
      a0[i]=a1[i];
    }


    IntegratedPath+=dl;
//    MPI_Bcast(&IntegratedPath,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  }


  //collect the data
  MPI_Allreduce(ResultVector,a0,ResultVectorLength,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  memcpy(ResultVector,a0,ResultVectorLength*sizeof(double));

//  MPI_Bcast(ResultVector,ResultVectorLength,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
}


void PIC::ColumnIntegration::GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)) {
  double IntegrationPathLength,xStart[3],xFinish[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *xStartNode,*xFinishNode;

  if (FindIntegrationLimits(x0,l,IntegrationPathLength,xStart,xStartNode,xFinish,xFinishNode)==false) {
    for (int i=0;i<ResultVectorLength;i++) ResultVector[i]=0.0;

    return;
  }
  if ((xStartNode==NULL)||(xFinishNode==NULL)) exit(__LINE__,__FILE__,"Error: mehs node is not determined");

  return GetCoulumnIntegral(ResultVector,ResultVectorLength,xStart,xStartNode,l,IntegrationPathLength,Integrand);
}


























