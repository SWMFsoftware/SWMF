//$Id$
//the set of function for processing of the cut-cells and cut-faces of the internal irregular surfaces

#include "pic.h"

void PIC::Mesh::IrregularSurface::InitExternalNormalVector() {
  //calculate external normals to the faces
  double xStart[3],xFinish[3],l,l0,*norm;
  long int nface;
  cTriangleFace *fcptr;
  int idim,iIntersections;

  const double angleCosMin=cos(85.0/180.0*3.141592654);

  int nStartFace,nFinishFace,nTotalThreads,ThisThread,nFaceThread;

  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

  nFaceThread=nBoundaryTriangleFaces/nTotalThreads;
  nStartFace=nFaceThread*ThisThread;
  nFinishFace=nStartFace+nFaceThread;
  if (ThisThread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

  for (nface=nStartFace;nface<nFinishFace;nface++) {
    fcptr=BoundaryTriangleFaces+nface;

    //get the center point and the normal of the face
    fcptr->GetCenterPosition(xStart);
    norm=fcptr->ExternalNormal;

    do {
      for (idim=0,l=0.0,l0=0.0;idim<3;idim++) {
        xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
        l+=pow(xFinish[idim],2);
        l0+=xFinish[idim]*norm[idim];
      }

      l=sqrt(l);
    }
    while (fabs(l0)/l<angleCosMin);

    if (l0<0.0) l*=-1.0;
    for (idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+xStart[idim];

    //count face intersections
    iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(xStart,xFinish,fcptr);

    if (iIntersections%2!=0) {
      //the norm has to be reversed
      for (idim=0;idim<3;idim++) fcptr->ExternalNormal[idim]*=-1.0;
    }
  }

  //collect the surface normals
  double sendBuffer[3*2*nFaceThread];
  int thread,cnt;

  for (thread=0;thread<nTotalThreads;thread++) {
    nStartFace=nFaceThread*thread;
    nFinishFace=nStartFace+nFaceThread;
    if (thread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

    if (thread==ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(sendBuffer+3*cnt,BoundaryTriangleFaces[nface].ExternalNormal,3*sizeof(double));
    }

    MPI_Bcast(sendBuffer,3*(nFinishFace-nStartFace),MPI_DOUBLE,thread,MPI_GLOBAL_COMMUNICATOR);

    if (thread!=ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(BoundaryTriangleFaces[nface].ExternalNormal,sendBuffer+3*cnt,3*sizeof(double));
    }
  }


}

//check weather a point (x0) in insed the domain:
//if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
//l -> is a random ray (intersection search) direction
bool PIC::Mesh::IrregularSurface::CheckPointInsideDomain(double *x,PIC::Mesh::IrregularSurface::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS) {
  int nface,nfaceStart,nfaceFinish,iIntersections;
  double l,xFinish[3];
  int idim;
  bool flag=true;

  if (SurfaceTriangulation==NULL) return true;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

  //distribute ditrction of the search
  for (l=0.0,idim=0;idim<3;idim++) {
    xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
    l+=pow(xFinish[idim],2);
  }

  for (l=sqrt(l),idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+x[idim];

  //xFinish is outside of the domain -> the point outside of the surface
  //calculate the number of the face intersections between 'x' and 'xFinish'

  iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(x,xFinish);

  return (iIntersections%2==0) ? true : false;

 /*   SearchDirection[idim]/=l;


  static bool initflag=false;
  static int ThisThread,nTotalThreads;

  if (initflag==false) {
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
    MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);
    initflag=true;
  }

  if (ParallelCheck==true) {
    nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
    nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
    if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;
  }
  else nfaceStart=0,nfaceFinish=nSurfaceTriangulation;

  bool flagbuffer[nTotalThreads];

  do {
    //distribute ditrction of the search
    for (l=0.0,idim=0;idim<3;idim++) {
      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
    }

    if (ParallelCheck==true) MPI_Bcast(SearchDirection,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (l=sqrt(l),idim=0;idim<3;idim++) SearchDirection[idim]/=l;
    iIntersections=0;
    flag=true;

    //find intersections with the faces on the mesh
    for (nface=nfaceStart;nface<nfaceFinish;nface++) {
      if (SurfaceTriangulation[nface].RayIntersection(x,SearchDirection,EPS)==true) iIntersections++;

      for (l=0.0,idim=0;idim<3;idim++) l+=pow(SurfaceTriangulation[nface].ExternalNormal[idim]*SearchDirection[idim],2);
      if (l<1.0E-10) {
        flag=false;
        break;
      }

    }

    if (ParallelCheck==true) {
      MPI_Gather(&flag,sizeof(bool),MPI_CHAR,flagbuffer,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
      if (ThisThread==0) for (int thread=1;thread<nTotalThreads;thread++) if (flagbuffer[thread]==false) flag=false;
      MPI_Bcast(&flag,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
    }
  }
  while (flag==false);

  if (ParallelCheck==true) {
    int t;

    MPI_Allreduce(&iIntersections,&t,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    iIntersections=t;
  }


  return (2*(iIntersections/2)==iIntersections) ? true : false;
*/}
