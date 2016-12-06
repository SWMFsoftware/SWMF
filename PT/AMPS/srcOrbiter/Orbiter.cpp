/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

//the object name and the names of the frames
/*char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Orbiter";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";*/

double Orbiter::UpstreamBC::Velocity[3]={-30.0E3,0.0,0.0};
double Orbiter::UpstreamBC::NumberDensity=1.5E21;
double Orbiter::UpstreamBC::Temperature=293.0;

double Orbiter::Sampling::DragCorfficient=0.0;
double Orbiter::ProjectionOrbiterSurfaceArea=0.0;

char Orbiter::SurfaceMeshName[_MAX_STRING_LENGTH_PIC_]="Orbiter";
double Orbiter::DomainSizeMultiplierX=1.0,Orbiter::DomainSizeMultiplierY=1.0,Orbiter::DomainSizeMultiplierZ=1.0;


//=======================================================================================
//calculate the projection area of the spacecraft
double Orbiter::CalculateProjectionArea() {
  unsigned char *map=NULL;
  double *x,de0,de1,c0,c1,e0[3],e1[3],e2[3],e0Min=0.0,e0Max=0.0,e1Min=0.0,e1Max=0.0;
  int idim,iSurfaceElement,iNode;
  int i,j,iOffset,jOffset;

  int nPixels=2000;
  int ByteMapLength=1+(nPixels*nPixels)/8;

  //get the frame of reference with e0 and e1 normal to the ambient gas flow
  for (idim=0;idim<3;idim++) e2[idim]=UpstreamBC::Velocity[idim];

  Vector3D::Normalize(e2);
  Vector3D::GetNormFrame(e0,e1,e2);

  //get the limits of the projected coordinates
  for (iNode=0;iNode<CutCell::nBoundaryTriangleNodes;iNode++) {
    x=CutCell::BoundaryTriangleNodes[iNode].x;

    c0=Vector3D::DotProduct(e0,x);
    c1=Vector3D::DotProduct(e1,x);

    if (iNode==0) e0Min=c0,e0Max=c0;
    else {
      if (c0>e0Max) e0Max=c0;
      if (c0<e0Min) e0Min=c0;
    }

    if (iNode==0) e1Min=c1,e1Max=c1;
    else {
      if (c1>e1Max) e1Max=c1;
      if (c1<e1Min) e1Min=c1;
    }
  }


  de0=(e0Max-e0Min)/nPixels;
  de1=(e1Max-e1Min)/nPixels;

  //allocate the projection map
  map=new unsigned char[ByteMapLength];

  for (i=0;i<ByteMapLength;i++) map[i]=0;

  //loop though all surface elements
  int iSurfaceElementStart,iSurfaceElementStartEnd,nSurfaceElementThread;

  nSurfaceElementThread=CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads;
  iSurfaceElementStart=nSurfaceElementThread*PIC::ThisThread;
  iSurfaceElementStartEnd=iSurfaceElementStart+nSurfaceElementThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iSurfaceElementStartEnd=CutCell::nBoundaryTriangleFaces;


  for (iSurfaceElement=iSurfaceElementStart;iSurfaceElement<iSurfaceElementStartEnd;iSurfaceElement++) {
    int nTestPoints;
    int t,iLocal,jLocal;
    double xLocal[2],dxLocal;

    double e0TriangleMax,e0TriangleMin,e1TriangleMax,e1TriangleMin;

    c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face);
    c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face);

    e0TriangleMin=c0,e0TriangleMax=c0;
    e1TriangleMin=c1,e1TriangleMax=c1;

    for (i=1;i<3;i++) {
      switch (i) {
      case 1:
        c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x1Face);
        c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x1Face);
        break;
      case 2:
        c0=Vector3D::DotProduct(e0,CutCell::BoundaryTriangleFaces[iSurfaceElement].x2Face);
        c1=Vector3D::DotProduct(e1,CutCell::BoundaryTriangleFaces[iSurfaceElement].x2Face);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: something is wrong with counting");
      }

      if (c0<e0TriangleMin) e0TriangleMin=c0;
      if (c0>e0TriangleMax) e0TriangleMax=c0;

      if (c1<e1TriangleMin) e1TriangleMin=c1;
      if (c1>e1TriangleMax) e1TriangleMax=c1;
    }

    nTestPoints=max(e0TriangleMax-e0TriangleMin,e1TriangleMax-e1TriangleMin)/min(de0,de1)*4;
    dxLocal=1.0/nTestPoints;

    for (iLocal=0;iLocal<nTestPoints;iLocal++) {
      xLocal[0]=iLocal*dxLocal;

      if (xLocal[0]<1.0) {
        for (jLocal=0;jLocal<nTestPoints;jLocal++) {
          xLocal[1]=jLocal*dxLocal;

          if ((xLocal[1]<1.0)&&(xLocal[0]+xLocal[1]<1.0)) {
            //the point belong to the surface triangle -> register it
            for (idim=0;idim<3;idim++) x[idim]=CutCell::BoundaryTriangleFaces[iSurfaceElement].x0Face[idim]+
                xLocal[0]*CutCell::BoundaryTriangleFaces[iSurfaceElement].e0[idim]+xLocal[1]*CutCell::BoundaryTriangleFaces[iSurfaceElement].e1[idim];

            c0=Vector3D::DotProduct(e0,x)-e0Min;
            c1=Vector3D::DotProduct(e1,x)-e1Min;

            int iPixel,jPixel,iGlobalPixelIndex,iByte,iBit;

            iPixel=c0/de0;
            jPixel=c1/de1;

            iGlobalPixelIndex=iPixel+jPixel*nPixels;

            iByte=iGlobalPixelIndex/8;
            iBit=iGlobalPixelIndex%8;

            //mark the shadow bit
            if ((iByte<ByteMapLength)&&(iBit<8)) {
              map[iByte]|=(1<<iBit);
            }
          }
        }
      }
    }
  }

  //combine maps from all processors, and calculate the total projection area
  int nShadowPixel=0;

  if (PIC::ThisThread==0) {
    MPI_Status status;
    int thread;
    unsigned char *buffer=new unsigned char[ByteMapLength];

    for (thread=1;thread<PIC::nTotalThreads;thread++) {
      MPI_Recv(buffer,ByteMapLength,MPI_UNSIGNED_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

      for (i=0;i<ByteMapLength;i++) map[i]|=buffer[i];
    }

    delete [] buffer;

    for (i=0;i<ByteMapLength;i++) for (j=0;j<8;j++) if ((map[i]&(1<<j)) != 0) nShadowPixel++;
  }
  else {
    MPI_Send(map,ByteMapLength,MPI_UNSIGNED_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
  }

  MPI_Bcast(&nShadowPixel,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

  ProjectionOrbiterSurfaceArea=(e0Max-e0Min)*(e1Max-e1Min)*double(nShadowPixel)/double(nPixels*nPixels);
  if (PIC::ThisThread==0) printf("$PREFIX: Orbiter Projection Surface Area: %e\n",ProjectionOrbiterSurfaceArea);

  //prepare the projection map of the orbiter
  if (PIC::ThisThread==0) {
    FILE *fout;
    char fname[100];

    sprintf(fname,"%s/OrbiterProjectionArea.dat",PIC::OutputDataFileDirectory);
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLE=\"X\",\"Y\",\"Surface Projection Flag\"\n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nPixels,nPixels);

    for (j=0;j<nPixels;j++) for (i=0;i<nPixels;i++) {
      int flag,iGlobalPixelIndex,iByte,iBit;

      iGlobalPixelIndex=i+j*nPixels;

      iByte=iGlobalPixelIndex/8;
      iBit=iGlobalPixelIndex%8;

      if ((iByte<ByteMapLength)&&(iBit<8)) {
        flag=((map[iByte]&(1<<iBit))!=0) ? 1 : 0;
      }
      else {
        flag=0;
      }

      fprintf(fout,"%e %e %i\n",e0Min+(i+0.5)*de0,e1Min+(j+0.5)*de1,flag);
    }

    fclose(fout);
  }

  delete [] map;
}


