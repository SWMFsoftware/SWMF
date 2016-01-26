//$Id$

/*
 * VIRTIS-M.cpp
 *
 *  Created on: Jan 20, 2016
 *      Author: vtenishe
 */


#include "VIRTIS-M.h"

//parameters of the field of view
const double cVirtisM::AngularFieldOfView=3.67/180.0*Pi;
const int cVirtisM::nFieldOfViewPixels=256;
const double cVirtisM::dAnglePixel=AngularFieldOfView/nFieldOfViewPixels;

//set the pixel limit
void cVirtisM::cBlockNucleus::SetPixelLimits(int dxPix,int dyPix) {
  nxVirtisBlockPixelRange=dxPix;
  nyVirtisBlockPixelRange=dyPix;
}


void cVirtisM::SetFrameAxis(SpiceDouble et) {
  SpiceDouble StateRosetta[6],lt;

  //get the location of the spacecraft
  spkezr_c("ROSETTA",et,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateRosetta,&lt);
  for (int i=0;i<3;i++) xRosetta[i]=1.0E3*StateRosetta[i];


  //get the coordinate frame related to the instrument
  int MAXBND=4;
  int WDSIZE=32;

  SpiceChar    frame  [WDSIZE];
  SpiceChar    shape  [WDSIZE];

  SpiceDouble  bounds [MAXBND][3];
  SpiceDouble  bsight [3];
  SpiceInt     n;

  getfov_c ( -226210, MAXBND, WDSIZE, WDSIZE, shape, frame, bsight, &n, bounds );

  SpiceDouble fmatrix[3][3];
  pxform_c(frame,"67P/C-G_CK",et,fmatrix);

  double e0V[3]={0.0,1.0,0.0},e1V[3]={1.0,0.0,0.0},e2V[3]={0.0,0.0,1};

  mxv_c(fmatrix,e0V,e0);
  mxv_c(fmatrix,e1V,e1);
  mxv_c(fmatrix,e2V,e2);
}

void cVirtisM::GetNucleusPixelCoordinate(int &e0PixelNucleusCoordinate,int &e1PixelNucleusCoordinate,SpiceDouble et) {
  double c0=0.0,c1=0.0,c2=0.0;
  double e1PixelPosition,e0PixelPosition;
  int idim;

  SetFrameAxis(et);

  for (idim=0,c0=0.0,c1=0.0;idim<3;idim++) c0-=xRosetta[idim]*e0[idim],c1-=xRosetta[idim]*e1[idim],c2-=xRosetta[idim]*e2[idim];

  e0PixelNucleusCoordinate=(int)(atan(c0/c2)/dAnglePixel)+nFieldOfViewPixels/2;
  e1PixelNucleusCoordinate=(int)(atan(c1/c2)/dAnglePixel)+nFieldOfViewPixels/2;
}

//init the mask
void cVirtisM::cBlockNucleus::SetBlock(SpiceDouble et,int nNucleusSurfaceFaces,CutCell::cTriangleFace *NucleusSurfaceFaces) {
  int i,j,idim,di,dj;
  double l[3];

  //get the total number of the processors used in the calculations
  int size,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  //set the limits of the search window
  int imin,imax,jmin,jmax;
  Virtis->GetNucleusPixelCoordinate(iCenter,jCenter,et);

  //the range of the veriation of the pixel indexes
  imin=iCenter-nxVirtisBlockPixelRange/2;
  if (imin<0) imin=0;

  imax=iCenter+nxVirtisBlockPixelRange/2;
  if (imax>=nFieldOfViewPixels) imax=nFieldOfViewPixels-1;

  jmin=jCenter-nyVirtisBlockPixelRange/2;
  if (jmin<0) jmin=0;

  jmax=jCenter+nyVirtisBlockPixelRange/2;
  if (jmax>=nFieldOfViewPixels) jmax=nFieldOfViewPixels-1;

  //reset the nucleus block mask
  for (i=0;i<nFieldOfViewPixels;i++) for (j=0;j<nFieldOfViewPixels;j++) NucleusMask[i][j]=0;


  //init the blocking by the nucleus
//#pragma omp parallel for
  for (i=imin;i<=imax;i++) {
    for (j=jmin;j<=jmax;j++) {
      di=i-nFieldOfViewPixels/2;
      dj=j-nFieldOfViewPixels/2;

      //get the pointing direction
      for (idim=0;idim<3;idim++) l[idim]=Virtis->e2[idim]+Virtis->e0[idim]*tan(di*dAnglePixel)+Virtis->e1[idim]*tan(dj*dAnglePixel);
      Vector3D::Normalize(l);
      memcpy(InstrumentPointing[i][j],l,3*sizeof(double));

      //check if the direction is intersected with the surface of the nucleus
      //determine the intersection time with the triangulated surface (if any)
      int iStartFace,iFinishFace,iFace,nFacesPerThread;
      double t;

      if (nNucleusSurfaceFaces!=0) {
        nFacesPerThread=nNucleusSurfaceFaces/size;

        iStartFace=rank*nFacesPerThread;
        iFinishFace=iStartFace+nFacesPerThread-1;

        if (rank==size-1) iFinishFace=nNucleusSurfaceFaces-1;

        //check intersection with the subset of the faces assigned to this processor
        for (iFace=iStartFace;iFace<=iFinishFace;iFace++) {
          if (NucleusSurfaceFaces[iFace].RayIntersection(Virtis->xRosetta,l,t,0.0)==true) {
            //there is intersection of the line of sight with the nucleus surface -> add it to the block mask
            NucleusMask[i][j]=1;
            break;
          }
        }
      }

    }
  }

  //collect NucleusMask.mask from all processors
  int tmpBuffer[nFieldOfViewPixels*nFieldOfViewPixels];

  MPI_Allreduce(NucleusMask[0],tmpBuffer,nFieldOfViewPixels*nFieldOfViewPixels,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  memcpy(NucleusMask[0],tmpBuffer,nFieldOfViewPixels*nFieldOfViewPixels*sizeof(int));

  //set mask with the additional blocking pixels
  for (i=0;i<nFieldOfViewPixels;i++) for (j=0;j<nFieldOfViewPixels;j++) {
    NucleusMask[i][j]=(NucleusMask[i][j]!=0) ? true : false;
    VirtisMask[i][j]=NucleusMask[i][j];
  }

  for (j=0;j<nFieldOfViewPixels;j++) {
    for (i=0;i<nFieldOfViewPixels;i++) {
      if (NucleusMask[i][j]==true) {
        //add additional blocking
        for (di=-AdditionalBlockedPixels;di<AdditionalBlockedPixels;di++) if ((i+di>=0)&&(i+di<nFieldOfViewPixels)) VirtisMask[i+di][j]=true;
      }
    }
  }

  //calculate the limb
  if (rank==0) {
    FILE *fNucleus=fopen("nucleus.block.dat","w");
    FILE *fVirtis=fopen("virtis.block.dat","w");

    fprintf(fNucleus,"VARIABLES=\"xPixel\", \"jPixel\"\n");
    fprintf(fVirtis,"VARIABLES=\"xPixel\", \"jPixel\"\n");

    for (i=0;i<nFieldOfViewPixels;i++) for (j=0;j<nFieldOfViewPixels;j++) {
      bool NucleusMaskBoundary=false,VirtisMaskBoundary=false;
      int di,dj;

      for (di=-1;di<=1;di+=2) if ((i+di>=0)&&(i+di<nFieldOfViewPixels)) {
        for (dj=-1;dj<=1;dj+=2) if ((j+dj>=0)&&(j+dj<nFieldOfViewPixels)) {
          if (NucleusMask[i][j]!=NucleusMask[i+di][j+dj]) NucleusMaskBoundary=true;
          if (VirtisMask[i][j]!=VirtisMask[i+di][j+dj]) VirtisMaskBoundary=true;
        }
      }

      if (NucleusMaskBoundary==true) fprintf(fNucleus,"%i %i\n",i,j);
      if (VirtisMaskBoundary==true) fprintf(fVirtis,"%i %i\n",i,j);
    }

    fclose(fNucleus);
    fclose(fVirtis);
  }


}

//======================================================================
//get the column integral map
void cVirtisM::cBlockNucleus::GetColumnIntegralMap(const char *fname,cPostProcess3D::cColumnIntegral::cColumnIntegrationSet* IntegrationSet,cPostProcess3D* PostProcessor) {
  int i,j,di,dj,idim;
  double l[3];

  //create the output file and print the variable list
  FILE *fout;
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  //init exis of the frame of the reference
  int imin,imax,jmin,jmax;

  //the range of the veriation of the pixel indexes
  imin=iCenter-nxVirtisBlockPixelRange/2;
  if (imin<0) imin=0;

  imax=iCenter+nxVirtisBlockPixelRange/2;
  if (imax>=nFieldOfViewPixels) imax=nFieldOfViewPixels-1;

  jmin=jCenter-nyVirtisBlockPixelRange/2;
  if (jmin<0) jmin=0;

  jmax=jCenter+nyVirtisBlockPixelRange/2;
  if (jmax>=nFieldOfViewPixels) jmax=nFieldOfViewPixels-1;

  //calculate the integral map
  int jLocalMin,jLocalMax;
  int StateVectorLength=IntegrationSet->IntegrantVectorLength();

  double DataBuffer[nFieldOfViewPixels*StateVectorLength];
  double GlobalDataBuffer[nFieldOfViewPixels*StateVectorLength];
  double StateVector[StateVectorLength];

  int nPointPerThread=(jmax-jmin+1)/size;

  jLocalMin=jmin+rank*nPointPerThread;
  jLocalMax=jLocalMin+nPointPerThread-1;

  if (rank==size-1) jLocalMax=jmax-1;

  if (rank==0) {
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"xPixel\", \"jPixel\", \"Nucleus Projection\", \"Additional Virtis Block\" ");

    if (PostProcessor!=NULL) {
      fprintf(fout,", ");
      IntegrationSet->PrintVariableList(fout);
    }

    fprintf(fout,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",jmax-jmin+1,imax-imin+1);
  }

  for (i=imin;i<=imax;i++) {
    for (j=0;j<nFieldOfViewPixels;j++) DataBuffer[j]=0.0,GlobalDataBuffer[j]=0.0;

    if (PostProcessor!=NULL) {
      for (j=jLocalMin;j<=jLocalMax;j++) {
        //calculate the integral
        memcpy(l,InstrumentPointing[i][j],3*sizeof(double));
        PostProcessor->ColumnIntegral.GetCoulumnIntegral(StateVector,StateVectorLength,Virtis->xRosetta,l,IntegrationSet->IntegrantVector);

        //save the state vector in the data buffer
        memcpy(DataBuffer+j*StateVectorLength,StateVector,StateVectorLength*sizeof(double));
      }

      //collect output value of the column integral
       MPI_Reduce(DataBuffer,GlobalDataBuffer,nFieldOfViewPixels*StateVectorLength,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }

     //output
     if (rank==0) {
       for (j=jmin;j<=jmax;j++) {
         fprintf(fout,"%i %i %i %i ",i,j,((NucleusMask[i][j]==false) ? 1 : 0),((VirtisMask[i][j]==false) ? 1 : 0));

         if (PostProcessor!=NULL) {
           IntegrationSet->PostProcessColumnIntegralVector(GlobalDataBuffer+j*StateVectorLength);
           for (int k=0;k<StateVectorLength;k++) fprintf(fout," %e ",GlobalDataBuffer[k+j*StateVectorLength]);
         }

         fprintf(fout,"\n");
       }
     }
    }

  //close the output file
  if (rank==0) fclose(fout);
}
