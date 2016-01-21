/*
 * VIRTIS-M.cpp
 *
 *  Created on: Jan 20, 2016
 *      Author: vtenishe
 */


#include "VIRTIS-M.h"

//parameters of the field of view
const double cVirtisM::cBlockNucleus::AngularFieldOfView=3.67/180.0*Pi;
const int cVirtisM::cBlockNucleus::nFieldOfViewPixels=256;
const double cVirtisM::cBlockNucleus::dAnglePixel=AngularFieldOfView/nFieldOfViewPixels;

//set the mask matrix
void cVirtisM::cBlockNucleus::cMaskPixel::SetMask(int nx,int ny) {
  int i,j;

  mask=new bool* [nx];
  mask[0]=new bool[nx*ny];

  for (i=1;i<nx;i++) mask[i]=mask[i-1]+ny;

  for (i=0;i<nx;i++) for (j=0;j<ny;j++) mask[i][j]=false;
}

//set the pixel limit
void cVirtisM::cBlockNucleus::SetPixelLimits(int dxPix,int dyPix) {
  nxVirtisBlockPixelRange=dxPix;
  nyVirtisBlockPixelRange=dyPix;
}


//init the mask
void cVirtisM::cBlockNucleus::SetBlock(SpiceDouble et,int nNucleusSurfaceFaces,CutCell::cTriangleFace *NucleusSurfaceFaces) {
  int i,j,idim,di,dj;
  double l[3];

  //init exis of the frame of the reference
  SpiceDouble StateRosetta[6],lt;

  spkezr_c("ROSETTA",et,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateRosetta,&lt);

  for (i=0;i<3;i++) {
    xRosetta[i]=1.0E3*StateRosetta[i];

    e0[i]=-xRosetta[i];
    e1[i]=-1.0E3*StateRosetta[i+3];
  }

  Vector3D::Normalize(e0);
  Vector3D::Orthogonalize(e0,e1);
  Vector3D::CrossProduct(e2,e0,e1);

  //init the blocking by the nucleus
  for (i=(nFieldOfViewPixels-nxVirtisBlockPixelRange)/2;i<(nFieldOfViewPixels+nxVirtisBlockPixelRange)/2;i++) {
    for (j=(nFieldOfViewPixels-nyVirtisBlockPixelRange)/2;j<(nFieldOfViewPixels+nyVirtisBlockPixelRange)/2;j++) {
      di=i-nFieldOfViewPixels/2;
      dj=j-nFieldOfViewPixels/2;

      //get the pointing direction
      for (idim=0;idim<3;idim++) l[idim]=e0[idim]+e1[idim]*sin(di*dAnglePixel)+e2[idim]*sin(dj*dAnglePixel);
      Vector3D::Normalize(l);

      //check if the direction is intersected with the surface of the nucleus
      //determine the intersection time with the triangulated surface (if any)
      int iStartFace,iFinishFace,iFace;
      double t;

      if (nNucleusSurfaceFaces!=0) {
        iStartFace=0;
        iFinishFace=nNucleusSurfaceFaces;

        for (iFace=iStartFace;iFace<iFinishFace;iFace++) {
          if (NucleusSurfaceFaces[iFace].RayIntersection(xRosetta,l,t,0.0)==true) {
            //there is intersection of the line of sight with the nucleus surface -> add it to the block mask
            NucleusMask.mask[i][j]=true;
          }
        }
      }

    }
  }

  //set mask with the additional blocking pixels
  for (i=0;i<nFieldOfViewPixels;i++) for (j=0;j<nFieldOfViewPixels;j++) VistisMask.mask[i][j]=NucleusMask.mask[i][j];

  for (i=0;i<nFieldOfViewPixels;i++) {
    for (j=0;j<nFieldOfViewPixels;j++) {
      if (NucleusMask.mask[i][j]==true) {
        //add additional blocking
        for (di=-AdditionalBlockedPixels;di<AdditionalBlockedPixels;di++) VistisMask.mask[i+di][j]=true;
      }
    }
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


  if (rank==0) {
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"xPixel\", \"jPixel\", \"Nucleus Projection\", \"Additional Virtis Block\" ");

    if (PostProcessor!=NULL) {
      fprintf(fout,", ");
      IntegrationSet->PrintVariableList(fout);
    }

    fprintf(fout,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nxVirtisBlockPixelRange+1,nyVirtisBlockPixelRange);
  }

  //calculate the integral map
  int jGlobalMin=(nFieldOfViewPixels-nyVirtisBlockPixelRange)/2,jGlobalMax=(nFieldOfViewPixels+nxVirtisBlockPixelRange)/2,jLocalMin,jLocalMax;
  int nHorizontalPoints=jGlobalMax-jGlobalMin+1;
  int StateVectorLength=IntegrationSet->IntegrantVectorLength();

  double DataBuffer[nHorizontalPoints*StateVectorLength];
  double GlobalDataBuffer[nHorizontalPoints*StateVectorLength];
  double StateVector[StateVectorLength];
  int nPointPerThread=nHorizontalPoints/size;

  jLocalMin=jGlobalMin+rank*nPointPerThread;
  jLocalMax=jGlobalMin+jLocalMin+nPointPerThread-1;

  if (rank==size-1) jLocalMax=jGlobalMax-1;


  for (i=(nFieldOfViewPixels-nxVirtisBlockPixelRange)/2;i<(nFieldOfViewPixels+nxVirtisBlockPixelRange)/2;i++) {
    for (j=0;j<nHorizontalPoints*StateVectorLength;j++) DataBuffer[j]=0.0,GlobalDataBuffer[j]=0.0;

    if (PostProcessor!=NULL) {
      for (j=jLocalMin;j<=jLocalMax;j++) {
        di=i-nxVirtisBlockPixelRange/2;
        dj=j-nyVirtisBlockPixelRange/2;

        //get the pointing direction
        for (idim=0;idim<3;idim++) l[idim]=e0[idim]+e1[idim]*sin(di*dAnglePixel)+e2[idim]*sin(dj*dAnglePixel);
        Vector3D::Normalize(l);

        //calculate the integral
        PostProcessor->ColumnIntegral.GetCoulumnIntegral(StateVector,StateVectorLength,xRosetta,l,IntegrationSet->IntegrantVector);

        //save the state vector in the data buffer
        memcpy(DataBuffer+(j-jGlobalMin)*StateVectorLength,StateVector,StateVectorLength*sizeof(double));
      }

      //collect output value of the column integral
       MPI_Reduce(DataBuffer,GlobalDataBuffer,nHorizontalPoints*StateVectorLength,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }

     //output
     if (rank==0) {
       for (j=jGlobalMin;j<=jGlobalMax;j++) {
         fprintf(fout,"%i %i %i %i ",i,j,((NucleusMask.mask[i][j]==false) ? 1 : 0),((VistisMask.mask[i][j]==false) ? 1 : 0));

         if (PostProcessor!=NULL) {
           IntegrationSet->PostProcessColumnIntegralVector(GlobalDataBuffer+(j-jGlobalMin)*StateVectorLength);
           for (int ii=0;ii<StateVectorLength;ii++) fprintf(fout," %e ",GlobalDataBuffer[ii+(j-jGlobalMin)*StateVectorLength]);
         }

         fprintf(fout,"\n");
       }
     }
    }

  //close the output file
  if (rank==0) fclose(fout);
}
