//$Id$

/*
 * PostProcess3D_ColumnIntegration.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include "PostProcess3D.h"

int cPostProcess3D::cColumnIntegral::x0PlaneNodeIndex[6][3]={ {0,0,0},{1,0,0},       {0,0,0},{0,1,0},           {0,0,0},{0,0,1}};
int cPostProcess3D::cColumnIntegral::x1PlaneNodeIndex[6][3]={ {0,1,0},{1,1,0},       {1,0,0},{1,1,0},           {1,0,0},{1,0,1}};
int cPostProcess3D::cColumnIntegral::x2PlaneNodeIndex[6][3]={ {0,0,1},{1,0,1},       {0,0,1},{0,1,1},           {0,1,0},{0,1,1}};
int cPostProcess3D::cColumnIntegral::PlaneNormal[6][3]=     { {1,0,0},{1,0,0},       {0,1,0},{0,1,0},           {0,0,1},{0,0,1}};


void cPostProcess3D::cColumnIntegral::Init(double *xGlobalMin,double *xGlobalMax) {
  double *xmin,*xmax,e0Length,e1Length;
  int nface,i;

  xmin=xGlobalMin;
  xmax=xGlobalMax;

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
bool cPostProcess3D::cColumnIntegral::FindIntegrationLimits(double *x0,double *l,
    double& IntegrationPathLength,double *xStart,double *xFinish,
    double *xGlobalMin,double *xGlobalMax) {
  double t,lPerp,c,c0,c1;
  int nface,idim;
  std::vector<double> IntersectionTime;

  //determine the intersection time for all boundary faces
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

         if ((c0>0.0)&&(c0<BoundingBoxFace[nface].e0Length+0.0) && (c1>-0.0)&&(c1<BoundingBoxFace[nface].e1Length+0.0)) {
           IntersectionTime.push_back(t);
         }

       }

     }
  }

  //determine the intersection time with the triangulated surface (if any)
  int iStartFace,iFinishFace;
  double CutCellIntersectionTime=-1.0;

  if (CutCell::nBoundaryTriangleFaces!=0) {
    iStartFace=0;
    iFinishFace=CutCell::nBoundaryTriangleFaces;

    for (int i=iStartFace;i<iFinishFace;i++) {
      if (CutCell::BoundaryTriangleFaces[i].RayIntersection(x0,l,t,0.0)==true) {
        if ((CutCellIntersectionTime<0.0) || (t<CutCellIntersectionTime)) CutCellIntersectionTime=t;
      }
    }

    if (CutCellIntersectionTime>0.0) IntersectionTime.push_back(CutCellIntersectionTime);
  }


  //check if any intersections with the boundary of the domain have found
  bool InsideDomainFlag=true;

  for (idim=0;idim<3;idim++) if ((x0[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<x0[idim])) InsideDomainFlag=false;

  //sort the intersection time
  if (IntersectionTime.size()==0) return false;
  std::sort(IntersectionTime.begin(),IntersectionTime.end());


  if (InsideDomainFlag==false) {
    //the point 'x0' is outside of the domain
    for (idim=0;idim<3;idim++) xStart[idim]=x0[idim]+IntersectionTime[0]*l[idim],xFinish[idim]=x0[idim]+IntersectionTime[1]*l[idim];

    //check if the 'start' and 'finish' nodes are within the domain
    for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xStart[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xStart[idim])) InsideDomainFlag=false;

    while (InsideDomainFlag==false) {
      for (idim=0;idim<3;idim++) xStart[idim]+=(xFinish[idim]-xStart[idim])/100000;
      for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xStart[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xStart[idim])) InsideDomainFlag=false;
    }

    for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xFinish[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xFinish[idim])) InsideDomainFlag=false;

    while (InsideDomainFlag==false) {
      for (idim=0;idim<3;idim++) xFinish[idim]-=(xFinish[idim]-xStart[idim])/100000;
      for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xFinish[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xFinish[idim])) InsideDomainFlag=false;
    }

    IntegrationPathLength=(IntersectionTime[1]-IntersectionTime[0])*sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
  }
  else {
    //the point 'x0' in within the domain
    for (idim=0;idim<3;idim++) xStart[idim]=x0[idim],xFinish[idim]=x0[idim]+IntersectionTime[0]*l[idim];

    for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xFinish[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xFinish[idim])) InsideDomainFlag=false;

    while (InsideDomainFlag==false) {
      for (idim=0;idim<3;idim++) xFinish[idim]-=(xFinish[idim]-xStart[idim])/100000;
      for (InsideDomainFlag=true,idim=0;idim<3;idim++) if ((xFinish[idim]<xGlobalMin[idim]) || (xGlobalMax[idim]<xFinish[idim])) InsideDomainFlag=false;
    }

    IntegrationPathLength=IntersectionTime[0]*sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
  }

  return true;
}

//====================================================
void cPostProcess3D::cColumnIntegral::GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,double *l,double IntegrationPathLength,void (*Integrand)(double*,double*)) {
  double lNormalized[3],c,x[3],dl=0.0,IntegratedPath=0.0,a0[ResultVectorLength],a1[ResultVectorLength];
  int idim,i;
  int IntegrationFinished=false;

  //the ratio between the step of the integration procedure and the local cell size
  static const double IntegrationStep2CellSizeRatio=0.3;

  for (i=0;i<ResultVectorLength;i++) ResultVector[i]=0.0,a0[i]=0.0,a1[i]=0.0;

  //normalize the pointing vector
  for (c=0.0,idim=0;idim<3;idim++) x[idim]=xStart[idim],c+=pow(l[idim],2);
  for (c=sqrt(c),idim=0;idim<3;idim++) lNormalized[idim]=l[idim]/c;

  //get the first value of the integrand function
  for (i=0;i<3;i++) {
    if (x[i]<=PostProcess3D->xmin[i]) x[i]=PostProcess3D->xmin[i]+1.0E-8*fabs(PostProcess3D->xmin[i]);
    if (x[i]>=PostProcess3D->xmax[i]) x[i]=PostProcess3D->xmax[i]-1.0E-8*fabs(PostProcess3D->xmax[i]);
  }

  Integrand(a0,x);

  //calculate the integral
  while ((IntegratedPath<IntegrationPathLength)&&(IntegrationFinished==false)) {
    dl=IntegrationStep2CellSizeRatio*PostProcess3D->CharacteristicCellSize(x);

    if (dl>IntegrationPathLength-IntegratedPath) {
      dl=IntegrationPathLength-IntegratedPath;
      IntegrationFinished=true;
    }

    for (idim=0;idim<3;idim++) x[idim]+=dl*lNormalized[idim];

 //   if (PostProcess3D->IsInsideDomain(x)==false) break;
    for (i=0;i<3;i++) {
      if (x[i]<=PostProcess3D->xmin[i]) x[i]=PostProcess3D->xmin[i]+1.0E-8*fabs(PostProcess3D->xmin[i]);
      if (x[i]>=PostProcess3D->xmax[i]) x[i]=PostProcess3D->xmax[i]-1.0E-8*fabs(PostProcess3D->xmax[i]);
    }

    Integrand(a1,x);

    for (i=0;i<ResultVectorLength;i++) {
      ResultVector[i]+=0.5*(a0[i]+a1[i])*dl;
      a0[i]=a1[i];
    }

    IntegratedPath+=dl;
  }

}


void cPostProcess3D::cColumnIntegral::GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,double*)) {
  double IntegrationPathLength,xStart[3],xFinish[3];

  if (FindIntegrationLimits(x0,l,IntegrationPathLength,xStart,xFinish,PostProcess3D->xmin,PostProcess3D->xmax)==false) {
    for (int i=0;i<ResultVectorLength;i++) ResultVector[i]=0.0;
    return;
  }

  return GetCoulumnIntegral(ResultVector,ResultVectorLength,xStart,l,IntegrationPathLength,Integrand);
}

//===========================================================
//generate a circular column integral map
void cPostProcess3D::cColumnIntegral::cMap::Circular(double *xObservation,double *PrimaryAxisDirection,double *SecondaryAxisDirection,double HalfMapAngle,int nAzimuthPoints,int nZenithPoints,const char *fname,cColumnIntegrationSet* IntegrationSet) {
  int iZenithPoint,iAzimuthPoint,idim,i;
  double e0[3],e1[3],e2[3],l[3],theta,phi;
  double ZenithAngle,AzimuthAngle;
  double cosPhi,sinPhi,sinTheta,cosTheta;

  HalfMapAngle*=Pi/180.0;

  const double dZenitAngle=2.0*HalfMapAngle/(nZenithPoints-1);
  const double dAzimuthAngle=2.0*HalfMapAngle/(nAzimuthPoints-1);

  //determine the frame of reference for the map capculation
  for (i=0;i<3;i++) {
    e0[i]=PrimaryAxisDirection[i]-xObservation[i];
    e1[i]=SecondaryAxisDirection[i]-xObservation[i];
  }

  Vector3D::Normalize(e0);
  Vector3D::Orthogonalize(e0,e1);
  Vector3D::Normalize(e1);
  Vector3D::CrossProduct(e2,e0,e1);

  //create the output file and print the variable list
  FILE *fout;

  fout=fopen(fname,"w");
  fprintf(fout,"VARIABLES=\"Lon\", \"Lat\", \"Cut-cell Projection\", ");
  IntegrationSet->PrintVariableList(fout);

  fprintf(fout,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthPoints,nZenithPoints);

  //create the column integraion map
  int StateVectorLength=IntegrationSet->IntegrantVectorLength();
  double StateVector[StateVectorLength];


  for (iZenithPoint=0;iZenithPoint<nZenithPoints;iZenithPoint++) {
    ZenithAngle=HalfMapAngle-dZenitAngle*iZenithPoint;

    for (iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints;iAzimuthPoint++) {
      AzimuthAngle=HalfMapAngle-dAzimuthAngle*iAzimuthPoint;

      cosPhi=cos(AzimuthAngle);
      sinPhi=sin(AzimuthAngle);
      sinTheta=sin(ZenithAngle);
      cosTheta=cos(ZenithAngle);

      for (idim=0;idim<3;idim++) l[idim]=e0[idim]+sinPhi*e1[idim]+sinTheta*e2[idim];

      //calculate and output of the column integral
      ColumnIntegral->GetCoulumnIntegral(StateVector,StateVectorLength,xObservation,l,IntegrationSet->IntegrantVector);

      //determine the cut-cell projection
      double t,ProjectionCode=-1.0;
      int iStartFace,iFinishFace;

      iStartFace=0;
      iFinishFace=CutCell::nBoundaryTriangleFaces;

      for (i=iStartFace;i<iFinishFace;i++) {
        if (CutCell::BoundaryTriangleFaces[i].RayIntersection(xObservation,l,t,0.0)==true) {
          ProjectionCode=1.0;
        }
      }

      //output value of the column integral
      fprintf(fout,"%e %e %e ",ZenithAngle*180.0/Pi,AzimuthAngle*180.0/Pi,ProjectionCode);
      for (i=0;i<StateVectorLength;i++) fprintf(fout," %e ",StateVector[i]);
      fprintf(fout,"\n");
    }
  }

  //close the output file
  fclose(fout);
}

