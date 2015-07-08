/*
 * Dust_sampling.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: vtenishe
 */

//$Id$
//sampling procedures for the dust model

#include "pic.h"
#include "pic__model__electrically_charged_dust.h"

vector<ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation> ElectricallyChargedDust::Sampling::FluxMap::SampleLocations;
int ElectricallyChargedDust::Sampling::FluxMap::nZenithSurfaceElements,ElectricallyChargedDust::Sampling::FluxMap::nAzimuthalSurfaceElements;



//init the flux sampling procedure
void ElectricallyChargedDust::Sampling::FluxMap::Init(int nZenithElements,int nAzimuthalElements) {
  nZenithSurfaceElements=nZenithElements,nAzimuthalSurfaceElements=nAzimuthalElements;
}

//set the new sampling location
void ElectricallyChargedDust::Sampling::FluxMap::SetSamplingLocation(double *xLocation,double *xPrimary,double *xSecondary) {
  cSampleLocation NewLocation;

  NewLocation.SetGeneralSurfaceMeshParameters(nZenithSurfaceElements,nAzimuthalSurfaceElements);
  NewLocation.SetLocation(xLocation,xPrimary,xSecondary);
  SampleLocations.push_back(NewLocation);

  SampleLocations[SampleLocations.size()-1].Allocate();
}

void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::SetLocation(double *xLocation,double *xPrimary,double *xSecondary) {
  //reconstract the frame of reference associated with the point of the observations
  //e0 is along xPrimary-xLocation
  //e1 is the in the direction to xSecondary-xLocation ortogonal to e0
  //e2 is the cross product of e0 and e1
  double l0,l1,l2;
  int idim;

  //calculate e0
  for (l0=0.0,idim=0;idim<3;idim++) {
    e0[idim]=xPrimary[idim]-xLocation[idim];
    l0+=pow(e0[idim],2);
  }

  //calculate e1
  for (l0=sqrt(l0),l1=0.0,idim=0;idim<3;idim++) {
    e0[idim]/=l0;
    e1[idim]=xSecondary[idim]-xLocation[idim];
    l1+=e0[idim]*e1[idim];
  }

  for (l2=0.0,idim=0;idim<3;idim++) {
    e1[idim]-=l1*e0[idim];
    l2+=pow(e1[idim],2);
  }

  l2=sqrt(l2);

  if (l2>1.0E-12) {
    for (idim=0;idim<3;idim++) e1[idim]/=l2;
  }
  else {
    //xPrimary-xLocation and xSecondary-xLocation are alinged -> chose a vector that is ortogonal to e0
    if (fabs(e0[0])>1.0E-5) e1[0]=-e0[1],e1[1]=e0[0],e1[2]=0.0;
    else e1[0]=0.0,e1[1]=e0[2],e1[2]=-e0[1];

    l2=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
    for (idim=0;idim<3;idim++) e1[idim]/=l2;
  }

  //calculate e2
  e2[0]=+(e0[1]*e1[2]-e1[1]*e0[2]);
  e2[1]=-(e0[0]*e1[2]-e1[0]*e0[2]);
  e2[2]=+(e0[0]*e1[1]-e1[0]*e0[1]);

  //save the location and determine the node
  for (idim=0;idim<3;idim++) x[idim]=xLocation[idim];

  node=PIC::Mesh::mesh.findTreeNode(x);
  PIC::Mesh::mesh.fingCellIndex(x,iCell,jCell,kCell,node);

  //determine lon and lat of the secondary direction
  int nZenithElement,nAzimuthalElement;
  double t[3];

  for (idim=0;idim<3;idim++) t[idim]=xSecondary[idim]-xLocation[idim];
  GetSpeed(t,nZenithElement,Lat_xSecondary,nAzimuthalElement,Lon_xSecondary);

  Lat_xSecondary=(Pi/2.0-Lat_xSecondary)*180.0/Pi;
  Lon_xSecondary*=180.0/Pi;
}


void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::Allocate() {
  //allocate the samplingdata buffer
  int nSizeGroup;

  SampleData_NumberDensityFlux=new double* [nDustSizeSamplingIntervals];

  for (nSizeGroup=0;nSizeGroup<nDustSizeSamplingIntervals;nSizeGroup++) {
    SampleData_NumberDensityFlux[nSizeGroup]=new double [nZenithSurfaceElements*nAzimuthalSurfaceElements];

    //initialize the data buffer
    for (int i=0;i<nZenithSurfaceElements*nAzimuthalSurfaceElements;i++) SampleData_NumberDensityFlux[nSizeGroup][i]=0.0;
  }
}


//sampling procedure
void ElectricallyChargedDust::Sampling::FluxMap::Sampling() {
  int nSamplePoint;

  for (nSamplePoint=0;nSamplePoint<SampleLocations.size();nSamplePoint++) SampleLocations[nSamplePoint].Sampling();
}

void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::Sampling() {
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_ //execute the body of the sampling procedure only when the dust is ON
  long int ptr;
  PIC::Mesh::cDataBlockAMR *block=NULL;
  double v[3],ParticleWeight;
  PIC::ParticleBuffer::byte *ParticleData;
  int spec;

  //data is sampled only of the currect processor
  if (node->Thread!=PIC::ThisThread) return;

  //sample the particle data
  block=node->block;
  if (block==NULL) return;

  //scroll through the particle list
  ptr=block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

  if (ptr!=-1) {
    ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    spec=PIC::ParticleBuffer::GetI(ParticleData);

    if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      ParticleWeight=block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
      PIC::ParticleBuffer::GetV(v,ParticleData);

      //determine the element for sampling of the particle
      int nZenithElement,nAzimuthalElement,el;
      double Speed,ZenithAngle,AzimuthalAngle;

      Speed=GetSpeed(v,nZenithElement,ZenithAngle,nAzimuthalElement,AzimuthalAngle);
      el=GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
      SampleData_NumberDensityFlux[spec-_DUST_SPEC_][el]+=Speed*ParticleWeight;
    }

    ptr=PIC::ParticleBuffer::GetNext(ParticleData);
  }
#endif //_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
}

//============================================================================
//output sampled data
void ElectricallyChargedDust::Sampling::FluxMap::PrintSurfaceData(int nDataSet, bool PrintStateVectorFlag) {
  char fname[_MAX_STRING_LENGTH_PIC_];
  int nSamplePoint;

  for (nSamplePoint=0;nSamplePoint<SampleLocations.size();nSamplePoint++) {
    sprintf(fname,"%s/amps.dust.flux-map.group=%i.nSampleLocation=%i.out=%i.dat",PIC::OutputDataFileDirectory,nDataSet-_DUST_SPEC_,nSamplePoint,PIC::DataOutputFileNumber);
    SampleLocations[nSamplePoint].PrintSurfaceData(fname,nDataSet,PrintStateVectorFlag);
  }
}

void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag) {
  long int iZenith,iAzimuthal;
  FILE *fout=NULL,*fout2d=NULL;
  double x[3];

  //only sampled data for the dust species should be outputed
  int nSizeGroup=nDataSet-_DUST_SPEC_;
  if ((nSizeGroup<0)||(nSizeGroup>=nDustSizeSamplingIntervals)) return;

  //collect all sampling data
  double TempBuffer[nZenithSurfaceElements*nAzimuthalSurfaceElements];

  MPI_Reduce(SampleData_NumberDensityFlux[nSizeGroup],TempBuffer,nZenithSurfaceElements*nAzimuthalSurfaceElements,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(SampleData_NumberDensityFlux,TempBuffer,nZenithSurfaceElements*nAzimuthalSurfaceElements*sizeof(double));

  //output the data file
  if (PIC::ThisThread==0) {
    char fname2d[300];

    sprintf(fname2d,"%s.2d.dat",fname);

    fout=fopen(fname,"w");
    fout2d=fopen(fname2d,"w");

    //print the output file title
    fprintf(fout,"TITLE=\"Primary direction:Lat=0,Log=0; Secondary direction: Lat=%e,Lon=%e\n",Lat_xSecondary,Lon_xSecondary);
    fprintf(fout2d,"TITLE=\"Primary direction:Lat=0,Log=0; Secondary direction: Lat=%e,Lon=%e\n",Lat_xSecondary,Lon_xSecondary);

    //print the variable list
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\", \"Flux [s^{-1}m^{-2}]\"");
    fprintf(fout2d,"VARIABLES=\"Lon\", \"Lat\"");

    //print the number of variables and blocks
    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",(nZenithSurfaceElements+1)*nAzimuthalSurfaceElements,nZenithSurfaceElements*nAzimuthalSurfaceElements);
    fprintf(fout2d,"ZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthalSurfaceElements,nZenithSurfaceElements+1);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      GetSurfaceCoordinate(x,iZenith,iAzimuthal);

      fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

      double lon,lat;
      GetSurfaceLonLatNormal(lon,lat,iZenith,iAzimuthal);
      fprintf(fout2d,"%e %e ",lon,lat);

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iZenith==0) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else if (iZenith==nZenithSurfaceElements) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else {
          int iA,iZ,A[2],Z[2];

          Z[0]=iZenith-1,Z[1]=iZenith;

          A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1;
          A[1]=iAzimuthal;

          for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
        }

        //calculate the flux interpolated into the node of the spherical mesh and output in into the data file
        int el,nInterpolationElement;
        double InterpolationCoefficient,summInterpolationCoefficient=0.0;
        double Flux=0.0;

        for (nInterpolationElement=0;nInterpolationElement<InterpolationListLength;nInterpolationElement++) {
          el=InterpolationList[nInterpolationElement];
          InterpolationCoefficient=GetSurfaceElementArea(el);
          summInterpolationCoefficient+=InterpolationCoefficient;

          //there is no need to multiply by the surface element surface because of the conversion from the flux
          //through the surface elemet into the flux per m^{-2}
          Flux+=SampleData_NumberDensityFlux[nSizeGroup][el];
        }

        //normalize the interpolated values
        if (PIC::LastSampleLength!=0) {
          Flux/=summInterpolationCoefficient*PIC::LastSampleLength;
        }

        //output the interpolated values into a file
        fprintf(fout," %e",Flux);
        fprintf(fout2d," %e",Flux);
      }

      fprintf(fout,"\n");
      fprintf(fout2d,"\n");
    }


    //print the connectivity list
    int iAzimuthalMax,iAzimuthalMin;
    int nd0,nd1,nd2,nd3;

    for (iZenith=0;iZenith<nZenithSurfaceElements;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
      iAzimuthalMin=iAzimuthal;

      nd0=1+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements;
      nd1=1+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements;
      nd2=1+iAzimuthalMax+(iZenith+1)*nAzimuthalSurfaceElements;
      nd3=1+iAzimuthalMin+(iZenith+1)*nAzimuthalSurfaceElements;

      fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
    }

    fclose(fout);
    fclose(fout2d);
  }

  //re-init the content of teh sampling buffer
  for (int i=0;i<nZenithSurfaceElements*nAzimuthalSurfaceElements;i++) SampleData_NumberDensityFlux[nSizeGroup][i]=0.0;
}


double ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::GetSpeed(double *v,int &nZenithElement,double &ZenithAngle,int &nAzimuthalElement,double &AzimuthalAngle) {
  double r;
  int idim;

  //Normalize the velocity vector and determine its coordinates in the frame related for the observation point
  double Speed=0.0,l[3]={0.0,0.0,0.0};

  for (idim=0;idim<3;idim++) Speed+=pow(v[idim],2);
  Speed=sqrt(Speed);

  for (idim=0;idim<3;idim++) {
    l[0]+=e0[idim]*v[idim]/Speed;
    l[1]+=e1[idim]*v[idim]/Speed;
    l[2]+=e2[idim]*v[idim]/Speed;
  }

  if ((r=pow(l[0],2)+pow(l[1],2))>0.0) {
    AzimuthalAngle=acos(l[0]/sqrt(r));
    if (l[1]<0.0) AzimuthalAngle=2.0*Pi-AzimuthalAngle;
  }
  else AzimuthalAngle=0.0;

  ZenithAngle=acos(l[2]);

  #if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
  nZenithElement=(long int)((1.0-l[2])/dCosZenithAngle);
  #elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
  nZenithElement=(long int)(ZenithAngle/dZenithAngle);
  #else
  exit(__LINE__,__FILE__,"Error: wrong option");
  #endif

  nAzimuthalElement=(long int)(AzimuthalAngle/dAzimuthalAngle);

  if (nZenithElement==nZenithSurfaceElements) --nZenithElement;
  if (nAzimuthalElement==nAzimuthalSurfaceElements) --nAzimuthalElement;

  return Speed;
}



