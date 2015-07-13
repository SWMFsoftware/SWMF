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
    e0[idim]=-(xPrimary[idim]-xLocation[idim]); //the minus is needed to have the particles comin from point 'xPrimary' have coordinated Lon=0, Lat=0 on the 2D map
    l0+=pow(e0[idim],2);
  }

  //calculate e1
  for (l0=sqrt(l0),l1=0.0,idim=0;idim<3;idim++) {
    e0[idim]/=l0;
    e1[idim]=-(xSecondary[idim]-xLocation[idim]); //similar reason as for 'e0' but for the point 'xSecondary'
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
  long int nZenithElement,nAzimuthalElement;
  double t[3];

  for (idim=0;idim<3;idim++) t[idim]=xSecondary[idim]-xLocation[idim];
  GetSpeed(t,Lat_xSecondary,nZenithElement,Lon_xSecondary,nAzimuthalElement);

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
  double v[3],ParticleWeight,GrainRadius;
  PIC::ParticleBuffer::byte *ParticleData;
  int spec,SizeGroup;

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
      GrainRadius=GetGrainRadius(ParticleData);

      //determine the element for sampling of the particle
      long int nZenithElement,nAzimuthalElement,el;
      double Speed,ZenithAngle,AzimuthalAngle;

      SizeGroup=(int)(log(GrainRadius/ElectricallyChargedDust::minDustRadius)/ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals);
      Speed=GetSpeed(v,ZenithAngle,nZenithElement,AzimuthalAngle,nAzimuthalElement);
      el=GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
      SampleData_NumberDensityFlux[SizeGroup][el]+=Speed*ParticleWeight;
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
    sprintf(fname,"%s/amps.dust.flux-map.nSampleLocation=%i.out=%i.dat",PIC::OutputDataFileDirectory,nSamplePoint,nDataSet);
    SampleLocations[nSamplePoint].PrintSurfaceData(fname,PrintStateVectorFlag);
  }
}

void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::PrintSurfaceData(const char *fname, bool PrintStateVectorFlag) {
  long int iZenith,iAzimuthal;
  int SizeGroup;
  FILE *fout=NULL;
  double x[3];

  //collect all sampling data
  for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
    double TempBuffer[nZenithSurfaceElements*nAzimuthalSurfaceElements];

    MPI_Reduce(SampleData_NumberDensityFlux[SizeGroup],TempBuffer,nZenithSurfaceElements*nAzimuthalSurfaceElements,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
    memcpy(SampleData_NumberDensityFlux[SizeGroup],TempBuffer,nZenithSurfaceElements*nAzimuthalSurfaceElements*sizeof(double));
  }

  //output the data file
  if (PIC::ThisThread==0) {
    char fname2d[300];

    sprintf(fname2d,"%s.2d.dat",fname);
    fout=fopen(fname,"w");

    //print the output file title
    fprintf(fout,"TITLE=\"Primary direction:Lat=0,Log=0; Secondary direction: Lat=%e,Lon=%e\"\n",Lat_xSecondary,Lon_xSecondary);

    //print the variable list
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");

    for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
      double r0,r1;

      r0=ElectricallyChargedDust::minDustRadius*exp(SizeGroup*ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals);
      r1=ElectricallyChargedDust::minDustRadius*exp((1+SizeGroup)*ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals);

      fprintf(fout,", \"Flux [s^{-1}m^{-2}] (%e<a<%e)\"",r0,r1);
    }

    fprintf(fout,", \"Total Flux [s^{-1}m^{-2}]\"");

    //print the number of variables and blocks
    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",(nZenithSurfaceElements+1)*nAzimuthalSurfaceElements,nZenithSurfaceElements*nAzimuthalSurfaceElements);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) {
      for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
        GetSurfaceCoordinate(x,iZenith,iAzimuthal);

        fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

        if (PrintStateVectorFlag==true) {
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
          //sampled data summed for all dust sizes
          double TotalFlux=0.0;

          for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
            int el,nInterpolationElement;
            double InterpolationCoefficient,summInterpolationCoefficient=0.0;
            double Flux=0.0;

            for (nInterpolationElement=0;nInterpolationElement<InterpolationListLength;nInterpolationElement++) {
              el=InterpolationList[nInterpolationElement];
              InterpolationCoefficient=GetSurfaceElementArea(el);
              summInterpolationCoefficient+=InterpolationCoefficient;

              //there is no need to multiply by the surface element surface because of the conversion from the flux
              //through the surface elemet into the flux per m^{-2}
              Flux+=InterpolationCoefficient*SampleData_NumberDensityFlux[SizeGroup][el];
            }

            //normalize the interpolated values
            if (PIC::LastSampleLength!=0) {
              Flux/=summInterpolationCoefficient*PIC::LastSampleLength;
            }

            //output the interpolated values into a file
            TotalFlux+=Flux;

            fprintf(fout," %e",Flux);
          }

          fprintf(fout," %e",TotalFlux);
        }

        fprintf(fout,"\n");
      }
   }

    //print the connectivity list
    int iAzimuthalMax,iAzimuthalMin;
    int nd0,nd1,nd2,nd3;


    for (iZenith=0;iZenith<nZenithSurfaceElements;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      //connectivity list for the spherical mesh
      iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
      iAzimuthalMin=iAzimuthal;

      nd0=1+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements;
      nd1=1+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements;
      nd2=1+iAzimuthalMax+(iZenith+1)*nAzimuthalSurfaceElements;
      nd3=1+iAzimuthalMin+(iZenith+1)*nAzimuthalSurfaceElements;

      fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
    }

    fclose(fout);
  }

  //output lon-lat map
  Print2dMap(fname,PrintStateVectorFlag);

  //re-init the content of teh sampling buffer
  for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
    for (int i=0;i<nZenithSurfaceElements*nAzimuthalSurfaceElements;i++) SampleData_NumberDensityFlux[SizeGroup][i]=0.0;
  }
}

void ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::Print2dMap(const char *fname, bool PrintStateVectorFlag) {
  long int iZenith,iAzimuthal;
  int SizeGroup;
  FILE *fout2d=NULL;

  //output the data file
  if (PIC::ThisThread==0) {
    char fname2d[300];

    sprintf(fname2d,"%s.2d.dat",fname);
    fout2d=fopen(fname2d,"w");

    //print the output file title
    fprintf(fout2d,"TITLE=\"Primary direction:Lat=0,Log=0; Secondary direction: Lat=%e,Lon=%e\"\n",Lat_xSecondary,Lon_xSecondary);

    //print the variable list
    fprintf(fout2d,"VARIABLES=\"Lon\", \"Lat\"");

    for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
      double r0,r1;

      r0=ElectricallyChargedDust::minDustRadius*exp(SizeGroup*ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals);
      r1=ElectricallyChargedDust::minDustRadius*exp((1+SizeGroup)*ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals);

      fprintf(fout2d,", \"Flux [s^{-1}m^{-2}] (%e<a<%e)\"",r0,r1);
    }

    fprintf(fout2d,", \"Total Flux [s^{-1}m^{-2}]\"");
    fprintf(fout2d,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthalSurfaceElements+2,nZenithSurfaceElements+1);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;
    double lon,lat;
    int iAzimuthalBreak=-1;

    do {
      iAzimuthalBreak++;
      GetSurfaceLonLatNormal(lon,lat,0,iAzimuthalBreak);
    }
    while (lon<=180.0);

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) {
      GetSurfaceLonLatNormal(lon,lat,iZenith,0);

      for (int iPass=0;iPass<2;iPass++) {
         int iAzimuthalMin,iAzimuthalMax;

         switch(iPass) {
         case 0:
           iAzimuthalMin=iAzimuthalBreak-1,iAzimuthalMax=nAzimuthalSurfaceElements-1;
           break;
         case 1:
           iAzimuthalMin=0,iAzimuthalMax=iAzimuthalBreak;
           break;
         default:
           exit(__LINE__,__FILE__,"Error: wrong option");
         }

        //first output the data for azimuthal angle -180 to 0 and than that from 0 to 180 degrees
        for (iAzimuthal=iAzimuthalMin;iAzimuthal<=iAzimuthalMax;iAzimuthal++) {
          lon=(iPass==1) ? iAzimuthal*dAzimuthalAngle*180.0/Pi : iAzimuthal*dAzimuthalAngle*180.0/Pi-360.0;
          fprintf(fout2d,"%e %e ",lon,lat);

          if (PrintStateVectorFlag==true) {
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
            //sampled data summed for all dust sizes
            double TotalFlux=0.0;

            for (SizeGroup=0;SizeGroup<nDustSizeSamplingIntervals;SizeGroup++) {
              int el,nInterpolationElement;
              double InterpolationCoefficient,summInterpolationCoefficient=0.0;
              double Flux=0.0;

              for (nInterpolationElement=0;nInterpolationElement<InterpolationListLength;nInterpolationElement++) {
                el=InterpolationList[nInterpolationElement];
                InterpolationCoefficient=GetSurfaceElementArea(el);
                summInterpolationCoefficient+=InterpolationCoefficient;

                //there is no need to multiply by the surface element surface because of the conversion from the flux
                //through the surface elemet into the flux per m^{-2}
                Flux+=InterpolationCoefficient*SampleData_NumberDensityFlux[SizeGroup][el];
              }

              //normalize the interpolated values
              if (PIC::LastSampleLength!=0) {
                Flux/=summInterpolationCoefficient*PIC::LastSampleLength;
              }

              //output the interpolated values into a file
              TotalFlux+=Flux;

              fprintf(fout2d," %e",Flux);
            }

            fprintf(fout2d," %e",TotalFlux);
          }

          fprintf(fout2d,"\n");
        }
      }
    }

    fclose(fout2d);
  }
}

double ElectricallyChargedDust::Sampling::FluxMap::cSampleLocation::GetSpeed(double *v,double &ZenithAngle,long int &nZenithElement, double &AzimuthalAngle,long int &nAzimuthalElement) {
  double r;
  int idim;

  //Normalize the velocity vector and determine its coordinates in the frame related for the observation point
  double Speed=0.0,l[3]={0.0,0.0,0.0};

  for (idim=0;idim<3;idim++) Speed+=pow(v[idim],2);
  Speed=sqrt(Speed);

  //convert 'v' into the reference frame related to the sa,mple location
  for (int idim=0;idim<3;idim++) {
    l[0]+=e0[idim]*v[idim]/Speed;
    l[1]+=e1[idim]*v[idim]/Speed;
    l[2]+=e2[idim]*v[idim]/Speed;
  }

  GetSurfaceElementProjectionIndex(l,ZenithAngle,nZenithElement,AzimuthalAngle,nAzimuthalElement);

  return Speed;
}



