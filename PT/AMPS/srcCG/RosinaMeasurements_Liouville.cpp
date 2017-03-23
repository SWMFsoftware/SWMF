//$Id$
//Evaluate density and pressure along the spacecraft trajectory using the Liouville theoreme

/*
 * RosinaMeasurements_Liouville.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "RosinaMeasurements.h"
#include "Comet.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

extern double *productionDistributionUniformNASTRAN;
extern double HeliocentricDistance,subSolarPointAzimuth,subSolarPointZenith;
extern bool *definedFluxBjorn;
extern double positionSun[3];
extern double **productionDistributionNASTRAN;

void RosinaSample::Liouville::EvaluateLocation(int spec,double& OriginalSourceRate, double& ModifiedSourceRate,double& NudeGaugePressure,double& NudeGaugeDensity,double& NudeGaugeFlux, double& RamGaugePressure,double& RamGaugeDensity,double& RamGaugeFlux,int iPoint) {
  int iTest,iSurfaceElement;
  double c,l[3],*xLocation,rLocation,xIntersection[3],beta;
  double SampledNudeGaugeDensity=0.0,SampledRamGaugeFlux=0.0;

  const int nTotalTests=_ROSINA_SAMPLE__LIOUVILLE__EVALUATE_LOCATION__TEST_NUMBER_;
  const double ThrehondSourceRate=1.8E18;

  //estimate the total source rate
  OriginalSourceRate=0.0,ModifiedSourceRate=0.0;

  for (int iSurfaceElement=0;iSurfaceElement<CutCell::nBoundaryTriangleFaces;iSurfaceElement++) {
    double SourceRate=productionDistributionNASTRAN[spec][iSurfaceElement]/CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    OriginalSourceRate+=SourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
    ModifiedSourceRate+=((SourceRate>ThrehondSourceRate) ? SourceRate : ThrehondSourceRate)*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    //evaluate the cosine between the external normal of teh face and the pointing to the spacecraft
    double x[3],l=0.0,c=0.0;
    CutCell::BoundaryTriangleFaces[iSurfaceElement].GetCenterPosition(x);

    for (int idim=0;idim<3;idim++) {
      c+=CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal[idim]*(Rosina[iPoint].x[idim]-x[idim]);
      l+=pow(Rosina[iPoint].x[idim]-x[idim],2);
    }

    CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.CrossProduct_FaceNormal_SpacecraftLocation=c/sqrt(l);
  }


  //set the angular limit for the ray direction generation
  xLocation=Rosina[iPoint].x;

  NudeGaugeDensity=0.0,NudeGaugeFlux=0.0;
  RamGaugeFlux=0.0,RamGaugeDensity=0.0;

  //loop through all surface elements
  int iStartSurfaceElement,iFinishSurfaceElement,nSurfaceElementThread;

  nSurfaceElementThread=CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads;
  iStartSurfaceElement=nSurfaceElementThread*PIC::ThisThread;
  iFinishSurfaceElement=iStartSurfaceElement+nSurfaceElementThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iFinishSurfaceElement=CutCell::nBoundaryTriangleFaces;

  double localNudeGaugeDensity=0.0,localNudeGaugeFlux=0.0,localRamGaugeFlux=0.0,localRamGaugeDensity=0.0;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic) default(none) shared(iPoint,iStartSurfaceElement,iFinishSurfaceElement,Rosina,xLocation,CutCell::BoundaryTriangleFaces,productionDistributionNASTRAN,spec,positionSun) \
  private(l,xIntersection,beta) reduction(+:localNudeGaugeDensity) reduction(+:localNudeGaugeFlux) reduction(+:localRamGaugeFlux) reduction(+:localRamGaugeDensity)
#endif
  for (iSurfaceElement=iStartSurfaceElement;iSurfaceElement<iFinishSurfaceElement;iSurfaceElement++) {
    int iTest,idim;
    double c=0.0,x[3],l[3],r,cosTheta,A;
    double tNudeGaugeDensity=0.0,tRamGaugeFlux=0.0,tRamGaugeDensity=0.0;

    CutCell::BoundaryTriangleFaces[iSurfaceElement].GetCenterPosition(x);

    for (idim=0;idim<3;idim++) c+=CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal[idim]*(xLocation[idim]-x[idim]);

    if (c<0.0) {
      //the normal is directed opposite to the point of data sampling -> check the next face
      continue;
    }
    else {
      //determine the surface temeprature, production rate, and the normalization coefficient
      //parameters of the primary species source at that surface element
      double SourceRate,Temperature,cosSubSolarAngle,x_LOCAL_SO_OBJECT[3]={0.0,0.0,0.0};

      SourceRate=productionDistributionNASTRAN[spec][iSurfaceElement]/CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

      if (SourceRate<ThrehondSourceRate) {
        productionDistributionNASTRAN[spec][iSurfaceElement]=ThrehondSourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
        SourceRate=ThrehondSourceRate;
      }

      cosSubSolarAngle=Vector3D::DotProduct(CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal,positionSun)/Vector3D::Length(positionSun);
      if (CutCell::BoundaryTriangleFaces[iSurfaceElement].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow

      Temperature=Comet::GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);

      beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*Temperature));
      A=2.0*SourceRate/Pi*pow(beta,4);
    }

    for (iTest=0;iTest<nTotalTests;iTest++) {
      CutCell::BoundaryTriangleFaces[iSurfaceElement].GetRandomPosition(x);
      for (idim=0;idim<3;idim++) l[idim]=xLocation[idim]-x[idim];

      if (PIC::RayTracing::FindFistIntersectedFace(x,l,xIntersection,CutCell::BoundaryTriangleFaces+iSurfaceElement)==-1) {
        //there is the direct access from the point on teh surface to the point of the observation ->  sample the number density and flux
        double ExternalNormal[3];

        switch (_ROSINA_SAMPLE__LIOUVILLE__INJECTED_GAS_DIRECTION_DISTRIBUTION_MODE_) {
        case _ROSINA_SAMPLE__LIOUVILLE__INJECTED_GAS_DIRECTION_DISTRIBUTION_MODE__NORMAL_:
          memcpy(ExternalNormal,CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal,3*sizeof(double));
          break;
        case _ROSINA_SAMPLE__LIOUVILLE__INJECTED_GAS_DIRECTION_DISTRIBUTION_MODE__RANDOM_:
          do {
            Vector3D::Distribution::Uniform(ExternalNormal);
          }
          while (Vector3D::DotProduct(ExternalNormal,CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal)<0.0);

          Vector3D::Normalize(ExternalNormal);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is not found");
        }

        r=Vector3D::Length(l);
        cosTheta=Vector3D::DotProduct(l,ExternalNormal)/r;

        if (Vector3D::DotProduct(l,Rosina[iPoint].NudeGauge.LineOfSight)<0.0) {
          //the particle flux can access the nude gauge
          tNudeGaugeDensity+=cosTheta/pow(r,2)/pow(beta,3);
        }

        if ((c=Vector3D::DotProduct(l,Rosina[iPoint].RamGauge.LineOfSight))<0.0)  {
           //the particle flux can bedetected by the ram gauge
          double cosLineOfSightAngle=-c/r;

          tRamGaugeFlux+=cosTheta/(2.0*pow(r,2)*pow(beta,4)) * cosLineOfSightAngle;
          tRamGaugeDensity+=cosTheta/pow(r,2)/pow(beta,3);
        }
      }
    }

    //sample the contribution of the surface element to the isntrument observation
    double t,*t_ptr;

    //nude gauge density
    t=tNudeGaugeDensity * A*sqrt(Pi)/4.0 / nTotalTests;
    localNudeGaugeDensity+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    t_ptr=CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.NudeGaugeDensityContribution;

    #pragma omp atomic
    t_ptr[spec]+=t;
    #else
    CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.NudeGaugeDensityContribution[spec]+=t;
    #endif

    //nude gauge flux
    localNudeGaugeFlux=0.0;

    //ram gauge density
    t=tRamGaugeDensity * A*sqrt(Pi)/4.0 / nTotalTests;
    localRamGaugeDensity+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    t_ptr=CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeDensityContribution;

    #pragma omp atomic
    t_ptr[spec]+=t;
    #else
    CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeDensityContribution[spec]+=t;
    #endif


    //ram gauge flux
    t=tRamGaugeFlux*A/2.0 / nTotalTests;
    localRamGaugeFlux+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    t_ptr=CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeFluxContribution;

    #pragma omp atomic
    t_ptr[spec]+=t;
    #else
    CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeFluxContribution[spec]+=t;
    #endif



  }

  //colles contribution from all processors
  MPI_Allreduce(&localNudeGaugeDensity,&NudeGaugeDensity,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localNudeGaugeFlux,&NudeGaugeFlux,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  MPI_Allreduce(&localRamGaugeFlux,&RamGaugeFlux,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localRamGaugeDensity,&RamGaugeDensity,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  //convert the sampled fluxes and density into the nude and ram gauges pressures
  const double rgTemperature=293.0;
  const double ngTemperature=293.0;

  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*rgTemperature));

//  RamGaugePressure=Kbol*rgTemperature*Pi*beta/2.0*RamGaugeFlux;
  RamGaugePressure=4.0*Kbol*rgTemperature*RamGaugeFlux/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(spec)));

  NudeGaugePressure=NudeGaugeDensity*Kbol*ngTemperature;


/* EXAMPLE OF CALCULATING THE PRESSURES
   if (_H2O_SPEC_>=0) {
    beta=sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/(2.0*Kbol*rgTemperature));
    rgTotalPressure+=Kbol*rgTemperature*Pi*beta/2.0*FluxRamGauge[_H2O_SPEC_+i*PIC::nTotalSpecies];

//        rgTotalPressure+=4.0*Kbol*rgTemperature*FluxRamGauge[_H2O_SPEC_+i*PIC::nTotalSpecies]/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(_H2O_SPEC_)));
    ngTotalPressure+=DensityNudeGauge[_H2O_SPEC_+i*PIC::nTotalSpecies]*Kbol*ngTemperature;
  }*/

}

//==========================================================================================
//estimate the solud angles that the nucleus is seen by the ram and nude gauges
void RosinaSample::Liouville::GetSolidAngle(double& NudeGaugeNucleusSolidAngle,double& RamGaugeNucleusSolidAngle,int iPoint) {
  NudeGaugeNucleusSolidAngle=0.0;
  RamGaugeNucleusSolidAngle=0.0;

  //solid angle occupied by the nucleus by the nude gauge
  int t,nTotalTests,iTest,iIntersectionFace,NucleusIntersectionCounter=0;
  double l[3],xIntersection[3];

  int iStartTest,iFinishTest,nTestThread;

  nTotalTests=_ROSINA_SAMPLE__LIOUVILLE__GET_SOLID_ANGLE__TEST_NUMBER_;
  iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
  iFinishTest=(nTotalTests/PIC::nTotalThreads)*(PIC::ThisThread+1);

  nTestThread=nTotalTests/PIC::nTotalThreads;
  iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
  iFinishTest=iStartTest+nTestThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iFinishTest=nTotalTests;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic) default(none) shared(CutCell::BoundaryTriangleFaces,iStartTest,iFinishTest,Rosina,iPoint) private(iTest,l,xIntersection,iIntersectionFace) reduction(+:NucleusIntersectionCounter)
#endif
  for (iTest=iStartTest;iTest<iFinishTest;iTest++) {
    //generate a ranfom direction
    do {
      Vector3D::Distribution::Uniform(l);
    }
    while (Vector3D::DotProduct(l,Rosina[iPoint].NudeGauge.LineOfSight)<0.0);

    //determine whether an intersection with the nucleus is found
    iIntersectionFace=PIC::RayTracing::FindFistIntersectedFace(Rosina[iPoint].x,l,xIntersection,NULL);

    if (iIntersectionFace!=-1) {
      NucleusIntersectionCounter++;
      CutCell::BoundaryTriangleFaces[iIntersectionFace].UserData.FieldOfView_NudeGauge=true;
    }
  }

  MPI_Allreduce(&NucleusIntersectionCounter,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  NudeGaugeNucleusSolidAngle=2.0*Pi*double(t)/double(nTotalTests);

  //solid angle occupied by the nucleus by the ram gauge
  NucleusIntersectionCounter=0;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic) default(none) shared(CutCell::BoundaryTriangleFaces,iStartTest,iFinishTest,Rosina,iPoint) private(iTest,l,xIntersection,iIntersectionFace) reduction(+:NucleusIntersectionCounter)
#endif
  for (iTest=iStartTest;iTest<iFinishTest;iTest++) {
    //generate a ranfom direction
    do {
      Vector3D::Distribution::Uniform(l);
    }
    while (Vector3D::DotProduct(l,Rosina[iPoint].RamGauge.LineOfSight)<0.0);

    //determine whether an intersection with the nucleus is found
    iIntersectionFace=PIC::RayTracing::FindFistIntersectedFace(Rosina[iPoint].x,l,xIntersection,NULL);

    if (iIntersectionFace!=-1) {
      NucleusIntersectionCounter++;
      CutCell::BoundaryTriangleFaces[iIntersectionFace].UserData.FieldOfView_RamGauge=true;
    }
  }

  MPI_Allreduce(&NucleusIntersectionCounter,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  RamGaugeNucleusSolidAngle=2.0*Pi*double(t)/double(nTotalTests);
}


//==========================================================================================
//evaluete the fulle set of the measuremetns
void RosinaSample::Liouville::Evaluate() {
#ifndef _NO_SPICE_CALLS_  //the function will be compiled only when SPICE is used (bacause SPICE is needed to adjust the model parameters for each point of the observations
  double NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux;
  double NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle;
  double OriginalSourceRate,ModifiedSourceRate;
  int iPoint,idim,spec;
  SpiceDouble lt,et,xRosetta[3],etStart;
  SpiceDouble       xform[6][6];

  FILE *fout[PIC::nTotalSpecies];
  FILE *fGroundTrack=NULL;
  FILE *fAllSpecies;

  const int Step=12;
  const int SurfaceOutputSter=2;

  if (PIC::ThisThread==0) {
    char fname[200];

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      sprintf(fname,"Liouville.spec=%s.dat",PIC::MolecularData::GetChemSymbol(spec));
      fout[spec]=fopen(fname,"w");
      fprintf(fout[spec],"VARIABLES=\"i\", \"Nude Gauge Pressure\", \"Nude Gauge Density\", \"Nude Gauge Flux\",  \"Ram Gauge Pressure\", \"Ram Gauge Density\", \"Ram Gauge Flux\", \"Seconds From The First Point\", \"Nude Guage Nucleus Solid angle\", \"Ram Gauge Nucleus Solid Angle\", \"Altitude\", \"Closest Surface Element Source Rate [m^-2 s^-1]\", \"Nude Gauge COPS Measurements\", \"Ram Gauge COPS Measurements\", \"Original Total Source Rate\", \"Modified Total Source Rate\" \n");
    }

    fAllSpecies=fopen("Liouville.spec=SUM.dat","w");
    fprintf(fAllSpecies,"VARIABLES=\"i\", \"Nude Gauge Pressure\", \"Nude Gauge Density\", \"Nude Gauge Flux\",  \"Ram Gauge Pressure\", \"Ram Gauge Density\", \"Ram Gauge Flux\", \"Seconds From The First Point\", \"Nude Guage Nucleus Solid angle\", \"Ram Gauge Nucleus Solid Angle\", \"Altitude\", \"Closest Surface Element Source Rate [m^-2 s^-1]\", \"Nude Gauge COPS Measurements\", \"Ram Gauge COPS Measurements\", \"Original Total Source Rate\", \"Modified Total Source Rate\" \n");

    fGroundTrack=fopen("GroundTracks.dat","w");
    printf("VARIABLES=\"i\", \"Nude Gauge Pressure\", \"Nude Gauge Density\", \"Nude Gauge Flux\",  \"Ram Gauge Pressure\", \"Ram Gauge Density\", \"Ram Gauge Flux\", \"Seconds From The First Point\", \"Nude Guage Nucleus Solid angle\", \"Ram Gauge Nucleus Solid Angle\", \"Altitude\", \"Closest Surface Element Source Rate [m^-2 s^-1]\", \"Nude Gauge COPS Measurements\", \"Ram Gauge COPS Measurements\", \"Original Total Source Rate\", \"Modified Total Source Rate\" \n");
  }

  for (iPoint=0;iPoint<RosinaSample::nPoints;iPoint+=Step) {
    //set the vectors, location of the Sub for the observation
    //init line-of-sight vectors
    //set the location of the spacecraft
    utc2et_c(RosinaSample::ObservationTime[iPoint],&et);
    spkpos_c("ROSETTA",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xRosetta,&lt);
    for (idim=0;idim<3;idim++) xRosetta[idim]*=1.0E3;

    if (iPoint==0) etStart=et;

    Rosina[iPoint].SecondsFromBegining=et-etStart;
    for (idim=0;idim<3;idim++) Rosina[iPoint].x[idim]=xRosetta[idim];

    sxform_c("ROS_SPACECRAFT","67P/C-G_CK",et,xform);

    //set pointing of COPS nude gauge: negative y
    SpiceDouble NudeGauge[6]={0.0,-1.0,0.0 ,0.0,0.0,0.0};
    SpiceDouble l[6];

    mxvg_c(xform,NudeGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[iPoint].NudeGauge.LineOfSight[idim]=l[idim];


    //set pointing of COPS ram gauge: positive z
    SpiceDouble RamGauge[6]={0.0,0.0,1.0 ,0.0,0.0,0.0};

    mxvg_c(xform,RamGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[iPoint].RamGauge.LineOfSight[idim]=l[idim];

    //update the location of the Sun and the nucleus shadowing
    SpiceDouble xSun[3];

    spkpos_c("SUN",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xSun,&lt);
    reclat_c(xSun,&HeliocentricDistance,&subSolarPointAzimuth,&subSolarPointZenith);

    HeliocentricDistance*=1.0E3;

    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

    PIC::RayTracing::SetCutCellShadowAttribute(positionSun,false);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      definedFluxBjorn[spec]=false;
      Comet::GetTotalProductionRateBjornNASTRAN(spec);
    }


    Comet::BjornNASTRAN::Init();


    //reset the fiew of view flags
    for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
      CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_NudeGauge=false;
      CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_RamGauge=false;
    }

    //simulate the solid angles and the measuremetns
    GetSolidAngle(NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,iPoint);

    //evaluate the location
    double TotalOriginalSourceRate=0.0,TotalModifiedSourceRate=0.0,TotalNudeGaugePressure=0.0;
    double TotalNudeGaugeDensity=0.0,TotalNudeGaugeFlux=0.0,TotalRamGaugePressure=0.0,TotalRamGaugeDensity=0.0,TotalRamGaugeFlux=0.0;



    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      //save the original source rate of the species
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++)  {
        CutCell::BoundaryTriangleFaces[iface].UserData.OriginalSourceRate[spec]=productionDistributionNASTRAN[spec][iface]/CutCell::BoundaryTriangleFaces[iface].SurfaceArea;
      }

      //evaluate the location
      EvaluateLocation(spec,OriginalSourceRate,ModifiedSourceRate,NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux,iPoint);

      TotalOriginalSourceRate+=OriginalSourceRate;
      TotalModifiedSourceRate+=ModifiedSourceRate;
      TotalNudeGaugePressure+=NudeGaugePressure;
      TotalNudeGaugeDensity+=NudeGaugeDensity;
      TotalNudeGaugeFlux+=NudeGaugeFlux;
      TotalRamGaugePressure+=RamGaugePressure;
      TotalRamGaugeDensity+=RamGaugeDensity;
      TotalRamGaugeFlux+=RamGaugeFlux;

      //save the modified source rate of the species
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++)  {
        CutCell::BoundaryTriangleFaces[iface].UserData.ModifiedSourceRate[spec]=productionDistributionNASTRAN[spec][iface]/CutCell::BoundaryTriangleFaces[iface].SurfaceArea;
      }


      if (PIC::ThisThread==0) {
        fprintf(fGroundTrack,"%e %e %e\n",Rosina[iPoint].xNucleusClosestPoint[0],Rosina[iPoint].xNucleusClosestPoint[1],Rosina[iPoint].xNucleusClosestPoint[2]);

        fprintf(fout[spec],"%i %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",iPoint,NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux,
            Rosina[iPoint].SecondsFromBegining,
            NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,
            Rosina[iPoint].Altitude,
            ((Rosina[iPoint].Altitude>0.0) ? productionDistributionNASTRAN[spec][Rosina[iPoint].iNucleusClosestFace]/CutCell::BoundaryTriangleFaces[Rosina[iPoint].iNucleusClosestFace].SurfaceArea : 0.0),
            RosinaSample::NudeGaugeReferenceData[iPoint],RosinaSample::RamGaugeReferenceData[iPoint],
            OriginalSourceRate,ModifiedSourceRate);


        printf("%i (%s) %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",iPoint,PIC::MolecularData::GetChemSymbol(spec),NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux,
            Rosina[iPoint].SecondsFromBegining,
            NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,
            Rosina[iPoint].Altitude,
            ((Rosina[iPoint].Altitude>0.0) ? productionDistributionNASTRAN[spec][Rosina[iPoint].iNucleusClosestFace]/CutCell::BoundaryTriangleFaces[Rosina[iPoint].iNucleusClosestFace].SurfaceArea : 0.0),
            RosinaSample::NudeGaugeReferenceData[iPoint],RosinaSample::RamGaugeReferenceData[iPoint],
            OriginalSourceRate,ModifiedSourceRate);
      }
    }

    if (PIC::ThisThread==0) {
      fprintf(fAllSpecies,"%i %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",iPoint,TotalNudeGaugePressure,TotalNudeGaugeDensity,TotalNudeGaugeFlux,TotalRamGaugePressure,TotalRamGaugeDensity,TotalRamGaugeFlux,
          Rosina[iPoint].SecondsFromBegining,
          NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,
          Rosina[iPoint].Altitude,
          ((Rosina[iPoint].Altitude>0.0) ? (productionDistributionNASTRAN[_H2O_SPEC_][Rosina[iPoint].iNucleusClosestFace]+productionDistributionNASTRAN[_CO2_SPEC_][Rosina[iPoint].iNucleusClosestFace])/CutCell::BoundaryTriangleFaces[Rosina[iPoint].iNucleusClosestFace].SurfaceArea : 0.0),
          RosinaSample::NudeGaugeReferenceData[iPoint],RosinaSample::RamGaugeReferenceData[iPoint],
          OriginalSourceRate,ModifiedSourceRate);
    }


    //save the surface properties
    if ((iPoint/Step)%SurfaceOutputSter==0) {
      //create the surface output file, and clean the buffers
      char fname[200];

      sprintf(fname,"%s/SurfaceContributionParameters.iPoint=%i.dat",PIC::OutputDataFileDirectory,iPoint);
      CutCell::PrintSurfaceData(fname);

      //clear the samplign buffers
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        CutCell::BoundaryTriangleFaces[iface].UserData.NudeGaugeDensityContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.NudeGaugeFluxContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.RamGaugeDensityContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.RamGaugeFluxContribution[spec]=0.0;
      }
    }
  }

  if (PIC::ThisThread==0) {
    fclose(fGroundTrack);
    fclose(fAllSpecies);
    for (spec=0;spec<PIC::nTotalSpecies;spec++) fclose(fout[spec]);
  }
#endif //_NO_SPICE_CALLS_
}



