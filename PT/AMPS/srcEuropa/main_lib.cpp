
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>


#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include <saturn.h>

//#define _MAIN_PROFILE_



//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "Europa.h"

//the modeling case
#define _EUROPA_MODEL_MODE__NEAR_PLANET_    0
#define _EUROPA_MODEL_MODE__TAIL_MODEL_     1

#define _EUROPA_MODEL_MODE_ _EUROPA_MODEL_MODE__TAIL_MODEL_

//Check the mesh consistency
//#define _CHECK_MESH_CONSISTENCY_

//create the nre list of points for ICES
#define _ICES_CREATE_COORDINATE_LIST_  _PIC_MODE_ON_


#define _ICES_LOAD_DATA_


//defiend the modes for included models
#define _EUROPA_MODE_ON_    0
#define _EUROPA_MODE_OFF_   1

#define _EUROPA_IMPACT_VAPORIZATION_MODE_ _EUROPA_MODE_ON_
#define _EUROPA_PSD_MODE_ _EUROPA_MODE_ON_
#define _EUROPA_THERMAL_DESORPTION_MODE_ _EUROPA_MODE_ON_

//sample particles that enter EPD
#define _GALILEO_EPD_SAMPLING_  _EUROPA_MODE_ON_



//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;


const double rSphere=_RADIUS_(_TARGET_);

//const double xMaxDomain=400.0; //modeling of the tail
const double xMaxDomain=5; //modeling the vicinity of the planet
const double yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


//const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*10.0/100,dxMaxSphere=DebugRunMultiplier*1.0/10.0;

const double dxMinGlobal=1,dxMaxGlobal=1;

//the species
//int NA=0;
//int NAPLUS=1;


//the total acceleration of the particles
#include "Na.h"

//photoionization lifetime of sodium
/*
double sodiumPhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {



  return 1.0E7;

	static const double LifeTime=3600.0*5.8;

	double res,r2=x[1]*x[1]+x[2]*x[2];

	if ((r2>rSphere*rSphere)||(x[0]<0.0)) {
		res=LifeTime,PhotolyticReactionAllowedFlag=true;
	}
	else {
		res=-1.0,PhotolyticReactionAllowedFlag=false;
	}

	return res;
}

int sodiumPhotoionizationReactionProcessor(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {
	spec=_O2_PLUS_SPEC_;

	PIC::ParticleBuffer::SetI(spec,ParticleData);
	  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;

	//return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
}
*/


//the codes for the soruces processes
#define _ALL_SOURCE_PROCESSES_                -1
#define _IMPACT_VAPORIZATION_SOURCE_PROCESS_   0
#define _PSD_SOURCE_PROCESS_                   1

//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
	double res=0.0;

#if _EUROPA_IMPACT_VAPORIZATION_MODE_ == _EUROPA_MODE_ON_
	if ((SourceProcessCode==-1)||(SourceProcessCode==_IMPACT_VAPORIZATION_SOURCE_PROCESS_)) { //the total impact vaporization flux
		res+=2.6E23;
	}
#endif

#if _EUROPA_PSD_MODE_ == _EUROPA_MODE_ON_
	if ((SourceProcessCode==-1)||(SourceProcessCode==_PSD_SOURCE_PROCESS_)) { //the photon stimulated desorption flux
		res+=3.6E24;
	}
#endif

	return res;
}





//impact vaporization
//const double ProductionRateIVsodium=2.6E23;
//const double TemperatureInjectionIVsodium=2000.0;



//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
        double res,r,l[3] = {1.0,0.0,0.0};
	int idim;
	double SubsolarAngle;

	for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
	if (r > 0.8 * rSphere)
	  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

	SubsolarAngle=acos(l[0]);


	SubsolarAngle=0.0;

	res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;


res/=2.1;

	return rSphere*res;
}

double localResolution(double *x) {
	int idim;
	double lnR,res,r=0.0;

	for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

	r=sqrt(r);

	if (r<2.0*rSphere) return localSphericalSurfaceResolution(x);

	if (r>dxMinGlobal*rSphere) {
		lnR=log(r);
		res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
	}
	else res=dxMinGlobal;

	//  if ((x[0]>0.0)&&(sqrt(x[1]*x[1]+x[2]*x[2])<3.0*rSphere)) res=(res<0.4) ? res : 0.4;  ///min(res,0.4);

	return rSphere*res;
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	double CellSize;

//	double CharacteristicSpeed_NA=5.0E3;

	const double maxInjectionEnergy=1E5*1000*ElectronCharge;
	double CharacteristicSpeed = sqrt(2*maxInjectionEnergy/_MASS_(_O_));


	//  CharacteristicSpeed*=sqrt(PIC::MolecularData::GetMass(NA)/PIC::MolecularData::GetMass(spec));

	CellSize=startNode->GetCharacteristicCellSize();




//	CharacteristicSpeed=1.0E6;

/*if (spec==_OPLUS_HIGH_SPEC_) CharacteristicSpeed=10.0*1.6E6;
if (spec==_OPLUS_THERMAL_SPEC_) CharacteristicSpeed=10.0*9.6E4;*/

  switch (spec) {
  case _OPLUS_HIGH_SPEC_:
    CharacteristicSpeed=10.0*1.6E6;
    break;
  case _OPLUS_THERMAL_SPEC_:
    CharacteristicSpeed=10.0*9.6E4;
    break;

  case _O2_SPEC_:
    CharacteristicSpeed=1.0e4;
    break;
  case _O2PLUS_SPEC_:
    CharacteristicSpeed=10*1.0e4;
    break;
  default:
    exit(__LINE__,__FILE__,"unknown species");
   }

	return 0.3*CellSize/CharacteristicSpeed;


}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

	/*
  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
	 */

	return 0.0;

	/*
  if (spec!=NA) return 0.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
	 */
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	bool ExternalFaces[6];
	double ExternalNormal[3],ModelParticlesInjectionRate;
	int nface;

	static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;

	if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
		for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
			startNode->GetExternalNormal(ExternalNormal,nface);
			ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_OPLUS_HIGH_SPEC_);

			if (ModelParticlesInjectionRate>0.0) return true;
		}
	}

	return false;
}



//injection of model particles through the faces of the bounding box
long int  BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	bool ExternalFaces[6];
	double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
	int nface,idim;
	//  long int nInjectedParticles;
	long int newParticle;
	PIC::ParticleBuffer::byte *newParticleData;
	long int nInjectedParticles=0;

	if ((spec!=_OPLUS_HIGH_SPEC_)&&(spec!=_OPLUS_THERMAL_SPEC_)) return 0; //inject only spec=0

	static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
	double v[3];


	double ModelParticlesInjectionRate;

	if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
		ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
		LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


		for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
			startNode->GetExternalNormal(ExternalNormal,nface);
			TimeCounter=0.0;

			ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_OPLUS_THERMAL_SPEC_);


			if (ModelParticlesInjectionRate>0.0) {
				ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

				PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

				while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
					//generate the new particle position on the face
					for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

					//generate a particle
					newParticle=PIC::ParticleBuffer::GetNewParticle();
					newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
					nInjectedParticles++;

					PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_OPLUS_THERMAL_SPEC_);

					PIC::ParticleBuffer::SetX(x,newParticleData);
					PIC::ParticleBuffer::SetV(v,newParticleData);
					PIC::ParticleBuffer::SetI(spec,newParticleData);

					PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);


					//inject the particle into the system
					//          PIC::Mover::MoveParticleTimeStep[spec](newParticle,LocalTimeStep-TimeCounter,startNode);

				}
			}


		}
	}

	return nInjectedParticles;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	long int nInjectedParticles=0;

	for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

	return nInjectedParticles;
}






double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
	double res=1.0;

	// for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

	return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
	double radiusSphere,*x0Sphere,l[3],r,vNorm,c;
	cInternalSphericalData *Sphere;
	cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
	int idim;

	//   long int newParticle;
	//   PIC::ParticleBuffer::byte *newParticleData;
	//   double ParticleStatWeight,WeightCorrection;


	Sphere=(cInternalSphericalData*)SphereDataPointer;
	startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

	Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

	for (r=0.0,idim=0;idim<DIM;idim++) {
		l[idim]=x[idim]-x0Sphere[idim];
		r+=pow(l[idim],2);
	}

	for (r=sqrt(r),vNorm=0.0,idim=0;idim<DIM;idim++) vNorm+=v[idim]*l[idim]/r;
	if (vNorm<0.0) for (c=2.0*vNorm/r,idim=0;idim<DIM;idim++) v[idim]-=c*l[idim];

	//sample the particle data
	double *SampleData;
	long int nSurfaceElement,nZenithElement,nAzimuthalElement;

	Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
	nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
	SampleData=Sphere->SamplingBuffer+PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


	SampleData[PIC::BC::InternalBoundary::Sphere::sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);


	//   if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<0.9*rSphere*rSphere) exit(__LINE__,__FILE__,"Particle inside the sphere");


	r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);


	//particle-surface interaction
	if (false) { /////(spec!=NA) { //surface reactiona
		exit(__LINE__,__FILE__,"no BC for the space is implemented");

		/*
     ParticleStatWeight=startNode->block->GetLocalParticleWeight(0);
     ParticleStatWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

     //model the rejected H+ and neutralized H
     //1. model rejected SW protons (SPEC=1)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(1,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.01*ParticleStatWeight/startNode->block->GetLocalParticleWeight(1);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[1](newParticle,dtTotal,startNode);

     //1. model rejected SW NEUTRALIZED protons (SPEC=2)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(2,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.2*ParticleStatWeight/startNode->block->GetLocalParticleWeight(2);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[2](newParticle,dtTotal,startNode);

     //remove the oroginal particle (s=0)
     PIC::ParticleBuffer::DeleteParticle(ptr);
     return _PARTICLE_DELETED_ON_THE_FACE_;
		 */
	}


	//delete all particles that was not reflected on the surface
//	PIC::ParticleBuffer::DeleteParticle(ptr);
	return _PARTICLE_DELETED_ON_THE_FACE_;
}




double sphereInjectionRate(int spec,void *SphereDataPointer) {
	double res=0.0;

	/*
	if (spec==NA) res=sodiumTotalProductionRate();
	else if (spec==NAPLUS) res=0.0;
	else exit(__LINE__,__FILE__,"Error: the source rate for the species is not determined");
*/

	return res;
}



void amps_init() {
//	MPI_Init(&argc,&argv);


PIC::InitMPI();


	//SetUp the alarm
	//  PIC::Alarm::SetAlarm(2000);



	//  VT_OFF();
	//  VT_traceoff();


#ifdef _MAIN_PROFILE_
	initSaturn ("");
#endif


	rnd_seed();



	char inputFile[]="Europa.input";



	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


	//init the Europa model
	//Exosphere::Init_BeforeParser();
	Europa::Init_BeforeParser();

	//init the particle solver
	PIC::Init_BeforeParser();



//	PIC::Parser::Run(inputFile);


	const int InitialSampleLength=100;

	//PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=1600;
	//PIC::RequiredSampleLength=InitialSampleLength; //0;

	Europa::Init_AfterParser();



	//output the PDS energy distribution function
	if (PIC::ThisThread==0) {
		//Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintCumulativeDistributionFunction("CumulativeEnergyDistribution-PSD.dat");
		//Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintDistributionFunction("EnergyDistribution-PSD.dat");

		//cout << Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
	}



	//register the sphere
	static const bool SphereInsideDomain=true;

	if (SphereInsideDomain==true) {
		double sx0[3]={0.0,0.0,0.0};
		cInternalBoundaryConditionsDescriptor SphereDescriptor;
		cInternalSphericalData *Sphere;


		//reserve memory for sampling of the surface balance of sticking species
		long int ReserveSamplingSpace[PIC::nTotalSpecies];

		for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_EUROPA_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


		cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);



		PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
		SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
		Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
		Sphere->SetSphereGeometricalParameters(sx0,rSphere);


		//init the object for distribution of the injection surface elements
    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
      Europa::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


    #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      Europa::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
    #endif

    #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      Europa::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
    #endif
    }






    Sphere->Radius=_RADIUS_(_TARGET_);
		Sphere->PrintSurfaceMesh("Sphere.dat");
		Sphere->PrintSurfaceData("SpheraData.dat",0);
		Sphere->localResolution=localSphericalSurfaceResolution;
		Sphere->InjectionRate=Europa::SourceProcesses::totalProductionRate;
		Sphere->faceat=0;
		Sphere->ParticleSphereInteraction=Europa::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
		Sphere->InjectionBoundaryCondition=Europa::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;

		Sphere->PrintTitle=Europa::Sampling::OutputSurfaceDataFile::PrintTitle;
		Sphere->PrintVariableList=Europa::Sampling::OutputSurfaceDataFile::PrintVariableList;
		Sphere->PrintDataStateVector=Europa::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

		//set up the planet pointer in Europa model
		Europa::Planet=Sphere;

		Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);

/*
    Sphere->SurfaceElementDesorptionFluxUP=new double*[PIC::nTotalSpecies];
    Sphere->SurfaceElementAdsorptionFluxDOWN=new double*[PIC::nTotalSpecies];
    Sphere->SurfaceElementPopulation=new double*[PIC::nTotalSpecies];

    Sphere->SurfaceElementDesorptionFluxUP[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
    Sphere->SurfaceElementAdsorptionFluxDOWN[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
    Sphere->SurfaceElementPopulation[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
      Sphere->SurfaceElementDesorptionFluxUP[spec]=Sphere->SurfaceElementDesorptionFluxUP[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
      Sphere->SurfaceElementAdsorptionFluxDOWN[spec]=Sphere->SurfaceElementAdsorptionFluxDOWN[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
      Sphere->SurfaceElementPopulation[spec]=Sphere->SurfaceElementPopulation[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
    }

    Sphere->SurfaceElementExternalNormal=new cInternalSphericalData_UserDefined::cSurfaceElementExternalNormal[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
    Sphere->SurfaceElementArea=new double[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];


    Sphere->ElementSourceRate=new cInternalSphericalData_UserDefined::cElementSourceRate*[PIC::nTotalSpecies];
    Sphere->ElementSourceRate[0]=new cInternalSphericalData_UserDefined::cElementSourceRate[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
      Sphere->ElementSourceRate[spec]=Sphere->ElementSourceRate[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
    }


    Sphere->SolarWindSurfaceFlux=new double[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    for (int el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        Sphere->SurfaceElementDesorptionFluxUP[spec][el]=0.0;
        Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el]=0.0;
        Sphere->SurfaceElementPopulation[spec][el]=0.0;
      }

      Sphere->SolarWindSurfaceFlux[el]=-1.0;

      Sphere->SurfaceElementArea[el]=Sphere->GetSurfaceElementArea(el);
      Sphere->GetSurfaceElementNormal((Sphere->SurfaceElementExternalNormal+el)->norm,el);
    }

    //allocate buffers for sampling surface sodium source rates and sodikum surface content
    int offsetSpecie,offsetElement,s,el,i;

    Sphere->SampleSpeciesSurfaceSourceRate=new double** [PIC::nTotalSpecies];
    Sphere->SampleSpeciesSurfaceSourceRate[0]=new double *[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
    Sphere->SampleSpeciesSurfaceSourceRate[0][0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];


    for (offsetSpecie=0,s=0,offsetElement=0;s<PIC::nTotalSpecies;s++) {
      Sphere->SampleSpeciesSurfaceSourceRate[s]=Sphere->SampleSpeciesSurfaceSourceRate[0]+offsetSpecie;
      offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        Sphere->SampleSpeciesSurfaceSourceRate[s][el]=Sphere->SampleSpeciesSurfaceSourceRate[0][0]+offsetElement;
        offsetElement+=_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;

        for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) Sphere->SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;
      }
    }

//    exit(__LINE__,__FILE__,"the above lines is commented");


    Sphere->SampleSpeciesSurfaceAreaDensity=new double* [PIC::nTotalSpecies];
    Sphere->SampleSpeciesSurfaceAreaDensity[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
      Sphere->SampleSpeciesSurfaceAreaDensity[s]=Sphere->SampleSpeciesSurfaceAreaDensity[0]+offsetSpecie;
      offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        Sphere->SampleSpeciesSurfaceAreaDensity[s][el]=0.0;
      }
    }

    Sphere->SampleSpeciesSurfaceReturnFlux=new double* [PIC::nTotalSpecies];
    Sphere->SampleSpeciesSurfaceReturnFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    Sphere->SampleSpeciesSurfaceInjectionFlux=new double* [PIC::nTotalSpecies];
    Sphere->SampleSpeciesSurfaceInjectionFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    Sphere->SampleReturnFluxBulkSpeed=new double* [PIC::nTotalSpecies];
    Sphere->SampleReturnFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    Sphere->SampleInjectedFluxBulkSpeed=new double* [PIC::nTotalSpecies];
    Sphere->SampleInjectedFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
      Sphere->SampleSpeciesSurfaceReturnFlux[s]=Sphere->SampleSpeciesSurfaceReturnFlux[0]+offsetSpecie;
      Sphere->SampleSpeciesSurfaceInjectionFlux[s]=Sphere->SampleSpeciesSurfaceInjectionFlux[0]+offsetSpecie;

      Sphere->SampleReturnFluxBulkSpeed[s]=Sphere->SampleReturnFluxBulkSpeed[0]+offsetSpecie;
      Sphere->SampleInjectedFluxBulkSpeed[s]=Sphere->SampleInjectedFluxBulkSpeed[0]+offsetSpecie;

      offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        Sphere->SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
        Sphere->SampleSpeciesSurfaceInjectionFlux[s][el]=0.0;

        Sphere->SampleReturnFluxBulkSpeed[s][el]=0.0;
        Sphere->SampleInjectedFluxBulkSpeed[s][el]=0.0;
      }
    }
*/

	}

	//=============   DEBUG ==============

	/*
  double xProbeMin[3]={-2172234,-1.236913E-10,2000742};
  double dx=-2172234-(-2229398),v=1.0,t;
  double xProbeMax[3];
  int IntersectionStatus;

  for (int i=0;i<3;i++) {
    xProbeMax[i]=xProbeMin[i]+dx;
    v*=dx;
  }

  t=Sphere->GetRemainedBlockVolume(xProbeMin,xProbeMax,0.0,&IntersectionStatus);
	 */

	//============= END DEBUG =============

	//init the solver
	PIC::Mesh::initCellSamplingDataBuffer();

	//init the mesh
	cout << "Init the mesh" << endl;

	int maxBlockCellsnumber,minBlockCellsnumber,idim;

	maxBlockCellsnumber=_BLOCK_CELLS_X_;
	if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
	if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);

	minBlockCellsnumber=_BLOCK_CELLS_X_;
	if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
	if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);

	double DomainLength[3],DomainCenterOffset[3],xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

	if (maxBlockCellsnumber==minBlockCellsnumber) {
		for (idim=0;idim<DIM;idim++) {
			DomainLength[idim]=2.0*xMaxDomain*rSphere;
			DomainCenterOffset[idim]=-xMaxDomain*rSphere;
		}
	}
	else {
		if (maxBlockCellsnumber!=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__);
		if (minBlockCellsnumber!=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__);
		if (minBlockCellsnumber!=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__);

		DomainLength[0]=xMaxDomain*rSphere*(1.0+double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_);
		DomainLength[1]=DomainLength[0]*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
		DomainLength[2]=DomainLength[0]*double(_BLOCK_CELLS_Z_)/_BLOCK_CELLS_X_;

		if (DomainLength[1]<2.01*yMaxDomain*_RADIUS_(_TARGET_)) {
			double r;

			printf("Size of the domain is smaller that the radius of the body: the size of the domain is readjusted\n");
			r=2.01*yMaxDomain*_RADIUS_(_TARGET_)/DomainLength[1];

			for (idim=0;idim<DIM;idim++) DomainLength[idim]*=r;
		}

		DomainCenterOffset[0]=-yMaxDomain*rSphere;////-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
		DomainCenterOffset[1]=-DomainLength[1]/2.0;
		DomainCenterOffset[2]=-DomainLength[2]/2.0;
	}

	/*
  for (idim=0;idim<DIM;idim++) {
    xmin[idim]=DomainCenterOffset[idim];
    xmax[idim]=DomainLength[idim]+DomainCenterOffset[idim];
  }
	 */

	for (idim=0;idim<DIM;idim++) {
		xmax[idim]=-DomainCenterOffset[idim];
		xmin[idim]=-(DomainLength[idim]+DomainCenterOffset[idim]);
	}

	//  double xmin[3]={-xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};
	//  double xmax[3]={xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};


	//generate only the tree
	PIC::Mesh::mesh.AllowBlockAllocation=false;
	PIC::Mesh::mesh.init(xmin,xmax,localResolution);
	PIC::Mesh::mesh.memoryAllocationReport();


	//  VT_ON();
	//  VT_USER_START("name");
	//  VT_ON();

	//  {
	//    VT_TRACER("name");

	if (PIC::Mesh::mesh.ThisThread==0) {
		PIC::Mesh::mesh.buildMesh();
		PIC::Mesh::mesh.saveMeshFile("mesh.msh");
		MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
	}
	else {
		MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
		PIC::Mesh::mesh.readMeshFile("mesh.msh");
	}

	//  }
	//  VT_USER_END("name");
	//  MPI_Finalize();
	//  return 1;

	cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

	PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");

	PIC::Mesh::mesh.memoryAllocationReport();
	PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
	PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

	PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
	PIC::Mesh::mesh.CreateNewParallelDistributionLists();

	//initialize the blocks
	PIC::Mesh::mesh.AllowBlockAllocation=true;
	PIC::Mesh::mesh.AllocateTreeBlocks();

	PIC::Mesh::mesh.memoryAllocationReport();
	PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
	PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

	//init the volume of the cells'
	PIC::Mesh::mesh.InitCellMeasure();





	//init the PIC solver
	PIC::Init_AfterParser ();
	PIC::Mover::Init();

	//Exosphere::Init_AfterParser();

	//	PIC::Mover::TotalParticleAcceleration=TotalParticleAcceleration;

	//	for (int s=0;s<PIC::nTotalSpecies;s++) {
	//	  PIC::Mover::MoveParticleTimeStep[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder; ///UniformWeight_UniformTimeStep_SecondOrder;
	//	  PIC::Mover::MoveParticleBoundaryInjection[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder;
	//	}


	//set up the time step
	PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
	PIC::ParticleWeightTimeStep::initTimeStep();

	//set up the particle weight
	PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjectionRate;
	PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_OPLUS_HIGH_SPEC_);
	PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_OPLUS_THERMAL_SPEC_);

	//copy the weight and time step from Na neutra to Na ions
	PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O2_SPEC_,_OPLUS_THERMAL_SPEC_,1.0e4);
//	PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(_O2_SPEC_,_OPLUS_THERMAL_SPEC_,1.0e4);

  PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O2PLUS_SPEC_,_O2_SPEC_,1.0E10*1.0E-7);
  PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(_O2PLUS_SPEC_,_OPLUS_THERMAL_SPEC_,1.0);

	//set photolytic reactions
	//PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(sodiumPhotoionizationReactionProcessor,_O2_SPEC_);
	//PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(sodiumPhotoionizationLifeTime,_O2_SPEC_);


	PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
	PIC::Mesh::mesh.outputMeshDataTECPLOT("mesh.data.dat",0);

	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
	if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;




	//output final data
	//  PIC::Mesh::mesh.outputMeshDataTECPLOT("final.data.dat",0);

	//create the list of mesh nodes where the injection boundary conditions are applied
	PIC::BC::BlockInjectionBCindicatior=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxParticleInjectionIndicator;
	PIC::BC::userDefinedBoundingBlockInjectionFunction=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjection;
	PIC::BC::InitBoundingBoxInjectionBlockList();




	//init the particle buffer
	PIC::ParticleBuffer::Init(10000000);
	//  double TimeCounter=time(NULL);
	int LastDataOutputFileNumber=-1;


	//init the sampling of the particls' distribution functions
	const int nSamplePoints=3;
	double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

	PIC::DistributionFunctionSample::vMin=-40.0E3;
	PIC::DistributionFunctionSample::vMax=40.0E3;
	PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

	PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);


#if _GALILEO_EPD_SAMPLING_ == _EUROPA_MODE_ON_
  //init sampling points along the s/c trajectory
  const int nFlybySamplePasses=1;
  const double FlybySamplingInterval=120.0*60.0,FlybySamplingIntervalStep=120.0; //in seconds
  const int nSampleSteps=(int)(FlybySamplingInterval/FlybySamplingIntervalStep);

  const char *FlybySamplePassesUTC[nFlybySamplePasses]={"1996-12-19T06:00:00"};

  SpiceDouble et,lt;
  SpiceDouble state[6];
  int nFlybyPass,n;

  /* FIPS POINTING
      INS-236720_FOV_FRAME       = 'MSGR_EPPS_FIPS'
      INS-236720_FOV_SHAPE       = 'CIRCLE'
      INS-236720_BORESIGHT       = ( 0.0, 0.0, 1.0 )
   */

  SpiceDouble pointing[3],bsight[3],bsight_INIT[3]={0.0,0.0,1.0};
  SpiceDouble rotate[3][3];

  const SpiceInt lenout = 35;
  SpiceChar utcstr[lenout+2];


  int nFluxSamplePoint=0;

  int nTotalFluxSamplePoints=nFlybySamplePasses*nSampleSteps;
  double FluxSampleLocations[nTotalFluxSamplePoints][3];
  double FluxSampleDirections[nTotalFluxSamplePoints][3];
  double Dist;


  for (nFlybyPass=0;nFlybyPass<nFlybySamplePasses;nFlybyPass++) {
    utc2et_c(FlybySamplePassesUTC[nFlybyPass],&et);

    if (PIC::ThisThread==0) {
      cout << "S/C Flyby Sampling: Pass=" << nFlybyPass << ":" << endl;
      //???      cout << "Flux Sample Point\tUTS\t\t\t x[km]\t\ty[km]\t\tz[km]\t\t\t lx\t\tly\t\tlz\t" << endl;
      cout << "Flux Sample Point\tUTS\t\t\tx[km]\t\ty[km]\t\tz[km]\tr[rTarget]\t" << endl;
    }

    for (n=0;n<nSampleSteps;n++) {
      //position of the s/c
      spkezr_c("GALILEO ORBITER",et,"GALL_EPHIOD","NONE","EUROPA",state,&lt);


      //get the pointing vector in the 'MSO' frame
      //memcpy(bsight,bsight_INIT,3*sizeof(double));

      //???      pxform_c ("MSGR_EPPS_FIPS","GALL_EPHIOD",et,rotate);
      //???      mxv_c(rotate,bsight,pointing);

      //print the pointing information
      if (PIC::ThisThread==0) {
        et2utc_c(et,"C",0,lenout,utcstr);
        printf("%i\t\t\t%s\t",nFluxSamplePoint,utcstr);
        for (idim=0;idim<3;idim++) printf("%e\t",state[idim]);
        Dist = 0.0;
        for (idim=0;idim<3;idim++) Dist+=pow(state[idim],2);
        printf("%.2f\t",sqrt(Dist)/_RADIUS_(_TARGET_)*1E3);
        cout << "\t";

        //???        for (idim=0;idim<3;idim++) printf("%e\t",pointing[idim]);
        cout << endl;
      }

      //save the sample pointing information
      for (idim=0;idim<3;idim++) {
        FluxSampleLocations[nFluxSamplePoint][idim]=state[idim]*1.0E3;
        //???  FluxSampleDirections[nFluxSamplePoint][idim]=pointing[idim];
      }

      //increment the flyby time
      et+=FlybySamplingIntervalStep;
      ++nFluxSamplePoint;
    }

    cout << endl;
  }

 // PIC::ParticleFluxDistributionSample::Init(FluxSampleLocations,FluxSampleDirections,30.0/180.0*Pi,nTotalFluxSamplePoints);

#elif _GALILEO_EPD_SAMPLING_ == _EUROPA_MODE_OFF_
  printf("No Galileo EPD sampling\n");
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif


	//init ICES

#if _ICES_CREATE_COORDINATE_LIST_ == _PIC_MODE_ON_
  const char IcesLocationPath[]="";//"/Users/dborovik/MyICES/ICES";
  const char IcesModelCase[]="";//"Europa09";


	PIC::CPLR::ICES::createCellCenterCoordinateList();
	PIC::CPLR::ICES::SetLocationICES(IcesLocationPath);
	PIC::CPLR::ICES::retriveSWMFdata(IcesModelCase);  ////("EUROPA_RESTART_n070001");
#endif


//#ifdef _ICES_LOAD_DATA_

#if _PIC_COUPLER_MODE_ ==	_PIC_COUPLER_MODE__ICES_
	PIC::CPLR::ICES::readSWMFdata(1.0);
	//  PIC::Mesh::mesh.outputMeshDataTECPLOT("ices.data.dat",0);


	//output the solar wind ion flux at the palnet's surface
	PIC::CPLR::ICES::PrintSphereSurfaceIonFlux("SurfaceIonFlux.dat",1.05*_RADIUS_(_TARGET_));
	PIC::CPLR::ICES::EvaluateSurfaceIonFlux(1.05);
	
	PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test.dat",0);
	
	Exosphere::SourceProcesses::Init();
	PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test-1.dat",0);


	//create the map of the solar wind flux
	int el;

	for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
		int i,j,k;
		long int nd;
		double FaceCenterPoint[3],PlasmaVelocity[3],PlasmaNumberDensity,FaceElementNormal[3],c;
		cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
		PIC::Mesh::cDataCenterNode *CenterNode;
		char *offset;

		Europa::Planet->GetSurfaceElementMiddlePoint(FaceCenterPoint,el);
		Europa::Planet->GetSurfaceElementNormal(FaceElementNormal,el);

		if (FaceElementNormal[0]<0.0) continue;

		node=PIC::Mesh::mesh.findTreeNode(FaceCenterPoint);

		if ((nd=PIC::Mesh::mesh.fingCellIndex(FaceCenterPoint,i,j,k,node,false))==-1) {
			exit(__LINE__,__FILE__,"Error: the cell is not found");
		}

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
		if ((node->Thread==PIC::ThisThread)&&(node->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif

		if (node->Thread==PIC::ThisThread) {
			CenterNode=node->block->GetCenterNode(nd);
			offset=CenterNode->GetAssociatedDataBufferPointer();

			if (*((int*)(offset+PIC::CPLR::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
				memcpy(PlasmaVelocity,offset+PIC::CPLR::ICES::PlasmaBulkVelocityOffset,3*sizeof(double));
				memcpy(&PlasmaNumberDensity,offset+PIC::CPLR::ICES::PlasmaNumberDensityOffset,sizeof(double));
			}
			else {
//				exit(__LINE__,__FILE__,"Error: PIC_ICES__STATUS is not _PIC_ICES__STATUS_OK_");

			  PlasmaVelocity[0]=0.0,PlasmaVelocity[1]=0.0,PlasmaVelocity[2]=0.0;
			  PlasmaNumberDensity=0.0;
			}

			c=-(PlasmaVelocity[0]*FaceElementNormal[0]+PlasmaVelocity[1]*FaceElementNormal[1]+PlasmaVelocity[2]*FaceElementNormal[2]);
			if (c<0.0) c=0.0;

			Europa::Planet->SolarWindSurfaceFlux[el]=c*PlasmaNumberDensity;
		}

		MPI_Bcast(Europa::Planet->SolarWindSurfaceFlux+el,1,MPI_DOUBLE,node->Thread,MPI_GLOBAL_COMMUNICATOR);
	}


#endif

	//prepopulate the solar wind protons
	//  prePopulateSWprotons(PIC::Mesh::mesh.rootTree);


	/*
  //check the volume of local cells
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  while (node!=NULL) {
    int i,j,k,di,dj,dk,idim;
    long int LocalCellNumber;
    double r,rmin=0.0,rmax=0.0,rprobe[3]={0.0,0.0,0.0};
    PIC::Mesh::cDataCenterNode *cell;

    if (node->Temp_ID==1862) {
      cout << __FILE__ << "@" << __LINE__ << endl;
    }

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            cell=node->block->GetCenterNode(LocalCellNumber);
            rmin=-1,rmax=-1;

            if (cell->Measure<=0.0) {
              for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
                node->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

                for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
                r=sqrt(r);

                if ((rmin<0.0)||(rmin>r)) rmin=r;
                if ((rmax<0.0)||(rmax<r)) rmax=r;
              }

              if ((rmin<rSphere)&&(rmax>rSphere)) {
                cout << "Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << r << endl;
              }

              PIC::Mesh::mesh.InitCellMeasure(node);

            }


          }
       }
    }

    node=node->nextNodeThisThread;
  }
	 */

	//  VT_TRACER("main");

}

	//time step

void amps_time_step () {
  //  cout <<'============= EUROPA PARAMS ==================' << endl;
  // cout << Europa::xEuropa[0] << endl;
  //  cout << Europa::xEuropa[1] << endl;
  //  cout << Europa::xEuropa[2] << endl;
  //  cout << Europa::vEuropa[0] << endl;
  //  cout << Europa::vEuropa[1] << endl;
  //  cout << Europa::vEuropa[2] << endl;
  //  cout << Europa::OrbitalMotion::et << endl;
  //  cout << Europa::OrbitalMotion::TAA << endl;
  //  cout << Europa::OrbitalMotion::CoordinateFrameRotationRate << endl;
  //  cout << Europa::OrbitalMotion::PlanetAxisToSunRotationAngle << endl;
  //  cout << Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0] << endl;
  //  cout << Europa::OrbitalMotion::SunDirection_IAU_EUROPA[1] << endl;
  //  cout << Europa::OrbitalMotion::SunDirection_IAU_EUROPA[2] << endl;
  //  int iTmp,jTmp;
  //  for(iTmp=0;iTmp<6;iTmp++)
  //    for(jTmp=0;jTmp<6;jTmp++)
  //      cout << Europa::OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix[iTmp][jTmp] << endl;
  //  for(iTmp=0;iTmp<6;iTmp++)
  //    for(jTmp=0;jTmp<6;jTmp++)
  //      cout << Europa::OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[iTmp][jTmp] << endl;
  //  for(iTmp=0;iTmp<6;iTmp++)
  //    for(jTmp=0;jTmp<6;jTmp++)
  //      cout << Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[iTmp][jTmp] << endl;
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int idim;


		//determine the parameters of the orbital motion of Europa
		SpiceDouble StateBegin[6],StateEnd[6],lt;
		double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;
		//    int idim;

		spkezr_c("Europa",Europa::OrbitalMotion::et,"GALL_EPHIOD","none","Jupiter",StateBegin,&lt);
		Europa::OrbitalMotion::et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_];
		spkezr_c("Europa",Europa::OrbitalMotion::et,"GALL_EPHIOD","none","Jupiter",StateEnd,&lt);

		for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
			StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
			StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;

			rBegin+=pow(StateBegin[idim],2);
			rEnd+=pow(StateEnd[idim],2);

			Europa::xEuropa[idim]=StateBegin[idim];
			Europa::vEuropa[idim]=StateBegin[3+idim];
		}

		rBegin=sqrt(rBegin);
		rEnd=sqrt(rEnd);

		for (idim=0;idim<3;idim++) {
			lBegin[idim]=StateBegin[idim]/rBegin;
			lEnd[idim]=StateEnd[idim]/rEnd;

			c0+=StateBegin[3+idim]*lBegin[idim];
			c1+=StateEnd[3+idim]*lEnd[idim];
		}

		Europa::xEuropaRadial=rBegin;
		Europa::vEuropaRadial=0.5*(c0+c1);

		//calculate TAA
		{
			double EccentricityVector[3];
			double Speed2,a,c,absEccentricity,cosTAA;

			const double GravitationalParameter=GravityConstant*_MASS_(_SUN_);

			Speed2=Europa::vEuropa[0]*Europa::vEuropa[0]+Europa::vEuropa[1]*Europa::vEuropa[1]+Europa::vEuropa[2]*Europa::vEuropa[2];
			c=Europa::xEuropa[0]*Europa::vEuropa[0]+Europa::xEuropa[1]*Europa::vEuropa[1]+Europa::xEuropa[2]*Europa::vEuropa[2];

			for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
				EccentricityVector[idim]=Speed2/GravitationalParameter*Europa::xEuropa[idim] - c/GravitationalParameter*Europa::vEuropa[idim] - Europa::xEuropa[idim]/rBegin;
				absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
				a+=EccentricityVector[idim]*Europa::xEuropa[idim];
			}

			absEccentricity=sqrt(absEccentricity);
			cosTAA=a/(absEccentricity*rBegin);

			Europa::OrbitalMotion::TAA=acos(cosTAA);
			if (c<0.0) Europa::OrbitalMotion::TAA=2.0*Pi-Europa::OrbitalMotion::TAA;
		}


		//the index varying from 0 to 1!!!!!! because the velocity must be perpendicular to z-axis, which is the axis of rotation
		for (idim=0;idim<2;idim++) {
			vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
			vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
		}

		vTangentialBegin=sqrt(vTangentialBegin);
		vTangentialEnd=sqrt(vTangentialEnd);

		Europa::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);


		//determine direction to the Sun and rotation angle in the coordinate frame related to Europa
		SpiceDouble state[6],l=0.0;

		spkezr_c("SUN",Europa::OrbitalMotion::et-0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_],"IAU_EUROPA","none","EUROPA",state,&lt);

		for (idim=0;idim<3;idim++) l+=pow(state[idim],2);

		for (l=sqrt(l),idim=0;idim<3;idim++) {
			Europa::OrbitalMotion::SunDirection_IAU_EUROPA[idim]=state[idim]/l;
		}

		Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=acos(Europa::OrbitalMotion::SunDirection_IAU_EUROPA[1]);
		if (Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0]<0.0) Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=2.0*Pi-Europa::OrbitalMotion::PlanetAxisToSunRotationAngle;


		//matrixes for transformation GALL_EPHIOD->IAU and IAU->GALL_EPHIOD coordinate frames
		sxform_c("GALL_EPHIOD","IAU_EUROPA",Europa::OrbitalMotion::et-0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_],Europa::OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix);
		sxform_c("IAU_EUROPA","GALL_EPHIOD",Europa::OrbitalMotion::et-0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_],Europa::OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix);

		sxform_c("IAU_EUROPA","GALL_ESOM",Europa::OrbitalMotion::et-0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_],Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix);

#else
		double xEuropaTest[3] = {-7.27596e-09, - 6.76112e+08, - 2.76134e+06}; 
		std::copy(&xEuropaTest[0], &xEuropaTest[0]+3, &Europa::xEuropa[0]);
		
		double vEuropaTest[3] = { 0,            77.5556,       92.079};
		std::copy(&vEuropaTest[0], &vEuropaTest[0]+3, &Europa::vEuropa[0]);

		double SunDirTest[3]  = {-0.364833, 0.931007, -0.0111228};
		std::copy(&SunDirTest[0], &SunDirTest[0]+3, &Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0]);
		
		Europa::OrbitalMotion::et  = -9.57527e+07;
		Europa::OrbitalMotion::TAA =  3.14159;
		Europa::OrbitalMotion::CoordinateFrameRotationRate = 5.5426e-10;
		Europa::OrbitalMotion::PlanetAxisToSunRotationAngle= 5.90955;

		SpiceDouble TransMatrix1[6][6] = {
		  {-0.0307892,   0.999518,    0.00388983,  0,         0,         0}, 
		  {-0.999503,   -0.0307619,  -0.00689729,  0,         0,         0},
		  {-0.0067743,  -0.00410026,  0.999969,    0,         0,         0},
		  {-3.05046e-07,-8.84681e-09,-1.41289e-07,-0.0307892, 0.999518,  0.00388983},
		  { 9.95762e-09,-3.05672e-07,-7.9686e-08, -0.999503, -0.0307619,-0.00689729},
		  {-8.27439e-08, 1.367e-07,  -2.56459e-14,-0.0067743,-0.00410026,0.999969}};
		std::copy(&TransMatrix1[0][0],&TransMatrix1[0][0]+36, &Europa::OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix[0][0]);

		SpiceDouble TransMatrix2[6][6] = {
		  {-0.0307892,  -0.999503,   -0.0067743,   0,          0,         0},
		  { 0.999518,   -0.0307619,  -0.00410026,  0,          0,         0},
		  { 0.00388983, -0.00689729,  0.999969,    0,          0,         0},
		  {-3.05046e-07, 9.95762e-09,-8.27439e-08,-0.0307892, -0.999503, -0.0067743},
		  {-8.84681e-09,-3.05672e-07, 1.367e-07,   0.999518,  -0.0307619,-0.00410026},
		  {-1.41289e-07,-7.9686e-08, -2.56459e-14, 0.00388983,-0.00689729,0.999969}};
		std::copy(&TransMatrix2[0][0],&TransMatrix2[0][0]+36,&Europa::OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]);

		SpiceDouble TransMatrix3[6][6] = {
		  {-0.364833,   0.931007,  -0.0111228,   0,         0,        0},
		  {-0.931064,  -0.364856,  -2.77556e-17, 0,         0,        0},
		  {-0.00405822, 0.0103561,  0.999938,    0,         0,        0},
		  { 1.90558e-05,7.46742e-06,1.08988e-09,-0.364833,  0.931007,-0.0111228},
		  {-7.46787e-06,1.9057e-05,-1.07033e-24,-0.931064, -0.364856,-2.77556e-17},
		  { 2.12366e-07,8.2049e-08, 1.21233e-11,-0.00405822,0.0103561,0.999938}};
		std::copy(&TransMatrix3[0][0],&TransMatrix3[0][0]+36,&Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[0][0]);

#endif

		//make the time advance
		static int LastDataOutputFileNumber=0;


		PIC::TimeStep();

		//     cell=(node->block!=NULL) ? node->block->GetCenterNode(nd) : NULL;


		if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
			PIC::RequiredSampleLength*=2;
			if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;


			LastDataOutputFileNumber=PIC::DataOutputFileNumber;
			if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
		}



}
