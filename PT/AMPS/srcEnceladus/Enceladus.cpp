/*
 * Enceladus.cpp
 *
 *  Created on: Jan 17, 2012
 *      Author: vtenishe
 */

//$Id$

#include "Enceladus.h"

//parameters of the injecting plumes
double Enceladus::sphericalSourceProductionRateConst;
double Enceladus::sphericalSourceExpansionVelocityConst;
double Enceladus::sphericalSource_Const0;
double Enceladus::sphericalSource_Const1;
Cplume Enceladus::Plume[nTotalPlumes];

/*---------------------------------    SURFACE BOUNDARY CONDITIONS     ------------------------------------*/
//production rate of the grains
double Enceladus::sphereInjectionRate(int spec,void *SphereDataPointer) {
   static bool InitFlag=false;
   static double MeanDustGrainMass=0.0;

   if (InitFlag==false) {
     InitFlag=true;

     //get mean mass of a dust grain
     const int nTotalTest=10000;
     double r,c,cTotal=0.0,m=0.0;


     for (int ntest=0;ntest<nTotalTest;ntest++) {
       SizeDistribution::GenerateGrainRandomRadius(r,c);
       m+=pow(r,3.0)*c;
       cTotal+=c;
     }

     MeanDustGrainMass=m*4.0/3.0*Pi*MeanDustDensity/cTotal;
   }


   if ((spec<_DUST_SPEC_)||(spec>=_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) exit(__LINE__,__FILE__,"Error: unknown speceis");

   return TotalMassDustProductionRate/MeanDustGrainMass;
}

int Enceladus::ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
  /*double r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  double c=v[0]*x[0]/r+v[1]*x[1]/r+v[2]*x[2]/r;

  v[0]-=2.0*c*x[0]/r;
  v[1]-=2.0*c*x[1]/r;
  v[2]-=2.0*c*x[2]/r;

  return _PARTICLE_REJECTED_ON_THE_FACE_;*/
  

  PIC::ParticleBuffer::byte *ParticleData;
  cInternalSphericalData *Sphere;
  double GrainMass,GrainWeightCorrection,ParticleWeight,LocalTimeStep;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  Sphere=(cInternalSphericalData*)SphereDataPointer;

  //sample the injection flux
  //sample the particle data
  double *SampleData;
  long int nSurfaceElement,nZenithElement,nAzimuthalElement;

  Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
  nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
  SampleData=Sphere->SamplingBuffer+PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(_DUST_SPEC_,nSurfaceElement);

  GrainMass=GetGrainMass(ParticleData);
  GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

  LocalTimeStep=((cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter)->block->GetLocalTimeStep(spec);

  SampleData[PIC::BC::InternalBoundary::Sphere::sampledFluxDownRelativeOffset]+=GrainMass*ParticleWeight*GrainWeightCorrection/LocalTimeStep/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);

  //  PIC::ParticleBuffer::DeleteParticle(ptr);
  return _PARTICLE_DELETED_ON_THE_FACE_;
}





/*---------------------------------  END SURFACE BOUNDARY CONDITIONS     ------------------------------------*/



/*---------------------------------   INJECT DUST GRAINS  ---------------------------------------------*/
bool Enceladus::GenerateInitialGrainParameters(double *x,double *v,double& GrainRadius, double& GrainWeightCorrection,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode) {
  int idim;
  int nInjectionPlume;

  static bool InitFlag=false;
  static double plumeCoordinateSystem[nTotalPlumes][3][3];
  static double maxPlumeProductionRate=0.0;
  static double plumeProductionRateTable[nTotalPlumes];

  //characteristic size of the exit vent
  const double characteristicVentSize=100.0;
  double Speed=2.0;

  if (InitFlag==false) {
    double l;

    InitFlag=true;

    //determine the clocal coordinate system related to the position of the origin of each plume
    for (int nplume=0;nplume<nTotalPlumes;nplume++) {
      for (idim=0;idim<3;idim++) plumeCoordinateSystem[nplume][2][idim]=Plume[nplume].plumeDirection[idim];

      if (fabs(Plume[nplume].plumeDirection[0])>1.0E-5) {
        plumeCoordinateSystem[nplume][0][0]=Plume[nplume].plumeDirection[1];
        plumeCoordinateSystem[nplume][0][1]=-Plume[nplume].plumeDirection[0];
        plumeCoordinateSystem[nplume][0][2]=0.0;
      }
      else {
        plumeCoordinateSystem[nplume][0][1]=Plume[nplume].plumeDirection[2];
        plumeCoordinateSystem[nplume][0][2]=-Plume[nplume].plumeDirection[1];
        plumeCoordinateSystem[nplume][0][0]=0.0;
      }

      for (l=0.0,idim=0;idim<3;idim++) l+=pow(plumeCoordinateSystem[nplume][0][idim],2.0);
      for (l=sqrt(l),idim=0;idim<3;idim++) plumeCoordinateSystem[nplume][0][idim]/=l;

      plumeCoordinateSystem[nplume][1][0]=plumeCoordinateSystem[nplume][0][1]*plumeCoordinateSystem[nplume][2][2]-plumeCoordinateSystem[nplume][2][1]*plumeCoordinateSystem[nplume][0][2];
      plumeCoordinateSystem[nplume][1][1]=-(plumeCoordinateSystem[nplume][0][0]*plumeCoordinateSystem[nplume][2][2]-plumeCoordinateSystem[nplume][2][0]*plumeCoordinateSystem[nplume][0][2]);
      plumeCoordinateSystem[nplume][1][2]=plumeCoordinateSystem[nplume][0][0]*plumeCoordinateSystem[nplume][2][1]-plumeCoordinateSystem[nplume][2][0]*plumeCoordinateSystem[nplume][0][1];

      //get the maximum plume production rate
      if (maxPlumeProductionRate<(plumeProductionRateTable[nplume]=Plume[nplume].SourceRate)) maxPlumeProductionRate=Plume[nplume].SourceRate;
    }

  }

  //determine the vent for injection of the dust grain
  double p,r,phi,e0,e1,x_LOCAL[3],v_LOCAL[3],plumeProductionRateTable_LOCAL[nTotalPlumes];

  memcpy(plumeProductionRateTable_LOCAL,plumeProductionRateTable,nTotalPlumes*sizeof(double));

  do {
    nInjectionPlume=(int)(rnd()*nTotalPlumes);
    p=plumeProductionRateTable_LOCAL[nInjectionPlume]/maxPlumeProductionRate;
  }
  while (p<rnd());

  //copy the coordinate vectors into a local variable
  double plumeCoordinateSystem_LOCAL[9];

  memcpy(plumeCoordinateSystem_LOCAL,plumeCoordinateSystem[nInjectionPlume],9*sizeof(double));

  //shift initial position of the injected particles within the vent
  memcpy(x_LOCAL,Plume[nInjectionPlume].plumePosition,3*sizeof(double));

  phi=PiTimes2*rnd();
  r=characteristicVentSize*sqrt(rnd());
  e0=r*sin(phi),e1=r*cos(phi);

  for (idim=0;idim<3;idim++) x_LOCAL[idim]+=e0*plumeCoordinateSystem_LOCAL[0*3+idim]+e1*plumeCoordinateSystem_LOCAL[1*3+idim];




  double dv,dv2,VentBulkFlowNumberDensity;
  VentBulkFlowNumberDensity=plumeProductionRateTable_LOCAL[nInjectionPlume]/(Pi*pow(characteristicVentSize,2))/Plume[nInjectionPlume].BulkVelocity;

  do {
    //get the radius of the grain
    SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);

    //calculate the terminal velocity of the dust grain
    dv2=4.0/3.0*GrainRadius*MeanDustDensity*GravityConstant*_MASS_(_TARGET_)/pow(_RADIUS_(_TARGET_),2.0)/(GrainDragCoefficient/2.0*_MASS_(_H2O_)*VentBulkFlowNumberDensity);
  }
  while (dv2<0.0);

  dv=sqrt(dv2);
  Speed=Plume[nInjectionPlume].BulkVelocity-dv;

  //generate the particle velocity
  for (idim=0;idim<3;idim++) v_LOCAL[idim]=Speed*plumeCoordinateSystem_LOCAL[2*3+idim];


  //check if the particle position is outside of Enceladus
  register double c=(_RADIUS_(_TARGET_)+PIC::Mesh::mesh.EPS)/sqrt(x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2]);
  for (idim=0;idim<3;idim++) x_LOCAL[idim]*=c;

  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL,startNode);

  //copy the position and velocity vectors to the global variables
  memcpy(v,v_LOCAL,3*sizeof(double));
  memcpy(x,x_LOCAL,3*sizeof(double));

  return (startNode->Thread==PIC::Mesh::mesh.ThisThread) ? true : false;
}




/*---------------------------------   END INJECT DUST GRAINS  ---------------------------------------------*/


//functions that output the data file
void Enceladus::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Local Atmosphere number density [m^{-3}]\"");
}


void Enceladus::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
}


void Enceladus::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {


  if (pipe->ThisThread==0) {
    double xNode[3],r,BackgroundDensity=0.0;
    int nplume;

    CenterNode->GetX(xNode);
    r=sqrt(pow(xNode[0],2)+pow(xNode[1],2)+pow(xNode[2],2));

    //print the local number density of the atmosphere
    if (r>_RADIUS_(_TARGET_)) {
      BackgroundDensity=sphericalSource_Const0/pow(r,2)+sphericalSource_Const1;
      for (nplume=0;nplume<nTotalPlumes;nplume++) BackgroundDensity+=Plume[nplume].GetDensity(xNode);
    }

    fprintf(fout,"%e  ",BackgroundDensity);

  }
}
