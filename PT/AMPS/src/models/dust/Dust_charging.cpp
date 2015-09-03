
//$Id$

/*
 * Dust_charging.cpp
 *
 *  Created on: Aug 30, 2015
 *      Author: vtenishe
 */

#include "pic.h"
#include "Dust.h"

//update charge of the dust grain
int ElectricallyChargedDust::Charging::UpdateGrainCharge__EQUILIBRIUM_POTENTIAL(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {
  return UpdateGrainCharge(xInit,xFinal,v,spec,ptr,ParticleData,dt,initNode,CHARGE_INTEGRATION_MODE__EQUILIBRIUM_POTENTIAL);
}

int ElectricallyChargedDust::Charging::UpdateGrainCharge__TIME_DEPENDENT(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {
  return UpdateGrainCharge(xInit,xFinal,v,spec,ptr,ParticleData,dt,initNode,CHARGE_INTEGRATION_MODE__TIME_DEPENDENT);
}


//calcualte the grain's currents and current devivatives
void ElectricallyChargedDust::Charging::GetGrainCurrent(
  double GrainRadius, double &GrainElectricCharge,double *GrainVelocity, //radius, the electric charge, and velocity of the grain
  double HeliocentricDistance, //heliocentric distance of the grain in AU
  double ni,double Ti, double *Vi,  //ion number density, temeprature, and velocity
  double ne, double Te,  //electron number density and temeprature
  double &Je, double &dJe, //the electron collection current
  double &Ji, double &dJi, //ion collection current
  double &Jpe, double &dJpe, //photo-electron current
  double &Jse, double &dJse //secondary-electron current
  ) {

  double DustPotential=GetDustGrainPotential(GrainRadius,GrainElectricCharge);

  //ELECRON COLLECTION CURRENT
  #if _DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_
  double J0e=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*ne*sqrt(Kbol*Te/(PiTimes2*ElectronMass));
  double XiElectron=-ElectronCharge*DustPotential/(Kbol*Te);

  if (XiElectron>=0.0) {
    double t2 = 0.1e1 / 0.3141592654e1;
    double t3 = 0.1e1 / VacuumPermittivity;
    double t6 = 0.1e1 / GrainRadius;
    double t7 = 0.1e1 / Kbol;
    double t9 = 0.1e1 / Te;
    double t17 = exp(ElectronCharge * GrainElectricCharge * t2 * t3 * t6 * t7 * t9 / 0.4e1);

    Je=J0e*exp(-XiElectron);
    dJe = J0e * ElectronCharge * t2 * t3 * t6 * t7 * t9 * t17 / 0.4e1;
  }
  else {
    Je=J0e*(1.0-XiElectron);
    dJe=J0e * ElectronCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1;
  }
  #else  //_DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_
  Je=0.0,dJe=0.0;
  #endif //_DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_

  //ION COLLECTION CURRENT
  #if _DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_
  double J0i=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*ni*sqrt(Kbol*Ti/(PiTimes2*ProtonMass));
  double M=sqrt((pow(Vi[0]-GrainVelocity[0],2)+pow(Vi[1]-GrainVelocity[1],2)+pow(Vi[2]-GrainVelocity[2],2))/(2.0*Kbol*Ti/ProtonMass));

  if (DustPotential<=0.0) {
    {
      double t1 = M * M;
      double t15 = sqrt(0.3141592654e1);
      double t17 = erf(M);
      double t21 = exp(-t1);
      Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - ElectronCharge * GrainElectricCharge / 0.3141592654e1 / VacuumPermittivity / Kbol / Ti / GrainRadius / 0.4e1) * t15 * t17 / M + t21) / 0.2e1;
    }

    {
      double t2 = sqrt(0.3141592654e1);
      double t11 = erf(M);
      dJi = -J0i * ElectronCharge / t2 / VacuumPermittivity / Kbol / Ti / GrainRadius * t11 / M / 0.8e1;
    }
  }
  else {
    {
      double t1 = M * M;
      double t2 = ElectronCharge * GrainElectricCharge;
      double t3 = 0.1e1 / 0.3141592654e1;
      double t5 = 0.1e1 / VacuumPermittivity;
      double t6 = 0.1e1 / Kbol;
      double t8 = 0.1e1 / Ti;
      double t9 = 0.1e1 / GrainRadius;
      double t12 = t2 * t3 * t5 * t6 * t8 * t9;
      double t15 = sqrt(0.3141592654e1);
      double t17 = sqrt(t12);
      double t18 = t17 / 0.2e1;
      double t19 = M + t18;
      double t20 = erf(t19);
      double t21 = M - t18;
      double t22 = erf(t21);
      double t24 = 0.1e1 / M;
      double t33 = sqrt(t2 * t3 * t5 * t6 * t8 * t9 * t24);
      double t34 = t33 / 0.2e1;
      double t36 = t21 * t21;
      double t37 = exp(-t36);
      double t40 = t19 * t19;
      double t41 = exp(-t40);
      Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - t12 / 0.4e1) * t15 * (t20 + t22) * t24 + (t34 + 0.1e1) * t37 - (t34 - 0.1e1) * t41) / 0.2e1;
    }

    {
      double t1 = sqrt(0.3141592654e1);
      double t4 = 0.1e1 / VacuumPermittivity;
      double t5 = 0.1e1 / Kbol;
      double t6 = t4 * t5;
      double t8 = 0.1e1 / Ti;
      double t9 = 0.1e1 / GrainRadius;
      double t10 = t8 * t9;
      double t11 = ElectronCharge * GrainElectricCharge;
      double t12 = 0.1e1 / 0.3141592654e1;
      double t14 = t6 * t10;
      double t15 = t11 * t12 * t14;
      double t16 = sqrt(t15);
      double t17 = t16 / 0.2e1;
      double t18 = M + t17;
      double t19 = erf(t18);
      double t20 = M - t17;
      double t21 = erf(t20);
      double t23 = 0.1e1 / M;
      double t28 = M * M;
      double t33 = 0.1e1 / t1 / 0.3141592654e1;
      double t34 = t18 * t18;
      double t35 = exp(-t34);
      double t38 = 0.1e1 / t16 * ElectronCharge;
      double t41 = t20 * t20;
      double t42 = exp(-t41);
      double t49 = t12 * t4;
      double t51 = t5 * t8;
      double t52 = t9 * t23;
      double t55 = sqrt(t11 * t49 * t51 * t52);
      double t58 = 0.1e1 / t55 * ElectronCharge * t49;
      double t63 = t55 / 0.2e1;
      double t66 = t38 * t12;
      dJi = J0i * (-ElectronCharge / t1 * t6 * t10 * (t19 + t21) * t23 / 0.4e1 + (t28 + 0.1e1 / 0.2e1 - t15 / 0.4e1) * t1 * (t33 * t35 * t38 * t14 - t33 * t42 * t38 * t14) * t23 / 0.2e1 + t58 * t51 * t52 * t42 / 0.4e1 + (t63 + 0.1e1) * t20 * t66 * t6 * t10 * t42 / 0.2e1 - t58 * t51 * t52 * t35 / 0.4e1 + (t63 - 0.1e1) * t18 * t66 * t6 * t10 * t35 / 0.2e1) / 0.2e1;

    }
  }
  #else //_DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_
   Ji=0.0,dJi=0.0;
  #endif //_DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_

  //PHOTO-ELECTRON CURRENT
  #if _DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  double J0pe=Pi*GrainRadius*GrainRadius*ElectronCharge*ElectonPhotoEmission::PhotoElectronEfficiency*ElectonPhotoEmission::PhotoElectronEfficiencyMaterialConstant/pow(HeliocentricDistance,2);

  if (DustPotential<=0.0) {
    Jpe=J0pe,dJpe=0.0;
  }
  else {
    {
      double t10 = exp(-GrainElectricCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / ElectonPhotoEmission::PhotoElectronEvergy / 0.4e1);
      Jpe = J0pe * t10;
    }
    {
      double t1 = 0.1e1 / 0.3141592654e1;
      double t3 = 0.1e1 / VacuumPermittivity;
      double t5 = 0.1e1 / GrainRadius;
      double t6 = 0.1e1 / ElectonPhotoEmission::PhotoElectronEvergy;
      double t13 = exp(-GrainElectricCharge * t1 * t6 * t3 * t5 / 0.4e1);
      dJpe = -J0pe * t1 * t3 * t5 * t6 * t13 / 0.4e1;

    }
  }
  #else //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  Jpe=0.0,dJpe=0.0;
  #endif //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_

  //SECONDARY ELECTON CURRENT
  #if _DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  exit(__LINE__,__FILE__,"Error: not implemented");
  #else //_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  Jse=0.0,dJse=0.0;
  #endif //_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
}




//the generic procedure of the grain charge update
int ElectricallyChargedDust::Charging::UpdateGrainCharge(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode,int CHARGE_INTEGRATION_MODE) {
  double PlasmaTemperature,PlasmaNumberDensity;
  double GrainElectricCharge,GrainElectricCharge_NEW,InitGrainCharge;

  PIC::Mesh::cDataCenterNode* cell;
  int i,j,k;
  long int LocalCellNumber;
  double swVel[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *finalNode=PIC::Mesh::mesh.findTreeNode(xFinal,initNode);

  //the procesure is applied only to dust
  if ((spec<_DUST_SPEC_) || (spec>=_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,finalNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell where the particle is located");
  cell=finalNode->block->GetCenterNode(LocalCellNumber);

  //get the grain electric potential
  char localParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double M,GrainRadius,dustPotential;

  memcpy((void*)localParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
  GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)localParticleData);
  GrainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)localParticleData);
  InitGrainCharge=GrainElectricCharge;


  //reserve space for different elecgtron and ion temepratures
  double Ti,Te,Je,dJe,Ji,dJi,Jpe,dJpe,Jse,dJse,pe;


  PIC::CPLR::InitInterpolationStencil(xInit,finalNode);
  PlasmaTemperature=PIC::CPLR::GetBackgroundPlasmaTemperature();
  PIC::CPLR::GetBackgroundPlasmaVelocity(swVel);
  PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
  pe=PIC::CPLR::GetBackgroundElectronPlasmaPressure();


  if (PlasmaNumberDensity<1.0E2) {
    PlasmaTemperature=200.0;
    PlasmaNumberDensity=1.0E2;
    swVel[0]=100.0,swVel[1]=0.0,swVel[2]=0.0;
    Ti=PlasmaTemperature;
    Te=PlasmaTemperature;
  }
  else{
    Ti=PlasmaTemperature;
    Te=pe/(Kbol*PlasmaNumberDensity);
  }


//==========  DEBUG  BEGIN =================
  static long int nFunctionCalls=0;

  nFunctionCalls++;
//==========   DEBUG END ===================

  //Evaluate the initial value of the current and devivative of the current
  double Jtotal=0.0,dJtotal=0.0;
  static double HeliocentricDistance=sqrt(Exosphere::xSun_SO[0]*Exosphere::xSun_SO[0]+
      Exosphere::xSun_SO[1]*Exosphere::xSun_SO[1]+
      Exosphere::xSun_SO[2]*Exosphere::xSun_SO[2])/_AU_;


  //equation that will be solved iteratively: f(q)=0, where
  //f(q)=q-qInit-dt*Jtotal (for time dependent integration)
  //f(q)=Jtotal (for the steady state)

  double ChargeIncrement;
  int Counter=0;

  static const double ChargeIncrementCup=100.0*ElectronCharge;
  static const double ChargeIncrementMaxFractionCup=3.0;
  static const double ChargeIncrementMinFractionCup=1.0/ChargeIncrementMaxFractionCup;


  do {
    GetGrainCurrent(
      GrainRadius,GrainElectricCharge,v, //radius, the electric charge, and velocity of the grain
      HeliocentricDistance, //heliocentric distance of the grain in AU
      PlasmaNumberDensity,Ti,swVel,  //ion number density, temeprature, and velocity
      PlasmaNumberDensity,Te,  //electron number density and temeprature
      Je, dJe, //the electron collection current
      Ji, dJi, //ion collection current
      Jpe, dJpe, //photo-electron current
      Jse, dJse //secondary-electron current
      );

    Jtotal=Je+Ji+Jpe+Jse;
    dJtotal=dJe+dJi+dJpe+dJse;


    switch (CHARGE_INTEGRATION_MODE) {
    case CHARGE_INTEGRATION_MODE__EQUILIBRIUM_POTENTIAL :
      ChargeIncrement=-Jtotal/dJtotal;
      break;
    case CHARGE_INTEGRATION_MODE__TIME_DEPENDENT:
      ChargeIncrement=InitGrainCharge-dt*Jtotal/dJtotal;
      break;
    default:
      exit(__LINE__,__FILE__,"Unknown option");
    }

    //limit variation of the grain's charge
    if (fabs(GrainElectricCharge)<ChargeIncrementCup) {
      double t=ChargeIncrement+GrainElectricCharge;

      if (t<-ChargeIncrementCup) t=-ChargeIncrementCup;
      else if (t>ChargeIncrementCup) t=ChargeIncrementCup;

      GrainElectricCharge_NEW=t;
    }
    else {
      double minCharge,maxCharge;
      double t=ChargeIncrement+GrainElectricCharge;

      if (GrainElectricCharge>0.0) {
        minCharge=ChargeIncrementMinFractionCup*GrainElectricCharge;
        maxCharge=ChargeIncrementMaxFractionCup*GrainElectricCharge;
      }
      else {
        minCharge=ChargeIncrementMaxFractionCup*GrainElectricCharge;
        maxCharge=ChargeIncrementMinFractionCup*GrainElectricCharge;
      }

      if (t<minCharge) t=minCharge;
      else if (t>maxCharge) t=maxCharge;

      GrainElectricCharge_NEW=t;
    }

    //evaluate the convergance of the iterations
    if (fabs(GrainElectricCharge-GrainElectricCharge_NEW)<1.0E-5*fabs(GrainElectricCharge+GrainElectricCharge_NEW)) {
      GrainElectricCharge=GrainElectricCharge_NEW;
      break;
    }

    GrainElectricCharge=GrainElectricCharge_NEW;
  }
  while (++Counter<10000);


  SetGrainCharge(GrainElectricCharge,ParticleData);

  //move the particle into diferent velocity group if needed
  int oldVelocityGroup,newVelocityGroup;

  oldVelocityGroup=spec-_DUST_SPEC_;
  newVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);

  if (oldVelocityGroup!=newVelocityGroup) {
    //move the particle into different velocity group
    double GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

    GrainWeightCorrection*=finalNode->block->GetLocalTimeStep(_DUST_SPEC_+newVelocityGroup)/finalNode->block->GetLocalTimeStep(_DUST_SPEC_+oldVelocityGroup);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+oldVelocityGroup]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+newVelocityGroup];
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,ParticleData);
    spec=_DUST_SPEC_+newVelocityGroup;
    PIC::ParticleBuffer::SetI(spec,ParticleData);
  }



  return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
}
