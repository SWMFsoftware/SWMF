
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
  double GrainRadius, double GrainElectricCharge,double *GrainVelocity, //radius, the electric charge, and velocity of the grain
  double HeliocentricDistance, //heliocentric distance of the grain in AU
  double ni,double Ti, double *Vi,  //ion number density, temeprature, and velocity
  double ne, double Te,  //electron number density and temeprature
  double &Je, double &dJe, //the electron collection current
  double &Ji, double &dJi, //ion collection current
  double &Jpe, double &dJpe, //photo-electron current
  double &Jse, double &dJse //secondary-electron current
  ) {

  double  DustPotential=GetDustGrainPotential(GrainRadius,GrainElectricCharge);
  // derivative of grain potential if potential = k*charge/r
  double dDustPotential=GetDustGrainPotential(GrainRadius,1.0);

  double  XiIon     = ElectronCharge* DustPotential/(Kbol*Ti);
  double dXiIon     = ElectronCharge*dDustPotential/(Kbol*Ti);
  double  XiElectron=-ElectronCharge* DustPotential/(Kbol*Te);
  double dXiElectron=-ElectronCharge*dDustPotential/(Kbol*Te);

  // ELECRON COLLECTION CURRENT------------------------------------------------
  // Horanyi-1996-ARAA Eq. 4
#if _DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_
  double J0e=-ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*ne*sqrt(Kbol*Te/(PiTimes2*ElectronMass));

  if (XiElectron>=0.0) {
    Je  = J0e * exp(-XiElectron);
    dJe =-Je  * dXiElectron;
  }
  else {
    Je  = J0e * (1.0-XiElectron);
    dJe =-J0e * dXiElectron;
  }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
  if (!isfinite(Je)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
  if (!isfinite(dJe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif//_PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
#endif//_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_

#else  //_DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_
  Je=0.0,dJe=0.0;
#endif //_DUST__CHARGING__ELECTRON_COLLECTION__MODE_ == _PIC_MODE_ON_

  //ION COLLECTION CURRENT-----------------------------------------------------
  // Horanyi-1996-ARAA Eq. 14,15
  // NOTE: typo in Eq.15, need to divide by 2????
#if _DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_
  double J0i=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*ni*sqrt(Kbol*Ti/(PiTimes2*ProtonMass));

  // relative Mach number: dust-to-plasma rel velocity over ion thermal speed
  double M2=(pow(Vi[0]-GrainVelocity[0],2)+pow(Vi[1]-GrainVelocity[1],2)+pow(Vi[2]-GrainVelocity[2],2))/(2.0*Kbol*Ti/ProtonMass);

  //avoid division by zero
  M2 = max(1E-8,M2);
  double M = sqrt(M2);

  if (DustPotential<=0.0) {

      double erfM      = erf(M);
      double expMinM2  = exp(-M2);
      double misc      = M2 + 0.5 - XiIon;

      Ji = 0.5*J0i*( misc*sqrtPi/M*erfM + expMinM2 );
      dJi= 0.5*J0i*(-dXiIon *sqrtPi/M*erfM);
  }
  else {

      double sqrtXiIon    = sqrt(XiIon);
      double  expMinMPls2 = exp(-pow(M+sqrtXiIon,2));
      double  expMinMMin2 = exp(-pow(M-sqrtXiIon,2));
      double dexpMinMPls2 =-expMinMPls2*(M/sqrtXiIon+1);
      double dexpMinMMin2 = expMinMMin2*(M/sqrtXiIon-1);
      double  erfMPls     = erf(M+sqrtXiIon);
      double  erfMMin     = erf(M-sqrtXiIon);
      double derfMPls     = erfMPls/(sqrtPi*sqrtXiIon);
      double derfMMin     =-erfMMin/(sqrtPi*sqrtXiIon);

      double sqrtM       = sqrt(M);
      double sqrtXiIonOverM= sqrtXiIon / sqrtM;
      double sqrtXiIonTimeM= sqrtXiIon * sqrtM;
      double misc        = M*M + 0.5 - XiIon;

      Ji = 0.25*J0i*( misc*sqrtPi/M * (erfMPls+erfMMin) + 
		      (sqrtXiIonOverM+1) * expMinMMin2 - 
		      (sqrtXiIonOverM-1) * expMinMPls2 );
      dJi= 0.25*J0i*(-    sqrtPi/M*( erfMPls+ erfMMin) + 
		     misc*sqrtPi/M*(derfMPls+derfMMin) + 
		     0.5/(sqrtXiIonTimeM  )* expMinMMin2 + 
		     (    sqrtXiIonOverM+1)*dexpMinMMin2 -
		     0.5/(sqrtXiIonTimeM  )* expMinMPls2 - 
		     (    sqrtXiIonOverM-1)*dexpMinMPls2
		     )*dXiIon;
	}

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
  if (!isfinite(Ji)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJi)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

  #else //_DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_
   Ji=0.0,dJi=0.0;
  #endif //_DUST__CHARGING__ION_COLLECTION__MODE_ == _PIC_MODE_ON_

  //PHOTO-ELECTRON CURRENT-----------------------------------------------------
  // Horanyi-1996-ARAA Eq. 13
  #if _DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  double J0pe=Pi*GrainRadius*GrainRadius*ElectronCharge*ElectonPhotoEmission::PhotoElectronEfficiency*ElectonPhotoEmission::PhotoElectronEfficiencyMaterialConstant/pow(HeliocentricDistance,2);

  if (DustPotential<=0.0) {
    Jpe=J0pe,dJpe=0.0;
  }
  else {
    double  misc = -ElectronCharge* DustPotential / ElectonPhotoEmission::PhotoElectronEnergy;
    double dmisc = -ElectronCharge*dDustPotential / ElectonPhotoEmission::PhotoElectronEnergy;
    double expmisc= exp(misc); 
    Jpe = J0pe * expmisc;
    dJpe= Jpe  * dmisc;
  }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(Jpe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJpe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

  #else //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  Jpe=0.0,dJpe=0.0;
  #endif //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_


  //SECONDARY ELECTON CURRENT--------------------------------------------------
  // Horanyi-1996-ARAA Eq. 10-12
  #if _DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  if (DustPotential<=0.0) {

    double misc     = 0.25*SecondaryEmissionPeakYieldEnergy/(Kbol*Te);
    double sqrtmisc = sqrt(misc);
    double exp4bymisc  = exp(0.25/misc);
    double F5 = (misc>1E-3)?
      ( 
       2*sqrtmisc*(1+2*misc)*(1+16*misc) - 
       exp4bymisc*sqrtPi*(1+20*misc*(1+3*misc))*erfc(0.5/sqrtmisc)
	) / (64*pow(sqrtmisc,7))
      :
      (120-5040*misc)*misc*misc;

    Jse =-3.7*SecondaryEmissionPeakYield*J0e*exp(-XiElectron)*F5;
    dJse=-Jse * dXiElectron;
  }
  else {

    double  XiS     = -ElectronCharge * DustPotential / (Kbol*SecondaryElectronTemperature);
    double dXiS     = -ElectronCharge *dDustPotential / (Kbol*SecondaryElectronTemperature);
    double expDifXi = exp(XiS-XiElectron);
    double misc     = 0.25*SecondaryEmissionPeakYieldEnergy/(Kbol*Te);
    double sqrtmisc = sqrt(misc);
    double exp4bymisc  = exp(0.25/misc);

    double B = sqrt(-XiElectron / misc);
    double expMinB = exp(-B);
    double expMinB2misc = exp(-B*B*misc);
    double F5B=(misc>1E-3)?
      (expMinB*expMinB2misc*2*sqrtmisc*
       (1+2*misc*(
		  9+16*misc+B*(
			       -1+2*misc*(
					  -7+B-2*(-4+B)*B*misc + 
					  4*B*B*B*misc*misc
					  ))))
       +
       exp4bymisc*sqrtPi*(1+20*misc*(1+3*misc))*erfc((0.5+B*misc)/sqrtmisc)
       )/(64*pow(sqrtmisc,7))
      :
      misc*misc*expMinB*
      (( 120+B*( 120+B*(  60+B*( 20+B*(  5+B))))) - 
       (5040+B*(5040+B*(2520+B*(840+B*(210+B*(42+B*(7+B)))))))*misc);

    Jse =-3.7*SecondaryEmissionPeakYield*J0e*F5B*(1-XiS)*expDifXi;

    dJse=-3.7*SecondaryEmissionPeakYield*J0e*F5B*expDifXi*
      (-dXiS+(1-XiS)*(dXiS-dXiElectron));
  }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(Jse)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJse)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

  #else //_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  Jse=0.0,dJse=0.0;
  #endif //_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
}


//the generic procedure of the grain charge update
int ElectricallyChargedDust::Charging::UpdateGrainCharge(
                double  dt,
		double* ParticleVelocity,
		double* PlasmaVelocity,
		double  PlasmaPressure,
                double  PlasmaTemperature,
		double  PlasmaNumberDensity,
		double  GrainRadius, 
		double& GrainElectricCharge,
		int     CHARGE_INTEGRATION_MODE) {

  double GrainElectricCharge_NEW,InitGrainCharge;
  
  //initial grain's charge
  InitGrainCharge=GrainElectricCharge;


  //reserve space for different electron and ion temepratures
  double Ti,Te,Je,dJe,Ji,dJi,Jpe,dJpe,Jse,dJse,pe;

  if (PlasmaNumberDensity<1.0E2) {
    PlasmaTemperature=200.0;
    PlasmaNumberDensity=1.0E2;
    PlasmaVelocity[0]=100.0,PlasmaVelocity[1]=0.0,PlasmaVelocity[2]=0.0;
    Ti=PlasmaTemperature;
    Te=PlasmaTemperature;
  }
  else{
    Ti=PlasmaTemperature;
    Te=PlasmaPressure/(Kbol*PlasmaNumberDensity);
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
  int nIter=10000;
  static const double ChargeIncrementCup=100.0*ElectronCharge;
  static const double ChargeIncrementMaxFractionCup=3.0;
  static const double ChargeIncrementMinFractionCup=1.0/ChargeIncrementMaxFractionCup;


  do {
    GetGrainCurrent(
      GrainRadius,GrainElectricCharge,ParticleVelocity, //radius, the electric charge, and velocity of the grain
      HeliocentricDistance, //heliocentric distance of the grain in AU
      PlasmaNumberDensity,Ti,PlasmaVelocity,  //ion number density, temeprature, and velocity
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
      //ChargeIncrement=InitGrainCharge-dt*Jtotal/dJtotal;
      ChargeIncrement= 
	-((GrainElectricCharge-InitGrainCharge) - Jtotal*dt) / (1-dJtotal*dt);
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
  while (++Counter<nIter);

}



//the generic procedure of the grain charge update
int ElectricallyChargedDust::Charging::UpdateGrainCharge(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode,int CHARGE_INTEGRATION_MODE) {
  double PlasmaTemperature,PlasmaNumberDensity, pe;
  double GrainElectricCharge,GrainElectricCharge_NEW,InitGrainCharge;

  PIC::Mesh::cDataCenterNode* cell;
  int i,j,k;
  long int LocalCellNumber;
  double swVel[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *finalNode=PIC::Mesh::mesh.findTreeNode(xFinal,initNode);

  if (_DUST__CHARGING_MODE_ == _DUST__CHARGING_MODE__OFF_) {
    //dust charing modeling is turned off
    return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
  }

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


  PIC::CPLR::InitInterpolationStencil(xInit,initNode);
  PlasmaTemperature=PIC::CPLR::GetBackgroundPlasmaTemperature();
  PIC::CPLR::GetBackgroundPlasmaVelocity(swVel);
  PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
  pe=PIC::CPLR::GetBackgroundElectronPlasmaPressure();


  UpdateGrainCharge(dt, v, 
		    swVel,
		    pe, PlasmaTemperature, PlasmaNumberDensity,
		    GrainRadius, GrainElectricCharge, CHARGE_INTEGRATION_MODE);


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
