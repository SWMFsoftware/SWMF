
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

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(Je)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

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

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(Ji)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJi)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

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

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(Jpe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      if (!isfinite(dJpe)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

  #else //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  Jpe=0.0,dJpe=0.0;
  #endif //_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_


  //SECONDARY ELECTON CURRENT
  #if _DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_ == _PIC_MODE_ON_
  if (DustPotential<=0.0) {

    {
      double t13 = exp(ElectronCharge*GrainElectricCharge/0.3141592653589793E1/
VacuumPermittivity/GrainRadius/Kbol/Te/4.0);
      double t16 = SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYieldEnergy;
      double t18 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t22 = sqrt(0.3141592653589793E1);
      double t27 = exp(1/SecondaryEmissionPeakYieldEnergy*Kbol*Te);
      double t30 = sqrt(Kbol);
      double t32 = sqrt(Te);
      double t34 = erf(t32*t30/t18);
      double t39 = Te*Te;
      double t41 = Kbol*Kbol;
      double t43 = t30*t41*t32*t39;
      double t44 = t27*t22;
      double t59 = t30*Kbol*t32*Te;
      Jse = 0.185E1*(15.0*t32*t30*t34*t27*t22*t16+4.0*t34*t44*t43+8.0*t18*t16+
4.0*t18*t39*t41+18.0*Kbol*Te*t18*SecondaryEmissionPeakYieldEnergy+20.0*
SecondaryEmissionPeakYieldEnergy*t34*t27*t22*t59-4.0*t44*t43-20.0*
SecondaryEmissionPeakYieldEnergy*t44*t59-15.0*t32*t30*t16*t44)*Te*Kbol/t18/t16/
SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYield*t13*J0e;
    }

    {
      double t2 = 1/0.3141592653589793E1;
      double t3 = 1/VacuumPermittivity;
      double t6 = 1/GrainRadius;
      double t16 = exp(1/Kbol/Te*t6*t3*t2*GrainElectricCharge*ElectronCharge/4.0);
      double t18 = SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYieldEnergy;
      double t20 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t24 = sqrt(0.3141592653589793E1);
      double t29 = exp(1/SecondaryEmissionPeakYieldEnergy*Kbol*Te);
      double t32 = sqrt(Kbol);
      double t34 = sqrt(Te);
      double t36 = erf(t34*t32/t20);
      double t41 = Te*Te;
      double t43 = Kbol*Kbol;
      double t45 = t32*t43*t34*t41;
      double t46 = t29*t24;
      double t61 = t32*Kbol*t34*Te;
      dJse = 0.4625*(15.0*t34*t32*t36*t29*t24*t18+4.0*t36*t46*t45+8.0*t20*t18+
4.0*t20*t43*t41+18.0*Kbol*Te*t20*SecondaryEmissionPeakYieldEnergy+20.0*
SecondaryEmissionPeakYieldEnergy*t36*t29*t24*t61-4.0*t46*t45-20.0*
SecondaryEmissionPeakYieldEnergy*t46*t61-15.0*t34*t32*t18*t46)/t20/t18/
SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYield*t16*t6*t3*t2*J0e*
ElectronCharge;
    }

  }
  else {

    {
      double t1 = ElectronCharge*GrainElectricCharge;
      double t6 = 1/VacuumPermittivity/GrainRadius;
      double t7 = 1/Kbol;
      double t8 = 1/Te;
      double t12 = t8*t7*t6/0.3141592653589793E1*t1/4.0;
      double t15 = SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYieldEnergy;
      double t17 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t23 = Kbol*Kbol;
      double t24 = sqrt(Kbol);
      double t26 = Te*Te;
      double t27 = sqrt(Te);
      double t29 = t27*t26*t24*t23;
      double t30 = 0.3141592653589793E1*0.3141592653589793E1;
      double t31 = sqrt(0.3141592653589793E1);
      double t32 = t31*t30;
      double t33 = Kbol*Te;
      double t36 = t6*t1;
      double t37 = sqrt(t36);
      double t38 = SecondaryEmissionPeakYieldEnergy*t37;
      double t40 = t17*SecondaryEmissionPeakYieldEnergy;
      double t43 = 1/t31;
      double t45 = exp(t43/t40*(t17*t31*t33+t38));
      double t46 = t45*t32;
      double t56 = 1/t17;
      double t60 = erf(t56/t27/t24*t43*(t37*t17+2.0*t31*t33)/2.0);
      double t62 = VacuumPermittivity*VacuumPermittivity;
      double t63 = GrainRadius*GrainRadius;
      double t64 = t63*t62;
      double t68 = t23*Kbol;
      double t69 = t26*Te;
      double t70 = t69*t68;
      double t71 = exp(-t12);
      double t73 = t31*0.3141592653589793E1*t71;
      double t78 = t26*t23;
      double t84 = t23*t23;
      double t86 = t26*t26;
      double t89 = t32*t27*t86*t24*t84;
      double t94 = ElectronCharge*t71;
      double t98 = 0.3141592653589793E1*VacuumPermittivity*GrainRadius;
      double t115 = t17*t15;
      double t122 = t27*t69*t24*t68;
      double t139 = ElectronCharge*ElectronCharge;
      double t141 = GrainElectricCharge*GrainElectricCharge;
      double t144 = t62*t45;
      double t158 = 60.0*t64*t15*t60*t46*t29-8.0*t64*t38*t73*t70-28.0*t64*t37*t15*t73*
t78+16.0*t64*t60*t45*t89+4.0*t98*t40*GrainElectricCharge*t94*t78+72.0*t64*t30*
t40*t71*t70-2.0*t64*t37*t36*t15*t31*t71*t33+32.0*t64*t30*t115*t71*t78+80.0*t64*
SecondaryEmissionPeakYieldEnergy*t60*t46*t122+16.0*t64*t30*t17*t71*t86*t84+8.0*
t98*t115*GrainElectricCharge*t94*t33+t115*t141*t139*t71-16.0*t63*t144*t89-80.0*
SecondaryEmissionPeakYieldEnergy*t63*t144*t32*t122-60.0*t15*t63*t144*t32*t29;
      double t162 = exp(-t37*t56*t43);
      Jse = 0.4625/t63/t62/t30*t162*t158*t8*t7/t17/t15/
SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYield*(1.0+t12)*J0e;

    }

    {
      double t2 = 0.3141592653589793E1*0.3141592653589793E1;
      double t5 = VacuumPermittivity*VacuumPermittivity;
      double t7 = 1/t5/VacuumPermittivity;
      double t9 = GrainRadius*GrainRadius;
      double t11 = 1/t9/GrainRadius;
      double t14 = Kbol*Kbol;
      double t16 = Te*Te;
      double t20 = SecondaryEmissionPeakYieldEnergy*SecondaryEmissionPeakYieldEnergy;
      double t22 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t24 = 1/t22/t20/SecondaryEmissionPeakYieldEnergy;
      double t25 = sqrt(Kbol);
      double t27 = sqrt(Te);
      double t29 = t27*t16*t25*t14;
      double t30 = sqrt(0.3141592653589793E1);
      double t31 = t30*t2;
      double t32 = Kbol*Te;
      double t35 = ElectronCharge*GrainElectricCharge;
      double t38 = 1/VacuumPermittivity/GrainRadius;
      double t39 = t38*t35;
      double t40 = sqrt(t39);
      double t41 = SecondaryEmissionPeakYieldEnergy*t40;
      double t43 = t22*SecondaryEmissionPeakYieldEnergy;
      double t46 = 1/t30;
      double t48 = exp(t46/t43*(t22*t30*t32+t41));
      double t49 = t48*t31;
      double t54 = t40*t22+2.0*t30*t32;
      double t59 = 1/t22;
      double t63 = erf(t59/t27/t25*t46*t54/2.0);
      double t65 = t9*t5;
      double t69 = t14*Kbol;
      double t70 = t16*Te;
      double t71 = t70*t69;
      double t72 = 1/0.3141592653589793E1;
      double t74 = 1/Kbol;
      double t75 = 1/Te;
      double t76 = t75*t74;
      double t79 = t76*t38*t72*t35/4.0;
      double t80 = exp(-t79);
      double t81 = t30*0.3141592653589793E1;
      double t82 = t81*t80;
      double t83 = t82*t71;
      double t87 = t16*t14;
      double t88 = t82*t87;
      double t93 = t14*t14;
      double t95 = t16*t16;
      double t97 = t27*t95*t25*t93;
      double t98 = t31*t97;
      double t103 = ElectronCharge*t80;
      double t104 = t103*t87;
      double t107 = 0.3141592653589793E1*VacuumPermittivity*GrainRadius;
      double t119 = t40*t39*t20;
      double t124 = t22*t20;
      double t131 = t27*t70*t25*t69;
      double t137 = t95*t93;
      double t148 = ElectronCharge*ElectronCharge;
      double t150 = GrainElectricCharge*GrainElectricCharge;
      double t153 = t5*t48;
      double t167 = 60.0*t65*t20*t63*t49*t29-8.0*t65*t41*t83-28.0*t65*t40*t20*t88+16.0
*t65*t63*t48*t98+4.0*t107*t43*GrainElectricCharge*t104+72.0*t65*t2*t43*t80*t71
-2.0*t65*t119*t30*t80*t32+32.0*t65*t2*t124*t80*t87+80.0*t65*
SecondaryEmissionPeakYieldEnergy*t63*t49*t131+16.0*t65*t2*t22*t80*t137+8.0*t107
*t124*GrainElectricCharge*t103*t32+t124*t150*t148*t80-16.0*t9*t153*t98-80.0*
SecondaryEmissionPeakYieldEnergy*t9*t153*t31*t131-60.0*t20*t9*t153*t31*t29;
      double t171 = exp(-t40*t59*t46);
      double t177 = (1.0+t79)*J0e;
      double t181 = 1/t40;
      double t182 = t181*t2;
      double t188 = t63*t48*GrainRadius*ElectronCharge*VacuumPermittivity;
      double t191 = t48*t81;
      double t192 = t54*t54;
      double t198 = exp(-1/SecondaryEmissionPeakYieldEnergy*t76*t72*t192/4.0);
      double t199 = t198*t191;
      double t202 = VacuumPermittivity*GrainRadius;
      double t204 = t20*t202*ElectronCharge*t181;
      double t207 = t30*ElectronCharge;
      double t215 = ElectronCharge*t202;
      double t216 = t215*t181*SecondaryEmissionPeakYieldEnergy;
      double t273 = t48*t202;
      double t287 = 30.0*t188*t43*t182*t29+30.0*t204*t199*t87+2.0*t40*
SecondaryEmissionPeakYieldEnergy*t80*t202*t207*t87-4.0*t216*t83+4.0*t40*t20*t80
*t202*t207*t32-14.0*t204*t88+8.0*t188*t59*t182*t97+8.0*t215*t181*t198*t191*t137
-t43*GrainElectricCharge*t80*t148*t32-14.0*t202*0.3141592653589793E1*t43*t104+
t119*t80*GrainRadius*VacuumPermittivity*t46*ElectronCharge/2.0+40.0*t188*t22*
t182*t131+40.0*t216*t199*t71-4.0*t22*t80*t202*ElectronCharge*
0.3141592653589793E1*t71-t124*t150*t80*t76*t38*t72*t148*ElectronCharge/4.0-8.0*
t273*ElectronCharge*t59*t182*t97-40.0*t273*ElectronCharge*t22*t182*t131-30.0*
t273*ElectronCharge*t43*t182*t29;
      double t299 = t20*t20;
      dJse = 0.115625*t171*t167*t24*SecondaryEmissionPeakYield/t16/t14*t11*t7/
t2/0.3141592653589793E1*J0e*ElectronCharge+0.4625/t9/t5/t2*t171*t287*t75*t74*
t24*SecondaryEmissionPeakYield*t177-0.23125*t171/t31*t11*t7*ElectronCharge*t181
*t167*t75*t74/t299*SecondaryEmissionPeakYield*t177;

    }

  }

  Jse*=-1.0;
  dJse*=-1.0;

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
