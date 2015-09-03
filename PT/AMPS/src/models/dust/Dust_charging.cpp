
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
  if (DustPotential<=0.0) {
    {
      double t13 = exp(ElectronCharge * GrainElectricCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1);
      double t16 = SecondaryEmissionPeakYieldEnergy * SecondaryEmissionPeakYieldEnergy;
      double t18 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t22 = Kbol * Kbol;
      double t23 = sqrt(Kbol);
      double t25 = Te * Te;
      double t26 = sqrt(Te);
      double t28 = t23 * t22 * t26 * t25;
      double t29 = sqrt(0.3141592654e1);
      double t33 = exp(0.1e1 / SecondaryEmissionPeakYieldEnergy * Kbol * Te);
      double t34 = t29 * t33;
      double t38 = erf(0.1e1 / t18 * t26 * t23);
      double t54 = t23 * Kbol * t26 * Te;
      Jse = 0.1850000000e1 * J0e * t13 * SecondaryEmissionPeakYield / t18 / t16 / SecondaryEmissionPeakYieldEnergy * Kbol * Te * (0.4e1 * t28 * t34 * t38 + 0.18e2 * Kbol * Te * t18 * SecondaryEmissionPeakYieldEnergy +
          0.15e2 * t16 * t29 * t33 * t38 * t26 * t23 + 0.20e2 * t54 * t29 * t33 * t38 * SecondaryEmissionPeakYieldEnergy + 0.4e1 * t22 * t25 * t18 + 0.8e1 * t18 * t16 -
          0.4e1 * t28 * t34 - 0.15e2 * t34 * t16 * t26 * t23 - 0.20e2 * t54 * t34 * SecondaryEmissionPeakYieldEnergy);
    }

    {
      double t2 = 0.1e1 / 0.3141592654e1;
      double t3 = 0.1e1 / VacuumPermittivity;
      double t6 = 0.1e1 / GrainRadius;
      double t16 = exp(ElectronCharge * GrainElectricCharge * t2 * t3 * t6 / Kbol / Te / 0.4e1);
      double t18 = SecondaryEmissionPeakYieldEnergy * SecondaryEmissionPeakYieldEnergy;
      double t20 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t24 = Kbol * Kbol;
      double t25 = sqrt(Kbol);
      double t27 = Te * Te;
      double t28 = sqrt(Te);
      double t30 = t25 * t24 * t28 * t27;
      double t31 = sqrt(0.3141592654e1);
      double t35 = exp(0.1e1 / SecondaryEmissionPeakYieldEnergy * Kbol * Te);
      double t36 = t31 * t35;
      double t40 = erf(0.1e1 / t20 * t25 * t28);
      double t56 = t25 * Kbol * t28 * Te;
      dJse = 0.4625000000e0 * J0e * ElectronCharge * t2 * t3 * t6 * t16 * SecondaryEmissionPeakYield / t20 / t18 / SecondaryEmissionPeakYieldEnergy * (0.4e1 * t30 * t36 * t40 + 0.18e2 * Kbol * Te * t20 * SecondaryEmissionPeakYieldEnergy +
          0.15e2 * t18 * t31 * t35 * t40 * t25 * t28 + 0.20e2 * t56 * t31 * t35 * t40 * SecondaryEmissionPeakYieldEnergy + 0.4e1 * t24 * t27 * t20 + 0.8e1 * t20 * t18 - 0.4e1 * t30 * t36 - 0.20e2 * t56 * t36 * SecondaryEmissionPeakYieldEnergy -
          0.15e2 * t36 * t18 * t25 * t28);
    }
  }
  else {

    {
      double t1 = ElectronCharge * GrainElectricCharge;
      double t2 = 0.1e1 / 0.3141592654e1;
      double t3 = t1 * t2;
      double t4 = 0.1e1 / VacuumPermittivity;
      double t5 = 0.1e1 / GrainRadius;
      double t6 = t4 * t5;
      double t7 = 0.1e1 / Kbol;
      double t8 = 0.1e1 / Te;
      double t21 = exp(-t3 * t6 * t7 * (0.1e1 / SecondaryElectronTemeprature - t8) / 0.4e1);
      double t23 = SecondaryEmissionPeakYieldEnergy * SecondaryEmissionPeakYieldEnergy;
      double t25 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t31 = Kbol * Kbol;
      double t32 = t31 * t31;
      double t33 = sqrt(Kbol);
      double t35 = Te * Te;
      double t36 = t35 * t35;
      double t37 = sqrt(Te);
      double t40 = sqrt(0.3141592654e1);
      double t41 = t40 * 0.3141592654e1;
      double t42 = t33 * t32 * t37 * t36 * t41;
      double t43 = t31 * t35;
      double t44 = pow(SecondaryEmissionPeakYieldEnergy, 0.1e1 / 0.4e1);
      double t45 = t44 * SecondaryEmissionPeakYieldEnergy;
      double t46 = pow(0.3141592654e1, 0.1e1 / 0.4e1);
      double t47 = t46 * t46;
      double t48 = t46 * t47;
      double t51 = t1 * t6;
      double t52 = pow(t51, 0.1e1 / 0.4e1);
      double t53 = t25 * SecondaryEmissionPeakYieldEnergy;
      double t56 = t40 * Te * Kbol;
      double t59 = sqrt(t51);
      double t60 = t44 * t44;
      double t61 = t60 * t44;
      double t62 = t61 * SecondaryEmissionPeakYieldEnergy;
      double t71 = 0.1e1 / t48 * t8 * t7;
      double t73 = exp((t43 * t45 * t48 + 0.2e1 * t52 * t53 * t56 + t59 * t62 * t46) / t44 / t23 * t71);
      double t75 = Kbol * Te;
      double t83 = 0.1e1 / t25;
      double t86 = erf((t44 * t52 + t75 * t46) / t46 / t33 / t37 * t83);
      double t88 = VacuumPermittivity * GrainRadius;
      double t92 = t31 * Kbol;
      double t93 = t35 * Te;
      double t94 = t92 * t93;
      double t97 = t48 * t52 * t88;
      double t108 = t52 * t52;
      double t117 = 0.3141592654e1 * VacuumPermittivity * GrainRadius;
      double t129 = t33 * t31 * t37 * t35;
      double t130 = t41 * t73;
      double t141 = t33 * t92 * t37 * t93;
      double t151 = t73 * VacuumPermittivity;
      double t165 = 0.4e1 * t42 * t73 * t86 * t88 - 0.4e1 * t94 * t61 * t97 + 0.4e1 * t1 * t53 + 0.4e1 * t43 * t40 * SecondaryEmissionPeakYieldEnergy * t59 * t88 - 0.4e1 * t75 * t45 * t46 * t108 * t52 * t88 +
          0.8e1 * t43 * t25 * t23 * t117 - 0.14e2 * t43 * t62 * t97 + 0.4e1 * t32 * t36 * t25 * t117 + 0.15e2 * t129 * t130 * t86 * t23 * t88 + 0.18e2 * t94 * t53 * t117 + 0.20e2 * t141 * t130 * t86 * SecondaryEmissionPeakYieldEnergy * t88 +
          0.8e1 * t56 * t23 * t59 * t88 - 0.4e1 * t42 * t151 * GrainRadius - 0.20e2 * t141 * t41 * t151 * GrainRadius * SecondaryEmissionPeakYieldEnergy - 0.15e2 * t129 * t41 * t151 * GrainRadius * t23;
      double t178 = exp(-(0.2e1 * t52 * t40 * t25 * Te * Kbol + t59 * t61 * t46) / t45 * t71);
      double t182 = 0.1e1 / t40 * t83;
      double t186 = exp(t182 * t59 * t8 * t7);
      double t187 = 0.1e1 / t186;
      Jse=0.1850000000e1 * J0e * (0.1e1 + t3 * t6 * t7 * t8 / 0.4e1) * t21 * SecondaryEmissionPeakYield / t25 / t23 / SecondaryEmissionPeakYieldEnergy * Kbol * Te * t165 * t178 * t2 * t4 * t5 / (t43 * t187 + t75 * t187 * t182 * t59);
    }

    {
      double t2 = 0.3141592654e1 * 0.3141592654e1;
      double t4 = VacuumPermittivity * VacuumPermittivity;
      double t6 = 0.1e1 / t2 / t4;
      double t7 = GrainRadius * GrainRadius;
      double t8 = 0.1e1 / t7;
      double t11 = ElectronCharge * GrainElectricCharge;
      double t12 = 0.1e1 / 0.3141592654e1;
      double t13 = t11 * t12;
      double t14 = 0.1e1 / VacuumPermittivity;
      double t15 = 0.1e1 / GrainRadius;
      double t16 = t14 * t15;
      double t17 = 0.1e1 / Kbol;
      double t19 = 0.1e1 / Te;
      double t20 = 0.1e1 / SecondaryElectronTemeprature - t19;
      double t25 = exp(-t13 * t16 * t17 * t20 / 0.4e1);
      double t27 = SecondaryEmissionPeakYieldEnergy * SecondaryEmissionPeakYieldEnergy;
      double t29 = sqrt(SecondaryEmissionPeakYieldEnergy);
      double t31 = 0.1e1 / t29 / t27 / SecondaryEmissionPeakYieldEnergy;
      double t32 = t25 * SecondaryEmissionPeakYield * t31;
      double t33 = Kbol * Kbol;
      double t34 = t33 * t33;
      double t35 = sqrt(Kbol);
      double t37 = Te * Te;
      double t38 = t37 * t37;
      double t39 = sqrt(Te);
      double t42 = sqrt(0.3141592654e1);
      double t43 = t42 * 0.3141592654e1;
      double t44 = t35 * t34 * t39 * t38 * t43;
      double t45 = t33 * t37;
      double t46 = pow(SecondaryEmissionPeakYieldEnergy, 0.1e1 / 0.4e1);
      double t47 = t46 * SecondaryEmissionPeakYieldEnergy;
      double t48 = pow(0.3141592654e1, 0.1e1 / 0.4e1);
      double t49 = t48 * t48;
      double t50 = t49 * t48;
      double t53 = t11 * t16;
      double t54 = pow(t53, 0.1e1 / 0.4e1);
      double t55 = t29 * SecondaryEmissionPeakYieldEnergy;
      double t57 = t42 * Te;
      double t58 = t57 * Kbol;
      double t61 = sqrt(t53);
      double t62 = t46 * t46;
      double t63 = t62 * t46;
      double t64 = t63 * SecondaryEmissionPeakYieldEnergy;
      double t69 = 0.1e1 / t46 / t27;
      double t73 = 0.1e1 / t50 * t19 * t17;
      double t75 = exp((t45 * t47 * t50 + 0.2e1 * t54 * t55 * t58 + t61 * t64 * t48) * t69 * t73);
      double t77 = Kbol * Te;
      double t79 = t46 * t54 + t77 * t48;
      double t85 = 0.1e1 / t29;
      double t88 = erf(t79 / t48 / t35 / t39 * t85);
      double t90 = VacuumPermittivity * GrainRadius;
      double t94 = t33 * Kbol;
      double t95 = t37 * Te;
      double t96 = t94 * t95;
      double t97 = t96 * t63;
      double t99 = t50 * t54 * t90;
      double t104 = t45 * t42;
      double t109 = t77 * t47;
      double t110 = t54 * t54;
      double t111 = t110 * t54;
      double t119 = 0.3141592654e1 * VacuumPermittivity * GrainRadius;
      double t122 = t45 * t64;
      double t125 = t34 * t38;
      double t131 = t35 * t33 * t39 * t37;
      double t132 = t43 * t75;
      double t143 = t35 * t94 * t39 * t95;
      double t153 = t75 * VacuumPermittivity;
      double t167 = 0.4e1 * t44 * t75 * t88 * t90 - 0.4e1 * t97 * t99 + 0.4e1 * t11 * t55 + 0.4e1 * t104 * SecondaryEmissionPeakYieldEnergy * t61 * t90 -
          0.4e1 * t109 * t48 * t111 * t90 + 0.8e1 * t45 * t29 * t27 * t119 - 0.14e2 * t122 * t99 + 0.4e1 * t125 * t29 * t119 + 0.15e2 * t131 * t132 * t88 * t27 * t90 +
          0.18e2 * t96 * t55 * t119 + 0.20e2 * t143 * t132 * t88 * SecondaryEmissionPeakYieldEnergy * t90 + 0.8e1 * t58 * t27 * t61 * t90 - 0.4e1 * t44 * t153 * GrainRadius -
          0.20e2 * t143 * t43 * t153 * GrainRadius * SecondaryEmissionPeakYieldEnergy - 0.15e2 * t131 * t43 * t153 * GrainRadius * t27;
      double t169 = t29 * Te;
      double t176 = 0.1e1 / t47;
      double t179 = exp(-(0.2e1 * t54 * t42 * t169 * Kbol + t61 * t63 * t48) * t176 * t73);
      double t180 = t167 * t179;
      double t181 = 0.1e1 / t42;
      double t182 = t181 * t85;
      double t186 = exp(t182 * t61 * t19 * t17);
      double t187 = 0.1e1 / t186;
      double t192 = t45 * t187 + t77 * t187 * t182 * t61;
      double t193 = 0.1e1 / t192;
      double t198 = t17 * t19;
      double t203 = J0e * (0.1e1 + t13 * t16 * t198 / 0.4e1);
      double t214 = t203 * t25;
      double t218 = 0.1e1 / t111;
      double t222 = Kbol * ElectronCharge * t16;
      double t224 = 0.1e1 / t61;
      double t228 = ElectronCharge * t14 * t15;
      double t231 = t50 * (t218 * t55 * t57 * t222 + t224 * t64 * t48 * t228) / 0.2e1;
      double t232 = t143 * t231;
      double t233 = t69 * t75;
      double t235 = t88 * VacuumPermittivity * GrainRadius;
      double t239 = t50 * t75;
      double t241 = t79 * t79;
      double t246 = exp(-t241 * t181 * t198 / SecondaryEmissionPeakYieldEnergy);
      double t247 = 0.1e1 / t46;
      double t249 = t218 * ElectronCharge;
      double t254 = t50 * t218 * ElectronCharge;
      double t272 = t35 * Kbol * t39 * Te * t231;
      double t273 = t247 * t75;
      double t282 = t131 * t231;
      double t283 = t176 * t75;
      double t305 = 0.4e1 * t232 * t233 * t235 + 0.2e1 * t125 * t239 * t246 * t247 * t249 - t97 * t254 + 0.4e1 * ElectronCharge * t55 +
          0.2e1 * t104 * SecondaryEmissionPeakYieldEnergy * t224 * ElectronCharge - 0.3e1 * t109 * t48 / t54 * ElectronCharge -
          0.7e1 / 0.2e1 * t122 * t254 + 0.15e2 * t272 * t273 * t235 + 0.15e2 / 0.2e1 * t45 * t239 * t246 * t64 * t249 +
          0.20e2 * t282 * t283 * t235 + 0.10e2 * t96 * t239 * t246 * t63 * t249 + 0.4e1 * t58 * t27 * t224 * ElectronCharge -
          0.4e1 * t232 * t233 * t90 - 0.20e2 * t282 * t283 * t90 - 0.15e2 * t272 * t273 * t90;
      double t314 = t27 * t27;
      double t341 = t192 * t192;
      dJse = 0.4625000000e0 * J0e * ElectronCharge * t6 * t8 * t32 * t180 * t193 -
          0.4625000000e0 * t203 * ElectronCharge * t6 * t8 * t20 * t32 * Te * t167 * t179 * t193 +
          0.1850000000e1 * t214 * SecondaryEmissionPeakYield * t31 * Kbol * Te * t305 * t179 * t12 * t14 * t15 * t193 -
          0.9250000000e0 * t214 * SecondaryEmissionPeakYield / t63 / t314 * t167 * (t218 * t42 * t169 * t222 + t224 * t63 * t48 * t228) / t50 / 0.3141592654e1 * t179 * t16 * t193 +
          0.9250000000e0 * t214 * SecondaryEmissionPeakYield / t29 / t314 * t77 * t180 * t6 * t8 / t341 * t187 * ElectronCharge;

    }

  }
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
