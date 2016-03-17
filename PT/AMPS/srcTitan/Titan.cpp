/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */


#include "pic.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Titan";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_TITAN";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="MSGR_MSO";

char Titan::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";


//sticking probability
double Exosphere::SurfaceInteraction::StickingProbability(int spec, double& ReemissionParticleFraction,double Temp) {



  ReemissionParticleFraction=1.0;
  return 1.0;

  if (Temp<300.0) return 1.0;
  if (Temp<650.0) return 1.0-(Temp-300.0)/350.0;

  return 0.0;
}

//surface temeprature of the planet
double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
  static const double Tn=110.0;
  static const double Td0_Aphelion=590.0,Td0_Perihelion=725.0;
  static const double TAA_Aphelion=Pi,TAA_Perihelion=0.0;
  static const double alpha=(Td0_Aphelion-Td0_Perihelion)/(TAA_Aphelion-TAA_Perihelion);

  double Td,Angle;

  Angle=(OrbitalMotion::TAA<Pi) ? OrbitalMotion::TAA : 2.0*Pi-OrbitalMotion::TAA;
  Td=Td0_Perihelion+alpha*(Angle-TAA_Perihelion);

  return (CosSubSolarAngle>0.0) ? Tn+(Td-Tn)*pow(CosSubSolarAngle,0.25) : Tn;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\",  \"Mean Speed Along the Line of Sight(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec),PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=2;
  }

  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (res[cnt]>0.0) res[cnt+1]/=res[cnt];
    cnt+=2;
  }
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;

  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double BulkVelocity_SO[3],v_LOCAL_SO_FROZEN[3],rHeliocentric,vHeliocentric;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity_SO,spec);

    v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+BulkVelocity_SO[0]+
        Exosphere::RotationVector_SO_FROZEN[1]*x[2]-Exosphere::RotationVector_SO_FROZEN[2]*x[1];

    v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+BulkVelocity_SO[1]-
        Exosphere::RotationVector_SO_FROZEN[0]*x[2]+Exosphere::RotationVector_SO_FROZEN[2]*x[0];

    v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+BulkVelocity_SO[2]+
        Exosphere::RotationVector_SO_FROZEN[0]*x[1]-Exosphere::RotationVector_SO_FROZEN[1]*x[0];

    rHeliocentric=sqrt(pow(x[0]-Exosphere::xObjectRadial,2)+(x[1]*x[1])+(x[2]*x[2]));
    vHeliocentric=(
        (v_LOCAL_SO_FROZEN[0]*(x[0]-Exosphere::xObjectRadial))+
        (v_LOCAL_SO_FROZEN[1]*x[1])+(v_LOCAL_SO_FROZEN[2]*x[2]))/rHeliocentric;


    //brightness of the exospheric sodium
    if (spec==_NA_SPEC_) {
      //check if the point is outside of Mercury's shadow
      if ( /*(Moon::EarthShadowCheck(x)==false) &&*/ ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))) ) {
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
      else res[cnt++]=0.0,res[cnt++]=0.0;
    }

  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  return GetTAA("TITAN","Sun",_MASS_(_SUN_),EphemerisTime);
}

//extract "hot" particles from the "thermal" population
void Titan::SpeciesEnergySeparation::Process(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int spec;

  spec=PIC::ParticleBuffer::GetI(ptr);

  if (spec==_N2_SPEC_) {
    //process only N2 thermal particles;
    //when a N2 model particle that is too energyzed is determined move it to the "hot" species

    double v[3],speed;

    PIC::ParticleBuffer::GetV(v,ptr);
    speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    //the "hot" criterion: the particle speed is 3 times of the escape speed
    if (speed>=3.0*2.0E3) {
      //the particle is "hot" -> move it to the "hot" category
      double ThermalWeight,ThermalTimeStep,HotTimeStep,HotWeight;
      double nHotParticles;

      #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
      ThermalWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
      HotWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_N2_HOT_SPEC_];
      #else
      ThermalWeight=0.0;
      HotWeight=0.0;
      exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
      #endif

      #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      ThermalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
      HotTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_N2_HOT_SPEC_];
      #else
      ThermalTimeStep=0.0;
      HotTimeStep=0.0;
      exit(__LINE__,__FILE__,"Error: the time step node is not defined");
      #endif

      nHotParticles=ThermalWeight/ThermalTimeStep*HotTimeStep/ThermalWeight;

      if ((nHotParticles<1.0)&&(_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_)) {
        //change the species ID of the particle and update the weoght correction factor
        PIC::ParticleBuffer::SetI(_N2_HOT_SPEC_,ptr);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(nHotParticles,ptr);

        //add the paricle to the new particle list
        PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
        PIC::ParticleBuffer::SetPrev(-1,ptr);

        if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
        FirstParticleCell=ptr;

      } else {
        //generate new "hot" paritlces and remove the original "thernal" particle
        int nGeneratedParticles,np;
        long int newptr;

        nGeneratedParticles=(int)nHotParticles;
        if (rnd()<nHotParticles-nGeneratedParticles) nGeneratedParticles++;

        for (np=0;np<nGeneratedParticles;np++) {
          PIC::ParticleBuffer::GetNewParticle(FirstParticleCell);
          PIC::ParticleBuffer::CloneParticle(newptr,ptr);

          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newptr);
          PIC::ParticleBuffer::SetI(_N2_HOT_SPEC_,ptr);
        }

        //remove the original thermal particle
        PIC::ParticleBuffer::DeleteParticle(ptr);
      }
    }
    else {
      //no transformation with the particle is parformed -> just add to the list of the new particles
      PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
      PIC::ParticleBuffer::SetPrev(-1,ptr);

      if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
      FirstParticleCell=ptr;
    }
  }
  else {
    //particle is not _N2_SPEC_ -> all it the list of the new particles
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
}
