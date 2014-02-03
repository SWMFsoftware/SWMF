//====================================================
//$Id$
//====================================================
//define specific functions for modeling collisions with the backgound atmosphere


#include "MTGCM.h"
#undef _MTGCM_







//collison cross section of the model particles with the background atmosphere

inline double GetCollisionCrossSectionBackgoundAtmosphereParticle1(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData,PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,double TranslationalEnergy,double cr2) {
#if _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ == _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__HARD_SPHERE_
  /* const static double BackgroundSpecieDiameter[2]={1.95E-10,4.0E-10};

  const static double RefDiameter_O=2.52E-10;

  return Pi*pow(0.5*(RefDiameter_O+BackgroundSpecieDiameter[BackgroundSpecieNumber]),2);*/

  return 3.0E-19;
#elif _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ == _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__FORWARD_SCATTERING_

  return ForwardCollisionCrossSection.GetTotalCrossSection(TranslationalEnergy);
#endif
}


inline double GetSigmaCrMax1(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData) {
#if _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ == _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__HARD_SPHERE_
  return 10000.0*3.0E-19;

#elif _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ == _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__FORWARD_SCATTERING_
  return 10000.0*ForwardCollisionCrossSection.GetMaxTotalCrossSection();
#else
  exit(__LINE__,__FILE__,"Error: option is not found");
#endif

}

//distribute the direction of the relative velocity of particles after collision
inline double UserDefinedVelocityScatteringAngle1(double* Vrel,double TranslationalEnergy,int spec,int BackgroundSpecieNumber) {
  return ForwardCollisionCrossSection.GetRandomScatteringAngle(TranslationalEnergy);
}




inline void GenerateBackgoundAtmosphereParticle1(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int idim;
  double *xmin,*xmax,*xMiddle,x[3],v[3],beta;

  static const double GlobalNeutalTemeprature=179.0;

  //generate positions of the background particle in the cell
  xmin=node->xmin;
  xmax=node->xmax;
  xMiddle=cell->GetX();

  x[0]=xMiddle[0]+(xmax[0]-xmin[0])/_BLOCK_CELLS_X_*(rnd()-0.5);
  x[1]=xMiddle[1]+(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_*(rnd()-0.5);
  x[2]=xMiddle[2]+(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_*(rnd()-0.5);

  PIC::ParticleBuffer::SetX(x,BackgroundAtmosphereParticleData);

  //generate velocity vector for a particle representing the bacground atmosphere
  beta=GetBackgroundMolecularMass(BackgroundSpecieNumber)/(2*Kbol*GlobalNeutalTemeprature);

  for (idim=0;idim<3;idim++) v[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())/beta);

  PIC::ParticleBuffer::SetV(v,BackgroundAtmosphereParticleData);
  PIC::ParticleBuffer::SetI(BackgroundSpecieNumber,BackgroundAtmosphereParticleData);
}

inline double GetBackgroundNumberDensity1(int BackgroundSpecieNumber,double *x) {
  double res=0.0;


  //static MTGSM file readers
  static cDataSetMTGCM O,CO2;
  static bool InitFlag=false;

  //init the readers
  if (InitFlag==false) {
    InitFlag=true;

    O.PlanetRadius=3376.2E3;
    O.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    O.ReadDataFile("O.h");

    CO2.PlanetRadius=3376.2E3;
    CO2.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
    CO2.ReadDataFile("CO2.h");
  }


#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
  switch (BackgroundSpecieNumber) {
  case _O_BACKGROUND_SPEC_ :
    res=(O.DataValueDefined(x)==true) ? O.Interpolate(x) : 0.0;
    break;
  case _CO2_BACKGROUND_SPEC_ :
    res=(CO2.DataValueDefined(x)==true) ? CO2.Interpolate(x) : 0.0;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the optino is not found");
  }

  res*=1E6;

#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
  switch (BackgroundSpecieNumber) {
  case _O_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_O(x);
    break;
  case _CO2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_CO2(x);
    break;
  case _N2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_N2(x);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the optino is not found");
  }
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

  return res;
}

inline double GetCellMeanBackgroundNumberDensity1(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double x[3];
  cell->GetX(x);

#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
  return GetBackgroundNumberDensity(BackgroundSpecieNumber,x);
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
  double res;

  switch (BackgroundSpecieNumber) {
  case _O_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_O(x);
    break;
  case _CO2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_CO2(x);
    break;
  case _N2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_N2(x);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the optino is not found");
  }

  return res;
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

}

inline double GetCellMaximumBackgroundNumberDensity1(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return *(BackgroundSpecieNumber+(double*)(maxLocalBackdroundDensityOffset+cell->GetAssociatedDataBufferPointer()));
}

inline double GetCellLocalBackgroundNumberDensity(double x[3],int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
  return GetBackgroundNumberDensity(BackgroundSpecieNumber,x);
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
  double res;

  switch (BackgroundSpecieNumber) {
  case _O_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_O(x);
    break;
  case _CO2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_CO2(x);
    break;
  case _N2_BACKGROUND_SPEC_ :
    res=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetBackgroundDensity_N2(x);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the optino is not found");
  }

  return res;
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif
}

//check if the particle should remain in the simulation
#if _THERMALIZED_PARTICLE_REMOVE_CRITERION_ == _THERMALIZED_PARTICLE_REMOVE_CRITERION__ESCAPE_SPEED_
inline bool KeepBackgroundAtmosphereParticle1(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData) {
  double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0};

  PIC::ParticleBuffer::GetV(v,BackgroundAtmosphereParticleData);
  PIC::ParticleBuffer::GetX(x,BackgroundAtmosphereParticleData);

  return (0.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])>GravityConstant*_MASS_(_TARGET_)/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])) ? true : false;
}
#elif _THERMALIZED_PARTICLE_REMOVE_CRITERION_ == _THERMALIZED_PARTICLE_REMOVE_CRITERION__LOCAL_BACKGROUND_THERMAL_SPEED_
inline bool KeepBackgroundAtmosphereParticle(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData) {

  //static MTGSM file readers
  static cDataSetMTGCM Tn;
  static bool InitFlag=false;

  //init the readers
  if (InitFlag==false) {
    InitFlag=true;

    Tn.PlanetRadius=3376.2E3;
    Tn.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
    Tn.ReadDataFile("Tn.h");
  }

  //only oxigen atoms can be keeped in the system
  if (PIC::ParticleBuffer::GetI(BackgroundAtmosphereParticleData)!=_O_SPEC_) return false;

  static const double massOxigen=26.56E-27;
  double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0},vThermal2;

  PIC::ParticleBuffer::GetV(v,BackgroundAtmosphereParticleData);
  PIC::ParticleBuffer::GetX(x,BackgroundAtmosphereParticleData);

  //the neutral temeprature
  double temp;


#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
  const static double GlobalNeutalTemeprature=179.0;

  temp=(Tn.DataValueDefined(x)==true) ? Tn.Interpolate(x) : GlobalNeutalTemeprature;
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
  temp=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetNeutralTemeprature(x);
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif


  vThermal2=3.0*Kbol*temp/massOxigen;

  return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>4*vThermal2) ? true : false;
}
#endif



