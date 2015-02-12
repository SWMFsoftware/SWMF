//$Id$
//data se tof the photolytic reaction model (Huebner 1992, ASS)


#include "PhotolyticReactions.h"



//-----------------------------------   H2O ----------------------------------------------
int PhotolyticReactions::H2O::Huebner1992ASS::ReactionProducts[nReactionChannels][nMaxReactionProducts]={
    {_OH_SPEC_,_H_SPEC_,-1},
    {_H2_SPEC_,_O_SPEC_,-1},
    {_H_SPEC_,_H_SPEC_,_O_SPEC_},
    {_H2O_PLUS_SPEC_,_ELECTRON_SPEC_,-1},
    {_H_SPEC_,_OH_PLUS_SPEC_,_ELECTRON_SPEC_},
    {_H2_SPEC_,_O_PLUS_SPEC_,_ELECTRON_SPEC_},
    {_OH_SPEC_,_H_PLUS_SPEC_,_ELECTRON_SPEC_}
};

double PhotolyticReactions::H2O::Huebner1992ASS::EmissionWaveLength[nReactionChannels]={
    -1.0,-1.0,-1.0,-1.0,684.4*_A_,684.8*_A_,662.3*_A_
};

double PhotolyticReactions::H2O::Huebner1992ASS::ReactionRateTable_QuietSun[nReactionChannels]={
  1.03E-5,5.97E-7,7.55E-7,3.31E-7,5.54E-8,5.85E-9,1.31E-8
};

double PhotolyticReactions::H2O::Huebner1992ASS::ReactionRateTable_ActiveSun[nReactionChannels]={
    1.76E-5,1.48E-6,1.91E-6,8.28E-7,1.51E-7,2.21E-8,4.07E-8
};

double *PhotolyticReactions::H2O::Huebner1992ASS::ReactionRateTable=NULL;
double PhotolyticReactions::H2O::Huebner1992ASS::TotalReactionRate=0.0;
double *PhotolyticReactions::H2O::Huebner1992ASS::ExcessEnergyTable=NULL;
int PhotolyticReactions::H2O::Huebner1992ASS::ReturnReactionProductList[nMaxReactionProducts];
double PhotolyticReactions::H2O::Huebner1992ASS::ReturnReactionProductVelocity[3*nMaxReactionProducts];

void PhotolyticReactions::H2O::Huebner1992ASS::Init() {
  //init the tables
  if (_PHOTOLYTIC_REACTIONS__SUN_ == _PHOTOLYTIC_REACTIONS__SUN__ACTIVE_) ReactionRateTable=ReactionRateTable_ActiveSun;
  else if  (_PHOTOLYTIC_REACTIONS__SUN_ == _PHOTOLYTIC_REACTIONS__SUN__QUIET_) ReactionRateTable=ReactionRateTable_QuietSun;
  else exit(__LINE__,__FILE__,"Error: the option is unknown");

  //calculate the total rate
  for (int nChannel=0;nChannel<nReactionChannels;nChannel++) TotalReactionRate+=ReactionRateTable[nChannel];
}


//-----------------------------------   O2 ----------------------------------------------
int PhotolyticReactions::O2::Huebner1992ASS::ReactionProducts[nReactionChannels][nMaxReactionProducts]={
    {_O_SPEC_,_O_SPEC_,-1},
    {_O_SPEC_,_O_SPEC_,-1},
    {_O_SPEC_,_O_SPEC_,-1},
    {_O2_PLUS_SPEC_,_ELECTRON_SPEC_,-1},
    {_O_SPEC_,_O_PLUS_SPEC_,_ELECTRON_SPEC_}
};

double PhotolyticReactions::O2::Huebner1992ASS::ReactionRateTable_QuietSun[nReactionChannels]={
  1.45E-7,4.05E-6,3.90E-8,4.46E-7,1.10E-7
};

double PhotolyticReactions::O2::Huebner1992ASS::ReactionRateTable_ActiveSun[nReactionChannels]={
  2.15E-7,6.47E-6,9.35E-8,1.18E-6,3.47E-7
};

double PhotolyticReactions::O2::Huebner1992ASS::EmissionWaveLength[nReactionChannels]={
    2423.7*_A_,1759.0*_A_,923.0*_A_,-1.0,-1.0
};

double PhotolyticReactions::O2::Huebner1992ASS::ExcessEnergyTable_QuietSun[nReactionChannels]={
    4.39*eV2J,1.33*eV2J,0.74*eV2J,15.9*eV2J,23.8*eV2J
};

double PhotolyticReactions::O2::Huebner1992ASS::ExcessEnergyTable_ActiveSun[nReactionChannels]={
    5.85*eV2J,1.55*eV2J,0.74*eV2J,19.3*eV2J,27.3*eV2J
};


double *PhotolyticReactions::O2::Huebner1992ASS::ReactionRateTable=NULL;
double *PhotolyticReactions::O2::Huebner1992ASS::ExcessEnergyTable=NULL;
double PhotolyticReactions::O2::Huebner1992ASS::TotalReactionRate=0.0;
int PhotolyticReactions::O2::Huebner1992ASS::ReturnReactionProductList[nMaxReactionProducts];
double PhotolyticReactions::O2::Huebner1992ASS::ReturnReactionProductVelocity[3*nMaxReactionProducts];

void PhotolyticReactions::O2::Huebner1992ASS::Init() {
  //init the tables
  if (_PHOTOLYTIC_REACTIONS__SUN_ == _PHOTOLYTIC_REACTIONS__SUN__ACTIVE_) {
    ReactionRateTable=ReactionRateTable_ActiveSun,ExcessEnergyTable=ExcessEnergyTable_ActiveSun;
  }
  else if (_PHOTOLYTIC_REACTIONS__SUN_ == _PHOTOLYTIC_REACTIONS__SUN__QUIET_) {
    ReactionRateTable=ReactionRateTable_QuietSun,ExcessEnergyTable=ExcessEnergyTable_QuietSun;
  }
  else exit(__LINE__,__FILE__,"Error: the option is unknown");

  //calculate the total rate
  for (int nChannel=0;nChannel<nReactionChannels;nChannel++) TotalReactionRate+=ReactionRateTable[nChannel];
}

//-----------------------------------  Chemical model ----------------------------------------
void PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* ReturnReactionProductTable,double *ReturnReactionProductVelocityTable,
    double *ReactionRateTable, int nReactionChannels,int* TotalReactionProductTable,int nMaxReactionProducts,
    double TotalReactionRate,double *ExcessEnergyTable) {
  int i;
  double summ=0.0;

  //1. Determine the reaction channel
  for (TotalReactionRate*=rnd(),i=0,summ=0.0;i<nReactionChannels;i++) {
    summ+=ReactionRateTable[i];
    if (summ>TotalReactionRate) break;
  }

  if (i==nReactionChannels) i=-1;
  ReactionChannel=i;

  //2. Init the list of the reaction products
  int t,iElectronProduct=-1;

  for (i=0,nReactionProducts=0;i<nMaxReactionProducts;i++) {
    if ((t=TotalReactionProductTable[i+ReactionChannel*nMaxReactionProducts])>=0) {
      ReturnReactionProductTable[nReactionProducts++]=t;

      if (t==_ELECTRON_SPEC_) {
        if (iElectronProduct!=-1) {
          exit(__LINE__,__FILE__,"Error: two electrons in the product list is not allowed");
        }

        iElectronProduct=nReactionProducts-1;
      }

    }
  }

  //3. Distribute velocity of the reaction components
  if ( ((iElectronProduct==-1)&&(nReactionProducts>2)) && ((iElectronProduct!=-1)&&(nReactionProducts>3)) ) {
    exit(__LINE__,__FILE__,"Too many reaction products. The models is not implemented for this mode");
  }

  int i0Product=0,i1Product=1; //the location of the products in the reaction product vector
  double TotalEnergy=0.0;
  double HeavyProductsCenterMassVelocity[3]={0.0,0.0,0.0};
  double l[3],ll=0.0; //the direction of the relative velocity of the "heavy" part of the products and the electron
  double i0ProductMass,i1ProductMass,MassHeavyProducts;

  if (ExcessEnergyTable!=NULL) TotalEnergy=ExcessEnergyTable[ReactionChannel];

  //determine the position of the heavy reaction products in the reacion product list
  if (iElectronProduct!=-1) {
    switch (iElectronProduct) {
    case 0:
      i0Product=2;
      break;
    case 1:
      i1Product=2;
      break;
    default:
      exit(__LINE__,__FILE__,"Too many species. The model is not implementes for a case where there are more that 2 products + 1 electron");
    }
  }

  //get the mass of the heavy reaction products
  i0ProductMass=PIC::MolecularData::GetMass(ReturnReactionProductTable[i0Product]);
  i1ProductMass=PIC::MolecularData::GetMass(ReturnReactionProductTable[i1Product]);
  MassHeavyProducts=i0ProductMass+i1ProductMass;

  //distribution the electron velocity (if used)
  if (iElectronProduct!=1) {
    //generate velocity of the electron and substract it from the total excsess energy
    double VelocityElectron;

    for (ll=0.0,i=0;i<3;i++) {
      l[i]=sqrt(-2.0*log(rnd()))*cos(2.0*Pi*rnd());
      ll+=pow(l[i],2);
    }

    ll=sqrt(ll);
    l[0]/=ll,l[1]/=ll,l[2]/=ll;

    VelocityElectron=sqrt(2.0*TotalEnergy/(_MASS_(_ELECTRON_)*(1.0+_MASS_(_ELECTRON_)/MassHeavyProducts)));

    for (i=0;i<3;i++) {
      ReturnReactionProductVelocityTable[i+3*iElectronProduct]=l[i]*VelocityElectron;
      HeavyProductsCenterMassVelocity[i]=-l[i]*VelocityElectron*_MASS_(_ELECTRON_)/MassHeavyProducts;
    }

    TotalEnergy-=_MASS_(_ELECTRON_)/2.0*pow(VelocityElectron,2);
    if (TotalEnergy<0.0) TotalEnergy=0.0;
  }

  //generate velocity of the heavy reaction products
  double i0Velocity,i1Velocity;

  for (ll=0.0,i=0;i<3;i++) {
    l[i]=sqrt(-2.0*log(rnd()))*cos(2.0*Pi*rnd());
    ll+=pow(l[i],2);
  }

  ll=sqrt(ll);
  l[0]/=ll,l[1]/=ll,l[2]/=ll;

  i0Velocity=sqrt(2.0*TotalEnergy/(i0ProductMass*(1.0+i0ProductMass/i1ProductMass)));
  i1Velocity=-i0Velocity*i0ProductMass/i1ProductMass;

  for (i=0;i<3;i++) {
    ReturnReactionProductVelocityTable[i+3*i0Product]=l[i]*i0Velocity+HeavyProductsCenterMassVelocity[i];
    ReturnReactionProductVelocityTable[i+3*i1Product]=l[i]*i1Velocity+HeavyProductsCenterMassVelocity[i];
  }
}

double PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(int spec,double *ReactionRateTable, int nReactionChannels, int* TotalReactionProductTable, int nMaxReactionProducts) {
  double Yield=0.0,TotalReactionRate=0.0;
  int ReactionChannel,i;

  for (ReactionChannel=0;ReactionChannel<nReactionChannels;ReactionChannel++) {
    TotalReactionRate+=ReactionRateTable[ReactionChannel];

    for (i=0;i<nMaxReactionProducts;i++) if (TotalReactionProductTable[i+ReactionChannel*nMaxReactionProducts]==spec) Yield+=TotalReactionRate;
  }

  return (TotalReactionRate>0.0) ? Yield/TotalReactionRate : 0.0;
}
