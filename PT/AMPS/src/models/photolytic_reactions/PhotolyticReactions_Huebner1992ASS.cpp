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
