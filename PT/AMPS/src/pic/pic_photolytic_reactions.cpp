//====================================================
//$Id$
//====================================================
//the model of photolytic reactions

#include "pic.h"

double *PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime=NULL;
PIC::ChemicalReactions::PhotolyticReactions::fReactionProcessor *PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable=NULL;
PIC::ChemicalReactions::PhotolyticReactions::fTotalLifeTime *PIC::ChemicalReactions::PhotolyticReactions::TotalLifeTime=NULL;


void PIC::ChemicalReactions::PhotolyticReactions::Init() {

  //only one particle transformation model can be used
  if (_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_) exit(__LINE__,__FILE__,"Error: only one particle transformation model can be used");

  if (TotalLifeTime==NULL) {
    ConstantTotalLifeTime=new double [PIC::nTotalSpecies];
    ReactionProcessorTable=new fReactionProcessor[PIC::nTotalSpecies];
    TotalLifeTime=new fTotalLifeTime [PIC::nTotalSpecies];

    for (int i=0;i<PIC::nTotalSpecies;i++) {
      ConstantTotalLifeTime[i]=-1.0;
      ReactionProcessorTable[i]=NULL;
      TotalLifeTime[i]=TotalLifeTime_default;
    }
  }
}


double PIC::ChemicalReactions::PhotolyticReactions::TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag) {
  ReactionAllowedFlag=(ConstantTotalLifeTime[spec]>0.0) ? true : false;
  return ConstantTotalLifeTime[spec];
}

void PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(fReactionProcessor f,int spec) {
  if (TotalLifeTime==NULL) Init();
  if ((spec<0)||(spec>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: out of range");

  ReactionProcessorTable[spec]=f;
}
void PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(fTotalLifeTime f,int spec) {
  if (TotalLifeTime==NULL) Init();
  if ((spec<0)||(spec>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: out of range");

  TotalLifeTime[spec]=f;
}

