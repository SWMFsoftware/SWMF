//source file for general functions using physical model

#include "Sputtering.h"

// initialization of sputtering model
void Sputtering::Init(){
  
  // number and list of target species  to be used for initialization
  int nSpecTargetInit;
  const int *SpecTargetInit;
  
  // determine surface type
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_
  
  // icy surface, determine ice model
#if   _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__RUBIN_
  // use Martin Rubin's model for ice sputtering
  SpecTargetInit  = Sputtering::Ice::Rubin::SpecTarget;
  nSpecTargetInit = Sputtering::Ice::Rubin::nSpecTarget;
#elif _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__TEOLIS_
  // use Benjamin Teolis's model for ice sputtering
  SpecTargetInit  = Sputtering::Ice::Teolis::SpecTarget;
  nSpecTargetInit = Sputtering::Ice::Teolis::nSpecTarget;
#else
  exit(__LINE__,__FILE__,"Error: ice sputtering model is not recognized");
#endif
  // icy surface done
  
#else
  exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif
  
  int count = 0; // count species actually used in the simulation     
  for(int iSpecTarget = 0; iSpecTarget < nSpecTargetInit; iSpecTarget++)
    if(SpecTargetInit[iSpecTarget] >= 0) ++count;
  // set number and list of species in simulation
  Sputtering::nSpecTarget = count;
  Sputtering::SpecTarget  = new int[count];
  Sputtering::Yield       = new double[count];
  for(int iSpecTarget=0, i=0; iSpecTarget < nSpecTargetInit; iSpecTarget++)
    if(SpecTargetInit[iSpecTarget] >= 0){
      Sputtering::Yield[i]        = 0.0;
      Sputtering::SpecTarget[i++] = SpecTargetInit[iSpecTarget];
    } 
}



// 1st Master function: 
  //   input: speed and species of a projectile
  //  output: list of species and vector of corresponding yields,
  //          returns number of species sputtered
int Sputtering::GetSputteringYield(double speed, int spec_projectile,
		       const int* &spec_target, double* &yield){
  // calculate yield for every species in the list SpecTarget
  for(int iSpec=0; iSpec < Sputtering::nSpecTarget; iSpec++){
    
    // determine surface type
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_
    
    // icy surface, determine ice model      
#if   _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__RUBIN_
    //use Martin Rubin's model
    Sputtering::Yield[iSpec]=Sputtering::Ice::Rubin::SputteringYield( speed,spec_projectile,SpecTarget[iSpec]);
#elif _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__TEOLIS_
    //use Benjamin Teolis's model
    Sputtering::Yield[iSpec]=Sputtering::Ice::Teolis::SputteringYield(speed,spec_projectile,SpecTarget[iSpec]);
#else
    exit(__LINE__,__FILE__,"Error: ice sputtering model is not recognized");
#endif
    // icy surface done

#else
    exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif
  }
  
  yield       = Sputtering::Yield;
  spec_target = Sputtering::SpecTarget;
  return Sputtering::nSpecTarget;
}

// 2nd Master function
//   input: species of a sputtered particle
//  output: its speed and statistical weight correction
void Sputtering::GetSputteringSpeed(int spec, double& speed, double& weight_correction){
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_
  Sputtering::Ice::GetSputteringSpeed(spec,speed, weight_correction);
#else
  exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif
}
