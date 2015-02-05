// Physical model for surface sputtering
// HOW TO ADD A NEW MODEL
// Two additions are required:
// 1) add a new namespace(s):
//    - namespace SurfaceType (if needed)
//    - namespace ModelName   (within namespace SurfaceType)
//    source file Sputtering.cpp contains yield table for all models
// 2) add a corresponding branch in:
//    - Sputtering::Init,
//    - Sputtering::GetSputteringYield,
//    - Sputtering::SurfaceType::GetSputteringSpeed
//------------------------------------------------------------------------

#ifndef _SPUTTERING_ 
#define _SPUTTERING_

#include "Sputtering.dfn"

#include "pic.h"


namespace Sputtering {
  // number and list of sputtered species, to be set by call of Init()
  static int nSpecTarget;
  static int *SpecTarget;
  // list of Yields corresponding to list of SpecTarget
  double *Yield;

  // icy surface ==============================================================
  namespace Ice {

    // yield model by Martin Rubin --------------------------------------------
    namespace Rubin {

      // sputtered species implemented in the model
      static const int nSpecTarget=1;
      static const int SpecTarget[nSpecTarget] = {_H2O_SPEC_};

      inline double SputteringYield(double speed, 
				    int spec_projectile,
				    int spec_target);
    } 
    // namespace Rubin --------------------------------------------------------

    // yield model by Benjamin Teolis -----------------------------------------
    namespace Teolis {

      // sputtered species implemented in the model
      static const int nSpecTarget = 2;
      static const int SpecTarget[nSpecTarget] = {_O2_SPEC_, _H2O_SPEC_};

      inline double SputteringYield(double speed, 
				       int spec_projectile,
				       int spec_target);
    } 
    // namespace Teolis -------------------------------------------------------

    double GetSputteringSpeed(int spec, double* v);
  }
  // namespace Ice ============================================================


  // initialization of sputtering model
  void Init(){

    // number and list of target species  to be used for initialization
    int nSpecTargetInit;
    const int *SpecTargetInit;
    
    // determine surface type
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_

    // icy surface, determine ice model
#if   _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__RUBIN_
    // use Martin Rubin's model for ice sputtering
    SpecTargetInit  = Ice::Rubin::SpecTarget;
    nSpecTargetInit = Ice::Rubin::nSpecTarget;
#elif _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__TEOLIS_
    // use Benjamin Teolis's model for ice sputtering
    SpecTargetInit  = Ice::Teolis::SpecTarget;
    nSpecTargetInit = Ice::Teolis::nSpecTarget;
#else
    exit(__LINE__,__FILE__,"Error: ice sputtering model is not recognized");
#endif
    // icy surface done

#else
    exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif

    int count = 0; // count species actually used in the simulation     
    for(int iSpecTarget = 0; iSpecTarget < nSpecTargetInit; iSpecTarget++)
      if(SpecTargetInit[iSpecTarget] > 0) ++count;
    if(count == 0)
      std::cout << "WARNING: no species will be sputtered";
    else { // set number and list of species in simulation
      nSpecTarget = count;
      SpecTarget  = new int[count];
      Yield       = new double[count];
      for(int iSpecTarget=0, i=0; iSpecTarget < nSpecTargetInit; iSpecTarget++)
	if(SpecTargetInit[iSpecTarget] > 0){
	  Yield[i] = 0.0;
	  SpecTarget[i++] = SpecTargetInit[iSpecTarget];
	} 
    }
  }
  

  // 1st Master function: 
  //   input: speed and species of a projectile
  //  output: list of species and vector of corresponding yields,
  //          returns number of species sputtered
  int GetSputteringYield(double speed, int spec_projectile,
			 const int* spec_target, double* yield){
    // calculate yield for every species in the list SpecTarget
    for(int iSpec=0; iSpec < nSpecTarget; iSpec++){

      // determine surface type
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_
      
      // icy surface, determine ice model      
#if   _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__RUBIN_
      //use Martin Rubin's model
      yield[iSpec]=Ice::Rubin::SputteringYield( speed,spec_projectile,SpecTarget[iSpec]);
#elif _SPUTTERING__ICE__MODEL_ == _SPUTTERING__ICE__TEOLIS_
      //use Benjamin Teolis's model
      yield[iSpec]=Ice::Teolis::SputteringYield(speed,spec_projectile,SpecTarget[iSpec]);
#else
      exit(__LINE__,__FILE__,"Error: ice sputtering model is not recognized");
#endif
    // icy surface done

#else
    exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif
    }

    spec_target = SpecTarget;
    return nSpecTarget;
  }
  
  // 2nd Master function
  //   input: species of a sputtered particle
  //  output: its velocity vector
  void GetSputteringSpeed(int spec, double* v){
#if _SPUTTERING__SURFACE_ == _SPUTTERING__ICE_
    Ice::GetSputteringSpeed(spec,v);
#else
    exit(__LINE__,__FILE__,"Error: sputtering surface type is not recognized");
#endif
  }

}

#endif//_SPUTTERING_
