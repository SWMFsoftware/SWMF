// ----------------------------------------------- //
//$Id$  
// Header file to read chemical reaction input
// Definition : To be continued

// ----------------------------------------------- //

#ifndef _ChemicalSolver_
#define _ChemicalSolver_

namespace ChemicalSolver {
  // Number of reactions read from the input file 
  const int nReactions = 0;
 
  // Number of sources and products
  const int nProductTotal = 0; 
  const int nSourceTotal = 0;
 
  // The species participate in the reaction
  const int Products[nProductTotal] = {};
  const int Sources[nSourceTotal] = {};

  // The offset for the products and sources
  double ReactionProductOffset[nReactions] = {};
  double ReactionSourceOffset[nReactions] = {};

  // The reaction rate
  double ReactionRate[nReactions] = {};

  // The coefficients of the sources and products
  double CoefProduct[nProductTotal] = {};
  double CoefSource[nSourceTotal] = {};

}
#endif
