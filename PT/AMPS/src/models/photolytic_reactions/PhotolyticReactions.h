/*
 * PhotolyticReactions.h
 *
 *  Created on: Feb 5, 2015
 *      Author: vtenishe
 */

//$Id$
//the data base of the photolytic reaction parameters for different species

#ifndef _PHOTOLYTICREACTIONS_
#define _PHOTOLYTICREACTIONS_

#include "pic.h"
#include "PhotolyticReactions.dfn"

namespace PhotolyticReactions {

  //the general functions and variables used by particular model
  namespace Huebner1992ASS {
    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* ReturnReactionProductTable,
        double *ReactionRateTable, int nReactionChannels,int* TotalReactionProductTable,int nMaxReactionProducts,double TotalReactionRate) {
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
      for (i=0,nReactionProducts=0;i<nMaxReactionProducts;i++) {
        if (TotalReactionProductTable[i+ReactionChannel*nMaxReactionProducts]>=0) {
          ReturnReactionProductTable[nReactionProducts++]=TotalReactionProductTable[i+ReactionChannel*nMaxReactionProducts];
        }
      }
    }
  }

  //reaction rates for H2O
  namespace H2O {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=3;
      const int nReactionChannels=7;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate);

        ReactionProductTable=ReturnReactionProductList;
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable);
    }

  }

  //reaction rates for O2
  namespace O2 {
    namespace Huebner1992ASS {
      using namespace PhotolyticReactions::H2O::Huebner1992ASS;

      const int nMaxReactionProducts=3;
      const int nReactionChannels=5;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate);

        ReactionProductTable=ReturnReactionProductList;
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable);
    }
  }

  inline void Init() {
    H2O::Huebner1992ASS::Init();
    O2::Huebner1992ASS::Init();
  }
}



#endif /* PHOTOLYTICREACTIONS_H_ */
