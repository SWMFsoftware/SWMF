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

    //generate the channel of the reaction, determin its products and velocities
    void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* ReturnReactionProductTable,double *ReturnReactionProductVelocityTable,
        double *ReactionRateTable, int nReactionChannels,int* TotalReactionProductTable,int nMaxReactionProducts,
        double TotalReactionRate,double *ExcessEnergyTable);


    //calculate the yield fpr a particular specie
    double GetSpeciesReactionYield(int spec,double *ReactionRateTable, int nReactionChannels, int* TotalReactionProductTable, int nMaxReactionProducts);
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
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];


      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }


    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }

  }

  //reaction rates for O2
  namespace O2 {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=3;
      const int nReactionChannels=5;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }
  }


  namespace H {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=2;
      const int nReactionChannels=1;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }
  }

  namespace H2 {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=3;
      const int nReactionChannels=4;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }
  }

  namespace OH {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=2;
      const int nReactionChannels=4;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }
  }

  namespace O {
    namespace Huebner1992ASS {
      const int nMaxReactionProducts=2;
      const int nReactionChannels=1;
      const double TableHeliocentricDistance=1.0*_AU_;

      extern int ReactionProducts[nReactionChannels][nMaxReactionProducts];
      extern double EmissionWaveLength[nReactionChannels];

      //the table of the rates for each channel
      extern double ReactionRateTable_QuietSun[nReactionChannels],ReactionRateTable_ActiveSun[nReactionChannels];
      extern double ExcessEnergyTable_QuietSun[nReactionChannels],ExcessEnergyTable_ActiveSun[nReactionChannels];
      extern double *ReactionRateTable,*ExcessEnergyTable;
      extern double TotalReactionRate;

      //the buffer to return the list of the reaction products
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern double ReturnReactionProductVelocity[3*nMaxReactionProducts];

      void Init();

      inline double GetTotalReactionRate(double HeliocentricDistance) {
        return TotalReactionRate*pow(TableHeliocentricDistance/HeliocentricDistance,2);
      }

      inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        PhotolyticReactions::Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReturnReactionProductList,ReturnReactionProductVelocity,
            ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts,TotalReactionRate,ExcessEnergyTable);

        ReactionProductTable=ReturnReactionProductList;
        ReactionProductVelocityTable=ReturnReactionProductVelocity;
      }

      inline double GetSpeciesReactionYield(int spec) {
        return PhotolyticReactions::Huebner1992ASS::GetSpeciesReactionYield(spec,ReactionRateTable,nReactionChannels,&ReactionProducts[0][0],nMaxReactionProducts);
      }
    }

    inline double GetTotalReactionRate(double HeliocentricDistance) {
      return Huebner1992ASS::GetTotalReactionRate(HeliocentricDistance);
    }

    inline void GenerateReactionProducts(int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Huebner1992ASS::GenerateReactionProducts(ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline double GetSpeciesReactionYield(int spec) {
      return Huebner1992ASS::GetSpeciesReactionYield(spec);
    }
  }

  inline void Init() {
    H2O::Huebner1992ASS::Init();
    O2::Huebner1992ASS::Init();
    H2::Huebner1992ASS::Init();
    H::Huebner1992ASS::Init();
    OH::Huebner1992ASS::Init();
    O::Huebner1992ASS::Init();
  }
}



#endif /* PHOTOLYTICREACTIONS_H_ */
