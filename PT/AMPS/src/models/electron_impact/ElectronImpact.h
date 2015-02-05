//$Id$
//The electron impact reactions


#ifndef _ELECTRON_IMPACT_
#define _ELECTRON_IMPACT_

#include "pic.h"

//Below is calcualtion of the parameters of the electron impact processes
//Units:
//   1. Electron Temeprature is in eV
//   2. The rate coefficient is in m^3 s^{-1}

namespace ElectronImpact {

  //the models of the electron impact of H2O
  namespace H2O {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review

      namespace General { //the functions and constant common for all species of the Burger2010SSR
        static const double minRateCoefficientTableElectronTemeprature_LOG10=0.0;
        static const double maxRateCoefficientTableElectronTemeprature_LOG10=2.0;
        static const int RateCoefficientTableLength=51;
        static const double dRateCoefficientTableElectronTemeprature_LOG10=(maxRateCoefficientTableElectronTemeprature_LOG10-minRateCoefficientTableElectronTemeprature_LOG10)/(RateCoefficientTableLength-1);

        //input parameters
        extern int nReactionChannels,nMaxReactionProducts;
        extern int *ReactionChannelProducts;
        extern double *RateCoefficientTable;
        extern char *ReactionChannelProductsString;


        double GetTotalRateCoefficient(double ElectronTemeprature);
        double GetTotalRateCoefficient(double ElectronTemeprature,double *RateTable);
        void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts);

        void Print(const char* fname,const char* OutputDirectory);


      }


      static const int nReactionChannels=5;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double RateCoefficientTable[nReactionChannels][General::RateCoefficientTableLength];
      extern int ReturnReactionProductList[nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        return General::GetTotalRateCoefficient(ElectronTemeprature);
      }

      inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        General::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];
        General::ReactionChannelProductsString=&ReactionChannelProductsString[0][0];

        General::Print(fname,OutputDirectory);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
      Burger2010SSR::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }
  }

  //the model of the electron impact of O2
  namespace O2 {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      using namespace ElectronImpact::H2O::Burger2010SSR;

      static const int nReactionChannels=3;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern double RateCoefficientTable[nReactionChannels][General::RateCoefficientTableLength];
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        return General::GetTotalRateCoefficient(ElectronTemeprature);
      }

      inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        General::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
      }

      inline void Print(const char* fname,const char* OutputDirectory) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];
        General::ReactionChannelProductsString=&ReactionChannelProductsString[0][0];

        General::Print(fname,OutputDirectory);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
      Burger2010SSR::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }
  }


  namespace H2 {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      using namespace ElectronImpact::H2O::Burger2010SSR;

      static const int nReactionChannels=3;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern double RateCoefficientTable[nReactionChannels][General::RateCoefficientTableLength];
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        return General::GetTotalRateCoefficient(ElectronTemeprature);
      }

      inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        General::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
      }

      inline void Print(const char* fname,const char* OutputDirectory) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];
        General::ReactionChannelProductsString=&ReactionChannelProductsString[0][0];

        General::Print(fname,OutputDirectory);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
      Burger2010SSR::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }
  }

  namespace H {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      using namespace ElectronImpact::H2O::Burger2010SSR;

      static const int nReactionChannels=1;
      static const int nMaxReactionProducts=2;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern double RateCoefficientTable[nReactionChannels][General::RateCoefficientTableLength];
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        return General::GetTotalRateCoefficient(ElectronTemeprature);
      }

      inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        General::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
      }

      inline void Print(const char* fname,const char* OutputDirectory) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];
        General::ReactionChannelProductsString=&ReactionChannelProductsString[0][0];

        General::Print(fname,OutputDirectory);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
      Burger2010SSR::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }
  }


  namespace O {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      using namespace ElectronImpact::H2O::Burger2010SSR;

      static const int nReactionChannels=1;
      static const int nMaxReactionProducts=2;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern double RateCoefficientTable[nReactionChannels][General::RateCoefficientTableLength];
      extern int ReturnReactionProductList[nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        return General::GetTotalRateCoefficient(ElectronTemeprature);
      }

      inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];

        General::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
      }

      inline void Print(const char* fname,const char* OutputDirectory) {
        //set up the species specific parameters
        General::nReactionChannels=nReactionChannels;
        General::nMaxReactionProducts=nMaxReactionProducts;
        General::ReactionChannelProducts=&ReactionChannelProducts[0][0];
        General::RateCoefficientTable=&RateCoefficientTable[0][0];
        General::ReactionChannelProductsString=&ReactionChannelProductsString[0][0];

        General::Print(fname,OutputDirectory);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
      Burger2010SSR::GenerateReactionChannel(ElectronTemeprature,ReactionProductsList,nReactionProducts);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }
  }


}






#endif
