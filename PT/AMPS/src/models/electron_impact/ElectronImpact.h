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

  //the ElectronImpact::Burger2010SSR functions and variables used by particular model
  namespace Burger2010SSR {
    static const double minlog10RateCoefficientTableElectronTemeprature_LOG10=0.0;
    static const double maxlog10RateCoefficientTableElectronTemeprature_LOG10=2.0;
    static const int log10RateCoefficientTableLength=51;
    static const double dlog10RateCoefficientTableElectronTemeprature_LOG10=(maxlog10RateCoefficientTableElectronTemeprature_LOG10-minlog10RateCoefficientTableElectronTemeprature_LOG10)/(log10RateCoefficientTableLength-1);

    void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int* ReturnReactionProductTable,double *ReturnReactionProductVelocityTable,int &nReactionProducts,
        int *ReactionChannelProducts,int nMaxReactionProducts,double *log10RateCoefficientTable,int nReactionChannels);

    double GetTotalRateCoefficient(double *RateCoefficientTable,double ElectronTemeprature,int nReactionChannels,double *log10RateCoefficientTable);
    void Print(const char* fname,const char* OutputDirectory,int nReactionChannels,int nMaxReactionProducts,int *ReactionChannelProducts,double *log10RateCoefficientTable,char *ReactionChannelProductsString);

    //calculate the yield fpr a particular specie
    double GetSpeciesReactionYield(int spec,double ElectronTemeprature,double *log10RateCoefficientTable, int nReactionChannels, int* ReactionChannelProducts, int nMaxReactionProducts);
  }


  //the models of the electron impact of H2O
  namespace H2O {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      static const int nReactionChannels=5;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double log10RateCoefficientTable[nReactionChannels][ElectronImpact::Burger2010SSR::log10RateCoefficientTableLength];

      extern int ReturnReactionProductTable[nMaxReactionProducts];
      extern double ReturnReactionProductVelocityTable[3*nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        double RateCoefficientTable[nReactionChannels];

        return ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(RateCoefficientTable,ElectronTemeprature,nReactionChannels,&log10RateCoefficientTable[0][0]);
      }

      inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        ElectronImpact::Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,ReturnReactionProductTable,ReturnReactionProductVelocityTable,nReactionProducts,
            &ReactionChannelProducts[0][0],nMaxReactionProducts,
            &log10RateCoefficientTable[0][0],nReactionChannels);

        ReactionProductTable=ReturnReactionProductTable;
        ReactionProductVelocityTable=ReturnReactionProductVelocityTable;
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        ElectronImpact::Burger2010SSR::Print(fname,OutputDirectory,nReactionChannels,nMaxReactionProducts,&ReactionChannelProducts[0][0],&log10RateCoefficientTable[0][0],&ReactionChannelProductsString[0][0]);
      }

      inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
        return ElectronImpact::Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature,&log10RateCoefficientTable[0][0],nReactionChannels,&ReactionChannelProducts[0][0],nMaxReactionProducts);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }

    inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
      return Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature);
    }
  }

  //the model of the electron impact of O2
  namespace O2 {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      static const int nReactionChannels=3;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double log10RateCoefficientTable[nReactionChannels][ElectronImpact::Burger2010SSR::log10RateCoefficientTableLength];

      extern int ReturnReactionProductTable[nMaxReactionProducts];
      extern double ReturnReactionProductVelocityTable[3*nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        double RateCoefficientTable[nReactionChannels];

        return ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(RateCoefficientTable,ElectronTemeprature,nReactionChannels,&log10RateCoefficientTable[0][0]);
      }

      inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        ElectronImpact::Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,ReturnReactionProductTable,ReturnReactionProductVelocityTable,nReactionProducts,
            &ReactionChannelProducts[0][0],nMaxReactionProducts,
            &log10RateCoefficientTable[0][0],nReactionChannels);

        ReactionProductTable=ReturnReactionProductTable;
        ReactionProductVelocityTable=ReturnReactionProductVelocityTable;
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        ElectronImpact::Burger2010SSR::Print(fname,OutputDirectory,nReactionChannels,nMaxReactionProducts,&ReactionChannelProducts[0][0],&log10RateCoefficientTable[0][0],&ReactionChannelProductsString[0][0]);
      }

      inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
        return ElectronImpact::Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature,&log10RateCoefficientTable[0][0],nReactionChannels,&ReactionChannelProducts[0][0],nMaxReactionProducts);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }

    inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
      return Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature);
    }
  }


  namespace H2 {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      static const int nReactionChannels=3;
      static const int nMaxReactionProducts=3;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double log10RateCoefficientTable[nReactionChannels][ElectronImpact::Burger2010SSR::log10RateCoefficientTableLength];

      extern int ReturnReactionProductTable[nMaxReactionProducts];
      extern double ReturnReactionProductVelocityTable[3*nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        double RateCoefficientTable[nReactionChannels];

        return ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(RateCoefficientTable,ElectronTemeprature,nReactionChannels,&log10RateCoefficientTable[0][0]);
      }

      inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        ElectronImpact::Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,ReturnReactionProductTable,ReturnReactionProductVelocityTable,nReactionProducts,
            &ReactionChannelProducts[0][0],nMaxReactionProducts,
            &log10RateCoefficientTable[0][0],nReactionChannels);

        ReactionProductTable=ReturnReactionProductTable;
        ReactionProductVelocityTable=ReturnReactionProductVelocityTable;
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        ElectronImpact::Burger2010SSR::Print(fname,OutputDirectory,nReactionChannels,nMaxReactionProducts,&ReactionChannelProducts[0][0],&log10RateCoefficientTable[0][0],&ReactionChannelProductsString[0][0]);
      }

      inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
        return ElectronImpact::Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature,&log10RateCoefficientTable[0][0],nReactionChannels,&ReactionChannelProducts[0][0],nMaxReactionProducts);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }

    inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
      return Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature);
    }
  }

  namespace H {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      static const int nReactionChannels=1;
      static const int nMaxReactionProducts=2;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double log10RateCoefficientTable[nReactionChannels][ElectronImpact::Burger2010SSR::log10RateCoefficientTableLength];

      extern int ReturnReactionProductTable[nMaxReactionProducts];
      extern double ReturnReactionProductVelocityTable[3*nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        double RateCoefficientTable[nReactionChannels];

        return ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(RateCoefficientTable,ElectronTemeprature,nReactionChannels,&log10RateCoefficientTable[0][0]);
      }

      inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        ElectronImpact::Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,ReturnReactionProductTable,ReturnReactionProductVelocityTable,nReactionProducts,
            &ReactionChannelProducts[0][0],nMaxReactionProducts,
            &log10RateCoefficientTable[0][0],nReactionChannels);

        ReactionProductTable=ReturnReactionProductTable;
        ReactionProductVelocityTable=ReturnReactionProductVelocityTable;
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        ElectronImpact::Burger2010SSR::Print(fname,OutputDirectory,nReactionChannels,nMaxReactionProducts,&ReactionChannelProducts[0][0],&log10RateCoefficientTable[0][0],&ReactionChannelProductsString[0][0]);
      }

      inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
        return ElectronImpact::Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature,&log10RateCoefficientTable[0][0],nReactionChannels,&ReactionChannelProducts[0][0],nMaxReactionProducts);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }

    inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
      return Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature);
    }
  }


  namespace O {

    namespace Burger2010SSR { //the model data are digitized from Burger et al., 2014 in Space Science Review
      static const int nReactionChannels=1;
      static const int nMaxReactionProducts=2;

      extern int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts];
      extern char ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_];
      extern double log10RateCoefficientTable[nReactionChannels][ElectronImpact::Burger2010SSR::log10RateCoefficientTableLength];

      extern int ReturnReactionProductTable[nMaxReactionProducts];
      extern double ReturnReactionProductVelocityTable[3*nMaxReactionProducts];

      inline double GetTotalRateCoefficient(double ElectronTemeprature) {
        double RateCoefficientTable[nReactionChannels];

        return ElectronImpact::Burger2010SSR::GetTotalRateCoefficient(RateCoefficientTable,ElectronTemeprature,nReactionChannels,&log10RateCoefficientTable[0][0]);
      }

      inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
        ElectronImpact::Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,ReturnReactionProductTable,ReturnReactionProductVelocityTable,nReactionProducts,
            &ReactionChannelProducts[0][0],nMaxReactionProducts,
            &log10RateCoefficientTable[0][0],nReactionChannels);

        ReactionProductTable=ReturnReactionProductTable;
        ReactionProductVelocityTable=ReturnReactionProductVelocityTable;
      }


      inline void Print(const char* fname,const char* OutputDirectory) {
        ElectronImpact::Burger2010SSR::Print(fname,OutputDirectory,nReactionChannels,nMaxReactionProducts,&ReactionChannelProducts[0][0],&log10RateCoefficientTable[0][0],&ReactionChannelProductsString[0][0]);
      }

      inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
        return ElectronImpact::Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature,&log10RateCoefficientTable[0][0],nReactionChannels,&ReactionChannelProducts[0][0],nMaxReactionProducts);
      }

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline void GenerateReactionProducts(double ElectronTemeprature,int &ReactionChannel,int &nReactionProducts, int* &ReactionProductTable,double* &ReactionProductVelocityTable) {
      Burger2010SSR::GenerateReactionProducts(ElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductTable,ReactionProductVelocityTable);
    }

    inline void Print(const char* fname,const char* OutputDirectory) {
      Burger2010SSR::Print(fname,OutputDirectory);
    }

    inline double GetSpeciesReactionYield(int spec,double ElectronTemeprature) {
      return Burger2010SSR::GetSpeciesReactionYield(spec,ElectronTemeprature);
    }
  }


}






#endif
