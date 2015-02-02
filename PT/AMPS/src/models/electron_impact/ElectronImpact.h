//$Id
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
      static const int nReactionChannels=5;
      static const int nMaxReactionProducts=3;

      static const int ReactionChannelProducts[nReactionChannels][nMaxReactionProducts]={
          {_H_PLUS_SPEC_,_OH_SPEC_,_ELECTRON_SPEC_},{_O_PLUS_SPEC_,_H2_SPEC_,_ELECTRON_SPEC_},{_OH_PLUS_SPEC_,_H_SPEC_,_ELECTRON_SPEC_},
          {_H2O_PLUS_SPEC_,_ELECTRON_SPEC_,-1},{_OH_SPEC_,_H_SPEC_,_ELECTRON_SPEC_}};

      static const double minRateCoefficientTableElectronTemeprature_LOG10=0.0;
      static const double maxRateCoefficientTableElectronTemeprature_LOG10=2.0;
      static const int RateCoefficientTableLength=51;
      static const double dRateCoefficientTableElectronTemeprature_LOG10=(maxRateCoefficientTableElectronTemeprature_LOG10-minRateCoefficientTableElectronTemeprature_LOG10)/(RateCoefficientTableLength-1);

      static const double RateCoefficientTable[RateCoefficientTableLength][nReactionChannels]={ //log10 (rate coefficient [cm^3 s{-1}]); the electron temperature grid is uniformly distributed in the logarithmic space
          {-12.00418,-12.00233,-12.00901,-12.00586,-11.67707},
          {-12.00364,-12.00233,-12.00836,-12.00586,-11.42510},
          {-12.00364,-12.00233,-12.00836,-12.00586,-11.09881},
          {-12.00364,-12.00233,-12.00836,-12.00586,-10.83414},
          {-12.00423,-12.01683,-12.01319,-12.00646,-10.57012},
          {-12.00000,-12.00293,-11.99929,-11.68864,-10.34106},
          {-11.99940,-12.01260,-11.99929,-11.34242,-10.13562},
          {-12.00423,-12.00896,-12.00478,-11.03002,-9.93865},
          {-11.99940,-12.00896,-12.00537,-10.72002,-9.74228},
          {-12.00000,-12.00896,-12.00537,-10.45904,-9.56040},
          {-12.00000,-12.00836,-11.68332,-10.22276,-9.40393},
          {-12.00000,-12.00836,-11.36967,-9.99316,-9.25045},
          {-11.75406,-12.00836,-11.07183,-9.77681,-9.10060},
          {-11.46224,-12.00836,-10.81682,-9.58228,-8.97611},
          {-11.17097,-12.00836,-10.57756,-9.38591,-8.85770},
          {-10.94017,-12.00836,-10.35702,-9.22336,-8.74228},
          {-10.70514,-11.79266,-10.14974,-9.08497,-8.62805},
          {-10.49427,-11.56062,-9.95157,-8.94598,-8.53684},
          {-10.30751,-11.32016,-9.76302,-8.80764,-8.42320},
          {-10.13958,-11.10864,-9.59987,-8.68256,-8.33379},
          {-9.94316,-10.92676,-9.45670,-8.55991,-8.23774},
          {-9.78490,-10.73522,-9.31592,-8.45535,-8.15251},
          {-9.63082,-10.57207,-9.17873,-8.35865,-8.07156},
          {-9.48879,-10.40958,-9.06091,-8.27227,-7.99419},
          {-9.34014,-10.25669,-8.94734,-8.18285,-7.92231},
          {-9.23020,-10.11950,-8.83007,-8.11152,-7.84315},
          {-9.11173,-9.97931,-8.74369,-8.02514,-7.78213},
          {-8.98730,-9.85488,-8.64944,-7.96352,-7.72110},
          {-8.88213,-9.74548,-8.56360,-7.88251,-7.65943},
          {-8.78728,-9.61616,-8.48808,-7.83055,-7.60204},
          {-8.69721,-9.50562,-8.40290,-7.76893,-7.55915},
          {-8.60600,-9.40952,-8.33579,-7.71638,-7.48119},
          {-8.52142,-9.31408,-8.28384,-7.66860,-7.44373},
          {-8.45551,-9.23188,-8.21255,-7.61605,-7.39117},
          {-8.36913,-9.14186,-8.16477,-7.56469,-7.35371},
          {-8.31109,-9.06091,-8.09772,-7.52663,-7.31565},
          {-8.24464,-9.00348,-8.05364,-7.48374,-7.29149},
          {-8.17938,-8.93159,-8.01010,-7.45469,-7.24915},
          {-8.13225,-8.87355,-7.96780,-7.41664,-7.22075},
          {-8.07910,-8.81134,-7.92611,-7.38764,-7.17846},
          {-8.03198,-8.76904,-7.90255,-7.35805,-7.15489},
          {-7.98480,-8.71589,-7.85906,-7.31636,-7.12590},
          {-7.93707,-8.67295,-7.83549,-7.29279,-7.09870},
          {-7.89418,-8.63489,-7.80710,-7.26500,-7.07090},
          {-7.87003,-8.59200,-7.76964,-7.25533,-7.05760},
          {-7.82713,-8.55877,-7.75574,-7.22754,-7.04311},
          {-7.80357,-8.51105,-7.72735,-7.20397,-7.03285},
          {-7.77029,-8.48808,-7.70313,-7.19008,-7.03285},
          {-7.75096,-8.46327,-7.67957,-7.17552,-7.01895},
          {-7.72740,-8.43971,-7.66024,-7.16163,-7.01895},
          {-7.71350,-8.43428,-7.66507,-7.15620,-7.00505}
      };

      extern int ReturnReactionProductList[nMaxReactionProducts];

      double GetTotalRateCoefficient(double ElectronTemeprature);
      double GetTotalRateCoefficient(double ElectronTemeprature,double *RateTable);
      void GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts);
      void Print(const char* fname,const char* OutputDirectory);

    }

    inline double RateCoefficient(double ElectronTemeprature) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature);
    }

    inline double GetTotalRateCoefficient(double ElectronTemeprature,double *RateTable) {
      return Burger2010SSR::GetTotalRateCoefficient(ElectronTemeprature,RateTable);
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
