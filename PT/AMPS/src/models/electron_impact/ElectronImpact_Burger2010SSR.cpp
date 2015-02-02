//$Id$
//interpolation of the electron impact data from Burger 2010 SSR

#include "ElectronImpact.h"

int ElectronImpact::H2O::Burger2010SSR::ReturnReactionProductList[ElectronImpact::H2O::Burger2010SSR::nMaxReactionProducts];

double ElectronImpact::H2O::Burger2010SSR::GetTotalRateCoefficient(double ElectronTemeprature,double *RateTable) {
  int i0,i1;
  double t,el,c0,c1,TempLOG10,res=0.0;
  int nChannel;

  //determine the element number in the data table and the corresponding interpolation coefficients;
  TempLOG10=log10(ElectronTemeprature);
  el=(TempLOG10-minRateCoefficientTableElectronTemeprature_LOG10)/dRateCoefficientTableElectronTemeprature_LOG10;

  if (el<0.0) {
    // take the data from the first elecemnt of the data array
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E6*pow(10.0,RateCoefficientTable[nChannel][0]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }
  else if (el>=RateCoefficientTableLength) {
    //take the last element of the data array
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E6*pow(10.0,RateCoefficientTable[nChannel][RateCoefficientTableLength-1]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }
  else {
    i0=(int)el;
    c0=i0+1-el;  //the interpolation coefficient for the point i0

    i1=i0+1,c1=1.0-c0;
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E6*pow(10.0,c0*RateCoefficientTable[nChannel][i0]+c1*RateCoefficientTable[nChannel][i1]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }

  return res;
}

double ElectronImpact::H2O::Burger2010SSR::GetTotalRateCoefficient(double ElectronTemeprature) {
  double res,RateTable[nReactionChannels];

  return GetTotalRateCoefficient(ElectronTemeprature,RateTable);
}

void ElectronImpact::H2O::Burger2010SSR::GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
  double RateTable[nReactionChannels],TotalRate,summ;
  int nChannel,i;

  TotalRate=GetTotalRateCoefficient(ElectronTemeprature,RateTable);

  //determine the channel of the impact reaction
  for (TotalRate*=rnd(),nChannel=0,summ=0.0;nChannel<nReactionChannels;nChannel++) {
    summ+=RateTable[nChannel];
    if (summ>TotalRate) break;
  }

  if (nChannel==nReactionChannels) nChannel-=1; //keep the reaction channel index within the range

  //prepare the list of the reaction products
  for (nReactionProducts=0,i=0;i<nMaxReactionProducts;i++) if (ReactionChannelProducts[nChannel][i]>=0) ReturnReactionProductList[nReactionProducts++]=ReactionChannelProducts[nChannel][i];

  //return the list of the reaction products
  ReactionProductsList=ReturnReactionProductList;
}

void ElectronImpact::H2O::Burger2010SSR::Print(const char* fname,const char* OutputDirectory) {
  FILE *fout;
  int nChannel,nprod,i;
  char str[500];

  sprintf(str,"%s/%s",OutputDirectory,fname);
  fout=fopen(str,"w");

  //print the variable list
  fprintf(fout,"VARIABLES=\"Te[eV],\"");

  for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
    sprintf(str,"Channel=%i(",nChannel);

    for (nprod=0,i=0;i<nMaxReactionProducts;i++) if (ReactionChannelProducts[nChannel][i]>=0) {
      if (nprod!=0) sprintf(str,"%s+",str);
      sprintf(str,"%s%s",str,PIC::MolecularData::GetChemSymbol(ReactionChannelProducts[nChannel][i]));
    }

    if (nChannel!=nReactionChannels-1) sprintf(str,"%s)[m^3s^{-1}],",str);
    else sprintf(str,"%s)[m^3s^{-1}]",str);

    fprintf(fout,"\"%s\"",str);
  }

  fprintf(fout,"\n");

  //print the data
  for (int i=0;i<RateCoefficientTableLength;i++) {
    fprintf(fout,"%e",pow(10.0,minRateCoefficientTableElectronTemeprature_LOG10+i*dRateCoefficientTableElectronTemeprature_LOG10));

    for (nChannel=0;nChannel<nReactionChannels;nChannel++) fprintf(fout,", %e",1.0E6*pow(10.0,RateCoefficientTable[i][nChannel]));
    fprintf(fout,"\n");
  }

  //close the file
  fclose(fout);
}






