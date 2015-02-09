//$Id$
//interpolation of the electron impact data from Burger 2010 SSR

#include "ElectronImpact.h"


int ElectronImpact::H2O::Burger2010SSR::General::nReactionChannels=0;
int ElectronImpact::H2O::Burger2010SSR::General::nMaxReactionProducts=0;
int *ElectronImpact::H2O::Burger2010SSR::General::ReactionChannelProducts=NULL;
double *ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTable=NULL;
char *ElectronImpact::H2O::Burger2010SSR::General::ReactionChannelProductsString=NULL;

//-------------------------------------   Burger 2010 SSR:  H2O ---------------------------------
int ElectronImpact::H2O::Burger2010SSR::ReturnReactionProductList[ElectronImpact::H2O::Burger2010SSR::nMaxReactionProducts];

int ElectronImpact::H2O::Burger2010SSR::ReactionChannelProducts[ElectronImpact::H2O::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::nMaxReactionProducts]={
    {_O_PLUS_SPEC_,_H2_SPEC_,_ELECTRON_SPEC_},
    {_H_PLUS_SPEC_,_OH_SPEC_,_ELECTRON_SPEC_},
    {_OH_PLUS_SPEC_,_H_SPEC_,_ELECTRON_SPEC_},
    {_H2O_PLUS_SPEC_,_ELECTRON_SPEC_,-1},
    {_OH_SPEC_,_H_SPEC_,-1}
};

char ElectronImpact::H2O::Burger2010SSR::ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_]={
    "H<sub>2</sub>O+e->O<sup>+</sup>+H<sub>2</sub>+2e",
    "H<sub>2</sub>O+e->H<sup>+</sup>+OH+2e",
    "H<sub>2</sub>O+e->OH<sup>+</sup>+H+2e",
    "H<sub>2</sub>O+e->H<sub>2</sub>O<sup>+</sup>+2e",
    "H<sub>2</sub>O+e->OH+H+e"
};


double ElectronImpact::H2O::Burger2010SSR::RateCoefficientTable[ElectronImpact::H2O::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTableLength]={ //log10 (rate coefficient [cm^3 s{-1}]); the electron temperature grid is uniformly distributed in the logarithmic space
    {-11.98831, -11.98831, -11.98831, -11.98831, -11.98038, -11.98038, -11.98038, -11.98038, -11.98095, -11.98095, -11.98095, -11.98095, -11.98095, -11.98095, -11.98147, -11.97406, -11.76574, -11.53353, -11.30455, -11.09776, -10.90475, -10.71967, -10.55521, -10.39128, -10.23633, -10.10200, -9.96663, -9.83600, -9.73128, -9.60541, -9.51289, -9.39443, -9.30558, -9.21725, -9.12312, -9.05650, -9.00152, -8.92640, -8.86349, -8.80898, -8.77040, -8.70430, -8.66781, -8.61759, -8.58423, -8.54722, -8.51709, -8.48697, -8.46154, -8.44938, -8.44040},
    {-11.97729, -11.97729, -11.97729, -11.97729, -11.97729, -11.97781, -11.97781, -11.97781, -11.97781, -11.97358, -11.97358, -11.98204, -11.70811, -11.44368, -11.18614, -10.95345, -10.70594, -10.49919, -10.29715, -10.12106, -9.94021, -9.77257, -9.62185, -9.47906, -9.33998, -9.22257, -9.10097, -8.99995, -8.88254, -8.78522, -8.68848, -8.60756, -8.51395, -8.44733, -8.37278, -8.30507, -8.23370, -8.18291, -8.12897, -8.07395, -8.02796, -7.98937, -7.93915, -7.91002, -7.86616, -7.83651, -7.81062, -7.78472, -7.75616, -7.73074, -7.72228},
    {-11.97619, -11.97619, -11.98413, -11.97196, -11.97196, -11.97196, -11.97196, -11.97196, -11.97672, -11.98042, -11.67266, -11.37335, -11.06611, -10.82178, -10.58064, -10.33949, -10.15018, -9.95608, -9.77152, -9.59861, -9.44789, -9.30829, -9.17819, -9.06025, -8.94602, -8.83602, -8.74298, -8.64249, -8.55415, -8.48226, -8.40291, -8.34000, -8.27338, -8.21835, -8.15544, -8.10465, -8.06288, -8.01584, -7.96979, -7.92855, -7.89049, -7.86089, -7.82701, -7.80639, -7.77308, -7.75616, -7.73502, -7.71705, -7.70066, -7.68902, -7.65143},
    {-11.97615, -11.97615, -11.97615, -11.97562, -11.97562, -11.68848, -11.33363, -10.99997, -10.72604, -10.43889, -10.21623, -9.97775, -9.76620, -9.56948, -9.39657, -9.22043, -9.07769, -8.93010, -8.80742, -8.67470, -8.56470, -8.46420, -8.36162, -8.26539, -8.18077, -8.10094, -8.02529, -7.94917, -7.88992, -7.83547, -7.76880, -7.71858, -7.67206, -7.62601, -7.57104, -7.52927, -7.49910, -7.45733, -7.42454, -7.38701, -7.35688, -7.33146, -7.30713, -7.28280, -7.25268, -7.23576, -7.21884, -7.20616, -7.18444, -7.16753, -7.15907},
    {-11.66634, -11.36808, -11.07614, -10.79693, -10.55260, -10.31145, -10.11897, -9.91274, -9.72819, -9.56050, -9.40346, -9.25379, -9.09883, -8.96508, -8.84239, -8.73291, -8.61127, -8.51395, -8.41298, -8.31196, -8.23631, -8.14907, -8.07291, -7.98462, -7.91743, -7.85134, -7.78419, -7.72176, -7.65410, -7.60387, -7.55255, -7.50233, -7.46422, -7.40606, -7.35902, -7.32942, -7.29554, -7.25377, -7.21200, -7.19086, -7.16491, -7.13213, -7.10675, -7.08138, -7.07292, -7.05601, -7.03481, -7.02213, -7.02213, -7.02213, -7.01367}
};

//-------------------------------------   Burger 2010 SSR:  O2 ---------------------------------
int ElectronImpact::O2::Burger2010SSR::ReturnReactionProductList[ElectronImpact::O2::Burger2010SSR::nMaxReactionProducts];


int ElectronImpact::O2::Burger2010SSR::ReactionChannelProducts[ElectronImpact::O2::Burger2010SSR::nReactionChannels][ElectronImpact::O2::Burger2010SSR::nMaxReactionProducts]={
    {_O_PLUS_SPEC_,_O_SPEC_,_ELECTRON_SPEC_},
    {_O2_PLUS_SPEC_,_ELECTRON_SPEC_,-1},
    {_O_SPEC_,_O_SPEC_,-1}};

char ElectronImpact::O2::Burger2010SSR::ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_]={
    "O<sub>2</sub>+e->O<sup>+</sup>+O+2e",
    "O<sub>2</sub>+e->O<sub>2</sub><sup>+</sup>+2e",
    "O<sub>2</sub>+e->O+O+e"
};

double ElectronImpact::O2::Burger2010SSR::RateCoefficientTable[ElectronImpact::O2::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTableLength]={
  {-11.00049, -11.00341, -11.00341, -11.00341, -10.99966, -10.99966, -10.99966, -10.99966, -11.00007, -11.00049, -11.00716, -10.99715, -11.00498, -10.99715, -10.81074, -10.54762, -10.30365, -10.08055, -9.89038, -9.70063, -9.50549, -9.35118, -9.19227, -9.05258, -8.90622, -8.79029, -8.66431, -8.56801, -8.44541, -8.35241, -8.26315, -8.17974, -8.09094, -8.01502, -7.93787, -7.86739, -7.81737, -7.75187, -7.69934, -7.64594, -7.59925, -7.55545, -7.52292, -7.47706, -7.45121, -7.41868, -7.38196, -7.35903, -7.32943, -7.30900, -7.29859},
  {-10.99869, -10.99869, -10.99869, -10.99869, -10.99869, -10.99869, -10.99869, -10.99869, -10.99869, -10.65089, -10.38818, -10.16838, -9.94528, -9.75178, -9.56454, -9.38858, -9.22551, -9.07207, -8.92567, -8.78890, -8.67672, -8.55288, -8.43735, -8.33768, -8.24468, -8.15209, -8.05662, -7.98321, -7.89980, -7.83307, -7.77012, -7.70005, -7.64710, -7.58704, -7.52742, -7.47526, -7.43607, -7.39354, -7.35348, -7.31388, -7.28176, -7.25508, -7.22923, -7.19584, -7.16624, -7.15208, -7.12206, -7.10202, -7.07826, -7.04866, -7.04866},
  {-10.99539, -10.99580, -11.00288, -11.00247, -11.00247, -10.91363, -10.63965, -10.37319, -10.15590, -9.93906, -9.75306, -9.56623, -9.41695, -9.24679, -9.10751, -8.98116, -8.86399, -8.75806, -8.66214, -8.56370, -8.48033, -8.40816, -8.33101, -8.27012, -8.21302, -8.15085, -8.11125, -8.06456, -8.01866, -7.98823, -7.94529, -7.91902, -7.88980, -7.87646, -7.85686, -7.84019, -7.82018, -7.80017, -7.79058, -7.78098, -7.77431, -7.76761, -7.76761, -7.76427, -7.77053, -7.77806, -7.78012, -7.79638, -7.80643, -7.82351, -7.85270}
};


//-------------------------------------   Burger 2010 SSR:  H2 ---------------------------------
int ElectronImpact::H2::Burger2010SSR::ReturnReactionProductList[ElectronImpact::H2::Burger2010SSR::nMaxReactionProducts];


int ElectronImpact::H2::Burger2010SSR::ReactionChannelProducts[ElectronImpact::H2::Burger2010SSR::nReactionChannels][ElectronImpact::H2::Burger2010SSR::nMaxReactionProducts]={
    {_H_PLUS_SPEC_,_H_SPEC_,_ELECTRON_SPEC_},
    {_H2_PLUS_SPEC_,_ELECTRON_SPEC_,-1},
    {_H_SPEC_,_H_SPEC_,-1}
};

char ElectronImpact::H2::Burger2010SSR::ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_]={
    "H<sub>2</sub>+e->H<sup>+</sup>+H+2e",
    "H<sub>2</sub>+e->H<sub>2</sub><sup>+</sup>+2e",
    "H<sub>2</sub>+e->H+H+e"
};

double ElectronImpact::H2::Burger2010SSR::RateCoefficientTable[ElectronImpact::H2::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTableLength]={
    {-10.99627, -10.99627, -10.99709, -10.99709, -10.99790, -10.99790, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -11.00365, -10.78689, -10.58655, -10.42808, -10.26138, -10.10168, -9.96331, -9.82660, -9.70878, -9.58399, -9.48627, -9.38527, -9.29043, -9.19931, -9.12129, -9.05932, -8.98170, -8.93530, -8.86635, -8.82116, -8.77151, -8.72879, -8.68364, -8.65040, -8.62126, -8.58592, -8.56294, -8.53708, -8.51451, -8.50137, -8.48536, -8.47592, -8.47592},
    {-11.00004, -11.00004, -11.00004, -11.00004, -11.00004, -11.00004, -10.99716, -10.99716, -10.99797, -10.92611, -10.67158, -10.37352, -10.12064, -9.89525, -9.70679, -9.51382, -9.34388, -9.18087, -9.03800, -8.88118, -8.76089, -8.64306, -8.53181, -8.44065, -8.34910, -8.25713, -8.18243, -8.09744, -8.03543, -7.96649, -7.91720, -7.85523, -7.80964, -7.76448, -7.71893, -7.68894, -7.64954, -7.61383, -7.58096, -7.56126, -7.53212, -7.50626, -7.48037, -7.46026, -7.44097, -7.42496, -7.40854, -7.39869, -7.38268, -7.37283, -7.37283},
    {-10.99753, -10.99753, -11.00409, -11.00125, -11.00125, -11.00125, -11.00125, -11.00125, -11.00125, -11.00741, -10.73315, -10.52993, -10.37599, -10.23596, -10.09475, -9.97405, -9.87511, -9.79053, -9.73143, -9.65957, -9.61114, -9.56923, -9.53640, -9.51054, -9.49413, -9.48428, -9.48055, -9.48715, -9.48715, -9.50029, -9.50726, -9.53190, -9.55695, -9.57336, -9.60619, -9.63534, -9.66861, -9.70801, -9.75401, -9.79340, -9.84103, -9.87758, -9.91658, -9.96542, -9.99828, -10.05329, -10.10829, -10.15469, -10.21504, -10.25776, -10.32055}
};

//-------------------------------------   Burger 2010 SSR:  H ---------------------------------
int ElectronImpact::H::Burger2010SSR::ReturnReactionProductList[ElectronImpact::H::Burger2010SSR::nMaxReactionProducts];


int ElectronImpact::H::Burger2010SSR::ReactionChannelProducts[ElectronImpact::H::Burger2010SSR::nReactionChannels][ElectronImpact::H::Burger2010SSR::nMaxReactionProducts]={
    {_H_PLUS_SPEC_,_ELECTRON_SPEC_}
};

char ElectronImpact::H::Burger2010SSR::ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_]={
    "H+e->H<sup>+</sup>+2e"
};


double ElectronImpact::H::Burger2010SSR::RateCoefficientTable[ElectronImpact::H::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTableLength]={
    {-10.00324, -9.99392, -10.00040, -10.00040, -10.00080, -10.00080, -10.00080, -10.00120, -10.00120, -9.99756, -9.99796, -9.99796, -9.99796, -9.78015, -9.56837, -9.40396, -9.25174, -9.10312, -8.97477, -8.85533, -8.73953, -8.63951, -8.54314, -8.45568, -8.36094, -8.30222, -8.23219, -8.16779, -8.11275, -8.04551, -8.00422, -7.95929, -7.91432, -7.88234, -7.84429, -7.80299, -7.77705, -7.75115, -7.73496, -7.69651, -7.67464, -7.65561, -7.63939, -7.62683, -7.60457, -7.59526, -7.57947, -7.57947, -7.56975, -7.56004, -7.56004}
};

//-------------------------------------   Burger 2010 SSR:  O ---------------------------------
int ElectronImpact::O::Burger2010SSR::ReturnReactionProductList[ElectronImpact::O::Burger2010SSR::nMaxReactionProducts];


int ElectronImpact::O::Burger2010SSR::ReactionChannelProducts[ElectronImpact::O::Burger2010SSR::nReactionChannels][ElectronImpact::O::Burger2010SSR::nMaxReactionProducts]={
    {_O_PLUS_SPEC_,_ELECTRON_SPEC_}
};

char ElectronImpact::O::Burger2010SSR::ReactionChannelProductsString[nReactionChannels][_MAX_STRING_LENGTH_PIC_]={
    "O+e->O<sup>+</sup>+2e"
};


double ElectronImpact::O::Burger2010SSR::RateCoefficientTable[ElectronImpact::O::Burger2010SSR::nReactionChannels][ElectronImpact::H2O::Burger2010SSR::General::RateCoefficientTableLength]={
    {-10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -10.00124, -9.82468, -9.64168, -9.45420, -9.28940, -9.11775, -8.95698, -8.81607, -8.70027, -8.57475, -8.47230, -8.37553, -8.27552, -8.19213, -8.10547, -8.03420, -7.94554, -7.88805, -7.81758, -7.75362, -7.70298, -7.64186, -7.58718, -7.53901, -7.49080, -7.45519, -7.41997, -7.38068, -7.34262, -7.32279, -7.28150, -7.25232, -7.23653, -7.21143, -7.18877, -7.17865, -7.14991, -7.14016, -7.13692}
};


//Calcualtion of the rate corfficient
double ElectronImpact::H2O::Burger2010SSR::General::GetTotalRateCoefficient(double ElectronTemeprature,double *RateTable) {
  int i0,i1;
  double t,el,c0,c1,TempLOG10,res=0.0;
  int nChannel;

  //determine the element number in the data table and the corresponding interpolation coefficients;
  TempLOG10=log10(ElectronTemeprature);
  el=(TempLOG10-minRateCoefficientTableElectronTemeprature_LOG10)/dRateCoefficientTableElectronTemeprature_LOG10;

  if (el<0.0) {
    // take the data from the first elecemnt of the data array
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E-6*pow(10.0,RateCoefficientTable[0+nChannel*RateCoefficientTableLength]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }
  else if (el>=RateCoefficientTableLength) {
    //take the last element of the data array
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E-6*pow(10.0,RateCoefficientTable[RateCoefficientTableLength-1+nChannel*RateCoefficientTableLength]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }
  else {
    i0=(int)el;
    c0=i0+1-el;  //the interpolation coefficient for the point i0

    i1=i0+1,c1=1.0-c0;
    for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
      t=1.0E-6*pow(10.0,c0*RateCoefficientTable[i0+nChannel*RateCoefficientTableLength]+c1*RateCoefficientTable[i1+nChannel*RateCoefficientTableLength]);

      RateTable[nChannel]=t;
      res+=t;
    }
  }

  return res;
}

double ElectronImpact::H2O::Burger2010SSR::General::GetTotalRateCoefficient(double ElectronTemeprature) {
  double res,RateTable[nReactionChannels];

  return GetTotalRateCoefficient(ElectronTemeprature,RateTable);
}

void ElectronImpact::H2O::Burger2010SSR::General::GenerateReactionChannel(double ElectronTemeprature,int* &ReactionProductsList,int &nReactionProducts) {
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
  for (nReactionProducts=0,i=0;i<nMaxReactionProducts;i++) if (ReactionChannelProducts[i+nChannel*nMaxReactionProducts]>=0) ReturnReactionProductList[nReactionProducts++]=ReactionChannelProducts[i+nChannel*nMaxReactionProducts];

  //return the list of the reaction products
  ReactionProductsList=ReturnReactionProductList;
}


void ElectronImpact::H2O::Burger2010SSR::General::Print(const char* fname,const char* OutputDirectory) {
  FILE *fout;
  int nChannel,nprod,i;
  char str[500];

  sprintf(str,"%s/%s",OutputDirectory,fname);
  fout=fopen(str,"w");

  //print the variable list
  fprintf(fout,"VARIABLES=\"Te[eV]\"");

  for (nChannel=0;nChannel<nReactionChannels;nChannel++) {
    sprintf(str,"Channel=%i(%s)[m^3s^{-1}]",nChannel,ReactionChannelProductsString+nChannel*_MAX_STRING_LENGTH_PIC_);

    fprintf(fout,", \"%s\"",str);
  }

  fprintf(fout,"\n");

  //print the data
  for (int i=0;i<RateCoefficientTableLength;i++) {
    fprintf(fout,"%e",pow(10.0,minRateCoefficientTableElectronTemeprature_LOG10+i*dRateCoefficientTableElectronTemeprature_LOG10));

    for (nChannel=0;nChannel<nReactionChannels;nChannel++) fprintf(fout,", %e",1.0E-6*pow(10.0,RateCoefficientTable[i+nChannel*RateCoefficientTableLength]));
    fprintf(fout,"\n");
  }

  //close the file
  fclose(fout);
}






