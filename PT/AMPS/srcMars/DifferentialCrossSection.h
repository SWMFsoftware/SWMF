//=====================================================
//$Id$
//=====================================================
//the class definition for the differectial cross section model

#ifndef DIFFERENTIALCROSSSECTION_H_
#define DIFFERENTIALCROSSSECTION_H_


#include <vector>
#include <algorithm>
#include "specfunc.h"
#include "constants.h"




struct cDiffCrossSectionAngularProfile {
  double ScatteringAngle,DifferentialCrossSection;
};

struct cDiffCrossSectionAngularProfile_ComparisonKey {
  bool operator () (const cDiffCrossSectionAngularProfile& t1, const cDiffCrossSectionAngularProfile& t2) const {
    return t1.ScatteringAngle < t2.ScatteringAngle;
  }
};

struct cDiffCrossSection_AngularProfile {
  double Energy,TotalIntegralCrossSection;
  int nPoints;
  vector<cDiffCrossSectionAngularProfile> data;
};

struct cDiffCrossSection_AngularProfile_ComarisonKey {
  bool operator () (const cDiffCrossSection_AngularProfile& t1, const cDiffCrossSection_AngularProfile& t2) const {
    return t1.Energy < t2.Energy;
  }
};


#define _DIFF_CROSS_SECTION__nANGULAR_POINTS_DEFAULT_        400
#define _DIFF_CROSS_SECTION__nENERGY_LAYERS_DEFAULT_         500
#define _DIFF_CROSS_SECTION__nCUMULATIVE_DISTRIBUTION_BINS_  400

class cDiffCrossSection {
private:
  bool Initialized;
  int nAngularPoints,nEnergyLayers,nCumulativeDistributionBins,ThisThread;

  double dLodEnergy,dScatteringAngle,**DifferentialCrossSectionTable,**cosCumulativeDistributionFunctionScatteringAngleBins;   //DifferentialCrossSectionTable,CumulativeDistributionTable[Energy Level][Scatterting Angle]
  double *TotalCrossSectionTable,ScatteringAngleMinimumIntegratedValue,maxTotalCrossSection;

  std::vector<cDiffCrossSection_AngularProfile> CrossSectionData;

public:

  double minEnergy,maxEnergy;
  double minLogEnergy,maxLogEnergy;
  double minScatteringAngle,maxScatteringAngle;

  cDiffCrossSection() {
    Initialized=false,ThisThread=0;
    nAngularPoints=_DIFF_CROSS_SECTION__nANGULAR_POINTS_DEFAULT_;
    nEnergyLayers=_DIFF_CROSS_SECTION__nENERGY_LAYERS_DEFAULT_;
    nCumulativeDistributionBins=_DIFF_CROSS_SECTION__nCUMULATIVE_DISTRIBUTION_BINS_;
    minEnergy=-1.0,maxEnergy=-1.0,minLogEnergy=-1.0,maxLogEnergy=-1.0;
    dLodEnergy=-1.0,dScatteringAngle=-1.0,DifferentialCrossSectionTable=NULL,TotalCrossSectionTable=NULL,cosCumulativeDistributionFunctionScatteringAngleBins=NULL;
    ScatteringAngleMinimumIntegratedValue=-1.0,minScatteringAngle=-1.0,maxScatteringAngle=-1.0,maxTotalCrossSection=-1.0;
  }

  void Set_nAngularPoints(int t) {nAngularPoints=t;}
  void Set_nEnergyLayers(int t) {nEnergyLayers=t;}
  void Set_nCumulativeDistributionBins(int t) {nCumulativeDistributionBins=t;}

  void Set_ScatteringAngleMinimumIntegratedValue_Degrees(double t) {
    //determine the number of the processor
    int mpiInitFlag;
    MPI_Initialized(&mpiInitFlag);
    if (mpiInitFlag==true)  MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);

    if (ThisThread==0) cout << "Diff. Cross Section: set lower limit on the scattering angles: min Scattering Angle=" << t << " degree" << endl;
    ScatteringAngleMinimumIntegratedValue=t;
  }


  void AddProfile(double Energy,int nPoints,const cDiffCrossSectionAngularProfile *CrossSectionProfile) {
    cDiffCrossSection_AngularProfile t;
    int i;
    double TotalIntegralCrossSection=0.0;

    double test=0.0;

    //determine the number of the processor
    int mpiInitFlag;
    MPI_Initialized(&mpiInitFlag);
    if (mpiInitFlag==true)  MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);

    //no additional profiles can be added after the cross section model is initalized
    if (Initialized==true) exit(__LINE__,__FILE__,"Error: the model is already initialied, no additional profileds can be added");

    //copy the buffer
    for (i=0,t.nPoints=0;i<nPoints;i++) if (CrossSectionProfile[i].ScatteringAngle>=ScatteringAngleMinimumIntegratedValue) {
      t.data.push_back(CrossSectionProfile[i]);
      t.data[t.nPoints].ScatteringAngle*=Pi/180.0;
      ++t.nPoints;
    }

    //sort the newly loaded cross section profile to be sure that the cross section angle is always increasing
    std::sort(t.data.begin(),t.data.end(),cDiffCrossSectionAngularProfile_ComparisonKey());


    //calcualte the total integral cross section
    for (i=0;i<t.nPoints-1;i++) {
      double alfa,a0,a1;

      a0=t.data[i].ScatteringAngle;
      a1=t.data[i+1].ScatteringAngle;
      alfa=log(t.data[i+1].DifferentialCrossSection/t.data[i].DifferentialCrossSection)/(a1-a0);

      if (alfa>-50000.0) TotalIntegralCrossSection+=2.0*Pi*(exp(alfa*a1)*(alfa*sin(a1)-cos(a1))-exp(alfa*a0)*(alfa*sin(a0)-cos(a0)))/(1.0+alfa*alfa)*t.data[i].DifferentialCrossSection/exp(alfa*a0);

      test+=2.0*Pi*0.5*(t.data[i+1].DifferentialCrossSection*sin(a1)+t.data[i].DifferentialCrossSection*sin(a0))*(a1-a0);
    }


    t.Energy=Energy;
    t.TotalIntegralCrossSection=TotalIntegralCrossSection;
    CrossSectionData.push_back(t);

    if ((minEnergy<0.0)||(minEnergy>Energy)) minEnergy=Energy;
    if (maxEnergy<Energy) maxEnergy=Energy;

    //find the minimum and maximum values of the scattering angle that present in ALL data sets
    if ((minScatteringAngle<0.0)||(minScatteringAngle<t.data[0].ScatteringAngle)) minScatteringAngle=t.data[0].ScatteringAngle;
    if ((maxScatteringAngle<0.0)||(maxScatteringAngle>t.data[t.nPoints-1].ScatteringAngle)) maxScatteringAngle=t.data[t.nPoints-1].ScatteringAngle;

    if (ThisThread==0) std::cout << "Diff. Cross Section: Collision Energy=" << t.Energy << "[J] (" << t.Energy/eV2J << "[eV]), Total Collision Cross Section=" << t.TotalIntegralCrossSection <<
        ", Applied range of scattering angles=" << t.data[0].ScatteringAngle/Pi*180.0 << " to " << t.data[t.nPoints-1].ScatteringAngle/Pi*180.0 << " [degrees]" << std::endl;
  }

  void InitInternalTable() {
    int iEnergyLayer,iScatteringAnglePoint,l0,l1,il0,il1,i;
    double E,ScattertingAngle,CrossSectionL0,CrossSectionL1,TotalIntegralCrossSection,alfa;

    //the model cannnot be re-initialized
    if (Initialized==true) exit(__LINE__,__FILE__,"Error: the model is already initialied");
    Initialized=true;

    if (ThisThread==0) {
      cout << "Scattering Angle Limit: " << minScatteringAngle/Pi*180.0 << " to " <<  maxScatteringAngle/Pi*180.0 << "[degrees]" << endl;
    }

    if (CrossSectionData.size()==0) exit(__LINE__,__FILE__,"Error: not differential cross section profile present");

    //If only one prifile for the collision crossection is added to the model => limit the number of the nEnergyLayers to 2: this will be the minimum number of the layers that is consistent with the rest of the differential cross section model
    if (CrossSectionData.size()==1) {
      nEnergyLayers=2;
    }


    std::sort(CrossSectionData.begin(),CrossSectionData.end(),cDiffCrossSection_AngularProfile_ComarisonKey());

    minLogEnergy=log(minEnergy);
    maxLogEnergy=log(maxEnergy);
    dLodEnergy=(maxLogEnergy-minLogEnergy)/(nEnergyLayers-1);
    dScatteringAngle=(maxScatteringAngle-minScatteringAngle)/(nAngularPoints-1);

    DifferentialCrossSectionTable=new double* [nEnergyLayers];
    DifferentialCrossSectionTable[0]=new double [nEnergyLayers*nAngularPoints];

    cosCumulativeDistributionFunctionScatteringAngleBins=new double* [nEnergyLayers];
    cosCumulativeDistributionFunctionScatteringAngleBins[0]=new double [nEnergyLayers*(nCumulativeDistributionBins+1)];


    TotalCrossSectionTable=new double [nEnergyLayers];

    for (iEnergyLayer=1;iEnergyLayer<nEnergyLayers;iEnergyLayer++) {
      DifferentialCrossSectionTable[iEnergyLayer]=DifferentialCrossSectionTable[iEnergyLayer-1]+nAngularPoints;
      cosCumulativeDistributionFunctionScatteringAngleBins[iEnergyLayer]=cosCumulativeDistributionFunctionScatteringAngleBins[iEnergyLayer-1]+nCumulativeDistributionBins+1;
    }

    //initialize the cross section table
    for (iEnergyLayer=0;iEnergyLayer<nEnergyLayers;iEnergyLayer++) {
      E=minEnergy*exp(dLodEnergy*iEnergyLayer);
      il0=1,il1=1;

      for (l0=0;l0<(int)CrossSectionData.size();l0++) if (CrossSectionData[l0].Energy<=E) break;

      if (l0>=(int)CrossSectionData.size()) l0=CrossSectionData.size()-1;
      l1=l0+1;
      if (l1>=(int)CrossSectionData.size()) l1=l0;

      for (iScatteringAnglePoint=0;iScatteringAnglePoint<nAngularPoints;iScatteringAnglePoint++) {
        ScattertingAngle=minScatteringAngle+iScatteringAnglePoint*dScatteringAngle;

        //find the value of the dufferential cross section at each energy levels used for the interpolation
        //the energy level L0
        for (;il0<CrossSectionData[l0].nPoints;il0++) if (ScattertingAngle<=CrossSectionData[l0].data[il0].ScatteringAngle) break;
        alfa=log(CrossSectionData[l0].data[il0].DifferentialCrossSection/CrossSectionData[l0].data[il0-1].DifferentialCrossSection)/(CrossSectionData[l0].data[il0].ScatteringAngle-CrossSectionData[l0].data[il0-1].ScatteringAngle);
        CrossSectionL0=(alfa>-50000.0) ? CrossSectionData[l0].data[il0-1].DifferentialCrossSection*exp(alfa*(ScattertingAngle-CrossSectionData[l0].data[il0-1].ScatteringAngle)) : CrossSectionData[l0].data[il0-1].DifferentialCrossSection;

        //the energy level L1
        for (;il1<CrossSectionData[l1].nPoints;il1++) if (ScattertingAngle<=CrossSectionData[l1].data[il1].ScatteringAngle) break;
        alfa=log(CrossSectionData[l1].data[il1].DifferentialCrossSection/CrossSectionData[l1].data[il1-1].DifferentialCrossSection)/(CrossSectionData[l1].data[il1].ScatteringAngle-CrossSectionData[l1].data[il1-1].ScatteringAngle);
        CrossSectionL1=(alfa>-50000.0) ? CrossSectionData[l1].data[il1-1].DifferentialCrossSection*exp(alfa*(ScattertingAngle-CrossSectionData[l1].data[il1-1].ScatteringAngle)) : CrossSectionData[l1].data[il1-1].DifferentialCrossSection;


        //interpolate the cross section value over the energy
        if (l0!=l1) {
          DifferentialCrossSectionTable[iEnergyLayer][iScatteringAnglePoint]=CrossSectionL0*exp(log(CrossSectionL1/CrossSectionL0)/log(CrossSectionData[l1].Energy/CrossSectionData[l0].Energy)*
            log(E/CrossSectionData[l0].Energy));
        } else DifferentialCrossSectionTable[iEnergyLayer][iScatteringAnglePoint]=CrossSectionL0;

      }

      //get the total collision integral for the current energy level
      for (i=0,TotalIntegralCrossSection=0.0;i<nAngularPoints-1;i++) {
        double a0,a1;

        a0=minScatteringAngle+i*dScatteringAngle;
        a1=minScatteringAngle+(i+1)*dScatteringAngle;
        alfa=log(DifferentialCrossSectionTable[iEnergyLayer][i+1]/DifferentialCrossSectionTable[iEnergyLayer][i])/(a1-a0);

        if (alfa>-50000.0) TotalIntegralCrossSection+=2.0*Pi*(exp(alfa*a1)*(alfa*sin(a1)-cos(a1))-exp(alfa*a0)*(alfa*sin(a0)-cos(a0)))/(1.0+alfa*alfa)*DifferentialCrossSectionTable[iEnergyLayer][i]/exp(alfa*a0);
      }

      TotalCrossSectionTable[iEnergyLayer]=TotalIntegralCrossSection;
      if (maxTotalCrossSection<TotalIntegralCrossSection) maxTotalCrossSection=TotalIntegralCrossSection;

      //calculate the cumulative distrivution function
      double dCumulativeDistribution,CumulativeDistribution=0.0;
      int iCummulativeDistributionIntervalFilled=0;
      double a0,a1,aHalf;
      int SearchSegment=1;
      double CumulativeDistributionBinIntegralFraction=TotalIntegralCrossSection/nCumulativeDistributionBins;

      for (i=0,a0=minScatteringAngle,cosCumulativeDistributionFunctionScatteringAngleBins[iEnergyLayer][0]=cos(minScatteringAngle);i<nAngularPoints-1;i++) {
        a1=minScatteringAngle+(i+1)*dScatteringAngle;
        SearchSegment=1;

        alfa=log(DifferentialCrossSectionTable[iEnergyLayer][i+1]/DifferentialCrossSectionTable[iEnergyLayer][i])/dScatteringAngle;

        dCumulativeDistribution=(alfa>-50000.0) ? 2.0*Pi*(exp(alfa*a1)*(alfa*sin(a1)-cos(a1))-exp(alfa*a0)*(alfa*sin(a0)-cos(a0)))/(1.0+alfa*alfa)*DifferentialCrossSectionTable[iEnergyLayer][i]/exp(alfa*(minScatteringAngle+i*dScatteringAngle)) : 0.0;

        if ((CumulativeDistribution+dCumulativeDistribution)>CumulativeDistributionBinIntegralFraction*(iCummulativeDistributionIntervalFilled+1.0)) {

          for (int nTry=0;(nTry<10)||((CumulativeDistribution+dCumulativeDistribution)>=CumulativeDistributionBinIntegralFraction*(iCummulativeDistributionIntervalFilled+2.0));nTry++) {
            aHalf=0.5*(a0+a1);
            dCumulativeDistribution=(alfa>-50000.0) ? 2.0*Pi*(exp(alfa*aHalf)*(alfa*sin(aHalf)-cos(aHalf))-exp(alfa*a0)*(alfa*sin(a0)-cos(a0)))/(1.0+alfa*alfa)*DifferentialCrossSectionTable[iEnergyLayer][i]/exp(alfa*(minScatteringAngle+i*dScatteringAngle)) : 0.0;

            if ((CumulativeDistribution+dCumulativeDistribution)>CumulativeDistributionBinIntegralFraction*(iCummulativeDistributionIntervalFilled+1.0)) {
              a1=aHalf;
              SearchSegment=0;
            }
            else {
              CumulativeDistribution+=dCumulativeDistribution;
              a0=aHalf;
            }
          }

          if (alfa>-50000.0) CumulativeDistribution+=2.0*Pi*(exp(alfa*a1)*(alfa*sin(a1)-cos(a1))-exp(alfa*a0)*(alfa*sin(a0)-cos(a0)))/(1.0+alfa*alfa)*DifferentialCrossSectionTable[iEnergyLayer][i]/exp(alfa*(minScatteringAngle+i*dScatteringAngle));
          cosCumulativeDistributionFunctionScatteringAngleBins[iEnergyLayer][++iCummulativeDistributionIntervalFilled]=cos(a1);

          if (iCummulativeDistributionIntervalFilled>nCumulativeDistributionBins) exit(__LINE__,__FILE__,"Error: iCummulativeDistributionIntervalFilled is out of range");

          a0=a1;
          if (SearchSegment==0) --i;
          continue;
        }
        else {
          CumulativeDistribution+=dCumulativeDistribution;
          a0=a1;
        }

      }

      if (iCummulativeDistributionIntervalFilled<nCumulativeDistributionBins-1) exit(__LINE__,__FILE__,"Error: something wrong in the integration of the cumulative distrtibution function");
      else if (iCummulativeDistributionIntervalFilled==nCumulativeDistributionBins-1) cosCumulativeDistributionFunctionScatteringAngleBins[iEnergyLayer][nCumulativeDistributionBins]=cos(minScatteringAngle+nAngularPoints*dScatteringAngle);

//      cout << CumulativeDistribution << endl;


    }
  }

  void PrintCrossSectionPlot() {
    int i,j;
    FILE *fout;
    char str[10000];

    if (ThisThread!=0) return;

    //print original data sets
    for (i=0;i<(int)CrossSectionData.size();i++) {
      sprintf(str,"DiffCrossSection.profile=%i.E=%e[J].E=%e[eV].dat",i,CrossSectionData[i].Energy,CrossSectionData[i].Energy/eV2J);
      fout=fopen(str,"w");
      fprintf(fout,"VARIABLES=\"Scattering Angle [degrees]\", \"Differential Cross Section [m^2]\"");

      for (j=0;j<(int)CrossSectionData[i].data.size();j++) fprintf(fout,"%e %e\n",CrossSectionData[i].data[j].ScatteringAngle/Pi*180.0,CrossSectionData[i].data[j].DifferentialCrossSection);

      fclose(fout);
    }

    //output the interpolated differential cross section
    fout=fopen("DifferentialCrossSection.dat","w");

    fprintf(fout,"VARIABLES=\"Scattering Angle\"");
    for (i=0;i<nEnergyLayers;i++) fprintf(fout," ,\"Sigma[E=%eV]\"",minEnergy*exp(dLodEnergy*i)/eV2J);
    fprintf(fout,"\n");

    for (i=0;i<nAngularPoints;i++) {
      fprintf(fout,"%e ",(minScatteringAngle+i*dScatteringAngle)/Pi*180.0);

      for (j=0;j<nEnergyLayers;j++) fprintf(fout,"%e ",DifferentialCrossSectionTable[j][i]);
      fprintf(fout,"\n");
    }

    fclose(fout);
  }

  //return the total integrate cross section for a particular energy level
  int GetEnergyLayerNumber(double CollisionEnergy) {
    int EnergyLayer;

    if (CollisionEnergy<=minEnergy) EnergyLayer=0;
    else if (CollisionEnergy>=maxEnergy) EnergyLayer=nEnergyLayers-1;
    else EnergyLayer=(int)(log(CollisionEnergy/minEnergy)/dLodEnergy);

    return EnergyLayer;
  }

  double GetTotalCrossSection(double CollisionEnergy) {
    return TotalCrossSectionTable[GetEnergyLayerNumber(CollisionEnergy)];
  }

  double GetRandomScatteringAngle(double CollisionEnergy) {
    int EnergyLayer,CumulativeDistributionBin;
    double cosA0,cosA1;

    EnergyLayer=GetEnergyLayerNumber(CollisionEnergy);
    CumulativeDistributionBin=(int)(rnd()*nCumulativeDistributionBins);

    cosA0=cosCumulativeDistributionFunctionScatteringAngleBins[EnergyLayer][CumulativeDistributionBin];
    cosA1=cosCumulativeDistributionFunctionScatteringAngleBins[EnergyLayer][CumulativeDistributionBin+1];

    return acos(cosA1+rnd()*(cosA0-cosA1));
  }

  double GetMaxTotalCrossSection() {return maxTotalCrossSection;}
};

#endif /* DIFFERENTIALCROSSSECTION_H_ */
