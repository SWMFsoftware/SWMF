//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * ElectricallyChargedDust.h
 *
 *  Created on: Jan 20, 2012
 *      Author: vtenishe
 */

//$Id$

//#include "pic.h"
#include "constants.h"



#ifndef ELECTRICALLYCHARGEDDUST_H_



using namespace std;



/*--------------------------------  Definitions of the model ---------------------------*/
//#define _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_      0
//#define _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__OFF_     1

#define _DUST__TIME_STEP_SUBGROUPING_MODE__ON_       0
#define _DUST__TIME_STEP_SUBGROUPING_MODE__OFF_      1

//depending the velocity of the dust grains change the time step sub-group of the dust grains
#define _DUST__TIME_STEP_SUBGROUP_MIGRATION_MODE__ON_   0
#define _DUST__TIME_STEP_SUBGROUP_MIGRATION_MODE__OFF_  1

//use a hard-wired constant or a user-provided function for calcualtion of dust grain properties
#define _DUST__CALCULATION_GRAIN_INTERNAL_PROPERTIES_MODE__CONSTANT_VALUE_         0
#define _DUST__CALCULATION_GRAIN_INTERNAL_PROPERTIES_MODE__USER_DEFINED_FUNCTION_  1


/*--------------------------------  Initial parameters of the model --------------------*/


#define _DUST__TIME_STEP_SUBGROUPING_MODE_ _DUST__TIME_STEP_SUBGROUPING_MODE__ON_
#define _DUST__TIME_STEP_SUBGROUP_MIGRATION_MODE_ _DUST__TIME_STEP_SUBGROUP_MIGRATION_MODE__ON_
#define _DUST__CALCULATION_GRAIN_INTERNAL_PROPERTIES_MODE_ _DUST__CALCULATION_GRAIN_INTERNAL_PROPERTIES_MODE__CONSTANT_VALUE_

/*-------------------------------- The model -------------------------------------------*/

namespace ElectricallyChargedDust {




  //the dust spece groups
  //  extern int _DUST_SPEC_;
  extern int nDustTimeStepGroups;
  extern double minDustRadius,maxDustRadius;

  //characteristic limits of the dust grain's velocityes
  extern double minDustVelocity,maxDustVelocty;

  //the locariphmic incremernt in the dust grains's size for each time step groups
//  extern double dLogGrainRadius_TimeStepGroup;


  //sample the mass injection of the dust grains
  extern double SampledDustMassInjectionRate;

  //Drag Coefficient
  typedef double (*fGrainDragCoeffiient)(double *GrainVelocity,double GrainRadius);
  const double GrainDragCoefficient=2.0;
  extern fGrainDragCoeffiient GrainDragCoeffiient_UserDefinedFunction;

  //the mean desntiy of the dust grains
  typedef double (*fGrainMassDensity)(double GrainRadius);
  static const double MeanDustDensity=1.0E3;
  extern fGrainMassDensity GrainMassDensity_UserDefinedFunction;

  //the total mass production rate of the dust grains
  extern double TotalMassDustProductionRate;

  //set ingternal proporties of the dust grains
  inline double GetGrainCharge(PIC::ParticleBuffer::byte *ParticleData) {
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    return *((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_));
#else
    exit(__LINE__,__FILE__,"Error: in the model configuration no charing of the dust grains is allowed");
#endif
  }

  inline void SetGrainCharge(double Charge,PIC::ParticleBuffer::byte *ParticleData) {
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    *((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_))=Charge;
#else
    exit(__LINE__,__FILE__,"Error: in the model configuration no charing of the dust grains is allowed");
#endif
  }

  inline double GetGrainMass(PIC::ParticleBuffer::byte *ParticleData) {return *((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_));}
  inline void SetGrainMass(double Mass,PIC::ParticleBuffer::byte *ParticleData) {*((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_))=Mass;}

  inline double GetGrainRadius(PIC::ParticleBuffer::byte *ParticleData) {return *((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_));}
  inline void SetGrainRadius(double Radius,PIC::ParticleBuffer::byte *ParticleData) {
    *((double*)(ParticleData+_PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_))=Radius;

    if ((Radius>maxDustRadius)||(Radius<minDustRadius)) exit(__LINE__,__FILE__,"Error: grains' radius is out of range");
  }


  //time step grouping of the grains
  namespace GrainVelocityGroup {
    extern double minGrainVelocity,maxGrainVelocity;
    extern double logMinGrainVelocity,logMaxGrainVelocity,dLogGrainVelocity;
    extern int nGroups;

    inline void GetGroupVelocityRange(double &vmin,double &vmax,int ngroup) {
      vmin=minGrainVelocity*exp(ngroup*dLogGrainVelocity);
      vmax=minGrainVelocity*exp((ngroup+1)*dLogGrainVelocity);
    }

    inline int GetGroupNumber(double *v) {
      double Speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

      if (Speed<=minGrainVelocity) return 0;
      else if (Speed>=maxGrainVelocity) return nGroups-1;
      else return (int)((log(Speed)-logMinGrainVelocity)/dLogGrainVelocity);
    }

    inline void Init() {
      logMinGrainVelocity=log(minGrainVelocity);
      logMaxGrainVelocity=log(maxGrainVelocity);
      dLogGrainVelocity=(logMaxGrainVelocity-logMinGrainVelocity)/nGroups;
    }
  }

  //sampling disfferent size ranges of dust macrospchpic parameters
  namespace Sampling {
    extern int nDustSizeSamplingIntervals;
    extern double dLogDustSamplingIntervals;
    extern int CellSamplingDataOffset,DustSizeSamplingInterval_NumberDensity_Offset,DustSizeSamplingInterval_Velocity_Offset,DustSizeSamplingInterval_Speed_Offset;
    extern int DustSizeSamplingInterval_ElectricCherge_Offset,DustSizeSamplingInterval_ElectricCurrent_Offset;
    extern int TotalDustElectricChargeDensitySamplingOffset,TotalDustElectricCurrentDensitySamplingOffset;

    void SetDustSamplingIntervals(int);
    int RequestSamplingData(int);

    void inline SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec) {
      double GrainRadius,v[3]={0.0,0.0,0.0},Speed=0.0;
      int dustSamplingInterval,idim;

      //determine if the particle is a dust grain
      if (!((_DUST_SPEC_<=spec)&&(spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups))) return;

      GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);
      dustSamplingInterval=(int)(log(GrainRadius/minDustRadius)/dLogDustSamplingIntervals);
      PIC::ParticleBuffer::GetV(v,(PIC::ParticleBuffer::byte*)ParticleData);

      //number density
      *(dustSamplingInterval+(double*)(SamplingBuffer+DustSizeSamplingInterval_NumberDensity_Offset))+=LocalParticleWeight;

      //velocity
      for (idim=0;idim<3;idim++) {
        v[idim]*=LocalParticleWeight;
        Speed+=v[idim]*v[idim];

        *(idim+3*dustSamplingInterval+(double*)(SamplingBuffer+DustSizeSamplingInterval_Velocity_Offset))+=v[idim];
      }

      //Speed
      Speed=sqrt(Speed);
      *(dustSamplingInterval+(double*)(SamplingBuffer+DustSizeSamplingInterval_Speed_Offset))+=Speed;

      //electric charge and currents
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
      double grainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);

      *(dustSamplingInterval+(double*)(SamplingBuffer+DustSizeSamplingInterval_ElectricCherge_Offset))+=grainElectricCharge*LocalParticleWeight;
      *((double*)(SamplingBuffer+Sampling::TotalDustElectricChargeDensitySamplingOffset))+=grainElectricCharge*LocalParticleWeight;

      for (idim=0;idim<3;idim++) {
        v[idim]*=grainElectricCharge;

        *(idim+3*dustSamplingInterval+(double*)(SamplingBuffer+DustSizeSamplingInterval_ElectricCurrent_Offset))+=v[idim];
        *(idim+(double*)(SamplingBuffer+Sampling::TotalDustElectricCurrentDensitySamplingOffset))+=v[idim];
      }
#endif
    }


    namespace SampleSizeDistributionFucntion {
      extern double dLog10DustRadius;

      //the number of the sampling size intervals
      extern int nSamplingIntervals;

      //sampling buffers
      extern double **SamplingBuffer;

      //total sampling offsets
      extern int SamplingIntervalDataLength;
      extern int NumberDensitySamplingOffset,VelocitySamplingOffset,SpeedSamplingOffset;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
      extern int ElectricChargeSamplingOffset,ElectricCurentSamplingOffset,AbsoluteElectricCurrentValueSamplingOffset;
#endif

      //the number of physical locations where the sampling will be performed; the list of the samplng locations; sanpling nodes; sampling local cells's number
      extern int nSamplingLocations;
      extern double **SamplingLocations;
      extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
      extern long int *SampleLocalCellNumber;

      void Init(double ProbeLocations[][DIM],int nProbeLocations,int nIntervals);

      void SampleDistributionFnction();
      void flushSamplingBuffers();

      void printDistributionFunction(int DataOutputFileNumber);
    }

  }

  //side ditribution of the dust grains: The distribution is calcualted as the LOGARIPHM!!!! of the grain radius
  namespace SizeDistribution {
    extern double PowerIndex,NormalizationFactor,LogMinDustGrainRadius,LogMaxDustGrainRadius,maxValue__non_normalized;

    inline double GetNormalizedValue(double GrainRadius) {
      return pow(GrainRadius,1.0-PowerIndex)/NormalizationFactor;
    }

    inline void PrintDistributionFucntion() {
      FILE *fout;

      const int nOutputPoints=400;

      if (PIC::ThisThread==0) {
        char fname[_MAX_STRING_LENGTH_PIC_];

        sprintf(fname,"%s/DustGrainSizeDistribution.dat",PIC::OutputDataFileDirectory);
        fout=fopen(fname,"w");
        fprintf(fout,"VARIABLE=\"Grain Radius\", \"Normalized Distribution Function\" \n");

        double r,dLogR=(LogMaxDustGrainRadius-LogMinDustGrainRadius)/(nOutputPoints-1);

        for (int i=0;i<nOutputPoints;i++) {
          r=exp(LogMinDustGrainRadius+dLogR*i);

          fprintf(fout,"%e  %e\n",r,GetNormalizedValue(r)/r); //transfer the distribution from f(ln R) -> f(R)
        }

        fclose(fout);
      }
    }

    inline void Init() {
      LogMinDustGrainRadius=log(minDustRadius);
      LogMaxDustGrainRadius=log(maxDustRadius);
      NormalizationFactor=(exp((1.0-PowerIndex)*LogMaxDustGrainRadius)-exp((1.0-PowerIndex)*LogMinDustGrainRadius))/(1.0-PowerIndex);
      maxValue__non_normalized=pow(minDustRadius,1.0-PowerIndex);

      if (PIC::ThisThread==0) {
        cout << "$PREFIX:Dust Grains Distribution:" << endl;
        cout << "$PREFIX:minDustRadius=" << minDustRadius << endl << "$PREFIX:maxDustRadius=" << maxDustRadius << endl << "$PREFIX:PowerIndex=" << PowerIndex << endl << "$PREFIX:NormalizationFactor=" << NormalizationFactor << endl << endl;
      }

      PrintDistributionFucntion();
    }

    inline void GenerateGrainRandomRadius(double& newRadius,double &WeightCorrectionFactor) {
      newRadius=minDustRadius*exp(rnd()*(LogMaxDustGrainRadius-LogMinDustGrainRadius));
      WeightCorrectionFactor=pow(newRadius,1.0-PowerIndex)/maxValue__non_normalized;
    }
  }





  //sampling procedues
  void OutputSampledModelData(int DataOutputFileNumber);
  void SampleModelData();
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);




  //evaluate the local time step: the function returns 'false' when the evaluation of the time step is requested for a NON-DUST specie
  bool EvaluateLocalTimeStep(int spec,double &dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

  //init the dsut model
  void Init_AfterParser();
  void Init_BeforeParser();

  //acceleration of the dust grains
  void inline TotalGrainAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    char ParticleData[PIC::ParticleBuffer::ParticleDataLength];
    double GrainCharge,GrainMass;

    memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

    GrainCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);
    GrainMass=GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);



    char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
    PIC::Mesh::cDataCenterNode *CenterNode;
    double *E,*B;
    int i,j,k;

    //copy to local variables
    double accl_LOCAL[3]={0.0,0.0,0.0},x_LOCAL[3],v_LOCAL[3];

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));


  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    long int nd;

    if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
      exit(__LINE__,__FILE__,"Error: the cell is not found");
    }

    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  #endif

    CenterNode=startNode->block->GetCenterNode(nd);
    memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+CenterNode->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);

    //Lorentz force
    E=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::ElectricFieldOffset);
    B=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::MagneticFieldOffset);

    accl_LOCAL[0]+=GrainCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/GrainMass;
    accl_LOCAL[1]+=GrainCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/GrainMass;
    accl_LOCAL[2]+=GrainCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/GrainMass;

    //the gravity force
    double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
    double c=GravityConstant*_MASS_(_TARGET_)/pow(r2,3.0/2.0);

    accl_LOCAL[0]-=c*x_LOCAL[0];
    accl_LOCAL[1]-=c*x_LOCAL[1];
    accl_LOCAL[2]-=c*x_LOCAL[2];

    //Drag force
    /*
    double A,cr2,GrainRadius;
    double *BackgroundAtmosphereBulkVelocity=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::NeutralBullVelocityOffset);
    double BackgroundAtmosphereNumberDensity=*((double*)(ICES_AssociatedData+PIC::CPLR::ICES::NeutralNumberDensityOffset));
    double GrainRadius=GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);

    cr2=(v_LOCAL[0]-BackgroundAtmosphereBulkVelocity[0])*(v_LOCAL[0]-BackgroundAtmosphereBulkVelocity[0])+
      (v_LOCAL[1]-BackgroundAtmosphereBulkVelocity[1])*(v_LOCAL[1]-BackgroundAtmosphereBulkVelocity[1])+
      (v_LOCAL[2]-BackgroundAtmosphereBulkVelocity[2])*(v_LOCAL[2]-BackgroundAtmosphereBulkVelocity[2]);


    A=Pi*pow(GrainRadius,2)/2.0*GrainDragCoefficient*sqrt(cr2)/GrainMass*BackgroundAtmosphereNumberDensity*_MASS_(_H2O_);

    accl_LOCAL[0]+=A*(BackgroundAtmosphereBulkVelocity[0]-v_LOCAL[0]);
    accl_LOCAL[1]+=A*(BackgroundAtmosphereBulkVelocity[1]-v_LOCAL[1]);
    accl_LOCAL[2]+=A*(BackgroundAtmosphereBulkVelocity[2]-v_LOCAL[2]);
    */

    //copy the accaleration vector from the internal buffer
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }


  //charging of the dust grains (as a generic transformation in the PIC code)
  inline int DustChargingProcessorIndicator(double *x,double *v,int spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double &dt,bool &TransformationTimeStepLimitFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
  }



//integrate the equation of charging of a duts grain

/*--------------------------------------     First Order Integration Scheme ----------------------------------------------------
inline int DustChargingProcessor(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {


  double electronCurrent,ionCurrent,plasmaTemperature,plasmaNumberDensity,basisElectronCurrent,basisIonCurrent;
  double GrainElectricCharge,GrainElectricCharge_FIRST_ORDER; //,GrainElectricCharge_SECOND_ORDER,GrainElectricCharge_SECOND_ORDER_MIDDLE;

  PIC::Mesh::cDataCenterNode* cell;
  int i,j,k;
  long int LocalCellNumber;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,initNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
  cell=initNode->block->GetCenterNode(LocalCellNumber);

  //ICES plasma data
  char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
  double *swVel;

  memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+cell->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);
  plasmaTemperature=*((double*)(PIC::CPLR::ICES::PlasmaTemperatureOffset+ICES_AssociatedData));
  swVel=(double*)(PIC::CPLR::ICES::PlasmaBulkVelocityOffset+ICES_AssociatedData);
  plasmaNumberDensity=*((double*)(PIC::CPLR::ICES::PlasmaNumberDensityOffset+ICES_AssociatedData));

  if (plasmaTemperature<=0.0) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;

  //get the grain electric potential
  char localParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double M,GrainRadius,dustPotential,xi;

  memcpy((void*)localParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
  GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)localParticleData);




//==========  DEBUG  BEGIN =================
  static long int nFunctionCalls=0;

  nFunctionCalls++;

  if ((PIC::nTotalThreads==7)&&(nFunctionCalls==1709)) {
    cout << __FILE__ << "@" << __LINE__ << endl;
  }

//==========   DEBUG END ===================


  for (int nTimeStepSplit=0;nTimeStepSplit<5;nTimeStepSplit++) {
    int nTotalSubSteps=(1<<nTimeStepSplit);
    double dtSubStepMax=dt/nTotalSubSteps;
    double dtSubStep=0.0,dtTotalCounter=0.0;
    double dtSubStepElectronCurrent,dtSubStepIonCurrent,IonCurrentMultiplier,ElectronCurrentMultiplier,newGrainElectricCharge,newDustPotential,newElectronCurrentMultiplier,newIonCurrentMultiplier;
    double dtSubStepMaxAbsolute=dtSubStepMax;



    GrainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)localParticleData);
    GrainElectricCharge_FIRST_ORDER=GrainElectricCharge;
//    GrainElectricCharge_SECOND_ORDER=GrainElectricCharge;

    basisElectronCurrent=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ElectronMass));
    basisIonCurrent=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ProtonMass));



    dtSubStep=-1.0;
    bool TimeStepIntegrationFinished=false;


//    while (TimeStepIntegrationFinished==false) {   ////(dtTotalCounter<dt) {
      dtSubStepElectronCurrent=dtSubStepMax;
      dtSubStepIonCurrent=dtSubStepMax;

      //FIRST ORDER INTEGRATION
      dustPotential=GrainElectricCharge_FIRST_ORDER/(4.0*Pi*VacuumPermittivity*GrainRadius);

      //get the electron current (Horanyi-1996-ARAA, Eq. 4)
      xi=-ElectronCharge*dustPotential/(Kbol*plasmaTemperature);

//      electronCurrent=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ElectronMass));
      ElectronCurrentMultiplier=(xi>=0.0) ? exp(-xi) : 1.0-xi;
      electronCurrent=ElectronCurrentMultiplier*basisElectronCurrent;

      if (ElectronCurrentMultiplier<1.0) {
        if (dtSubStepElectronCurrent*ElectronCurrentMultiplier>0.1) dtSubStepElectronCurrent=0.1/ElectronCurrentMultiplier;
      }
      else dtSubStepElectronCurrent=0.1/ElectronCurrentMultiplier;


      //get the ion current (Horanyi-1996-ARAA, Eq. 4,14)

      xi=ElectronCharge*dustPotential/(Kbol*plasmaTemperature);
      M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*plasmaTemperature/ProtonMass));

//      ionCurrent=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ProtonMass));

      if (dustPotential<0.0) {
        IonCurrentMultiplier=((pow(M,2)+0.5-xi)*sqrt(Pi)/M*erf(M)+exp(-pow(M,2)))/2.0;
      }
      else {
        double sqrt_xi=sqrt(xi);

        IonCurrentMultiplier=(
            (pow(M,2)+0.5-xi)*sqrt(Pi)/M*(erf(M+sqrt_xi)+erf(M-sqrt_xi))  +
            (sqrt(xi/M)+1.0)*exp(-pow(M-sqrt_xi,2.0)) -
            (sqrt(xi/M)-1.0)*exp(-pow(M+sqrt_xi,2.0))
            )/2.0;
      }

      ionCurrent=IonCurrentMultiplier*basisIonCurrent;

      if (IonCurrentMultiplier<1.0) {
        if (dtSubStepIonCurrent*IonCurrentMultiplier>0.1) dtSubStepIonCurrent=0.1/IonCurrentMultiplier;
      }
      else dtSubStepIonCurrent=0.1/IonCurrentMultiplier;



      //increase the value of 'dtSubStepMax' -> in case if the integration procedure will allow increasing if the integration step
//      dtSubStepMax*=2.0;

//      dtSubStep=min(min(min(dt-dtTotalCounter,dtSubStepMax),min(dtSubStepElectronCurrent,dtSubStepIonCurrent)),dtSubStepMaxAbsolute);

      if (dtSubStep<0.0) dtSubStep=min(min(min(dt-dtTotalCounter,dtSubStepMax),min(dtSubStepElectronCurrent,dtSubStepIonCurrent)),dtSubStepMaxAbsolute);


      dtSubStep/=4.0;

      while (TimeStepIntegrationFinished==false) {


      dtSubStep*=4.0;
      if (dtSubStep+dtTotalCounter>=dt) dtSubStep=dt-dtTotalCounter;

      do {
      //update the grain charge (Horanyi-1996-ARAA, Eq.1)
      newGrainElectricCharge=GrainElectricCharge_FIRST_ORDER+(ionCurrent+electronCurrent)*dtSubStep;
      newDustPotential=newGrainElectricCharge/(4.0*Pi*VacuumPermittivity*GrainRadius);

      xi=-ElectronCharge*newDustPotential/(Kbol*plasmaTemperature);
      newElectronCurrentMultiplier=(xi>=0.0) ? exp(-xi) : 1.0-xi;

      xi=ElectronCharge*newDustPotential/(Kbol*plasmaTemperature);
      if (newDustPotential<0.0) {
        newIonCurrentMultiplier=((pow(M,2)+0.5-xi)*sqrt(Pi)/M*erf(M)+exp(-pow(M,2)))/2.0;
      }
      else {
        double sqrt_xi=sqrt(xi);

        newIonCurrentMultiplier=(
            (pow(M,2)+0.5-xi)*sqrt(Pi)/M*(erf(M+sqrt_xi)+erf(M-sqrt_xi))  +
            (sqrt(xi/M)+1.0)*exp(-pow(M-sqrt_xi,2.0)) -
            (sqrt(xi/M)-1.0)*exp(-pow(M+sqrt_xi,2.0))
            )/2.0;

      }

      if (  ((fabs(newElectronCurrentMultiplier-ElectronCurrentMultiplier)<=0.1*fabs(newElectronCurrentMultiplier+ElectronCurrentMultiplier)) &&
          (fabs(newIonCurrentMultiplier-IonCurrentMultiplier)<=0.1*fabs(newIonCurrentMultiplier+IonCurrentMultiplier)))  ||
          (dtSubStep<1.0E-5*dt) ) break;


      dtSubStep/=4.0;
 //     dtSubStepMax=dtSubStep;

      }
      while (true); // (dtSubStep>1.0E-4*dt);




//==================   DEBUG ==================

if (fabs(newGrainElectricCharge)>1.0E-3) {
  cout << __FILE__ << "%" << __LINE__ << endl;
}

//==================   END DEBUG ==============




      dustPotential=newDustPotential;

      ElectronCurrentMultiplier=newElectronCurrentMultiplier;
      IonCurrentMultiplier=newIonCurrentMultiplier;

      electronCurrent=ElectronCurrentMultiplier*basisElectronCurrent;
      ionCurrent=IonCurrentMultiplier*basisIonCurrent;

      dtTotalCounter+=dtSubStep;
      GrainElectricCharge_FIRST_ORDER=newGrainElectricCharge;

      if (dtSubStep+dtTotalCounter>=dt) TimeStepIntegrationFinished=true;

      }
  break;

  }


//  if (GrainElectricCharge_SECOND_ORDER>1.0) exit(__LINE__,__FILE__,"Error: out of range");


  SetGrainCharge(GrainElectricCharge_FIRST_ORDER,ParticleData);

  //move the particle into diferent velocity group if needed
  int oldVelocityGroup,newVelocityGroup;

  oldVelocityGroup=spec-_DUST_SPEC_;
  newVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);

  if (oldVelocityGroup!=newVelocityGroup) {
    //move the particle into different velocity group
    double GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

    GrainWeightCorrection*=initNode->block->GetLocalTimeStep(_DUST_SPEC_+newVelocityGroup)/initNode->block->GetLocalTimeStep(_DUST_SPEC_+oldVelocityGroup);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+oldVelocityGroup]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+newVelocityGroup];
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,ParticleData);
    spec=_DUST_SPEC_+newVelocityGroup;
    PIC::ParticleBuffer::SetI(spec,ParticleData);
  }



  return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
}

*/

  /*---------------------------------------------------    Second Order Integration Scheme -----------------------------------*/
  inline int DustChargingProcessor_Implicit_SecondOrder(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {
  double plasmaTemperature,plasmaNumberDensity;
  double GrainElectricCharge,GrainElectricCharge_NEW_Iteration,GrainElectricCharge_LAST_Iteration;

  PIC::Mesh::cDataCenterNode* cell;
  int i,j,k;
  long int LocalCellNumber;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,initNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
  cell=initNode->block->GetCenterNode(LocalCellNumber);

  //ICES plasma data
  char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
  double *swVel;

  memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+cell->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);
  plasmaTemperature=*((double*)(PIC::CPLR::ICES::PlasmaTemperatureOffset+ICES_AssociatedData));
  swVel=(double*)(PIC::CPLR::ICES::PlasmaBulkVelocityOffset+ICES_AssociatedData);
  plasmaNumberDensity=*((double*)(PIC::CPLR::ICES::PlasmaNumberDensityOffset+ICES_AssociatedData));

  if (plasmaTemperature<=0.0) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;

  //get the grain electric potential
  char localParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double M,GrainRadius,dustPotential;

  memcpy((void*)localParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
  GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)localParticleData);
  GrainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)localParticleData);


  //reserve space for different elecgtron and ion temepratures
  double Ti,Te,J0i,J0e,Je,dJe,Ji,dJi,XiElectron;

  Ti=plasmaTemperature;
  Te=plasmaTemperature;


//==========  DEBUG  BEGIN =================
  static long int nFunctionCalls=0;

  nFunctionCalls++;
/*
  if ((PIC::nTotalThreads==7)&&(nFunctionCalls==1709)) {
    cout << __FILE__ << "@" << __LINE__ << endl;
  }
  */

//==========   DEBUG END ===================


  J0e=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*Te/(PiTimes2*ElectronMass));
  J0i=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*Ti/(PiTimes2*ProtonMass));

  //evaluate the first interation step
  dustPotential=GrainElectricCharge/(4.0*Pi*VacuumPermittivity*GrainRadius);

  //get the electon current and the jacobian of the electron current
  XiElectron=-ElectronCharge*dustPotential/(Kbol*Te);

  if (XiElectron>=0.0) {
    double t2 = 0.1e1 / 0.3141592654e1;
    double t3 = 0.1e1 / VacuumPermittivity;
    double t6 = 0.1e1 / GrainRadius;
    double t7 = 0.1e1 / Kbol;
    double t9 = 0.1e1 / Te;
    double t17 = exp(ElectronCharge * GrainElectricCharge * t2 * t3 * t6 * t7 * t9 / 0.4e1);

    Je=J0e*exp(-XiElectron);
    dJe = J0e * ElectronCharge * t2 * t3 * t6 * t7 * t9 * t17 / 0.4e1;
  }
  else {
    Je=J0e*(1.0-XiElectron);
    dJe=J0e * ElectronCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1;
  }

  //get the current and hacobian of the ion current
  M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*Ti/ProtonMass));


  if (dustPotential<=0.0) {
    {
      double t1 = M * M;
      double t15 = sqrt(0.3141592654e1);
      double t17 = erf(M);
      double t21 = exp(-t1);
      Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - ElectronCharge * GrainElectricCharge / 0.3141592654e1 / VacuumPermittivity / Kbol / Ti / GrainRadius / 0.4e1) * t15 * t17 / M + t21) / 0.2e1;
    }

    {
      double t2 = sqrt(0.3141592654e1);
      double t11 = erf(M);
      dJi = -J0i * ElectronCharge / t2 / VacuumPermittivity / Kbol / Ti / GrainRadius * t11 / M / 0.8e1;
    }
  }
  else {
    {
      double t1 = M * M;
      double t2 = ElectronCharge * GrainElectricCharge;
      double t3 = 0.1e1 / 0.3141592654e1;
      double t5 = 0.1e1 / VacuumPermittivity;
      double t6 = 0.1e1 / Kbol;
      double t8 = 0.1e1 / Ti;
      double t9 = 0.1e1 / GrainRadius;
      double t12 = t2 * t3 * t5 * t6 * t8 * t9;
      double t15 = sqrt(0.3141592654e1);
      double t17 = sqrt(t12);
      double t18 = t17 / 0.2e1;
      double t19 = M + t18;
      double t20 = erf(t19);
      double t21 = M - t18;
      double t22 = erf(t21);
      double t24 = 0.1e1 / M;
      double t33 = sqrt(t2 * t3 * t5 * t6 * t8 * t9 * t24);
      double t34 = t33 / 0.2e1;
      double t36 = t21 * t21;
      double t37 = exp(-t36);
      double t40 = t19 * t19;
      double t41 = exp(-t40);
      Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - t12 / 0.4e1) * t15 * (t20 + t22) * t24 + (t34 + 0.1e1) * t37 - (t34 - 0.1e1) * t41) / 0.2e1;
    }

    {
      double t1 = sqrt(0.3141592654e1);
      double t4 = 0.1e1 / VacuumPermittivity;
      double t5 = 0.1e1 / Kbol;
      double t6 = t4 * t5;
      double t8 = 0.1e1 / Ti;
      double t9 = 0.1e1 / GrainRadius;
      double t10 = t8 * t9;
      double t11 = ElectronCharge * GrainElectricCharge;
      double t12 = 0.1e1 / 0.3141592654e1;
      double t14 = t6 * t10;
      double t15 = t11 * t12 * t14;
      double t16 = sqrt(t15);
      double t17 = t16 / 0.2e1;
      double t18 = M + t17;
      double t19 = erf(t18);
      double t20 = M - t17;
      double t21 = erf(t20);
      double t23 = 0.1e1 / M;
      double t28 = M * M;
      double t33 = 0.1e1 / t1 / 0.3141592654e1;
      double t34 = t18 * t18;
      double t35 = exp(-t34);
      double t38 = 0.1e1 / t16 * ElectronCharge;
      double t41 = t20 * t20;
      double t42 = exp(-t41);
      double t49 = t12 * t4;
      double t51 = t5 * t8;
      double t52 = t9 * t23;
      double t55 = sqrt(t11 * t49 * t51 * t52);
      double t58 = 0.1e1 / t55 * ElectronCharge * t49;
      double t63 = t55 / 0.2e1;
      double t66 = t38 * t12;
      dJi = J0i * (-ElectronCharge / t1 * t6 * t10 * (t19 + t21) * t23 / 0.4e1 + (t28 + 0.1e1 / 0.2e1 - t15 / 0.4e1) * t1 * (t33 * t35 * t38 * t14 - t33 * t42 * t38 * t14) * t23 / 0.2e1 + t58 * t51 * t52 * t42 / 0.4e1 + (t63 + 0.1e1) * t20 * t66 * t6 * t10 * t42 / 0.2e1 - t58 * t51 * t52 * t35 / 0.4e1 + (t63 - 0.1e1) * t18 * t66 * t6 * t10 * t35 / 0.2e1) / 0.2e1;

    }
  }


  double IterationIncrement;
  double dtSubStep=0.0,dtCounter=0.0,dtIons,dtElectrons,totalJ;
  int Counter=0;

  while (dtCounter<dt) {
    dtIons=(fabs(dJi*(Ji+Je)/Ji)*(dt-dtCounter)>0.1) ? 0.1/fabs(dJi) : dt-dtCounter;
    dtElectrons=(fabs(dJe*(Ji+Je)/Je*(dt-dtCounter))>0.1) ? 0.1/fabs(dJe) : dt-dtCounter;

    dtSubStep=min(dt-dtCounter,min(dtIons,dtElectrons));
    dtCounter+=dtSubStep;
    totalJ=Ji+Je;

    GrainElectricCharge_LAST_Iteration=GrainElectricCharge;
    Counter=0;

    //do the iteration loop
    do {
      dustPotential=GrainElectricCharge_LAST_Iteration/(4.0*Pi*VacuumPermittivity*GrainRadius);

      //get the electon current and the jacobian of the electron current
      XiElectron=-ElectronCharge*dustPotential/(Kbol*Te);

      if (XiElectron>=0.0) {
        double t2 = 0.1e1 / 0.3141592654e1;
        double t3 = 0.1e1 / VacuumPermittivity;
        double t6 = 0.1e1 / GrainRadius;
        double t7 = 0.1e1 / Kbol;
        double t9 = 0.1e1 / Te;
        double t17 = exp(ElectronCharge * GrainElectricCharge_LAST_Iteration * t2 * t3 * t6 * t7 * t9 / 0.4e1);

        Je=J0e*exp(-XiElectron);
        dJe = J0e * ElectronCharge * t2 * t3 * t6 * t7 * t9 * t17 / 0.4e1;
      }
      else {
        Je=J0e*(1.0-XiElectron);
        dJe=J0e * ElectronCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1;
      }

      //get the current and hacobian of the ion current
      M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*Ti/ProtonMass));


      if (dustPotential<=0.0) {
        {
          double t1 = M * M;
          double t15 = sqrt(0.3141592654e1);
          double t17 = erf(M);
          double t21 = exp(-t1);
          Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - ElectronCharge * GrainElectricCharge_LAST_Iteration / 0.3141592654e1 / VacuumPermittivity / Kbol / Ti / GrainRadius / 0.4e1) * t15 * t17 / M + t21) / 0.2e1;
        }

        {
          double t2 = sqrt(0.3141592654e1);
          double t11 = erf(M);
          dJi = -J0i * ElectronCharge / t2 / VacuumPermittivity / Kbol / Ti / GrainRadius * t11 / M / 0.8e1;
        }
      }
      else {
        {
          double t1 = M * M;
          double t2 = ElectronCharge * GrainElectricCharge_LAST_Iteration;
          double t3 = 0.1e1 / 0.3141592654e1;
          double t5 = 0.1e1 / VacuumPermittivity;
          double t6 = 0.1e1 / Kbol;
          double t8 = 0.1e1 / Ti;
          double t9 = 0.1e1 / GrainRadius;
          double t12 = t2 * t3 * t5 * t6 * t8 * t9;
          double t15 = sqrt(0.3141592654e1);
          double t17 = sqrt(t12);
          double t18 = t17 / 0.2e1;
          double t19 = M + t18;
          double t20 = erf(t19);
          double t21 = M - t18;
          double t22 = erf(t21);
          double t24 = 0.1e1 / M;
          double t33 = sqrt(t2 * t3 * t5 * t6 * t8 * t9 * t24);
          double t34 = t33 / 0.2e1;
          double t36 = t21 * t21;
          double t37 = exp(-t36);
          double t40 = t19 * t19;
          double t41 = exp(-t40);
          Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - t12 / 0.4e1) * t15 * (t20 + t22) * t24 + (t34 + 0.1e1) * t37 - (t34 - 0.1e1) * t41) / 0.2e1;
        }

        {
          double t1 = sqrt(0.3141592654e1);
          double t4 = 0.1e1 / VacuumPermittivity;
          double t5 = 0.1e1 / Kbol;
          double t6 = t4 * t5;
          double t8 = 0.1e1 / Ti;
          double t9 = 0.1e1 / GrainRadius;
          double t10 = t8 * t9;
          double t11 = ElectronCharge * GrainElectricCharge_LAST_Iteration;
          double t12 = 0.1e1 / 0.3141592654e1;
          double t14 = t6 * t10;
          double t15 = t11 * t12 * t14;
          double t16 = sqrt(t15);
          double t17 = t16 / 0.2e1;
          double t18 = M + t17;
          double t19 = erf(t18);
          double t20 = M - t17;
          double t21 = erf(t20);
          double t23 = 0.1e1 / M;
          double t28 = M * M;
          double t33 = 0.1e1 / t1 / 0.3141592654e1;
          double t34 = t18 * t18;
          double t35 = exp(-t34);
          double t38 = 0.1e1 / t16 * ElectronCharge;
          double t41 = t20 * t20;
          double t42 = exp(-t41);
          double t49 = t12 * t4;
          double t51 = t5 * t8;
          double t52 = t9 * t23;
          double t55 = sqrt(t11 * t49 * t51 * t52);
          double t58 = 0.1e1 / t55 * ElectronCharge * t49;
          double t63 = t55 / 0.2e1;
          double t66 = t38 * t12;
          dJi = J0i * (-ElectronCharge / t1 * t6 * t10 * (t19 + t21) * t23 / 0.4e1 + (t28 + 0.1e1 / 0.2e1 - t15 / 0.4e1) * t1 * (t33 * t35 * t38 * t14 - t33 * t42 * t38 * t14) * t23 / 0.2e1 + t58 * t51 * t52 * t42 / 0.4e1 + (t63 + 0.1e1) * t20 * t66 * t6 * t10 * t42 / 0.2e1 - t58 * t51 * t52 * t35 / 0.4e1 + (t63 - 0.1e1) * t18 * t66 * t6 * t10 * t35 / 0.2e1) / 0.2e1;
        }
      }

      IterationIncrement=((GrainElectricCharge+totalJ*dtSubStep/2.0)-(GrainElectricCharge_LAST_Iteration-(Ji+Je)*dtSubStep/2.0))/(1.0-(dJi+dJe)*dtSubStep/2.0);
      GrainElectricCharge_NEW_Iteration=GrainElectricCharge_LAST_Iteration+IterationIncrement;


      //evaluate the convergance of the iterations
      if (fabs(IterationIncrement)<1.0E-4*fabs(GrainElectricCharge_LAST_Iteration+GrainElectricCharge_NEW_Iteration)) {
        GrainElectricCharge=GrainElectricCharge_NEW_Iteration;
        break;
      }

      GrainElectricCharge_LAST_Iteration=GrainElectricCharge_NEW_Iteration;
    }
    while (++Counter<10000);
  }


  SetGrainCharge(GrainElectricCharge,ParticleData);

  //move the particle into diferent velocity group if needed
  int oldVelocityGroup,newVelocityGroup;

  oldVelocityGroup=spec-_DUST_SPEC_;
  newVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);

  if (oldVelocityGroup!=newVelocityGroup) {
    //move the particle into different velocity group
    double GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

    GrainWeightCorrection*=initNode->block->GetLocalTimeStep(_DUST_SPEC_+newVelocityGroup)/initNode->block->GetLocalTimeStep(_DUST_SPEC_+oldVelocityGroup);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+oldVelocityGroup]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+newVelocityGroup];
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,ParticleData);
    spec=_DUST_SPEC_+newVelocityGroup;
    PIC::ParticleBuffer::SetI(spec,ParticleData);
  }



  return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
}




  /*-----------------------  Steady State Dust Charging ------------------------------- */
  inline int DustChargingProcessor_SteadyState(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {
    double plasmaTemperature,plasmaNumberDensity;
    double GrainElectricCharge,GrainElectricCharge_NEW;

    PIC::Mesh::cDataCenterNode* cell;
    int i,j,k;
    long int LocalCellNumber;

    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,initNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
    cell=initNode->block->GetCenterNode(LocalCellNumber);

    //ICES plasma data
    char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
    double *swVel;

    memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+cell->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);
    plasmaTemperature=*((double*)(PIC::CPLR::ICES::PlasmaTemperatureOffset+ICES_AssociatedData));
    swVel=(double*)(PIC::CPLR::ICES::PlasmaBulkVelocityOffset+ICES_AssociatedData);
    plasmaNumberDensity=*((double*)(PIC::CPLR::ICES::PlasmaNumberDensityOffset+ICES_AssociatedData));

    if (plasmaTemperature<=0.0) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;

    //get the grain electric potential
    char localParticleData[PIC::ParticleBuffer::ParticleDataLength];
    double M,GrainRadius,dustPotential;

    memcpy((void*)localParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
    GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)localParticleData);
    GrainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)localParticleData);


    //reserve space for different elecgtron and ion temepratures
    double Ti,Te,J0i,J0e,Je,dJe,Ji,dJi,XiElectron;

    Ti=plasmaTemperature;
    Te=plasmaTemperature;


  //==========  DEBUG  BEGIN =================
    static long int nFunctionCalls=0;

    nFunctionCalls++;

    /*
    if ((PIC::nTotalThreads==7)&&(nFunctionCalls==1709)) {
      cout << __FILE__ << "@" << __LINE__ << endl;
    }
    */

  //==========   DEBUG END ===================


    J0e=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*Te/(PiTimes2*ElectronMass));
    J0i=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*Ti/(PiTimes2*ProtonMass));

    //evaluate the first interation step
    dustPotential=GrainElectricCharge/(4.0*Pi*VacuumPermittivity*GrainRadius);

    //get the electon current and the jacobian of the electron current
    XiElectron=-ElectronCharge*dustPotential/(Kbol*Te);

    if (XiElectron>=0.0) {
      double t2 = 0.1e1 / 0.3141592654e1;
      double t3 = 0.1e1 / VacuumPermittivity;
      double t6 = 0.1e1 / GrainRadius;
      double t7 = 0.1e1 / Kbol;
      double t9 = 0.1e1 / Te;
      double t17 = exp(ElectronCharge * GrainElectricCharge * t2 * t3 * t6 * t7 * t9 / 0.4e1);

      Je=J0e*exp(-XiElectron);
      dJe = J0e * ElectronCharge * t2 * t3 * t6 * t7 * t9 * t17 / 0.4e1;
    }
    else {
      Je=J0e*(1.0-XiElectron);
      dJe=J0e * ElectronCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1;
    }

    //get the current and hacobian of the ion current
    M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*Ti/ProtonMass));


    if (dustPotential<=0.0) {
      {
        double t1 = M * M;
        double t15 = sqrt(0.3141592654e1);
        double t17 = erf(M);
        double t21 = exp(-t1);
        Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - ElectronCharge * GrainElectricCharge / 0.3141592654e1 / VacuumPermittivity / Kbol / Ti / GrainRadius / 0.4e1) * t15 * t17 / M + t21) / 0.2e1;
      }

      {
        double t2 = sqrt(0.3141592654e1);
        double t11 = erf(M);
        dJi = -J0i * ElectronCharge / t2 / VacuumPermittivity / Kbol / Ti / GrainRadius * t11 / M / 0.8e1;
      }
    }
    else {
      {
        double t1 = M * M;
        double t2 = ElectronCharge * GrainElectricCharge;
        double t3 = 0.1e1 / 0.3141592654e1;
        double t5 = 0.1e1 / VacuumPermittivity;
        double t6 = 0.1e1 / Kbol;
        double t8 = 0.1e1 / Ti;
        double t9 = 0.1e1 / GrainRadius;
        double t12 = t2 * t3 * t5 * t6 * t8 * t9;
        double t15 = sqrt(0.3141592654e1);
        double t17 = sqrt(t12);
        double t18 = t17 / 0.2e1;
        double t19 = M + t18;
        double t20 = erf(t19);
        double t21 = M - t18;
        double t22 = erf(t21);
        double t24 = 0.1e1 / M;
        double t33 = sqrt(t2 * t3 * t5 * t6 * t8 * t9 * t24);
        double t34 = t33 / 0.2e1;
        double t36 = t21 * t21;
        double t37 = exp(-t36);
        double t40 = t19 * t19;
        double t41 = exp(-t40);
        Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - t12 / 0.4e1) * t15 * (t20 + t22) * t24 + (t34 + 0.1e1) * t37 - (t34 - 0.1e1) * t41) / 0.2e1;
      }

      {
        double t1 = sqrt(0.3141592654e1);
        double t4 = 0.1e1 / VacuumPermittivity;
        double t5 = 0.1e1 / Kbol;
        double t6 = t4 * t5;
        double t8 = 0.1e1 / Ti;
        double t9 = 0.1e1 / GrainRadius;
        double t10 = t8 * t9;
        double t11 = ElectronCharge * GrainElectricCharge;
        double t12 = 0.1e1 / 0.3141592654e1;
        double t14 = t6 * t10;
        double t15 = t11 * t12 * t14;
        double t16 = sqrt(t15);
        double t17 = t16 / 0.2e1;
        double t18 = M + t17;
        double t19 = erf(t18);
        double t20 = M - t17;
        double t21 = erf(t20);
        double t23 = 0.1e1 / M;
        double t28 = M * M;
        double t33 = 0.1e1 / t1 / 0.3141592654e1;
        double t34 = t18 * t18;
        double t35 = exp(-t34);
        double t38 = 0.1e1 / t16 * ElectronCharge;
        double t41 = t20 * t20;
        double t42 = exp(-t41);
        double t49 = t12 * t4;
        double t51 = t5 * t8;
        double t52 = t9 * t23;
        double t55 = sqrt(t11 * t49 * t51 * t52);
        double t58 = 0.1e1 / t55 * ElectronCharge * t49;
        double t63 = t55 / 0.2e1;
        double t66 = t38 * t12;
        dJi = J0i * (-ElectronCharge / t1 * t6 * t10 * (t19 + t21) * t23 / 0.4e1 + (t28 + 0.1e1 / 0.2e1 - t15 / 0.4e1) * t1 * (t33 * t35 * t38 * t14 - t33 * t42 * t38 * t14) * t23 / 0.2e1 + t58 * t51 * t52 * t42 / 0.4e1 + (t63 + 0.1e1) * t20 * t66 * t6 * t10 * t42 / 0.2e1 - t58 * t51 * t52 * t35 / 0.4e1 + (t63 - 0.1e1) * t18 * t66 * t6 * t10 * t35 / 0.2e1) / 0.2e1;

      }
    }


    double IterationIncrement,InitIterationIncrement=-(Ji+Je)/(dJi+dJe),IterationParameter=1.0;
    int Counter=0;

    GrainElectricCharge_NEW=GrainElectricCharge+InitIterationIncrement;


    //do the iteration loop
    do {

      GrainElectricCharge_NEW=GrainElectricCharge+IterationParameter*InitIterationIncrement;

      dustPotential=GrainElectricCharge_NEW/(4.0*Pi*VacuumPermittivity*GrainRadius);

      //get the electon current and the jacobian of the electron current
      XiElectron=-ElectronCharge*dustPotential/(Kbol*Te);

      if (XiElectron>=0.0) {
        double t2 = 0.1e1 / 0.3141592654e1;
        double t3 = 0.1e1 / VacuumPermittivity;
        double t6 = 0.1e1 / GrainRadius;
        double t7 = 0.1e1 / Kbol;
        double t9 = 0.1e1 / Te;
        double t17 = exp(ElectronCharge * GrainElectricCharge_NEW * t2 * t3 * t6 * t7 * t9 / 0.4e1);

        Je=J0e*exp(-XiElectron);
        dJe = J0e * ElectronCharge * t2 * t3 * t6 * t7 * t9 * t17 / 0.4e1;
      }
      else {
        Je=J0e*(1.0-XiElectron);
        dJe=J0e * ElectronCharge / 0.3141592654e1 / VacuumPermittivity / GrainRadius / Kbol / Te / 0.4e1;
      }

      //get the current and hacobian of the ion current
      M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*Ti/ProtonMass));


      if (dustPotential<=0.0) {
        {
          double t1 = M * M;
          double t15 = sqrt(0.3141592654e1);
          double t17 = erf(M);
          double t21 = exp(-t1);
          Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - ElectronCharge * GrainElectricCharge_NEW / 0.3141592654e1 / VacuumPermittivity / Kbol / Ti / GrainRadius / 0.4e1) * t15 * t17 / M + t21) / 0.2e1;
        }

        {
          double t2 = sqrt(0.3141592654e1);
          double t11 = erf(M);
          dJi = -J0i * ElectronCharge / t2 / VacuumPermittivity / Kbol / Ti / GrainRadius * t11 / M / 0.8e1;
        }
      }
      else {
        {
          double t1 = M * M;
          double t2 = ElectronCharge * GrainElectricCharge_NEW;
          double t3 = 0.1e1 / 0.3141592654e1;
          double t5 = 0.1e1 / VacuumPermittivity;
          double t6 = 0.1e1 / Kbol;
          double t8 = 0.1e1 / Ti;
          double t9 = 0.1e1 / GrainRadius;
          double t12 = t2 * t3 * t5 * t6 * t8 * t9;
          double t15 = sqrt(0.3141592654e1);
          double t17 = sqrt(t12);
          double t18 = t17 / 0.2e1;
          double t19 = M + t18;
          double t20 = erf(t19);
          double t21 = M - t18;
          double t22 = erf(t21);
          double t24 = 0.1e1 / M;
          double t33 = sqrt(t2 * t3 * t5 * t6 * t8 * t9 * t24);
          double t34 = t33 / 0.2e1;
          double t36 = t21 * t21;
          double t37 = exp(-t36);
          double t40 = t19 * t19;
          double t41 = exp(-t40);
          Ji = J0i * ((t1 + 0.1e1 / 0.2e1 - t12 / 0.4e1) * t15 * (t20 + t22) * t24 + (t34 + 0.1e1) * t37 - (t34 - 0.1e1) * t41) / 0.2e1;
        }

        {
          double t1 = sqrt(0.3141592654e1);
          double t4 = 0.1e1 / VacuumPermittivity;
          double t5 = 0.1e1 / Kbol;
          double t6 = t4 * t5;
          double t8 = 0.1e1 / Ti;
          double t9 = 0.1e1 / GrainRadius;
          double t10 = t8 * t9;
          double t11 = ElectronCharge * GrainElectricCharge_NEW;
          double t12 = 0.1e1 / 0.3141592654e1;
          double t14 = t6 * t10;
          double t15 = t11 * t12 * t14;
          double t16 = sqrt(t15);
          double t17 = t16 / 0.2e1;
          double t18 = M + t17;
          double t19 = erf(t18);
          double t20 = M - t17;
          double t21 = erf(t20);
          double t23 = 0.1e1 / M;
          double t28 = M * M;
          double t33 = 0.1e1 / t1 / 0.3141592654e1;
          double t34 = t18 * t18;
          double t35 = exp(-t34);
          double t38 = 0.1e1 / t16 * ElectronCharge;
          double t41 = t20 * t20;
          double t42 = exp(-t41);
          double t49 = t12 * t4;
          double t51 = t5 * t8;
          double t52 = t9 * t23;
          double t55 = sqrt(t11 * t49 * t51 * t52);
          double t58 = 0.1e1 / t55 * ElectronCharge * t49;
          double t63 = t55 / 0.2e1;
          double t66 = t38 * t12;
          dJi = J0i * (-ElectronCharge / t1 * t6 * t10 * (t19 + t21) * t23 / 0.4e1 + (t28 + 0.1e1 / 0.2e1 - t15 / 0.4e1) * t1 * (t33 * t35 * t38 * t14 - t33 * t42 * t38 * t14) * t23 / 0.2e1 + t58 * t51 * t52 * t42 / 0.4e1 + (t63 + 0.1e1) * t20 * t66 * t6 * t10 * t42 / 0.2e1 - t58 * t51 * t52 * t35 / 0.4e1 + (t63 - 0.1e1) * t18 * t66 * t6 * t10 * t35 / 0.2e1) / 0.2e1;
        }
      }

      IterationIncrement=-(Ji+Je)/(dJi+dJe);

      if (IterationIncrement*InitIterationIncrement<0.0) {
        //the iteration resutled in overshoot of the dust charge
        IterationParameter/=2.0;
        continue;
      }
      else {

        //evaluate the convergance of the iterations
        if (fabs(GrainElectricCharge-GrainElectricCharge_NEW)<1.0E-4*fabs(GrainElectricCharge+GrainElectricCharge_NEW)) {
          GrainElectricCharge=GrainElectricCharge_NEW;
          break;
        }

        GrainElectricCharge=GrainElectricCharge_NEW;
        InitIterationIncrement=IterationIncrement;

        IterationParameter*=2.0;
        if (IterationParameter>1.0) IterationParameter=1.0;
      }

    }
    while (++Counter<10000);


    SetGrainCharge(GrainElectricCharge,ParticleData);

    //move the particle into diferent velocity group if needed
    int oldVelocityGroup,newVelocityGroup;

    oldVelocityGroup=spec-_DUST_SPEC_;
    newVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);

    if (oldVelocityGroup!=newVelocityGroup) {
      //move the particle into different velocity group
      double GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

      GrainWeightCorrection*=initNode->block->GetLocalTimeStep(_DUST_SPEC_+newVelocityGroup)/initNode->block->GetLocalTimeStep(_DUST_SPEC_+oldVelocityGroup);

  #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
      GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+oldVelocityGroup]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+newVelocityGroup];
  #else
      exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
  #endif


      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,ParticleData);
      spec=_DUST_SPEC_+newVelocityGroup;
      PIC::ParticleBuffer::SetI(spec,ParticleData);
    }



    return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
  }


//injection of the dust grains
typedef bool (*fGenerateNewDustGrainInternalProperties)(double *x,double *v,double& GrainRadius, double& GrainWeightCorrection,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode);
extern fGenerateNewDustGrainInternalProperties GenerateNewDustGrainInternalProperties;

long int DustInjection__Sphere(void *SphereDataPointer);
}





#define ELECTRICALLYCHARGEDDUST_H_


#endif /* ELECTRICALLYCHARGEDDUST_H_ */
