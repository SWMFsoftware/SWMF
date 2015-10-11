/*
 * Comet.h
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and  vtenishe
 */

#ifndef COMET_H_
#define COMET_H_

#include <vector>

#include "PhotolyticReactions.h"
#include "Exosphere.h"
#include "CG-BC.h"
#include "Gravity.h"

#include "Comet.dfn"

namespace Comet {
  using namespace Exosphere;

  unsigned int GetParticleSurfaceElement(PIC::ParticleBuffer::byte *ParticleDataStart);
  void SetParticleSurfaceElement(int SurfaceElement,PIC::ParticleBuffer::byte *ParticleDataStart);

  void PrintMaxLiftableSizeSurfaceTriangulationMesh(const char *fname);

  double GetTotalProduction(int spec,void *BoundaryElement);
  double GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement);

  //init the model
  void Init_BeforeParser();
  void Init_AfterParser(const char*);
  void InitGravityData();

  static int ndist=4;
  static double Bjorn_SourceRate[]={0.0};
  static double Jet_SourceRate[]={0.0};
  static double Uniform_SourceRate[]={0.0};

  extern int  GravityFieldOffset;

  static cInternalNastranSurfaceData *CG;
  void GetNucleusNastranInfo(cInternalNastranSurfaceData *CG);

  void PrintSurfaceTriangulationMesh(const char *fname,CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS);
  
  double GetTotalProductionRateBjornNASTRAN(int spec);
  bool GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData);

  double GetTotalProductionRateUniformNASTRAN(int spec);
  bool GenerateParticlePropertiesUniformNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData);

  double GetTotalProductionRateJetNASTRAN(int spec);
  bool GenerateParticlePropertiesJetNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData);

  long int InjectionBoundaryModel_Limited();
  long int InjectionBoundaryModel_Limited(int spec);

  double radiativeCoolingRate_Crovisier(PIC::Mesh::cDataCenterNode *CenterNode);
  void StepOverTime();

  int RequestDataBuffer(int offset);
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

  void GetGravityAcceleration(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  namespace LossProcesses {
   extern double PhotolyticReactionRate,ElectronImpactRate,ElectronTemeprature;

   //the constant of the artificial increase of the primary species loss
   //the modification of the rate is compensated by the appropricate particle weight of the daugher products and
   //probabilities of the destruction of the primary species
   const double NumericalLossRateIncrease=1000.0; //1000.0

   double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
   int ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
  }

  //the condition for begining of the dust trajectory tracking
  namespace TrajectoryTracking {
    const int nZenithSurfaceElements=300,nAzimuthalSurfaceElements=300;
    const int nTotalTracedTrajectories=25000;
    const double TracingSurfaceRadius=2.7e3;

    extern int nTracedTrajectoriesPerElement;
    extern int TracedParticleNumber_GLOBAL[nZenithSurfaceElements*nAzimuthalSurfaceElements],TracedParticleNumber_LOCAL[nZenithSurfaceElements*nAzimuthalSurfaceElements];
    extern cInternalSphericalData Sphere;

    void Init();
    bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData);
    void UpdateParticleCounter();
  }

  namespace CometData {
    extern int NeutralsFromBinaryOffset;
    extern int nNeutrals;
    extern int iSpecies;

    void WriteBinaryOutput(const char *fNameBase,int s,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,FILE *fout);
    void LoadBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,FILE *fout);
    int LoadSpeciesNumberBinaryFile(const char *fNameBase,FILE *fout);
    int RequestDataBuffer(int offset);
    void PrintVariableList(FILE* fout,int DataSetNumber);
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

    double GetNeutralsMassDensity(int s,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    void GetNeutralsVelocity(double *x, int s,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    int GetnNeutrals();
    void SetnNeutrals(int n);

    int GetiSpecies();
    void SetiSpecies(int s);
  }

  //Sampling: sample quantaty a^2*Weight for the further integration of the brightness of the dust grains
  namespace Sampling {
     //sample particle data
     extern int SamplingDataOffset;

     //sample particle data
     void SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec);
     int RequestSamplingData(int offset);

     //empty model sampling function and function output of the column integrals
     void SampleModelData();
     void PrintBrightnessMap(double halfAngleRange,int iTestPoint,int DataOutputFileNumber);
     void PrintBrightnessMap(int DataOutputFileNumber);
  }

}
    //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};
    int idim;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      int nd,i,j,k;

      memcpy(x_LOCAL,x,3*sizeof(double));
      memcpy(v_LOCAL,v,3*sizeof(double));
   
      nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode);
      
      //the gravity force
      if (_CG_DUST_FORCE_MODE__GRAVITY_ == _PIC_MODE_ON_) {
        #if _3DGRAVITY__MODE_ == _3DGRAVITY__MODE__ON_
        //the gravity force non spherical case
        Comet::GetGravityAcceleration(accl_LOCAL,nd,startNode);

        #else
        //the gravity force spherical case
        double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
        double r=sqrt(r2);
        double mass=1.0e13;

         for (idim=0;idim<DIM;idim++) accl_LOCAL[idim]-=GravityConstant*mass/r2*x_LOCAL[idim]/r;
         #endif
      }
      
      //Drag force
      char ParticleData[PIC::ParticleBuffer::ParticleDataLength];
      
      memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);
      
      double GrainDragCoefficient=2.0;
      double A,cr2;
      double GasBulkVelocity[3];
      double GasNumberDensity;
      double GrainRadius=ElectricallyChargedDust::GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);
      double GrainMass=ElectricallyChargedDust::GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);
     
      if (_CG_DUST_FORCE_MODE__DRAG_FORCE_ == _PIC_MODE_ON_) for (int s=0;s<PIC::nTotalSpecies;s++) if (_DUST_SPEC_>s || s>=_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
	      startNode->block->GetCenterNode(nd)->GetBulkVelocity(GasBulkVelocity,s);
	      GasNumberDensity=startNode->block->GetCenterNode(nd)->GetNumberDensity(s);
	  
	      cr2=(v_LOCAL[0]-GasBulkVelocity[0])*(v_LOCAL[0]-GasBulkVelocity[0])+
	       (v_LOCAL[1]-GasBulkVelocity[1])*(v_LOCAL[1]-GasBulkVelocity[1])+
	       (v_LOCAL[2]-GasBulkVelocity[2])*(v_LOCAL[2]-GasBulkVelocity[2]);
	  
	      A=Pi*pow(GrainRadius,2)/2.0*GrainDragCoefficient*sqrt(cr2)/GrainMass*GasNumberDensity*PIC::MolecularData::GetMass(s); //_MASS_(_H2O_);
	  
        accl_LOCAL[0]+=A*(GasBulkVelocity[0]-v_LOCAL[0]);
        accl_LOCAL[1]+=A*(GasBulkVelocity[1]-v_LOCAL[1]);
        accl_LOCAL[2]+=A*(GasBulkVelocity[2]-v_LOCAL[2]);

        //Generation of the tangential force corresponding to diffusion
        double e0[3],e1[3]={0.0,0.0,0.0},e2[3]={0.0,0.0,0.0};
        double norm=0,drag=0,fraction=0.02;
        double theta=2*Pi*rnd();

        for (idim=0;idim<DIM;idim++) {
          e0[idim]=GasBulkVelocity[idim]-v_LOCAL[idim];
          norm+=pow(e0[idim],2.0);
          drag+=pow(A*(GasBulkVelocity[idim]-v_LOCAL[idim]),2.0);
        }
        norm=sqrt(norm);
        drag=sqrt(drag);

        for (idim=0;idim<DIM;idim++) e0[idim]/=norm;

        norm=0.0;

        if (fabs(e0[0])>1.0E-5) e1[0]=-e0[1],e1[1]=e0[0],e1[2]=0.0;
        else e1[0]=0.0,e1[1]=e0[2],e1[2]=-e0[1];

        double l=sqrt((e1[0]*e1[0])+(e1[1]*e1[1])+(e1[2]*e1[2]));
        e1[0]/=l,e1[1]/=l,e1[2]/=l;

        e2[0]=+(e0[1]*e1[2]-e1[1]*e0[2]);
        e2[1]=-(e0[0]*e1[2]-e1[0]*e0[2]);
        e2[2]=+(e0[0]*e1[1]-e1[0]*e0[1]);

        if (_CG_DUST_FORCE_MODE__DRAG_FORCE__TANGENTIAL_COMPONENT_ == _PIC_MODE_ON_) for (idim=0;idim<DIM;idim++) accl_LOCAL[idim]+=drag*fraction*(cos(theta)*e1[idim]+sin(theta)*e2[idim]);
      }
      
      //the Lorentz force
      double E[3],B[3];
      double GrainCharge=ElectricallyChargedDust::GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);

      if (_CG_DUST_FORCE_MODE__LORENTZ_FORCE_ == _PIC_MODE_ON_) if (fabs(GrainCharge)>0.0) {
        PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
        PIC::CPLR::GetBackgroundMagneticField(B);
        PIC::CPLR::GetBackgroundElectricField(E);

        accl_LOCAL[0]+=GrainCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/GrainMass;
        accl_LOCAL[1]+=GrainCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/GrainMass;
        accl_LOCAL[2]+=GrainCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/GrainMass;
      }

      //Radiation Pressure: beta is taken from fig3-Gombosi-2015-AA
      if (_CG_DUST_FORCE_MODE__RADIATION_PRESSURE_ == _PIC_MODE_ON_) {
        static const int n_Gombosi2015AA=101;
        static const double dR_Gombosi2015AA=0.01220*1.0E-9;
        static const double Rmin_Gombosi2015AA=0.10000*1.0E-9;
        static const double Rmax_Gombosi2015AA=10000.00000*1.0E-9;
        static const double beta_Gombosi2015AA[n_Gombosi2015AA]={0.00191, 0.00191, 0.00191, 0.00387, 0.00447, 0.00584, 0.00584, 0.00721, 0.00721, 0.00858, 0.01055, 0.01261, 0.01269,
            0.01474, 0.01671, 0.01671, 0.01868, 0.02064, 0.02261, 0.02603, 0.02800, 0.03167, 0.03698, 0.03894, 0.04366, 0.04913, 0.05520, 0.06204, 0.06812, 0.07556, 0.08420, 0.09506,
            0.10507, 0.11585, 0.13031, 0.14280, 0.16246, 0.18180, 0.20079, 0.21995, 0.24518, 0.27101, 0.30276, 0.33637, 0.37367, 0.40754, 0.44491, 0.48520, 0.53225, 0.58169, 0.62814,
            0.67348, 0.72729, 0.77827, 0.82257, 0.87236, 0.90785, 0.94814, 0.97628, 0.99528, 0.99998, 0.99339, 0.96816, 0.93565, 0.88398, 0.81264, 0.72685, 0.65542, 0.58947, 0.51343,
            0.46800, 0.41677, 0.37006, 0.32601, 0.29487, 0.26305, 0.23585, 0.20942, 0.18572, 0.16323, 0.14629, 0.13201, 0.11600, 0.10394, 0.09325, 0.08513, 0.07409, 0.06537, 0.05920,
            0.05245, 0.04706, 0.04107, 0.03825, 0.03354, 0.03132, 0.02713, 0.02448, 0.02106, 0.01901, 0.01695, 0.01413};

        static const double log10Rmin_Gombosi2015AA=log10(Rmin_Gombosi2015AA);
        static const double log10Rmax_Gombosi2015AA=log10(Rmax_Gombosi2015AA);
        static const double dLog10R=(log10Rmax_Gombosi2015AA-log10Rmin_Gombosi2015AA)/(n_Gombosi2015AA-1);

        double beta;

        if (GrainRadius<Rmin_Gombosi2015AA) beta=0.00191;
        else if (GrainRadius<Rmax_Gombosi2015AA) {
          int i;
          double c;

  //        c=(GrainRadius-Rmin_Gombosi2015AA)/dR_Gombosi2015AA;

          c=(log10(GrainRadius)-log10Rmin_Gombosi2015AA)/dLog10R;
          i=(int)c;
          c-=i;

          beta=(1.0-c)*beta_Gombosi2015AA[i]+c*beta_Gombosi2015AA[i+1];
        }
        else beta=0.01413;

        double SunDirection[3]={1.0,0.0,0.0};
        double HeliocentricDistance=1.5*_AU_;
        double RadiationPressure;

        RadiationPressure=beta*GravityConstant*_SUN__MASS_/pow(HeliocentricDistance,2);

        //the total acceleration due to Sun's radiation pressure and gravity is RadiationPressure + SunGravity - Acceleration of the frame of the refecence (Acceleration of the nucleus due to Sun's gravity
        for (int idim=0;idim<3;idim++) accl_LOCAL[idim]-=SunDirection[idim]*RadiationPressure;
      }


     //Centrifugal and Centripital forces
      if (_CG_DUST_FORCE_MODE__FRAME_ROTATION_ == _PIC_MODE_ON_) {
        double aCen[3],aCorr[3],t3,t7,t12,RotationVector_SO_FROZEN[3];

        static const double RotationAxis[3]={0.0,0.0,1.0};
        static const double RotationRate=2.0*Pi/(12.0*3600.0);

        for (int idim=0;idim<3;idim++) RotationVector_SO_FROZEN[idim]=RotationRate*RotationAxis[idim];

        t3 = RotationVector_SO_FROZEN[0] * x_LOCAL[1] - RotationVector_SO_FROZEN[1] * x_LOCAL[0];
        t7 = RotationVector_SO_FROZEN[2] * x_LOCAL[0] - RotationVector_SO_FROZEN[0] * x_LOCAL[2];
        t12 = RotationVector_SO_FROZEN[1] * x_LOCAL[2] - RotationVector_SO_FROZEN[2] * x_LOCAL[1];

        aCen[0] = -RotationVector_SO_FROZEN[1] * t3 + RotationVector_SO_FROZEN[2] * t7;
        aCen[1] = -RotationVector_SO_FROZEN[2] * t12 + RotationVector_SO_FROZEN[0] * t3;
        aCen[2] = -RotationVector_SO_FROZEN[0] * t7 + RotationVector_SO_FROZEN[1] * t12;


        aCorr[0] = -2.0*(RotationVector_SO_FROZEN[1] * v_LOCAL[2] - RotationVector_SO_FROZEN[2] * v_LOCAL[1]);
        aCorr[1] = -2.0*(RotationVector_SO_FROZEN[2] * v_LOCAL[0] - RotationVector_SO_FROZEN[0] * v_LOCAL[2]);
        aCorr[2] = -2.0*(RotationVector_SO_FROZEN[0] * v_LOCAL[1] - RotationVector_SO_FROZEN[1] * v_LOCAL[0]);

        accl_LOCAL[0]+=aCen[0]+aCorr[0];
        accl_LOCAL[1]+=aCen[1]+aCorr[1];
        accl_LOCAL[2]+=aCen[2]+aCorr[2];
      }

      #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      #if _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_
      PIC::Debugger::CatchOutLimitValue(accl_LOCAL,DIM,__LINE__,__FILE__);
      #endif
      #endif

    }  
#endif
           
  //Test: no acceleration:
    memcpy(accl,accl_LOCAL,3*sizeof(double));
}
      
#endif

