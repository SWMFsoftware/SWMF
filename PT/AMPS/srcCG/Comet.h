/*
 * Mercury.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

#ifndef COMET_H_
#define COMET_H_

#include <vector>

#include "Exosphere.h"
#include "CG-BC.h"
#include "Gravity.h"

#include "Comet.dfn"

namespace Comet {
  using namespace Exosphere;

  void PrintMaxLiftableSizeSurfaceTriangulationMesh(const char *fname);

  double GetTotalProduction(int spec,void *BoundaryElement);
  double GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement);

  //init the model
  void Init_BeforeParser();
  void Init_AfterParser();
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
#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
      //the gravity force non spherical case
      Comet::GetGravityAcceleration(accl_LOCAL,nd,startNode);

#else
      //the gravity force spherical case
      double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
      double r=sqrt(r2);
      double mass=1.0e13;

      for (idim=0;idim<DIM;idim++) {
	accl_LOCAL[idim]-=GravityConstant*mass/r2*x_LOCAL[idim]/r;
      }
#endif
      
      //Drag force
      char ParticleData[PIC::ParticleBuffer::ParticleDataLength];
      
      memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);
      
      double GrainDragCoefficient=2.0;
      double A,cr2;
      double GasBulkVelocity[3];
      double GasNumberDensity;
      double GrainRadius=ElectricallyChargedDust::GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);
      double GrainMass=ElectricallyChargedDust::GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);
     
      for (int s=0;s<PIC::nTotalSpecies;s++) if (_DUST_SPEC_>s || s>=_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
	  startNode->block->GetCenterNode(nd)->GetBulkVelocity(GasBulkVelocity,s);
	  GasNumberDensity=startNode->block->GetCenterNode(nd)->GetNumberDensity(s);
	  
	  cr2=(v_LOCAL[0]-GasBulkVelocity[0])*(v_LOCAL[0]-GasBulkVelocity[0])+                                                                                         
	    (v_LOCAL[1]-GasBulkVelocity[1])*(v_LOCAL[1]-GasBulkVelocity[1])+                                                                                           
	    (v_LOCAL[2]-GasBulkVelocity[2])*(v_LOCAL[2]-GasBulkVelocity[2]);                                                                                           
	  
	  A=Pi*pow(GrainRadius,2)/2.0*GrainDragCoefficient*sqrt(cr2)/GrainMass*GasNumberDensity*PIC::MolecularData::GetMass(s); //_MASS_(_H2O_);
	  

	  accl_LOCAL[0]+=A*(GasBulkVelocity[0]-v_LOCAL[0]);                                                                                                                             
	  accl_LOCAL[1]+=A*(GasBulkVelocity[1]-v_LOCAL[1]);                                                                                                                             
	  accl_LOCAL[2]+=A*(GasBulkVelocity[2]-v_LOCAL[2]);                                                                                                                             
	}
      //Generation of the tangential force corresponding to diffusion                                                                                                                               
      double e0[3],e1[3]={0.0,0.0,0.0},e2[3]={0.0,0.0,0.0};
      double norm=0,drag=0,fraction=0.2;
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

      e1[0]=-e0[1];
      e1[1]=e0[0];

      for (idim=0;idim<DIM;idim++) norm+=pow(e1[idim],2.0);
      norm=sqrt(norm);
      for (idim=0;idim<DIM;idim++) e1[idim]/=norm;

      e2[0]=-e0[2]*e1[1];
      e2[1]=e0[2]*e1[0];
      e2[2]=e0[0]*e1[1]-e0[1]*e1[0];

      for (idim=0;idim<DIM;idim++) accl_LOCAL[idim]+=drag*fraction*(cos(theta)*e1[idim]+sin(theta)*e2[idim]);
      
      /*      //the Lorentz force
      char *offset;
      double E[3],B[3];
      double GrainCharge=ElectricallyChargedDust::GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);

      offset=startNode->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer();
      
      if (*((int*)(offset+PIC::CPLR::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
        memcpy(E,offset+PIC::CPLR::ICES::ElectricFieldOffset,3*sizeof(double));
        memcpy(B,offset+PIC::CPLR::ICES::MagneticFieldOffset,3*sizeof(double));
      }
      else {
        memcpy(E,swE_Typical,3*sizeof(double));
        memcpy(B,Exosphere_swB_Typical,3*sizeof(double));
      }
            
      accl_LOCAL[0]+=GrainCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/GrainMass;
      accl_LOCAL[1]+=GrainCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/GrainMass;
      accl_LOCAL[2]+=GrainCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/GrainMass;
      */
    }  
#endif
           
  //Test: no acceleration:
    memcpy(accl,accl_LOCAL,3*sizeof(double));
    return;}
      
#endif

