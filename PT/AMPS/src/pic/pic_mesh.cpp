//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the function for particle's data sampling

#include "pic.h"
namespace PIC{
  namespace Mesh{
    //basic macroscopic parameters sampled in the simulation
    cDatumTimed    DatumParticleWeight(1,"\"Particle Weight\"",      false);
    cDatumTimed    DatumParticleNumber(1,"\"Particle Number\"",       true);
    cDatumTimed    DatumNumberDensity( 1,"\"Number Density[1/m^3]\"", true);

    cDatumWeighted DatumParticleVelocity(3, "\"Vx [m/s]\", \"Vy [m/s]\", \"Vz [m/s]\"", true);
    cDatumWeighted DatumParticleVelocity2(3,"\"Vx^2 [(m/s)^2]\", \"Vy^2 [(m/s)^2]\", \"Vz^2 [(m/s)^2]\"", false);

    cDatumWeighted DatumParticleSpeed(1,            "\"|V| [m/s]\"",     true);
    cDatumWeighted DatumParticleParallelVelocity(1, "\"Vpar [m/s]\"",   false);
    cDatumWeighted DatumParticleParallelVelocity2(1,"\"Vpar^2 [(m/s)^2]\"",false);

    //-------------------------------------------------------------------------
    // IMPORTANT: some data may be oncluded only for certain species!!!!
    //if(GetSpecieType(DataSetNumber)==_PIC_SPECIE_TYPE__GAS_)
    //-------------------------------------------------------------------------
    cDatumDerived DatumTranslationalTemperature(1, "\"Translational Temperature [K]\"", true);
    cDatumDerived DatumParallelTranslationalTemperature(1, "\"Parallel Translational Temperature [K]\"", true);
    cDatumDerived DatumTangentialTranslationalTemperature(1, "\"Tangential Translational Temperature [K]\"", true);

    // vector of active sampling data
    vector<PIC::Datum::cDatumSampled*> DataSampledCenterNodeActive;
    // vector of active derived data
    vector<cDatumDerived*> DataDerivedCenterNodeActive;
  }
}

//init the block's global data
int PIC::Mesh::cDataBlockAMR::LocalTimeStepOffset=0;
int PIC::Mesh::cDataBlockAMR::LocalParticleWeightOffset=0;
int PIC::Mesh::cDataBlockAMR::totalAssociatedDataLength=0;
bool PIC::Mesh::cDataBlockAMR::InternalDataInitFlag=false;
int PIC::Mesh::cDataBlockAMR::UserAssociatedDataOffset=0;

//init the cells' global data
int PIC::Mesh::cDataCenterNode::totalAssociatedDataLength=0;
int PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset=0;

//in case OpenMP is used: tempParticleMovingListTableThreadOffset is the offset in the associatedDataPointer vector to the position when the temporary particle list begins
int PIC::Mesh::cDataBlockAMR::tempParticleMovingListTableThreadOffset=-1;
int PIC::Mesh::cDataBlockAMR::tempParticleMovingListTableThreadLength=0;

int PIC::Mesh::cDataBlockAMR::LoadBalancingMeasureOffset=0;

//the offsets to the sampled data stored in 'center nodes'
int PIC::Mesh::completedCellSampleDataPointerOffset=0,PIC::Mesh::collectingCellSampleDataPointerOffset=0;

int PIC::Mesh::sampledParticleWeghtRelativeOffset=0,PIC::Mesh::sampledParticleNumberRelativeOffset=0,PIC::Mesh::sampledParticleNumberDensityRelativeOffset=0;
int PIC::Mesh::sampledParticleVelocityRelativeOffset=0,PIC::Mesh::sampledParticleVelocity2RelativeOffset=0,PIC::Mesh::sampledParticleSpeedRelativeOffset=0;
int PIC::Mesh::sampledParticleNormalParallelVelocityRelativeOffset=0,PIC::Mesh::sampledParticleNormalParallelVelocity2RelativeOffset=0;
int PIC::Mesh::sampledExternalDataRelativeOffset=0;
int PIC::Mesh::sampleSetDataLength=0;

//domain block decomposition used in OpenMP loops
unsigned int PIC::DomainBlockDecomposition::nLocalBlocks=0;
int PIC::DomainBlockDecomposition::LastMeshModificationID=-1;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> **PIC::DomainBlockDecomposition::BlockTable=NULL;

//the mesh parameters
double PIC::Mesh::xmin[3]={0.0,0.0,0.0},PIC::Mesh::xmax[3]={0.0,0.0,0.0};
PIC::Mesh::fLocalMeshResolution PIC::Mesh::LocalMeshResolution=NULL;
#if DIM == 3
cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR> PIC::Mesh::mesh;
#elif DIM == 2
cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  PIC::Mesh::mesh;
#else
cMeshAMR1d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  PIC::Mesh::mesh;
#endif


//the user defined functions for output of the 'ceneter node' data into a data file
vector<PIC::Mesh::fPrintVariableListCenterNode> PIC::Mesh::PrintVariableListCenterNode;
vector<PIC::Mesh::fPrintDataCenterNode> PIC::Mesh::PrintDataCenterNode;
vector<PIC::Mesh::fInterpolateCenterNode> PIC::Mesh::InterpolateCenterNode;

void PIC::Mesh::AddVaraibleListFunction(fPrintVariableListCenterNode f) {
  PrintVariableListCenterNode.push_back(f);
} 



void PIC::Mesh::SetCellSamplingDataRequest() {
  exit(__LINE__,__FILE__,"not implemented yet");
}

void PIC::Mesh::cDataCenterNode::PrintVariableList(FILE* fout,int DataSetNumber) {
  // printe sampled data names
  vector<PIC::Datum::cDatumSampled*>::iterator ptrDatumSampled;

  for (ptrDatumSampled = DataSampledCenterNodeActive.begin(); ptrDatumSampled!= DataSampledCenterNodeActive.end(); ptrDatumSampled++) {
    if ((*ptrDatumSampled)->doPrint==true) (*ptrDatumSampled)->PrintName(fout);
  }

  // print derived data names
  vector<cDatumDerived*>::iterator ptrDatumDerived;

  for(ptrDatumDerived = DataDerivedCenterNodeActive.begin(); ptrDatumDerived!= DataDerivedCenterNodeActive.end(); ptrDatumDerived++) {
    if ((*ptrDatumDerived)->doPrint==true) (*ptrDatumDerived)->PrintName(fout);
  }

  //print the user defind 'center node' data
  vector<fPrintVariableListCenterNode>::iterator fptr;
  for (fptr=PrintVariableListCenterNode.begin();fptr!=PrintVariableListCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber);
  
  //if drift velocity is output -> print the variable name here
  if (_PIC_OUTPUT__DRIFT_VELOCITY__MODE_==_PIC_MODE_ON_) fprintf(fout, ", \"vxDrift\", \"vyDrift\", \"vzDrift\"");

  //print varialbes sampled by the user defined sampling procedures
  if (PIC::IndividualModelSampling::PrintVariableList.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintVariableList.size();i++) PIC::IndividualModelSampling::PrintVariableList[i](fout,DataSetNumber);
  }
}

void PIC::Mesh::cDataCenterNode::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread) {
  int idim;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  static unsigned long int nCallCounter=0;
  ++nCallCounter;
#endif

  static int  nOutput=0;
  static bool IsFirstCall=true;
  static double* OutputData;

  static vector<cDatumTimed*>    DataTimedPrint;
  static vector<cDatumWeighted*> DataWeightedPrint;
  static vector<cDatumDerived*>  DataDerivedPrint;

  // find size of message at the first call
  if (IsFirstCall==true) {
    // include sampled data
    vector<PIC::Datum::cDatumSampled*>::iterator itrDatumSampled;
    cDatumTimed*    ptrDatumTimed;
    cDatumWeighted* ptrDatumWeighted;

    for(itrDatumSampled = DataSampledCenterNodeActive.begin();itrDatumSampled!= DataSampledCenterNodeActive.end();itrDatumSampled++) {
      if ((*itrDatumSampled)->doPrint==true) {
	      nOutput+=(*itrDatumSampled)->length;

	      if((*itrDatumSampled)->type == PIC::Datum::cDatumSampled::Timed_) {
	        ptrDatumTimed = static_cast<cDatumTimed*> ((*itrDatumSampled));
	        DataTimedPrint.push_back(ptrDatumTimed);
	      }
	      else {
	        ptrDatumWeighted = static_cast<cDatumWeighted*> ((*itrDatumSampled));
	        DataWeightedPrint.push_back(ptrDatumWeighted);
	      }
      }
    }

    // include derived data
    vector<cDatumDerived*>::iterator itrDatumDerived;

    for(itrDatumDerived = DataDerivedCenterNodeActive.begin();itrDatumDerived!= DataDerivedCenterNodeActive.end(); itrDatumDerived++) if ((*itrDatumDerived)->doPrint) {
	    nOutput+=(*itrDatumDerived)->length;
	    DataDerivedPrint.push_back((*itrDatumDerived));
    }

    // allocate memory for output
    OutputData = new double[nOutput];

    // mark exit from the first (initializing) call
    IsFirstCall=false;
  }
  
  if (pipe->ThisThread==CenterNodeThread) {
    // compose a message
    int iOutput=0;
    // timed data values
    vector<cDatumTimed*>::iterator itrDatumTimed;
    vector<cDatumWeighted*>::iterator itrDatumWeighted;
    vector<cDatumDerived*>::iterator itrDatumDerived;

    for (itrDatumTimed = DataTimedPrint.begin();itrDatumTimed!= DataTimedPrint.end(); itrDatumTimed++) {
      GetDatumAverage(*(*itrDatumTimed),&OutputData[iOutput], DataSetNumber);
      iOutput += (*itrDatumTimed)->length;
    }

    for(itrDatumWeighted = DataWeightedPrint.begin();itrDatumWeighted!= DataWeightedPrint.end(); itrDatumWeighted++) {
      GetDatumAverage(*(*itrDatumWeighted),&OutputData[iOutput],DataSetNumber);
      iOutput += (*itrDatumWeighted)->length;
    }

    for(itrDatumDerived = DataDerivedPrint.begin();itrDatumDerived!= DataDerivedPrint.end(); itrDatumDerived++) {
      GetDatumAverage(*(*itrDatumDerived),&OutputData[iOutput], DataSetNumber);
      iOutput += (*itrDatumDerived)->length;
    }
  }
  
  if (pipe->ThisThread==0) {
    //print values to the output file
    if (CenterNodeThread!=0) pipe->recv((char*)OutputData,nOutput*sizeof(double),CenterNodeThread);

    for(int iOutput=0; iOutput<nOutput; iOutput++) fprintf(fout, "%e ", OutputData[iOutput]);
  }
  else pipe->send((char*)OutputData,nOutput*sizeof(double));

  //print the user defind 'center node' data
  vector<fPrintDataCenterNode>::iterator fptr;

  for (fptr=PrintDataCenterNode.begin();fptr!=PrintDataCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber,pipe,CenterNodeThread,this);

  //if drift velocity is output -> print the variable name here
  if (_PIC_OUTPUT__DRIFT_VELOCITY__MODE_==_PIC_MODE_ON_) {
    double vDrift[3];

    if (pipe->ThisThread==CenterNodeThread) {
      //calculate the drift velocity
      double BulkVelocity[3],ParticleMass,ParticleCharge;

      ParticleMass=PIC::MolecularData::GetMass(DataSetNumber);
      ParticleCharge=PIC::MolecularData::GetElectricCharge(DataSetNumber);
      GetBulkVelocity(BulkVelocity,DataSetNumber);

      //set up the points of the interpolation. move the point inside domain if on the boundary
      double xTest[3];

      for (int idim=0;idim<3;idim++) {
        xTest[idim]=x[idim];
        if (xTest[idim]==PIC::Mesh::mesh.xGlobalMax[idim]) xTest[idim]-=1.0E-10*(PIC::Mesh::mesh.xGlobalMax[idim]-PIC::Mesh::mesh.xGlobalMin[idim]);
      }

      PIC::CPLR::InitInterpolationStencil(xTest,NULL);
      PIC::CPLR::GetDriftVelocity(vDrift,BulkVelocity,ParticleMass,ParticleCharge);
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv((char*)vDrift,3*sizeof(double),CenterNodeThread);
      fprintf(fout," %e %e %e ",vDrift[0],vDrift[1],vDrift[2]);
    }
    else pipe->send((char*)vDrift,3*sizeof(double));
  }

  //print data sampled by the user defined sampling functions
  if (PIC::IndividualModelSampling::PrintSampledData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintSampledData.size();i++) PIC::IndividualModelSampling::PrintSampledData[i](fout,DataSetNumber,pipe,CenterNodeThread,this);
  }

}

void PIC::Mesh::cDataCenterNode::Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoefficients,int nInterpolationCoefficients) {
  int i,s,idim;
  double c;



  //==============================  DEBUGGER ===============
           if (nInterpolationCoefficients!=0) Measure=InterpolationList[0]->Measure;

           static long int nCallCounter=0;
           nCallCounter++;

  //============================== END DEBUGGER ============


  #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if (associatedDataPointer==NULL) exit(__LINE__,__FILE__,"Error: The associated data buffer is not initialized");
  #endif

  double InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedBulkVelocity[3]={0.0,0.0,0.0},InterpolatedBulk2Velocity[3]={0.0,0.0,0.0};
  double InterpolatedParticleSpeed=0.0;
  double pWeight;

#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
#else
  double InterpolatedBulkParallelVelocity,InterpolatedBulk2ParallelVelocity;
#endif

  for (s=0;s<PIC::nTotalSpecies;s++) {
    InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedParticleSpeed=0.0;
    for (idim=0;idim<3;idim++) InterpolatedBulkVelocity[idim]=0.0,InterpolatedBulk2Velocity[idim]=0.0;

#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
#else
    InterpolatedBulkParallelVelocity=0.0,InterpolatedBulk2ParallelVelocity=0.0;
#endif

    // interpolate all sampled data
    vector<PIC::Datum::cDatumSampled*>::iterator ptrDatum;

    for(ptrDatum = DataSampledCenterNodeActive.begin();ptrDatum!= DataSampledCenterNodeActive.end(); ptrDatum++) {
      InterpolateDatum(**ptrDatum,InterpolationList,InterpolationCoefficients,nInterpolationCoefficients, s);
    }
  }

  //print the user defind 'center node' data
  vector<fInterpolateCenterNode>::iterator fptr;

  for (fptr=InterpolateCenterNode.begin();fptr!=InterpolateCenterNode.end();fptr++) (*fptr)(InterpolationList,InterpolationCoefficients,nInterpolationCoefficients,this);

  //interpolate data sampled by user defiend sampling procedures
  if (PIC::IndividualModelSampling::InterpolateCenterNodeData.size()!=0) {
    for (unsigned int ifunc=0;ifunc<PIC::IndividualModelSampling::PrintVariableList.size();ifunc++) PIC::IndividualModelSampling::InterpolateCenterNodeData[ifunc](InterpolationList,InterpolationCoefficients,nInterpolationCoefficients,this);
  }
}

void PIC::Mesh::initCellSamplingDataBuffer() {

//  if (cDataBlockAMR::totalAssociatedDataLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of the blocks associated data offsets");

  //local time step
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  if (PIC::ThisThread==0) cout << "$PREFIX:Time step mode: specie dependent local time step" << endl;
  cDataBlockAMR::LocalTimeStepOffset=cDataBlockAMR::RequestInternalBlockData(sizeof(double)*PIC::nTotalSpecies);
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  //do nothing for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  #elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  //do nothing for _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //local particle weight
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  if (PIC::ThisThread==0) cout << "$PREFIX:Particle weight mode: specie dependent local weight" << endl;
  cDataBlockAMR::LocalParticleWeightOffset=cDataBlockAMR::RequestInternalBlockData(sizeof(double)*PIC::nTotalSpecies);
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  //do nothing for _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SINGLE_GLOBAL_PARTICLE_WEIGHT_
  //do nothing for _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //set up the offsets for 'center node' sampled data
  long int offset=0;
  DatumParticleWeight.activate(offset, &DataSampledCenterNodeActive);
  DatumNumberDensity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleNumber.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleVelocity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleVelocity2.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleSpeed.activate(offset, &DataSampledCenterNodeActive);
  DatumTranslationalTemperature.activate(&cDataCenterNode::GetTranslationalTemperature, &DataDerivedCenterNodeActive);

  //sampling the 'parallel' and 'tangential' kinetic temperatures
#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
  //do nothing
#else
  DatumParticleParallelVelocity.activate(offset, &DataSampledCenterNodeActive);
  DatumParticleParallelVelocity2.activate(offset,&DataSampledCenterNodeActive);
  DatumParallelTranslationalTemperature.activate(&cDataCenterNode::GetParallelTranslationalTemperature, &DataDerivedCenterNodeActive);
  DatumTangentialTranslationalTemperature.activate(&cDataCenterNode::GetTangentialTranslationalTemperature, &DataDerivedCenterNodeActive);
#endif

  //check if user defined sampling data is requested
  sampledExternalDataRelativeOffset=offset;

  if (PIC::IndividualModelSampling::DataSampledList.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::DataSampledList.size();i++) PIC::IndividualModelSampling::DataSampledList[i]->activate(offset, &DataSampledCenterNodeActive);
  }

  if (PIC::IndividualModelSampling::RequestSamplingData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::RequestSamplingData.size();i++) offset+=PIC::IndividualModelSampling::RequestSamplingData[i](offset);
  }




  PIC::Mesh::sampleSetDataLength=offset;

  PIC::Mesh::completedCellSampleDataPointerOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::collectingCellSampleDataPointerOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+PIC::Mesh::sampleSetDataLength;

  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=2*PIC::Mesh::sampleSetDataLength;

  //the volume partilce injection: save the volume particle injection rate
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

  //allocate the model requested static (not sampling) cell data
  if (PIC::IndividualModelSampling::RequestStaticCellData.size()!=0) {
    for (unsigned int i=0;i<PIC::IndividualModelSampling::RequestStaticCellData.size();i++) {
      PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=PIC::IndividualModelSampling::RequestStaticCellData[i](PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
    }
  }
}


//flush and switch the sampling buffers in 'center' nodes
void PIC::Mesh::flushCompletedSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::flushCollectingSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::switchSamplingBuffers() {
  int tempOffset;

  //exchange the CellSampleData offsets
  tempOffset=PIC::Mesh::completedCellSampleDataPointerOffset;
  PIC::Mesh::completedCellSampleDataPointerOffset=PIC::Mesh::collectingCellSampleDataPointerOffset;
  PIC::Mesh::collectingCellSampleDataPointerOffset=tempOffset;

  //switch the offsets for the internal spherical surfaces installed into the mesh
  PIC::BC::InternalBoundary::Sphere::switchSamplingBuffers();

}

//==============================================================
//init set and set up the computational mesh
void PIC::Mesh::Init(double* xMin,double* xMax,fLocalMeshResolution ResolutionFunction) {

  for (int idim=0;idim<DIM;idim++) xmin[idim]=xMin[idim],xmax[idim]=xMax[idim];

  LocalMeshResolution=ResolutionFunction;
  mesh.init(xMin,xMax,LocalMeshResolution);
}

void PIC::Mesh::buildMesh() {
  mesh.buildMesh();
}



void PIC::Mesh::cDataBlockAMR::sendBoundaryLayerBlockData(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->send(cell->associatedDataPointer,cell->totalAssociatedDataLength);
    pipe->send(cell->Measure);
  }

  pipe->send(associatedDataPointer+UserAssociatedDataOffset,totalAssociatedDataLength-UserAssociatedDataOffset);
}

void PIC::Mesh::cDataBlockAMR::sendMoveBlockAnotherProcessor(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
//  PIC::Mesh::cDataCenterNode *cell=NULL;

  sendBoundaryLayerBlockData(pipe);

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  //send all blocks' data when the blocks is moved to another processor
  long int Particle,NextParticle;
  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
//    LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
//    cell=GetCenterNode(LocalCellNumber);

    //    Particle=cell->FirstCellParticle;
    Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

    if  (Particle!=-1) {
      LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      pipe->send(_CENTRAL_NODE_NUMBER_SIGNAL_);
      pipe->send(LocalCellNumber);

      while (Particle!=-1) {
        PIC::ParticleBuffer::PackParticleData(buffer,Particle);
        pipe->send(_NEW_PARTICLE_SIGNAL_);
        pipe->send(buffer,PIC::ParticleBuffer::ParticleDataLength);

        NextParticle=PIC::ParticleBuffer::GetNext(Particle);
//        PIC::ParticleBuffer::DeleteParticle(Particle);
        PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(Particle);

        Particle=NextParticle;
      }

     // cell->FirstCellParticle=-1;
      FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=-1;
    }
  }

  pipe->send(_END_COMMUNICATION_SIGNAL_);
  delete [] buffer;
}

void PIC::Mesh::cDataBlockAMR::recvBoundaryLayerBlockData(CMPI_channel *pipe,int From) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->recv(cell->associatedDataPointer,cell->totalAssociatedDataLength,From);
    pipe->recv(cell->Measure,From);
  }

  pipe->recv(associatedDataPointer+UserAssociatedDataOffset,totalAssociatedDataLength-UserAssociatedDataOffset,From);
}

//recieve all blocks' data when the blocks is moved to another processo
void PIC::Mesh::cDataBlockAMR::recvMoveBlockAnotherProcessor(CMPI_channel *pipe,int From) {
  long int LocalCellNumber=-1;
//  PIC::Mesh::cDataCenterNode *cell=NULL;
  int i=-10,j=-10,k=-10;

  recvBoundaryLayerBlockData(pipe,From);

  long int Particle;
  char buffer[PIC::ParticleBuffer::ParticleDataLength];

//  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  int Signal;
  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;


  pipe->recv(Signal,From);
  //LocalCellNumber=-1,cell=NULL;

  while (Signal!=_END_COMMUNICATION_SIGNAL_) {
    switch (Signal) {
    case _CENTRAL_NODE_NUMBER_SIGNAL_ :
      pipe->recv(LocalCellNumber,From);
//      cell=GetCenterNode(LocalCellNumber);

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(LocalCellNumber,i,j,k);
      break;
    case _NEW_PARTICLE_SIGNAL_ :
      pipe->recv(buffer,PIC::ParticleBuffer::ParticleDataLength,From);

//      Particle=PIC::ParticleBuffer::GetNewParticle(cell->FirstCellParticle);

      Particle=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)],true);
      PIC::ParticleBuffer::UnPackParticleData(buffer,Particle);
      break;
    default :
      exit(__LINE__,__FILE__,"Error: unknown option");
    }

    pipe->recv(Signal,From);
  }

//  delete [] buffer;
}

//===============================================================================================================
//update the domain decomposition table
void PIC::DomainBlockDecomposition::UpdateBlockTable() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  if (LastMeshModificationID==PIC::Mesh::mesh.nMeshModificationCounter) return; //no modification of the mesh is made

  //calculate the new number of the blocks
  for (nLocalBlocks=0,node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    nLocalBlocks++;
  }

  //deallocate and allocat the block pointe buffer
  if (BlockTable!=NULL) delete [] BlockTable;
  BlockTable=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nLocalBlocks];

  //populate the block pointer buffer
  for (nLocalBlocks=0,node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    BlockTable[nLocalBlocks++]=node;
  }

  //update the domain decomposition table ID
  LastMeshModificationID=PIC::Mesh::mesh.nMeshModificationCounter;
}


//===============================================================================================================
//get the interpolation stencil for visualization of the model results (used only when the linear interpolation routine is set)
int PIC::Mesh::GetCenterNodesInterpolationCoefficients(double *x,double *CoefficientsList,PIC::Mesh::cDataCenterNode **InterpolationStencil,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,int nMaxCoefficients) {
  int iCell,cnt=0;
  double SumWeight=0.0;

  if (_PIC_COUPLER__INTERPOLATION_MODE_ != _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    exit(__LINE__,__FILE__,"Error: the function should be used only when the linear interpolation routine is set");
  }

  //construct the interpolation stencil
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  PIC::CPLR::InitInterpolationStencil(x,startNode);

  //if the length of the coefficient list is not enough -> exist with an error message
  if (nMaxCoefficients<PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Length) {
    exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
    return -1;
  }

  for (iCell=0;iCell<PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Length;iCell++) {
    if (PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].cell[iCell]->Measure>0.0) {
      CoefficientsList[cnt]=PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Weight[iCell];
      InterpolationStencil[cnt]=PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].cell[iCell];

      SumWeight+=CoefficientsList[cnt];
      cnt++;
    }
  }

  if (cnt!=0) for (int ii=0;ii<cnt;ii++) CoefficientsList[ii]/=SumWeight;

  return cnt;
}



//===============================================================================================================
//cell dat aaccess routines
void PIC::Mesh::cDataCenterNode::SampleDatum(Datum::cDatumSampled Datum, double* In, int spec, double weight) {
  for(int i=0; i<Datum.length; i++) {
    *(i + Datum.length * spec + (double*)(associatedDataPointer + collectingCellSampleDataPointerOffset+Datum.offset))+= In[i] * weight;
  }
}

void PIC::Mesh::cDataCenterNode::SampleDatum(Datum::cDatumSampled Datum, double In, int spec,  double weight) {
  *(spec + (double*)(associatedDataPointer + collectingCellSampleDataPointerOffset+Datum.offset))+= In * weight;
}

//.......................................................................
void PIC::Mesh::cDataCenterNode::SetDatum(Datum::cDatumSampled Datum, double* In, int spec) {
  for(int i=0; i<Datum.length; i++) *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) = In[i];
}

//get accumulated data
//.......................................................................
void PIC::Mesh::cDataCenterNode::GetDatumCumulative(Datum::cDatumSampled Datum, double* Out, int spec) {
  for(int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset));
}

double PIC::Mesh::cDataCenterNode::GetDatumCumulative(Datum::cDatumSampled Datum, int spec) {
  return *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset));
}

//get data averaged over time
void PIC::Mesh::cDataCenterNode::GetDatumAverage(cDatumTimed Datum, double* Out, int spec) {
  if (PIC::LastSampleLength > 0) for (int i=0; i<Datum.length; i++) {
    Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / PIC::LastSampleLength;
  }
  else for (int i=0; i<Datum.length; i++) Out[i] = 0.0;
}

double PIC::Mesh::cDataCenterNode::GetDatumAverage(cDatumTimed Datum, int spec) {
  return (PIC::LastSampleLength > 0) ?  *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / PIC::LastSampleLength : 0.0;
}

//get data averaged over sampled weight
void PIC::Mesh::cDataCenterNode::GetDatumAverage(cDatumWeighted Datum, double* Out, int spec) {
  double TotalWeight=0.0;

  GetDatumCumulative(DatumParticleWeight, &TotalWeight, spec);

  if (TotalWeight > 0) {
    for(int i=0; i<Datum.length; i++) Out[i] = *(i + Datum.length * spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / TotalWeight;
  }
  else for(int i=0; i<Datum.length; i++) Out[i] = 0.0;
}

double PIC::Mesh::cDataCenterNode::GetDatumAverage(cDatumWeighted Datum, int spec) {
  double TotalWeight=0.0;

  GetDatumCumulative(DatumParticleWeight, &TotalWeight, spec);
  return (TotalWeight > 0) ? *(spec + (double*)(associatedDataPointer + completedCellSampleDataPointerOffset+Datum.offset)) / TotalWeight : 0.0;
}

//get average for derived data
void PIC::Mesh::cDataCenterNode::GetDatumAverage(cDatumDerived Datum, double* Out, int spec) {
  (this->*Datum.GetAverage)(Out,spec);
}
//-----------------------------------------------------------------------

//backward compatible access
//-----------------------------------------------------------------------
double PIC::Mesh::cDataCenterNode::GetNumberDensity(int spec) {
  return GetDatumAverage(DatumNumberDensity, spec);
}

void PIC::Mesh::cDataCenterNode::GetBulkVelocity(double* vOut, int spec) {
  GetDatumAverage(DatumParticleVelocity, vOut, spec);
}

double PIC::Mesh::cDataCenterNode::GetMeanParticleSpeed(int spec) {
  return GetDatumAverage(DatumParticleSpeed, spec);
}

double PIC::Mesh::cDataCenterNode::GetCompleteSampleCellParticleWeight(int spec) {
  return GetDatumAverage(DatumParticleWeight, spec);
}
//-----------------------------------------------------------------------

// data interpolation
//-----------------------------------------------------------------------
void PIC::Mesh::cDataCenterNode::InterpolateDatum(Datum::cDatumSampled Datum, cDataCenterNode** InterpolationList,double *InterpolationCoefficients,int nInterpolationCoefficients, int spec) {
  // container for the interpolated value; set it to be zero
  double value[Datum.length], interpolated[Datum.length];

  for (int i=0; i<Datum.length; i++) {
    value[i]=0.0; interpolated[i]=0.0;
  }

  // interpolation loop
  for (int i=0; i<nInterpolationCoefficients; i++) {
    InterpolationList[i]->GetDatumCumulative(Datum, value, spec);

    for(int j=0; j<Datum.length; j++) interpolated[j] += InterpolationCoefficients[i] * value[j];
  }

  SetDatum(Datum, interpolated, spec);

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    for(int i=0; i<Datum.length; i++) if (isfinite(interpolated[i])==false) exit(__LINE__,__FILE__,"Error: Floating Point Exception");
  #endif
  #endif
}

