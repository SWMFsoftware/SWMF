//$Id$
//sampling function for the Earth model


/*
 * Earth_Sampling.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "Earth.h"

int Earth::Sampling::ParticleData::_GYRO_RADIUS_SAMPLING_OFFSET_=-1;
int Earth::Sampling::ParticleData::_GYRO_FRECUENCY_SAMPLING_OFFSET_=-1;
bool Earth::Sampling::ParticleData::SamplingMode=false;


//init the sampling module
void Earth::Sampling::ParticleData::Init() {
  if (SamplingMode==true) {
    //request sampling data
    PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);

    //print out of the otuput file
    PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
    PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
    PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
  }
}

//request sampling data
int Earth::Sampling::ParticleData::RequestSamplingData(int offset) {
  int DataLength=0;

  _GYRO_RADIUS_SAMPLING_OFFSET_=offset;
  DataLength+=PIC::nTotalSpecies*sizeof(double);

  _GYRO_FRECUENCY_SAMPLING_OFFSET_=offset+DataLength;
  DataLength+=PIC::nTotalSpecies*sizeof(double);

  return DataLength;
}

//sample particle data
void Earth::Sampling::ParticleData::SampleParticleData(char* tempParticleData,double LocalParticleWeight,char *SamplingData,int s,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double B[3],x[3],v[3];
  double GyroRadius=0.0,GyroFrequency=0.0;
  double mass,ElectricCharge;

  ElectricCharge=fabs(PIC::MolecularData::GetElectricCharge(s));

  if (ElectricCharge>0.0) {
    PIC::ParticleBuffer::GetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::GetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);

    PIC::CPLR::InitInterpolationStencil(x,node);
    PIC::CPLR::GetBackgroundMagneticField(B);

    mass=PIC::MolecularData::GetMass(s);

    GyroRadius=Relativistic::GetGyroRadius(v,mass,ElectricCharge,B);
    GyroFrequency=Relativistic::GetGyroFrequency(v,mass,ElectricCharge,B);
  }

  *(s+(double*)(SamplingData+_GYRO_RADIUS_SAMPLING_OFFSET_))+=GyroRadius*LocalParticleWeight;
  *(s+(double*)(SamplingData+_GYRO_FRECUENCY_SAMPLING_OFFSET_))+=GyroFrequency*LocalParticleWeight;
}

//interpolate the cell aberaged data on a corner of a cell
void Earth::Sampling::ParticleData::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int spec,i;
  double GyroRadiiTable[PIC::nTotalSpecies];
  double GyroFrequencyTable[PIC::nTotalSpecies];
  double TotalMeasure=0.0,Measure=0.0;
  char *SamplingBuffer,*CellNodeSamplingBuffer;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) GyroRadiiTable[spec]=0.0,GyroFrequencyTable[spec]=0.0;

  CellNodeSamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

    Measure=InterpolationCoeficients[i];
    TotalMeasure+=Measure;

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      GyroRadiiTable[spec]+=(*(spec+(double*)(SamplingBuffer+_GYRO_RADIUS_SAMPLING_OFFSET_)))*Measure;
      GyroFrequencyTable[spec]+=(*(spec+(double*)(SamplingBuffer+_GYRO_FRECUENCY_SAMPLING_OFFSET_)))*Measure;
    }
  }


  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    *(spec+(double*)(CellNodeSamplingBuffer+_GYRO_RADIUS_SAMPLING_OFFSET_))=GyroRadiiTable[spec]/TotalMeasure;
    *(spec+(double*)(CellNodeSamplingBuffer+_GYRO_FRECUENCY_SAMPLING_OFFSET_))=GyroFrequencyTable[spec]/TotalMeasure;
  }
}


//output the variable list
void Earth::Sampling::ParticleData::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Mean Gyro Radius\", \"Mean Gyro Frequency\"");
}

//print sampled data into a file
void Earth::Sampling::ParticleData::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  char *SamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
  double TotalWeight,GyroRadius,GyroFrequency;


  if (pipe->ThisThread==CenterNodeThread) {
    TotalWeight=CenterNode->GetDatumCumulative(PIC::Mesh::DatumParticleWeight,DataSetNumber);

    if (TotalWeight>0.0) {
      GyroRadius=(*(DataSetNumber+(double*)(SamplingBuffer+_GYRO_RADIUS_SAMPLING_OFFSET_)))/TotalWeight;
      GyroFrequency=(*(DataSetNumber+(double*)(SamplingBuffer+_GYRO_FRECUENCY_SAMPLING_OFFSET_)))/TotalWeight;
    }
    else GyroRadius=0.0,GyroFrequency=0.0;
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      GyroRadius=pipe->recv<double>(CenterNodeThread);
      GyroFrequency=pipe->recv<double>(CenterNodeThread);
    }

    fprintf(fout," %e  %e ",GyroRadius,GyroFrequency);
  }
  else {
    pipe->send(GyroRadius);
    pipe->send(GyroFrequency);
  }
}
















