//$Id$
//functions and variables that are user to process errors generated during the particle moving procedure

/*
 * pic_mover_error.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include <sstream>
#include <iostream>
#include <string>

std::vector<int> PIC::Mover::Sampling::Errors::RemovedModelParticles;
std::vector<double> PIC::Mover::Sampling::Errors::ModelParticleRemovingRate;
std::vector<std::string> PIC::Mover::Sampling::Errors::ErrorLineID;

int PIC::Mover::Sampling::Errors::ErrorDetectionFlag=0;

void PIC::Mover::Sampling::Errors::Init() {
  int spec;

  ErrorDetectionFlag=0;

  ErrorLineID.clear();
  RemovedModelParticles.clear();
  ModelParticleRemovingRate.clear();
}

void PIC::Mover::Sampling::Errors::PrintData() {

  //combine the error data from all processors
  int ErrorDetected=0;
  MPI_Allreduce(&ErrorDetectionFlag,&ErrorDetected,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  //combine and output the error data
  if (ErrorDetected!=0) {

    if (PIC::ThisThread==0) {
      int spec,thread,nErrorIDs,nerror;
      int RecvBufferNumber[PIC::nTotalSpecies];
      double RecvBufferRate[PIC::nTotalSpecies];
      MPI_Status status;
      char ErrorID[_MAX_STRING_LENGTH_PIC_];
      std::string ErrorIDstring;

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        MPI_Recv(&nErrorIDs,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (nerror=0;nerror<nErrorIDs;nerror++) {
          MPI_Recv(ErrorID,_MAX_STRING_LENGTH_PIC_,MPI_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          ErrorIDstring=ErrorID;

          MPI_Recv(RecvBufferNumber,PIC::nTotalSpecies,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          MPI_Recv(RecvBufferRate,PIC::nTotalSpecies,MPI_DOUBLE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

          for (spec=0;spec<PIC::nTotalSpecies;spec++) AddRemovedParticleData(RecvBufferRate[spec],RecvBufferNumber[spec],spec,ErrorIDstring);
        }
      }

    }
    else {
      int nErrorIDs=ErrorLineID.size();
      int spec,nerror;

      MPI_Send(&nErrorIDs,1,MPI_INT,0,0,MPI_GLOBAL_COMMUNICATOR);

      for (nerror=0;nerror<nErrorIDs;nerror++) {
        char ErrorID[_MAX_STRING_LENGTH_PIC_];

        strcpy(ErrorID, ErrorLineID[nerror].c_str());
        MPI_Send(ErrorID,_MAX_STRING_LENGTH_PIC_,MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);

        int SendBufferNumber[PIC::nTotalSpecies];
        double SendBufferRate[PIC::nTotalSpecies];

        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          SendBufferNumber[spec]=RemovedModelParticles[spec+nerror*PIC::nTotalSpecies];
          SendBufferRate[spec]=ModelParticleRemovingRate[spec+nerror*PIC::nTotalSpecies];
        }

        MPI_Send(SendBufferNumber,PIC::nTotalSpecies,MPI_INT,0,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Send(SendBufferRate,PIC::nTotalSpecies,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);
      }
    }



    //print the sampled data
    if (PIC::ThisThread==0) {
      printf("$PREFIX: MODEL PARTICLES ARE REMOVED BY THE PARTICLE MOVING PROCEDURE DUE TO AN ERROR\n");
      printf("$PREFIX: SPEC\tNUMBER OF THE REMOVED PARTICLE PER ITERATION\tREMOVING RATE PER SEC\n");

      for (int nerror=0;nerror<ErrorLineID.size();nerror++) {
        printf("$PREFIX: ERROR ID: %s\n",ErrorLineID[nerror].c_str());

        for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
          printf("$PREFIX: %i\t%i\t%e\n",spec,
              RemovedModelParticles[spec+nerror*PIC::nTotalSpecies]/PIC::LastSampleLength,
              ModelParticleRemovingRate[spec+nerror*PIC::nTotalSpecies]/PIC::LastSampleLength);
        }
      }
    }

  }


  //reset the error buffers
  Init();
}

void PIC::Mover::Sampling::Errors::AddRemovedParticleData(double Rate, int spec, int line,const char *fname) {
  std::string FullLineID;

  //create the full erro line ID and comapre with the previously registered line IDs
  std::ostringstream f1,f2;

  f1 << line;
  f2 << fname;

  FullLineID=f1.str() + "@" + f2.str();

  AddRemovedParticleData(Rate,spec,FullLineID);
}

void PIC::Mover::Sampling::Errors::AddRemovedParticleData(double Rate, int spec, std::string &FullLineID) {
  AddRemovedParticleData(Rate,1,spec,FullLineID);
}

void PIC::Mover::Sampling::Errors::AddRemovedParticleData(double Rate, int nRemovedParticles, int spec, std::string &FullLineID) {
  ErrorDetectionFlag=1;

  //search for the previously registered error
  bool foundflag=false;
  int ErrorID=0;

  for (std::vector<std::string>::iterator ptr=ErrorLineID.begin();ptr!=ErrorLineID.end();ptr++) {
    if (ptr->compare(FullLineID)==0) {
      foundflag=true;
      break;
    }

    ErrorID++;
  }

  if (foundflag==false) {
    //the error position is not found -> append the error list
    vector<double> Rate(PIC::nTotalSpecies);
    vector<int> Number(PIC::nTotalSpecies);

    for (int i=0;i<PIC::nTotalSpecies;i++) Rate[i]=0.0,Number[i]=0;

    RemovedModelParticles.insert(RemovedModelParticles.end(),Number.begin(),Number.end());
    ModelParticleRemovingRate.insert(ModelParticleRemovingRate.end(),Rate.begin(),Rate.end());
    ErrorLineID.push_back(FullLineID);
  }

  //save the error data
  RemovedModelParticles[spec+ErrorID*PIC::nTotalSpecies]+=nRemovedParticles;
  ModelParticleRemovingRate[spec+ErrorID*PIC::nTotalSpecies]+=Rate;
}
