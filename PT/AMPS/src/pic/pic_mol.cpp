//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//==========================================================
//$Id$
//==========================================================
//molecular data of model particles 

#include "pic.h"

//init the global variables
//double* PIC::MolecularData::MolMass=NULL;
//double* PIC::MolecularData::ElectricCharge=NULL;
//char** PIC::MolecularData::ChemTable=NULL;
//int *PIC::MolecularData::SpcecieTypeTable=NULL;
char** PIC::MolecularData::LoadingSpeciesList=NULL;
int PIC::MolecularData::MolModelCode=_HS_MOL_MODEL_;
bool PIC::MolecularData::ExternalSpeciesModelingFlag=_EXTERNAL_SPECIES_OFF_;
bool PIC::MolecularData::UnimolecularReactionFlag=_UNIMOLECULAR_REACTIONS_OFF_;
bool PIC::MolecularData::InternalDegreesOfFreedomModelingFlag=_INTERNAL_DEGRESS_OF_FREEDOM_OFF_;

//the functions:
//==========================================================
//set and get current molecular model
void PIC::MolecularData::SetMolType(int Code) {MolModelCode=Code;}
int PIC::MolecularData::GetMolType() {return MolModelCode;}

//==========================================================
//get the counting number of the specie
int PIC::MolecularData::GetSpecieNumber(char* sname) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"The chemical table is not defined");
  for (int s=0;s<PIC::nTotalSpecies;s++) if (strcmp(ChemTable[s],sname)==0) return s;

  return -1;
}

//==========================================================
//get the specie BASE symbol
void PIC::MolecularData::GetChemBaseSymbol(char* sym,int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  sprintf(sym,"%s",GetChemBaseSymbol(spec));
}

char* PIC::MolecularData::GetChemBaseSymbol(int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  static char TempString[_MAX_STRING_LENGTH_PIC_];

  for (int i=0;i<_MAX_STRING_LENGTH_PIC_;i++) {
    TempString[i]=ChemTable[spec][i];

    if ((TempString[i]=='%')||(TempString[i]=='\n')||(TempString[i]=='\r')||(TempString[i]==0)) {
      TempString[i]=0;
      break;
    }
  }

  return TempString;
}

//==========================================================
//get the specie symbol
void PIC::MolecularData::GetChemSymbol(char* sym,int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  sprintf(sym,"%s",ChemTable[spec]);
}

const char* PIC::MolecularData::GetChemSymbol(int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  return ChemTable[spec];
}






















