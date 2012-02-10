//==========================================================
//$Id$
//==========================================================
//molecular data of model particles 

#include "pic.h"

//init the global variables
double* PIC::MolecularData::MolMass=NULL;
double* PIC::MolecularData::ElectricCharge=NULL;
char** PIC::MolecularData::ChemTable=NULL;
int *PIC::MolecularData::SpcecieTypeTable=NULL;
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
//set and get molecular mass of the species
double PIC::MolecularData::GetMass(int s) {return MolMass[s];}

void PIC::MolecularData::SetMass(double mass,int spec) {
  if (MolMass==NULL) {
    MolMass=new double[PIC::nTotalSpecies];
    for (int i=0;i<PIC::nTotalSpecies;i++) MolMass[i]=0.0;
  }

  MolMass[spec]=mass;
}

/*
extern double *MolMass;
void SetMass(double,int);
double GetMass(int);
*/

//==========================================================
//get the counting number of the specie
int PIC::MolecularData::GetSpecieNumber(char* sname) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"The chemical table is not defined");
  for (int s=0;s<PIC::nTotalSpecies;s++) if (strcmp(ChemTable[s],sname)==0) return s;

  return -1;
}

//==========================================================
//set and get the specie symbol
void PIC::MolecularData::SetChemSymbol(char* sym,int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  sprintf(ChemTable[spec],"%s",sym);
}


void PIC::MolecularData::GetChemSymbol(char* sym,int spec) {
  if (ChemTable==NULL) exit(__LINE__,__FILE__,"ChemTable is not initialized");

  sprintf(sym,"%s",ChemTable[spec]);
}

//===========================================================
//set and get the species type
void PIC::MolecularData::SetSpecieType(int SpcecieType,int spec) {
  if (SpcecieTypeTable==NULL) exit(__LINE__,__FILE__,"Error: SpcecieTypeTable is not initialized");

  SpcecieTypeTable[spec]=SpcecieType;
}

int PIC::MolecularData::GetSpecieType(int spec) {
  if (SpcecieTypeTable==NULL) exit(__LINE__,__FILE__,"Error: SpcecieTypeTable is not initialized");

  return SpcecieTypeTable[spec];
}
























