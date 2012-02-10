//====================================================
//$Id$
//====================================================
//reads the input file for the pic solver 

#include "pic.h"
//====================================================
void PIC::Parser::Run(char* InputFile) {
  CiFileOperations ifile;
//  char nonstandardBlock_endname[_MAX_STRING_LENGTH_PIC_];
//  bool nonstandardBlock_flag=false;
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::ThisThread==0) printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
	printf("Cannot find the input file:%s\n",InputFile);
	exit(__LINE__,__FILE__);
  }

  ifile.openfile(InputFile);

  while (ifile.eof()==false) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(str1,str);

    //check if the block is used defined
    bool usedDefinedBlock=false;

/*
     for (list<TExternalInputFileReader>::iterator ptr=ExternalInputFileReaderList.begin();ptr!=ExternalInputFileReaderList.end();ptr++) if ((*ptr)(str1,ifile)==true) {
	      usedDefinedBlock=true;
	      break;

	    }
*/

    if (usedDefinedBlock==true) {
      if (ThisThread==0) printf("Read a user defined block \"%s\"\n",str1);
	}
	else if (strcmp("#MAIN",str1)==0) readMain(ifile);
	else if (strcmp("#SPECIES",str1)==0) PIC::MolecularData::Parser::run(ifile);
	else if (strcmp("#END",str1)==0) return;
	else if (strcmp("",str1)!=0) ifile.error();
  }
}

//====================================================
//read the main block of the input file
void PIC::Parser::readMain(CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_];
  char *endptr;

  while (ifile.eof()==false) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(str1,str);

    if (strcmp("MOLECULARMODEL",str1)==0) {
	  ifile.CutInputStr(str1,str);
	  if (strcmp("HS",str1)==0) PIC::MolecularData::SetMolType(_HS_MOL_MODEL_);
	  else if (strcmp("VHS",str1)==0) PIC::MolecularData::SetMolType(_VHS_MOL_MODEL_);
	  else if (strcmp("VSS",str1)==0) PIC::MolecularData::SetMolType(_VSS_MOL_MODEL_);
	  else ifile.error();}
	else if (strcmp("EXTERNALSPECIES",str1)==0) {
	  ifile.CutInputStr(str1,str);
	  if (strcmp("ON",str1)==0) PIC::MolecularData::ExternalSpeciesModelingFlag=_EXTERNAL_SPECIES_ON_;
	  else if (strcmp("OFF",str1)==0) PIC::MolecularData::ExternalSpeciesModelingFlag=_EXTERNAL_SPECIES_OFF_;
	  else ifile.error();}
	else if (strcmp("UNIMOLECULARTRANSFORMATIONS",str1)==0) {
	  ifile.CutInputStr(str1,str);
	  if (strcmp("ON",str1)==0) PIC::MolecularData::UnimolecularReactionFlag=_UNIMOLECULAR_REACTIONS_ON_;
	  else if (strcmp("OFF",str1)==0) PIC::MolecularData::UnimolecularReactionFlag=_UNIMOLECULAR_REACTIONS_OFF_;
	  else ifile.error();}
	else if (strcmp("INTERNALDEGREESOFFREEDOM",str1)==0) {
	  ifile.CutInputStr(str1,str);
	  if (strcmp("ON",str1)==0) PIC::MolecularData::InternalDegreesOfFreedomModelingFlag=_INTERNAL_DEGRESS_OF_FREEDOM_OFF_;
	  else if (strcmp("OFF",str1)==0) PIC::MolecularData::InternalDegreesOfFreedomModelingFlag=_INTERNAL_DEGRESS_OF_FREEDOM_OFF_;
	  else ifile.error();}
    else if (strcmp("NS",str1)==0) {
	  ifile.CutInputStr(str1,str);
	  PIC::nTotalSpecies=(unsigned char)strtol(str1,&endptr,10);
	  if (PIC::nTotalSpecies<=0) exit(__LINE__,__FILE__,"NS is out of the range");
	  if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();}
	else if (strcmp("SPECIESLIST",str1)==0) {
	  int i,spec;

	  //initialize the table of chemical species
	  if (PIC::MolecularData::ChemTable!=NULL) exit(__LINE__,__FILE__"The chemical table is already defined");
	  if (PIC::nTotalSpecies==0) exit(__LINE__,__FILE__,"The value of NS is nor defined yet");

	  PIC::MolecularData::ChemTable=new char* [PIC::nTotalSpecies];
	  PIC::MolecularData::ChemTable[0]=new char[PIC::nTotalSpecies*_MAX_STRING_LENGTH_PIC_];
	  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
	    PIC::MolecularData::ChemTable[spec]=PIC::MolecularData::ChemTable[0]+spec*_MAX_STRING_LENGTH_PIC_;
	    PIC::MolecularData::ChemTable[spec][0]='\0';
	  }

    PIC::MolecularData::LoadingSpeciesList=new char* [PIC::nTotalSpecies];
    PIC::MolecularData::LoadingSpeciesList[0]=new char[PIC::nTotalSpecies*_MAX_STRING_LENGTH_PIC_];
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      PIC::MolecularData::LoadingSpeciesList[spec]=PIC::MolecularData::LoadingSpeciesList[0]+spec*_MAX_STRING_LENGTH_PIC_;
      PIC::MolecularData::LoadingSpeciesList[spec][0]='\0';
    }

    PIC::MolecularData::SpcecieTypeTable=new int [PIC::nTotalSpecies];
    for (spec=0;spec<PIC::nTotalSpecies;spec++) PIC::MolecularData::SpcecieTypeTable[spec]=-1;

	  spec=0;

	   while (str[0]!='\0') {
	     ifile.CutInputStr(str1,str);
	     for (i=0;str1[i]!='\0';i++) PIC::MolecularData::LoadingSpeciesList[spec][i]=str1[i];
	     PIC::MolecularData::LoadingSpeciesList[spec][i]='\0';

	     spec+=1;
 	   }
	 }
	 else if (strcmp("#ENDMAIN",str1)==0) {
//     mol.init(NS);
       return;}
     else ifile.error();
   }
}

