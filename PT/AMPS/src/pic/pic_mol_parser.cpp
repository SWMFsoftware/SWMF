//====================================================
//$Id$
//====================================================
//reads the file with the molecules data

#include "pic.h"



//====================================================
//read the part of the input file related to a particular specie
void PIC::MolecularData::Parser::SpeciesBlock(char *ChemSymbol,int SpecieNumber,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_];
  char *endptr;
  double a;

  PIC::MolecularData::SetChemSymbol(ChemSymbol,SpecieNumber);
  PIC::MolecularData::SetSpecieType(_PIC_SPECIE_TYPE__GAS_,SpecieNumber);

  while (ifile.eof()==false) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(str1,str);

    if ((strcmp("#COMPONENT",str1)==0)||(strcmp("#ENDSPECIES",str1)==0)) {
      ifile.moveLineBack();
      return;
    }
    else if (strcmp("MASS",str1)==0) {
      ifile.CutInputStr(str1,str);
      a=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();
      PIC::MolecularData::SetMass(a,SpecieNumber);}

    else ifile.error();
  }
}
//====================================================
//control the parsing of the molecule data's part of the input file
void PIC::MolecularData::Parser::run(CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_];

  //read the data file
  int s,SpecieType;
  bool SpecieFound;

  int LoadedSpecieCountingNumber=0;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);
    if (strcmp("#COMPONENT",str1)==0) {
      ifile.CutInputStr(str1,str);
      SpecieFound=false;
      SpecieType=_PIC_SPECIE_TYPE__GAS_;

      if (strcmp("EXTERNALSPEC",str1)==0) {
        SpecieType=_PIC_SPECIE_TYPE__EXTERNAL_;
        ifile.CutInputStr(str1,str);
      }
      else if (strcmp("BACKGROUNDSPEC",str1)==0) {
        SpecieType=_PIC_SPECIE_TYPE__BACKGROUND_;
        ifile.CutInputStr(str1,str);
      }

      for (s=0;s<PIC::nTotalSpecies;s++) if (strcmp(LoadingSpeciesList[s],str1)==0) {
        SpecieFound=true;
        break;
      }

      if (SpecieFound==false) {
        do {
          ifile.GetInputStr(str,sizeof(str));
          ifile.CutInputStr(str1,str);
        }
        while ((strcmp("#COMPONENT",str1)!=0)&&(strcmp("#ENDSPECIES",str1)!=0));

        ifile.moveLineBack();
      }
      else if (SpecieType==_PIC_SPECIE_TYPE__GAS_) {
        SpeciesBlock(str1,LoadedSpecieCountingNumber,ifile);
        LoadedSpecieCountingNumber++;
      }
      else if ((SpecieType==_PIC_SPECIE_TYPE__EXTERNAL_)||(SpecieType==_PIC_SPECIE_TYPE__BACKGROUND_)) {
        PIC::MolecularData::SetChemSymbol(str1,LoadedSpecieCountingNumber);
        PIC::MolecularData::SetSpecieType(SpecieType,LoadedSpecieCountingNumber);
        LoadedSpecieCountingNumber++;
      }
    }
    else if (strcmp("#ENDSPECIES",str1)==0) {
//      PIC::MolecularData::Init();
      return;
    }
    else ifile.error();
  }

  if (LoadedSpecieCountingNumber!=PIC::nTotalSpecies) exit(__LINE__,__FILE__,"Error: the number of loaded species is not consistent with the requested number of species");

}
