//====================================================
//$Id$
//====================================================
//reads the file with the molecules data

#include "pic.h"



//====================================================
//read the part of the input file related to a particular specie
void PIC::MolecularData::Parser::SpeciesBlock(int SpecieNumber,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_];
  char *endptr;
  double a;

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
  int s;
  bool SpecieFound;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);
    if (strcmp("#COMPONENT",str1)==0) {
      ifile.CutInputStr(str1,str);
      SpecieFound=false;

      for (s=0;s<PIC::nTotalSpecies;s++) if (strcmp(ChemTable[s],str1)==0) {
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
      else if ((strcmp("EXTERNALSPEC",str1)!=0)&&(strcmp("BACKGROUNDSPEC",str1)!=0)) SpeciesBlock(PIC::MolecularData::GetSpecieNumber(str1),ifile);
    }
    else if (strcmp("#ENDSPECIES",str1)==0) {
//      PIC::MolecularData::Init();
      return;
    }
    else ifile.error();
  }


}
