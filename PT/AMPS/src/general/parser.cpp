//===================================================
//$Id$
//===================================================

#ifdef MPI_ON
#include "mpi.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "data.h"
#include "parser.h"

#if CompilationTarget==DSMCTARGET
  #include "dsmc.h"
  #include "mol.h"
#elif CompilationTarget==PICTARGET
  #include "dsmc.h"
  #include "mol.h"
  #include "pic.h"
#elif CompilationTarget==EULERTARGET
  #include "euler.h"
#elif CompilationTarget==HYBRIDTARGET
  #include "dsmc.h"
  #include "mol.h"
  #include "euler.h"
  #include "hybrid.h"
#endif

list <PARSER::TExternalInputFileReader> PARSER::ExternalInputFileReaderList;

//===================================================
void PARSER::GeneralBlock(CiFileOperations& ifile) {
  static bool GeneralBlockLoadedFlag=false;
  char str1[300],str[300];
  char *endptr;

  if (GeneralBlockLoadedFlag!=false) {
    do {
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(str1,str);
    }
    while((ifile.eof()==false)&&(strcmp("#ENDGENERAL",str1)!=0)); 

    return;
  }

  GeneralBlockLoadedFlag=true; 

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("TAU",str1)==0) {
      ifile.CutInputStr(str1,str);
      tau=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();}
#if (CompilationTarget==DSMCTARGET)||(CompilationTarget==PICTARGET)||(CompilationTarget==HYBRIDTARGET) 
    else if (strcmp("MOLECULARMODEL",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("HS",str1)==0) mol.SetMolType(HS_MOL_TYPE); 
      else if (strcmp("VHS",str1)==0) mol.SetMolType(VHS_MOL_TYPE);
      else if (strcmp("VSS",str1)==0) mol.SetMolType(VSS_MOL_TYPE);
      else ifile.error();}
#endif
    else if (strcmp("TMAX",str1)==0) {
      ifile.CutInputStr(str1,str);
      tmax=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();}
    else if (strcmp("DSMC",str1)==0) { 
      ifile.CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) dsmc_flag=true;
      else if (strcmp("OFF",str1)==0) dsmc_flag=false;
      else ifile.error();}
    else if (strcmp("EXTERNALSPECIES",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) ExternalSpeciesUsingFlag=true;
      else if (strcmp("OFF",str1)==0) ExternalSpeciesUsingFlag=false;
      else ifile.error();} 
    else if (strcmp("UNIMOLECULARTRANSFORMATIONS",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) mol.UniMolTransformations=true;
      else if (strcmp("OFF",str1)==0) mol.UniMolTransformations=false;
      else ifile.error();}
    else if (strcmp("CHEMISTRY",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) chem_flag=true;
      else if (strcmp("OFF",str1)==0) chem_flag=false;
      else ifile.error();}
    else if (strcmp("IDF",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) idf_flag=true;
      else if (strcmp("OFF",str1)==0) idf_flag=false;
      else ifile.error();}
    else if (strcmp("NS",str1)==0) {
      ifile.CutInputStr(str1,str);
      NS=(unsigned char)strtol(str1,&endptr,10);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();}
    else if (strcmp("DIM",str1)==0) {
      ifile.CutInputStr(str1,str);
      DIM=(int)strtol(str1,&endptr,10);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) ifile.error();
      if ((DIM<0)||(DIM>3)) ifile.error();}
    else if (strcmp("SYMMETRY",str1)==0) {
      ifile.CutInputStr(str1,str);
      if (strcmp("NONE",str1)==0)
        SymmetryMode=no_symmetry;
      else if (strcmp("CYLINDRICAL",str1)==0)
        SymmetryMode=cylindrical_symmetry;
      else if (strcmp("SPHERICAL",str1)==0)
        SymmetryMode=spherical_symmetry;
      else ifile.error(); 
    }
    else if (strcmp("SPECIESLIST",str1)==0) {
      while (str[0]!='\0') {
        char *newSpecieSymbol=new char[100];
        int i;

        ifile.CutInputStr(str1,str);
        for (i=0;str1[i]!='\0';i++) newSpecieSymbol[i]=str1[i];
        newSpecieSymbol[i]='\0';

        mol.SpeciesList.push_back(newSpecieSymbol);
      }
    } 
    else if (strcmp("#ENDGENERAL",str1)==0) {
#if (CompilationTarget==DSMCTARGET)||(CompilationTarget==PICTARGET)||(CompilationTarget==HYBRIDTARGET) 
      mol.init(NS);
#endif
      return;}
    else ifile.error();
  }
}
//===================================================
void PARSER::parser(char* InputFile) {
  CiFileOperations ifile;
  char nonstandardBlock_endname[100]; 
  bool nonstandardBlock_flag=false;

  char str1[300],str[300];

  if (ThisThread==0) printf("InputFile: %s\n",InputFile);

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
    for (list<TExternalInputFileReader>::iterator ptr=ExternalInputFileReaderList.begin();ptr!=ExternalInputFileReaderList.end();ptr++) if ((*ptr)(str1,ifile)==true) { 
      usedDefinedBlock=true;
      break;
    }
   
    if (usedDefinedBlock==true) { 
      if (ThisThread==0) printf("Read a user defined block \"%s\"\n",str1);
    } 
    else if (strcmp("#GENERAL",str1)==0) PARSER::GeneralBlock(ifile);
    else if (strcmp("#INCLUDE",str1)==0) {
      ifile.CutInputStr(str1,ifile.init_str);
      ifile.CutInputStr(str1,ifile.init_str);
      parser(str1);
    }

#if CompilationTarget==DSMCTARGET 
    else if (strcmp("#SPECIES",str1)==0) mol.parser(ifile);
    else if (strcmp("#DSMC",str1)==0) dsmc.parser(ifile);
#endif

#if CompilationTarget==PICTARGET
    else if (strcmp("#SPECIES",str1)==0) mol.parser(ifile);
    else if (strcmp("#DSMC",str1)==0) pic.parser(ifile);
#endif

#if CompilationTarget==HYBRIDTARGET
    else if (strcmp("#SPECIES",str1)==0) mol.parser(ifile);
    else if (strcmp("#DSMC",str1)==0) hybrid.dsmc.parser(ifile);
    else if (strcmp("#EULER",str1)==0) hybrid.euler.parser(ifile);
    else if (strcmp("#HYBRID",str1)==0) hybrid.parser(ifile);
#endif

#if CompilationTarget==EULERTARGET 
    else if (strcmp("#EULER",str1)==0) euler.parser(ifile);
#endif

    else if (strcmp("#END",str1)==0) return;
    else if (strcmp("",str1)!=0) ifile.error();
  }
}  

void PARSER::readGeneralBlock(CiFileOperations& ifile) {
  GeneralBlock(ifile);
}

void PARSER::readGeneralBlock(char* InputFile) {
  CiFileOperations ifile;
  char str1[300],str[300];

  if (ThisThread==0) printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
    printf("Cannot find the input file:%s\n",InputFile);
    exit(__LINE__,__FILE__);
  }

  ifile.openfile(InputFile);
  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("#GENERAL",str1)==0) {
      readGeneralBlock(ifile);
      ifile.closefile();
      return;
    }
  }

  printf("Error: void parser_readGeneralBlock(char* InputFile), cannot find #GENERAL block \n");
  exit(__LINE__,__FILE__);
}











