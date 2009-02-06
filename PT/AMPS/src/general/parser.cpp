//===================================================
//$Id$
//===================================================

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


char init_str[100];
long int line=0;
//===================================================
//Separators:' ', ',', '=', ';', ':' '"' 
void CutInputStr(char* dest, char* src) {
  int i,j;

  for (i=0;(src[i]!='\0')&&(src[i]!=' ')&&
    (src[i]!=',')&&(src[i]!='=')&&(src[i]!=';')&&(src[i]!=':');i++) dest[i]=src[i];
    dest[i]='\0';

  for (;(src[i]!='\0')&&((src[i]==' ')||
    (src[i]==',')||(src[i]=='=')||(src[i]==';')||(src[i]==':'));i++);

  if (src[i]=='\0')
    src[0]='\0';
  else {
    for (j=0;src[j+i]!='\0';j++) src[j]=src[j+i];
    src[j]='\0';
  } 
}
//===================================================
void GetInputStr(char* str,int n,FILE* fd){
  int i,j; 
 
  do {
    line++; 
    fgets(str,n,fd);

    for (i=0;(str[i]!='\0')&&(str[i]!='\n')&&(str[i]==' ');i++);
    for (j=0;(str[i+j]!='\0')&&(str[i+j]!='\n');j++) str[j]=str[i+j];
    str[j]='\0';

    for (i=0;str[i]!='\0';i++) if (str[i]=='!') str[i]='\0'; 
  } while ((str[0]=='\0')&&(!feof(fd)));

  for(i=0;str[i]!='\0';i++) { 
    if (str[i]=='"') str[i]=' '; 
    init_str[i]=str[i];
    if ((str[i]>='a')&&(str[i]<='z')) str[i]=str[i]-32;
  }
  init_str[i]='\0';
} 

//===================================================
void error() {
  printf("Error in input file (line=%ld)\n",line);
  exit(__LINE__,__FILE__);
}

//===================================================
//===================================================
void GeneralBlock(FILE* fd) {
  static bool GeneralBlockLoadedFlag=false;
  char str1[100],str[100];
  char *endptr;

  if (GeneralBlockLoadedFlag!=false) {
    do {
      GetInputStr(str,sizeof(str),fd);
      CutInputStr(str1,str);
    }
    while((!feof(fd))&&(strcmp("#ENDGENERAL",str1)!=0)); 

    return;
  }

  GeneralBlockLoadedFlag=true; 

  while (!feof(fd)) {
    GetInputStr(str,sizeof(str),fd);
    CutInputStr(str1,str);

    if (strcmp("TAU",str1)==0) {
      CutInputStr(str1,str);
      tau=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) error();}
#if (CompilationTarget==DSMCTARGET)||(CompilationTarget==PICTARGET)||(CompilationTarget==HYBRIDTARGET) 
    else if (strcmp("MOLECULARMODEL",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("HS",str1)==0) mol.SetMolType(HS_MOL_TYPE); 
      else if (strcmp("VHS",str1)==0) mol.SetMolType(VHS_MOL_TYPE);
      else if (strcmp("VSS",str1)==0) mol.SetMolType(VSS_MOL_TYPE);
      else error();}
#endif
    else if (strcmp("TMAX",str1)==0) {
      CutInputStr(str1,str);
      tmax=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) error();}
    else if (strcmp("DSMC",str1)==0) { 
      CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) dsmc_flag=true;
      else if (strcmp("OFF",str1)==0) dsmc_flag=false;
      else error();}
    else if (strcmp("EXTERNALSPECIES",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) ExternalSpeciesUsingFlag=true;
      else if (strcmp("OFF",str1)==0) ExternalSpeciesUsingFlag=false;
      else error();} 
    else if (strcmp("UNIMOLECULARTRANSFORMATIONS",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) mol.UniMolTransformations=true;
      else if (strcmp("OFF",str1)==0) mol.UniMolTransformations=false;
      else error();}
    else if (strcmp("CHEMISTRY",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) chem_flag=true;
      else if (strcmp("OFF",str1)==0) chem_flag=false;
      else error();}
    else if (strcmp("IDF",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("ON",str1)==0) idf_flag=true;
      else if (strcmp("OFF",str1)==0) idf_flag=false;
      else error();}
    else if (strcmp("NS",str1)==0) {
      CutInputStr(str1,str);
      NS=(unsigned char)strtol(str1,&endptr,10);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) error();}
    else if (strcmp("DIM",str1)==0) {
      CutInputStr(str1,str);
      DIM=(int)strtol(str1,&endptr,10);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) error();
      if ((DIM<0)||(DIM>3)) error();}
    else if (strcmp("SYMMETRY",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("NONE",str1)==0)
        SymmetryMode=no_symmetry;
      else if (strcmp("CYLINDRICAL",str1)==0)
        SymmetryMode=cylindrical_symmetry;
      else if (strcmp("SPHERICAL",str1)==0)
        SymmetryMode=spherical_symmetry;
      else error(); 
    }
    else if (strcmp("SPECIESLIST",str1)==0) {
      while (str[0]!='\0') {
        char *newSpecieSymbol=new char[100];
        CutInputStr(str1,str);
        for (int i=0;str1[i]!='\0';i++) newSpecieSymbol[i]=str1[i];

        mol.SpeciesList.push_back(newSpecieSymbol);
      }
    } 
    else if (strcmp("#ENDGENERAL",str1)==0) {
#if (CompilationTarget==DSMCTARGET)||(CompilationTarget==PICTARGET)||(CompilationTarget==HYBRIDTARGET) 
      mol.init(NS);
#endif
      return;}
    else error();
  }
}
//===================================================
void parser(char* InputFile) {
  char nonstandardBlock_endname[100]; 
  bool nonstandardBlock_flag=false;
  FILE* fd;

  char str1[100],str[100];

  if (ThisThread==0) printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
    printf("Cannot find the input file:%s\n",InputFile);
    exit(__LINE__,__FILE__);
  }

  fd=fopen(InputFile,"r");
  while (!feof(fd)) {
    GetInputStr(str,sizeof(str),fd);
    CutInputStr(str1,str);
    if (strcmp("#GENERAL",str1)==0) GeneralBlock(fd);
    else if (strcmp("#INCLUDE",str1)==0) {
      CutInputStr(str1,init_str);
      CutInputStr(str1,init_str);
      parser(str1);
    }

#if CompilationTarget==DSMCTARGET 
    else if (strcmp("#SPECIES",str1)==0) mol.parser(fd,line,InputFile);
    else if (strcmp("#DSMC",str1)==0) dsmc.parser(fd,line,InputFile);
#endif

#if CompilationTarget==PICTARGET
    else if (strcmp("#SPECIES",str1)==0) mol.parser(fd,line,InputFile);
    else if (strcmp("#DSMC",str1)==0) pic.parser(fd,line,InputFile);
#endif

#if CompilationTarget==HYBRIDTARGET
    else if (strcmp("#SPECIES",str1)==0) mol.parser(fd,line,InputFile);
    else if (strcmp("#DSMC",str1)==0) hybrid.dsmc.parser(fd,line,InputFile);
    else if (strcmp("#EULER",str1)==0) hybrid.euler.parser(fd,line,InputFile);
    else if (strcmp("#HYBRID",str1)==0) hybrid.parser(fd,line,InputFile);
#endif

#if CompilationTarget==EULERTARGET 
    else if (strcmp("#EULER",str1)==0) euler.parser(fd,line,InputFile);
#endif

    else if (strcmp("#END",str1)==0) return;
    else if (strcmp("",str1)!=0) error();
  }
}  

void parser_readGeneralBlock(FILE* fin,long int& line_in) {
  line=line_in;

  GeneralBlock(fin);
  line_in=line;
}

void parser_readGeneralBlock(char* InputFile) {
  CiFileOperations ifile;
  FILE* fd;
  char str1[100],str[100];

  if (ThisThread==0) printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
    printf("Cannot find the input file:%s\n",InputFile);
    exit(__LINE__,__FILE__);
  }

  fd=ifile.openfile(InputFile);
  while (!feof(fd)) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("#GENERAL",str1)==0) {
      parser_readGeneralBlock(fd,ifile.CurrentLine());
      ifile.closefile();
      return;
    }
  }

  printf("Error: void parser_readGeneralBlock(char* InputFile), cannot find #GENERAL block \n");
  exit(__LINE__,__FILE__);
}











