//===================================================
//$Id$
//===================================================

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "dsmc.h"
#include "mol.h"
#include "data.h"

FILE* fd;
char init_str[100];
long int line=0;
//===================================================
//Separators:' ', ',', '=', ';', ':' 
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
void GetInputStr(char* str,int n){
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
    init_str[i]=str[i];
    if ((str[i]>='a')&&(str[i]<='z')) str[i]=str[i]-32;
  }
  init_str[i]='\0';
} 

//===================================================
void error() {
  printf("Error in input file (line=%i)\n",line);
  exit(0);
}

//===================================================
//===================================================
void GeneralBlock() {
  char str1[100],str[100];
  char *endptr;

  while (!feof(fd)) {
    GetInputStr(str,sizeof(str));
    CutInputStr(str1,str);
    if (strcmp("MOLECULARMODEL",str1)==0) {
      CutInputStr(str1,str);
      if (strcmp("HS",str1)==0) mol.SetMolType(HS_MOL_TYPE); 
      else if (strcmp("VHS",str1)==0) mol.SetMolType(VHS_MOL_TYPE);
      else if (strcmp("VSS",str1)==0) mol.SetMolType(VSS_MOL_TYPE);
      else error();}
    else if (strcmp("TAU",str1)==0) {
      CutInputStr(str1,str);
      tau=strtod(str1,&endptr);
      if ((str1[0]=='\0')||(endptr[0]!='\0')) error();} 
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
      else error();}
    else if (strcmp("#ENDGENERAL",str1)==0) {
      mol.init(NS);
      return;}
    else error();
  }
}
//===================================================
void parser(char* InputFile) {
  char str1[100],str[100];

  printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
    printf("Cannot find the input file:%s\n",InputFile);
    exit(0);
  }

  fd=fopen(InputFile,"r");
  while (!feof(fd)) {
    GetInputStr(str,sizeof(str));
    CutInputStr(str1,str);
    if (strcmp("#GENERAL",str1)==0) GeneralBlock();
    else if (strcmp("#SPECIES",str1)==0) mol.parser(fd,line,InputFile); 
    else if (strcmp("#DSMC",str1)==0) dsmc.parser(fd,line,InputFile);
    else if (strcmp("#END",str1)==0) return;
    else error();
  }
}  

void parser_readGeneralBlock(FILE* fin,long int& line_in) {
  fd=fin;
  line=line_in;

  GeneralBlock();
  line_in=line;
}
