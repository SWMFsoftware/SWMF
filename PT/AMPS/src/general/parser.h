//===================================================
//$Id$
//===================================================

#ifndef PARSER
#define PARSER

#include "ifileopr.h"

void parser(char*);
void parser_readGeneralBlock(FILE*,long int&);

template<class T>
void parser(char* BlockName,char* InputFile,T& v) {
  CiFileOperations ifile;
  FILE* fd;
  char str1[100],str[100];

  printf("InputFile: %s\n",InputFile);

  if (access(InputFile,R_OK)!=0) {
    printf("Cannot find the input file:%s\n",InputFile);
    exit(0);
  }

  fd=ifile.openfile(InputFile);
  while (!feof(fd)) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("#GENERAL",str1)==0) parser_readGeneralBlock(fd,ifile.CurrentLine());
    else if (strcmp(BlockName,str1)==0) {
      v.parser(fd,ifile.CurrentLine(),InputFile,BlockName);
      ifile.closefile();
      return;
    }
  }

  printf("Error: parser(char* BlockName,char* InputFile,T& v), cannot fine block \"%s\"\n",BlockName);
  exit(0);
}



#endif
