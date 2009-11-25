//===================================================
//$Id$
//===================================================

#ifndef _PARSER_ 
#define _PARSER_ 

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>

#include "ifileopr.h"
#include "specfunc.h" 

using namespace std;

namespace PARSER {


  void parser(char*);
  void readGeneralBlock(char*);
  void readGeneralBlock(CiFileOperations&);
  void GeneralBlock(CiFileOperations&);

  //the list of external functions for reading user defined blockes of an input file
  typedef bool (*TExternalInputFileReader)(char*,CiFileOperations&);
  extern list <TExternalInputFileReader> ExternalInputFileReaderList; 

  template<class T>
  void parser(char* BlockName,char* InputFile,T& v) {
    CiFileOperations ifile;
    FILE* fd;
    char str1[100],str[100];

    if (ThisThread==0) printf("InputFile: %s\n",InputFile);

    if (access(InputFile,R_OK)!=0) {
      if (ThisThread==0) printf("Cannot find the input file:%s\n",InputFile);
      exit(__LINE__,__FILE__);
    }

    fd=ifile.openfile(InputFile);
    while (!feof(fd)) {
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(str1,str);

      if (strcmp("#GENERAL",str1)==0) readGeneralBlock(ifile);
      else if (strcmp(BlockName,str1)==0) {
        v.parser(fd,ifile.CurrentLine(),InputFile,BlockName);
        ifile.closefile();
        return;
      }
    }

    printf("Error: parser(char* BlockName,char* InputFile,T& v), cannot find block \"%s\"\n",BlockName);
    exit(__LINE__,__FILE__);
  }
}

#endif
