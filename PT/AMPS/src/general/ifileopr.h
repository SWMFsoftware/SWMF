//===================================================
//$Id$
//===================================================

#ifndef IFILEOPR 
#define IFILEOPR


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

class CiFileOperations {
public:
  FILE* fd;
  char fname[100],init_str[1000];
  long int line;

  CiFileOperations() {
    line=-1;
  };

  FILE* openfile(char* ifile) {
    sprintf(fname,"%s",ifile); 
    line=0;
    if ((fd=fopen(ifile,"r"))==NULL) {
      printf ("FILE %s is not found\n",ifile);
      exit(0);
    } 
     
    return fd;
  }; 

  void closefile() {
    fclose(fd);
  };

  void setfile(FILE* input_fd,long int input_line,char* InputFile) {
    fd=input_fd;
    line=input_line;
    sprintf(fname,"%s",InputFile);
  };

  long int& CurrentLine() {
    return line;
  }; 

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
  }; 

  bool eof() {
    return (!feof(fd)) ? false : true;
  };

  bool GetInputStr(char* str,int n){
    int i,j;

    if (!feof(fd)) do {
      line++;
      fgets(str,n,fd);

      for (i=0;(str[i]!='\0')&&(str[i]!='\n')&&(str[i]==' ');i++);
      for (j=0;(str[i+j]!='\0')&&(str[i+j]!='\n');j++) str[j]=str[i+j];
      str[j]='\0';

      for (i=0;str[i]!='\0';i++) if ((str[i]=='!')||(str[i]=='\r')) str[i]='\0';
    } while ((str[0]=='\0')&&(!feof(fd)));

    if (feof(fd)!=0) {
      str[0]='\0',init_str[0]='\0';  
      return false; 
    }

    for(i=0;str[i]!='\0';i++) {
      init_str[i]=str[i];
      if ((str[i]>='a')&&(str[i]<='z')) str[i]=str[i]-32;
    }
    init_str[i]='\0';
 
    return true;
  }; 

  void error() {
    printf("Error in file %s, line=%ld\n",fname,line);
    printf("%s\n",init_str);
    exit(0);
  }; 
};

#endif
