//$Id$
//compare files in two directories including subdirectories 

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <dirent.h>
#include <errno.h>


#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include <iostream>
#include <fstream>


using namespace std; 


//the list of default extensions to be considered 
vector<string> ExtensionList;

void initExtensionList() {
  ExtensionList.push_back(".cpp");
  ExtensionList.push_back(".h");
  ExtensionList.push_back(".pl");
}

void help() {
  cout << endl;
  cout << "The utility looks for files that have the same name but different content. Subdirectories are included in the search." << endl;
  cout << "The format of the argument line:" << endl;
  cout << "compareSources.exe DIR1 DIR2 -e .ext" << endl;
  cout << "compareSources.exe is the executible name" << endl;
  cout << "DIR1, DIR2 are the names of directories where the search is performed" << endl;
  cout << "-e .ext allows to use additional extension in the search. In front of each additional extension should be '-e'." << endl;
  cout << "Tons of additional extension can be included in the search" << endl;
} 
 

bool CompareFiles(string strfile1,string strfile2) {
  ifstream file1stream,file2stream;
  char fname[10000];
  string fileLine1,fileLine2;

  //open file streams
  strfile1.copy(fname,strfile1.length(),0);
  fname[strfile1.length()]='\0';
  file1stream.open(fname); 
  if (file1stream.fail()==true) {
    cout << "Error: cannot open file " << strfile1 << endl;
    return false;
  }

  strfile2.copy(fname,strfile2.length(),0);
  fname[strfile2.length()]='\0';
  file2stream.open(fname);
  if (file2stream.fail()==true) {
    cout << "Error: cannot open file " << strfile2 << endl;
    file1stream.close();
    return false;
  }

  //compare files 
  bool res=true;

  while (!file1stream.eof()) {
    file1stream >> fileLine1;

    if (!file2stream.eof()) {
      file2stream >> fileLine2;  
      if (fileLine1.compare(fileLine2)!=0) {
        res=false;
        break;
      }
    }
    else {
      res=false;
      break;
    }
  }

  if (!file2stream.eof()) res=false;
  file1stream.close();
  file2stream.close();

  return res;
}


void CompareDirectories(string strDir1,string strDir2) {
  char fullname[10000];
  DIR *ptrDir[2]; 
  struct dirent *ep;
  int nList;

  list<string> subDirList[2],fileList[2]; 
  string str;

  int eLength,pos;
  vector<string>::iterator eptr;
  list<string>::iterator dir1ptr,dir2ptr,dirtmpptr,filetmpptr,file1ptr,file2ptr;

  //open the directories 
  for (nList=0;nList<2;nList++) {
    if (nList==0) {
      strDir1.copy(fullname,strDir1.length(),0); 
      fullname[strDir1.length()]='\0';
    }
    else {
     strDir2.copy(fullname,strDir2.length(),0); 
     fullname[strDir2.length()]='\0';
    }

    if ((ptrDir[nList]=opendir(fullname)) ==NULL) {

      if (nList==0) {
        cout << "ERROR: cannot open directory " << strDir1 << endl;
      }
      else {
        cout << "ERROR: cannot open directory " << strDir2 << endl;
        closedir(ptrDir[0]);
      }

      return;
    }
  }

  //read and sort the directory strDir1
  for (nList=0;nList<2;nList++) while (ep==readdir(ptrDir[nList])) if (ep->d_name[0]!='.') {
    str.assign(ep->d_name);

    switch (ep->d_type) {
    case DT_REG:

      for (eptr=ExtensionList.begin();eptr!=ExtensionList.end();eptr++) {
        eLength=eptr->length();
        pos=str.length()-eLength;

        if (pos>0) if (eptr->compare(str.substr(pos))==0) {
          fileList[nList].push_back(str);
          break;
        }
      }

      break;
    case DT_DIR: 
      subDirList[nList].push_back(str); 
    }
  } 

  //close current directories
  for (nList=0;nList<2;nList++) closedir(ptrDir[nList]);

  //check subdirectories
  for (dir1ptr=subDirList[0].begin();dir1ptr!=subDirList[0].end();) {
    dir2ptr=find(subDirList[1].begin(),subDirList[1].end(),*dir1ptr);

    if (dir2ptr!=subDirList[1].end()) {
      string dir1fullpath,dir2fullpath;

      dir1fullpath=strDir1;
      dir1fullpath.append("/");
      dir1fullpath.append(*dir1ptr); 
      
      dir2fullpath=strDir2;
      dir2fullpath.append("/");
      dir2fullpath.append(*dir1ptr);

      CompareDirectories(dir1fullpath,dir2fullpath);    

      dirtmpptr=dir1ptr;
      dirtmpptr++;
      subDirList[0].erase(dir1ptr);
      dir1ptr=dirtmpptr;

      subDirList[1].erase(dir2ptr);
    }
    else dir1ptr++;
  }

  //print the list of directories that do not match
  for (nList=0;nList<2;nList++) for (dirtmpptr=subDirList[nList].begin();dirtmpptr!=subDirList[nList].end();dirtmpptr++) {
    string dirfullpath;

    if (nList==0) dirfullpath=strDir1; else dirfullpath=strDir2;
    dirfullpath.append("/");
    dirfullpath.append(*dirtmpptr);
    
    cout << "DIRECTORY: the directory " << dirfullpath << " has not match" << endl;
  } 
 

  //compare files
  for (file1ptr=fileList[0].begin();file1ptr!=fileList[0].end();) {
    file2ptr=find(fileList[1].begin(),fileList[1].end(),*file1ptr);

    if (file2ptr!=fileList[1].end()) {
      string file1fullpath,file2fullpath;

      file1fullpath=strDir1;
      file1fullpath.append("/");
      file1fullpath.append(*file1ptr);

      file2fullpath=strDir2;
      file2fullpath.append("/");
      file2fullpath.append(*file1ptr);

      if (CompareFiles(file1fullpath,file2fullpath)==true) {    
        filetmpptr=file1ptr;
        filetmpptr++;
        fileList[0].erase(file1ptr);
        file1ptr=filetmpptr;
      }
      else file1ptr++;

      fileList[1].erase(file2ptr);
    }
    else file1ptr++;
  }


  //print the list of files that do not match
  for (nList=0;nList<2;nList++) for (filetmpptr=fileList[nList].begin();filetmpptr!=fileList[nList].end();filetmpptr++) {
    string filefullpath;

    if (nList==0) filefullpath=strDir1; else filefullpath=strDir2;
    filefullpath.append("/");
    filefullpath.append(*filetmpptr);

    cout << "FILE: the file " << filefullpath << " does not match" << endl;
  }


  
}


int main(int argc,char **argv) {
  string rootDir[2],str;
  int nDir=0,nArg=0; 

  initExtensionList();

  //read the input line
  for (nArg=1;nArg<argc;nArg++) { 
    str.assign(argv[nArg]);

    if (str.compare("-e")==0) {
      if (nArg==argc-1) {
        cout << "Error: not enought arguments in the argument line" << endl;
        help();
        return 0;
      }
      else {
        str.assign(argv[++nArg]);
        ExtensionList.push_back(str);
      }
    }
    else if (str.compare("-help")==0) {
      help();
      return 0;
    }
    else {
      if (nDir==2) {
        cout << "Too many seach directories in the argument line" << endl;
        help();
        return 0; 
      }
      else rootDir[nDir++].assign(str); 
    } 
  }

  if (nDir!=2) {
    cout << "ERROR: two directories must be present in the argument line" << endl;
    help();
    return 0;
  } 

  cout << "The following extensions are considered:" << endl;
  for (vector<string>::iterator ptr=ExtensionList.begin();ptr!=ExtensionList.end();ptr++) cout << *ptr << "  ";
  cout << endl; 

  

  CompareDirectories(rootDir[0],rootDir[1]); 

  return 1;
}

