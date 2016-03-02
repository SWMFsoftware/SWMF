/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef READPARAM_H
#define READPARAM_H

#include <iomanip>
#include <iostream>
#include <sstream>
#include <climits>
#include <cstdlib>

extern bool isProc0;

inline std::stringstream *char_to_stringstream(char* chararray, int length, int linelength, int iProc ){

  isProc0 = iProc == 0;

  std::stringstream  *ss;
  for(int i = linelength; i < length; i += linelength)
    chararray[i-1] = '\n';
  ss = new std::stringstream;
  ss->write(chararray, length);
  return(ss);
}

template <class T> 
void read_var(std::stringstream *ss, std::string description, T *var){
  if(*ss) {
    *ss >> *var;
    ss->ignore(INT_MAX, '\n');
      if(isProc0){
          std::cout<<"PC: "<<std::left<<std::setw(40)<<*var<<description<<std::endl;
          return;
      }
  }  else {
    std::cout<<" Can not find input parameter. Abort!"<<std::endl;
    std::abort();
  }
}

inline void read_var(std::stringstream *ss, std::string description, bool *var){
  std::string text;
  if(*ss) {
    *ss >> text;
    ss->ignore(INT_MAX, '\n');
    *var = false;
    if(text == "T") *var = true;
      if(isProc0){
          std::cout<<"PC: "<<std::left<<std::setw(40)<<text<<description<<std::endl;
          return;
      }
  } else {
    std::cout<<" Can not find input parameter. Abort!"<<std::endl;
    std::abort();
  }
}

inline void read_var(std::stringstream *ss, std::string description,
		     string *var){
  // Read input until \t or next line.
  std::string text;
  if(*ss) {
    string temp;
    std::getline(*ss,temp);
    std::string delimiter = "\t";
    (*var) = temp.substr(0, temp.find(delimiter));    
    if(isProc0){
      std::cout<<"PC: "<<std::left<<std::setw(40)<<(*var)<<description<<std::endl;
      return;
    }
  } else {
    std::cout<<" Can not find input parameter. Abort!"<<std::endl;
    std::abort();
  }
}

inline void get_next_command(std::stringstream *ss, std::string *id){
  *id = "";
  ss->ignore(INT_MAX, '#');
  ss->unget();
  *ss >> *id;
  if(isProc0 && id->find('#')!=std::string::npos)
    std::cout<<"\n"<<"PC: "<<*id<<std::endl;
  ss->ignore(INT_MAX, '\n');
}

inline void skip_lines(std::stringstream *ss, int nlines){
    for (int i = 0; i < nlines; i++){
    ss->ignore(INT_MAX, '\n');
  }
}

#endif
