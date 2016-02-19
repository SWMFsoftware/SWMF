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


#ifndef __PARALLELIO_H__
#define __PARALLELIO_H__

#ifdef USEH5HUT
#  include "H5hut-io.h"
#endif

#ifdef PHDF5
#  include "phdf5.h"
#endif

#include "ipicfwd.h"
#include "arraysfwd.h"
#include <string>
using std::string;

void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle);
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);

void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid);
void ReadFieldsH5hut(int nspec, EMfields3D *EMf,       Collective *col, VCtopology3D *vct, Grid3DCU *grid);

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);



/**************MPI_IO*********************/
int WriteFieldsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
			float**** fieldwritebuffer,MPI_Request requestArr[4],MPI_File fhArr[4]);

int WriteMomentsVTKNonblk(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct,int cycle,
			float*** momentswritebuffer,MPI_Request requestArr[14],MPI_File fhArr[14]);

void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle);
void WriteFieldsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float**** fieldwritebuffer);
void WriteMomentsVTK(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle,float***  momentswritebuffer);
void WriteTestPclsVTK(int nspec, Grid3DCU *grid, Particles3D *part, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle, MPI_Request *testpartMPIReq, MPI_File *fh);
void ByteSwap(unsigned char * b, int n);
#endif
