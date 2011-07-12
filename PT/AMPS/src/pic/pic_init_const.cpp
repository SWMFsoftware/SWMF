

#include "pic.h"

int PIC::nTotalSpecies=0;
int PIC::ThisThread=0,PIC::nTotalThreads=1;

//the mesh parameters
double PIC::Mesh::xmin[3]={0.0,0.0,0.0},PIC::Mesh::xmax[3]={0.0,0.0,0.0};
PIC::Mesh::fLocalMeshResolution PIC::Mesh::LocalMeshResolution=NULL;
#if DIM == 3
cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR> PIC::Mesh::mesh;
#elif DIM == 2
cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR>  mesh;
#else
Not implemented yet
#endif

