//$Id$
//the library for reading and post-processing of AMPS' output files

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <algorithm>

#include "meshAMRcutcell.h"
#include "specfunc.h"
#include "global.h"

#define _POST_PROCESS_MAP_SIZE_ (1<<8)

class cPostProcess3D {
public:

  //methods for reading and reconstruction of the output file
  struct cCell {
    int n[8];
  };

  class cBlock {
  public:
    cCell ***cell;
    double xmin[3],xmax[3],dx[3];
    int nCellX,nCellY,nCellZ,id;

    void Init(int nx,int ny,int nz);
  };


  cCell *ConnectivityList;
  cBlock *Block;
  vector<cBlock*> Mesh[_POST_PROCESS_MAP_SIZE_][_POST_PROCESS_MAP_SIZE_][_POST_PROCESS_MAP_SIZE_]; //contains pointers to the blocks
  double xmin[3],xmax[3],dxBlockMin[3];
  int nNodes,nCells,nBlocks;
  int nBlockCellX,nBlockCellY,nBlockCellZ;

  //data structure
  class sDataStructure {
  public:
    int nNodes,nVariables;
    std::vector<std::string> VariableList;
    double** data;

    void CopyVariable(int iTarget,sDataStructure* Source,int iSource);
    void InitDataBuffer(int nnodes,int nvars);
  };

  sDataStructure X,data;

  //interpolation stencil
  struct cStencil {
    int Node[8];
    double Weight[8];
  };


  //trajectories of the infividual particles
  class cParticleTrajectory {
  public:
    cPostProcess3D *PostProcess3D;

    //data for all trajectories
    std::vector<std::string> VariableList;
    int TrajectoryStartingFaceOffset;
    int nTrajectoryVariables;
    int nTotalTrajectories;

    //Individual Particle Trajectory
    class cIndividualTrajectoryData {
    public:
      double *Data;
      int nDataPoints;
    };

    std::vector<cIndividualTrajectoryData> IndividualTrajectories;

    //Load and output trajectory files
    void LoadDataFile(const char* fname,const char* path);
    void AddIndividualTrajectoryData(int& nDataPoints,std::vector<double>& data);

    void PrintDataFileHeader(const char* fname);
    void AddTrajectoryDataFile(cIndividualTrajectoryData* Trajectory,int TrajectoryNumber,const char* fname);

    //plot the surface mesh showing averaged propserties of the trajectories
    void PrintSurfaceData(const char *fname);

    cParticleTrajectory() {
      nTotalTrajectories=0,nTrajectoryVariables=0,TrajectoryStartingFaceOffset=-1;
    }

  };


  cParticleTrajectory ParticleTrajectory;

  //column integration
  class cColumnIntegral {
  public:
    cPostProcess3D *PostProcess3D;

    //class that contains the user-defined integtated parametes
    struct cColumnIntegrationSet {
      int (*IntegrantVectorLength)();
      void (*IntegrantVector)(double* Data,double *Location);
      void (*PrintVariableList)(FILE* fout);
    };

    class cMap {
    public:
      cColumnIntegral *ColumnIntegral;

      //Types of the column maps integration
      void Circular(double *xObservation,double *PrimaryAxisDirection,double *SecondaryAxisDirection,double HalfMapAngle,int nAzimuthPoints,int nZenithPoints,const char *fname,cColumnIntegrationSet* IntegrationSet);
    };

    cMap Map;

  public:
    //define 3 nodes on the surface of a bounding plane; index value: 0 -> xmin component of the coordinate, 1 -> xmax component of the coordinate
    static int x0PlaneNodeIndex[6][3]; //={ {0,0,0},{1,0,0},       {0,0,0},{0,1,0},           {0,0,0},{0,0,1}};
    static int x1PlaneNodeIndex[6][3]; //={ {0,1,0},{1,1,0},       {1,0,0},{1,1,0},           {1,0,0},{1,0,1}};
    static int x2PlaneNodeIndex[6][3]; //={ {0,0,1},{1,0,1},       {0,0,1},{0,1,1},           {0,1,0},{0,1,1}};
    static int PlaneNormal[6][3]; //=     { {1,0,0},{1,0,0},       {0,1,0},{0,1,0},           {0,0,1},{0,0,1}};

    struct cBoundingBoxFace {
      double x[3];
      double e0[3];
      double e1[3];
      double Normal[3];
      double e0Length;
      double e1Length;
    };

    cBoundingBoxFace BoundingBoxFace[6];

    //operations with the column integrals
    void Init(double *xGlobalMin,double *xGlobalMax);
    bool FindIntegrationLimits(double *x0,double *l,double& IntegrationPathLength,double *xStart,double *xFinish,double *xGlobalMin,double *xGlobalMax);

    //calculate the integral
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,double *l,double IntegrationPathLength,void (*Integrand)(double*,double*));
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,double*));


    //constructor
    cColumnIntegral() {
    }
  };

  cColumnIntegral ColumnIntegral;

  //mesh operations
  cBlock* GetBlock(double *x);
  double CharacteristicCellSize(double *x);
  bool IsInsideDomain(double *x);

  //data file operations
  void LoadDataFile(const char *fname,const char* path);
  void SaveDataFile(const char *fname,sDataStructure &dat);
  void PrintVariableList();

  void SetBlockSize(int nx,int ny,int nz) {nBlockCellX=nx,nBlockCellY=ny,nBlockCellZ=nz;}
  void GetInterpolationStencil(double *x,cStencil* stencil);

  //init MPI
  int size,rank;

  inline void InitMPI() {
    MPI_Init(NULL,NULL);

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
  }

  inline void FinalizeMPI() {
    MPI_Finalize();
  }



  //constructor
  cPostProcess3D() {
    data.nVariables=0,data.nNodes=0;

    ColumnIntegral.PostProcess3D=this;
    ColumnIntegral.Map.ColumnIntegral=&ColumnIntegral;

    ParticleTrajectory.PostProcess3D=this;

  }
};
