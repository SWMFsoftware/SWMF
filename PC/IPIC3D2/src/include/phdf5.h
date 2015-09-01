
#ifndef __PHDF5_H__
#define __PHDF5_H__

#ifndef NO_HDF5
#include "ipichdf5.h"
#endif

#include "mpi.h"
#include "arraysfwd.h"
#include <string>
using namespace std;

class PHDF5fileClass{
#ifndef NO_HDF5
  public:

    PHDF5fileClass(string filestr, int nd, const int *coord, MPI_Comm mpicomm);
    PHDF5fileClass(string filestr);

    void SetDefaultGroups(void);
    void CreatePHDF5file(double *L, int *dglob, int *dlocl, bool bp);
    void ClosePHDF5file();
    void OpenPHDF5file();
    void ReadPHDF5dataset_double(string dataset, arr3_double data);
    void ReadPHDF5param();
    int  WritePHDF5dataset(string grpname, string datasetname, const_arr3_double data, int nx, int ny, int nz);

    int  getPHDF5ndim();
    int  getPHDF5ncx();
    int  getPHDF5ncy();
    int  getPHDF5ncz();

  private:

    /* The group names are fixed and must be initialized in the constructor */
    string grpnames[3];
    static const int ngrp = 3;

    /* Private variables */
    MPI_Comm comm;

    int      ndim;
    hid_t    file_id;
    double   LxLyLz  [3];  // Using dynamic allocation of these vectors caused segfault and mpi problems
    int      mpicoord[3];  // as the vectors were corrupted at the second call of the Write function.
    hsize_t  dim     [3];  // Is there a way to avoid this problem?
    hsize_t  chdim   [3];  //

    string   filename;
    bool     bparticles;
#endif
};

#endif
