//$Id$


#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"


void InitBackgroundFields() {
  //init the nackground magnetic field
  switch (_PIC_COUPLER_MODE_) {
  case _PIC_COUPLER_MODE__DATAFILE_ :

    if (PIC::CPLR::DATAFILE::BinaryFileExists("MOP-BATSRUS")==true)  {
      PIC::CPLR::DATAFILE::LoadBinaryFile("MOP-BATSRUS");
    }
    else if (_PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_) {
      // BATL reader

      if (PIC::CPLR::DATAFILE::BinaryFileExists("MOP-BATSRUS")==true)  {
        PIC::CPLR::DATAFILE::LoadBinaryFile("MOP-BATSRUS");
      }
      else {
        // initialize the reader
        #if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_
        PIC::CPLR::DATAFILE::BATSRUS::Init("3d__ful_2_t00000000_n00020000.idl");
        PIC::CPLR::DATAFILE::BATSRUS::LoadDataFile();
        #endif

        //initialize derived data
        if (PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate==true) {
          #if _PIC_COUPLER__INTERPOLATION_MODE_==_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_
          exit(__LINE__,__FILE__,"ERROR: magnetic field gradient can't be computed with 0th order interpolation method");
          #endif

          for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
            PIC::CPLR::DATAFILE::GenerateMagneticFieldGradient(node);
          }

          //Exchange derived data betwenn the boundary nodes
          PIC::Mesh::mesh.ParallelBlockDataExchange();
        }

        PIC::CPLR::DATAFILE::SaveBinaryFile("MOP-BATSRUS");
      }
    }
    else {
      exit(__LINE__,__FILE__,"ERROR: the background importing procedure is not defined");
    }

    break;
  case _PIC_COUPLER_MODE__KMAG_ :
    if (PIC::CPLR::DATAFILE::BinaryFileExists("MOP-KMAG")==true)  {
      PIC::CPLR::DATAFILE::LoadBinaryFile("MOP-KMAG");
    }
    else {
      //calculate the geomegnetic filed

      //init the model of the background magnetic field
//      T96::Init(Exosphere::SimulationStartTimeString,NULL);

      class cSetBackgroundMagneticField {
      public:
        void Set(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

          const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
          const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
          const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

          if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
            int ii,S=(kMax-kMin+1)*(jMax-jMin+1)*(iMax-iMin+1);

            if (startNode->block!=NULL) {
              #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
              #pragma omp parallel for schedule(dynamic,1) default (none) shared (PIC::Mesh::mesh,iMin,jMin,kMin,S,PIC::CPLR::DATAFILE::Offset::MagneticField,PIC::CPLR::DATAFILE::Offset::ElectricField,startNode)
              #endif

              for (ii=0;ii<S;ii++) {
                int i,j,k;
                double *xNodeMin=startNode->xmin;
                double *xNodeMax=startNode->xmax;
                double x[3],B[3],E[3]={0.0,0.0,0.0},xCell[3];
                PIC::Mesh::cDataCenterNode *CenterNode;

                //set the value of the geomagnetic field calculated at the centers of the cells
                int nd,idim;
                char *offset;

                //determine the coordinates of the cell
                int S1=ii;

                i=iMin+S1/((kMax-kMin+1)*(jMax-jMin+1));
                S1=S1%((kMax-kMin+1)*(jMax-jMin+1));

                j=jMin+S1/(kMax-kMin+1);
                k=kMin+S1%(kMax-kMin+1);

                //locate the cell
                nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
                if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
                offset=CenterNode->GetAssociatedDataBufferPointer();

                //the interpolation location
                xCell[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i));
                xCell[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j));
                xCell[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k));

                //calculate the geomagnetic field
                MOP::SaturninanSystem::Magnetosphere::GetMagneticField(B,xCell);

                if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==true) {
                  MOP::SaturninanSystem::Magnetosphere::GetElectricField(E,xCell);
                }

                //save E and B
                for (idim=0;idim<3;idim++) {
                  if (PIC::CPLR::DATAFILE::Offset::MagneticField.active==true) {
                    *((double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+idim*sizeof(double)))=B[idim];
                  }

                  if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==true) {
                    *((double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+idim*sizeof(double)))=E[idim];
                  }
                }
              }

            }
          }
          else {
            int i;
            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

            for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) Set(downNode);
          }


        }
      } SetBackgroundMagneticField;

      SetBackgroundMagneticField.Set(PIC::Mesh::mesh.rootTree);
      PIC::CPLR::DATAFILE::SaveBinaryFile("MOP-KMAG");
    }


    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }



   //init particle weight of neutral species that primary source is sputtering



 if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) {
   PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);
 }
}
