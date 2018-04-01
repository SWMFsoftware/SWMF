/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta.
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010):
 *1509-1519.'
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

#include "Com3DNonblk.h"
#include "my_mpi.h"
#include "Alloc.h"

// isCenterFlag: 1 communicateCenter; 0 communicateNode
void NBDerivedHaloComm(int nx, int ny, int nz, int n4, double ****vectorIn_GI,
                       const VirtualTopology3D *vct, EMfields3D *EMf,
                       bool isCenterFlag, bool isFaceOnlyFlag, bool needInterp,
                       bool isParticle, bool doSwapGhostCell = false,
                       bool isMassMatrix = false) {

  MPI_Errhandler_set(MPI_COMM_MYSIM, MPI_ERRORS_RETURN);
  MPI_Status stat[12];
  MPI_Request reqList[12]; // at most 6 requests x 2 (send recv)
  int communicationCnt[6] = { 0, 0, 0, 0, 0, 0 }; // 1 if there communication on
                                                  // that dir
  int recvcnt = 0, sendcnt = 0;                   // request counter

  const int tag_XL = 1, tag_YL = 2, tag_ZL = 3, tag_XR = 4, tag_YR = 5,
            tag_ZR = 6; // To address same rank as left and right neighbour in
                        // periodic case
  const int myrank = vct->getCartesian_rank();
  const MPI_Comm comm =
      isParticle ? vct->getParticleComm() : vct->getFieldComm();
  const int right_neighborX =
      isParticle ? vct->getXright_neighbor_P() : vct->getXright_neighbor();
  const int left_neighborX =
      isParticle ? vct->getXleft_neighbor_P() : vct->getXleft_neighbor();
  const int right_neighborY =
      isParticle ? vct->getYright_neighbor_P() : vct->getYright_neighbor();
  const int left_neighborY =
      isParticle ? vct->getYleft_neighbor_P() : vct->getYleft_neighbor();
  const int right_neighborZ =
      isParticle ? vct->getZright_neighbor_P() : vct->getZright_neighbor();
  const int left_neighborZ =
      isParticle ? vct->getZleft_neighbor_P() : vct->getZleft_neighbor();
  bool isCenterDim =
      (isCenterFlag && !needInterp); // This is to address moment interp use
                                     // nodes dimension but exchange center face



  // The communicating grid is: (nxc+iLay)*(nyc+iLay)*(nzc*iLay), and 
  // there are n4 variables on each node/cell. The corresponding 
  // data types are already defined in advanced. 
  double iLayAdd;
  if(isCenterDim)
    iLayAdd = 0; 
  else
    iLayAdd = 1; 

  // Added one extra ghost cell layer at each side to receive ghost cell
  // values from nearyby processors. 
  if(doSwapGhostCell) iLayAdd = 2; 

  const MPI_Datatype yzFacetype = EMf->getYZFacetype(iLayAdd, n4);
  const MPI_Datatype xzFacetype = EMf->getXZFacetype(iLayAdd, n4);
  const MPI_Datatype xyFacetype = EMf->getXYFacetype(iLayAdd, n4);
  const MPI_Datatype xEdgetype = EMf->getXEdgetype(iLayAdd, n4);
  const MPI_Datatype yEdgetype = EMf->getYEdgetype(iLayAdd, n4);
  const MPI_Datatype zEdgetype = EMf->getZEdgetype(iLayAdd, n4);
  const MPI_Datatype xEdgetype2 = EMf->getXEdgetype2(iLayAdd, n4);
  const MPI_Datatype yEdgetype2 = EMf->getYEdgetype2(iLayAdd, n4);
  const MPI_Datatype zEdgetype2 = EMf->getZEdgetype2(iLayAdd, n4);
  const MPI_Datatype cornertype = EMf->getCornertype(iLayAdd, n4);

  double ****vector_GI = NULL;
  if (!doSwapGhostCell) {
    vector_GI = vectorIn_GI;
  } else {
    // Creating a new tempary array so that there is space to store
    // ghost cells from nearby processors. 
    nx +=2; ny +=2; nz +=2; 
    vector_GI = newArr4(double, nx, ny, nz,n4);   
    for (int i = 1; i < nx-1; i++)
      for (int j = 1; j < ny-1; j++)
        for (int k = 1; k < nz-1; k++)
          for (int i4 = 0; i4 < n4; i4++) {
            vector_GI[i][j][k][i4] = vectorIn_GI[i-1][j-1][k-1][i4];
          }
  }


  int offset = (isCenterFlag? 0:1);

  void *vp = NULL;


  // Face Exchange as long as neighbor exits on that direction
  // Tag is based on the sender: XL for sender = XR for receiver
  // Post Recv before Send
  if (left_neighborX != MPI_PROC_NULL && left_neighborX != myrank) {
    vp = &vector_GI[0][1][1][0];
    MPI_Irecv(vp, 1, yzFacetype, left_neighborX, tag_XR, comm,
              &reqList[recvcnt++]);
    communicationCnt[0] = 1;
  }
  if (right_neighborX != MPI_PROC_NULL && right_neighborX != myrank) {
    vp = &vector_GI[nx - 1][1][1][0];
    MPI_Irecv(vp, 1, yzFacetype, right_neighborX, tag_XL, comm,
              &reqList[recvcnt++]);
    communicationCnt[1] = 1;
  }
  if (left_neighborY != MPI_PROC_NULL && left_neighborY != myrank) {
    vp = &vector_GI[1][0][1][0];
    MPI_Irecv(vp, 1, xzFacetype, left_neighborY, tag_YR, comm,
              &reqList[recvcnt++]);
    communicationCnt[2] = 1;
  }
  if (right_neighborY != MPI_PROC_NULL && right_neighborY != myrank) {
    vp = &vector_GI[1][ny - 1][1][0];
    MPI_Irecv(vp, 1, xzFacetype, right_neighborY, tag_YL, comm,
              &reqList[recvcnt++]);
    communicationCnt[3] = 1;
  }
  if (left_neighborZ != MPI_PROC_NULL && left_neighborZ != myrank) {
    vp = &vector_GI[1][1][0][0];
    MPI_Irecv(vp, 1, xyFacetype, left_neighborZ, tag_ZR, comm,
              &reqList[recvcnt++]);
    communicationCnt[4] = 1;
  }
  if (right_neighborZ != MPI_PROC_NULL && right_neighborZ != myrank) {
    vp = &vector_GI[1][1][nz - 1][0];
    MPI_Irecv(vp, 1, xyFacetype, right_neighborZ, tag_ZL, comm,
              &reqList[recvcnt++]);
    communicationCnt[5] = 1;
  }

  sendcnt = recvcnt;

  if (communicationCnt[0] == 1) {
    vp = &vector_GI[1 + offset][1][1][0];
    MPI_Isend(vp, 1, yzFacetype, left_neighborX, tag_XL, comm,
              &reqList[sendcnt++]);
  }
  if (communicationCnt[1] == 1) {
    vp = &vector_GI[nx - 2 - offset][1][1][0];
    MPI_Isend(vp, 1, yzFacetype, right_neighborX, tag_XR, comm,
              &reqList[sendcnt++]);
  }
  if (communicationCnt[2] == 1) {
    vp = &vector_GI[1][1 + offset][1][0];
    MPI_Isend(vp, 1, xzFacetype, left_neighborY, tag_YL, comm,
              &reqList[sendcnt++]);
  }
  if (communicationCnt[3] == 1) {
    vp = &vector_GI[1][ny - 2 - offset][1][0];
    MPI_Isend(vp, 1, xzFacetype, right_neighborY, tag_YR, comm,
              &reqList[sendcnt++]);
  }
  if (communicationCnt[4] == 1) {
    vp = &vector_GI[1][1][1 + offset][0];
    MPI_Isend(vp, 1, xyFacetype, left_neighborZ, tag_ZL, comm,
              &reqList[sendcnt++]);
  }
  if (communicationCnt[5] == 1) {
    vp = &vector_GI[1][1][nz - 2 - offset][0];
    MPI_Isend(vp, 1, xyFacetype, right_neighborZ, tag_ZR, comm,
              &reqList[sendcnt++]);
  }
  assert_eq(recvcnt, sendcnt - recvcnt);

  // Buffer swap if any (done before waiting for receiving msg done to delay
  // sync)
  if (right_neighborX == myrank && left_neighborX == myrank) {
    for (int iy = 1; iy < ny - 1; iy++)
      for (int iz = 1; iz < nz - 1; iz++)
        for (int i4 = 0; i4 < n4; i4++) {
          vector_GI[0][iy][iz][i4] = vector_GI[nx - 2 - offset][iy][iz][i4];	  
          vector_GI[nx - 1][iy][iz][i4] = vector_GI[1 + offset][iy][iz][i4];
        }
  }
  if (right_neighborY == myrank && left_neighborY == myrank) {
    for (int ix = 1; ix < nx - 1; ix++)
      for (int iz = 1; iz < nz - 1; iz++)
        for (int i4 = 0; i4 < n4; i4++) {
          vector_GI[ix][0][iz][i4] = vector_GI[ix][ny - 2 - offset][iz][i4];
          vector_GI[ix][ny - 1][iz][i4] = vector_GI[ix][1 + offset][iz][i4];
        }
  }
  if (right_neighborZ == myrank && left_neighborZ == myrank) {
    for (int ix = 1; ix < nx - 1; ix++)
      for (int iy = 1; iy < ny - 1; iy++)
        for (int i4 = 0; i4 < n4; i4++) {
          vector_GI[ix][iy][0][i4] = vector_GI[ix][iy][nz - 2 - offset][i4];
          vector_GI[ix][iy][nz - 1][i4] = vector_GI[ix][iy][1 + offset][i4];
        }
  }
  // Need to finish receiving + sending before edge exchange
  if (sendcnt > 0) {
    MPI_Waitall(sendcnt, &reqList[0], &stat[0]);
    bool stopFlag = false;
    // The following code is useful for debug, but error_code
    // will be always a random number on Pleiades (with mpt).
    // for(int si=0;si< sendcnt;si++){
    //   int error_code = stat[si].MPI_ERROR;
    // if (error_code != MPI_SUCCESS) {
    //     stopFlag = true;
    //     char error_string[100];
    //     int length_of_error_string, error_class;
    //     MPI_Error_class(error_code, &error_class);
    //     MPI_Error_string(error_class, error_string, &length_of_error_string);
    //     dprintf("MPI_Waitall error at %d  %s\n",si, error_string);
    // }
    //}
    if (stopFlag)
      exit(EXIT_FAILURE);
  }

  double tmp; 

  if (!isFaceOnlyFlag) {

    // Exchange yEdge only when Z X neighbours exist
    // if Zleft + Zright, use merged yEdgeType2 to send two edges in one msg
    // Otherwise, only send one yEdge
    recvcnt = 0, sendcnt = 0;
    if (communicationCnt[0] == 1) {
      if (communicationCnt[4] == 1 && communicationCnt[5] == 1) {
        vp = &vector_GI[0][1][0][0];
        MPI_Irecv(vp, 1, yEdgetype2, left_neighborX, tag_XR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[4] == 1) {
        vp = &vector_GI[0][1][0][0];
        MPI_Irecv(vp, 1, yEdgetype, left_neighborX, tag_XR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[5] == 1) {
        vp = &vector_GI[0][1][nz - 1][0];
        MPI_Irecv(vp, 1, yEdgetype, left_neighborX, tag_XR, comm,
                  &reqList[recvcnt++]);
      }
    }
    if (communicationCnt[1] == 1) {
      if (communicationCnt[4] == 1 && communicationCnt[5] == 1) {
        vp = &vector_GI[nx - 1][1][0][0];
        MPI_Irecv(vp, 1, yEdgetype2, right_neighborX, tag_XL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[4] == 1) {
        vp = &vector_GI[nx - 1][1][0][0];
        MPI_Irecv(vp, 1, yEdgetype, right_neighborX, tag_XL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[5] == 1) {
        vp = &vector_GI[nx - 1][1][nz - 1][0];
        MPI_Irecv(vp, 1, yEdgetype, right_neighborX, tag_XL, comm,
                  &reqList[recvcnt++]);
      }
    }
    // Exchange zEdge only when X Y neighbours exist
    // if Xleft + Xright, use merged zEdgeType2
    // Otherwise, only send one zEdge
    if (communicationCnt[2] == 1) {
      if (communicationCnt[0] == 1 && communicationCnt[1] == 1) {
        vp = &vector_GI[0][0][1][0];
        MPI_Irecv(vp, 1, zEdgetype2, left_neighborY, tag_YR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[0] == 1) {
        vp = &vector_GI[0][0][1][0];
        MPI_Irecv(vp, 1, zEdgetype, left_neighborY, tag_YR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 1][0][1][0];
        MPI_Irecv(vp, 1, zEdgetype, left_neighborY, tag_YR, comm,
                  &reqList[recvcnt++]);
      }
    }
    if (communicationCnt[3] == 1) {
      if (communicationCnt[0] == 1 && communicationCnt[1] == 1) {
        vp = &vector_GI[0][ny - 1][1][0];
        MPI_Irecv(vp, 1, zEdgetype2, right_neighborY, tag_YL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[0] == 1) {
        vp = &vector_GI[0][ny - 1][1][0];
        MPI_Irecv(vp, 1, zEdgetype, right_neighborY, tag_YL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 1][ny - 1][1][0];
        MPI_Irecv(vp, 1, zEdgetype, right_neighborY, tag_YL, comm,
                  &reqList[recvcnt++]);
      }
    }
    // Exchange xEdge only when Y Z neighbours exist
    // if Yleft + Yright exist, use merged xEdgeType2
    // Otherwise, only send one xEdge
    if (communicationCnt[4] == 1) {
      if (communicationCnt[2] == 1 && communicationCnt[3] == 1) {
        vp = &vector_GI[1][0][0][0];
        MPI_Irecv(vp, 1, xEdgetype2, left_neighborZ, tag_ZR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[2] == 1) {
        vp = &vector_GI[1][0][0][0];
        MPI_Irecv(vp, 1, xEdgetype, left_neighborZ, tag_ZR, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[3] == 1) {
        vp = &vector_GI[1][ny - 1][0][0];
        MPI_Irecv(vp, 1, xEdgetype, left_neighborZ, tag_ZR, comm,
                  &reqList[recvcnt++]);
      }
    }
    if (communicationCnt[5] == 1) {
      if (communicationCnt[2] == 1 && communicationCnt[3] == 1) {
        vp = &vector_GI[1][0][nz - 1][0];
        MPI_Irecv(vp, 1, xEdgetype2, right_neighborZ, tag_ZL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[2] == 1) {
        vp = &vector_GI[1][0][nz - 1][0];
        MPI_Irecv(vp, 1, xEdgetype, right_neighborZ, tag_ZL, comm,
                  &reqList[recvcnt++]);
      } else if (communicationCnt[3] == 1) {
        vp = &vector_GI[1][ny - 1][nz - 1][0];
        MPI_Irecv(vp, 1, xEdgetype, right_neighborZ, tag_ZL, comm,
                  &reqList[recvcnt++]);
      }
    }

    sendcnt = recvcnt;

    if (communicationCnt[0] == 1) {
      if (communicationCnt[4] == 1 && communicationCnt[5] == 1) {
        vp = &vector_GI[1 + offset][1][0][0];
        MPI_Isend(vp, 1, yEdgetype2, left_neighborX, tag_XL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[4] == 1) {
        vp = &vector_GI[1 + offset][1][0][0];
        MPI_Isend(vp, 1, yEdgetype, left_neighborX, tag_XL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[5] == 1) {
        vp = &vector_GI[1 + offset][1][nz - 1][0];
        MPI_Isend(vp, 1, yEdgetype, left_neighborX, tag_XL, comm,
                  &reqList[sendcnt++]);
      }
    }
    if (communicationCnt[1] == 1) {
      if (communicationCnt[4] == 1 && communicationCnt[5] == 1) {
        vp = &vector_GI[nx - 2 - offset][1][0][0];
        MPI_Isend(vp, 1, yEdgetype2, right_neighborX, tag_XR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[4] == 1) {
        vp = &vector_GI[nx - 2 - offset][1][0][0];
        MPI_Isend(vp, 1, yEdgetype, right_neighborX, tag_XR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[5] == 1) {
        vp = &vector_GI[nx - 2 - offset][1][nz - 1][0];
        MPI_Isend(vp, 1, yEdgetype, right_neighborX, tag_XR, comm,
                  &reqList[sendcnt++]);
      }
    }
    if (communicationCnt[2] == 1) {
      if (communicationCnt[0] == 1 && communicationCnt[1] == 1) {
        vp = &vector_GI[0][1 + offset][1][0];
        MPI_Isend(vp, 1, zEdgetype2, left_neighborY, tag_YL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[0] == 1) {
        vp = &vector_GI[0][1 + offset][1][0];
        MPI_Isend(vp, 1, zEdgetype, left_neighborY, tag_YL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 1][1 + offset][1][0];
        MPI_Isend(vp, 1, zEdgetype, left_neighborY, tag_YL, comm,
                  &reqList[sendcnt++]);
      }
    }
    if (communicationCnt[3] == 1) {
      if (communicationCnt[0] == 1 && communicationCnt[1] == 1) {
        vp = &vector_GI[0][ny - 2 - offset][1][0];
        MPI_Isend(vp, 1, zEdgetype2, right_neighborY, tag_YR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[0] == 1) {
        vp = &vector_GI[0][ny - 2 - offset][1][0];
        MPI_Isend(vp, 1, zEdgetype, right_neighborY, tag_YR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 1][ny - 2 - offset][1][0];
        MPI_Isend(vp, 1, zEdgetype, right_neighborY, tag_YR, comm,
                  &reqList[sendcnt++]);
      }
    }
    if (communicationCnt[4] == 1) {
      if (communicationCnt[2] == 1 && communicationCnt[3] == 1) {
        vp = &vector_GI[1][0][1 + offset][0];
        MPI_Isend(vp, 1, xEdgetype2, left_neighborZ, tag_ZL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[2] == 1) {
        vp = &vector_GI[1][0][1 + offset][0];
        MPI_Isend(vp, 1, xEdgetype, left_neighborZ, tag_ZL, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[3] == 1) {
        vp = &vector_GI[1][ny - 1][1 + offset][0];
        MPI_Isend(vp, 1, xEdgetype, left_neighborZ, tag_ZL, comm,
                  &reqList[sendcnt++]);
      }
    }
    if (communicationCnt[5] == 1) {
      if (communicationCnt[2] == 1 && communicationCnt[3] == 1) {
        vp = &vector_GI[1][0][nz - 2 - offset][0];
        MPI_Isend(vp, 1, xEdgetype2, right_neighborZ, tag_ZR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[2] == 1) {
        vp = &vector_GI[1][0][nz - 2 - offset][0];
        MPI_Isend(vp, 1, xEdgetype, right_neighborZ, tag_ZR, comm,
                  &reqList[sendcnt++]);
      } else if (communicationCnt[3] == 1) {
        vp = &vector_GI[1][ny - 1][nz - 2 - offset][0];
        MPI_Isend(vp, 1, xEdgetype, right_neighborZ, tag_ZR, comm,
                  &reqList[sendcnt++]);
      }
    }

    assert_eq(recvcnt, sendcnt - recvcnt);

    // Swap Local Edges

    if (right_neighborX == myrank && left_neighborX == myrank) {
      if (right_neighborZ != MPI_PROC_NULL) {
        for (int iy = 1; iy < ny - 1; iy++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[nx - 2 - offset][iy][nz - 1][i4];
            vector_GI[nx - 1][iy][nz - 1][i4] =
	      vector_GI[1 + offset][iy][nz - 1][i4];
	    vector_GI[0][iy][nz - 1][i4] = tmp; 
          }
      }
      if (left_neighborZ != MPI_PROC_NULL) {
        for (int iy = 1; iy < ny - 1; iy++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = vector_GI[nx - 2 - offset][iy][0][i4];
            vector_GI[nx - 1][iy][0][i4] = vector_GI[1 + offset][iy][0][i4];
	    vector_GI[0][iy][0][i4] = tmp; 
          }
      }
      if (right_neighborY != MPI_PROC_NULL) {
        for (int iz = 1; iz < nz - 1; iz++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
                vector_GI[nx - 2 - offset][ny - 1][iz][i4];
            vector_GI[nx - 1][ny - 1][iz][i4] =
                vector_GI[1 + offset][ny - 1][iz][i4];
	    vector_GI[0][ny - 1][iz][i4] = tmp; 
          }
      }
      if (left_neighborY != MPI_PROC_NULL) {
        for (int iz = 1; iz < nz - 1; iz++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[nx - 2 - offset][0][iz][i4];
            vector_GI[nx - 1][0][iz][i4] = vector_GI[1 + offset][0][iz][i4];
	    vector_GI[0][0][iz][i4] = tmp; 
          }
      }
    }

    if (right_neighborY == myrank && left_neighborY == myrank) {
      if (right_neighborX != MPI_PROC_NULL) {
        for (int iz = 1; iz < nz - 1; iz++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
                vector_GI[nx - 1][ny - 2 - offset][iz][i4];
            vector_GI[nx - 1][ny - 1][iz][i4] =
                vector_GI[nx - 1][1 + offset][iz][i4];
	    vector_GI[nx - 1][0][iz][i4] = tmp; 
          }
      }  
      if (left_neighborX != MPI_PROC_NULL) {
        for (int iz = 1; iz < nz - 1; iz++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[0][ny - 2 - offset][iz][i4];
            vector_GI[0][ny - 1][iz][i4] = vector_GI[0][1 + offset][iz][i4];
	    vector_GI[0][0][iz][i4]  = tmp; 
          }
      }
      if (right_neighborZ != MPI_PROC_NULL) {
        for (int ix = 1; ix < nx - 1; ix++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
                vector_GI[ix][ny - 2 - offset][nz - 1][i4];
            vector_GI[ix][ny - 1][nz - 1][i4] =
                vector_GI[ix][1 + offset][nz - 1][i4];
	    vector_GI[ix][0][nz - 1][i4] = tmp; 
          }
      }
      if (left_neighborZ != MPI_PROC_NULL) {
        for (int ix = 1; ix < nx - 1; ix++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[ix][ny - 2 - offset][0][i4];
            vector_GI[ix][ny - 1][0][i4] = vector_GI[ix][1 + offset][0][i4];
	    vector_GI[ix][0][0][i4] = tmp; 
          }
      }
    }

    if (right_neighborZ == myrank && left_neighborZ == myrank) {
      if (right_neighborY != MPI_PROC_NULL) {
        for (int ix = 1; ix < nx - 1; ix++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
                vector_GI[ix][ny - 1][nz - 2 - offset][i4];
            vector_GI[ix][ny - 1][nz - 1][i4] =
                vector_GI[ix][ny - 1][1 + offset][i4];
	    vector_GI[ix][ny - 1][0][i4] = tmp; 
          }
      }
      if (left_neighborY != MPI_PROC_NULL) {
        for (int ix = 1; ix < nx - 1; ix++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[ix][0][nz - 2 - offset][i4];
            vector_GI[ix][0][nz - 1][i4] = vector_GI[ix][0][1 + offset][i4];
	    vector_GI[ix][0][0][i4] = tmp; 
          }
      }
      if (right_neighborX != MPI_PROC_NULL) {
        for (int iy = 1; iy < ny - 1; iy++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
                vector_GI[nx - 1][iy][nz - 2 - offset][i4];
            vector_GI[nx - 1][iy][nz - 1][i4] =
                vector_GI[nx - 1][iy][1 + offset][i4];
	    vector_GI[nx - 1][iy][0][i4] = tmp; 
          }
      }
      if (left_neighborX != MPI_PROC_NULL) {
        for (int iy = 1; iy < ny - 1; iy++)
          for (int i4 = 0; i4 < n4; i4++) {
            tmp = 
	      vector_GI[0][iy][nz - 2 - offset][i4];
            vector_GI[0][iy][nz - 1][i4] = vector_GI[0][iy][1 + offset][i4];
	    vector_GI[0][iy][0][i4] = tmp; 
          }
      }
    }
    // Need to finish receiving edges for corner exchange
    if (sendcnt > 0) {
      MPI_Waitall(sendcnt, &reqList[0], &stat[0]);
      bool stopFlag = false;
      // The following code is useful for debug, but error_code
      // will be always a random number on Pleiades (with mpt).
      // for(int si=0;si< sendcnt;si++){
      // 	int error_code = stat[si].MPI_ERROR;
      // 	if (error_code != MPI_SUCCESS) {
      // 		stopFlag = true;
      // 		char error_string[100];
      // 		int length_of_error_string, error_class;

      // 		MPI_Error_class(error_code, &error_class);
      // 		MPI_Error_string(error_class, error_string,
      // &length_of_error_string);
      // 		dprintf("MPI_Waitall error at %d  %s\n",si,
      // error_string);
      // 	}
      // }

      if (stopFlag)
        exit(EXIT_FAILURE);
      ;
    }

    // Corner Exchange only needed if XYZ neighbours all exist
    // 4 corners communicated in one message
    // Assume Non-periodic will be handled in BC
    // Define corner types for X communication
    recvcnt = 0, sendcnt = 0;
    if ((communicationCnt[2] == 1 || communicationCnt[3] == 1) &&
        (communicationCnt[4] == 1 || communicationCnt[5] == 1)) {
      // if XLeft exists, send 4 corners to XLeft
      if (communicationCnt[0] == 1) {
        vp = &vector_GI[0][0][0][0];
        MPI_Irecv(vp, 1, cornertype, left_neighborX, tag_XR, comm,
                  &reqList[recvcnt++]);
      }
      // if XRight exist
      if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 1][0][0][0];
        MPI_Irecv(vp, 1, cornertype, right_neighborX, tag_XL, comm,
                  &reqList[recvcnt++]);
      }

      sendcnt = recvcnt;

      if (communicationCnt[0] == 1) {
        vp = &vector_GI[1 + offset][0][0][0];
        MPI_Isend(vp, 1, cornertype, left_neighborX, tag_XL, comm,
                  &reqList[sendcnt++]);
      }
      if (communicationCnt[1] == 1) {
        vp = &vector_GI[nx - 2 - offset][0][0][0];
        MPI_Isend(vp, 1, cornertype, right_neighborX, tag_XR, comm,
                  &reqList[sendcnt++]);
      }
    }

    assert_eq(recvcnt, sendcnt - recvcnt);

    // Delay local data copy
    if (left_neighborX == myrank && right_neighborX == myrank) {
      if ((left_neighborY != MPI_PROC_NULL) &&
          (left_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          vector_GI[0][0][0][i4] = vector_GI[nx - 2 - offset][0][0][i4];
          vector_GI[nx - 1][0][0][i4] = vector_GI[1 + offset][0][0][i4];
        }
      }

      if ((left_neighborY != MPI_PROC_NULL) &&
          (right_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 2 - offset][0][nz - 1][i4];
          vector_GI[nx - 1][0][nz - 1][i4] =
              vector_GI[1 + offset][0][nz - 1][i4];
	  vector_GI[0][0][nz - 1][i4] = tmp; 
        }
      }

      if ((right_neighborY != MPI_PROC_NULL) &&
          (left_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 2 - offset][ny - 1][0][i4];
          vector_GI[nx - 1][ny - 1][0][i4] =
              vector_GI[1 + offset][ny - 1][0][i4];
	  vector_GI[0][ny - 1][0][i4] = tmp; 
        }
      }

      if ((right_neighborY != MPI_PROC_NULL) &&
          (right_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 2 - offset][ny - 1][nz - 1][i4];
          vector_GI[nx - 1][ny - 1][nz - 1][i4] =
              vector_GI[1 + offset][ny - 1][nz - 1][i4];
	  vector_GI[0][ny - 1][nz - 1][i4] = tmp; 
        }
      }
    } else if (left_neighborY == myrank && right_neighborY == myrank) {
      if ((left_neighborX != MPI_PROC_NULL) &&
          (left_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
	    vector_GI[0][ny - 2 - offset][0][i4];
          vector_GI[0][ny - 1][0][i4] = vector_GI[0][1 + offset][0][i4];
	  vector_GI[0][0][0][i4] = tmp; 
        }
      }

      if ((left_neighborX != MPI_PROC_NULL) &&
          (right_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[0][ny - 2 - offset][nz - 1][i4];
          vector_GI[0][ny - 1][nz - 1][i4] =
              vector_GI[0][1 + offset][nz - 1][i4];
	  vector_GI[0][0][nz - 1][i4] = tmp; 
        }
      }
      if ((right_neighborX != MPI_PROC_NULL) &&
          (left_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 1][ny - 2 - offset][0][i4];
          vector_GI[nx - 1][ny - 1][0][i4] =
              vector_GI[nx - 1][1 + offset][0][i4];
	  vector_GI[nx - 1][0][0][i4] = tmp; 
        }
      }
      if ((right_neighborX != MPI_PROC_NULL) &&
          (right_neighborZ != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {	  
          tmp = 
              vector_GI[nx - 1][ny - 2 - offset][nz - 1][i4];
          vector_GI[nx - 1][ny - 1][nz - 1][i4] =
              vector_GI[nx - 1][1 + offset][nz - 1][i4];
	  vector_GI[nx - 1][0][nz - 1][i4] = tmp; 
        }
      }
    } else if (left_neighborZ == myrank && right_neighborZ == myrank) {
      if ((left_neighborY != MPI_PROC_NULL) &&
          (left_neighborX != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
	    vector_GI[0][0][nz - 2 - offset][i4];
          vector_GI[0][0][nz - 1][i4] = vector_GI[0][0][1 + offset][i4];
	  vector_GI[0][0][0][i4] = tmp; 
        }
      }

      if ((left_neighborY != MPI_PROC_NULL) &&
          (right_neighborX != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 1][0][nz - 2 - offset][i4];
          vector_GI[nx - 1][0][nz - 1][i4] =
              vector_GI[nx - 1][0][1 + offset][i4];
	  vector_GI[nx - 1][0][0][i4] = tmp; 
        }
      }
      if ((right_neighborY != MPI_PROC_NULL) &&
          (left_neighborX != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[0][ny - 1][nz - 2 - offset][i4];
          vector_GI[0][ny - 1][nz - 1][i4] =
              vector_GI[0][ny - 1][1 + offset][i4];
	  vector_GI[0][ny - 1][0][i4] = tmp; 
        }
      }

      if ((right_neighborY != MPI_PROC_NULL) &&
          (right_neighborX != MPI_PROC_NULL)) {
        for (int i4 = 0; i4 < n4; i4++) {
          tmp = 
              vector_GI[nx - 1][ny - 1][nz - 2 - offset][i4];
          vector_GI[nx - 1][ny - 1][nz - 1][i4] =
              vector_GI[nx - 1][ny - 1][1 + offset][i4];
	  vector_GI[nx - 1][ny - 1][0][i4] = tmp; 
        }
      }
    }

    if (sendcnt > 0) {
      MPI_Waitall(sendcnt, &reqList[0], &stat[0]);
      bool stopFlag = false;
      // for(int si=0;si< sendcnt;si++){
      // 	int error_code = stat[si].MPI_ERROR;
      // 	if (error_code != MPI_SUCCESS) {
      // 		stopFlag = true;
      // 		char error_string[100];
      // 		int length_of_error_string, error_class;

      // 		MPI_Error_class(error_code, &error_class);
      // 		MPI_Error_string(error_class, error_string,
      // &length_of_error_string);
      // 		dprintf("MPI_Waitall error at %d  %s\n",si,
      // error_string);
      // 	}
      // }

      if (stopFlag)
        exit(EXIT_FAILURE);
    }
  }

  // if this is NodeInterpolation operation
  if (needInterp || doSwapGhostCell) {
    addBoundary(nx, ny, nz, n4, vector_GI, vct, doSwapGhostCell);
  }
  
  if (doSwapGhostCell){
    for (int i = 1; i < nx-1; i++)
      for (int j = 1; j < ny-1; j++)
        for (int k = 1; k < nz-1; k++)
          for (int i4 = 0; i4 < n4; i4++) {
            vectorIn_GI[i-1][j-1][k-1][i4] = vector_GI[i][j][k][i4];
          }
    
    delArr4(vector_GI, nx, ny, nz);
  }

}

void NBDerivedHaloComm(int nx, int ny, int nz, double ***vector_G,
                       const VirtualTopology3D *vct, EMfields3D *EMf,
                       bool isCenterFlag, bool isFaceOnlyFlag, bool needInterp,
                       bool isParticle, bool doSwapGhostCell = false) {

  // Create a set of new pointers point to the input array vector_G. 
  double ****vector_GI;
  vector_GI = new double ***[nx];
  for (int i = 0; i < nx; i++) {
    vector_GI[i] = new double **[ny];
    for (int j = 0; j < ny; j++)
      vector_GI[i][j] = new double *[nz];
  } // i

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vector_GI[i][j][k] = &(vector_G[i][j][k]);

  NBDerivedHaloComm(nx, ny, nz, 1, vector_GI, vct, EMf, isCenterFlag,
                    isFaceOnlyFlag, needInterp, isParticle, doSwapGhostCell);

  // Delete the local pointers, but not the memory store the values. 
  for(int i = 0; i<nx; i++){
    for(int j = 0; j<ny; j++){
      delete [] vector_GI[i][j];
    }
    delete [] vector_GI[i];
  }	
  delete [] vector_GI;
}

void NBDerivedHaloComm(int nx, int ny, int nz, int n4, int n5,
                       double *****vector_GII, const VirtualTopology3D *vct,
                       EMfields3D *EMf, bool isCenterFlag, bool isFaceOnlyFlag,
                       bool needInterp, bool isParticle,
                       bool doSwapGhostCell = false) {

  // Create a set of new pointers point to the input array vector_G. 
  double ****vector_GI;
  vector_GI = new double ***[nx];
  for (int i = 0; i < nx; i++) {
    vector_GI[i] = new double **[ny];
    for (int j = 0; j < ny; j++)
      vector_GI[i][j] = new double *[nz];
  } // i

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vector_GI[i][j][k] = &(vector_GII[i][j][k][0][0]);

  bool isMassMatrix = true;
  NBDerivedHaloComm(nx, ny, nz, n4 * n5, vector_GI, vct, EMf, isCenterFlag,
                    isFaceOnlyFlag, needInterp, isParticle, doSwapGhostCell,
                    isMassMatrix);

  // Delete the local pointers, but not the memory store the values. 
  for(int i = 0; i<nx; i++){
    for(int j = 0; j<ny; j++){
      delete [] vector_GI[i][j];
    }
    delete [] vector_GI[i];
  }	
  delete [] vector_GI;

}

void communicateNodeBC(int nx, int ny, int nz, arr3_double _vector,
                       int bcFaceXrght, int bcFaceXleft, int bcFaceYrght,
                       int bcFaceYleft, int bcFaceZrght, int bcFaceZleft,
                       const VirtualTopology3D *vct, EMfields3D *EMf) {
  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false, false, false, false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft,
         bcFaceZrght, bcFaceZleft, vct);
}

void communicateNodeBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
                                 int bcFaceXrght, int bcFaceXleft,
                                 int bcFaceYrght, int bcFaceYleft,
                                 int bcFaceZrght, int bcFaceZleft,
                                 const VirtualTopology3D *vct,
                                 EMfields3D *EMf) {

  double ***vector = _vector.fetch_arr3();

  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false, true, false, false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft,
         bcFaceZrght, bcFaceZleft, vct);
}

void communicateNodeBoxStencilBC_P(int nx, int ny, int nz, arr3_double _vector,
                                   int bcFaceXrght, int bcFaceXleft,
                                   int bcFaceYrght, int bcFaceYleft,
                                   int bcFaceZrght, int bcFaceZleft,
                                   const VirtualTopology3D *vct,
                                   EMfields3D *EMf) {

  double ***vector = _vector.fetch_arr3();

  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false, true, false, true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght,
           bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateNodeBC_P(int nx, int ny, int nz, arr3_double _vector,
                         int bcFaceXrght, int bcFaceXleft, int bcFaceYrght,
                         int bcFaceYleft, int bcFaceZrght, int bcFaceZleft,
                         const VirtualTopology3D *vct, EMfields3D *EMf) {
  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false, false, false, true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght,
           bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateNode_P(int nx, int ny, int nz, double ***vector,
                       const VirtualTopology3D *vct, EMfields3D *EMf) {
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, false, false, false, true);
}

void communicateCenterBC(int nx, int ny, int nz, arr3_double _vector,
                         int bcFaceXrght, int bcFaceXleft, int bcFaceYrght,
                         int bcFaceYleft, int bcFaceZrght, int bcFaceZleft,
                         const VirtualTopology3D *vct, EMfields3D *EMf) {
  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false, false, false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft,
         bcFaceZrght, bcFaceZleft, vct);
}

void communicateSwapAddGhostCell(int nx, int ny, int nz, double ***vector_C,
				 const VirtualTopology3D *vct, EMfields3D *EMf) {
  bool isCenter = true;
  bool isFaceOnly = false;
  bool doInterp = false;
  bool isParticle = true;
  bool doSwapGhostCell = true;
  NBDerivedHaloComm(nx, ny, nz, vector_C, vct, EMf, isCenter, isFaceOnly,
                    doInterp, isParticle, doSwapGhostCell);
}

void communicateSwapAddGhostCell(int nx, int ny, int nz, int n4, 
			      double ****vector_CI,
                              const VirtualTopology3D *vct, EMfields3D *EMf) {
  bool isCenter = true;
  bool isFaceOnly = false;
  bool doInterp = false;
  bool isParticle = true;
  bool doSwapGhostCell = true;
  NBDerivedHaloComm(nx, ny, nz, n4, vector_CI, vct, EMf, isCenter, isFaceOnly,
                    doInterp, isParticle, doSwapGhostCell);
}

void communicateCenter_P(int nx, int ny, int nz, int n4, double ****vector_CI,
                           const VirtualTopology3D *vct, EMfields3D *EMf) {

  NBDerivedHaloComm(nx, ny, nz, n4, vector_CI, vct, EMf, true, false, false, true);
}



void communicateCenterBC_P(int nx, int ny, int nz, arr3_double _vector,
                           int bcFaceXrght, int bcFaceXleft, int bcFaceYrght,
                           int bcFaceYleft, int bcFaceZrght, int bcFaceZleft,
                           const VirtualTopology3D *vct, EMfields3D *EMf) {

  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false, false, true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght,
           bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateCenterBC_P(int nx, int ny, int nz, double ***vector,
                           int bcFaceXrght, int bcFaceXleft, int bcFaceYrght,
                           int bcFaceYleft, int bcFaceZrght, int bcFaceZleft,
                           const VirtualTopology3D *vct, EMfields3D *EMf) {
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false, false, true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght,
           bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}

void communicateCenterBoxStencilBC(int nx, int ny, int nz, arr3_double _vector,
                                   int bcFaceXrght, int bcFaceXleft,
                                   int bcFaceYrght, int bcFaceYleft,
                                   int bcFaceZrght, int bcFaceZleft,
                                   const VirtualTopology3D *vct,
                                   EMfields3D *EMf) {
  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, true, false, false);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght, bcFaceYleft,
         bcFaceZrght, bcFaceZleft, vct);
}

void communicateCenterBoxStencilBC_P(
    int nx, int ny, int nz, arr3_double _vector, int bcFaceXrght,
    int bcFaceXleft, int bcFaceYrght, int bcFaceYleft, int bcFaceZrght,
    int bcFaceZleft, const VirtualTopology3D *vct, EMfields3D *EMf) {

  double ***vector = _vector.fetch_arr3();
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, true, false, true);

  // ////////////////////////////////////////////////////////////////////////
  // ///////////////// APPLY the boundary conditions ////////////////////////
  // ////////////////////////////////////////////////////////////////////////
  BCface_P(nx, ny, nz, vector, bcFaceXrght, bcFaceXleft, bcFaceYrght,
           bcFaceYleft, bcFaceZrght, bcFaceZleft, vct);
}


/** communicate and sum shared ghost cells */
void communicateInterp(int nx, int ny, int nz, double ***vector,
                       const VirtualTopology3D *vct, EMfields3D *EMf) {
  NBDerivedHaloComm(nx, ny, nz, vector, vct, EMf, true, false, true, true);
}

/** communicate and sum shared ghost cells for 5D matrix*/
void communicateInterp(int nx, int ny, int nz, int n4, int n5,
                       double *****vector, const VirtualTopology3D *vct,
                       EMfields3D *EMf) {
  NBDerivedHaloComm(nx, ny, nz, n4, n5, vector, vct, EMf, true, false, true,
                    false);
}

void communicateNode_P(int nx, int ny, int nz, int n4, int n5,
                       double *****vector, const VirtualTopology3D *vct,
                       EMfields3D *EMf) {
  NBDerivedHaloComm(nx, ny, nz, n4, n5, vector, vct, EMf, false, false, false,
                    true);
}

void addBoundary(int nx, int ny, int nz, int n4, double ****vector_GI,
		 const VirtualTopology3D *vct, bool doSwapGhostCell) {
  const int iShift = (doSwapGhostCell? 1:0);

  const int nxr = nx - 2 - iShift; // nx_receive
  const int nyr = ny - 2 - iShift;
  const int nzr = nz - 2 - iShift;
  const int ir = 1 + iShift; 


  const int nxs = nx - 1; // nx_send
  const int nys = ny - 1;
  const int nzs = nz - 1;
  const int myrank = vct->getCartesian_rank();

  // Xright
  if (vct->hasXrghtNeighbor_P()) {
    for (int j = ir; j <= nyr; j++)
      for (int k = ir; k <= nzr; k++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[nxr][j][k][i4] += vector_GI[nxs][j][k][i4];
  }

  // XLEFT
  if (vct->hasXleftNeighbor_P()) {
    for (int j = ir; j <= nyr; j++)
      for (int k = ir; k <= nzr; k++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[ir][j][k][i4] += vector_GI[0][j][k][i4];
  }

  // Yright
  if (vct->hasYrghtNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int k = ir; k <= nzr; k++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[i][nyr][k][i4] += vector_GI[i][nys][k][i4];
  }

  // Yleft
  if (vct->hasYleftNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int k = ir; k <= nzr; k++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[i][ir][k][i4] += vector_GI[i][0][k][i4];
  }

  // Zright
  if (vct->hasZrghtNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int j = ir; j <= nyr; j++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[i][j][nzr][i4] += vector_GI[i][j][nzs][i4];
  }

  // ZLEFT
  if (vct->hasZleftNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int j = ir; j <= nyr; j++)
        for (int i4 = 0; i4 < n4; i4++)
          vector_GI[i][j][ir][i4] += vector_GI[i][j][0][i4];
  }

  // Add edge z
  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P()) {
    for (int k = ir; k <= nzr; k++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[nxr][nyr][k][i4] += vector_GI[nxs][nys][k][i4];
  }

  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P()) {
    for (int k = ir; k <= nzr; k++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[ir][ir][k][i4] += vector_GI[0][0][k][i4];
  }

  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P()) {
    for (int i = ir; i <= nzr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[nxr][ir][i][i4] += vector_GI[nxs][0][i][i4];
  }

  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P()) {
    for (int i = ir; i <= nzr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[ir][nyr][i][i4] += vector_GI[0][nys][i][i4];
  }

  // Add edge y
  if (vct->hasXrghtNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = ir; i <= nyr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[nxr][i][nzr][i4] += vector_GI[nxs][i][nzs][i4];
  }

  if (vct->hasXleftNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = ir; i <= nyr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[ir][i][ir][i4] += vector_GI[0][i][0][i4];
  }

  if (vct->hasXleftNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = ir; i <= nyr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[ir][i][nzr][i4] += vector_GI[0][i][nzs][i4];
  }

  if (vct->hasXrghtNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = ir; i <= nyr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[nxr][i][ir][i4] += vector_GI[nxs][i][0][i4];
  }

  // Add edge x
  if (vct->hasYrghtNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[i][nyr][nzr][i4] += vector_GI[i][nys][nzs][i4];
  }

  if (vct->hasYleftNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[i][ir][ir][i4] += vector_GI[i][0][0][i4];
  }

  if (vct->hasYleftNeighbor_P() && vct->hasZrghtNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[i][ir][nzr][i4] += vector_GI[i][0][nzs][i4];
  }

  if (vct->hasYrghtNeighbor_P() && vct->hasZleftNeighbor_P()) {
    for (int i = ir; i <= nxr; i++)
      for (int i4 = 0; i4 < n4; i4++)
        vector_GI[i][nyr][ir][i4] += vector_GI[i][nys][0][i4];
  }

  // Add corner.
  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P() &&
      vct->hasZrghtNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[nxr][nyr][nzr][i4] += vector_GI[nxs][nys][nzs][i4];

  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P() &&
      vct->hasZrghtNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[ir][nyr][nzr][i4] += vector_GI[0][nys][nzs][i4];

  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P() &&
      vct->hasZrghtNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[nxr][ir][nzr][i4] += vector_GI[nxs][0][nzs][i4];

  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P() &&
      vct->hasZrghtNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[ir][ir][nzr][i4] += vector_GI[0][0][nzs][i4];

  if (vct->hasXrghtNeighbor_P() && vct->hasYrghtNeighbor_P() &&
      vct->hasZleftNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[nxr][nyr][ir][i4] += vector_GI[nxs][nys][0][i4];

  if (vct->hasXleftNeighbor_P() && vct->hasYrghtNeighbor_P() &&
      vct->hasZleftNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[ir][nyr][ir][i4] += vector_GI[0][nys][0][i4];

  if (vct->hasXrghtNeighbor_P() && vct->hasYleftNeighbor_P() &&
      vct->hasZleftNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[nxr][ir][ir][i4] += vector_GI[nxs][0][0][i4];

  if (vct->hasXleftNeighbor_P() && vct->hasYleftNeighbor_P() &&
      vct->hasZleftNeighbor_P())
    for (int i4 = 0; i4 < n4; i4++)
      vector_GI[ir][ir][ir][i4] += vector_GI[0][0][0][i4];
}
