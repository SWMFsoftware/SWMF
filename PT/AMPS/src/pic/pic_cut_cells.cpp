//$Id$
//the set of function for processing of the cut-cells and cut-faces of the internal irregular surfaces

#include "pic.h"

int PIC::Mesh::IrregularSurface::nCutFaceInformationCopyAttempts=0; 

void PIC::Mesh::IrregularSurface::InitExternalNormalVector() {
  //calculate external normals to the faces
  double xStart[3],xFinish[3],l,l0,*norm;
  long int nface;
  cTriangleFace *fcptr;
  int idim,iIntersections;

  const double angleCosMin=cos(85.0/180.0*3.141592654);

  int nStartFace,nFinishFace,nTotalThreads,ThisThread,nFaceThread;

  //check whether external normal vectors are already have been determined for the surface trianguletion
  unsigned long int TriangulationSignature;
  char fname[_MAX_STRING_LENGTH_PIC_];
  FILE *fExternalVectorFile=NULL;

  TriangulationSignature=CutCell::GetTriangulationSignature();
  sprintf(fname,"amr.sig=0x%lx.TriangulationExternalNormalVector.bin",TriangulationSignature);

  fExternalVectorFile=fopen(fname,"r");

  if (fExternalVectorFile!=NULL) {
    //the binary file containing the external normal information exists -> load it
    if (PIC::ThisThread==0) printf("$PREFIX: Binary file with the external normals (%s) is found: loading.... ",fname);

    for (nface=0;nface<nBoundaryTriangleFaces;nface++) fread(BoundaryTriangleFaces[nface].ExternalNormal,sizeof(double),3,fExternalVectorFile);

    fclose(fExternalVectorFile);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) printf("done\n");

    return;
  }

  if (PIC::ThisThread==0) printf("$PREFIX: Binary file with the external normals is not found: generating.... "); 

  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

  nFaceThread=nBoundaryTriangleFaces/nTotalThreads;
  nStartFace=nFaceThread*ThisThread;
  nFinishFace=nStartFace+nFaceThread;
  if (ThisThread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none) shared(nStartFace,nFinishFace,BoundaryTriangleFaces,lDomain) \
  private (nface,fcptr,norm,idim,xFinish,xStart,l,l0,iIntersections)
#endif
  for (nface=nStartFace;nface<nFinishFace;nface++) {
    fcptr=BoundaryTriangleFaces+nface;

    //get the center point and the normal of the face
    fcptr->GetCenterPosition(xStart);
    norm=fcptr->ExternalNormal;

    do {
      for (idim=0,l=0.0,l0=0.0;idim<3;idim++) {
        xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
        l+=pow(xFinish[idim],2);
        l0+=xFinish[idim]*norm[idim];
      }

      l=sqrt(l);
    }
    while (fabs(l0)/l<angleCosMin);

    if (l0<0.0) l*=-1.0;
    for (idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+xStart[idim];

    //count face intersections
    iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(xStart,xFinish,fcptr->MeshFileID,fcptr);

    if (iIntersections%2!=0) {
      //the norm has to be reversed
      for (idim=0;idim<3;idim++) fcptr->ExternalNormal[idim]*=-1.0;
    }
  }

  //collect the surface normals
  double *sendBuffer=new double[3*2*nFaceThread];
  int thread,cnt;

  for (thread=0;thread<nTotalThreads;thread++) {
    nStartFace=nFaceThread*thread;
    nFinishFace=nStartFace+nFaceThread;
    if (thread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

    if (thread==ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(sendBuffer+3*cnt,BoundaryTriangleFaces[nface].ExternalNormal,3*sizeof(double));
    }

    MPI_Bcast(sendBuffer,3*(nFinishFace-nStartFace),MPI_DOUBLE,thread,MPI_GLOBAL_COMMUNICATOR);

    if (thread!=ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(BoundaryTriangleFaces[nface].ExternalNormal,sendBuffer+3*cnt,3*sizeof(double));
    }
  }

  delete [] sendBuffer;

  //save a file with the external normals for the future use
  if (PIC::ThisThread==0) {
    fExternalVectorFile=fopen(fname,"w");

    for (nface=0;nface<nBoundaryTriangleFaces;nface++) fwrite(BoundaryTriangleFaces[nface].ExternalNormal,sizeof(double),3,fExternalVectorFile);

    fclose(fExternalVectorFile);
  }


  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) printf("done\n");
}

//check weather a point (x0) in insed the domain:
//if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
//l -> is a random ray (intersection search) direction
bool PIC::Mesh::IrregularSurface::CheckPointInsideDomain_default(double *x,PIC::Mesh::IrregularSurface::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS) {
  int nface,nfaceStart,nfaceFinish,iIntersections;
  double l,xFinish[3];
  int idim;
  bool flag=true;

  if (SurfaceTriangulation==NULL) return true;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

  //distribute ditrction of the search
  for (l=0.0,idim=0;idim<3;idim++) {
    xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
    l+=pow(xFinish[idim],2);
  }

  for (l=sqrt(l),idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+x[idim];

  //xFinish is outside of the domain -> the point outside of the surface
  //calculate the number of the face intersections between 'x' and 'xFinish'

  iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(x,xFinish,-1);

  return (iIntersections%2==0) ? true : false;

 /*   SearchDirection[idim]/=l;


  static bool initflag=false;
  static int ThisThread,nTotalThreads;

  if (initflag==false) {
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);
    initflag=true;
  }

  if (ParallelCheck==true) {
    nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
    nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
    if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;
  }
  else nfaceStart=0,nfaceFinish=nSurfaceTriangulation;

  bool flagbuffer[nTotalThreads];

  do {
    //distribute ditrction of the search
    for (l=0.0,idim=0;idim<3;idim++) {
      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
    }

    if (ParallelCheck==true) MPI_Bcast(SearchDirection,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    for (l=sqrt(l),idim=0;idim<3;idim++) SearchDirection[idim]/=l;
    iIntersections=0;
    flag=true;

    //find intersections with the faces on the mesh
    for (nface=nfaceStart;nface<nfaceFinish;nface++) {
      if (SurfaceTriangulation[nface].RayIntersection(x,SearchDirection,EPS)==true) iIntersections++;

      for (l=0.0,idim=0;idim<3;idim++) l+=pow(SurfaceTriangulation[nface].ExternalNormal[idim]*SearchDirection[idim],2);
      if (l<1.0E-10) {
        flag=false;
        break;
      }

    }

    if (ParallelCheck==true) {
      MPI_Gather(&flag,sizeof(bool),MPI_CHAR,flagbuffer,sizeof(bool),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
      if (ThisThread==0) for (int thread=1;thread<nTotalThreads;thread++) if (flagbuffer[thread]==false) flag=false;
      MPI_Bcast(&flag,sizeof(bool),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
    }
  }
  while (flag==false);

  if (ParallelCheck==true) {
    int t;

    MPI_Allreduce(&iIntersections,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    iIntersections=t;
  }


  return (2*(iIntersections/2)==iIntersections) ? true : false;
*/}


//copy information about the triangular cut faces into the nighbouring blocks 
void PIC::Mesh::IrregularSurface::CopyCutFaceInformation(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  static unsigned int *BoundaryTriangleMap=NULL;
  static const int nResetBoundaryTriangleMapTestCounter=100000000;
  static unsigned int BoundaryTriangleMapTestCounter=0;

  static int CutFaceDescriptorTablePointer=-1;
  static CutCell::cTriangleFaceDescriptor *CutFaceDescriptorTable=NULL;
  static int CutFaceDescriptorTableLength=(int)(0.5*CutCell::nBoundaryTriangleFaces);     

  struct CutFaceData {
    void ResetTempPointerTable(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *nd) {
      if (nd->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
        int i;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

        for (i=0;i<(1<<DIM);i++) if ((downNode=nd->downNode[i])!=NULL) {
          ResetTempPointerTable(downNode);
        }
      }
      else nd->neibFirstTriangleCutFace_temp=NULL;
    } 

    void SetTempCutFacePointers(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *nd) {
      if (nd->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
        int i;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

        for (i=0;i<(1<<DIM);i++) if ((downNode=nd->downNode[i])!=NULL) {
          SetTempCutFacePointers(downNode);
        }
      }
      else if (nd->neibFirstTriangleCutFace_temp!=NULL) {
        //add the pointers to the neibFirstTriangleCutFace list
        CutCell::cTriangleFaceDescriptor *t;

        if (nd->neibFirstTriangleCutFace==NULL) {
          nd->neibFirstTriangleCutFace=nd->neibFirstTriangleCutFace_temp;
          nd->neibFirstTriangleCutFace_temp=NULL;
        }
        else {
          t=nd->neibFirstTriangleCutFace;
          while (t->next!=NULL) t=t->next; 

          //link the lists
          t->next=nd->neibFirstTriangleCutFace_temp;
          nd->neibFirstTriangleCutFace_temp->prev=t;

          nd->neibFirstTriangleCutFace_temp=NULL;
        }
      }
    }

    //add cut faces from block NodeFrom to NodeTo 
    void AddCutFaceData(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *To,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *From) { 
      int iTestNeib,iTestNode;
      CutCell::cTriangleFaceDescriptor *t;

      //save all cut-faces from the neib node
      if (++BoundaryTriangleMapTestCounter==nResetBoundaryTriangleMapTestCounter) {
        //reset the counters
        BoundaryTriangleMapTestCounter=1;
        for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) BoundaryTriangleMap[i]=0;
      }

      for (iTestNode=0;iTestNode<3;iTestNode++) {
        switch (iTestNode) {
        case 0: 
          t=To->FirstTriangleCutFace;
          break;
        case 1: 
          t=To->neibFirstTriangleCutFace;
          break;
        case 2:
          t=To->neibFirstTriangleCutFace_temp;
          break;
        default:
          exit(__LINE__,__FILE__,"Error: something is wrong");
        }

        for (;t!=NULL;t=t->next) BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
      }

      //check whether all cut-faces from startNode are regiteres in neibNode
      for (iTestNode=0;iTestNode<2;iTestNode++) {
        t=(iTestNode==0) ? From->FirstTriangleCutFace : From->neibFirstTriangleCutFace;

        for (;t!=NULL;t=t->next) if (BoundaryTriangleMap[t->TriangleFace->Temp_ID]!=BoundaryTriangleMapTestCounter) {
          CutCell::cTriangleFaceDescriptor *NewDescriptor;

          if (CutFaceDescriptorTablePointer==-1) {
            CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
            CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;
          }

          NewDescriptor=CutFaceDescriptorTable+CutFaceDescriptorTablePointer;
          CutFaceDescriptorTablePointer--;

          NewDescriptor->TriangleFace=t->TriangleFace;
          NewDescriptor->prev=NULL;
          NewDescriptor->next=To->neibFirstTriangleCutFace_temp;

          if (To->neibFirstTriangleCutFace_temp!=NULL) To->neibFirstTriangleCutFace_temp->prev=NewDescriptor;
          To->neibFirstTriangleCutFace_temp=NewDescriptor;

          BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
        }
      }
    }

 
  } ProcessCutFaceData;
     

  if (startNode==PIC::Mesh::mesh.rootTree) {
    BoundaryTriangleMapTestCounter=1;
    BoundaryTriangleMap=new unsigned int [CutCell::nBoundaryTriangleFaces];  
         
    for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
      BoundaryTriangleMap[i]=0,CutCell::BoundaryTriangleFaces[i].Temp_ID=i;
    }

    //allocate the buffer with the new descriptors if needed 
    if (CutFaceDescriptorTable==NULL) {
      CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
      CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;
    }


    ProcessCutFaceData.ResetTempPointerTable(PIC::Mesh::mesh.rootTree); 
  }

  //move to the botton of the tree
  if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;
    double c;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      CopyCutFaceInformation(downNode); 
    }
  }
  else if ((startNode->FirstTriangleCutFace!=NULL)||(startNode->neibFirstTriangleCutFace!=NULL)) {
    const int ProcessedNeibBlockTableLength=60;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ProcessedNeibBlockTable[ProcessedNeibBlockTableLength];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
    int i,j,k,iFace,iEdge,iCorner,iNeib,ProcessedNodeCounter=0;
    bool found;

    //connection through corners
    for (iCorner=0;iCorner<(1<<DIM);iCorner++) if ((neibNode=startNode->GetNeibCorner(iCorner))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information 
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);  
      }
    }
  
    //connection through edges 
    for (iEdge=0;iEdge<12;iEdge++) for (i=0;i<2;i++) if ((neibNode=startNode->GetNeibEdge(iEdge,i))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);
      }
    }

    //connection through faces   
    for (iFace=0;iFace<6;iFace++) for (i=0;i<2;i++) for (j=0;j<2;j++) if ((neibNode=startNode->GetNeibFace(iFace,i,j))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);
      }
    }

    if (ProcessedNodeCounter>ProcessedNeibBlockTableLength) exit(__LINE__,__FILE__,"Error: the counting is wrong"); 
  } 

 





/*
    int iNeib,jNeib,kNeib;
    CutCell::cTriangleFaceDescriptor *t,*tNeib;

    //scan through all neibours
    int iNeibNode,jNeibNode,kNeibNode;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;

    static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=-1,kNeibNodeMax=1;

    for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) { 
      for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) {
        for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  {
          if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) {
            if ((neibNode=PIC::Mesh::mesh.getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
              int iTestNeib,iTestNode;

              //save all cut-faces from the neib node
              if (++BoundaryTriangleMapTestCounter==nResetBoundaryTriangleMapTestCounter) {
                //reset the counters 
                BoundaryTriangleMapTestCounter=1;
                for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) BoundaryTriangleMap[i]=0; 
              }

              for (iTestNode=0;iTestNode<2;iTestNode++) {
                t=(iTestNode==0) ? neibNode->FirstTriangleCutFace : neibNode->neibFirstTriangleCutFace;

                for (;t!=NULL;t=t->next) BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
              }     

              //check whether all cut-faces from startNode are regiteres in neibNode
              for (iTestNode=0;iTestNode<2;iTestNode++) {
                t=(iTestNode==0) ? startNode->FirstTriangleCutFace : startNode->neibFirstTriangleCutFace; 

                for (;t!=NULL;t=t->next) if (BoundaryTriangleMap[t->TriangleFace->Temp_ID]!=BoundaryTriangleMapTestCounter) {
                  CutCell::cTriangleFaceDescriptor *NewDescriptor;

                  if (CutFaceDescriptorTablePointer==-1) {
                    CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
                    CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;  
                  }

                  NewDescriptor=CutFaceDescriptorTable+CutFaceDescriptorTablePointer;
                  CutFaceDescriptorTablePointer--;

                  NewDescriptor->TriangleFace=t->TriangleFace;
                  NewDescriptor->prev=NULL;
                  NewDescriptor->next=neibNode->neibFirstTriangleCutFace;

                  if (neibNode->neibFirstTriangleCutFace!=NULL) neibNode->neibFirstTriangleCutFace->prev=NewDescriptor;
                  neibNode->neibFirstTriangleCutFace=NewDescriptor;   

                  BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter; 
                }
              }

            }
          }
        }
      }
    } 
  }
*/

  if (startNode==PIC::Mesh::mesh.rootTree) {
    delete [] BoundaryTriangleMap;

    BoundaryTriangleMapTestCounter=0;
    BoundaryTriangleMap=NULL;

    ProcessCutFaceData.SetTempCutFacePointers(PIC::Mesh::mesh.rootTree);
  }
}
 





