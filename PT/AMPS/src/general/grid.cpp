//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//==================================================
//$Id$
//==================================================

#define grid_cpp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <algorithm>

#include "dsmc.h"
#include "cell.h"
#include "face.h"
#include "node.h"
#include "array_1d.h"

extern int DIM;
//==================================================
void cutstr(char* dest, char* src)
{
  int i,p;

  for(p=0;(src[p]==' ')&&(src[p]!='\0')&&(src[p]!='\n');p++);
  for(i=0;(src[i+p]!='\0')&&(src[i+p]!='\n');i++) src[i]=src[i+p];
  src[i]='\0';

  p=0;
  while ((src[p]!='\0')&&(src[p]!='\n')&&(src[p]!=' ')) {
    dest[p]=src[p];
    p++;
    }
  dest[p]='\0';

  for(;(src[p]!='\0')&&(src[p]!='\n')&&(src[p]==' ');p++);
  
  for (i=0;(src[i+p]!='\0')&&(src[i+p]!='\n');i++) src[i]=src[i+p];
  src[i]='\0'; 
}

//==================================================
void error(long int line) {
  printf("Error in reading of grid's file (line=%i)\n",line);
  exit(__LINE__,__FILE__);
}

//==================================================
void Cgrid::OutputInFile(char *fname) {
  FILE* fout;
  int idim;
  long int n,f,c;

  fout=fopen(fname,"w");

  fprintf(fout," < Nnodes >< Nfaces >< Ncells >\n");
  fprintf(fout,"  %i  %i  %i\n",nnodes,nfaces,ncells);

  fprintf(fout,"node's coordinates\n");
  for (n=0;n<nnodes;n++) fprintf(fout,node[n]);

  fprintf(fout,"faces: nodes, faceat\n");
  for (f=0;f<nfaces;f++) {
    fprintf(fout,"face=%i   nodes: ",f);
    for (idim=0;idim<DIM+1;idim++) {
      for (n=0;n<nnodes;n++)  
        if (&(node[n])==face[f].node[idim]) fprintf(fout,"%i  ",n); 
    } 
    fprintf(fout,"faceat=%i\n",(int)face[f].faceat);
  }

  fprintf(fout,"cells: nodes,faces,neighbours\n");
  for (c=0;c<ncells;c++) {
    fprintf(fout,"cell=%i nodes: ",c);
    for (idim=0;idim<DIM+1;idim++) 
      for (n=0;n<nnodes;n++)
        if (&(node[n])==cell[c].node[idim]) fprintf(fout,"%i  ",n);
 
    fprintf(fout,"faces:  ");
    for (idim=0;idim<DIM+1;idim++) 
      for (f=0;f<nfaces;f++)
        if (&(face[f])==cell[c].face[idim]) fprintf(fout,"%i  ",f);
  
    fprintf(fout,"neighbours:  ");
    for (idim=0;idim<DIM+1;idim++) 
       fprintf(fout,"%i  ",(long int)cell[c].neighbour_cellno[idim]);

    fprintf(fout,"\n"); 
  } 


  long int b; 
  int e1,e2;
  fprintf(fout","\n");
  fprintf(fout,"============ connection data ===========\n");
  for (c=0;c<ncells;c++)   
    for (f=0;f<DIM+1;f++) {
      fprintf(fout,"\ncell=%i, face=%i\n",c,f);
 
      for (b=0;b<DIM+1;b++) { 
        if (b==f) continue;
        fprintf(fout,"basis=%i  => basis=%i\n",b,cells_connection_data[c].nbasis[b][f]);
        fprintf(fout,"ne=%i\n",cells_connection_data[c].ne[b][f]); 
      
        for (idim=0;idim<DIM-1;idim++) {
          e1=cells_connection_data[c].ne1[idim][b][f]; 
          e2=cells_connection_data[c].ne2[idim][b][f];
          fprintf(fout,"e1=%i  =>  e2=%i\n",e1,e2);
        }
      }
    }
  fclose(fout);
}

 
  
//==================================================
void Cgrid::LoadFromFile(char* file_name)
{
  FILE* fd;
  char str1[200],str[200];
  long int i;
  int idim;
  float x[3];
  long int line;

  if (DIM==0) {
    cell=new Ccell[1];
    ncells=1; 
    return;
  }

  fd=fopen(file_name,"r");
  fgets(str,200,fd);

  fgets(str,200,fd);
  cutstr(str1,str);
  nnodes=atoi(str1);
  cutstr(str1,str);
  nfaces=atoi(str1);
  cutstr(str1,str);
  ncells=atoi(str1);  

  line=2;
  if (errno!=0) error(line);

  node=new Cnode[nnodes];
  face=new Cface[nfaces];
  cell=new Ccell[ncells];

  fgets(str,200,fd);line++;
  for (i=0;i<nnodes;i++) {
    fgets(str,200,fd);
    for (idim=0;idim<DIM;idim++) {
      cutstr(str1,str);
      x[idim]=atof(str1);
    } 
    node[i].SetX(x);
    line++;if (errno!=0) error(line);
  }

  fgets(str,200,fd);line++;
  for (i=0;i<nfaces;i++) {
    fgets(str,200,fd);
    for (idim=0;idim<DIM;idim++) {
      cutstr(str1,str);
      face[i].nodeno[idim]=atoi(str1)-1;
      face[i].node[idim]=&node[face[i].nodeno[idim]];
    }
    cutstr(str1,str);
    face[i].faceat=atoi(str1);
    line++;if (errno!=0) error(line);
  }

  fgets(str,200,fd);line++;
  for (i=0;i<ncells;i++) {
    fgets(str,200,fd);
    for (idim=0;idim<DIM+1;idim++) {
      cutstr(str1,str);
      cell[i].nodeno[idim]=atoi(str1)-1;
      cell[i].node[idim]=&node[cell[i].nodeno[idim]];
    }
    
    for (idim=0;idim<DIM+1;idim++) {
      cutstr(str1,str);
      cell[i].faceno[idim]=atoi(str1)-1;
      cell[i].face[idim]=&face[cell[i].faceno[idim]];
    }
     
    for (idim=0;idim<DIM+1;idim++) {
      cutstr(str1,str);
      cell[i].neighbour_cellno[idim]=atoi(str1)-1;
    }
    line++;if (errno!=0) error(line);
  }
}

//==================================================
double Cgrid::Measure()
{
  double measure;
  long int ncell;

  measure=0.0;
  for (ncell=0;ncell<ncells;ncell++) 
    measure+=cell[ncell].Measure();

  return measure;
}

//==================================================
void Cgrid::InitGridData() {

  if (DIM==0) {
    cell[0].first_ptr=-1;
    cell[0].thread=0;
    cell[0].sbdm=0;
    return;
  }

  for (long int ncell=0;ncell<ncells;ncell++) {
    cell[ncell].first_ptr=-1;
    cell[ncell].thread=0;
    cell[ncell].sbdm=0;
    cell[ncell].InitCellNormal();
  }

  switch(DIM) {
  case 3 :
    GetTMatrix3D();
    break;
  case 2 :
    GetTMatrix2D();
    break;
  case 1 :
    GetTMatrix1D();
    break;
  default :
    printf("Error: Cgrid::InitGridData()\n");
    printf("value of variable DIM is out range (DIM=%i)\n",DIM);
    exit(__LINE__,__FILE__);
  }  

  InitInterpolationData();
  InitCellConnectionData();
}

//==================================================
void Cgrid::InitInterpolationData() {
  long int k,ncell,nnode;
  long int* nnn=new long int[nnodes];
  float* sum=new float[nnodes];
  float measure;
  int idim;

  for (nnode=0;nnode<nnodes;nnode++) {
    nnn[nnode]=0;
    sum[nnode]=0.0;
  } 

  for (ncell=0;ncell<ncells;ncell++) {
    measure=cell[ncell].Measure();
    for (idim=0;idim<DIM+1;idim++) {
      k=cell[ncell].nodeno[idim];
      sum[k]+=measure;
      nnn[k]+=1;
    }
  }

  for (nnode=0;nnode<nnodes;nnode++) {    
    node[nnode].InterpolationWeight=new float[nnn[nnode]+1];
    node[nnode].InterpolationMask=new long int [nnn[nnode]+1];
    nnn[nnode]=0;
  } 

  for (ncell=0;ncell<ncells;ncell++) {
    measure=cell[ncell].Measure();
    for (idim=0;idim<DIM+1;idim++) {
      k=cell[ncell].nodeno[idim];
      node[k].InterpolationMask[nnn[k]]=ncell;
      node[k].InterpolationWeight[nnn[k]]=measure/sum[k];
      nnn[k]++;
    }
  }

  for (nnode=0;nnode<nnodes;nnode++) 
    node[nnode].InterpolationMask[nnn[nnode]]=-1;
} 
       
//==================================================
//Get "global" number of cell, which contains point x 
long int Cgrid::GetNCell(float* x) {
  int idim,i;
  long int nnode,ncell;
  float xmin[3],xmax[3],locx[3],summ;
  array_1d<float> x_node(DIM);
  bool flag;

  for (ncell=0;ncell<ncells;ncell++) {
    x_node=cell[ncell].node[0]->X();
    for (idim=0;idim<DIM;idim++) {
      xmax[idim]=x_node(idim);
      xmin[idim]=x_node(idim);
    }
   
    for(nnode=1;nnode<DIM+1;nnode++) {
      x_node=cell[ncell].node[nnode]->X();
      for (idim=0;idim<DIM;idim++) {
        if (xmax[idim]<x_node(idim)) xmax[idim]=x_node(idim);
        if (xmin[idim]>x_node(idim)) xmin[idim]=x_node(idim);
      }
    }

    flag=true;
    for (idim=0;idim<DIM;idim++) 
      if ((x[idim]<xmin[idim])||(x[idim]>xmax[idim])) flag=false;

    if (flag==false) continue;

    for (i=0;i<DIM;i++) locx[i]=0.0; 
    x_node=cell[ncell].node[0]->X();

    for (idim=0;idim<DIM;idim++)
      for (i=0;i<DIM;i++)  
        locx[idim]+=TMatrix(idim,i,0,ncell)*(x[i]-x_node(i)); 

    summ=0.0;
    flag=true;
    for (idim=0;idim<DIM;idim++) {
      summ+=locx[idim];
      if ((locx[idim]<-1.0E-5)||(locx[idim]>1.00001)) flag=false;
    }
    if ((summ<-1.0E-5)||(summ>1.00001)) flag=false;

    if (flag==true) break;
  }

  if (flag==false) {
    printf("Error: proc. Cgrid::GetNCell\n");
    printf("Cannot find cell which contain point x=");
    for (idim=0;idim<DIM;idim++) printf("  %e",x[idim]);
    printf("\n");
    exit(__LINE__,__FILE__);
  }  

  return ncell;     
}  

//==================================================
void Cgrid::GetTMatrix2D() {
  array_1d<float> a0(DIM),a1(DIM),b(DIM);
  int i,n[3];
  double detA,det;

  TMatrix.init(DIM,DIM,DIM+1,ncells);

  for(long int ncell=0;ncell<ncells;ncell++) 
    for (int nbasis=0;nbasis<DIM+1;nbasis++) {  
      for (i=0;i<DIM;i++) {
        n[i]=nbasis+i+1;
        if (n[i]>=DIM+1) n[i]-=(DIM+1);
      }

    a0=cell[ncell].node[n[0]]->X()-cell[ncell].node[nbasis]->X();
    a1=cell[ncell].node[n[1]]->X()-cell[ncell].node[nbasis]->X();

    detA=a0(0)*a1(1)-a1(0)*a0(1);
    
    for (i=0;i<DIM;i++) {
      b=0.0;
      b(i)=1.0;

      det=b(0)*a1(1)-b(1)*a1(0);
      TMatrix(0,i,nbasis,ncell)=det/detA;

      det=a0(0)*b(1)-b(0)*a0(1);
      TMatrix(1,i,nbasis,ncell)=det/detA;
    }
  }
}

//==================================================
void Cgrid::InitCellConnectionData() {
  long int ncell,nbr_ncell; 
  int i,nface,nbasis,nbr_nbasis;

  cells_connection_data=new cells_connection_data_type[ncells]; 
  bvector.init(ncells,DIM,DIM+1,DIM+1);

  for (ncell=0;ncell<ncells;ncell++) {
    for (nface=0;nface<=DIM;nface++) {
      nbr_ncell=cell[ncell].neighbour_cellno[nface]; 

      for (nbasis=0;nbasis<DIM+1;nbasis++) {
        if (nbasis==nface) continue;

        if (nbr_ncell>=0) {
          for (nbr_nbasis=0;nbr_nbasis<DIM+1;nbr_nbasis++) 
            if (cell[ncell].node[nbasis]==cell[nbr_ncell].node[nbr_nbasis]) {
              cells_connection_data[ncell].nbasis[nbasis][nface]=nbr_nbasis;
              break;
            }
        }

        char ne,ne1,ne_node;
        for (ne=0;ne<DIM;ne++) {
          ne_node=nbasis+ne+1;
          if (ne_node>DIM) ne_node-=DIM+1;  
          if (ne_node==nface) {
            cells_connection_data[ncell].ne[nbasis][nface]=ne;
            break;
          }
        }
      
        ne1=0;
        for (ne=0;ne<DIM;ne++) {
          if (ne==cells_connection_data[ncell].ne[nbasis][nface]) continue;

          cells_connection_data[ncell].ne1[ne1][nbasis][nface]=ne;
          if (nbr_ncell>=0) {
            i=nbasis+ne+1; 
            if (i>DIM) i-=DIM+1;
    
            char nbr_nface,nbr_ne,nbr_i;
            for (nbr_ne=0;nbr_ne<DIM;nbr_ne++) {
              nbr_i=cells_connection_data[ncell].nbasis[nbasis][nface]+1+nbr_ne;
              if (nbr_i>DIM) nbr_i-=DIM+1;
              if (cell[ncell].node[i]==cell[nbr_ncell].node[nbr_i]) {
                cells_connection_data[ncell].ne2[ne1][nbasis][nface]=nbr_ne;
                break;
              }
            }
          }
          ne1++;
        }

//get bvector
        if (nbr_ncell>=0) {
          int idim;
          array_1d<float> e(DIM);
          double a;
          int n,n_common_nodes;

          n=cells_connection_data[ncell].ne[nbasis][nface]+nbasis+1; 
          if (n>DIM) n-=DIM+1;
          e=cell[ncell].node[n]->X()-cell[ncell].node[nbasis]->X();

          for (idim=0;idim<DIM;idim++) {
            a=0.0;
            for (i=0;i<DIM;i++) a+=e(i)*TMatrix(idim,i,nbr_nbasis,nbr_ncell); 
            bvector(ncell,idim,nbasis,nface)=a; 
          }

          for (char nbr_nface=0;nbr_nface<DIM+1;nbr_nface++) {
            n_common_nodes=0;
            for (char nnode=0;nnode<DIM;nnode++) 
              for (char nbr_nnode=0;nbr_nnode<DIM;nbr_nnode++) 
                if (cell[ncell].face[nface]->node[nnode]==
                  cell[nbr_ncell].face[nbr_nface]->node[nbr_nnode])
                    n_common_nodes++;
              
            if (n_common_nodes==DIM) {
              cells_connection_data[ncell].pface[nbasis][nface]=nbr_nface;
              break;  
            }
          }
        }
      }
    }
  }
}

//==================================================
void Cgrid::ChangeLocalVector2D(float* a,char nb1,char nb2) {
  int idim,dn;
  float b[2];

  for (idim=0;idim<DIM;idim++) b[idim]=0.0;
  dn=nb2-nb1;
  if (dn<0) dn+=3;

  switch(dn) {
  case 1:
    b[0]=a[1];
    b[1]=-(a[0]+a[1]);
    break;
  case 2:
    b[0]=-(a[0]+a[1]);
    b[1]=a[0];
    break;
  default :
    printf("Error: proc. grid.cpp::Cdsmc::ChangeLocalVector2D\n");
    exit(__LINE__,__FILE__);
  }

  for (idim=0;idim<DIM;idim++) a[idim]=b[idim];
}
//==================================================
void Cgrid::ChangeLocalPositionVector2D(float* a,char nb1, char nb2) {
  float b[2];
  int idim,dn;

  for (idim=0;idim<DIM;idim++) b[idim]=0.0;
  dn=nb2-nb1;
  if (dn<0) dn+=3;
  
  switch(dn) {
  case 1:
    b[0]=a[1];
    b[1]=1.0-(a[0]+a[1]);
    break;
  case 2:
    b[0]=1.0-(a[0]+a[1]);
    b[1]=a[0];
    break;
  default :
    printf("Error: proc. grid.cpp::Cdsmc::ChangeLocalPositionVector2D\n");
    exit(__LINE__,__FILE__);
  }  

  for (idim=0;idim<DIM;idim++) a[idim]=b[idim];
}

//==================================================
void Cgrid::GetTMatrix1D() {
  int i,n;
  double l;

  TMatrix.init(DIM,DIM,DIM+1,ncells);

  for(long int ncell=0;ncell<ncells;ncell++) 
    for (int nbasis=0;nbasis<DIM+1;nbasis++) {
      n=nbasis+1;
      if (n>=DIM+1) n-=(DIM+1);

      l=(cell[ncell].node[n]->X()-cell[ncell].node[nbasis]->X())(0);

      TMatrix(0,0,nbasis,ncell)=1.0/l;
    }
}

//==================================================
void Cgrid::ChangeLocalVector1D(float* a,char nb1,char nb2) {
  int idim,dn;

  dn=nb2-nb1;
  if ((dn!=1)&&(dn!=-1)) {
    printf("Error: proc. Cgrid::ChangeLocalVector1D\n");
    printf("dn=%i\n",dn);
    exit(__LINE__,__FILE__);
  }

  a[0]=-a[0];
}

//==================================================
void Cgrid::ChangeLocalPositionVector1D(float* a,char nb1,char nb2) {
  int idim,dn;

  dn=nb2-nb1;
  if ((dn!=1)&&(dn!=-1)) {
    printf("Error: proc. Cgrid::ChangeLocalPositionVector1D\n");
    printf("dn=%i\n",dn);
    exit(__LINE__,__FILE__);
  }

  a[0]=1.0-a[0];
}

//==================================================
void Cgrid::SaveImageFile(int fd) {
  long int ncell;

  for (ncell=0;ncell<ncells;ncell++) {
    write(fd,&cell[ncell].first_ptr,sizeof(cell[ncell].first_ptr));
    write(fd,&cell[ncell].thread,sizeof(cell[ncell].thread));
    write(fd,&cell[ncell].sbdm,sizeof(cell[ncell].sbdm));
  }
}
  
//==================================================
void Cgrid::LoadImageFile(int fd) {
  long int ncell;

  for (ncell=0;ncell<ncells;ncell++) {
    read(fd,&cell[ncell].first_ptr,sizeof(cell[ncell].first_ptr));
    read(fd,&cell[ncell].thread,sizeof(cell[ncell].thread));
    read(fd,&cell[ncell].sbdm,sizeof(cell[ncell].sbdm));
  }
}

//==================================================
void Cgrid::GetTMatrix3D() {
  array_1d<float> a0(DIM),a1(DIM),a2(DIM),b(DIM);
  int i,n[3];
  double detA,det;

  TMatrix.init(DIM,DIM,DIM+1,ncells);

  for(long int ncell=0;ncell<ncells;ncell++)
    for (int nbasis=0;nbasis<DIM+1;nbasis++) {
      for (i=0;i<DIM;i++) {
        n[i]=nbasis+i+1;
        if (n[i]>=DIM+1) n[i]-=(DIM+1);
      }

    a0=cell[ncell].node[n[0]]->X()-cell[ncell].node[nbasis]->X();
    a1=cell[ncell].node[n[1]]->X()-cell[ncell].node[nbasis]->X();
    a2=cell[ncell].node[n[2]]->X()-cell[ncell].node[nbasis]->X();

    detA=a0(0)*(a1(1)*a2(2)-a1(2)*a2(1))-
      a1(0)*(a0(1)*a2(2)-a0(2)*a2(1))+
      a2(0)*(a0(1)*a1(2)-a0(2)*a1(1));

    for (i=0;i<DIM;i++) {
      b=0.0;
      b(i)=1.0;

      det=b(0)*(a1(1)*a2(2)-a1(2)*a2(1))-
        a1(0)*(b(1)*a2(2)-b(2)*a2(1))+
        a2(0)*(b(1)*a1(2)-b(2)*a1(1));
      TMatrix(0,i,nbasis,ncell)=det/detA;

      det=a0(0)*(b(1)*a2(2)-b(2)*a2(1))-
        b(0)*(a0(1)*a2(2)-a0(2)*a2(1))+
        a2(0)*(a0(1)*b(2)-a0(2)*b(1));
      TMatrix(1,i,nbasis,ncell)=det/detA;

      det=a0(0)*(a1(1)*b(2)-a1(2)*b(1))-
        a1(0)*(a0(1)*b(2)-a0(2)*b(1))+
        b(0)*(a0(1)*a1(2)-a0(2)*a1(1));
      TMatrix(2,i,nbasis,ncell)=det/detA;
    }
  }
}

//==================================================
void Cgrid::ChangeLocalVector3D(float* a,char nb1,char nb2) {
  int idim,dn;
  float b[3];

  for (idim=0;idim<DIM;idim++) b[idim]=0.0;
  dn=nb2-nb1;
  if (dn<0) dn+=4;

  switch (dn) {
  case 1:
    b[0]=a[1];
    b[1]=a[2];
    b[2]=-(a[0]+a[1]+a[2]);
    break;
  case 2:
    b[0]=a[2];
    b[1]=-(a[0]+a[1]+a[2]);
    b[2]=a[0];
    break;
  case 3:
    b[0]=-(a[0]+a[1]+a[2]);
    b[1]=a[0];
    b[2]=a[1];
    break;
  default :
    printf("Error: proc. grid.cpp,Cdsmc::ChangeLocalVector3D\n");
    printf("dn=%i\n",dn);
    exit(__LINE__,__FILE__);
  }

  for (idim=0;idim<DIM;idim++) a[idim]=b[idim];
}

//==================================================
void Cgrid::ChangeLocalPositionVector3D(float* a,char nb1,char nb2) {
  int idim,dn;
  float b[3];

  for (idim=0;idim<DIM;idim++) b[idim]=0.0;
  dn=nb2-nb1;
  if (dn<0) dn+=4;

  switch (dn) {
  case 1:
    b[0]=a[1];
    b[1]=a[2];
    b[2]=1.0-(a[0]+a[1]+a[2]);
    break;
  case 2:
    b[0]=a[2];
    b[1]=1.0-(a[0]+a[1]+a[2]);
    b[2]=a[0];
    break;
  case 3:
    b[0]=1.0-(a[0]+a[1]+a[2]);
    b[1]=a[0];
    b[2]=a[1];
    break;
  default :
    printf("Error: proc. grid.cpp,Cdsmc::ChangeLocalPositionVector3D\n");
    printf("dn=%i\n",dn);
    exit(__LINE__,__FILE__);
  }

  for (idim=0;idim<DIM;idim++) a[idim]=b[idim];
}

//==================================================
void Cgrid::InitSurfaceInterpolationData(vector< vector<long int> >&SurfaceSDataGroups) {
  long int ngroup;
  long int nnode,nface,face_type,faceat;
  int idim;
  bool* use_flag;
  double measure, *sum;
  long int n,faces_in_the_group; 

  if (SurfaceSDataGroups.size()==0) return;
  surface_interpolation_data.resize(SurfaceSDataGroups.size()); 

  use_flag=new bool[nfaces];
  sum=new double[nnodes]; 

  for (ngroup=0;ngroup<SurfaceSDataGroups.size();ngroup++) {
    faces_in_the_group=0;
    for (nface=0;nface<nfaces;nface++) use_flag[nface]=false;
    for (nnode=0;nnode<nnodes;nnode++) sum[nnode]=0.0;
    
    for (nface=0;nface<nfaces;nface++) {
      faceat=face[nface].faceat;
      n=0;
      count(SurfaceSDataGroups[ngroup].begin(),SurfaceSDataGroups[ngroup].end(),faceat,n);
      if (n!=0) {
        faces_in_the_group++;
        measure=face[nface].Measure();
        for (idim=0;idim<DIM;idim++) {
          use_flag[face[nface].nodeno[idim]]=true;
          sum[face[nface].nodeno[idim]]+=measure;
        } 
      }
    }

    long int nodes_in_the_group=0,current_node=0;

    for (nnode=0;nnode<nnodes;nnode++) 
      if (use_flag[nnode]==true) nodes_in_the_group++;

    surface_interpolation_data[ngroup].node.resize(nodes_in_the_group);

    for (nnode=0;nnode<nnodes;nnode++) if (use_flag[nnode]==true) {
      surface_interpolation_data[ngroup].node[current_node].nnode=nnode;

      for (nface=0;nface<nfaces;nface++) {
        faceat=face[nface].faceat;
        n=0;
        count(SurfaceSDataGroups[ngroup].begin(),SurfaceSDataGroups[ngroup].end(),faceat,n);
        
        if (n!=0) for (idim=0;idim<DIM;idim++) if (face[nface].nodeno[idim]==nnode) { 
          surface_interpolation_data[ngroup].node[current_node].face.push_back(nface);
          surface_interpolation_data[ngroup].node[current_node].
            weight.push_back(face[nface].Measure()/sum[nnode]);
          break;
        }
      }
    
      surface_interpolation_data[ngroup].node[current_node].face.push_back(-1);
      current_node++;
    }

    surface_interpolation_data[ngroup].nnodes=nodes_in_the_group;
    surface_interpolation_data[ngroup].nfaces=faces_in_the_group;
  }

  delete [] use_flag;
  delete [] sum;
}
      
//==================================================
