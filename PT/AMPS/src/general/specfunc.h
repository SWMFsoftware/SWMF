#ifndef SPECFUNC
#define SPECFUNC

long int nint(double);

void rnd_seed();
double rnd();
double erf(double);
double gam(double);

void PrintErrorLog(char*);


template<class TMesh>
bool GetGradient(double* gradQ,double cellQ,double* Q,long int ncell,TMesh &grid) {
  int counter,idim,pface;
  long int neib;
  double dx2,x[4][3],x0[3],df[4];
  double A[3][3],aa[3][3],af[3],detaa;

  switch (DIM) {
  case 1:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,dx2=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      gradQ[0]+=(Q[pface]-cellQ)*(x[pface][0]-x0[0]);
      dx2=pow(x[pface][0]-x0[0],2);
    }

    gradQ[0]/=(counter!=0) ? dx2 : 1.0; 

    if (counter!=DIM+1) return false;
    break;
  case 2:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,gradQ[1]=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      A[pface][0]=x[pface][0]-x0[0],A[pface][1]=x[pface][1]-x0[1];
      df[pface]=Q[pface]-cellQ;
    }

    if (counter!=DIM+1) return false;

    aa[0][0]=pow(A[0][0],2)+pow(A[1][0],2)+pow(A[2][0],2);
    aa[0][1]=A[0][0]*A[0][1]+A[1][0]*A[1][1]+A[2][0]*A[2][1];
    aa[1][0]=aa[0][1];
    aa[1][1]=pow(A[0][1],2)+pow(A[1][1],2)+pow(A[2][1],2);

    af[0]=A[0][0]*df[0]+A[1][0]*df[1]+A[2][0]*df[2];
    af[1]=A[0][1]*df[0]+A[1][1]*df[1]+A[2][1]*df[2];

    detaa=aa[0][0]*aa[1][1]-aa[1][0]*aa[0][1];

    gradQ[0]=(af[0]*aa[1][1]-af[1]*aa[0][1])/detaa;
    gradQ[1]=(aa[0][0]*af[1]-aa[1][0]*af[0])/detaa;

    break;
  default:
    printf("Error: GetGradient. DIM=%i is not implemented\n",DIM);
    exit(0);
  }

  return true;
}

#endif
   
