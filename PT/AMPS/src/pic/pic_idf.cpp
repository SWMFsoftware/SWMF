//===================================================
//$Id$
//===================================================

#ifdef MPI_ON
#include "mpi.h"
#endif

#include <stdlib.h>
#include <stdio.h>

#include "data.h"
#include "dsmc.h"
#include "specfunc.h"
#include "dsmc.dfn"
#include "const.dfn"

//===================================================
void Cdsmc::Cidf::init(Cdsmc* this_dsmc_obj) {

  this_dsmc=this_dsmc_obj;
  this_grid=&this_dsmc_obj->grid;
  this_pbuffer=&this_dsmc_obj->pbuffer;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.init(this_dsmc);
    break;
  case QuantumLarsenBorgnakke :
    qlb.init(this_dsmc);
    break;
  default :
    exit(__LINE__,__FILE__,"cannot find the idf model");
  }
}


//===================================================
double Cdsmc::Cidf::GetRotE(ParticlePtr ptr) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetRotE(ptr);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetRotE(ptr);
    break;
  default :
    exit(__LINE__,__FILE__,"cannot find the idf model");
  }

  return res;
}

//===================================================
void Cdsmc::Cidf::SetRotE(double e,ParticlePtr ptr) {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.SetRotE(e,ptr);
    break;
  default : 
    exit(__LINE__,__FILE__,"cannot find the idf model");
  }  
}


//===================================================
double Cdsmc::Cidf::GetVibE(long int nmode,ParticlePtr ptr) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetVibE(nmode,ptr);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetVibE(nmode,ptr);
    break;
  }

  return res;
}

//===================================================
void Cdsmc::Cidf::RedistEnergy(ParticlePtr ptr0,ParticlePtr ptr1,double& vrel,bool* change_internal_properties_flag,long int ncell,long int nsubcl) {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.RedistEnergy(ptr0,ptr1,vrel,change_internal_properties_flag,ncell,nsubcl); 
    break;
  case QuantumLarsenBorgnakke :
    qlb.RedistEnergy(ptr0,ptr1,vrel,change_internal_properties_flag,ncell,nsubcl);
    break;
  }
}

//===================================================
double Cdsmc::Cidf::GetRotTemp(int s,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetRotTemp(s,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetRotTemp(s,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetRotTemp(int s,long int nsubcl,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetRotTemp(s,nsubcl,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetRotTemp(s,nsubcl,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetVibTemp(int s,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetVibTemp(s,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetVibTemp(s,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetVibTemp(int s,long int nsubcl,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetVibTemp(s,nsubcl,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetVibTemp(s,nsubcl,ncell);
    break;
  }

  return res;
}

//===================================================
void Cdsmc::Cidf::InitRotTemp(double temp,ParticlePtr ptr) {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.InitRotTemp(temp,ptr);
    break;
  case QuantumLarsenBorgnakke :
    qlb.InitRotTemp(temp,ptr);
    break;
  }
}

//===================================================
void Cdsmc::Cidf::InitVibTemp(double temp,ParticlePtr ptr) {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.InitVibTemp(temp,ptr);
    break;
  case QuantumLarsenBorgnakke :
    qlb.InitVibTemp(temp,ptr);
    break;
  }
}

//===================================================
double Cdsmc::Cidf::GetMeanVibE(long int nmode,int s,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetMeanVibE(nmode,s,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetMeanVibE(nmode,s,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetMeanVibE(long int nmode,int s,long int nsubcl,long int ncell) {
 double res=0.0;  

 switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetMeanVibE(nmode,s,nsubcl,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetMeanVibE(nmode,s,nsubcl,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetMeanRotE(int s,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetMeanRotE(s,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetMeanRotE(s,ncell);
    break;
  }

  return res;
}

//===================================================
double Cdsmc::Cidf::GetMeanRotE(int s,long int nsubcl,long int ncell) {
  double res=0.0;

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    res=lb.GetMeanRotE(s,nsubcl,ncell);
    break;
  case QuantumLarsenBorgnakke :
    res=qlb.GetMeanRotE(s,nsubcl,ncell);
    break;
  }

  return res;
}

//===================================================
void Cdsmc::Cidf::InitSampling() {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.Sampling();
    break;
  case QuantumLarsenBorgnakke :
    qlb.Sampling();
    break;
  default :
    exit(__LINE__,__FILE__,"cannot find the idf model");
  }
  
//  InitEnergyDistributionSampling();
}

//===================================================
void Cdsmc::Cidf::Sampling(bool breakSampling) {

  switch (this_dsmc->ModelOfInternalDegreesOfFreedom) {
  case None_internal_degrees_of_freedom :
    break;
  case ClassicalLarsenBorgnakke :
    lb.Sampling();
    break;
  case QuantumLarsenBorgnakke :
    qlb.Sampling();
    break; 
  default :
    exit(__LINE__,__FILE__,"cannot find the idf model");
  }
  
//  EnergyDistributionSampling(breakSampling);
}
  
//===================================================
void Cdsmc::Cidf::InitEnergyDistributionSampling() {
  double de,e;
  long int n,dataset,ncell,nsubcl,NSubCells;
  ParticlePtr ptr;
  unsigned char s;

  if (IDF_EDistrSData.size()<1) return;

  if (IDF_EDistrSData[0].data.size()==0) {
    long int nvibmodes_max;
 
    nvibmodes_max=0;
    for (s=0;s<NS;s++) if (nvibmodes_max<mol.GetNVibModes(s)) nvibmodes_max=mol.GetNVibModes(s); 
    
    for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
      IDF_EDistrSData[dataset].data.init(NI_IDF_EDistrSData,2+nvibmodes_max,NS);
      IDF_EDistrSData[dataset].data=0.0;

      IDF_EDistrSDataInt[dataset].data.init(NI_IDF_EDistrSData,2+nvibmodes_max,NS);
      IDF_EDistrSDataInt[dataset].data=0.0;
    }
  }

  de=Kbol*Tmax_IDF_EDistrSData/NI_IDF_EDistrSData;

  for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
    ncell=IDF_EDistrSDataInt[dataset].ncell;

    NSubCells=1<<(DIM*this_dsmc->grid.cell[ncell].GetSubcellLayer());
    for (nsubcl=0;nsubcl<NSubCells;nsubcl++) {
      ptr=this_dsmc->grid.GetFirstPtr(nsubcl,ncell);

      while (ptr>=0) {
        s=this_pbuffer->GetI(ptr);

//energy destribution (energy of internal degrees of freedom)
//rotational energy destribution
        e=GetRotE(ptr);
        n=(long int)(e/de);
        if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,0,s)+=1.0;

//vidrational energy distribution (over all vibational modes)
        if (mol.GetNVibModes(s)>0) {
          e=GetVibE(-1,ptr);
          n=(long int)(e/de);
          if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,1,s)+=1.0; 
        }

//vidrational energy distribution (for each vibational modes)
        if (mol.GetNVibModes(s)>0) 
          for (long int nmode=0;nmode<mol.GetNVibModes(s);nmode++) {
            e=GetVibE(nmode,ptr); 
            n=(long int)(e/de);
            if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,2+nmode,s)+=1.0;
          }
        
        ptr=this_pbuffer->GetNext(ptr);
      }
    }
  }

  if ((ThisThread==0)&&((this_dsmc->SamplingMode!=UnsteadySM)||
  ((this_dsmc->SamplingMode==UnsteadySM)&&(this_dsmc->SampleCounter==this_dsmc->SampleLength)))) { 
    for (dataset=0;dataset<IDF_EDistrSData.size();dataset++)
      IDF_EDistrSData[dataset].data=IDF_EDistrSDataInt[dataset].data;

    //OutputIDF_EDistibutionDataIntoTECPLOT();
  }
}

//===================================================
void Cdsmc::Cidf::EnergyDistributionSampling(bool breakSampling) {
  double de,e;
  long int n,dataset,ncell,nsubcl,NSubCells;
  ParticlePtr ptr;
  unsigned char s;
  static long int SampleCounterInt=0;

  if (IDF_EDistrSData.size()<1) return;

  de=Kbol*Tmax_IDF_EDistrSData/NI_IDF_EDistrSData;

  for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
    ncell=IDF_EDistrSDataInt[dataset].ncell;

    NSubCells=1<<(DIM*this_dsmc->grid.cell[ncell].GetSubcellLayer());
    for (nsubcl=0;nsubcl<NSubCells;nsubcl++) {
      ptr=this_dsmc->grid.GetFirstPtr(nsubcl,ncell);

      while (ptr>=0) {
        s=this_pbuffer->GetI(ptr);

//energy destribution (energy of internal degrees of freedom)
//rotational energy destribution
        e=GetRotE(ptr);
        n=(long int)(e/de);
        if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,0,s)+=1.0;

//vidrational energy distribution (over all vibational modes)
        if (mol.GetNVibModes(s)>0) {
          e=GetVibE(-1,ptr);
          n=(long int)(e/de);
          if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,1,s)+=1.0;
        }

//vidrational energy distribution (for each vibational modes)
        if (mol.GetNVibModes(s)>0)
          for (long int nmode=0;nmode<mol.GetNVibModes(s);nmode++) {
            e=GetVibE(nmode,ptr);
            n=(long int)(e/de);
            if ((n<NI_IDF_EDistrSData)&&(n>=0)) IDF_EDistrSDataInt[dataset].data(n,2+nmode,s)+=1.0;
          }

        ptr=this_pbuffer->GetNext(ptr);
      }
    }
  }

  SampleCounterInt++;

  if ((breakSampling==true)||(SampleCounterInt==this_dsmc->SampleLength)) {
    switch (this_dsmc->SamplingMode) {
    case UnsteadySM :
      for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
        IDF_EDistrSData[dataset].data=IDF_EDistrSDataInt[dataset].data;
        IDF_EDistrSDataInt[dataset].data=0.0;
      }
      break;
    case SteadySM :
      for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
        IDF_EDistrSData[dataset].data=IDF_EDistrSDataInt[dataset].data;
        IDF_EDistrSDataInt[dataset].data=0.0;
      }
      break;
    case AccumulationSM :
      if (this_dsmc->FlowDataOutputFileNumber<=2) {
        for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {
          IDF_EDistrSData[dataset].data=IDF_EDistrSDataInt[dataset].data;
          IDF_EDistrSDataInt[dataset].data=0.0;
        }
      }
      else {
        for (dataset=0;dataset<IDF_EDistrSData.size();dataset++)
          IDF_EDistrSData[dataset].data=IDF_EDistrSDataInt[dataset].data;
      }
      break;
    default :
      printf("Error: proc Cdsmc::Cidf::EnergyDistributionSamplng\n");
      printf("Unrecognizable sampling mode (SamplingMode=%i)\n",this_dsmc->SamplingMode);
      exit(__LINE__,__FILE__);
    }

    SampleCounterInt=0;
  }
}

//===================================================
void Cdsmc::Cidf::OutputEnergyDistibutionDataIntoTECPLOT() {

/*
  double de,norm;
  long int i,j,dataset,ncell;
  float x[3];
  unsigned char s;
  char fname[50],str[50];
  ofstream fout;
*/

  long int dataset;

  for (dataset=0;dataset<IDF_EDistrSData.size();dataset++) {


printf("void Cdsmc::Cidf::OutputEnergyDistibutionDataIntoTECPLOT is not implemented for subcells\n");
exit(__LINE__,__FILE__);

/*

    ncell=IDF_EDistrSData[dataset].ncell;
    for (int idim=0;idim<DIM;idim++) x[idim]=IDF_EDistrSData[dataset].x[idim];

//normalization of distribution functions
    for (s=0;s<NS;s++) for (j=0;j<2+mol.GetNVibModes(s);j++) {
      norm=0.0;
      for (i=0;i<NI_IDF_EDistrSData;i++)
        norm+=IDF_EDistrSData[dataset].data(i,j,s);
      if (norm>1.0E-8)
        for (i=0;i<NI_IDF_EDistrSData;i++)
          IDF_EDistrSData[dataset].data(i,j,s)/=norm;
    }

//output into file
    for (s=0;s<NS;s++) {
      if (NS==1)
        sprintf(fname,"E%i.%i.dsmc.dat",this_dsmc->FlowDataOutputFileNumber,dataset);
      else
        sprintf(fname,"E%i.%i.%i.dsmc.dat",this_dsmc->FlowDataOutputFileNumber,
          dataset,s+1);

      fout.open(fname);
      fout.precision(5);
      fout.setf(ios::scientific,ios::floatfield);
      fout << "TITLE=\"Internal energy distribution function, cell="<<ncell<<" (x= ";
      for (int idim=0;idim<DIM;idim++) {
        fout <<x[idim];
        if (idim!=DIM-1)
          fout << ",";
        else
          fout << ") " ;
      }
      mol.GetChemSymbol(str,s);
      fout << ",componenet="<<str<<", over "<<this_dsmc->SampleCounter<<" itrs.\""<<endl;
      fout << "VARIABLES=\"e\",\"Erot\",\"Evib\"";
      for (long int nmode=0;(nmode<mol.GetNVibModes(s))&&(mol.GetNVibModes(s)!=1);nmode++)
        fout << ",\"Evib(mode"<< nmode << ")\"";
      fout << endl;

      de=Kbol*Tmax_IDF_EDistrSData/NI_IDF_EDistrSData;

      for (i=0;i<NI_IDF_EDistrSData;i++) {
        long int j_max;
        j_max=(mol.GetNVibModes(s)==1) ? 2 : 2+mol.GetNVibModes(s);
 
        fout << de*(i+0.5) << "  ";
        for (j=0;j<j_max;j++)
          fout << "   " << IDF_EDistrSData[dataset].data(i,j,s);
        fout << endl;
      }
      fout.close();
    }

*/
  }
}



