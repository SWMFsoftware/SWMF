//===================================================
//$Id$
//===================================================


#include "pic.h"

//Sampling
//Erot: this_dsmc->DoubleSDataInt(DoubleSDataOffset,s,ncell)
//Evib for given mode: this_dsmc->DoubleSDataInt(1+3*nmode+DoubleSDataOffset,s,ncell)+=GetVibE(nmode,ptr);
//population of ground vib. state:
//   this_dsmc->DoubleSDataInt(2+3*nmode+DoubleSDataOffset,s,ncell)+=GetVibE(nmode,ptr);
//population of first exited level
//   this_dsmc->DoubleSDataInt(3+3*nmode+DoubleSDataOffset,s,ncell)+=GetVibE(nmode,ptr); 

//int PIC::IDF::qLB::NmaxVibModes=0;
int PIC::IDF::qLB::_VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_=0,PIC::IDF::qLB::_VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_=0;

//===================================================
//request the data for the model
int PIC::IDF::qLB::RequestSamplingData(int offset) {
  int SampleDataLength=0;

  _TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  _ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    _VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s]=offset+SampleDataLength;
    SampleDataLength+=nTotalVibtationalModes[s]*sizeof(double);
  }

  _VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  _VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  return SampleDataLength;
}

//===================================================
void PIC::IDF::qLB::Init() {
  //request the additional particle data
  long int offset;
  int DataLength;

  DataLength=(1+PIC::IDF::nSpeciesMaxVibrationalModes)*sizeof(double);
  PIC::ParticleBuffer::RequestDataStorage(offset,DataLength);

  _ROTATIONAL_ENERGY_OFFSET_=offset;
  _VIBRATIONAL_ENERGY_OFFSET_=offset+sizeof(double);

  //request sampling data
  PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
}

//===================================================
double PIC::IDF::qLB::GetVibE(int nmode_in,PIC::ParticleBuffer::byte *ParticleDataStart) {
  double res=0.0;
  int s,nmode,nmodes_total;

  s=PIC::ParticleBuffer::GetI(ParticleDataStart);
  nmodes_total=PIC::IDF::nTotalVibtationalModes[s];

  for (nmode=((nmode_in>=0) ? nmode_in : 0);((nmode<nmodes_total)&&(nmode_in<0))||(nmode==nmode_in);nmode++) {
    res+=Kbol*PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes]*(*(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_)));
  }

  return res;
} 


//===================================================
void PIC::IDF::qLB::SetVibLevel(double VibQuantumNumber,int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
  *(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_))=VibQuantumNumber;
}

//===================================================
int PIC::IDF::qLB::GetVibLevel(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
  return (int)(*(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_)));
}

//===================================================
void PIC::IDF::qLB::InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
  int s,nmode;

  const int MaxVibQuantumLevel=10;

  s=PIC::ParticleBuffer::GetI(ParticleDataStart);

  if (VibTemp<1.0E-8) {
    for (nmode=0;nmode<PIC::IDF::nSpeciesMaxVibrationalModes;nmode++) SetVibLevel(0,nmode,ParticleDataStart);
    return;
  }
  else for (nmode=PIC::IDF::nTotalVibtationalModes[s];nmode<PIC::IDF::nSpeciesMaxVibrationalModes;nmode++) SetVibLevel(0,nmode,ParticleDataStart);

  for (nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[s];nmode++) {
    double n,p;

//Bird., eq 5.58,5.57
    do {
      n=(int)(rnd()*MaxVibQuantumLevel);
      p=exp(-n*PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes]/VibTemp);
    }
    while (rnd()>p);  

    SetVibLevel(n,nmode,ParticleDataStart);
  }
} 

//===================================================
double PIC::IDF::qLB::GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell) {
  return GetCellVibTemp(-1,s,cell);
}


//===================================================
double PIC::IDF::qLB::GetCellVibTemp(int nmode_in,int s,PIC::Mesh::cDataCenterNode* cell) {
  long int nmode; 
  double Tvib,ThetaVib,Ev,EtaVib,EtaVibSumm=0.0;
  double TempOfTheMode,VibTempSumm=0.0;

  for(nmode=((nmode_in<0) ? 0 : nmode_in);((nmode_in<0)&&(nmode<nTotalVibtationalModes[s]))||((nmode_in>=0)&&(nmode==nmode_in));nmode++) {
    double n0,n1;

    ThetaVib=PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes];
    n0=GetGroundVibLevelPopulationFraction(nmode,s,cell);
    n1=GetFirstExitedVibLevelPopulationFraction(nmode,s,cell);
    Ev=GetCellMeanVibE(nmode,s,cell);

    TempOfTheMode=((n0>0.0)&&(n1>0.0)) ? ThetaVib/log(n0/n1) : 0.0; 
    EtaVib=(TempOfTheMode>1.0-20) ? 2.0*Ev/(Kbol*TempOfTheMode) : 0.0; 

    EtaVibSumm+=EtaVib;
    VibTempSumm+=EtaVib*TempOfTheMode;
  }

  Tvib=(EtaVibSumm>1.0E-300) ? VibTempSumm/EtaVibSumm : 0.0;
  return Tvib;
} 

     
//===================================================
double PIC::IDF::qLB::GetGroundVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell) {
  double wtot,res;

  wtot=(*(s+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_)));
  res=(*(s+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_)));

  if (wtot>0.0) res/=wtot;

  return res;
}


double PIC::IDF::qLB::GetFirstExitedVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell) {
  double wtot,res;

  wtot=(*(s+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_)));
  res=(*(s+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_)));

  if (wtot>0.0) res/=wtot;

  return res;
} 


//===================================================
void PIC::IDF::qLB::RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell) {
/*  PIC::ParticleBuffer::byte ptr[2];
  long int nptr;
  unsigned char s[2];
  double TempIndex,Ec,TrDF;
  double m[2],mr;

//RedistributionEnergyFlag==false <= none of internal modes were relaxed
//RedistributionEnergyFlag==true  <= at least one of internal modes was relaxed 
  bool RedistributionEnergyFlag=false;

  ptr[0]=ptr0; 
  ptr[1]=ptr1;

  s[0]=PIC::ParticleBuffer::GetI(ptr0);
  s[1]=PIC::ParticleBuffer::GetI(ptr1);

  TempIndex=mol.GetTempIndex(s[0],s[1]);
  TrDF=2.5-TempIndex;

  m[0]=mol.GetMass(s[0]);
  m[1]=mol.GetMass(s[1]);
  mr=m[0]*m[1]/(m[0]+m[1]);

  Ec=0.5*mr*vrel*vrel;

  for (nptr=0;nptr<2;nptr++) {

//Vibrational-Translational (VT) relaxaton
  if ((this_dsmc->idf.VT_model_flag==true)&&(mol.GetNVibModes(s[nptr])!=0)) {
    double Evib,ThetaVib;

    Evib=GetVibE(0,ptr[nptr]);
    ThetaVib=mol.GetThetaVib(0,s[nptr]);

//Bird, Eq.5.52 , 1.31

    double Zv,Temp_coll;

exit(__LINE__,__FILE__,"QLB:VT - there are some problems");

/---*
    double a,EtaVib;  

    a=Evib/(Kbol*ThetaVib);
    if (a>1.0E-300) {
      Tvib=ThetaVib/log(1.0+1.0/a);
      EtaVib=2.0*Evib/(Tvib*Kbol);
    }
    else { 
      Tvib=0.0;
      EtaVib=0.0;
    }
*---/

    Temp_coll=this_dsmc->GetTrTemp(s[nptr],ncell);  


    switch (nptr) {
    case 0 :
      Zv=mol.GetZvib(s[0],Temp_coll,s[1]);
      break;
    case 1 :
      Zv=mol.GetZvib(s[1],Temp_coll,s[0]); 
      break;
    }
    
    if (Zv<1.0E-5) Zv=1.0E-5; 


    if (rnd()<1.0/Zv) { 
      double i,i_max,p;
     
      Ec+=Evib; 

//Bird, eq  5.62
      i_max=(int)(Ec/(Kbol*ThetaVib)); 

//Bird, eq 5.61  
      do { 
        i=(int)(rnd()*(i_max+1));
        p=pow(1.0-i*Kbol*ThetaVib/Ec,TrDF-1.0); 
      }
      while(rnd()>p);  
      
      Evib=ThetaVib*Kbol*i;
      Ec-=Evib;

      if (change_flag[nptr]==true) SetVibLevel(i,0,ptr[nptr]);

      RedistributionEnergyFlag=true;
      break;
    } 
  }


//Vibrational-Vibrational (VV) relaxation

//Rotational-Translational (TR) relaxation
    double Zrot,RotDF;
    RotDF=mol.GetRotDF(s[nptr]);
    Zrot=mol.GetZrot(s[nptr]);
    
    if ((this_dsmc->idf.RT_model_flag==true)&&(RotDF>0.0)&&(rnd()<1.0/Zrot)) {
      double a1,a2,Erot,p,c;

      Ec+=GetRotE(ptr[nptr]);
      a1=RotDF/2.0-1.0;
      a2=TrDF-1.0;


      if (fabs(a1)<1.0E-5) Erot=(1.0-pow(rnd(),1.0/TrDF))*Ec;
      else if (a1>0.0) {
        do {
          c=rnd();
          p=pow((a1+a2)/a1*c,a1)*pow((a1+a2)/a2*(1.0-c),a2);
        }
        while (rnd()>p);

        Erot=c*Ec;
      }
      else {
        double Pmax;
        c=0.0001;
        Pmax=pow(c,a1)*pow(1.0-c,a2);

        do {
          c=rnd();
          p=pow(c,a1)*pow(1.0-c,a2);
        }
        while (Pmax*rnd()>p);

        Erot=c*Ec;
      }

      if (change_flag[nptr]==true) SetRotE(Erot,ptr[nptr]);
      Ec-=Erot;

      RedistributionEnergyFlag=true;
      break;
    }
  }
  
  if (RedistributionEnergyFlag==true) vrel=sqrt(2.0*Ec/mr);*/
}

