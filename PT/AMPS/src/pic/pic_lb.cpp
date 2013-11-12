//$Id$
//the LB model of the internal degrees of freedom

#include "pic.h"

//the total number of the vibrational and rotational modes
//int *PIC::IDF::nTotalVibtationalModes;
//int *PIC::IDF::nTotalRotationalModes;
//double **PIC::IDF::CharacteristicVibrationalTemperature;

//the offset of the rotation and vibration energy in the particle's data state vector
int PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_=-1,PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_=-1;
int PIC::IDF::_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_=-1;
int PIC::IDF::_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_=-1;
int PIC::IDF::_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[PIC::nTotalSpecies];
//double *PIC::IDF::RotationZnumber;


//init the particle's vibrational energy
void PIC::IDF::LB::InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
  int s,nmode;
  double ThetaVib;

  s=PIC::ParticleBuffer::GetI(ParticleDataStart);

  if ((nTotalVibtationalModes[s]==0)||(VibTemp<1.0E-8)) {
    for (nmode=0;nmode<PIC::IDF::nTotalVibtationalModes[s];nmode++) SetVibE(0.0,nmode,ParticleDataStart);
    return;
  }

  for (nmode=0;nmode<nTotalVibtationalModes[s];nmode++) {
    ThetaVib=PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes];

    if (VibTemp<0.0001*ThetaVib) {
      SetVibE(0.0,nmode,ParticleDataStart);
      continue;
    }

    double c,p,EtaVib,Evib,dE,EvibMax;
    double F,Fmax;

//Bird, Eq 5.53
    c=ThetaVib/VibTemp;
    EtaVib=2.0*c/(exp(c)-1.0);

    EvibMax=10.0;
    dE=EvibMax/10000.0;

    Fmax=pow(dE,EtaVib/2.0-1.0)*exp(-EtaVib/2.0*dE);

    do {
      Evib=dE+rnd()*(EvibMax-dE);
      F=pow(Evib,EtaVib/2.0-1.0)*exp(-EtaVib/2.0*Evib);
      p=F/Fmax;
    }
    while (rnd()>p);

    Evib*=EtaVib/2.0*Kbol*VibTemp;
    SetVibE(Evib,nmode,ParticleDataStart);
  }
}

void PIC::IDF::LB::InitRotTemp(double RotTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
  int s;
  double RotDF,Er,p;

  s=PIC::ParticleBuffer::GetI(ParticleDataStart);
  RotDF=nTotalRotationalModes[s];

  if (RotDF<1.0E-8) {
    SetRotE(0.0,ParticleDataStart);
    return;
  }

  if (RotDF==2.0) Er=-log(rnd());
  else {
    //Bird: Eq. 11.23
    do {
      Er=10.0*rnd();
      p=pow(Er/(RotDF/2.0-1.0),RotDF/2.0-1.0)*exp(-Er+RotDF/2.0-1.0);
    }
    while (rnd()>p);
  }

  Er*=Kbol*RotTemp;
  SetRotE(Er,ParticleDataStart);
}

double PIC::IDF::LB::GetCellRotTemp(int s,PIC::Mesh::cDataCenterNode* cell) {
  double RotDF,Erot,w,wtot=0.0,res=0.0;
  int sStart,sStop;

  if (s==-1) sStart=0,sStop=PIC::nTotalSpecies;
  else sStart=s,sStop=s+1;

  for (s=sStart;s<sStop;s++) {
    RotDF=nTotalRotationalModes[s];
    Erot=GetCellMeanRotE(s,cell);
    w=cell->GetNumberDensity(s);

    if ((w<1.0E-8)||(Erot<1.0E-300)||(RotDF<1.0E-10)) continue;

    res+=w*2.0/(RotDF*Kbol)*Erot;
    wtot+=w;
  }

  if (wtot>0.0) res/=wtot;
  return res;
}

double PIC::IDF::LB::GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell) {
  return GetCellVibTemp(-1,s,cell);
}

double PIC::IDF::LB::GetCellVibTemp(int nmode,int s,PIC::Mesh::cDataCenterNode* cell) {
  int nmodeStart,nmodeStop,sStart,sStop;
  double ThetaVib,Ev,EtaVib,EtaVibSumm,TempOfTheMode,VibTempSumm,Tvib;
  double w,wtot=0.0,res=0.0;

  if (s==-1) nmode=-1,sStart=0,sStop=PIC::nTotalSpecies;
  else sStart=s,sStop=s+1;

  for (s=sStart;s<sStop;s++) {
    EtaVibSumm=0.0,VibTempSumm=0.0;

    if (nmode==-1) nmodeStart=0,nmodeStop=PIC::IDF::nTotalVibtationalModes[s];
    else nmodeStart=nmode,nmodeStop=nmode+1;

    for(nmode=nmodeStart;nmode<nmodeStop;nmode++) {
      ThetaVib=PIC::IDF::CharacteristicVibrationalTemperature[nmode+s*PIC::IDF::nSpeciesMaxVibrationalModes];

  //Get effective temperature of the mode: Bird. Eq 11.26
  //Get effective number of vibrational degrees of freedom for the mode: Bird. Eq 11.28, 5.53
      Ev=GetCellMeanVibE(nmode,s,cell);

      double a;
      a=Ev/(Kbol*ThetaVib);
      if (a>1.0E-300) {
        TempOfTheMode=ThetaVib/log(1.0+1.0/a);
        EtaVib=2.0*Ev/(Kbol*TempOfTheMode);
      }
      else {
        TempOfTheMode=0.0;
        EtaVib=0.0;
      }

      EtaVibSumm+=EtaVib;
      VibTempSumm+=EtaVib*TempOfTheMode;
    }

    Tvib=(EtaVibSumm>1.0E-300) ? VibTempSumm/EtaVibSumm : 0.0;

    w=cell->GetNumberDensity(s);
    res+=w*Tvib;
    wtot+=w;
  }

  if (wtot>0.0) res/=wtot;
  return res;
}

double PIC::IDF::LB::GetCellMeanRotE(int s,PIC::Mesh::cDataCenterNode* cell) {
  double wtot=0.0,res=0.0;
  int sStart,sStop;

  if (s==-1) sStart=0,sStop=PIC::nTotalSpecies;
  else sStart=s,sStop=s+1;

  for (s=sStart;s<sStop;s++) {
    wtot+=(*(s+(double*)(cell->associatedDataPointer+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_)));
    res+=(*(s+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_)));
  }

  if (wtot>0.0) res/=wtot;
  return res;
}

double PIC::IDF::LB::GetCellMeanVibE(int nmode,int s,PIC::Mesh::cDataCenterNode* cell) {
  double wtot=0.0,res=0.0;
  int sStart,sStop,nmodeStart,nmodeStop;

  if (s==-1) sStart=0,sStop=PIC::nTotalSpecies;
  else sStart=s,sStop=s+1;

  for (s=sStart;s<sStop;s++) {
    wtot+=(*(s+(double*)(cell->associatedDataPointer+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_)));

    if (nmode==-1) nmodeStart=0,nmodeStop=PIC::IDF::nTotalVibtationalModes[s];
    else nmodeStart=nmode,nmodeStop=nmode+1;

    for(nmode=nmodeStart;nmode<nmodeStop;nmode++) {
      res+=(*(nmode+(double*)(cell->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s])));
    }
  }

  if (wtot>0.0) res/=wtot;
  return res;
}

void PIC::IDF::LB::RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell) {
  PIC::ParticleBuffer::byte* ptr[2];
  int s[2],nptr;
  double TempIndex,Ec,TrDF;
  double m[2],mr;

//RedistributionEnergyFlag==false <= none of internal modes were relaxed
//RedistributionEnergyFlag==true  <= at least one of internal modes was relaxed
  bool RedistributionEnergyFlag=false;

  ptr[0]=ptr0;
  ptr[1]=ptr1;

  s[0]=PIC::ParticleBuffer::GetI(ptr0);
  s[1]=PIC::ParticleBuffer::GetI(ptr1);
  TempIndex=GetTempIndex(s[0],s[1]);
  TrDF=2.5-TempIndex;

  m[0]=PIC::MolecularData::GetMass(s[0]);
  m[1]=PIC::MolecularData::GetMass(s[1]);
  mr=m[0]*m[1]/(m[0]+m[1]);

  Ec=0.5*mr*vrel*vrel;

  for (nptr=0;nptr<2;nptr++) {

//Vibrational-Translational (VT) relaxaton
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_  == _PIC_MODE_ON_
    if ((this_dsmc->idf.VT_model_flag==true)&&(mol.GetNVibModes(s[nptr])!=0)) {
      double Evib,ThetaVib,EtaVib,VibDF;

      Evib=GetVibE(0,ptr[nptr]);
      ThetaVib=mol.GetThetaVib(0,s[nptr]);

  //Bird, Eq.5.52 , 1.31
      double a,Tvib; //,Temp_coll;
      a=Evib/(Kbol*ThetaVib);
      if (a>1.0E-300) {
        Tvib=ThetaVib/log(1.0+1.0/a);
        EtaVib=2.0*Evib/(Tvib*Kbol);
      }
      else {
        Tvib=0.0;
        EtaVib=0.0;
      }

      VibDF=EtaVib/2.0;


      printf("There are some pronbems in void Cdsmc::Cidf::Cqlb::RedistEnergy's VT relaxation\n");
      exit(__LINE__,__FILE__);

      //    Temp_coll=(Ec+Evib)/(TrDF+VibDF)/Kbol;

      //    Temp_coll=this_dsmc->GetTrTemp(0,0);
      //    a=ThetaVib/Temp_coll;
      //    VibDF=a/(exp(a)-1.0);

      /*
          switch (nptr) {
          case 0 :
            Zv=mol.GetZvib(s[0],Temp_coll,s[1]);
            break;
          case 1 :
            Zv=mol.GetZvib(s[1],Temp_coll,s[0]);
            break;
          }
      */


        if (rnd()<0.5) { //1.0/Zv) {
          double c,p,a1,a2;

          Ec+=Evib;
          a1=VibDF-1.0,a2=TrDF-1.0;

          if (fabs(a1)<1.0E-5) Evib=(1.0-pow(rnd(),1.0/TrDF))*Ec;
          else if (a1>0.0) {
            do {
              c=rnd();
              p=pow((a1+a2)/a1*c,a1)*pow((a1+a2)/a2*(1.0-c),a2);
            }
            while (rnd()>p);

            Evib=c*Ec;
          }
          else {
            double Pmax;
            c=0.0001*Ec;
            Pmax=pow(c,a1)*pow(1.0-c,a2);

            do {
              c=rnd()*Ec;
              p=pow(c,a1)*pow(1.0-c,a2);
            }
            while (Pmax*rnd()>p);

            Evib=c;
          }

          Ec-=Evib;
          if (change_flag[nptr]==true) SetVibE(Evib,0,ptr[nptr]);

          RedistributionEnergyFlag=true;
          break;
        }
      }
#endif

  //Vibrational-Vibrational (VV) relaxation
#if  _PIC_INTERNAL_DEGREES_OF_FREEDOM__VV_RELAXATION_MODE_  == _PIC_MODE_ON_
    exit(__LINE__,__FILE__,"not implemented");
#endif

  //Rotational-Translational (TR) relaxation
#if  _PIC_INTERNAL_DEGREES_OF_FREEDOM__TR_RELAXATION_MODE_  == _PIC_MODE_ON_
    double RotDF=PIC::IDF::nTotalRotationalModes[s[nptr]];
    double Zrot=PIC::IDF::RotationZnumber[s[nptr]];

    if (RotDF>0.0) if (rnd()<1.0/Zrot) {
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

      if (ChangeParticlePropertiesFlag[nptr]==true) SetRotE(Erot,ptr[nptr]);
      Ec-=Erot;

      RedistributionEnergyFlag=true;
      break;
    }
#endif

  }

  if (RedistributionEnergyFlag==true) vrel=sqrt(2.0*Ec/mr);
}

//request the data for the model
int PIC::IDF::LB::RequestSamplingData(int offset) {
  int SampleDataLength=0;

  _TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  _ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_=offset+SampleDataLength;
  SampleDataLength+=PIC::nTotalSpecies*sizeof(double);

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    _VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s]=offset+SampleDataLength;
    SampleDataLength+=nTotalVibtationalModes[s]*sizeof(double);
  }

  return SampleDataLength;
}

void PIC::IDF::LB::Init_BeforeParser() {
  //request the additional particle data
  long int offset;
  int DataLength;

  DataLength=(1+nSpeciesMaxVibrationalModes)*sizeof(double);
  PIC::ParticleBuffer::RequestDataStorage(offset,DataLength);

  _ROTATIONAL_ENERGY_OFFSET_=offset;
  _VIBRATIONAL_ENERGY_OFFSET_=offset+sizeof(double);

  //request sampling data
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestSamplingData);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);

}

//print the model data
void PIC::IDF::LB::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Trot\", \"Tvib\"");
}

void PIC::IDF::LB::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double Trot,Tvib;

  if (pipe->ThisThread==CenterNodeThread) {
    Trot=GetCellRotTemp(DataSetNumber,CenterNode);
    Tvib=GetCellVibTemp(DataSetNumber,CenterNode);
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv(Trot,CenterNodeThread);
      pipe->recv(Tvib,CenterNodeThread);
    }

    fprintf(fout,"%e  %e",Trot,Tvib);
  }
  else {
    pipe->send(Trot);
    pipe->send(Tvib);
  }
}

void PIC::IDF::LB::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i,s,nmode;

  for (s=0;s<PIC::nTotalSpecies;s++) for (i=0;i<nInterpolationCoeficients;i++) {
    *(s+(double*)(CenterNode->GetAssociatedDataBufferPointer()+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_))+=
      *(s+(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+_TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_));

    *(s+(double*)(CenterNode->GetAssociatedDataBufferPointer()+_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_))+=
      *(s+(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+_ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_));

    for (nmode=0;nmode<nTotalVibtationalModes[s];nmode++) {
      *(nmode+(double*)(CenterNode->GetAssociatedDataBufferPointer()+_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s]))+=
        *(nmode+(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+_VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[s]));
    }
  }
}


