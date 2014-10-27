/*
 * newMars.cpp
 *
 *  Created on: Dec 6, 2011
 *      Author: vtenishe
 */
/*
 * $Log$
 * Revision 1.2  2013/03/06 21:42:08  ylee5
 * Te&Ti eqns from Dr.Bougher Fox93
 *
 * Revision 1.1  2013/03/06 21:40:30  ylee5
 *
 */
//The placeholder for the glocal variables of the Mars model


int maxLocalBackdroundDensityOffset=-1;

#include "rnd.h"
#include "pic.h"
#include "newMars.h"

int newMars::maxLocalCellOxigenProductionRateOffset=-1;//,newMars::maxLocalBackdroundDensityOffset=-1;
int newMars::minLocalCellOxigenProductionRateOffset=-1,newMars::minLocalBackdroundDensityOffset=-1;
int newMars::sampledLocalInjectionRateOffset=-1;
double *newMars::SampledEscapeRate=NULL;
//spherical mesh for sampling macroscopic parameters
cSphericalVolumeMesh newMars::SphericalSamplingMesh;

cDataSetMTGCM newMars::Te,newMars::Ti,newMars::Tn,newMars::O,newMars::CO2,newMars::O2p,newMars::Un,newMars::Vn,newMars::Wn,newMars::COp,newMars::CO,newMars::E;

//read the forward collision cross section data
//cDiffCrossSection newMars::ForwardCollisionCrossSection;

cDiffCrossSection ForwardCollisionCrossSection;
//output sampled data to the spherical sampling mesh
int newMars::PrintVariableListSphericalVolumeMesh(FILE* fout) {
    fprintf(fout,"\"Theoretical O Production Rate\", \"O2p density\", \"Te\", \"Local Hot O Density\", \"Numerical Source Rate\"");
    
    return 5;
}

double newMars::GetSampledVariableSphericalSamplingMesh(int var,double* x,int el) {
    double res=0.0;
    
    switch (var) {
        case 0:
            res=ProductionRateCaluclation(x);
            break;
        case 1:
            res=O2p.Interpolate(x);
            break;
        case 2:
            res=Te.Interpolate(x);
            break;
            
        case 3:
            if (PIC::LastSampleLength!=0) {
                res=SphericalSamplingMesh.SamplingBuffer[LocalDensityOffsetSamplingSphericalVolumeMesh][el]/PIC::LastSampleLength/SphericalSamplingMesh.GetVolumeSamplingElement(el);
            }
            break;
        case 4:
            if (PIC::LastSampleLength!=0) {
                res=SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el]/PIC::LastSampleLength/SphericalSamplingMesh.GetVolumeSamplingElement(el);
            }
            break;
            
        default:
            exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    
    return res;
}



int newMars::RequestSamplingData(int offset) {
  sampledLocalInjectionRateOffset=offset;
  return PIC::nTotalSpecies*sizeof(double);
}

int newMars::RequestStaticCellData(int offset) {
  maxLocalCellOxigenProductionRateOffset=offset;
  offset+=sizeof(double);

  maxLocalBackdroundDensityOffset=offset;
  offset+=PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies()*sizeof(double);

  minLocalCellOxigenProductionRateOffset=offset;
  offset+=sizeof(double);

  minLocalBackdroundDensityOffset=offset;
  offset+=PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies()*sizeof(double);

  return 2*(1+PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies())*sizeof(double);
}

  void newMars::OutputSampledModelData(int DataOutputFileNumber) {
    double buffer[PIC::nTotalSpecies*PIC::nTotalThreads];
    int s;
      
      //output data sampled on the spherical mesh
      char fname[300];
      sprintf(fname,"pic.SphericalMesh.out=%i.dat",DataOutputFileNumber);
      SphericalSamplingMesh.plot(fname);
      


    MPI_Gather(SampledEscapeRate,PIC::nTotalSpecies,MPI_DOUBLE,buffer,PIC::nTotalSpecies,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (PIC::ThisThread==0) {
      char ChemSymbol[_MAX_STRING_LENGTH_PIC_];

      for (int thread=1;thread<PIC::nTotalThreads;thread++) for (s=0;s<PIC::nTotalSpecies;s++) buffer[s]+=buffer[thread*PIC::nTotalSpecies+s];
      printf("Total escape rate: \n Specie \t Escape rate [s^{-1}]\n");

      for (s=0;s<PIC::nTotalSpecies;s++) {
        PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
        printf("%s (s=%i):\t %e\n",ChemSymbol,s,buffer[s]/PIC::LastSampleLength);
      }

      printf("\n");
    }

    if (PIC::SamplingMode==_RESTART_SAMPLING_MODE_) for (s=0;s<PIC::nTotalSpecies;s++) SampledEscapeRate[s]=0.0;
  }

  void newMars::PrintVariableList(FILE* fout,int DataSetNumber) {
    fprintf(fout,", \"Maximum cell Theoretical Oxigen Production Rate\", \"Minimum cell Theoretical Oxigen Production Rate\", \"O Theoretical Production Rate\", \"O Mean Cell Numerical Production Rate\", \"Electron Temeprature\", \"Electron Density\", \"O2PLUS Density\"");

    for (int bspec=0;bspec<PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies();bspec++) {
      fprintf(fout,", \"Minimum Background Desnity s=%i\", \"Maximum Background Desnity s=%i\"",bspec,bspec);
    }
  }

void newMars::SampleModelData() {
  //do nothing
}

void newMars::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

  struct cDataExchengeBuffer {
    double NumericalLocalInjectionRate,TheoreticalLocalInjectionRate;
    double maxTheoreticalLocalInjectionRate,minTheoreticalLocalInjectionRate;
  };

  struct cBackgroundDensityBuffer {
    double minValue,maxValue;
  };

  int nBackgroundSpecies=PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies();

  cDataExchengeBuffer buffer;
  cBackgroundDensityBuffer BackgroundDensityBuffer[nBackgroundSpecies];

  if (pipe->ThisThread==CenterNodeThread) {
    buffer.TheoreticalLocalInjectionRate=PIC::VolumeParticleInjection::GetCellInjectionRate(_O_SPEC_,CenterNode);
    buffer.NumericalLocalInjectionRate=*(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));

    buffer.maxTheoreticalLocalInjectionRate=*((double*)(maxLocalCellOxigenProductionRateOffset+CenterNode->GetAssociatedDataBufferPointer()));
    buffer.minTheoreticalLocalInjectionRate=*((double*)(minLocalCellOxigenProductionRateOffset+CenterNode->GetAssociatedDataBufferPointer()));

    for (int bspec=0;bspec<nBackgroundSpecies;bspec++) {
      BackgroundDensityBuffer[bspec].minValue=*(bspec+(double*)(minLocalBackdroundDensityOffset+CenterNode->GetAssociatedDataBufferPointer()));
      BackgroundDensityBuffer[bspec].maxValue=*(bspec+(double*)(maxLocalBackdroundDensityOffset+CenterNode->GetAssociatedDataBufferPointer()));
    }
  }


  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv((char*)&buffer,sizeof(cDataExchengeBuffer),CenterNodeThread);
      pipe->recv((char*)&BackgroundDensityBuffer,nBackgroundSpecies*sizeof(cBackgroundDensityBuffer),CenterNodeThread);
    }

    double *x=CenterNode->GetX();

    if (PIC::LastSampleLength!=0) buffer.NumericalLocalInjectionRate/=PIC::LastSampleLength;

    if (Te.DataValueDefined(x)==true) fprintf(fout,"%e  %e  %e  %e  %e  %e  %e  ",buffer.maxTheoreticalLocalInjectionRate,buffer.minTheoreticalLocalInjectionRate,buffer.TheoreticalLocalInjectionRate,buffer.NumericalLocalInjectionRate,Te.Interpolate(x),1.0E6*E.Interpolate(x),1.0E6*O2p.Interpolate(x));
    else  fprintf(fout,"0.0  0.0  0.0  0.0  0.0  0.0  0.0 ");

    for (int bspec=0;bspec<nBackgroundSpecies;bspec++) fprintf(fout,"%e  %e  ",BackgroundDensityBuffer[bspec].minValue,BackgroundDensityBuffer[bspec].maxValue);
  }
  else {
    pipe->send((char*)&buffer,sizeof(cDataExchengeBuffer));
    pipe->send((char*)&BackgroundDensityBuffer,nBackgroundSpecies*sizeof(cBackgroundDensityBuffer));
  }


}

void newMars::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i;
  double c;


  //interpolate the local productions rate
  double InterpoaltedLocalOxigenProductionRate=0.0,TheoreticalOxigenInjectionRate=0.0,maxTheoreticalLocalInjectionRate=0.0;

  //interpolate the sampled data
  for (i=0;i<nInterpolationCoeficients;i++) {
    c=InterpolationCoeficients[i];

    InterpoaltedLocalOxigenProductionRate+=c*(*(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer())));
    TheoreticalOxigenInjectionRate+=c*PIC::VolumeParticleInjection::GetCellInjectionRate(_O_SPEC_,InterpolationList[i])/InterpolationList[i]->Measure;

    if (maxTheoreticalLocalInjectionRate<*((double*)(maxLocalCellOxigenProductionRateOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()))) {
      maxTheoreticalLocalInjectionRate=*((double*)(maxLocalCellOxigenProductionRateOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
    }
  }

  //stored the interpolated data in the associated data buffer
  *(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpoaltedLocalOxigenProductionRate;
  *(_O_SPEC_+(double*)(PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset+CenterNode->GetAssociatedDataBufferPointer()))=TheoreticalOxigenInjectionRate;

  *((double*)(maxLocalCellOxigenProductionRateOffset+CenterNode->GetAssociatedDataBufferPointer()))=maxTheoreticalLocalInjectionRate;
}

//PRODUCTION RATE
void newMars::ProductionRateCaluclation(bool *InjectionFlag,double *Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double xCenter[3]={0.0,0.0,0.0},dx[3]={0.0,0.0,0.0},xmin[3],xmax[3];
  int iLevel,i,LastLevel=0;
  double LastResult=0.0,res;

  const int iLevelMax=100;


  if (node==NULL) {
    InjectionFlag[0]=false;
    Rate[0]=0.0;
    return;
  }


  cell->GetX(xCenter);

  dx[0]=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  dx[1]=(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_;
  dx[2]=(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;

  for (i=0;i<DIM;i++) xmin[i]=xCenter[i]-0.5*dx[i],xmax[i]=xCenter[i]+0.5*dx[i];

  for (iLevel=1;iLevel<iLevelMax;iLevel*=2) {
#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
    res=Quadrature::Gauss::Cube::GaussLegendre(DIM,iLevel,ProductionRateCaluclation,xmin,xmax);
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
    res=Quadrature::Gauss::Cube::GaussLegendre(DIM,iLevel,MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate,xmin,xmax);
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

    if (res>10.0) {
      if (fabs(LastResult-res)/res<0.01) {
        LastResult=res;
        LastLevel=iLevel;
        break;
      }
      else {
        LastResult=res;
        LastLevel=iLevel;
      }
    }
    else {
      LastResult=res;
      LastLevel=iLevel;
      break;
    }
  }

  //evaluate the maximum and minimum value of the oxigen production rate and lacal background density in a cell
  int j,k,bspec,bspecMax;
  double x[3],c,*maxRate=NULL,*minRate=NULL,*maxBdens=NULL,*minBdens=NULL;

  bspecMax=PIC::MolecularCollisions::BackgroundAtmosphere::GetTotalNumberBackgroundSpecies();

  maxRate=(double*)(maxLocalCellOxigenProductionRateOffset+cell->GetAssociatedDataBufferPointer());
  minRate=(double*)(minLocalCellOxigenProductionRateOffset+cell->GetAssociatedDataBufferPointer());

  maxBdens=(double*)(maxLocalBackdroundDensityOffset+cell->GetAssociatedDataBufferPointer());
  minBdens=(double*)(minLocalBackdroundDensityOffset+cell->GetAssociatedDataBufferPointer());

  *maxRate=0.0;
  *minRate=-1.0;

  for (bspec=0;bspec<bspecMax;bspec++) maxBdens[bspec]=0.0,minBdens[bspec]=-1.0;

  iLevel=LastLevel+1; //increse the number of points to check for the min and max values

  for (i=0;i<iLevel;i++) {
    x[0]=xCenter[0]+(double(i)/double(iLevel-1)-0.5)*dx[0];

    for (j=0;j<iLevel;j++) {
      x[1]=xCenter[1]+(double(j)/double(iLevel-1)-0.5)*dx[1];

      for (k=0;k<iLevel;k++) {
        x[2]=xCenter[2]+(double(k)/double(iLevel-1)-0.5)*dx[2];

#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
        c=ProductionRateCaluclation(x);
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
        c=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate(x);
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif


        if (*maxRate<c) *maxRate=c;
        if ((*minRate<0.0)||(*minRate>c)) *minRate=c;

        for (bspec=0;bspec<bspecMax;bspec++) {
          c=PIC::MolecularCollisions::BackgroundAtmosphere::GetBackgroundNumberDensity(bspec,x);
          if (maxBdens[bspec]<c) maxBdens[bspec]=c;
          if ((minBdens[bspec]<0.0)||(minBdens[bspec]>c)) minBdens[bspec]=c;
        }
      }
    }
  }

//    LastResult/=dx[0]*dx[1]*dx[2];

  InjectionFlag[0]=true;
  Rate[0]=LastResult;
}

//BOUNDARY CONDITIONS
int newMars::ProcessOutsideDomainParticles(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int spec;
  double rate;
  PIC::ParticleBuffer::byte *ParticleData;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  rate=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)*startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec);

#if _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ == _MARS_ESCAPE_PARTICLES_COUNTING_MODE__COUNT_ALL_
  SampledEscapeRate[spec]+=rate;
#elif  _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ == _MARS_ESCAPE_PARTICLES_COUNTING_MODE__ESCAPE_SPEED_
  double x[3],v[3];

  PIC::ParticleBuffer::GetX(x,ParticleData);
  PIC::ParticleBuffer::GetV(v,ParticleData);

  if (0.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])>GravityConstant*_MASS_(_TARGET_)/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])) SampledEscapeRate[spec]+=rate;

#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

  return _PARTICLE_DELETED_ON_THE_FACE_;
}
//HOT OXIGEN BLOCK

     long int newMars::HotOxygen::HotOProduction(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       long int nInjectedParticles=0;


 //      return 0;


       //get the productino rate
//       bool InjectionFlag[PIC::nTotalSpecies];
//       double Rate[PIC::nTotalSpecies];
       double ModelParticleInjectionRate,TimeCounter=0.0,LocalTimeStep;
       long int newParticle;
       PIC::ParticleBuffer::byte *newParticleData;
       double LocalParticleWeight=node->block->GetLocalParticleWeight(_O_SPEC_);

       /*
       ProductionRateCaluclation(InjectionFlag,Rate,iCellIndex,jCellIndex,kCellIndex,cell,node);
       ModelParticleInjectionRate=Rate[0]*cell->Measure/node->block->GetLocalParticleWeight(0);
       */


       //two O particles will be produce ar the same time
       ModelParticleInjectionRate=0.5*PIC::VolumeParticleInjection::GetCellInjectionRate(_O_SPEC_,cell)/LocalParticleWeight;
       LocalTimeStep=node->block->GetLocalTimeStep(_O_SPEC_);



 //      return 0;



//=======================  DEBUG =============
       /*
//test the particle injection at a particular altitude
double xTest[3]={0.0,0.0,200.0E3+_RADIUS_(_TARGET_)};
static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *testNode=PIC::Mesh::mesh.findTreeNode(xTest);

if (testNode!=node) return 0;
ModelParticleInjectionRate=0.1/LocalTimeStep;
*/
//=======================  END DEBUG =============



       if (ModelParticleInjectionRate<=0.0) return 0;

       while ((TimeCounter+=-log(rnd())/ModelParticleInjectionRate)<LocalTimeStep) {

    //Calculate the velocity for each hot O (2 hot O's moving oppositely to each other)
      //Assigin random position of the new hot O on the sphere
    double r=0.0, Randomposition[3], x[3];
    int idim;

    //get particle position
#if _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE__UNIFORM_

    //distribute randomly the initial positino of a particle in a cell
#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM
    do {
      PIC::VolumeParticleInjection::GetRandomCellPosition(x,iCellIndex,jCellIndex,kCellIndex,node);
    } while (Ti.DataValueDefined(x)==false);
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
      PIC::VolumeParticleInjection::GetRandomCellPosition(x,iCellIndex,jCellIndex,kCellIndex,node);
#else
      exit(__LINE__,__FILE__,"Error: the option is net recognized");
#endif

#elif _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE__RATE_DEPENDENT_
    //distribute initial particle position according to the local production within the cell
    double p,maxInjectionRate=*((double*)(maxLocalCellOxigenProductionRateOffset+cell->GetAssociatedDataBufferPointer()));

    do {
#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
      do {
        PIC::VolumeParticleInjection::GetRandomCellPosition(x,iCellIndex,jCellIndex,kCellIndex,node);
      } while (Ti.DataValueDefined(x)==false);

      p=ProductionRateCaluclation(x)/maxInjectionRate;
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
      PIC::VolumeParticleInjection::GetRandomCellPosition(x,iCellIndex,jCellIndex,kCellIndex,node);
      p=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate(x)/maxInjectionRate;
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

    }
    while (p<rnd());
#endif



    for (idim=0;idim<3;idim++) {
      Randomposition[idim]=sqrt(-2.0*log(rnd()))*cos(2*Pi*rnd());
      r+=Randomposition[idim]*Randomposition[idim];}
      r=sqrt(r);
    for (idim=0;idim<3;idim++) {Randomposition[idim]/=r;}

    //Determine the energy for new hot O
    double BrRatio=rnd(); //Branch Ratio for 4 channels from dissociative recombination Kella et al. ['97]
    double speed;
/*
    if      (BrRatio<.22) {speed=sqrt(KineticEnergy4/massO);}
    else if (BrRatio<.64) {speed=sqrt(KineticEnergy3/massO);}
    else if (BrRatio<.95) {speed=sqrt(KineticEnergy2/massO);}
    else                  {speed=sqrt(KineticEnergy1/massO);}
*/
    //Determine the energy for new hot C, Branch Ration for 3 channels from dissociative recombination Rosen et al. ['98]
    if      (BrRatio<.761) {speed=sqrt((KineticEnergy1)/((massC*massO+massC*massC)/(2*massO)));}
    else if (BrRatio<.906) {speed=sqrt((KineticEnergy2)/((massC*massO+massC*massC)/(2*massO)));}
    else                   {speed=sqrt((KineticEnergy3)/((massC*massO+massC*massC)/(2*massO)));}

//=======================  DEBUG =============
//test the particle injection with particular energy
//speed=sqrt(3.0*eV2J*2.0/massO);
//=======================  END DEBUG =============

    //Speed of the parent molecule (O2+)
    double vparent[3],bulkVelocity[3]={0.0,0.0,0.0},c;
    double betaO2;
    double rr=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double Altitude=rr-_RADIUS_(_TARGET_);


#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
	if (Altitude>200.0E3 && Altitude<=300.0E3) {//Dr.Bougher provided Ti eqn, Fox93 Nitrogen paper
	   betaO2=massCO/(2*k*(pow(10,(2.243+((Altitude/1000)-180)/95))));
	   BGMeanFlowVelocity(bulkVelocity,x);
	   }
	else if (Altitude>300.0E3) { //Ti=Te above 300km
	   betaO2=massCO/(2*k*(4200-3750*exp((180-(Altitude/1000))/89.6)));
	   BGMeanFlowVelocity(bulkVelocity,x);
	   }
	else {
	   betaO2=massCO/(2*k*Ti.Interpolate(x));
	   BGMeanFlowVelocity(bulkVelocity,x);
	   }
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
    betaO2=massCO/(2*k*MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetNeutralTemeprature(x));
#endif



        for (idim=0;idim<3;idim++) {
      c=sqrt(-log(rnd())/betaO2);
//      vparent[idim]=c*Randomposition[idim]+bulkVelocity[idim];

      vparent[idim]=cos(2*Pi*rnd())*c+bulkVelocity[idim];


        }



//TEST!!!!!!
//        for (idim=0;idim<3;idim++) vparent[idim]=0.0; //,Randomposition[idim]=0.0;
//        Randomposition[2]=1.0;

        //=============


        /*--------- Take Energy distribution from J. Fox distribution -----------*/
//        memcpy(x,xTest,3*sizeof(double));



//        double KineticEnergy=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetInjectionEnergy(x)*eV2J;
//        speed=sqrt(2.0*KineticEnergy/massO);



    //Calculate velocity (vector) for new hot O
    double velocityO1[3],velocityO2[3];
    for (idim=0;idim<3;idim++) {velocityO1[idim]=  vparent[idim]+    speed*Randomposition[idim];}
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
           if (!isfinite(velocityO1[idim])) {
               exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
           }
#endif
       
   // for (idim=0;idim<3;idim++) {velocityO2[idim]=  /*vparent[idim]+*/    speed*(-Randomposition[idim]);}

    //generate new particles


    //particle ONE
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

    nInjectedParticles++;
    PIC::BC::nInjectedParticles[_O_SPEC_]++;
    PIC::BC::ParticleProductionRate[_O_SPEC_]+=LocalParticleWeight/LocalTimeStep;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(velocityO1,newParticleData);
    PIC::ParticleBuffer::SetI(_O_SPEC_,newParticleData);


    //TEST!!!!!!!
//    PIC::ParticleBuffer::SetX(xTest,newParticleData);
    //TEST END

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);


    //inject the particle into the system
    //PIC::Mover::MoveParticleBoundaryInjection[0](newParticle,LocalTimeStep-TimeCounter,node,true);
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,LocalTimeStep-TimeCounter,node,true);
/*
    //particle TWO
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

    nInjectedParticles++;
    PIC::BC::nInjectedParticles[_O_SPEC_]++;
    PIC::BC::ParticleProductionRate[_O_SPEC_]+=LocalParticleWeight/LocalTimeStep;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(velocityO2,newParticleData);
    PIC::ParticleBuffer::SetI(_O_SPEC_,newParticleData);


    //TEST!!!!!!!
//    PIC::ParticleBuffer::SetX(xTest,newParticleData);
    //TEST END



    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);


    //inject the particle into the system
    //PIC::Mover::MoveParticleBoundaryInjection[0](newParticle,LocalTimeStep-TimeCounter,node,true);
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,LocalTimeStep-TimeCounter,node,true);
*/
    //increment the source rate counter
    PIC::VolumeParticleInjection::SourceRate[_O_SPEC_]+=1.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep;
    *(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=1.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep/cell->Measure;
         
         //sample production rate on the spherical mesh
         int el=SphericalSamplingMesh.GetSamplingElementNumber(x);
         if (el!=-1) {
             SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el]+=1.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep;
             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
             if (!isfinite(SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el])) {
                 exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
             }
#endif
}


       }

    return nInjectedParticles;
  }

//  double LocalTimeStepSize(double m) {

   double newMars::HotOxygen::LocalTimeStep(int spec,bool& TimeStepLimitationImposed, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {


    //Find maximum velocity of new hot O (need for determining local time step size)
    double Vmax=sqrt((2*KineticEnergy1)/massC);

    TimeStepLimitationImposed=true;

    double dx=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;

    if (DIM>1) dx=min(dx,(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_);
    if (DIM>2) dx=min(dx,(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_);

    return dx/Vmax;

  }

   void newMars::TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {


     accl[0]=0.0,accl[1]=0.0,accl[2]=0.0;

   //  return;

     //the gravity force
     double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
     double r=sqrt(r2);
     int idim;

     for (idim=0;idim<DIM;idim++) {
       accl[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
     }



   }

