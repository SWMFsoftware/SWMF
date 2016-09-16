/*
 * newMars.cpp
 *
 *  Created on: Dec 6, 2011
 *      Author: vtenishe
 */
/*
 * $Log$
 * Revision 1.10  2016/03/24 21:16:56  yunilee
 * All the functions and parameters that have required manual chnages for different model species are now merged for a target-species switch in input file.
 * Mars codes (Mars.cpp Mars.h main.cpp) are now specie-dependent.
 * Model species is switchable by changing "SpeciesList" in input file (single/multispecies).
 * Modifications were made in the functions listed below:
 * 1. Production rate calculation function in Mars.h (double ProductionRateCaluclation_HotO, double ProductionRateCaluclation_HotC)
 * 2. Independent namespaces are created for target species in Mars.h
 * 3. In Mars.cpp: Former is a wrapper function for latter (see below):
 * void newMars::ProductionRateCaluclation --> void newMars::SpeciesProductionRateCaluclation
 * long int newMars::HotAtomProduction_wrapper --> long int newMars::HotOxygen::HotOProduction, long int newMars::HotCarbon::HotCProduction
 * 4. main.cpp will now initialize species assigned by input file
 *
 * Working cases for current verstion: Hot O (Dissociative recombination of O2+), Hot C (Photodissociation of CO)
 *
 * Revision 1.9  2015/10/23 18:51:23  dborovik
 * reverting to the previous versions
 *
 * Revision 1.7  2015/05/26 14:56:27  vtenishe
 * a bug in the loop limits is fixed
 *
 * Revision 1.6  2015/03/31 22:40:30  vtenishe
 * the particle allocation flag is set when the new particle is generated
 *
 * Revision 1.5  2015/03/13 19:17:12  yunilee
 * Changes for Mars code initialization
 *
 * Revision 1.4  2014/11/25 04:52:40  vtenishe
 * few function calls are changed
 *
 * Revision 1.1  2014/10/27 16:32:04  yunilee
 * source files for mars simulation
 *
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
#include "Mars.h"
#include "PhotolyticReactions.h"

int newMars::maxLocalCellOxigenProductionRateOffset=-1;//,newMars::maxLocalBackdroundDensityOffset=-1;
int newMars::minLocalCellOxigenProductionRateOffset=-1,newMars::minLocalBackdroundDensityOffset=-1;
int newMars::sampledLocalInjectionRateOffset=-1;
double *newMars::SampledEscapeRate=NULL;
double *newMars::IonizationRate=NULL;


//spherical mesh for sampling macroscopic parameters
cSphericalVolumeMesh newMars::SphericalSamplingMesh;

//cDataSetMTGCM newMars::Te,newMars::Ti,newMars::Tn,newMars::O,newMars::CO2,newMars::O2p,newMars::Un,newMars::Vn,newMars::Wn,newMars::COp,newMars::CO,newMars::E,newMars::TnEQU,newMars::OEQU,newMars::CO2EQU;
cDataSetMGITM newMars::Te,newMars::Ti,newMars::Tn,newMars::O,newMars::CO2,newMars::CO,newMars::N2,newMars::O2p,newMars::Un,newMars::Vn,newMars::Wn,newMars::E;

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
            res=0.0;//ProductionRateCaluclation(x);
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
    double TotalIonizationRate[PIC::nTotalSpecies*PIC::nTotalThreads];

      int s;
      
      //output data sampled on the spherical mesh
      char fname[300];
      sprintf(fname,"pic.SphericalMesh.out=%i.dat",DataOutputFileNumber);
      SphericalSamplingMesh.plot(fname);
      


    MPI_Gather(SampledEscapeRate,PIC::nTotalSpecies,MPI_DOUBLE,buffer,PIC::nTotalSpecies,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(IonizationRate,PIC::nTotalSpecies,MPI_DOUBLE,TotalIonizationRate,PIC::nTotalSpecies,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (PIC::ThisThread==0) {
      char ChemSymbol[_MAX_STRING_LENGTH_PIC_];

        
        
      for (int thread=1;thread<PIC::nTotalThreads;thread++) for (s=0;s<PIC::nTotalSpecies;s++) buffer[s]+=buffer[thread*PIC::nTotalSpecies+s];
      printf("Total escape rate: \n Specie \t Escape rate [s^{-1}]\n");

      for (s=0;s<PIC::nTotalSpecies;s++) {
        PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
        printf("%s (s=%i):\t %e\n",ChemSymbol,s,buffer[s]/PIC::LastSampleLength);
      }

      printf("\n");
        
        
        for (int thread=1;thread<PIC::nTotalThreads;thread++) for (s=0;s<PIC::nTotalSpecies;s++) TotalIonizationRate[s]+=TotalIonizationRate[thread*PIC::nTotalSpecies+s];
        printf("Total Inoization rate: \n Specie \t Ionization rate [s^{-1}]\n");
        
        for (s=0;s<PIC::nTotalSpecies;s++) {
            PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
            printf("%s (s=%i):\t %e\n",ChemSymbol,s,TotalIonizationRate[s]/PIC::LastSampleLength);
        }
        
        printf("\n");
        
    }

      if (PIC::SamplingMode==_RESTART_SAMPLING_MODE_) for (s=0;s<PIC::nTotalSpecies;s++) {
          
          
      SampledEscapeRate[s]=0.0;
      IonizationRate[s]=0.0;
          
      }
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
      if (_C_SPEC_>=0) {
          buffer.TheoreticalLocalInjectionRate=PIC::VolumeParticleInjection::GetCellInjectionRate(_C_SPEC_,CenterNode);
          buffer.NumericalLocalInjectionRate=*(_C_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
      }
      if (_O_SPEC_>=0) {
          buffer.TheoreticalLocalInjectionRate=PIC::VolumeParticleInjection::GetCellInjectionRate(_O_SPEC_,CenterNode);
          buffer.NumericalLocalInjectionRate=*(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
      }

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

      if (_C_SPEC_>=0) {
          InterpoaltedLocalOxigenProductionRate+=c*(*(_C_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer())));
          TheoreticalOxigenInjectionRate+=c*PIC::VolumeParticleInjection::GetCellInjectionRate(_C_SPEC_,InterpolationList[i])/InterpolationList[i]->Measure;
      }
      if (_O_SPEC_>=0) {
          InterpoaltedLocalOxigenProductionRate+=c*(*(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer())));
          TheoreticalOxigenInjectionRate+=c*PIC::VolumeParticleInjection::GetCellInjectionRate(_O_SPEC_,InterpolationList[i])/InterpolationList[i]->Measure;
      }
      
    if (maxTheoreticalLocalInjectionRate<*((double*)(maxLocalCellOxigenProductionRateOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()))) {
      maxTheoreticalLocalInjectionRate=*((double*)(maxLocalCellOxigenProductionRateOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
    }
  }

  //stored the interpolated data in the associated data buffer
    if (_C_SPEC_>=0) {
        *(_C_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpoaltedLocalOxigenProductionRate;
        *(_C_SPEC_+(double*)(PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset+CenterNode->GetAssociatedDataBufferPointer()))=TheoreticalOxigenInjectionRate;
    }
    if (_O_SPEC_>=0) {
        *(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpoaltedLocalOxigenProductionRate;
        *(_O_SPEC_+(double*)(PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset+CenterNode->GetAssociatedDataBufferPointer()))=TheoreticalOxigenInjectionRate;
    }

  *((double*)(maxLocalCellOxigenProductionRateOffset+CenterNode->GetAssociatedDataBufferPointer()))=maxTheoreticalLocalInjectionRate;
}

//PRODUCTION RATE
void newMars::ProductionRateCaluclation(bool *InjectionFlag,double *Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    bool SingleSpeciesInjectionFlag;
    double SingleSpeciesProductionRate;

    SpeciesProductionRateCaluclation(spec,SingleSpeciesInjectionFlag,SingleSpeciesProductionRate,iCellIndex,jCellIndex,kCellIndex,cell,node);

    InjectionFlag[spec]=SingleSpeciesInjectionFlag;
    Rate[spec]=SingleSpeciesProductionRate;
  }

}


void newMars::SpeciesProductionRateCaluclation(int spec,bool &InjectionFlag,double &Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double xCenter[3]={0.0,0.0,0.0},dx[3]={0.0,0.0,0.0},xmin[3],xmax[3];
  int iLevel,i,LastLevel=0;
  double LastResult=0.0,res;

  const int iLevelMax=10;


  if (node==NULL) {
    InjectionFlag=false;
    Rate=0.0;
    return;
  }


  cell->GetX(xCenter);

  dx[0]=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  dx[1]=(node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_;
  dx[2]=(node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;

  for (i=0;i<DIM;i++) xmin[i]=xCenter[i]-0.5*dx[i],xmax[i]=xCenter[i]+0.5*dx[i];

  for (iLevel=1;iLevel<iLevelMax;iLevel*=2) {
#if _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_
    switch (spec) {
    case _C_SPEC_:
      res=Quadrature::Gauss::Cube::GaussLegendre(DIM,iLevel,ProductionRateCaluclation_HotC,xmin,xmax);
      break;
    case _O_SPEC_:
      res=Quadrature::Gauss::Cube::GaussLegendre(DIM,iLevel,ProductionRateCaluclation_HotO,xmin,xmax);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: unknown species");
    }

#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
    res=Quadrature::Gauss::Cube::GaussLegendre(DIM,iLevel,MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate,xmin,xmax);
    exit(__LINE__,__FILE__,"Error: the block is not generalized for multimple species");
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
        switch (spec) {
        case _C_SPEC_:
	 c=ProductionRateCaluclation_HotC(x);
          break;
        case _O_SPEC_:
	 c=ProductionRateCaluclation_HotO(x);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: unknown species");
        }
#elif _MARS_BACKGROUND_ATMOSPHERE_MODEL_ == _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
        c=MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate(x);
        exit(__LINE__,__FILE__,"Error: the block is not generalized for multimple species");
#else
        exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif


        if (*maxRate<c) *maxRate=c;
        if ((*minRate<0.0)||(*minRate>c)) *minRate=c;

        for (bspec=0;bspec<bspecMax;bspec++) {
          c=PIC::MolecularCollisions::BackgroundAtmosphere::GetBackgroundLocalNumberDensity(bspec,x);
          if (maxBdens[bspec]<c) maxBdens[bspec]=c;
          if ((minBdens[bspec]<0.0)||(minBdens[bspec]>c)) minBdens[bspec]=c;
        }
      }
    }
  }

//    LastResult/=dx[0]*dx[1]*dx[2];

  InjectionFlag=true;
  Rate=LastResult;
}

//BOUNDARY CONDITIONS
int newMars::ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) { //(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int spec;
  double rate;
  PIC::ParticleBuffer::byte *ParticleData;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  rate=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)*startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec);

/*#if _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ == _MARS_ESCAPE_PARTICLES_COUNTING_MODE__COUNT_ALL_
  SampledEscapeRate[spec]+=rate;
#elif  _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ == _MARS_ESCAPE_PARTICLES_COUNTING_MODE__ESCAPE_SPEED_
*/
#if  _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ == _MARS_ESCAPE_PARTICLES_COUNTING_MODE__ESCAPE_SPEED_
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

      p=ProductionRateCaluclation_HotO(x)/maxInjectionRate;
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

    if      (BrRatio<.22) {speed=sqrt(KineticEnergy4/massO);}
    else if (BrRatio<.64) {speed=sqrt(KineticEnergy3/massO);}
    else if (BrRatio<.95) {speed=sqrt(KineticEnergy2/massO);}
    else                  {speed=sqrt(KineticEnergy1/massO);}
/*
    //Determine the energy for new hot C, Branch Ration for 3 channels from dissociative recombination Rosen et al. ['98]
    if      (BrRatio<.761) {speed=sqrt((KineticEnergy1)/((massC*massO+massC*massC)/(2*massO)));}
    else if (BrRatio<.906) {speed=sqrt((KineticEnergy2)/((massC*massO+massC*massC)/(2*massO)));}
    else                   {speed=sqrt((KineticEnergy3)/((massC*massO+massC*massC)/(2*massO)));}
*/
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
           betaO2=massO2/(2*k*Ti.Interpolate(x));
           BGMeanFlowVelocity(bulkVelocity,x);
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

    for (idim=0;idim<3;idim++) {
      velocityO1[idim]=vparent[idim]+speed*Randomposition[idim];
      velocityO2[idim]=vparent[idim]+speed*(-Randomposition[idim]);

      #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(velocityO1[idim])) {
        exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
      }
      #endif
    }
       
   // for (idim=0;idim<3;idim++) {velocityO2[idim]=  /*vparent[idim]+*/    speed*(-Randomposition[idim]);}

    //generate new particles


    //particle ONE
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)newParticleData);

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

    //particle TWO
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)newParticleData);

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

    //increment the source rate counter
    PIC::VolumeParticleInjection::SourceRate[_O_SPEC_]+=2.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep;
    *(_O_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=2.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep/cell->Measure;
         
         //sample production rate on the spherical mesh
         int el=SphericalSamplingMesh.GetSamplingElementNumber(x);
         if (el!=-1) {
             SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el]+=2.0*node->block->GetLocalParticleWeight(_O_SPEC_)/LocalTimeStep;
             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
             if (!isfinite(SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el])) {
                 exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
             }
#endif
}


       }

    return nInjectedParticles;
  }
long int newMars::HotCarbon::HotCProduction(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    long int nInjectedParticles=0;
    
    
    //      return 0;
    
    
    //get the productino rate
    //       bool InjectionFlag[PIC::nTotalSpecies];
    //       double Rate[PIC::nTotalSpecies];
    double ModelParticleInjectionRate,TimeCounter=0.0,LocalTimeStep;
    long int newParticle;
    PIC::ParticleBuffer::byte *newParticleData;
    double LocalParticleWeight=node->block->GetLocalParticleWeight(_C_SPEC_);
    
    /*
     ProductionRateCaluclation(InjectionFlag,Rate,iCellIndex,jCellIndex,kCellIndex,cell,node);
     ModelParticleInjectionRate=Rate[0]*cell->Measure/node->block->GetLocalParticleWeight(0);
     */
    
    
    //two O particles will be produce ar the same time
    ModelParticleInjectionRate=PIC::VolumeParticleInjection::GetCellInjectionRate(_C_SPEC_,cell)/LocalParticleWeight;
    LocalTimeStep=node->block->GetLocalTimeStep(_C_SPEC_);
    
    
    
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
            
            p=ProductionRateCaluclation_HotC(x)/maxInjectionRate;
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
        
#if _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Dissociative_Recombination_COp_
        //Determine the energy for new hot C, Branch Ration for 3 channels from dissociative recombination Rosen et al. ['98]
        if      (BrRatio<.761) {speed=sqrt((KineticEnergy1)/((massC*massO+massC*massC)/(2*massO)));}
        else if (BrRatio<.906) {speed=sqrt((KineticEnergy2)/((massC*massO+massC*massC)/(2*massO)));}
        else                   {speed=sqrt((KineticEnergy3)/((massC*massO+massC*massC)/(2*massO)));}
        
#elif _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Photodissociation_CO_
        //Determine the energy for new hot C from photodissociation of CO, Huebner[92]
        speed=sqrt((KineticEnergy)/((massC*massO+massC*massC)/(2*massO)));
#else
        exit(__LINE__,__FILE__,"Error: HotC source is not defined");
#endif
        
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
        
#if _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Dissociative_Recombination_COp_
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
#elif _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Photodissociation_CO_
        betaO2=massCO/(2*k*Tn.Interpolate(x));
        BGMeanFlowVelocity(bulkVelocity,x);
#else
        exit(__LINE__,__FILE__,"Error: HotC source is not defined");
#endif
        
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
        
        for (idim=0;idim<3;idim++) {
            velocityO1[idim]=vparent[idim]+speed*Randomposition[idim];
            
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
            if (!isfinite(velocityO1[idim])) {
                exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
            }
#endif
        }
        
        // for (idim=0;idim<3;idim++) {velocityO2[idim]=  /*vparent[idim]+*/    speed*(-Randomposition[idim]);}
        
        //generate new particles
        
        
        //particle ONE
        newParticle=PIC::ParticleBuffer::GetNewParticle();
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
        PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)newParticleData);
        
        nInjectedParticles++;
        PIC::BC::nInjectedParticles[_C_SPEC_]++;
        PIC::BC::ParticleProductionRate[_C_SPEC_]+=LocalParticleWeight/LocalTimeStep;
        
        PIC::ParticleBuffer::SetX(x,newParticleData);
        PIC::ParticleBuffer::SetV(velocityO1,newParticleData);
        PIC::ParticleBuffer::SetI(_C_SPEC_,newParticleData);
        
        
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
         PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)newParticleData);
         
         nInjectedParticles++;
         PIC::BC::nInjectedParticles[_C_SPEC_]++;
         PIC::BC::ParticleProductionRate[_C_SPEC_]+=LocalParticleWeight/LocalTimeStep;
         
         PIC::ParticleBuffer::SetX(x,newParticleData);
         PIC::ParticleBuffer::SetV(velocityO2,newParticleData);
         PIC::ParticleBuffer::SetI(_C_SPEC_,newParticleData);
         
         
         //TEST!!!!!!!
         //    PIC::ParticleBuffer::SetX(xTest,newParticleData);
         //TEST END
         
         
         
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);
         
         
         //inject the particle into the system
         //PIC::Mover::MoveParticleBoundaryInjection[0](newParticle,LocalTimeStep-TimeCounter,node,true);
         _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,LocalTimeStep-TimeCounter,node,true);
         */
        //increment the source rate counter
        PIC::VolumeParticleInjection::SourceRate[_C_SPEC_]+=1.0*node->block->GetLocalParticleWeight(_C_SPEC_)/LocalTimeStep;
        *(_C_SPEC_+(double*)(sampledLocalInjectionRateOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=1.0*node->block->GetLocalParticleWeight(_C_SPEC_)/LocalTimeStep/cell->Measure;
        
        //sample production rate on the spherical mesh
        int el=SphericalSamplingMesh.GetSamplingElementNumber(x);
        if (el!=-1) {
            SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el]+=1.0*node->block->GetLocalParticleWeight(_C_SPEC_)/LocalTimeStep;
            
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
            if (!isfinite(SphericalSamplingMesh.SamplingBuffer[LocalSourceRateOffsetSamplingSphericalVolumeMesh][el])) {
                exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
            }
#endif
        }
        
        
    }
    
    return nInjectedParticles;
}


/////////////////////////Ionization////////////////////////////////////////////////////////


double newMars::PhotoIonizationHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData) {
    double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0};
    double frequency=0.0;
    //const double SemiMajorDistance=1.523679; //227,939,100 km = 1.523679 AU, from Wikipedia
    //const double SemiMajorDistance=1.381497; //227,939,100 km = 1.523679 AU, from Wikipedia
    const double SemiMajorDistance=1.665861; //227,939,100 km = 1.523679 AU, from Wikipedia
    
    const double PhotoionizationTimeSolarMin=3.7E6*pow(SemiMajorDistance,2); // [s] at 1AU: 2.2E6 from Lammer 03 or 3.7E6 from Hodges 00.
    //[s-1] at Mars: 8.89E-8 from Schunk 00 or 10.7E-8 from Modolo 05 or 11.7E-8 from Hodges 00.
    const double PhotoionizationTimeSolarMax=1.7E6*pow(SemiMajorDistance,2); // [s] at 1AU: 1.7E6 from Hodges 00.
    //[s-1] at Mars: 2.73E-7 from Schunk 00 or 2.6E-7 from Hodges 00 or 31.25E-8 from Modolo 05.
    
    //const double PhotoionizationTime=(PhotoionizationTimeSolarMin+PhotoionizationTimeSolarMax)/2;//solar moderate
    const double PhotoionizationTime=PhotoionizationTimeSolarMin;//solar min
    //const double PhotoionizationTime=PhotoionizationTimeSolarMax;//solar max
    
    double IonopauseAltitude=300000.0; // Day-side ionopause altitude assumed constant and equal to 300 km from Zhang 93
    
    PIC::ParticleBuffer::GetV(v,modelParticleData);
    PIC::ParticleBuffer::GetX(x,modelParticleData);
    
double AltI=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_);
    //frequency=(1/PhotoionizationTime)*(1/pow(AltI,2));
    frequency=1/PhotoionizationTime;

  
    //determine the positon of the subsolar point (tilt angle = 25.19 deg)
    double Lat,Lon,xsun[3],xunit[3],x_proj;
    Lat=25.19*Pi/180; //Aphelion=25.19, Perihelion=-25.19, Equinox=0.0
    Lon=0.0*Pi/180;
    
    //unit vector for anti-subsolar point
    xunit[0]=-cos(Lat)*cos(Lon);
    xunit[1]=-cos(Lat)*sin(Lon);
    xunit[2]=-sin(Lat);
    
    //x_proj_to_anti-subsolar=unit_vec(unit_vec.x)
    xsun[0]=xunit[0]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[1]=xunit[1]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[2]=xunit[2]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    
    x_proj=sqrt(pow((x[0]-xsun[0]),2)+pow((x[1]-xsun[1]),2)+pow((x[2]-xsun[2]),2));
    
    if (((xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2])>0.0)&&(x_proj<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    
    //if ((x[0]<0.0)&&((sqrt(x[1]*x[1]+x[2]*x[2]))<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    if ((sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_))<IonopauseAltitude) frequency=0.0;
    
    return frequency;
}


double newMars::ChargeExchangeHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData) {
    double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0};
    double frequency=0.0;
    double CXcrossSection=8.0E-20; //[m2] from Stebbings 64.
    double IonopauseAltitude=300000.0; // Day-side ionopause altitude assumed constant and equal to 300 km from Zhang 93
    
    
    double SolarWindProtonDensityx1EUV=2.5E6; //[m-3] from Lammer 03 from Selsis 02 or Zhang 93.
    double SolarWindSpeedx1EUV=4E5; //[m.s-1] from Lammer 03 or Zhang 93 from Newkirk 80.
    //    const double SolarWindProtonDensityx3EUV=8E6; //[m-3] from Lammer 03 from Selsis 02 or 6E6 from Zhang 93.
    //    const double SolarWindSpeedx3EUV=4.8E5; //[m.s-1] from Lammer 03 or 5E5 from Zhang from Newkirk 80.
    //    const double SolarWindProtonDensityx6EUV=4E7; //[m-3] from Lammer 03 from Selsis 02 or 4.3E6 from Zhang 93.
    //    const double SolarWindSpeedx6EUV=7E5; //[m.s-1] from Lammer 03 or Zhang 93 from Newkirk 80.
    
    double SolarWindProtonDensity=SolarWindProtonDensityx1EUV;
    double SolarWindSpeed=SolarWindSpeedx1EUV;

 
double AltI=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_);
    //frequency=CXcrossSection*SolarWindProtonDensity*SolarWindSpeed*(1/pow(AltI,2)); //[s-1]
    frequency=CXcrossSection*SolarWindProtonDensity*SolarWindSpeed; //[s-1]
    
    PIC::ParticleBuffer::GetV(v,modelParticleData);
    PIC::ParticleBuffer::GetX(x,modelParticleData);
    
    //determine the positon of the subsolar point (tilt angle = 25.19 deg)
    double Lat,Lon,xsun[3],xunit[3],x_proj;
    Lat=25.19*Pi/180; //Aphelion=25.19, Perihelion=-25.19, Equinox=0.0
    Lon=0.0*Pi/180;
    
    //unit vector for anti-subsolar point
    xunit[0]=-cos(Lat)*cos(Lon);
    xunit[1]=-cos(Lat)*sin(Lon);
    xunit[2]=-sin(Lat);
    
    //x_proj_to_anti-subsolar=unit_vec(unit_vec.x)
    xsun[0]=xunit[0]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[1]=xunit[1]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[2]=xunit[2]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    
    x_proj=sqrt(pow((x[0]-xsun[0]),2)+pow((x[1]-xsun[1]),2)+pow((x[2]-xsun[2]),2));
    
    if (((xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2])>0.0)&&(x_proj<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    
    //if ((x[0]<0.0)&&((sqrt(x[1]*x[1]+x[2]*x[2]))<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    if ((sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_))<IonopauseAltitude) frequency=0.0;
   
    return frequency;
}


double newMars::ElectronImpactHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData) {
    double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0};
    double frequency=0.0;
    double IonopauseAltitude=300000.0; // Day-side ionopause altitude assumed constant and equal to 300 km from Zhang 93
    
    double SolarWindProtonDensityx1EUV=2.5E6; //[m-3] from Lammer 03 from Selsis 02 or Zhang 93.
    //    const double SolarWindProtonDensityx3EUV=8E6; //[m-3] from Lammer 03 from Selsis 02 or 6E6 from Zhang 93.
    //    const double SolarWindProtonDensityx6EUV=4E7; //[m-3] from Lammer 03 from Selsis 02 or 4.3E6 from Zhang 93.
    
    double SolarWindProtonDensity=SolarWindProtonDensityx1EUV;
    
    double SolarWindElectronDensityAboveBowShock=SolarWindProtonDensity; //[cm-3] from Zhang 93.
    double SolarWindElectronDensityBehindBowShock=4*SolarWindProtonDensity; //[cm-3] from Zhang 93.
    double frequencyPerIncidentElectronAboveBowShock=2E-14; //[m3.s-1] function of Te from Cravens 87.
    double frequencyPerIncidentElectronBehindBowShock=8E-14; //[m3.s-1] function of Te from Cravens 87.
    
    PIC::ParticleBuffer::GetV(v,modelParticleData);
    PIC::ParticleBuffer::GetX(x,modelParticleData);
    
    double Xf=0.6*_RADIUS_(_TARGET_);
    double L=2.081*_RADIUS_(_TARGET_);
    double epsi=1.026;
    double Ybs=sqrt(max((pow(epsi,2)-1)*pow((x[0]-Xf),2)-2*epsi*L*(x[0]-Xf)+pow(L,2),0.1));
    double YY=sqrt(x[1]*x[1]+x[2]*x[2]);
    
    if (YY<Ybs) {frequency=SolarWindElectronDensityBehindBowShock*frequencyPerIncidentElectronBehindBowShock;}
    else {frequency=SolarWindElectronDensityAboveBowShock*frequencyPerIncidentElectronAboveBowShock;}
   

double AltI=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_);
//frequency=frequency*(1/pow(AltI,2)); 
 


    //determine the positon of the subsolar point (tilt angle = 25.19 deg)
    double Lat,Lon,xsun[3],xunit[3],x_proj;
    Lat=25.19*Pi/180; //Aphelion=25.19, Perihelion=-25.19, Equinox=0.0
    Lon=0.0*Pi/180;
    
    //unit vector for anti-subsolar point
    xunit[0]=-cos(Lat)*cos(Lon);
    xunit[1]=-cos(Lat)*sin(Lon);
    xunit[2]=-sin(Lat);
    
    //x_proj_to_anti-subsolar=unit_vec(unit_vec.x)
    xsun[0]=xunit[0]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[1]=xunit[1]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    xsun[2]=xunit[2]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
    
    x_proj=sqrt(pow((x[0]-xsun[0]),2)+pow((x[1]-xsun[1]),2)+pow((x[2]-xsun[2]),2));
    
    if (((xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2])>0.0)&&(x_proj<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    
    //if ((x[0]<0.0)&&((sqrt(x[1]*x[1]+x[2]*x[2]))<(_RADIUS_(_TARGET_)+IonopauseAltitude))) frequency=0.0;
    if ((sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_))<IonopauseAltitude) frequency=0.0;

    return frequency;
}

/*
//Total ionization of hot O corona
double newMars::IonizationHotO(PIC::ParticleBuffer::byte *modelParticleData) {
    double x[3]={0.0,0.0,0.0},v[3]={0.0,0.0,0.0};
    double localTimeStep;
    int spec;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
    
    PIC::ParticleBuffer::GetV(v,modelParticleData);
    PIC::ParticleBuffer::GetX(x,modelParticleData);
    spec=PIC::ParticleBuffer::GetI(modelParticleData);
    localTimeStep=node->block->GetLocalTimeStep(spec);
    
    //    localTimeStep=node->block->GetLocalTimeStep(_O_SPEC_);
    
    double TotalIonizationFrequency=PhotoIonizationHotOFrequency(modelParticleData)+ChargeExchangeHotOFrequency(modelParticleData)+ElectronImpactHotOFrequency(modelParticleData);
    double p=1.0-exp(-localTimeStep*TotalIonizationFrequency);
    //    double p=1.0-exp(-1.0*TotalIonizationFrequency);
    
    //printf("%e, %i",localTimeStep,localTimeStep);
    //exit(__LINE__,__FILE__,"Timestep printing");
    
    return (rnd()<p)? true:false;
}*/


//default function for calculation of the lifetile
double newMars::TotalLifeTime(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    static bool initflag=false;
    static double ReactionLifeTimeTable[PIC::nTotalSpecies];
    static bool ReactionAllowedTable[PIC::nTotalSpecies];
    
    //init the lifetime table
    if (initflag==false) {
        initflag=true;
        
        for (int s=0;s<PIC::nTotalSpecies;s++) {
            if (PhotolyticReactions::ModelAvailable(s)==true) {
                ReactionAllowedTable[s]=true;
         //       ReactionLifeTimeTable[s]=1.0/PhotolyticReactions::GetTotalReactionRate(s,HeliocentricDistance);
            }
            else {
                ReactionAllowedTable[s]=false;
       //         ReactionLifeTimeTable[s]=-1.0;
            }
        }
    }
    
    //return the lifetime
    ReactionAllowedFlag=ReactionAllowedTable[spec];
   // return ReactionLifeTimeTable[spec];
   
    
    double v[3]={0.0,0.0,0.0};
    double localTimeStep;
    //int spec;
//    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
    PIC::ParticleBuffer::byte *modelParticleData;

    //get the particle data
    modelParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    //PIC::ParticleBuffer::GetV(v,modelParticleData);
   // PIC::ParticleBuffer::GetX(x,modelParticleData);
    //spec=PIC::ParticleBuffer::GetI(modelParticleData);
    localTimeStep=node->block->GetLocalTimeStep(spec);
    
    //    localTimeStep=node->block->GetLocalTimeStep(_O_SPEC_);
    
    double TotalIonizationFrequency=PhotoIonizationHotOFrequency(modelParticleData)+ChargeExchangeHotOFrequency(modelParticleData)+ElectronImpactHotOFrequency(modelParticleData);
   // double p=1.0-exp(-localTimeStep*TotalIonizationFrequency);
    //    double p=1.0-exp(-1.0*TotalIonizationFrequency);
    
    //return (rnd()<p)? true:false;
  /* 
if (TotalIonizationFrequency==0.0){
    TotalIonizationFrequency=TotalIonizationFrequency;
}
else{ 
    TotalIonizationFrequency=1/TotalIonizationFrequency;}
    */
return TotalIonizationFrequency;

}


//process particle chemical transformation
void newMars::PhotochemicalModelProcessor(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    int *ReactionProductsList,nReactionProducts;
    double *ReactionProductVelocity;
    int ReactionChannel,spec;
    bool PhotolyticReactionRoute;
    PIC::ParticleBuffer::byte *ParticleData;
    double vParent[3],xParent[3],ParentLifeTime;
    double rateI;

    
    
    //get the particle data
    ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    spec=PIC::ParticleBuffer::GetI(ParticleData);
    PIC::ParticleBuffer::GetV(vParent,ParticleData);
    PIC::ParticleBuffer::GetX(xParent,ParticleData);
    
    rateI=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)*node->block->GetLocalParticleWeight(spec)/node->block->GetLocalTimeStep(spec);
    
    bool PhotolyticReactionAllowedFlag;
    
    ParentLifeTime=TotalLifeTime(NULL,spec,ptr,PhotolyticReactionAllowedFlag,node);
    
//    exit(__LINE__,__FILE__,"test0");
    if (PhotolyticReactionAllowedFlag==false) {
        //no reaction occurs -> add the particle to the list
        //the particle is remain in the system
        PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
        PIC::ParticleBuffer::SetPrev(-1,ptr);
        
 //   exit(__LINE__,__FILE__,"test0.5");
        if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
        FirstParticleCell=ptr;
        return;
    }
    
    //exit(__LINE__,__FILE__,"test1");
  /*  //init the reaction tables
    static bool initflag=false;
    static double ProductionYieldTable[PIC::nTotalSpecies][PIC::nTotalSpecies];
    
    if (initflag==false) {
        int iParent,iProduct;
        initflag=true;
        
        for (iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
            ProductionYieldTable[iParent][iProduct]=0.0;
            
            if (PhotolyticReactions::ModelAvailable(iParent)==true) {
                ProductionYieldTable[iParent][iProduct]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
            }
        }
    }*/
    
    //inject the products of the reaction
    double ParentTimeStep,ParentParticleWeight;
    
#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
    ParentParticleWeight=0.0;
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif
    
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
    ParentTimeStep=node->block->GetLocalTimeStep(spec);
#endif
    
    
    //account for the parent particle correction factor
    ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
    
    //the particle buffer used to set-up the new particle data
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);
    
    //copy the state of the initial parent particle into the new-daugher particle (just in case....)
    PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*)tempParticleData,ParticleData);
    
    for (int specProduct=0;specProduct<PIC::nTotalSpecies;specProduct++) if (specProduct!=spec) {
        double ProductTimeStep,ProductParticleWeight;
        double ModelParticleInjectionRate,TimeCounter=0.0,TimeIncrement,ProductWeightCorrection=1.0;
        int iProduct;
        long int newParticle;
        PIC::ParticleBuffer::byte *newParticleData;
        
        
#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
        ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[specProduct];
#else
        ProductParticleWeight=0.0;
        exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif
        
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
        ProductTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[specProduct];
#else
        ProductTimeStep=node->block->GetLocalTimeStep(specProduct);
#endif
        
    /*
        double anpart=ProductionYieldTable[spec][specProduct]*(1.0-exp(-ProductTimeStep/ParentLifeTime))*ParentParticleWeight/ProductParticleWeight;
        int npart=(int)anpart;
        if (anpart-npart>rnd()) npart+=1;
        
        for (int n=0;n<npart;n++) {
            //generate model particle with spec=specProduct
            bool flag=false;
            
            do {
                //generate a reaction channel
                PhotolyticReactions::GenerateReactionProducts(spec,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
                
                //check whether the products contain species with spec=specProduct
                for (iProduct=0;iProduct<nReactionProducts;iProduct++) if (ReactionProductsList[iProduct]==specProduct) {
                    flag=true;
                    break;
                }
            }
            while (flag==false);
   
   
            //determine the velocity of the product specie
            double x[3],v[3],c=rnd();
            
            for (int idim=0;idim<3;idim++) {
                x[idim]=xParent[idim];
                v[idim]=vParent[idim]+ReactionProductVelocity[idim+3*iProduct];
            }
            
            //generate a particle
            PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
            PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
            PIC::ParticleBuffer::SetI(specProduct,(PIC::ParticleBuffer::byte*)tempParticleData);
            
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
            PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
#endif
            
            //apply condition of tracking the particle
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
            PIC::ParticleTracker::InitParticleID(tempParticleData);
            PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,specProduct,tempParticleData,(void*)node);
#endif
            
            
            //get and injection into the system the new model particle
            newParticle=PIC::ParticleBuffer::GetNewParticle(FirstParticleCell);
            newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
            
            PIC::ParticleBuffer::CloneParticle(newParticleData,(PIC::ParticleBuffer::byte *)tempParticleData);
        }
     */
    }
    
            
            
    double probI=exp(-ParentTimeStep*ParentLifeTime);
    double randomN=rnd();

 //determine whether the parent particle is removed
    //if (rnd()<exp(-ParentTimeStep*ParentLifeTime)) {
   
 if (randomN<probI) {
        //the particle is remain in the system
        PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
        PIC::ParticleBuffer::SetPrev(-1,ptr);
      
/* 
if (ParentLifeTime<3.0 && ParentLifeTime>0.0){ 
    printf("%e  %e  %e\n",ParentTimeStep,ParentLifeTime,probI);
    exit(__LINE__,__FILE__,"test2");
    }
*/
        if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
        FirstParticleCell=ptr;
    }
    else {
        //the particle is removed from the system
        PIC::ParticleBuffer::DeleteParticle(ptr);
        IonizationRate[spec]+=rateI;
   // exit(__LINE__,__FILE__,"test3");
        
    }
}
/////////////////////////Ionization////////////////////////////////////////////////////////



   double newMars::LocalTimeStep(int spec,bool& TimeStepLimitationImposed, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {


    //Find maximum velocity of new hot O (need for determining local time step size)
   //double Vmax=sqrt((2*KineticEnergy1)/massC);
       double Vmax;
       if (spec==_O_SPEC_) {
           Vmax=sqrt((newMars::HotOxygen::KineticEnergy)/massO);}
       
       
       else if (spec==_C_SPEC_) {
#if _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Dissociative_Recombination_COp_
           
           Vmax=sqrt((2*newMars::HotCarbon::KineticEnergy1)/massC);}
#elif _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Photodissociation_CO_

           Vmax=sqrt((2*newMars::HotCarbon::KineticEnergy)/massC);}
       
#else
exit(__LINE__,__FILE__,"Error: HotC source is not defined");
#endif

       else {
           exit(__LINE__,__FILE__,"Error: species is not defined");
       }
       

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

       if (spec==_O_SPEC_) {
     for (idim=0;idim<DIM;idim++) {
       accl[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
     }
	}

       else if (spec==_C_SPEC_) {
     for (idim=0;idim<DIM;idim++) {
       accl[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
     }
	}
       else {
           exit(__LINE__,__FILE__,"Error: species is not defined");
       }


   }

//wrapper function
long int newMars::HotAtomProduction_wrapper(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node){
    long int O_nInjectedParticles=0;
    long int C_nInjectedParticles=0;
    
    if (_O_SPEC_>=0) {
        
        O_nInjectedParticles=newMars::HotOxygen::HotOProduction(iCellIndex,jCellIndex,kCellIndex,cell,node);
        
    }
    
    if (_C_SPEC_>=0) {
        
        C_nInjectedParticles=newMars::HotCarbon::HotCProduction(iCellIndex,jCellIndex,kCellIndex,cell,node);
        
    }
   /* 
    if (_O_SPEC_>=0 &&_C_SPEC_>=0) {
        
        O_nInjectedParticles=newMars::HotOxygen::HotOProduction(iCellIndex,jCellIndex,kCellIndex,cell,node);
        C_nInjectedParticles=newMars::HotCarbon::HotCProduction(iCellIndex,jCellIndex,kCellIndex,cell,node);
        
    }*/
    
    //sum up the outputs from the hotXProduction functions
    return O_nInjectedParticles+C_nInjectedParticles;
    
    
}

