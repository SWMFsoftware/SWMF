//====================================================
//$Id$
//====================================================
//define specific functions for modeling collisions with the backgound atmosphere

#include "pic.h"
//#include "MTGCM.h"

//collison cross section of the model particles with the background atmosphere
double PIC::MolecularCollisions::BackgroundAtmosphere::GetCollisionCrossSectionBackgoundAtmosphereParticle(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData,PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,double TranslationalEnergy,double cr2) {
  // come up with realistic collision cross section
  return 3.0E-19;
}


double PIC::MolecularCollisions::BackgroundAtmosphere::GetSigmaCrMax(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData) {
  return 1000.0*3.0E-19;

}


void PIC::MolecularCollisions::BackgroundAtmosphere::GenerateBackgoundAtmosphereParticle(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int idim;
  double *xmin,*xmax,*xMiddle,x[3],v[3],beta;

  static const double GlobalNeutalTemeprature=179.0;

  //generate positions of the background particle in the cell
  xmin=node->xmin;
  xmax=node->xmax;
  xMiddle=cell->GetX();

  x[0]=xMiddle[0]+(xmax[0]-xmin[0])/_BLOCK_CELLS_X_*(rnd()-0.5);
  x[1]=xMiddle[1]+(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_*(rnd()-0.5);
  x[2]=xMiddle[2]+(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_*(rnd()-0.5);

  PIC::ParticleBuffer::SetX(x,BackgroundAtmosphereParticleData);

  //generate velocity vector for a particle representing the bacground atmosphere
  beta=GetBackgroundMolecularMass(BackgroundSpecieNumber)/(2*Kbol*GlobalNeutalTemeprature);

  for (idim=0;idim<3;idim++) v[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())/beta);

  PIC::ParticleBuffer::SetV(v,BackgroundAtmosphereParticleData);
  PIC::ParticleBuffer::SetI(BackgroundSpecieNumber,BackgroundAtmosphereParticleData);
}

double PIC::MolecularCollisions::BackgroundAtmosphere::GetBackgroundLocalNumberDensity(int BackgroundSpecieNumber,double *x) {
  double res=0.0;
  res = 1.0e17;
  return res;
}


double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellMeanBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double x[3];
  cell->GetX(x);
  return GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,x);
}

//double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellMaximumBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
// return *(BackgroundSpecieNumber+(double*)(maxLocalBackdroundDensityOffset+cell->GetAssociatedDataBufferPointer()));
//}

double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellLocalBackgroundNumberDensity(double x[3],int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,x);
}




