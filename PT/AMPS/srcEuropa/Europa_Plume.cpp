//$Id$
//model of the Euroap plume

#include "pic.h"
#include "Europa.h"

double Europa::Plume::xSource_SO[3]; //the location of the source
SpiceDouble Europa::Plume::xSource_IAU[3]; //the location of the source
double Europa::Plume::PlumeExternalNormal_SO[3]; //the normal to the Europa's surface at the location of the plume

//determine the location of the plume source in the "global" frame of reference
void Europa::Plume::SetPlumeLocation() {
  SpiceInt id,n;
  SpiceBoolean  found;
  int idim;

  bodn2c_c("EUROPA",&id,&found);

  //srfrec_c  converts planetocentric latitude and longitude of a surface
  //point on a specified body to rectangular coordinates.
  //planetocentic Longitude == east Longitude -> convert WLongitude to east Longitude before calling srfrec_c
  srfrec_c(id,(360.0-PlumeWLon)*rpd_c(),PlumeLat*rpd_c(),xSource_IAU);


  //transform position of the plume location into the frame where the simulation is performed (Exosphere::SO_FRAME)
  int iSO,iIAU;
  double l=0.0;

  for (iIAU=0;iIAU<3;iIAU++) xSource_IAU[iIAU]*=1.0E3;

  for (iSO=0;iSO<3;iSO++) {
    for (xSource_SO[iSO]=0.0,iIAU=0;iIAU<3;iIAU++) {
      xSource_SO[iSO]+=OrbitalMotion::IAU_to_SO_TransformationMartix[iSO][iIAU]*xSource_IAU[iIAU];
    }

    l+=pow(xSource_SO[iSO],2);
  }

  for (iSO=0,l=sqrt(l);iSO<3;iSO++) PlumeExternalNormal_SO[iSO]=-xSource_SO[iSO]/l;
}

double Europa::Plume::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
  return PlumeSourceRate[spec];
}

bool Europa::Plume::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
  int idim;

  if (PlumeSourceRate[spec]<0.0) return false;

  startNode=PIC::Mesh::mesh.findTreeNode(xSource_SO,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

  //calculate velocity of the injected particle in the IAU frame (relative to the surface of Europa)
  double BulkVelocity[3]={0.0,0.0,0.0};
  PIC::Distribution::InjectMaxwellianDistribution(v_IAU_OBJECT,BulkVelocity,PlumeSourceTemeprature,PlumeExternalNormal_SO,spec);

  //copy the local coordinate to the 'global' variable
  for (idim=0;idim<3;idim++) x_IAU_OBJECT[idim]=xSource_IAU[idim],x_SO_OBJECT[idim]=xSource_SO[idim];

  //transfer the velocity vector into the 'SO' frame where the simualtion is performed
  v_SO_OBJECT[0]=
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_IAU_OBJECT[2])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_IAU_OBJECT[2]);

  v_SO_OBJECT[1]=
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_IAU_OBJECT[2])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_IAU_OBJECT[2]);

  v_SO_OBJECT[2]=
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_IAU_OBJECT[2])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_IAU_OBJECT[0])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_IAU_OBJECT[1])+
      (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_IAU_OBJECT[2]);

  return true;
}
