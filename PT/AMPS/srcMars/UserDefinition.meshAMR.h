//the user defined classes used by the mesher
#ifndef _USER_DEFINED_AMR_MESH_CLASSES_
#define _USER_DEFINED_AMR_MESH_CLASSES_

class cInternalSphericalData_UserDefined {
public :
  int faceat;
  double *SamplingBuffer;

  double *maxIntersectedNodeTimeStep;


  //injection of particles from the sphere
  typedef double (*fInjectionRate)(int spec,void *InternalSphere);
  fInjectionRate InjectionRate;

  //injection of particles from the sphere
  typedef long int (*fInjectionBoundaryCondition)(void *InternalSphere);
  fInjectionBoundaryCondition InjectionBoundaryCondition;

  //interaction of aprticles with the sphere
  typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal,void *startNode,void *InternalSphere);
  fParticleSphereInteraction ParticleSphereInteraction;

  //the flag indicating that the boundary prosidure associated with the sphere is already executed
  bool ProcessedBCflag;

  cInternalSphericalData_UserDefined () {
    faceat=-1,SamplingBuffer=NULL,ProcessedBCflag=false;
    InjectionRate=NULL,InjectionBoundaryCondition=NULL,ParticleSphereInteraction=NULL;

    maxIntersectedNodeTimeStep=NULL;
  }
};

#endif
