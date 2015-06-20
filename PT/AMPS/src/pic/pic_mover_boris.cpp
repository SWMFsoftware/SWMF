//$Id$


#include "pic.h"

void PIC::Mover::BorisSplitAcceleration_default(double *accl, double *rotation, int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  /* function finds acceleration and splits into a simple and gyroscopic parts
   * for Lorentz force (considered here)
   * a_{total} = (q/m) * E + (q/m) V\times B = a_{simple} + \Omega\times V
   *             \_simple_/ \_ gyroscopic _/
   * acceleration due to Lorentz Force only is considered 
   **********************************************************************/
  double accl_LOCAL[3]={0.0}, rotation_LOCAL[3]={0.0};
  
#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
  // find electro-magnetic field
  double E[3],B[3];
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
  // if no coupler is used -> use characteristic values
  memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
  memcpy(B,Exosphere::swB_Typical,3*sizeof(double));
#else 
  // if coupler is used -> get values from it
  //......................................................................
  // find the cell based on particles' position x and block startNode
  // input: x, startNode; output: nd, i,j,k
  long int nd;  // cell's number in the block
  int i,j,k;    // cell's coordinates in the block
  
  // fail-safe check: if the block doesn't exist => exit
  if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");

  // flag: true - exit if a point is not found in the block / false: don't
  nd = PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false);

  // fail-safe check: if a point isn't found, try seacrhing in other blocks
  if (nd==-1) {
    // try to found the block the point is in;
    // starting point for search is block startNode
    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false);
    // if still not found => exit

    if (nd==-1) exit(__LINE__,__FILE__,"Error: the cell is not found");
  }
  
  //......................................................................
  // finally, get fields' values at the cell
  PIC::CPLR::GetBackgroundElectricField(E,x,nd,startNode);
  PIC::CPLR::GetBackgroundMagneticField(B,x,nd,startNode);
#endif//_PIC_COUPLER_MODE_
  
  //......................................................................
  // calculate acceleraton due to Lorentz force
  double ElectricCharge=0.0,mass,Charge2Mass;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if ((_DUST_SPEC_<=spec) && (spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) {
    char ParticleData[PIC::ParticleBuffer::ParticleDataLength];

    memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

    ElectricCharge=ElectricallyChargedDust::GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);
    mass=ElectricallyChargedDust::GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);
  }
  else {
    ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    mass=PIC::MolecularData::GetMass(spec);
  }
#else
  ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
  mass=PIC::MolecularData::GetMass(spec);
#endif //_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_

  if (ElectricCharge!=0.0) {
    Charge2Mass=ElectricCharge/mass;

    for (int idim = 0; idim<DIM; idim++){
      accl_LOCAL[idim]    += Charge2Mass*E[idim];
      rotation_LOCAL[idim]-= Charge2Mass*B[idim];
    }
  }
#endif//_FORCE_LORENTZ_MODE_


  //calculate the acceleration due to gravity where applied
#if _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_
  double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
  }
#endif //_TARGET_ != _TARGET_NONE_


  memcpy(accl,    accl_LOCAL,    3*sizeof(double));
  memcpy(rotation,rotation_LOCAL,3*sizeof(double));
}


int PIC::Mover::Boris(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode){
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  double dtTemp;
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0},vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
  int idim;

  long int LocalCellNumber;
  int i,j,k,spec;

  PIC::Mesh::cDataCenterNode *cell;
  bool MovingTimeFinished=false;

  double xmin,xmax;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

  static long int nCall=0;
  nCall++;

  memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
  memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));
  
  MovingTimeFinished=true;

  double acclInit[3],rotInit[3],vv;
  
#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
  _PIC_PARTICLE_MOVER__BORIS_SPLIT_ACCELERATION_(acclInit,rotInit,spec,ptr,xInit,vInit,startNode);
#endif
  

  // Integrate the equations of motion
  /************************ Boris method: *************************
   * dV/dt    = A        + \Omega \times V                        *
   * X_1      = X_0      + dt * V_{+1/2}                          *
   * V_{+1/2} = U        + 0.5*dt*A_0                             *
   * U        = u        + (u + (u\times h))\times s              *
   * u        = V_{-1/2} + 0.5*dt*A_0                             *
   * h        =-0.5*dt*\Omega                                     *
   * s        = 2*h/(1+|h|^2)                                     *
   ****************************************************************/
  double u[3]={0.0}, U[3]={0.0};
  dtTemp=dtTotal/2.0;
  u[0]=vInit[0]+dtTemp*acclInit[0];
  u[1]=vInit[1]+dtTemp*acclInit[1];
  u[2]=vInit[2]+dtTemp*acclInit[2];

  double h[3];
  h[0]=-dtTemp*rotInit[0];
  h[1]=-dtTemp*rotInit[1];
  h[2]=-dtTemp*rotInit[2];
  double h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  double uh = u[0]*h[0]+u[1]*h[1]+u[2]*h[2];

  U[0] = ( (1-h2)*u[0] + 2*(u[1]*h[2]-h[1]*u[2]+uh*h[0]) ) / (1+h2);
  U[1] = ( (1-h2)*u[1] + 2*(u[2]*h[0]-h[2]*u[0]+uh*h[1]) ) / (1+h2);
  U[2] = ( (1-h2)*u[2] + 2*(u[0]*h[1]-h[0]*u[1]+uh*h[2]) ) / (1+h2);


  vFinal[0]=U[0]+dtTemp*acclInit[0];
  vFinal[1]=U[1]+dtTemp*acclInit[1];
  vFinal[2]=U[2]+dtTemp*acclInit[2];

  xFinal[0]=xInit[0]+dtTotal*vFinal[0];
  xFinal[1]=xInit[1]+dtTotal*vFinal[1];
  xFinal[2]=xInit[2]+dtTotal*vFinal[2];
  
  

  //rotate the final position
#if _PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_
    // rotate to the y=0 plane
    double cosPhi, sinPhi, vTmpX, vTmpY;
    double xNormFinal = pow(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1], 0.5);
    cosPhi = xFinal[0] / xNormFinal;
    sinPhi = xFinal[1] / xNormFinal;
    xFinal[0] = xNormFinal;
    xFinal[1] = 0.0;
    vTmpX = vFinal[0]; vTmpY = vFinal[1];
    vFinal[0] = vTmpX*cosPhi + vTmpY*sinPhi;
    vFinal[1] =-vTmpX*sinPhi + vTmpY*cosPhi;
#endif //_PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_

  newNode=PIC::Mesh::mesh.findTreeNode(xFinal,startNode);

  //advance the particle's position and velocity
  //interaction with the faces of the block and internal surfaces
  
  //check whether the particle trajectory is intersected the spherical body
#if _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_
  double rFinal2;

  if ((rFinal2=xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_)) {
    //the particle is inside the sphere -> apply the boundary condition procedure
    int code;
    cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptor=startNode->InternalBoundaryDescriptorList;
    cInternalSphericalData *Sphere=(cInternalSphericalData*)(InternalBoundaryDescriptor->BoundaryElement);

    //move the particle location at the surface of the sphere
    double r=sqrt(rFinal2);

    for (int idim=0;idim<DIM;idim++) xFinal[idim]*=_RADIUS_(_TARGET_)/r;

    //apply the boundary condition
    code=Sphere->ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)startNode,InternalBoundaryDescriptor->BoundaryElement);

    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
  }
#endif //_TARGET_ == _TARGET_NONE_


  if (newNode==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;
    
    //call the function that process particles that leaved the coputational domain
    switch(code){
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }



  //adjust the value of 'startNode'
  startNode=newNode;
  memcpy(vInit,vFinal,3*sizeof(double));
  memcpy(xInit,xFinal,3*sizeof(double));
  
  //save the trajectory point
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::RecordTrajectoryPoint(xInit,vInit,spec,ParticleData);
#endif
  
  
  
  
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
#if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData);
#endif
#endif
  
  

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xFinal,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
  


  PIC::Mesh::cDataBlockAMR *block=startNode->block;
  long int tempFirstCellParticle=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);


  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ptr;


  return _PARTICLE_MOTION_FINISHED_;

  
}
