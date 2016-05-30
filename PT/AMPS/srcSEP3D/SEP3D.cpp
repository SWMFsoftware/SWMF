//$Id$

// source code for functions defined in SEPD3D.h

#include "SEP3D.h"

int SEP3D::Physics::Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode){
  return SEP3D::Physics::Mover_Axisymmetric_SecondOrder(ptr, dtTotal, startNode,false);
}


int SEP3D::Physics::Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, bool FirstBoundaryFlag){

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL,*middleNode=NULL;
  double dtTemp;
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0},vMiddle[3],xMiddle[3],vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
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

  double acclInit[3],acclMiddle[3],vv;
  
#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_
  
  _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclInit,spec,ptr,xInit,vInit,startNode);

#endif
  

  //Integrate the equations of motion
  //predictor
  dtTemp=dtTotal/2.0;
  xMiddle[0]=xInit[0]+dtTemp*vInit[0];
  xMiddle[1]=xInit[1]+dtTemp*vInit[1];
  xMiddle[2]=xInit[2]+dtTemp*vInit[2];
  
  vMiddle[0]=vInit[0]+dtTemp*acclInit[0];
  vMiddle[1]=vInit[1]+dtTemp*acclInit[1];
  vMiddle[2]=vInit[2]+dtTemp*acclInit[2];
  
  // rotate to the y=0 plane
  double xNormMiddle = pow(xMiddle[0]*xMiddle[0]+xMiddle[1]*xMiddle[1], 0.5);
  double cosPhi = xMiddle[0] / xNormMiddle;
  double sinPhi = xMiddle[1] / xNormMiddle;
  xMiddle[0] = xNormMiddle;
  xMiddle[1] = 0.0;
  double vTmpX = vMiddle[0], vTmpY = vMiddle[1];
  vMiddle[0] = vTmpX*cosPhi + vTmpY*sinPhi;
  vMiddle[1] =-vTmpX*sinPhi + vTmpY*cosPhi;
  
  middleNode=PIC::Mesh::mesh.findTreeNode(xMiddle,startNode);
  
  if(middleNode == NULL){
    // the particle has left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

#if _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);
    //rotate acceleration and velocity back
    double acclTmp = acclMiddle[0];
    acclMiddle[0] = acclTmp*cosPhi - acclMiddle[1]*sinPhi;
    acclMiddle[1] = acclTmp*sinPhi + acclMiddle[1]*cosPhi;
    vMiddle[0] = vTmpX;
    vMiddle[1] = vTmpY;
#endif

    
    //adjust the particle moving time
    xFinal[0]=xInit[0]+dtTotal*vMiddle[0];
    xFinal[1]=xInit[1]+dtTotal*vMiddle[1];
    xFinal[2]=xInit[2]+dtTotal*vMiddle[2];
    
    vFinal[0]=vInit[0]+dtTotal*acclMiddle[0];
    vFinal[1]=vInit[1]+dtTotal*acclMiddle[1];
    vFinal[2]=vInit[2]+dtTotal*acclMiddle[2];
    
    //rotate the final position

    // rotate to the y=0 plane
    double xNormFinal = pow(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1], 0.5);
    cosPhi = xFinal[0] / xNormFinal;
    sinPhi = xFinal[1] / xNormFinal;
    xFinal[0] = xNormFinal;
    xFinal[1] = 0.0;
    vTmpX = vFinal[0]; vTmpY = vFinal[1];
    vFinal[0] = vTmpX*cosPhi + vTmpY*sinPhi;
    vFinal[1] =-vTmpX*sinPhi + vTmpY*cosPhi;

    newNode=PIC::Mesh::mesh.findTreeNode(xFinal,middleNode);

    
    //advance the particle's position and velocity
    //interaction with the faces of the block and internal surfaces

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
    PIC::ParticleTracker::RecordTrajectoryPoint(xInit,vInit,spec,ParticleData,(void*)newNode);
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


double SEP3D::GlobalEnergyDistribution::Distribution[nBin];

void SEP3D::GlobalEnergyDistribution::sample() {
  // go through all particles in the domain and sample energy distribution
  
  //aliases
  namespace PB = PIC::ParticleBuffer;
  namespace MD = PIC::MolecularData;
  namespace MS = PIC::Mesh;
  
  // reset the distribution
  for(int iBin=0; iBin<nBin; iBin++) Distribution[iBin] = 0.0;
  
  // go through all particles
  for (long int ptr=0;ptr<PB::MaxNPart;ptr++) {
    PB::byte *PData = PB::GetParticleDataPointer(ptr); 
    
    // skip if particle is not used
    if (!PB::IsParticleAllocated(PData))
      continue;
    
    // the particle is used, sample;
    // first, find the weight of a particle
    double x[3];
    cTreeNodeAMR<MS::cDataBlockAMR> *node=NULL;
    PIC::ParticleBuffer::GetX(x,(PB::byte*)PData);
    node=MS::Search::FindBlock(x);
    int spec  = PB::GetI(ptr);
    double Weight=node->block->GetLocalParticleWeight(spec);
    Weight*=PB::GetIndividualStatWeightCorrection((PB::byte*)PData);
    
    // second, find its total kinetic energy
    double m0 = MD::GetMass(spec);
    
    double mu = PB::GetMagneticMoment(ptr);
    double v[3];PB::GetV(v, ptr);
    
    double B[3];
    PIC::CPLR::InitInterpolationStencil(x,node);
    PIC::CPLR::GetBackgroundMagneticField(B);
    double AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2],0.5);
    
    double E = mu*AbsB + 0.5 * m0 * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    
    
    
    //third, determine which bin it falls into; sample
    if(E >= EnergyMin && E < EnergyMax){
      int iBin = (int) (nBin * (E - EnergyMin) / (EnergyMax - EnergyMin));
      Distribution[iBin] += Weight;
    }
  } 
}

//---------------------------------------------------------------------------
    
//print the distribution function into a file
void SEP3D::GlobalEnergyDistribution::print(char *fname,int spec) {

  //aliases
  namespace MS = PIC::Mesh;
  
  long int idim,nProbe,i,nVariable,thread,offset;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,dInterval=0.0;
  char str[_MAX_STRING_LENGTH_PIC_];
  
  if (MS::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);
  
  
  if (MS::mesh.ThisThread==0) {
    sprintf(str,"%s.dat",fname,nProbe);
    fout=fopen(str,"w");
    
    fprintf(fout,"TITLE=\"Global energy distribution function\"\n");
    fprintf(fout,"VARIABLES=\"Energy\",\"f\"\n");
    
    //collect the sampled information from other processors
    for (thread=1;thread<MS::mesh.nTotalThreads;thread++) 
      for (int iBin=0; iBin<nBin; iBin++) 
	Distribution[iBin]+=pipe.recv<double>(thread);
    
    //normalize the distribution functions
    double norm = 0.0;
    for (int iBin=0; iBin<nBin; iBin++)
      norm += Distribution[iBin] * (EnergyMax-EnergyMin)/nBin;
    
    if (fabs(norm)>0.0) 
      for (int iBin=0; iBin<nBin; iBin++)
	Distribution[iBin]/=norm;
    
    //print the output file
    for (int iBin=0; iBin<nBin; iBin++) {
      double E = EnergyMin + iBin * ((EnergyMax - EnergyMin) / nBin);
      fprintf(fout,"%e ",E);
      fprintf(fout,"  %e\n",Distribution[iBin]);
    }
    
    //close the output file
    fclose(fout);
  }
  else {
    for (int iBin=0; iBin<nBin;iBin++) {
      pipe.send(Distribution[iBin]);
      Distribution[iBin]=0.0;
    }
    
  }
  
  
  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();
  
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}



