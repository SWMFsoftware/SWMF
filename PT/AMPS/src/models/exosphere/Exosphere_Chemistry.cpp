//$Id$


/*
 * Exosphere_Chemistry.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "Exosphere.h"
#include "PhotolyticReactions.h"

//heliocentric distance to the object
double Exosphere::ChemicalModel::HeliocentricDistance=1.0*_AU_;

//the total loss and source rates in the exospheric chemistry
double Exosphere::ChemicalModel::TotalSourceRate[PIC::nTotalSpecies];
double Exosphere::ChemicalModel::TotalLossRate[PIC::nTotalSpecies];

//set the helioncentric distance
void Exosphere::ChemicalModel::Init(double rHeliocentric) {
  HeliocentricDistance=rHeliocentric;

  //init the source and loss sampling buffers
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) TotalSourceRate[spec]=0.0,TotalLossRate[spec]=0.0;
}

//default function for calculation of the lifetile
double Exosphere::ChemicalModel::TotalLifeTime(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  static bool initflag=false;
  static double ReactionLifeTimeTable[PIC::nTotalSpecies];
  static bool ReactionAllowedTable[PIC::nTotalSpecies];

  //init the lifetime table
  if (initflag==false) {
    initflag=true;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      if (PhotolyticReactions::ModelAvailable(s)==true) {
        ReactionAllowedTable[s]=true;
        ReactionLifeTimeTable[s]=1.0/PhotolyticReactions::GetTotalReactionRate(s,HeliocentricDistance);
      }
      else {
        ReactionAllowedTable[s]=false;
        ReactionLifeTimeTable[s]=-1.0;
      }
    }
  }

  //return the lifetime
  ReactionAllowedFlag=ReactionAllowedTable[spec];
  return ReactionLifeTimeTable[spec];
}


//process particle chemical transformation
void Exosphere::ChemicalModel::PhotochemicalModelProcessor(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int *ReactionProductsList,nReactionProducts;
  double *ReactionProductVelocity;
  int ReactionChannel,spec;
  bool PhotolyticReactionRoute;
  PIC::ParticleBuffer::byte *ParticleData;
  double vParent[3],xParent[3],ParentLifeTime;

  //get the particle data
  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  PIC::ParticleBuffer::GetV(vParent,ParticleData);
  PIC::ParticleBuffer::GetX(xParent,ParticleData);

  bool PhotolyticReactionAllowedFlag;

  ParentLifeTime=TotalLifeTime(NULL,spec,ptr,PhotolyticReactionAllowedFlag,node);

  if (PhotolyticReactionAllowedFlag==false) {
    //no reaction occurs -> add the particle to the list
    //the particle is remain in the system
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
    return;
  }

  //init the reaction tables
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
  }

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
  }

  //determine whether the parent particle is removed
  if (rnd()<exp(-ParentTimeStep/ParentLifeTime)) {
    //the particle is remain in the system
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
  else {
    //the particle is removed from the system
    PIC::ParticleBuffer::DeleteParticle(ptr);
  }
}
