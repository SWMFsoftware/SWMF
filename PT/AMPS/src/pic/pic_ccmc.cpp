//$Id$


/*
 * pic_ccmc.cpp
 *
 *  Created on: Sep 28, 2015
 *      Author: vtenishe
 */
//the procedures for tracking model particles WITHOUT initialization of the weight, sampleing, ets.
//The functions use the standard particle integration procedure -> ALL particle transformation reactions MUST be turned OFF


#include "pic.h"


vector<PIC::CCMC::ParticleInjection::cInjectionDescriptor> PIC::CCMC::ParticleInjection::InjectionDescriptorList;
char PIC::CCMC::Parser::ControlFileName[_MAX_STRING_LENGTH_PIC_]="ccmc.dat";


//read the control file
void PIC::CCMC::Parser::LoadControlFile() {
  CiFileOperations ifile;
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr,fname[_MAX_STRING_LENGTH_PIC_];
  PIC::CCMC::ParticleInjection::cInjectionDescriptor InjectionBlock;

  sprintf(fname,"%s/%s",PIC::UserModelInputDataPath,ControlFileName);

  if (PIC::ThisThread==0) printf("$PREFIX:Trajectory Tracking Control File: %s\n",fname);

  if (access(fname,R_OK)!=0) {
  printf("Cannot find the input file:%s\n",fname);
  exit(__LINE__,__FILE__);
  }

  //read the file
  ifile.openfile(fname);

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);


    if (strcmp("#TRACKING",str1)==0) {
      ifile.CutInputStr(str1,str);
      InjectionBlock.nTestParticles=strtol(str1,&endptr,10);
    }
    else if (strcmp("#ENDTRACKING",str1)==0) {
      PIC::CCMC::ParticleInjection::InjectionDescriptorList.push_back(InjectionBlock);
    }
    else if (strcmp("STARTTIME",str1)==0) {
      ifile.CutInputStr(str1,str);
      InjectionBlock.StartTime=atof(str1);
    }

    else if (strcmp("#SOURCEREGION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("SPHERE",str1)==0) Read::SourceRegion::Sphere(InjectionBlock,ifile);
      else if (strcmp("CONSTANT",str1)==0) Read::SourceRegion::Constant(InjectionBlock,ifile);
      else exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    else if (strcmp("#VELOCITYDISTRIBUTION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("MAXWELLIAN",str1)==0) Read::VelocityDistribution::Maxwellian(InjectionBlock,ifile);
      else if (strcmp("CONSTANT",str1)==0) Read::VelocityDistribution::Constant(InjectionBlock,ifile);
      else exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    else if (strcmp("",str1)==0) {
      //do nothing -> move to the next line
    }
    else exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //close the file
  ifile.closefile();
}

void PIC::CCMC::Parser::Read::SourceRegion::Sphere(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  InjectionBlock.SpatialDistribution.Type=PIC::CCMC::DEF::SOURCE::TYPE::Sphere;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("CENTER",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Spherical.Origin[i]=atof(str1);
      }
    }
    else if (strcmp("RADIUS",str1)==0) {
      ifile.CutInputStr(str1,str);
      InjectionBlock.SpatialDistribution.Spherical.Radius=atof(str1);
    }
    else if (strcmp("SPATIALDISTRIBUTION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("UNIFORM",str1)==0) InjectionBlock.SpatialDistribution.Spherical.SpatialDistributionType=PIC::CCMC::DEF::SOURCE::SHPERE::TYPE::Uniform;
      else if (strcmp("UNIFORM",str1)==0) {
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
      else exit(__LINE__,__FILE__,"Error: the option is not found");
    }
    else if (strcmp("#ENDSOURCEREGION",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}

void PIC::CCMC::Parser::Read::SourceRegion::Constant(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;
  int i,idim;
  PIC::CCMC::ParticleInjection::cInjectionRegionConstant c;

  for (i=0;(i<InjectionBlock.nTestParticles)&&(ifile.eof()==false);i++) {
    ifile.GetInputStr(str,sizeof(str));

    for (idim=0;idim<3;idim++) {
      ifile.CutInputStr(str1,str);
      c.x[idim]=atof(str1);
    }

    InjectionBlock.SpatialDistribution.Constant.push_back(c);
  }


  InjectionBlock.SpatialDistribution.Type=PIC::CCMC::DEF::SOURCE::TYPE::Constant;

  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);

  if (strcmp("#ENDSOURCEREGION",str1)!=0) exit(__LINE__,__FILE__,"Error: the option is not found");
}

void PIC::CCMC::Parser::Read::VelocityDistribution::Constant(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;
  int i,idim;
  PIC::CCMC::ParticleInjection::cVelocityDistributionConstant c;

  for (i=0;(i<InjectionBlock.nTestParticles)&&(ifile.eof()==false);i++) {
    ifile.GetInputStr(str,sizeof(str));

    for (idim=0;idim<3;idim++) {
      ifile.CutInputStr(str1,str);
      c.v[idim]=atof(str1);
    }

    InjectionBlock.VelocityDistribution.Constant.push_back(c);
  }

  InjectionBlock.VelocityDistribution.Type=PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Constant;

  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);

  if (strcmp("#ENDVELOCITYDISTRIBUTION",str1)!=0) exit(__LINE__,__FILE__,"Error: the option is not found");
}

void PIC::CCMC::Parser::Read::VelocityDistribution::Maxwellian(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  InjectionBlock.VelocityDistribution.Type=PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Maxwellian;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("BULKVELOCITY",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.VelocityDistribution.Maxwellian.BulkVelocity[i]=atof(str1);
      }
    }
    else if (strcmp("TEMPERATURE",str1)==0) {
      ifile.CutInputStr(str1,str);
      InjectionBlock.VelocityDistribution.Maxwellian.Temeprature=atof(str1);
    }
    else if (strcmp("#ENDVELOCITYDISTRIBUTION",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}


//load particles that will be tracked in the simulation
void PIC::CCMC::LoadParticles() {
  int np,spec,idim;
  int iInjectionEntry;
  double x[3],v[3];
  double r,phi,sinTheta,cosTheta;

  PIC::ParticleBuffer::byte *newParticleData;
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iInjectionEntry=0;iInjectionEntry<ParticleInjection::InjectionDescriptorList.size();iInjectionEntry++) {
    for (np=0;np<ParticleInjection::InjectionDescriptorList[iInjectionEntry].nTestParticles;np++) {

      //generate the location
      switch (ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Type) {
      case PIC::CCMC::DEF::SOURCE::TYPE::Constant:
        for (idim=0;idim<3;idim++) x[idim]=ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Constant[np].x[idim];
        break;

      case PIC::CCMC::DEF::SOURCE::TYPE::Sphere:
        r=ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Radius*pow(rnd(),1.0/3.0);
        phi=2.0*Pi*rnd();
        sinTheta=2.0*rnd()-1;
        cosTheta=sqrt(1.0-sinTheta*sinTheta);

        x[0]=r*cos(phi)*cosTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[0];
        x[1]=r*sin(phi)*cosTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[1];
        x[2]=r*sinTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[2];
        break;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }

      //generate the velocity vector
      switch (ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Type) {
      case PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Constant:
        for (idim=0;idim<3;idim++) v[idim]=ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Constant[np].v[idim];
        break;

      case PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Maxwellian:
        PIC::Distribution::MaxwellianVelocityDistribution(v,
            ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Maxwellian.BulkVelocity,
            ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Maxwellian.Temeprature,
            spec);
        break;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }

      //determine the block and processor number for the particle
      static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

      node=PIC::Mesh::mesh.findTreeNode(x,node);
      if (node->Thread!=PIC::Mesh::mesh.ThisThread) continue;

      //create and inject the particle into the system
      long int newParticle;

      PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
      PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
      PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

      #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,(PIC::ParticleBuffer::byte*)tempParticleData);
      #endif

      newParticle=PIC::ParticleBuffer::GetNewParticle();
      newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
      memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

      //apply condition of tracking the particle
      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      PIC::ParticleTracker::InitParticleID(newParticleData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData);
      #endif

      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,0.0,node,true);
    }
  }
}

//the main tracking procedure
int PIC::CCMC::TraceParticles() {
  long int nTotalParticles;
  char fname[_MAX_STRING_LENGTH_PIC_];

  //load particles
  LoadParticles();

  //continue simulation untill particles are in the system
  MPI_Allreduce(&PIC::ParticleBuffer::NAllPart,&nTotalParticles,PIC::nTotalSpecies,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) printf("$PREFIX: Particle Tracking: The total number of tacked particle is %i\n",nTotalParticles);

  int niter=0;
  const int nOutputStep=1000;

  while (nTotalParticles!=0) {
    PIC::Mover::MoveParticles();
    PIC::Parallel::ExchangeParticleData();

    niter++;

    if (niter%nOutputStep==0) {
      sprintf(fname,"%s/amps.TrajectoryTracking.out=%ld",OutputDataFileDirectory,nOutputStep/nOutputStep);
      PIC::ParticleTracker::OutputTrajectory(fname);
    }

    MPI_Allreduce(&PIC::ParticleBuffer::NAllPart,&nTotalParticles,PIC::nTotalSpecies,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }

  //output sampled trajectories
  sprintf(fname,"%s/amps.TrajectoryTracking",OutputDataFileDirectory);
  PIC::ParticleTracker::OutputTrajectory(fname);

  //combine all trajectory files into a single reference file
  if ((_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_)&&(PIC::ThisThread==0)) {
    char cmd[_MAX_STRING_LENGTH_PIC_];

    sprintf(cmd,"cat %s/amps.TrajectoryTracking.s=*.dat > %s/test_CCMC-Individual_Trajectories.dat",OutputDataFileDirectory,OutputDataFileDirectory);
    system(cmd);
  }


  return _PIC_TIMESTEP_RETURN_CODE__SUCCESS_;
}
