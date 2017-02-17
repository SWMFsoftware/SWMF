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
char PIC::CCMC::Parser::ControlFileName[_MAX_STRING_LENGTH_PIC_]="ccmc.InjectionLocation.dat";
char PIC::CCMC::BackgroundDataFileName[_MAX_STRING_LENGTH_PIC_]="amps.Background.data.cdf";

//maximum trajectory integration time
double PIC::CCMC::MaxTrajectoryIntegrationTime=-1.0;

//characteristic speed of traced particles
double *PIC::CCMC::ParticleCharacteristicSpeedTable=NULL;

//computational domain limits
double PIC::CCMC::Domain::xmax[3]={2.070575000000000e+08, 2.994370000000000e+08, 2.994370000000000e+08};
double PIC::CCMC::Domain::xmin[3]={-1.599121000000000e+09,-2.994370000000000e+08,-2.994370000000000e+08};

//resoluton of the computational domain
int PIC::CCMC::Domain::Resolution::mode=-1;
double PIC::CCMC::Domain::Resolution::BackgroundParameterVariationLimit=0.0;
double PIC::CCMC::Domain::Resolution::dxMax=-1.0;
double PIC::CCMC::Domain::Resolution::dxMin=-1.0;
double PIC::CCMC::Domain::Resolution::rXmax=-1.0;
double PIC::CCMC::Domain::Resolution::rXmin=-1.0;

//pocessing particles that cross the internal boundary
int PIC::CCMC::InternalBoundary::ParticleProcessingMode=PIC::CCMC::InternalBoundary::ParticleBoundaryInteractionCode::NoBoundary;
double PIC::CCMC::InternalBoundary::Sphere::Radius=0.0;
double PIC::CCMC::InternalBoundary::Sphere::x0[3]={0.0,0.0,0.0};
int PIC::CCMC::InternalBoundary::Sphere::SamplingMode=_PIC_MODE_OFF_;

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
    else if (strcmp("#MAXTRAJECTORYINTEGRATIONTIME",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::MaxTrajectoryIntegrationTime=atof(str1);
    }
    else if (strcmp("#INTERNALBOUNDARYSPHERE",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("ON",str1)==0) {
        PIC::CCMC::InternalBoundary::ParticleProcessingMode=PIC::CCMC::InternalBoundary::ParticleBoundaryInteractionCode::Sphere;
        Read::InternalBoundarySphere(ifile);
      }
      else if (strcmp("OFF",str1)!=0) exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    else if (strcmp("#BACKGROUNDMODELSPATIALUNITS2METERSFACTOR",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CPLR::DATAFILE::KAMELEON::cdfDataFile2m_ConversionFactor=atof(str1);
    }
    else if (strcmp("#SOURCEREGION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("SPHERE",str1)==0) Read::SourceRegion::Sphere(InjectionBlock,ifile);
      else if (strcmp("CIRCLE",str1)==0) Read::SourceRegion::Circle(InjectionBlock,ifile);
      else if (strcmp("TABLE",str1)==0) Read::SourceRegion::Table(InjectionBlock,ifile);
      else if (strcmp("QUADRILATERAL",str1)==0) Read::SourceRegion::Quadrilateral(InjectionBlock,ifile);
      else exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    else if (strcmp("#VELOCITYDISTRIBUTION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("MAXWELLIAN",str1)==0) Read::VelocityDistribution::Maxwellian(InjectionBlock,ifile);
      else if (strcmp("TABLE",str1)==0) Read::VelocityDistribution::Table(InjectionBlock,ifile);
      else if (strcmp("CONSTANT",str1)==0) Read::VelocityDistribution::Constant(InjectionBlock,ifile);
      else exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    else if (strcmp("#DOMAIN",str1)==0) {
      Read::DomainLimits(ifile);
    }
    else if (strcmp("#CHARACTERISTICSPEED",str1)==0) {
      PIC::CCMC::Parser::Read::CharacteristicSpeedTable(ifile);
    }

    else if (strcmp("",str1)==0) {
      //do nothing -> move to the next line
    }
    else exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //close the file
  ifile.closefile();
}

//Read the internal boundary sphere section
void PIC::CCMC::Parser::Read::InternalBoundarySphere(CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("RADIUS",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::InternalBoundary::Sphere::Radius=atof(str1);
    }
    else if (strcmp("SAMPLING",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("ON",str1)==0) PIC::CCMC::InternalBoundary::Sphere::SamplingMode=_PIC_MODE_ON_;
      else if (strcmp("OFF",str1)==0) PIC::CCMC::InternalBoundary::Sphere::SamplingMode=_PIC_MODE_OFF_;
      else exit(__LINE__,__FILE__,"Error: the option is not recognized");
    }
    else if (strcmp("X0",str1)==0) {
       for (int idim=0;idim<3;idim++) {
         ifile.CutInputStr(str1,str);
         PIC::CCMC::InternalBoundary::Sphere::x0[idim]=atof(str1);
       }
    }
    else if (strcmp("#ENDINTERNALBOUNDARYSPHERE",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}

//Read individual section of the input file
void PIC::CCMC::Parser::Read::DomainLimits(CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("XMIN",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        PIC::CCMC::Domain::xmin[i]=atof(str1);
      }
    }
    else if (strcmp("XMAX",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        PIC::CCMC::Domain::xmax[i]=atof(str1);
      }
    }
    else if (strcmp("DXMAX",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::Domain::Resolution::dxMax=atof(str1);
    }
    else if (strcmp("DXMIN",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::Domain::Resolution::dxMin=atof(str1);
    }
    else if (strcmp("RXMIN",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::Domain::Resolution::rXmin=atof(str1);
    }
    else if (strcmp("RXMAX",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::Domain::Resolution::rXmax=atof(str1);
    }
    else if (strcmp("LOCALRESOLUTIONMODE",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("CONSTANT",str1)==0) {
        PIC::CCMC::Domain::Resolution::mode=PIC::CCMC::Domain::Resolution::TYPE::Constant;
      }
      else if (strcmp("BACKGROUNDFIELDVARIATION",str1)==0) {
        PIC::CCMC::Domain::Resolution::mode=PIC::CCMC::Domain::Resolution::TYPE::BackgroundFieldVariation;
      }
      else if (strcmp("LOGARITHMIC",str1)==0) {
        PIC::CCMC::Domain::Resolution::mode=PIC::CCMC::Domain::Resolution::TYPE::Logarithmic;
      }
      else exit(__LINE__,__FILE__,"Error: the option is not recognized");
    }
    else if (strcmp("BACKGROUNDFIELDVARIATIONLIMIT",str1)==0) {
      ifile.CutInputStr(str1,str);
      PIC::CCMC::Domain::Resolution::BackgroundParameterVariationLimit=atof(str1);
    }
    else if (strcmp("#ENDDOMAIN",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}

//calculate the requested local resolution
double PIC::CCMC::Domain::Resolution::GetlocalResolution(double *x) {
  double res=0.0;

  switch (mode) {
  case TYPE::Constant:
    res=0.5*(dxMin+dxMax);
    break;

  case TYPE::Logarithmic:
    res=log(dxMin)+(log(dxMax)-log(dxMin))/(rXmax-rXmin)*(sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2))-rXmin);
    if (res<dxMin) res=dxMin;
    if (res>dxMax) res=dxMax;
    break;

  case TYPE::BackgroundFieldVariation:
    exit(__LINE__,__FILE__,"Error: not implemented");
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the option is not found");
  }

  return res;
}

//return the characteristic particle speed
void PIC::CCMC::Parser::Read::CharacteristicSpeedTable(CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;
  int spec;

  if (PIC::CCMC::ParticleCharacteristicSpeedTable==NULL) {
    PIC::CCMC::ParticleCharacteristicSpeedTable=new double [PIC::nTotalSpecies];
    for (spec=0;spec<PIC::nTotalSpecies;spec++) PIC::CCMC::ParticleCharacteristicSpeedTable[spec]=-1.0;
  }

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("V",str1)==0) {
      ifile.CutInputStr(str1,str);
      spec=PIC::MolecularData::GetSpecieNumber(str1);

      if (spec==-1) {
        char msg[100];

        sprintf(msg,"Error: species %s is not defined in the input file",str1);
        exit(__LINE__,__FILE__,msg);
      }

      ifile.CutInputStr(str1,str);
      PIC::CCMC::ParticleCharacteristicSpeedTable[spec]=atof(str1);
    }

    else if (strcmp("#ENDCHARACTERISTICSPEED",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
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

void PIC::CCMC::Parser::Read::SourceRegion::Circle(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  InjectionBlock.SpatialDistribution.Type=PIC::CCMC::DEF::SOURCE::TYPE::Circle;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("CENTER",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Circle.Origin[i]=atof(str1);
      }
    }
    else if (strcmp("RADIUS",str1)==0) {
      ifile.CutInputStr(str1,str);
      InjectionBlock.SpatialDistribution.Circle.Radius=atof(str1);
    }
    else if (strcmp("NORMAL",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Circle.Normal[i]=atof(str1);
      }

      //normalize the vector
      double l=0.0;
      int i;

      for (i=0;i<3;i++) l+=pow(InjectionBlock.SpatialDistribution.Circle.Normal[i],2);
      for (l=sqrt(l),i=0;i<3;i++) InjectionBlock.SpatialDistribution.Circle.Normal[i]/=l;

      //determine the coordinate frame related to the sphere
      Vector3D::GetNormFrame(InjectionBlock.SpatialDistribution.Circle.e0,InjectionBlock.SpatialDistribution.Circle.e1,InjectionBlock.SpatialDistribution.Circle.Normal);
    }
    else if (strcmp("SPATIALDISTRIBUTION",str1)==0) {
      ifile.CutInputStr(str1,str);

      if (strcmp("UNIFORM",str1)==0) InjectionBlock.SpatialDistribution.Circle.SpatialDistributionType=PIC::CCMC::DEF::SOURCE::SHPERE::TYPE::Uniform;
      else if (strcmp("UNIFORM",str1)==0) {
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
      else exit(__LINE__,__FILE__,"Error: the option is not found");
    }
    else if (strcmp("#ENDSOURCEREGION",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}

void PIC::CCMC::Parser::Read::SourceRegion::Table(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;
  int i,idim;
  PIC::CCMC::ParticleInjection::cInjectionRegionTable c;

  for (i=0;(i<InjectionBlock.nTestParticles)&&(ifile.eof()==false);i++) {
    ifile.GetInputStr(str,sizeof(str));

    for (idim=0;idim<3;idim++) {
      ifile.CutInputStr(str1,str);
      c.x[idim]=atof(str1);
    }

    InjectionBlock.SpatialDistribution.Table.push_back(c);
  }


  InjectionBlock.SpatialDistribution.Type=PIC::CCMC::DEF::SOURCE::TYPE::Table;

  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);

  if (strcmp("#ENDSOURCEREGION",str1)!=0) exit(__LINE__,__FILE__,"Error: the option is not found");
}

void PIC::CCMC::Parser::Read::SourceRegion::Quadrilateral(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  InjectionBlock.SpatialDistribution.Type=PIC::CCMC::DEF::SOURCE::TYPE::Quadrilateral;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("XCENTER",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Quadrilateral.xCenter[i]=atof(str1);
      }
    }
    else if (strcmp("DX0",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Quadrilateral.dX0[i]=atof(str1);
      }
    }
    else if (strcmp("DX1",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.SpatialDistribution.Quadrilateral.dX1[i]=atof(str1);
      }
    }
    else if (strcmp("#ENDSOURCEREGION",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}

void PIC::CCMC::Parser::Read::VelocityDistribution::Table(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;
  int i,idim;
  PIC::CCMC::ParticleInjection::cVelocityDistributionTable c;

  for (i=0;(i<InjectionBlock.nTestParticles)&&(ifile.eof()==false);i++) {
    ifile.GetInputStr(str,sizeof(str));

    for (idim=0;idim<3;idim++) {
      ifile.CutInputStr(str1,str);
      c.v[idim]=atof(str1);
    }

    InjectionBlock.VelocityDistribution.Table.push_back(c);
  }

  InjectionBlock.VelocityDistribution.Type=PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Table;

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

void PIC::CCMC::Parser::Read::VelocityDistribution::Constant(PIC::CCMC::ParticleInjection::cInjectionDescriptor& InjectionBlock,CiFileOperations& ifile) {
  char str1[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],*endptr;

  InjectionBlock.VelocityDistribution.Type=PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Constant;

  while (ifile.eof()==false) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(str1,str);

    if (strcmp("VELOCITY",str1)==0) {
      for (int i=0;i<3;i++) {
        ifile.CutInputStr(str1,str);
        InjectionBlock.VelocityDistribution.Constant.v[i]=atof(str1);
      }
    }
    else if (strcmp("#ENDVELOCITYDISTRIBUTION",str1)==0) return;
    else exit(__LINE__,__FILE__,"Error: the option is not recognized");
  }
}


//load particles that will be tracked in the simulation
void PIC::CCMC::LoadParticles() {
  int np,spec,idim;
  int iInjectionEntry,iTriangle;
  double x[3],v[3];
  double r,phi,sinTheta,cosTheta;

  //variables that are used to generate new particle locations of a Quadrilateral
  double xLocal[2];
  PIC::CCMC::ParticleInjection::cInjectionRegionQuadrilateral* Quadrilateral;
  static const int e0SignTable[4]={1,-1,-1,1},e1SignTable[4]={1,1,-1,-1};

  PIC::ParticleBuffer::byte *newParticleData;
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iInjectionEntry=0;iInjectionEntry<ParticleInjection::InjectionDescriptorList.size();iInjectionEntry++) {
    for (np=0;np<ParticleInjection::InjectionDescriptorList[iInjectionEntry].nTestParticles;np++) {
      static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* t=NULL;

      //generate the location
      switch (ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Type) {
      case PIC::CCMC::DEF::SOURCE::TYPE::Table:
        for (idim=0;idim<3;idim++) x[idim]=ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Table[np].x[idim];
        break;

      case PIC::CCMC::DEF::SOURCE::TYPE::Sphere:

        //verify the the origin of the sphere is inside the domain
        for (idim=0;idim<3;idim++) {
          if ((ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[idim]<PIC::Mesh::mesh.rootTree->xmin[idim])  ||
          (ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[idim]>PIC::Mesh::mesh.rootTree->xmax[idim])) { 
            //the origin of the sphere is outside of thedomain => terminate execution

            char message[1000];
            double *x=ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin;
            double *xmin=PIC::Mesh::mesh.rootTree->xmin;
            double *xmax=PIC::Mesh::mesh.rootTree->xmax; 

            sprintf(message," Error: injection sphere is outside of the domain. Redefine the location of the injection spherical region.\nThe origin of the sphere: %e, %e, %e.\nThe domain limits: xmin= %e, %e, %e, and xmax=%e, %e, %e\n",x[0],x[1],x[2], xmin[0],xmin[1],xmin[2], xmax[0],xmax[1],xmax[2]); 

            exit(__LINE__,__FILE__,message);
          } 
        }

        do {
          r=ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Radius*pow(rnd(),1.0/3.0);
          phi=2.0*Pi*rnd();
          sinTheta=2.0*rnd()-1;
          cosTheta=sqrt(1.0-sinTheta*sinTheta);

          x[0]=r*cos(phi)*cosTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[0];
          x[1]=r*sin(phi)*cosTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[1];
          x[2]=r*sinTheta+ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Spherical.Origin[2];
        }
        while ((t=PIC::Mesh::mesh.findTreeNode(x,t))==NULL);

        break;
      case PIC::CCMC::DEF::SOURCE::TYPE::Circle:
        //verify the the origin of the sphere is inside the domain
        for (idim=0;idim<3;idim++) {
          if ((ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Circle.Origin[idim]<PIC::Mesh::mesh.rootTree->xmin[idim])  ||
          (ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Circle.Origin[idim]>PIC::Mesh::mesh.rootTree->xmax[idim])) {
            exit(__LINE__,__FILE__,"Error: injection Circle is outside of the domain. Redefine the location of the injection spherical region");
          }
        }

        do {
          PIC::CCMC::ParticleInjection::cInjectionRegionCircle *Circle=&ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Circle;

          r=Circle->Radius*pow(rnd(),1.0/2.0);
          phi=2.0*Pi*rnd();

          for (idim=0;idim<3;idim++) x[idim]=Circle->Origin[idim]+r*(sin(phi)*Circle->e0[idim]+cos(phi)*Circle->e1[idim]);
        }
        while ((t=PIC::Mesh::mesh.findTreeNode(x,t))==NULL);

        break;
      case PIC::CCMC::DEF::SOURCE::TYPE::Quadrilateral:
        //verify that at least one corner point of the Quadrilateral is within the domain
        {
          bool InsideDomainFlag=false;

          Quadrilateral=&ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Quadrilateral;

          for (iTriangle=0;(iTriangle<4)&&(InsideDomainFlag==false);iTriangle++) for (int np=0;np<3;np++) {
            InsideDomainFlag=true;

            switch (np) {
            case 0:
              xLocal[0]=0.0,xLocal[1]=0.0;
              break;
            case 1:
              xLocal[0]=1.0,xLocal[1]=0.0;
              break;
            case 2:
              xLocal[0]=0.0,xLocal[1]=1.0;
              break;
            }

            for (idim=0;idim<3;idim++) {
              x[idim]=Quadrilateral->xCenter[idim]+
                xLocal[0]*Quadrilateral->dX0[idim]*e0SignTable[iTriangle]+
                xLocal[1]*Quadrilateral->dX1[idim]*e1SignTable[iTriangle];

              if ((x[idim]<PIC::Mesh::mesh.rootTree->xmin[idim])||(x[idim]>PIC::Mesh::mesh.rootTree->xmax[idim])) {
                //the point is outside of the domain
                InsideDomainFlag=false;
                break;
              }
            }
          }

          if (InsideDomainFlag==false) exit(__LINE__,__FILE__,"Error: all points of the Quadrilateral are outside of the domain");
        }


        //determine the local coordinates of the injection point in a triangle
        do {
          xLocal[0]=1.0-sqrt(rnd());
          xLocal[1]=rnd()*(1.0-xLocal[0]);

          //determine the triangle at which the point will be generated
          iTriangle=(int)(4.0*rnd());
          Quadrilateral=&ParticleInjection::InjectionDescriptorList[iInjectionEntry].SpatialDistribution.Quadrilateral;

          for (idim=0;idim<3;idim++) x[idim]=Quadrilateral->xCenter[idim]+
            xLocal[0]*Quadrilateral->dX0[idim]*e0SignTable[iTriangle]+
            xLocal[1]*Quadrilateral->dX1[idim]*e1SignTable[iTriangle];
        }
        while ((t=PIC::Mesh::mesh.findTreeNode(x,t))==NULL);

        break;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }

      //generate the velocity vector
      switch (ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Type) {
      case PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Table:
        for (idim=0;idim<3;idim++) v[idim]=ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Table[np].v[idim];
        break;

      case PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Maxwellian:
        PIC::Distribution::MaxwellianVelocityDistribution(v,
            ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Maxwellian.BulkVelocity,
            ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Maxwellian.Temeprature,
            spec);
        break;

      case PIC::CCMC::DEF::VELOCITY_DISTRIBUTION::TYPE::Constant:
        for (idim=0;idim<3;idim++) v[idim]=ParticleInjection::InjectionDescriptorList[iInjectionEntry].VelocityDistribution.Constant.v[idim];
        break;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }

      //the new particle must be outside of the internal surface
      if (InternalBoundary::ParticleProcessingMode==InternalBoundary::ParticleBoundaryInteractionCode::Sphere) {
        double r=0.0;

        for (int idim=0;idim<3;idim++) r+=pow(x[idim]-InternalBoundary::Sphere::x0[idim],2);
        if (r<pow(InternalBoundary::Sphere::Radius,2)) continue;
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
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)node);
      #endif

      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,0.0,node,true);
    }
  }
}

//the main tracking procedure
int PIC::CCMC::TraceParticles() {
  long int nTotalParticles;
  char fname[_MAX_STRING_LENGTH_PIC_];
  double TimeCounter[PIC::nTotalSpecies];
  int spec;

  //load particles
  LoadParticles();

  //init the time counter
  for (spec=0;spec<PIC::nTotalSpecies;spec++) TimeCounter[spec]=0.0;

  //continue simulation untill particles are in the system
  MPI_Allreduce(&PIC::ParticleBuffer::NAllPart,&nTotalParticles,PIC::nTotalSpecies,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) printf("$PREFIX: Particle Tracking: The total number of tracked particle is %ld\n",nTotalParticles);

  int niter=0;
  const int nOutputStep=1000;

  while (nTotalParticles!=0) {
    PIC::Mover::MoveParticles();
    PIC::Parallel::ExchangeParticleData();

    niter++;

    if (niter%nOutputStep==0) {
      sprintf(fname,"%s/amps.TrajectoryTracking.out=%i",OutputDataFileDirectory,niter/nOutputStep);
      PIC::ParticleTracker::OutputTrajectory(fname);
    }

    MPI_Allreduce(&PIC::ParticleBuffer::NAllPart,&nTotalParticles,PIC::nTotalSpecies,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    //update the time counter
    PIC::SimulationTime::Update();

    //compare the TimeCounter with the requested max trajectory integraion time and quit if the time is up
    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
      if ((MaxTrajectoryIntegrationTime>0.0) && (MaxTrajectoryIntegrationTime<TimeCounter[spec])) {
        //the integration time for species 'spec' exseeds that of the requested -> terminate trajectory integration for species 'spec'
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
        PIC::Mesh::cDataBlockAMR *block;
        int i,j,k;
        long int ptr,ptrNext;

        for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
          block=node->block;

          if (block==NULL) continue;

          for (k=0;k<_BLOCK_CELLS_Z_;k++) {
            for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
              for (i=0;i<_BLOCK_CELLS_X_;i++) {
                ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

                while (ptr!=-1) {
                  ptrNext=PIC::ParticleBuffer::GetNext(ptr);
                  if (PIC::ParticleBuffer::GetI(ptr)==spec) PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
                  else {
                    //delete the particle when: 1) it is within a boundary sphere, and 2) sampling of the partiucle on the sphere if turned off
                    //(if the flux sampling is on then the particel is deleted during partucle motion step)
                    if ((InternalBoundary::ParticleProcessingMode==InternalBoundary::ParticleBoundaryInteractionCode::Sphere)&&(InternalBoundary::Sphere::SamplingMode==_PIC_MODE_OFF_)) {
                      double r=0.0,*x=PIC::ParticleBuffer::GetX(ptr);

                      for (int idim=0;idim<3;idim++) r+=pow(x[idim]-InternalBoundary::Sphere::x0[idim],2);

                      if (r<pow(InternalBoundary::Sphere::Radius,2)) {
                        PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
                      }
                    }
                  }

                  ptr=ptrNext;
                }
              }
            }
          }
        }
      }
      else TimeCounter[spec]+=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
    }
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

//================================================
//return characteristic particle speed
double PIC::CCMC::GetParticleCharacteristicSpeed(int spec) {
  if (ParticleCharacteristicSpeedTable==NULL) exit(__LINE__,__FILE__,"Error: PIC::CCMC::ParticleCharacteristicSpeedTable is not initialized");

  if (ParticleCharacteristicSpeedTable[spec]<0.0) {
    char msg[200];

    sprintf(msg,"Characteristic speed for species %s is not defined",PIC::MolecularData::GetChemSymbol(spec));
    exit(__LINE__,__FILE__,msg);
  }

  return ParticleCharacteristicSpeedTable[spec];
}
