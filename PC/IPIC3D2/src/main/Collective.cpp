
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "input_array.h"
#include "ipichdf5.h"
#include "Collective.h"
#include "ConfigFile.h"
#include "limits.h" // for INT_MAX
#include "MPIdata.h"
#include "debug.h"
#include "asserts.h" // for assert_ge
#include "string.h"

// order must agree with Enum in Collective.h
static const char *enumNames[] =
{
  "default",
  "initial",
  "final",
  // used by ImplSusceptMode
  "explPredict",
  "implPredict",
  // marker for last enumerated symbol of this class
  "NUMBER_OF_ENUMS",
  "INVALID_ENUM"
};

int Collective::read_enum_parameter(const char* option_name, const char* default_value,
  const ConfigFile& config)
{
  string enum_name = config.read < string >(option_name,default_value);
  // search the list (could use std::map)
  //
  for(int i=0;i<NUMBER_OF_ENUMS;i++)
  {
    if(!strcmp(enum_name.c_str(),enumNames[i]))
      return i;
  }
  // could not find enum, so issue error and quit.
  if(!MPIdata::get_rank())
  {
    eprintf("in input file %s there is an invalid option %s\n",
      inputfile.c_str(), enum_name.c_str());
  }
  MPIdata::exit(1);
  // this is a better way
  return INVALID_ENUM;
}

const char* Collective::get_name_of_enum(int in)
{
  assert_ge(in, 0);
  assert_lt(in, NUMBER_OF_ENUMS);
  return enumNames[in];
}

/*! Read the input file from text file and put the data in a collective wrapper: if it's a restart read from input file basic sim data and load particles and EM field from restart file */
void Collective::ReadInput(string inputfile) {
  using namespace std;
  int test_verbose;
  // Loading the input file 
  ConfigFile config(inputfile);
  // the following variables are ALWAYS taken from inputfile, even if restarting 
  {

#ifdef BATSRUS
    if(RESTART1)
    {
      cout<<" The fluid interface can not handle RESTART yet, aborting!\n"<<flush;
      abort();
    }
#endif

    dt = config.read < double >("dt");
    ncycles = config.read < int >("ncycles");
    th = config.read < double >("th",1.0);

    Smooth = config.read < double >("Smooth",1.0);
    SmoothNiter = config.read < int >("SmoothNiter",6);

    SaveDirName = config.read < string > ("SaveDirName","data");
    RestartDirName = config.read < string > ("RestartDirName","data");
    ns = config.read < int >("ns");
    nstestpart = config.read < int >("nsTestPart", 0);
    NpMaxNpRatio = config.read < double >("NpMaxNpRatio",1.5);
    assert_ge(NpMaxNpRatio, 1.);
    // mode parameters for second order in time
    PushWithBatTime = config.read < double >("PushWithBatTime",0);
    PushWithEatTime = config.read < double >("PushWithEatTime",1);
    ImplSusceptTime = config.read < double >("ImplSusceptTime",0);
    ImplSusceptMode = read_enum_parameter("ImplSusceptMode", "initial",config);
    switch(ImplSusceptMode)
    {
      // values not yet supported:
      case explPredict:
      case implPredict:
      default:
        unsupported_value_error(ImplSusceptMode);
      // supported values:
      case initial:
        ;
    }
    // GEM Challenge 
    B0x = config.read <double>("B0x",0.0);
    B0y = config.read <double>("B0y",0.0);
    B0z = config.read <double>("B0z",0.0);

    // Earth parameters
    B1x = 0.0;
    B1y = 0.0;
    B1z = 0.0;
    B1x = config.read <double>("B1x",0.0);
    B1y = config.read <double>("B1y",0.0);
    B1z = config.read <double>("B1z",0.0);

    delta = config.read < double >("delta",0.5);

    Case              = config.read<string>("Case");
    wmethod           = config.read<string>("WriteMethod");
    SimName           = config.read<string>("SimulationName");
    PoissonCorrection = config.read<string>("PoissonCorrection");

    rhoINIT = new double[ns];
    array_double rhoINIT0 = config.read < array_double > ("rhoINIT");
    rhoINIT[0] = rhoINIT0.a;
    if (ns > 1)
      rhoINIT[1] = rhoINIT0.b;
    if (ns > 2)
      rhoINIT[2] = rhoINIT0.c;
    if (ns > 3)
      rhoINIT[3] = rhoINIT0.d;
    if (ns > 4)
      rhoINIT[4] = rhoINIT0.e;
    if (ns > 5)
      rhoINIT[5] = rhoINIT0.f;

    rhoINJECT = new double[ns];
    array_double rhoINJECT0 = config.read<array_double>( "rhoINJECT" );
    rhoINJECT[0]=rhoINJECT0.a;
    if (ns > 1)
      rhoINJECT[1]=rhoINJECT0.b;
    if (ns > 2)
      rhoINJECT[2]=rhoINJECT0.c;
    if (ns > 3)
      rhoINJECT[3]=rhoINJECT0.d;
    if (ns > 4)
      rhoINJECT[4]=rhoINJECT0.e;
    if (ns > 5)
      rhoINJECT[5]=rhoINJECT0.f;

    // take the tolerance of the solvers
    CGtol = config.read < double >("CGtol",1e-3);
    GMREStol = config.read < double >("GMREStol",1e-3);
    NiterMover = config.read < int >("NiterMover",3);
    // take the injection of the particless
    Vinj = config.read < double >("Vinj",0.0);

    // take the output cycles
    FieldOutputCycle = config.read < int >("FieldOutputCycle",100);
    ParticlesOutputCycle = config.read < int >("ParticlesOutputCycle",0);
    FieldOutputTag     =   config.read <string>("FieldOutputTag","");
    ParticlesOutputTag =   config.read <string>("ParticlesOutputTag","");
    MomentsOutputTag   =   config.read <string>("MomentsOutputTag","");
    TestParticlesOutputCycle = config.read < int >("TestPartOutputCycle",0);
    testPartFlushCycle = config.read < int >("TestParticlesOutputCycle",10);
    RestartOutputCycle = config.read < int >("RestartOutputCycle",5000);
    DiagnosticsOutputCycle = config.read < int >("DiagnosticsOutputCycle", FieldOutputCycle);
    CallFinalize = config.read < bool >("CallFinalize", true);
  }

  //read everything from input file, if restart is true, overwrite the setting - bug fixing

  restart_status = 0;
  last_cycle = -1;
  c = config.read < double >("c",1.0);

#ifdef BATSRUS
  // set grid size and resolution based on the initial file from fluid code
  Lx =  getFluidLx();
  Ly =  getFluidLy();
  Lz =  getFluidLz();
  nxc = getFluidNxc();
  nyc = getFluidNyc();
  nzc = getFluidNzc();
#else
  Lx = config.read < double >("Lx",10.0);
  Ly = config.read < double >("Ly",10.0);
  Lz = config.read < double >("Lz",10.0);
  nxc = config.read < int >("nxc",64);
  nyc = config.read < int >("nyc",64);
  nzc = config.read < int >("nzc",64);
#endif
  XLEN = config.read < int >("XLEN",1);
  YLEN = config.read < int >("YLEN",1);
  ZLEN = config.read < int >("ZLEN",1);
  PERIODICX = config.read < bool >("PERIODICX",true);
  PERIODICY = config.read < bool >("PERIODICY",true);
  PERIODICZ = config.read < bool >("PERIODICZ",true);

  PERIODICX_P = config.read < bool >("PERIODICX_P",PERIODICX);
  PERIODICY_P = config.read < bool >("PERIODICY_P",PERIODICY);
  PERIODICZ_P = config.read < bool >("PERIODICZ_P",PERIODICZ);

  x_center = config.read < double >("x_center",5.0);
  y_center = config.read < double >("y_center",5.0);
  z_center = config.read < double >("z_center",5.0);
  L_square = config.read < double >("L_square",5.0);


  uth = new double[ns];
  vth = new double[ns];
  wth = new double[ns];
  u0 = new double[ns];
  v0 = new double[ns];
  w0 = new double[ns];

  array_double uth0 = config.read < array_double > ("uth");
  array_double vth0 = config.read < array_double > ("vth");
  array_double wth0 = config.read < array_double > ("wth");
  array_double u00 = config.read < array_double > ("u0");
  array_double v00 = config.read < array_double > ("v0");
  array_double w00 = config.read < array_double > ("w0");

  uth[0] = uth0.a;
  vth[0] = vth0.a;
  wth[0] = wth0.a;
  u0[0] = u00.a;
  v0[0] = v00.a;
  w0[0] = w00.a;
  if (ns > 1) {
    uth[1] = uth0.b;
    vth[1] = vth0.b;
    wth[1] = wth0.b;
    u0[1] = u00.b;
    v0[1] = v00.b;
    w0[1] = w00.b;
  }
  if (ns > 2) {
    uth[2] = uth0.c;
    vth[2] = vth0.c;
    wth[2] = wth0.c;
    u0[2] = u00.c;
    v0[2] = v00.c;
    w0[2] = w00.c;
  }
  if (ns > 3) {
    uth[3] = uth0.d;
    vth[3] = vth0.d;
    wth[3] = wth0.d;
    u0[3] = u00.d;
    v0[3] = v00.d;
    w0[3] = w00.d;
  }
  if (ns > 4) {
    uth[4] = uth0.e;
    vth[4] = vth0.e;
    wth[4] = wth0.e;
    u0[4] = u00.e;
    v0[4] = v00.e;
    w0[4] = w00.e;
  }
  if (ns > 5) {
    uth[5] = uth0.f;
    vth[5] = vth0.f;
    wth[5] = wth0.f;
    u0[5] = u00.f;
    v0[5] = v00.f;
    w0[1] = w00.f;
  }

  if (nstestpart > 0) {
		array_double pitch_angle0 = config.read < array_double > ("pitch_angle");
		array_double energy0 	  = config.read < array_double > ("energy");
		pitch_angle = new double[nstestpart];
		energy      = new double[nstestpart];
		if (nstestpart > 0) {
			pitch_angle[0] = pitch_angle0.a;
			energy[0] 	   = energy0.a;
		}
		if (nstestpart > 1) {
			pitch_angle[1] = pitch_angle0.b;
			energy[1] 	   = energy0.b;
		}
		if (nstestpart > 2) {
			pitch_angle[2] = pitch_angle0.c;
			energy[2] 	   = energy0.c;
		}
		if (nstestpart > 3) {
			pitch_angle[3] = pitch_angle0.d;
			energy[3] 	   = energy0.d;
		}
		if (nstestpart > 4) {
			pitch_angle[4] = pitch_angle0.e;
			energy[4] 	   = energy0.e;
		}
		if (nstestpart > 5) {
			pitch_angle[5] = pitch_angle0.f;
			energy[5] 	   = energy0.f;
		}
		if (nstestpart > 6) {
			pitch_angle[6] = pitch_angle0.g;
			energy[6] 	   = energy0.g;
		}
		if (nstestpart > 7) {
			pitch_angle[7] = pitch_angle0.h;
			energy[7] 	   = energy0.h;
		}
  }


  npcelx = new int[ns+nstestpart];
  npcely = new int[ns+nstestpart];
  npcelz = new int[ns+nstestpart];
  qom = new double[ns+nstestpart];
  array_int npcelx0 = config.read < array_int > ("npcelx");
  array_int npcely0 = config.read < array_int > ("npcely");
  array_int npcelz0 = config.read < array_int > ("npcelz");
  array_double qom0 = config.read < array_double > ("qom");
  npcelx[0] = npcelx0.a;
  npcely[0] = npcely0.a;
  npcelz[0] = npcelz0.a;
  qom[0]	  = qom0.a;
  int ns_tot =ns+nstestpart;
  if (ns_tot > 1) {
    npcelx[1] = npcelx0.b;
    npcely[1] = npcely0.b;
    npcelz[1] = npcelz0.b;
    qom[1]	= qom0.b;
  }
  if (ns_tot > 2) {
    npcelx[2] = npcelx0.c;
    npcely[2] = npcely0.c;
    npcelz[2] = npcelz0.c;
    qom[2] 	= qom0.c;
  }
  if (ns_tot > 3) {
    npcelx[3] = npcelx0.d;
    npcely[3] = npcely0.d;
    npcelz[3] = npcelz0.d;
    qom[3] 	= qom0.d;
  }
  if (ns_tot > 4) {
    npcelx[4] = npcelx0.e;
    npcely[4] = npcely0.e;
    npcelz[4] = npcelz0.e;
    qom[4] 	= qom0.e;
  }
  if (ns_tot > 5) {
    npcelx[5] = npcelx0.f;
    npcely[5] = npcely0.f;
    npcelz[5] = npcelz0.f;
    qom[5] 	= qom0.f;
  }
  if (ns_tot > 6) {
    npcelx[6] = npcelx0.g;
    npcely[6] = npcely0.g;
    npcelz[6] = npcelz0.g;
    qom[6] 	= qom0.g;
  }
  if (ns_tot > 7) {
    npcelx[7] = npcelx0.h;
    npcely[7] = npcely0.h;
    npcelz[7] = npcelz0.h;
    qom[7] 	= qom0.h;
  }
  if (ns_tot > 8) {
    npcelx[8] = npcelx0.i;
    npcely[8] = npcely0.i;
    npcelz[8] = npcelz0.i;
    qom[8] 	= qom0.i;
  }
  if (ns_tot > 9) {
    npcelx[9] = npcelx0.j;
    npcely[9] = npcely0.j;
    npcelz[9] = npcelz0.j;
    qom[9] 	= qom0.j;
  }
  if (ns_tot > 10) {
    npcelx[10] = npcelx0.k;
    npcely[10] = npcely0.k;
    npcelz[10] = npcelz0.k;
    qom[10] 	 = qom0.k;
  }
  if (ns_tot > 11) {
    npcelx[11] = npcelx0.l;
    npcely[11] = npcely0.l;
    npcelz[11] = npcelz0.l;
    qom[11] 	 = qom0.l;
  }



  //verbose = config.read < bool > ("verbose",false);

  // PHI Electrostatic Potential
  bcPHIfaceXright = config.read < int >("bcPHIfaceXright",1);
  bcPHIfaceXleft  = config.read < int >("bcPHIfaceXleft",1);
  bcPHIfaceYright = config.read < int >("bcPHIfaceYright",1);
  bcPHIfaceYleft  = config.read < int >("bcPHIfaceYleft",1);
  bcPHIfaceZright = config.read < int >("bcPHIfaceZright",1);
  bcPHIfaceZleft  = config.read < int >("bcPHIfaceZleft",1);

  // EM field boundary condition
  bcEMfaceXright = config.read < int >("bcEMfaceXright");
  bcEMfaceXleft  = config.read < int >("bcEMfaceXleft");
  bcEMfaceYright = config.read < int >("bcEMfaceYright");
  bcEMfaceYleft  = config.read < int >("bcEMfaceYleft");
  bcEMfaceZright = config.read < int >("bcEMfaceZright");
  bcEMfaceZleft  = config.read < int >("bcEMfaceZleft");

  /*  ---------------------------------------------------------- */
  /*  Electric and Magnetic field boundary conditions for BCface */
  /*  ---------------------------------------------------------- */
  // if bcEM* is 0: perfect conductor, if bcEM* is not 0: perfect mirror
  // perfect conductor: normal = free, perpendicular = 0
  // perfect mirror   : normal = 0,    perpendicular = free
  /*  ---------------------------------------------------------- */

  /* X component in faces Xright, Xleft, Yright, Yleft, Zright and Zleft (0, 1, 2, 3, 4, 5) */
  bcEx[0] = bcEMfaceXright == 0 ? 2 : 1;   bcBx[0] = bcEMfaceXright == 0 ? 1 : 2;
  bcEx[1] = bcEMfaceXleft  == 0 ? 2 : 1;   bcBx[1] = bcEMfaceXleft  == 0 ? 1 : 2;
  bcEx[2] = bcEMfaceYright == 0 ? 1 : 2;   bcBx[2] = bcEMfaceYright == 0 ? 2 : 1;
  bcEx[3] = bcEMfaceYleft  == 0 ? 1 : 2;   bcBx[3] = bcEMfaceYleft  == 0 ? 2 : 1;
  bcEx[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBx[4] = bcEMfaceZright == 0 ? 2 : 1;
  bcEx[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBx[5] = bcEMfaceZleft  == 0 ? 2 : 1;
  /* Y component */
  bcEy[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBy[0] = bcEMfaceXright == 0 ? 2 : 1;
  bcEy[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBy[1] = bcEMfaceXleft  == 0 ? 2 : 1;
  bcEy[2] = bcEMfaceYright == 0 ? 2 : 1;   bcBy[2] = bcEMfaceYright == 0 ? 1 : 2;
  bcEy[3] = bcEMfaceYleft  == 0 ? 2 : 1;   bcBy[3] = bcEMfaceYleft  == 0 ? 1 : 2;
  bcEy[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBy[4] = bcEMfaceZright == 0 ? 2 : 1;
  bcEy[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBy[5] = bcEMfaceZleft  == 0 ? 2 : 1;
  /* Z component */
  bcEz[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBz[0] = bcEMfaceXright == 0 ? 2 : 1;
  bcEz[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBz[1] = bcEMfaceXleft  == 0 ? 2 : 1;
  bcEz[2] = bcEMfaceYright == 0 ? 1 : 1;   bcBz[2] = bcEMfaceYright == 0 ? 2 : 1;
  bcEz[3] = bcEMfaceYleft  == 0 ? 1 : 1;   bcBz[3] = bcEMfaceYleft  == 0 ? 2 : 1;
  bcEz[4] = bcEMfaceZright == 0 ? 2 : 1;   bcBz[4] = bcEMfaceZright == 0 ? 1 : 2;
  bcEz[5] = bcEMfaceZleft  == 0 ? 2 : 1;   bcBz[5] = bcEMfaceZleft  == 0 ? 1 : 2;

  // Particles Boundary condition
  bcPfaceXright = config.read < int >("bcPfaceXright",1);
  bcPfaceXleft  = config.read < int >("bcPfaceXleft",1);
  bcPfaceYright = config.read < int >("bcPfaceYright",1);
  bcPfaceYleft  = config.read < int >("bcPfaceYleft",1);
  bcPfaceZright = config.read < int >("bcPfaceZright",1);
  bcPfaceZleft  = config.read < int >("bcPfaceZleft",1);


  if (RESTART1) {               // you are restarting
    RestartDirName = config.read < string > ("RestartDirName","data");
    ReadRestart(RestartDirName);
  }

  TrackParticleID = new bool[ns];
  array_bool TrackParticleID0 = config.read < array_bool > ("TrackParticleID");
  TrackParticleID[0] = TrackParticleID0.a;
  if (ns > 1)
    TrackParticleID[1] = TrackParticleID0.b;
  if (ns > 2)
    TrackParticleID[2] = TrackParticleID0.c;
  if (ns > 3)
    TrackParticleID[3] = TrackParticleID0.d;
  if (ns > 4)
    TrackParticleID[4] = TrackParticleID0.e;
  if (ns > 5)
    TrackParticleID[5] = TrackParticleID0.f;
}

bool Collective::field_output_is_off()const
{
  return (FieldOutputCycle <= 0);
}

bool Collective::particle_output_is_off()const
{
  return getParticlesOutputCycle() <= 0;
}
bool Collective::testparticle_output_is_off()const
{
  return getTestParticlesOutputCycle() <= 0;
}

/*! Read the collective information from the RESTART file in HDF5 format
 * There are three restart status: restart_status = 0 ---> new inputfile
 * restart_status = 1 ---> RESTART and restart and result directories does not coincide
 * restart_status = 2 ---> RESTART and restart and result directories coincide */
int Collective::ReadRestart(string inputfile) {
#ifdef NO_HDF5
  eprintf("restart requires compiling with HDF5");
#else
  restart_status = 1;
  // hdf stuff 
  hid_t file_id;
  hid_t dataset_id;
  herr_t status;
  // Open the setting file for the restart.
  file_id = H5Fopen((inputfile + "/settings.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << inputfile << endl;
    return -1;
  }

  // read c
  dataset_id = H5Dopen2(file_id, "/collective/c", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &c);
  status = H5Dclose(dataset_id);

  // read Lx 
  dataset_id = H5Dopen2(file_id, "/collective/Lx", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
  status = H5Dclose(dataset_id);
  // read Ly 
  dataset_id = H5Dopen2(file_id, "/collective/Ly", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
  status = H5Dclose(dataset_id);
  // read Lz 
  dataset_id = H5Dopen2(file_id, "/collective/Lz", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lz);
  status = H5Dclose(dataset_id);
  // read x_center
  dataset_id = H5Dopen2(file_id, "/collective/x_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x_center);
  status = H5Dclose(dataset_id);
  // read y_center
  dataset_id = H5Dopen2(file_id, "/collective/y_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y_center);
  status = H5Dclose(dataset_id);
  // read z_center
  dataset_id = H5Dopen2(file_id, "/collective/z_center", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z_center);
  status = H5Dclose(dataset_id);
  // read L_square
  dataset_id = H5Dopen2(file_id, "/collective/L_square", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L_square);
  status = H5Dclose(dataset_id);
  // read nxc
  dataset_id = H5Dopen2(file_id, "/collective/Nxc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nyc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyc);
  status = H5Dclose(dataset_id);
  // read nyc 
  dataset_id = H5Dopen2(file_id, "/collective/Nzc", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzc);
  status = H5Dclose(dataset_id);
  // read ns
  dataset_id = H5Dopen2(file_id, "/collective/Ns", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ns);
  status = H5Dclose(dataset_id);
  //read number of test particles species
  dataset_id = H5Dopen2(file_id, "/collective/NsTestPart", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nstestpart);
  status = H5Dclose(dataset_id);

  /*! Boundary condition information */
  // read EMfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXleft);
  status = H5Dclose(dataset_id);
  // read EMfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceXright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceXright);
  status = H5Dclose(dataset_id);
  // read EMfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYleft);
  status = H5Dclose(dataset_id);
  // read EMfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceYright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceYright);
  status = H5Dclose(dataset_id);
  // read EMfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZleft", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZleft);
  status = H5Dclose(dataset_id);
  // read EMfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/EMfaceZright", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcEMfaceZright);
  status = H5Dclose(dataset_id);

  // read PHIfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceXright);
  status = H5Dclose(dataset_id);
  // read PHIfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceYright);
  status = H5Dclose(dataset_id);
  // read PHIfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PHIfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPHIfaceZright);
  status = H5Dclose(dataset_id);

  // read PfaceXleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXleft);
  status = H5Dclose(dataset_id);
  // read PfaceXright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceXright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceXright);
  status = H5Dclose(dataset_id);
  // read PfaceYleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYleft);
  status = H5Dclose(dataset_id);
  // read PfaceYright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceYright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceYright);
  status = H5Dclose(dataset_id);
  // read PfaceZleft
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZleft", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZleft);
  status = H5Dclose(dataset_id);
  // read PfaceZright
  dataset_id = H5Dopen2(file_id, "/collective/bc/PfaceZright", H5P_DEFAULT); // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcPfaceZright);
  status = H5Dclose(dataset_id);
  // allocate fields depending on species
  npcelx = new int[ns+nstestpart];
  npcely = new int[ns+nstestpart];
  npcelz = new int[ns+nstestpart];
  qom = new double[ns+nstestpart];
  uth = new double[ns];
  vth = new double[ns];
  wth = new double[ns];
  u0 = new double[ns];
  v0 = new double[ns];
  w0 = new double[ns];
  // read data from species0, species 1, species2,...
  string *name_species = new string[ns];
  stringstream *ss = new stringstream[ns];
  string *name_testspecies;
  stringstream *testss;

  for (int i = 0; i < ns; i++) {
    ss[i] << i;
    name_species[i] = "/collective/species_" + ss[i].str() + "/";
  }
  if(nstestpart>0){
	  name_testspecies = new string[nstestpart];
	  testss = new stringstream[nstestpart];
	  for (int i = 0; i < nstestpart; i++) {
		  testss[i] << (i+ns);
		  name_testspecies[i] = "/collective/testspecies_" + testss[i].str() + "/";
	  }

	  pitch_angle = new double[nstestpart];
	  energy      = new double[nstestpart];
	  for (int i = 0; i < nstestpart; i++) {
	    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "pitch_angle").c_str(), H5P_DEFAULT); // HDF 1.8.8
	    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pitch_angle[i]);
	    status = H5Dclose(dataset_id);

	    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "energy").c_str(), H5P_DEFAULT); // HDF 1.8.8
	    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &energy[i]);
	    status = H5Dclose(dataset_id);
	  }
  }

  // npcelx for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelx").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  // npcely for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcely").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcely[i]);
    status = H5Dclose(dataset_id);
  }
  // npcelz for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "Npcelz").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelz[i]);
    status = H5Dclose(dataset_id);
  }
  // qom for different species
  for (int i = 0; i < ns; i++) {
    dataset_id = H5Dopen2(file_id, (name_species[i] + "qom").c_str(), H5P_DEFAULT);  // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &qom[i]);
    status = H5Dclose(dataset_id);
  }

  //Test Particle
  // npcelx for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcelx").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  // npcely for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcely").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcely[i]);
    status = H5Dclose(dataset_id);
  }
  // npcelz for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "Npcelz").c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npcelz[i]);
    status = H5Dclose(dataset_id);
  }
  // qom for different species
  for (int i = ns; i < (ns+nstestpart); i++) {
    dataset_id = H5Dopen2(file_id, (name_testspecies[i] + "qom").c_str(), H5P_DEFAULT);  // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &qom[i]);
    status = H5Dclose(dataset_id);
  }


  /*! not needed for restart * */
  for (int i = 0; i < ns; i++)
    uth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    vth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    wth[i] = 0.0;
  for (int i = 0; i < ns; i++)
    u0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    v0[i] = 0.0;
  for (int i = 0; i < ns; i++)
    w0[i] = 0.0;
  // verbose on
  //verbose = 1;


  // if RestartDirName == SaveDirName overwrite dt,Th,Smooth (append to old hdf files)
  if (RestartDirName == SaveDirName) {
    restart_status = 2;
    // read dt
    dataset_id = H5Dopen2(file_id, "/collective/Dt", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dt);
    status = H5Dclose(dataset_id);
    // read th 
    dataset_id = H5Dopen2(file_id, "/collective/Th", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &th);
    status = H5Dclose(dataset_id);
    // read Smooth
    dataset_id = H5Dopen2(file_id, "/collective/Smooth", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Smooth);
    status = H5Dclose(dataset_id);
    dataset_id = H5Dopen2(file_id, "/collective/SmoothNiter", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &SmoothNiter);
    status = H5Dclose(dataset_id);
  }

  status = H5Fclose(file_id);


  // read last cycle (not from settings, but from restart0.hdf)

  file_id = H5Fopen((inputfile + "/restart0.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: " << inputfile << endl;
    return -1;
  }

  dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT);  // HDF 1.8.8
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &last_cycle);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  // deallocate
  delete[]name_species;
  delete[]ss;
#endif
  return (0);
}



void Collective::read_field_restart(
    const VCtopology3D* vct,
    const Grid* grid,
    arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
    arr3_double Ex, arr3_double Ey, arr3_double Ez,
    array4_double* rhons_, int ns)const
{
#ifdef NO_HDF5
  eprintf("Require HDF5 to read from restart file.");
#else
    const int nxn = grid->getNXN();
    const int nyn = grid->getNYN();
    const int nzn = grid->getNZN();
    if (vct->getCartesian_rank() == 0)
    {
      printf("LOADING EM FIELD FROM RESTART FILE in %s/restart.hdf\n",getRestartDirName().c_str());
    }

    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = getRestartDirName() + "/restart" + ss.str() + ".hdf";

    // hdf stuff
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("Failed to open file: %s\n ", name_file.c_str());
    }

    //find the last cycle
    int lastcycle=0;
    dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &lastcycle);
    status = H5Dclose(dataset_id);

    // Bxn
    ss.str("");ss << "/fields/Bx/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    ss.str("");ss << "/fields/By/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    ss.str("");ss << "/fields/Bz/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    ss.str("");ss << "/fields/Ex/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey
    ss.str("");ss << "/fields/Ey/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez
    ss.str("");ss << "/fields/Ez/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species
    for (int is = 0; is < ns; is++) {
      ss.str("");ss << "/moments/species_" << is << "/rho/cycle_" << lastcycle;
      dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      status = H5Dclose(dataset_id);
      array4_double& rhons = *rhons_;
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int jj = 1; jj < nzn - 1; jj++)
            rhons[is][i][j][jj] = temp_storage[k++];
    }

    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
#endif
}

// extracted from Particles3Dcomm.cpp
//
void Collective::read_particles_restart(
    const VCtopology3D* vct,
    int species_number,
    vector_double& u,
    vector_double& v,
    vector_double& w,
    vector_double& q,
    vector_double& x,
    vector_double& y,
    vector_double& z,
    vector_double& t)const
{
#ifdef NO_HDF5
  eprintf("Require HDF5 to read from restart file.");
#else
    if (vct->getCartesian_rank() == 0 && species_number == 0)
    {
      printf("LOADING PARTICLES FROM RESTART FILE in %s/restart.hdf\n",
        getRestartDirName().c_str());
    }
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[1];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("couldn't open file: %s\n"
        "\tRESTART NOT POSSIBLE", name_file.c_str());
      //cout << "couldn't open file: " << name_file << endl;
      //cout << "RESTART NOT POSSIBLE" << endl;
    }


    //find the last cycle
    int lastcycle=0;
    dataset_id = H5Dopen2(file_id, "/last_cycle", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &lastcycle);
    status = H5Dclose(dataset_id);

    stringstream species_name;
    species_name << species_number;

    ss.str("");ss << "/particles/species_" << species_number << "/x/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id); /* dataspace handle */
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    // get how many particles there are on this processor for this species
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    const int nop = dims_out[0]; // number of particles in this process
    //Particles3Dcomm::resize_SoA(nop);
    {
      //
      // allocate space for particles including padding
      //
      const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
      u.reserve(padded_nop);
      v.reserve(padded_nop);
      w.reserve(padded_nop);
      q.reserve(padded_nop);
      x.reserve(padded_nop);
      y.reserve(padded_nop);
      z.reserve(padded_nop);
      t.reserve(padded_nop);
      //
      // define size of particle data
      //
      u.resize(nop);
      v.resize(nop);
      w.resize(nop);
      q.resize(nop);
      x.resize(nop);
      y.resize(nop);
      z.resize(nop);
      t.resize(nop);
    }
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x[0]);
    // close the data set
    status = H5Dclose(dataset_id);

    // get y
    ss.str("");ss << "/particles/species_" << species_number << "/y/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y[0]);
    status = H5Dclose(dataset_id);

    // get z
    ss.str("");ss << "/particles/species_" << species_number << "/z/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z[0]);
    status = H5Dclose(dataset_id);

    // get u
    ss.str("");ss << "/particles/species_" << species_number << "/u/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    status = H5Dclose(dataset_id);

    // get v
    ss.str("");ss << "/particles/species_" << species_number << "/v/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    status = H5Dclose(dataset_id);

    // get w
    ss.str("");ss << "/particles/species_" << species_number << "/w/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &w[0]);
    status = H5Dclose(dataset_id);

    // get q
    ss.str("");ss << "/particles/species_" << species_number << "/q/cycle_" << lastcycle;
    dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &q[0]);
    status = H5Dclose(dataset_id);

    // get ID
    if (false) {//TrackParticleID
		ss.str("");ss << "/particles/species_" << species_number << "/ID/cycle_" << lastcycle;
		dataset_id = H5Dopen2(file_id, ss.str().c_str(), H5P_DEFAULT); // HDF 1.8.8
		status = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t[0]);
		status = H5Dclose(dataset_id);
    }

    status = H5Fclose(file_id);
#endif
}



/*! constructor */
Collective::Collective(int argc, char **argv) {
  if (argc < 2) {
    inputfile = "inputfile";
    RESTART1 = false;
  }
  else if (argc < 3) {
    inputfile = argv[1];
    RESTART1 = false;
  }
  else {
    if (strcmp(argv[1], "restart") == 0) {
      inputfile = argv[2];
      RESTART1 = true;
    }
    else if (strcmp(argv[2], "restart") == 0) {
      inputfile = argv[1];
      RESTART1 = true;
    }
    else {
      cout << "Error: syntax error in mpirun arguments. Did you mean to 'restart' ?" << endl;
      return;
    }
  }
  ReadInput(inputfile);
  init_derived_parameters();
}

void Collective::init_derived_parameters()
{
  /*! fourpi = 4 greek pi */
  fourpi = 16.0 * atan(1.0);
  /*! dx = space step - X direction */
  dx = Lx / (double) nxc;
  /*! dy = space step - Y direction */
  dy = Ly / (double) nyc;
  /*! dz = space step - Z direction */
  dz = Lz / (double) nzc;
  /*! npcel = number of particles per cell */
  npcel = new int[ns+nstestpart];
  /*! np = number of particles of different species */
  //np = new int[ns];
  /*! npMax = maximum number of particles of different species */
  //npMax = new int[ns];

  /* quantities per process */

  // check that procs divides grid
  // (this restriction should be removed).
  //
  if(0==MPIdata::get_rank())
  {
    fflush(stdout);
    bool xerror = false;
    bool yerror = false;
    bool zerror = false;
    if(nxc % XLEN) xerror=true;
    if(nyc % YLEN) yerror=true;
    if(nzc % ZLEN) zerror=true;
    if(xerror) warning_printf("XLEN=%d does not divide nxc=%d\n", XLEN,nxc);
    if(yerror) warning_printf("YLEN=%d does not divide nyc=%d\n", YLEN,nyc);
    if(zerror) warning_printf("ZLEN=%d does not divide nzc=%d\n", ZLEN,nzc);
    fflush(stdout);
    bool error = (xerror||yerror||zerror);
    // Comment out this check if your postprocessing code does not
    // require the field output subarrays to be the same size.
    // Alternatively, you could modify the output routine to pad
    // with zeros...
    //if(error)
    //{
    //  eprintf("For WriteMethod=default processor dimensions "
    //          "must divide mesh cell dimensions");
    //}
  }

  int num_cells_r = nxc*nyc*nzc;
  //num_procs = XLEN*YLEN*ZLEN;
  //ncells_rs = nxc_rs*nyc_rs*nzc_rs;

  for (int i = 0; i < (ns+nstestpart); i++)
  {
    npcel[i] = npcelx[i] * npcely[i] * npcelz[i];
    //np[i] = npcel[i] * num_cells;
    //nop_rs[i] = npcel[i] * ncells_rs;
    //maxnop_rs[i] = NpMaxNpRatio * nop_rs[i];
    // INT_MAX is about 2 billion, surely enough
    // to index the particles in a single MPI process:
    //assert_le(NpMaxNpRatio * npcel[i] * ncells_proper_per_proc , double(INT_MAX));
    //double npMaxi = (NpMaxNpRatio * np[i]);
    //npMax[i] = (int) npMaxi;
  }
}

/*! destructor */
Collective::~Collective() {
  //delete[]np;
  delete[]npcel;
  delete[]npcelx;
  delete[]npcely;
  delete[]npcelz;
  //delete[]npMax;
  delete[]qom;

  delete[]uth;
  delete[]vth;
  delete[]wth;

  delete[]u0;
  delete[]v0;
  delete[]w0;

  delete[]TrackParticleID;

  delete[]rhoINIT;
  delete[]rhoINJECT;

  delete[]pitch_angle;
  delete[]energy;
}
/*! Print Simulation Parameters */
void Collective::Print() {
  cout << endl;
  cout << "Simulation Parameters" << endl;
  cout << "---------------------" << endl;
  cout << "Number of species    = " << ns << endl;
  for (int i = 0; i < ns; i++)
    cout << "qom[" << i << "] = " << qom[i] << endl;
  cout << "x-Length                 = " << Lx << endl;
  cout << "y-Length                 = " << Ly << endl;
  cout << "z-Length                 = " << Lz << endl;
  cout << "Number of cells (x)      = " << nxc << endl;
  cout << "Number of cells (y)      = " << nyc << endl;
  cout << "Number of cells (z)      = " << nzc << endl;
  cout << "Time step                = " << dt << endl;
  cout << "Number of cycles         = " << ncycles << endl;
  cout << "Results saved in  : " << SaveDirName << endl;
  cout << "Case type         : " << Case << endl;
  cout << "Simulation name   : " << SimName << endl;
  cout << "Poisson correction: " << PoissonCorrection << endl;
  cout << "---------------------" << endl;
  cout << "Check Simulation Constraints" << endl;
  cout << "---------------------" << endl;
  cout << "Accuracy Constraint:  " << endl;
  for (int i = 0; i < ns; i++) {
    cout << "u_th < dx/dt species " << i << ".....";
    if (uth[i] < (dx / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;

    cout << "v_th < dy/dt species " << i << "......";
    if (vth[i] < (dy / dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;
  }
  cout << endl;
  cout << "Finite Grid Stability Constraint:  ";
  cout << endl;
  for (int is = 0; is < ns; is++) {
    if (uth[is] * dt / dx > .1)
      cout << "OK u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " > .1" << endl;
    else
      cout << "WARNING. u_th*dt/dx (species " << is << ") = " << uth[is] * dt / dx << " < .1" << endl;

    if (vth[is] * dt / dy > .1)
      cout << "OK v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " > .1" << endl;
    else
      cout << "WARNING. v_th*dt/dy (species " << is << ") = " << vth[is] * dt / dy << " < .1"  << endl;

  }


}
/*! Print Simulation Parameters */
void Collective::save() {
  string temp;
  temp = SaveDirName + "/SimulationData.txt";
  ofstream my_file(temp.c_str());
  my_file << "---------------------------" << endl;
  my_file << "-  Simulation Parameters  -" << endl;
  my_file << "---------------------------" << endl;

  my_file << "Number of species    = " << ns << endl;
  for (int i = 0; i < ns; i++)
    my_file << "qom[%d] = " << qom[i] << endl;
  my_file << "---------------------------" << endl;
  my_file << "x-Length                 = " << Lx << endl;
  my_file << "y-Length                 = " << Ly << endl;
  my_file << "z-Length                 = " << Lz << endl;
  my_file << "Number of cells (x)      = " << nxc << endl;
  my_file << "Number of cells (y)      = " << nyc << endl;
  my_file << "Number of cells (z)      = " << nzc << endl;
  my_file << "---------------------------" << endl;
  my_file << "Time step                = " << dt << endl;
  my_file << "Number of cycles         = " << ncycles << endl;
  my_file << "---------------------------" << endl;
  for (int is = 0; is < ns; is++){
    my_file << "rho init species   " << is << " = " << rhoINIT[is] << endl;
    my_file << "rho inject species " << is << " = " << rhoINJECT[is]  << endl;
  }
  my_file << "current sheet thickness  = " << delta << endl;
  my_file << "B0x                      = " << B0x << endl;
  my_file << "BOy                      = " << B0y << endl;
  my_file << "B0z                      = " << B0z << endl;
  my_file << "---------------------------" << endl;
  my_file << "Smooth                   = " << Smooth << endl;
  my_file << "SmoothNiter              = " << SmoothNiter<< endl;
  my_file << "GMRES error tolerance    = " << GMREStol << endl;
  my_file << "CG error tolerance       = " << CGtol << endl;
  my_file << "Mover error tolerance    = " << NiterMover << endl;
  my_file << "---------------------------" << endl;
  my_file << "Results saved in: " << SaveDirName << endl;
  my_file << "Restart saved in: " << RestartDirName << endl;
  my_file << "---------------------" << endl;
  my_file.close();

}

