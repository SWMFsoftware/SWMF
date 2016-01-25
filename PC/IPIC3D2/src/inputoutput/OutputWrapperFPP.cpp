
#include "mpi.h"
#include "OutputWrapperFPP.h"

#include "CollectiveIO.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"

void OutputWrapperFPP::init_output_files(
	    Collective    *col,
	    VCtopology3D  *vct,
	    Grid3DCU      *grid,
	    EMfields3D    *EMf,
	    Particles3D   *part,
	    int 		  ns,
	    Particles3D   *testpart,
	    int 		  nstestpart)
{
#ifndef NO_HDF5
    cartesian_rank = vct->getCartesian_rank();
    stringstream num_proc_ss;
    num_proc_ss << cartesian_rank;
    string num_proc_str = num_proc_ss.str();
    SaveDirName = col->getSaveDirName();
    RestartDirName = col->getRestartDirName();
    int restart_status = col->getRestart_status();

    string iRegion;
    if(col->getCase()=="BATSRUS"){
      stringstream ss;
      ss<<col->getiRegion();
      iRegion = "_region"+ss.str();
      RestartDirName = "PC/restartOUT";
    }else{
      iRegion="";
    }
    
    output_file = SaveDirName + "/proc"   + num_proc_str + iRegion + ".hdf";
    if(col->getCase()=="BATSRUS"){
      restart_file= RestartDirName + "/restart"+ num_proc_str + iRegion + ".hdf";
    }else{
      restart_file= SaveDirName + "/restart"+ num_proc_str + iRegion + ".hdf";
    }
    
    // Initialize the output (simulation results and restart file)
    hdf5_agent.set_simulation_pointers(EMf, grid, vct, col);

    for (int i = 0; i < ns; ++i){
      hdf5_agent.set_simulation_pointers_part(&part[i]);
    }
    for (int i = 0; i < nstestpart; ++i){
      hdf5_agent.set_simulation_pointers_part(&testpart[i]);
    }

    // Add the HDF5 output agent to the Output Manager's list
    output_mgr.push_back(&hdf5_agent);

    if(col->getWriteMethod() == "shdf5"||(col->getWriteMethod()=="pvtk"&&!col->particle_output_is_off()) ){
        if (cartesian_rank == 0 && restart_status < 2) {
	  stringstream filename;

	  if(col->getCase()=="BATSRUS")
	    filename<<"/settings"<<"_region"<<col->getiRegion()<<".hdf";
	  else
	    filename<<"/settings.hdf";
	  
          hdf5_agent.open(SaveDirName + filename.str());
          output_mgr.output("collective + total_topology + proc_topology", 0);
          hdf5_agent.close();
        }

    	if (restart_status == 0) {
    	      hdf5_agent.open(output_file);
    	}else {
    	      hdf5_agent.open_append(output_file);
    	    }
		output_mgr.output("proc_topology ", 0);
		hdf5_agent.close();
    }

    if(col->getCallFinalize() || col->getRestartOutputCycle()>0){

        if (cartesian_rank == 0 && restart_status < 2) {
	  stringstream filename;
	  if(col->getCase()=="BATSRUS"){
	    filename<<"/settings"<<"_region"<<col->getiRegion()<<".hdf";
	  }else{
	    filename<<"/settings.hdf"<<endl;
	  }

	  hdf5_agent.open(RestartDirName + filename.str());
	  output_mgr.output("collective + total_topology + proc_topology", 0);
	  hdf5_agent.close();
        }

    	if (restart_status == 0) {
    		hdf5_agent.open(restart_file);
    		hdf5_agent.close();
    	}

    }


#endif
}

void OutputWrapperFPP::append_output(const char* tag, int cycle)
{
#ifndef NO_HDF5
    hdf5_agent.open_append(output_file);
    output_mgr.output(tag, cycle);
    hdf5_agent.close();
#endif
}

void OutputWrapperFPP::append_output(const char* tag, int cycle, int sample)
{
#ifndef NO_HDF5
    hdf5_agent.open_append(output_file);
    output_mgr.output(tag, cycle, sample);
    hdf5_agent.close();
#endif
}

void OutputWrapperFPP::append_restart(int cycle)
{
#ifndef NO_HDF5
  int cycle0=cycle;
#ifdef BATSRUS
  cycle0 = 0;
#endif

  hdf5_agent.open_append(restart_file);
  output_mgr.output("proc_topology ", cycle0);  
  output_mgr.output("Eall + Ball + rhos + Js + pressure", cycle0);
  #ifdef BATSRUS
  output_mgr.output("Bcall", cycle0);
  output_mgr.output("pseudo_random_seed", cycle0,0);
  #endif
  output_mgr.output("position + velocity + q ", cycle0, 0);
  output_mgr.output("testpartpos + testpartvel + testpartcharge", cycle0, 0);
  output_mgr.output("last_cycle", cycle);
  hdf5_agent.close();
#endif
}

#ifdef BATSRUS
void OutputWrapperFPP::getFluidState(int nDim, int nPoint, double *Xyz_I,
				     double *data_I,int nVar , string *NameVar){
  hdf5_agent.getFluidSIatPoint(nDim, nPoint, Xyz_I, data_I, nVar, NameVar);
}

#endif
