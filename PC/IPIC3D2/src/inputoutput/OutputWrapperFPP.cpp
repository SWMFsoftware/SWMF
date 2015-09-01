
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
    output_file = SaveDirName + "/proc"   + num_proc_str + ".hdf";
    restart_file= SaveDirName + "/restart"+ num_proc_str + ".hdf";

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
          hdf5_agent.open(SaveDirName + "/settings.hdf");
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
  		  hdf5_agent.open(RestartDirName + "/settings.hdf");
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
		hdf5_agent.open_append(restart_file);
		output_mgr.output("proc_topology ", cycle);
		output_mgr.output("Eall + Ball + rhos + Js + pressure", cycle);
		output_mgr.output("position + velocity + q ", cycle, 0);
		output_mgr.output("testpartpos + testpartvel + testpartcharge", cycle, 0);
		output_mgr.output("last_cycle", cycle);
		hdf5_agent.close();
#endif
}
