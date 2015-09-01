
#include <mpi.h>
#include <stdarg.h>
#include "TimeTasks.h"
#include "asserts.h"
#include "MPIdata.h" // for get_rank
#include "parallel.h"
#include "debug.h"
#include "errors.h"
#include "Collective.h"
#include "string.h" // for strcmp

/** implementation of declarations in utility/TimeTasks.h **/

TimeTasks timeTasks;

static const char *taskNames[] = // order must agree with Tasks in TimeTasks.h
{
  "none = 0",
  //
  "fields",
  "particles",
  "moments",
  "AFTER_EXCLUSIVE",
  //
  "communicating",
  "BEFORE_COMMUNICATION",
  "flds_comm",
  "flds_mpi_allreduce",
  "flds_mpi_sendrecv",
  "pcls_comm",
  "pcls_mpi_allreduce",
  "pcls_mpi_sendrecv",
  "moms_comm",
  "moms_mpi_allreduce",
  "moms_mpi_sendrecv",
  //
  "BEFORE_REPORT_LIST",
  "reduce_fields",
  "bfield",
  "moment_pcl_sorting",
  "moment_accumulation",
  "moment_reduction",
  "mover_pcl_sorting",
  "mover_pcl_moving",
  "transpose_pcls_to_aos",
  "transpose_pcls_to_soa",
  "write_fields",
  "write_particles",
  "PCLS_MPI_Isend",
  "PCLS_MPI_Irecv",
  "PCLS_MPI_Wait",
  "PCLS_MPI_Cancel",
  "PCLS_MPI_Request_free",
  "PCLS_MPI_Test",
  "PCLS_MPI_Waitany",
  //
  "number_of_tasks"
};

const char* TimeTasks::get_taskname(int arg)
{
  assert_le(arg,NUMBER_OF_TASKS);
  return taskNames[arg];
}

void TimeTasks::resetCycle()
{
  assert(!strcmp(taskNames[NUMBER_OF_TASKS],"number_of_tasks"));

  for(int e=0;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]=0.;
    //communicate[e]=0.;
    //sendrecv[e]=0.;
    //allreduce[e]=0.;
    active[e]=false;
    stack_depth[e]=0;
    start_times[e]=0.;
  }
  active_task=NONE;
  //communicating=false;
}
void TimeTasks::start_main_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(is_exclusive(taskid));
  assert_ne(active_task, taskid);
  active_task = taskid;
  assert(!active[taskid]);
  active[taskid]=true;
}
void TimeTasks::start_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(!is_exclusive(taskid));
  assert(!active[taskid]);
  active[taskid]=true;

  // special behavior per task
  switch(taskid)
  {
    default:
      break;
    // trigger appropriate communication task
    case COMMUNICATING:
      switch(active_task)
      {
        default:
          invalid_value_error(active_task);
        case FIELDS:
          timeTasks.start_task(FLDS_COMM);
          break;
        case PARTICLES:
          timeTasks.start_task(PCLS_COMM);
          break;
        case MOMENTS:
          timeTasks.start_task(MOMS_COMM);
          break;
      }
      break;
    case REDUCE_FIELDS:
      assert_eq(active_task, FIELDS);
      break;
      }
}
// have to manage the task stack explicitly
void TimeTasks::start_task(TimeTasks::Tasks taskid, double start_time)
{
  if(omp_get_thread_num()) return;
  if(stack_depth[taskid]==0)
  {
    start_times[taskid]=start_time;
    start_task(taskid);
  }
  stack_depth[taskid]++;
}
void TimeTasks::end_main_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  end_task(taskid, start_time);
  active_task = NONE;
}
void TimeTasks::end_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert(active[taskid]);
  double now = MPI_Wtime();
  // compute time spent on task
  task_duration[taskid] += now - start_time;
  active[taskid] = false;

  switch(taskid)
  {
    default:
      break;
    // trigger appropriate communication task
    case COMMUNICATING:
      switch(active_task)
      {
        default:
          invalid_value_error(active_task);
        case FIELDS:
          timeTasks.end_task(FLDS_COMM, start_time);
          break;
        case PARTICLES:
          timeTasks.end_task(PCLS_COMM, start_time);
          break;
        case MOMENTS:
          timeTasks.end_task(MOMS_COMM, start_time);
          break;
      }
      break;
    case REDUCE_FIELDS:
      assert_eq(active_task, FIELDS);
      break;
      }
}
// have to manage the task stack explicitly
void TimeTasks::end_task(TimeTasks::Tasks taskid)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  stack_depth[taskid]--;
  assert_ge(stack_depth[taskid],0);
  if(stack_depth[taskid]==0)
  {
    end_task(taskid, start_times[taskid]);
  }
}
// update appropriate XXXX_MPI_SENDRECV task
//
void TimeTasks::end_sendrecv(double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  int sendrecv_task;
  switch(active_task)
  {
    default:
      unsupported_value_error(get_taskname(active_task));
    case FIELDS:
      sendrecv_task = FLDS_MPI_SENDRECV;
      break;
    case PARTICLES:
      sendrecv_task = PCLS_MPI_SENDRECV;
      break;
    case MOMENTS:
      sendrecv_task = MOMS_MPI_SENDRECV;
      break;
  }
  const double additional_time = MPI_Wtime()-start_time;

  task_duration[sendrecv_task] += additional_time;
}

void TimeTasks::end_allreduce(double start_time)
{
  assert(!omp_get_thread_num()); //if(omp_get_thread_num()) return;
  assert_eq(active_task,FIELDS);
  double additional_communication_time = MPI_Wtime()-start_time;
  task_duration[FLDS_MPI_ALLREDUCE] += additional_communication_time;
}

void TimeTasks::print_cycle_times(int cycle,
  double* tskdur,
  const char* reduce_mode)
{
  // restrict output to master process
  //
  if(MPIdata::get_rank()) return;
  FILE* file = stdout;
  {
    fflush(file);
    double commun[NUMBER_OF_TASKS];
    double sndrcv[NUMBER_OF_TASKS];
    double allred[NUMBER_OF_TASKS];
    //commun[NUMBER_OF_TASKS];
    commun[FIELDS] = tskdur[FLDS_COMM];
    sndrcv[FIELDS] = tskdur[FLDS_MPI_SENDRECV];
    allred[FIELDS] = tskdur[FLDS_MPI_ALLREDUCE];
    commun[PARTICLES] = tskdur[PCLS_COMM];
    sndrcv[PARTICLES] = tskdur[PCLS_MPI_SENDRECV];
    allred[PARTICLES] = tskdur[PCLS_MPI_ALLREDUCE];
    commun[MOMENTS] = tskdur[MOMS_COMM];
    sndrcv[MOMENTS] = tskdur[MOMS_MPI_SENDRECV];
    allred[MOMENTS] = tskdur[MOMS_MPI_ALLREDUCE];
    sndrcv[PARTICLES] +=
      tskdur[PCLS_MPI_Isend]+
      tskdur[PCLS_MPI_Irecv]+
      tskdur[PCLS_MPI_Wait]+
      tskdur[PCLS_MPI_Cancel]+
      tskdur[PCLS_MPI_Request_free]+
      tskdur[PCLS_MPI_Test]+
      tskdur[PCLS_MPI_Waitany];
    double tskdurtot=0.;
    double computtot=0.;
    double communtot=0.;
    double allredtot=0.;
    double sndrcvtot=0.;
    fprintf(file, "%s_|total  comput commun task\n", reduce_mode);
    assert_eq(FIELDS+2,MOMENTS);
    for(int e=FIELDS; e<=MOMENTS; e++)
    {
      const double comput = tskdur[e]-commun[e];
      tskdurtot += tskdur[e];
      computtot += comput;
      communtot += commun[e];
      allredtot += allred[e];
      sndrcvtot += sndrcv[e];
      fprintf(file, "%s_|%6.3f %6.3f %6.3f %s\n",
        reduce_mode,
        tskdur[e],
        comput,
        commun[e],
        //loccom[e],
        get_taskname(e));
    }

    // report total times
    fprintf(file, "%s_|%6.3f %6.3f %6.3f %s\n",
      reduce_mode,
      tskdurtot,
      computtot,
      communtot,
      //loccomtot,
      "[total times]");

   
    fflush(file);
  }
}

static void reduce_doubles(int len, MPI_Op mpi_op, void* inbuff, void* outbuff)
{
    MPI_Allreduce(inbuff,outbuff,
      len,MPI_DOUBLE,mpi_op,MPI_COMM_WORLD);
}

void TimeTasks::print_cycle_times(int cycle, const char* reduce_mode)
{
  MPI_Op mpi_op;
  if(!strcmp(reduce_mode,"max"))
    mpi_op = MPI_MAX;
  else if(!strcmp(reduce_mode,"min"))
    mpi_op = MPI_MIN;
  else if(!strcmp(reduce_mode,"avg"))
    mpi_op = MPI_SUM;
  else
    invalid_value_error(reduce_mode);

  // assume that only main thread is active
  assert(!omp_get_thread_num());
  // perform all-reduce to get max times for all processes
  double reported_task_duration[NUMBER_OF_TASKS];
  reduce_doubles(NUMBER_OF_TASKS, mpi_op, task_duration, reported_task_duration);
  if(!strcmp(reduce_mode,"avg"))
  {
    const int nprocs = MPIdata::get_nprocs();
    for(int i=0;i<NUMBER_OF_TASKS;i++)
    {
      reported_task_duration[i] /=nprocs;
    }
  }
  print_cycle_times(cycle, reported_task_duration, reduce_mode);
}

void TimeTasks::print_cycle_times(int cycle)
{
  assert(!omp_get_thread_num());

  if(!MPIdata::get_rank()) fflush(stdout);
  //printf0("=== times for cycle %d (main process) ===\n", cycle);
  //print_cycle_times(cycle, task_duration, "main");
  //printf0("=== times for cycle %d (maximum over all processes) ===\n", cycle);
  //print_cycle_times(cycle, task_duration, "max");
  printf0("=== times for cycle %d (averaged over all processes) ===\n", cycle);
  print_cycle_times(cycle, task_duration, "avg");
  //printf0("=== times for cycle %d (minimum over all processes) ===\n", cycle);
  //print_cycle_times(cycle, task_duration, "min");
  if(!MPIdata::get_rank()) fflush(stdout);
}

// The following three methods provide for a hack by which
// the timeTasks copies of all threads are averaged.
// 
void TimeTasks::operator/=(int num)
{
  assert(false); // this method is not in use.
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]/=num;
    start_times[e]/=num;
  }
}
void TimeTasks::operator+=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    assert_eq(active[e], arg.active[e]);
    assert_eq(stack_depth[e], arg.stack_depth[e]);
    task_duration[e]+=arg.task_duration[e];
    start_times[e]+=arg.start_times[e];
  }
}
void TimeTasks::operator=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    active[e] = arg.active[e];
    task_duration[e]=arg.task_duration[e];
    stack_depth[e] = arg.stack_depth[e];
    start_times[e] = arg.start_times[e];
  }
}

TimeTasks_caller_to_set_main_task_for_scope::
TimeTasks_caller_to_set_main_task_for_scope(TimeTasks::Tasks _task) :
  task(_task)
{
  //if(omp_get_thread_num()) return;
  // assume that only one thread is active
  assert(!omp_get_thread_num());
  start_time = MPI_Wtime();
  timeTasks.start_main_task(task);
}
TimeTasks_caller_to_set_main_task_for_scope::
~TimeTasks_caller_to_set_main_task_for_scope()
{
  //if(omp_get_thread_num()) return;
  // assume that only one thread is active
  assert(!omp_get_thread_num());

  //MPI_Barrier(MPI_COMM_WORLD);
  timeTasks.end_main_task(task, start_time);
}

TimeTasks_caller_to_set_task_for_scope::
TimeTasks_caller_to_set_task_for_scope(TimeTasks::Tasks task_)
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  task = task_;
  already_active = timeTasks.is_active(task);
  // if the task is already active then
  // we cannot tell timeTasks to start it.
  if(already_active)
    return;
  //#pragma omp barrier
  start_time = MPI_Wtime();
  timeTasks.start_task(task);

}
TimeTasks_caller_to_set_task_for_scope::
~TimeTasks_caller_to_set_task_for_scope()
{
  assert(!omp_get_thread_num()); // if(omp_get_thread_num()) return;
  if(already_active)
  {
    assert(timeTasks.is_active(task));
    return;
  }
  //#pragma omp barrier
  timeTasks.end_task(task, start_time);
}

