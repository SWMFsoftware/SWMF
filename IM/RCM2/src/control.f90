program control_RCM
  use ModControl
  use rcm_variables
  use rcm_io
  implicit none
  include "mpif.h"

  integer :: ierror
  integer :: limit(3)
  integer :: time_i, time_f, idt2, rec_i, rec_f, istep, time_i_overall, time_f_overall
  REAL, dimension (isize,jsize) :: MHD_rho, MHD_p, MHD_vol, MHD_xmin, MHD_ymin,MHD_v
  REAL, dimension (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc) :: vm_tmp, xmin_tmp, ymin_tmp

  ! setup mpi
  call MPI_INIT(ierror)
  comm_WORLD = MPI_COMM_WORLD
  call MPI_COMM_RANK(comm_WORLD, pe_WORLD, ira)
  call MPI_COMM_SIZE(comm_WORLD, numPEs_WORLD, ira)

  ! Put everyone into one big group
  call MPI_COMM_GROUP(comm_WORLD, group_WORLD, ira)

  !\
  ! Create MPI subgroups
  !/

  ! RCM

  ! Put processors into RCM group
  limit(1) = 0                   !first
! limit(2) = min(0,numPEs_WORLD-1)   !last
  limit(2) = numPEs_WORLD-1
  limit(3) = 1                   !stride
  call MPI_Group_Range_Incl(group_WORLD, 1, limit, group_RCM, ira)
  pe0_RCM = limit(1)

  ! Get the rank of a PE withing this group
  call MPI_Group_Rank(group_RCM, pe_RCM, ira)

  ! Create the communicator for this group
  call MPI_Comm_Create(comm_WORLD, group_RCM, comm_RCM, ira)

  if (pe_RCM /= MPI_UNDEFINED) then
     ! Get the number of PEs in this group
     call MPI_COMM_SIZE(comm_RCM, numPEs_RCM, ira)
  end if
  call MPI_Bcast(numPEs_RCM,1,MPI_Integer,pe0_RCM,comm_WORLD,ira)

!
! Get magnetic field for RCM (constant):
  CALL Read_array ('rcmvm_inp', 1, label, ARRAY_2D=vm_tmp, ASCI=asci_flag)
  CALL Read_array ('rcmxmin_inp', 1, label, ARRAY_2D=xmin_tmp, ASCI=asci_flag)
  CALL Read_array ('rcmymin_inp', 1, label, ARRAY_2D=ymin_tmp, ASCI=asci_flag)
  MHD_vol = vm_tmp(1:isize,1:jsize)
  MHD_xmin = xmin_tmp(1:isize,1:jsize)
  MHD_ymin = ymin_tmp(1:isize,1:jsize)
  CALL Rcm_putfieldvolume (MHD_vol, 1, MHD_xmin, MHD_ymin)
  fstoff = label%real(13)
!
! Get plasma boundary conditions for RCM (constant):
! MHD_rho = 0.4E+6*xmass(2) ! in kg/m3
! MHD_p   = 5000.0*1.6E-19*0.4E+6 ! in Pa
  CALL Read_array ('rcmdensity_inp', 1, label, ARRAY_2D=density, ASCI=asci_flag)
  CALL Read_array ('rcmtemperature_inp', 1, label, ARRAY_2D=temperature, ASCI=asci_flag)
  MHD_rho = density(1:isize,1:jsize)*xmass(2)*1.0E+6 ! in kg/m3
  MHD_p   = temperature(1:isize,1:jsize)*1.6E-19*1.0E+6*density(1:isize,1:jsize) ! in Pa
  CALL Rcm_putrhop (MHD_rho, MHD_p)
!
! Get MHD potential:
  CALL Read_array ('rcmv_inp', 1, label, ARRAY_2D=v, ASCI=asci_flag)
  MHD_v = v (1:isize,1:jsize)
  CALL Rcm_putpotential (MHD_v)
!
! CALL RCM:
!
! Initialize (read in files, etc):
  call RCM_advec (1, 0, 14400, 2, 1, 1, 600)
!
!
! Run calculational loop:
  time_i_overall = 0
  time_f_overall = 14400
  idt2 = 600
  DO istep = 1, (time_f_overall-time_i_overall)/idt2
     time_i = (istep-1)*idt2 + time_i_overall
     time_f = time_i + idt2
     rec_i = 1 + (time_f-time_i)/idt2*(istep-1)
     rec_f = rec_i
     print*,'calling rcm:', time_i, time_f, rec_i, rec_f
     CALL RCM_advec (2, time_i, time_f, 2, rec_i, rec_f, idt2)
  END DO
!
! Finalize:
  call RCM_advec (3, 0, 14400, 2, 1, 1, 600)
!
  CALL MPI_FINALIZE (ira)
  STOP
end program control_RCM
