!$Id$
!the wrapper for reading BATSRUS data files into AMPS

! The substitute for the subroutine called by SWMF shared subroutines
! in the case of an error
!subroutine CON_stop(StringError)
!  implicit none
!  character (len=*), intent(in) :: StringError
!  !----------------------------------------------------------------------------
!  write(*,'(a)')StringError
!  write(*,'(a)')'!!! AMPS_ABORT !!!'
!  stop
!end subroutine CON_stop

!==============================================================================                                                   
! The substitute for the timing subroutines called by SWMF
subroutine timing_start(name)
  character (LEN=*), intent(in):: name
end subroutine timing_start

subroutine timing_stop(name)
  character (LEN=*), intent(in):: name
end subroutine timing_stop

!=============================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the number of the variables in the state vector that are returned by the interpolation routine
subroutine batsrus2amps_get_nvar(res)
  use ModReadAmr, ONLY: nVar

  implicit none
  integer,intent(out)::res
  !-----------------------------------------------------  
  res=nVar 
end subroutine batsrus2amps_get_nvar    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get the string that containes the names of the variables that are returned by the interpolation routine 
subroutine batsrus2amps_get_namevardata(varlist,varlistlength) 
  use ModReadAmr, ONLY:NameVarData
  
  implicit none
  integer,intent(in)::varlistlength
  character(len=varlistlength),intent(out)::varlist
  !-----------------------------------------------------  
  varlist(1:len(NameVarData))=NameVarData(1:len(NameVarData))
end subroutine batsrus2amps_get_namevardata  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get the string that containes defienition of the units of the 1. the point coordinate and 2. the interpolated variables
subroutine batsrus2amps_get_nameunitdata(unitlist,unitlistlength) 
  use ModReadAmr, ONLY:NameUnitData
  
  implicit none
  integer,                      intent(in )::unitlistlength
  character(len=unitlistlength),intent(out)::unitlist
  !-----------------------------------------------------
  unitlist(1:len(NameUnitData))=NameUnitData(1:len(NameUnitData))
end subroutine batsrus2amps_get_nameunitdata  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Determine the boundary of the computational domain in the normalized units
subroutine batsrus2amps_domain_limits(xmin,xmax) 
  use ModReadAmr, ONLY: CoordMin_D, CoordMax_D
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_),intent(out)::xmin(3),xmax(3)
  !-----------------------------------------------------
  xmin(:)=CoordMin_D(:)
  xmax(:)=CoordMax_D(:)
end subroutine batsrus2amps_domain_limits  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read the data file (.idl)
subroutine batsrus2amps_openfile(FileName,FileNameLength)
  use ModReadAmr, ONLY:readamr_read 
  
  implicit none
  integer,                      intent(in):: FileNameLength
  character(len=FileNameLength),intent(in):: FileName
  !-----------------------------------------------------
  call readamr_read(FileName, IsNewGridIn = .false., IsVerboseIn=.true.)
end subroutine batsrus2amps_openfile 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialize the reading procedure. Read the header file. After this step it is possible to find the list of the variables and units as well as the limits of the computational domain
subroutine batsrus2amps_read_file_header(FileName,FileNameLength)
  use ModReadAmr, ONLY:readamr_init
  
  implicit none
  integer,                      intent(in):: FileNameLength
  character(len=FileNameLength),intent(in):: FileName
  integer::l
  
  l = index(FileName,".",BACK=.true.)
  call readamr_init(FileName(1:l-1), .true.)
end subroutine batsrus2amps_read_file_header
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Close the data file and deallocate all memory buffers allocated by the interpolation routine
subroutine batsrus2amps_closefile()
  use ModReadAmr, ONLY:readamr_clean
  
  call readamr_clean
end subroutine batsrus2amps_closefile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get the interpolated values in a point
subroutine batsrus2amps_get_data_point(x,res,FoundFlag) 
  use ModReadAmr, ONLY:nVar,readamr_get
  use BATL_lib,   ONLY:nDim
  use ModKind,    ONLY: Real8_
  
  implicit none
  real(Real8_),dimension(3),     intent(in )::x
  real(Real8_),dimension(0:nVar),intent(out)::res
  integer,                       intent(out)::FoundFlag
  
  real(Real8_),dimension(0:nVar)::State
  real(Real8_),dimension(nDim)::Xyz_D
  logical::IsFound
  !-----------------------------------------------------  
  Xyz_D(:)=x(:)
  call readamr_get(Xyz_D, State, IsFound)
  
  res(:)=State(:)
  if(IsFound)then
     FoundFlag=1
  else
     FoundFlag=0
  end if
end subroutine batsrus2amps_get_data_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Export the MPI parameters into the interpolation routine
subroutine batsrus2amps_set_mpi_parameters(ThisThread,nTotalThreads,Communicator)
  use BATL_lib,  ONLY: iProc,nProc,iComm
  
  implicit none
  integer,intent(in)::ThisThread,nTotalThreads,Communicator
  !-----------------------------------------------------
  iProc=ThisThread
  nProc=nTotalThreads
  iComm=Communicator
end subroutine batsrus2amps_set_mpi_parameters















  
  
