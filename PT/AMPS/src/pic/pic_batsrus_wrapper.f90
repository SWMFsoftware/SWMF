!$Id$
!the wrapper for reading BATSRUS data files into AMPS

!use BATL_lib, ONLY: iProc, nProc, iComm, nI, nJ, nK, nDim, init_mpi, clean_mpi, coord_to_xyz
!use ModReadAmr, ONLY: nVar, CoordMin_D, CoordMax_D, readamr_read, readamr_get, readamr_clean
!use ModConst, ONLY: cPi
!use ModMpi, ONLY: MPI_REAL, MPI_SUM, MPI_allreduce


subroutine batsrus2amps_get_nvar(res)
  use ModReadAmr, ONLY: nVar

  implicit none
  integer,intent(out)::res
  
  res=nVar 
end subroutine batsrus2amps_get_nvar    

subroutine batsrus2amps_get_namevardata(varlist,varlistlength) 
  use ModReadAmr, ONLY:NameVarData
  
  implicit none
  integer,intent(in)::varlistlength
  character(len=varlistlength),intent(out)::varlist
  
  varlist(1:varlistlength)=NameVarData(1:varlistlength)
end subroutine batsrus2amps_get_namevardata  
  

subroutine batsrus2amps_domain_limits(xmin,xmax) 
  use ModReadAmr, ONLY: CoordMin_D, CoordMax_D
  
  implicit none
  real(8),intent(out)::xmin(:),xmax(:)
  
  xmin(:)=CoordMin_D(:)
  xmax(:)=CoordMax_D(:)
end subroutine batsrus2amps_domain_limits  

subroutine batsrus2amps_openfile(FileName,FileNameLength)
  use ModReadAmr, ONLY:readamr_read 
  
  implicit none
  integer,intent(in)::FileNameLength
  character(len=FileNameLength),intent(in):: FileName

  call readamr_read(FileName, IsNewGridIn = .true., IsVerboseIn=.true., UseCoordTest=.false.)
end subroutine batsrus2amps_openfile 

subroutine batsrus2amps_read_file_header(FileName,FileNameLength)
  use ModReadAmr, ONLY:readamr_init
  
  implicit none
  integer,intent(in)::FileNameLength
  character(len=FileNameLength),intent(in):: FileName
  integer::l
  
  l = index(FileName,".",BACK=.true.)
  call readamr_init(FileName(1:l-1), .true.)
end subroutine batsrus2amps_read_file_header
  

subroutine batsrus2amps_closefile()
  use ModReadAmr, ONLY:readamr_clean
  
  call readamr_clean
end subroutine batsrus2amps_closefile

subroutine batsrus2amps_get_data_point(x,res,FoundFlag) 
  use ModReadAmr, ONLY:nVar,readamr_get
  use BATL_lib, ONLY:nDim
  
  implicit none
  real(8),dimension(3),intent(in)::x
  real(8),dimension(0:nVar),intent(out)::res
  integer,intent(out)::FoundFlag
  
  real,dimension(0:nVar)::State
  real,dimension(nDim)::Xyz_D
  logical::IsFound
  
  Xyz_D(:)=x(:)
  call readamr_get(Xyz_D, State, IsFound)
  
  res(:)=State(:)
  FoundFlag=IsFound
  
  
end subroutine batsrus2amps_get_data_point


subroutine batsrus2amps_set_mpi_parameters(ThisThread,nTotalThreads,Communicator)
  use BATL_lib,  ONLY: iProc,nProc,iComm
  
  implicit none
  integer,intent(in)::ThisThread,nTotalThreads,Communicator
  
  iProc=ThisThread
  nProc=nTotalThreads
  iComm=Communicator
end subroutine batsrus2amps_set_mpi_parameters















  
  