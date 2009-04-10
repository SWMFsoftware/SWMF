program select_snapshot

  ! Read a *.outs IDL file of type ascii/real4/real8, select snapshots, 
  ! then write it into a new file, of any type of ascii/real4/real8 

  use ModIoUnit,    ONLY: io_unit_new
  use ModPlotFile,  ONLY: save_plot_file, read_plot_file
  use ModUtilities, ONLY: split_string

  implicit none

  character(len=100) :: NameFileIn
  character(len=100) :: NameFileOut
  character(len=10)  :: TypeFileIn
  character(len=10)  :: TypeFileOut
  character(len=10)  :: TypePosition
  character(len=100) :: StringSnapshot
  character(len=500) :: StringHeader
  character(len=500) :: NameVar
  integer            :: nStep, nDim, nParam, nVar
  real               :: Time
  real               :: Param_I(100)
  real,  allocatable :: Coord_DIII(:,:,:,:), Var_VIII(:,:,:,:)

  integer, parameter :: MaxSnapshot = 1000
  logical :: DoSave_I(MaxSnapshot)

  integer :: nSnapshot

  character(len=100) :: StringPart_I(MaxSnapshot)
  character(len=100) :: StringSubpart_I(3)

  integer :: n_D(0:3), n1, n2, n3 
  integer :: i, iPart, nPart, nSubpart
  integer :: iSnapshot, iSnapshot1, iSnapshot2, DiSnapshot, nSnapshotWritten

  integer :: iUnit
  character(len=*), parameter:: NameSub = 'select_snapshot'

  !-------------------------------------------------------------------------

  ! Name and type of the input file
  write(*,'(a)') 'Input filename:'
  read(*,'(a)') NameFileIn
  write(*,'(a)') 'Type of input file (ascii/real4/real8):'
  read(*,'(a)') TypeFileIn

  ! Input snapshots 
  write(*,*) '----------------------------------------------------------'  
  write(*,*) 'Example for snapshots selecting: 2, 3 , 6:12:2, 20:25'
  write(*,*) 'Please use "," to indicate snapshots as the example shows.'  
  write(*,*) 'Blanks are allowed. '   
  write(*,*) 'The maximum number of snapshots is 1000.'
  write(*,*) '----------------------------------------------------------'  
  write(*,'(a)') 'Select snapshots:'
  read(*,'(a)') StringSnapshot

  ! Input name and type of snapshot file
  write(*,'(a)') 'Output filename:'
  read(*,'(a)') NameFileOut 
  write(*,'(a)') 'Type of output file (ascii/real4/real8):'
  read(*,'(a)') TypeFileOut

  ! Get index of snapshots from StringSnapshot 

  DoSave_I = .false.
  iSnapshot = 0
  call split_string(StringSnapshot, MaxSnapshot, StringPart_I, nPart, ',')
  do iPart = 1, nPart
     call split_string(StringPart_I(iPart), 3, StringSubpart_I, nSubpart, ':')
     if(nSubpart == 1) StringSubpart_I(2) = StringSubpart_I(1)
     if(nSubpart < 3) StringSubpart_I(3) = '1'
     read(StringSubpart_I(1), *) iSnapshot1
     read(StringSubpart_I(2), *) iSnapshot2
     read(StringSubpart_I(3), *) DiSnapshot

     DoSave_I(iSnapshot1:iSnapshot2:DiSnapshot) = .true.
    
  end do


  ! Now we get the index of snapshots need to extract out.
  nSnapshot = count(DoSave_I)
  write(*,*) 'Number of snapshots selected:', nSnapshot


  ! Open a I/O unit for writing the new file 
  iUnit = io_unit_new()
  iSnapshot = 0  ! number of snapshots read from the *.outs file
  nSnapshotWritten = 0 ! number of snapshots have been written

  TypePosition = 'rewind'

  do i = 1, MaxSnapshot
     ! Read header info
     call read_plot_file(NameFileIn, iUnit, TypeFileIn, &
          StringHeader, nStep,  Time, nDim, nParamOut = nParam, &
          nVarOut = nVar, n1Out = n1, n2Out = n2, n3Out = n3, &
          nOut_D = n_D, ParamOut_I = Param_I, NameVarOut = NameVar)

     ! Determine the shape of arrays from the header
     allocate(Coord_DIII(nDim, n1, n2, n3), Var_VIII(nVar, n1, n2, n3))

     ! Read the coord and var arrays
     call read_plot_file(NameFileIn, iUnit, TypeFileIn, &
          CoordOut_DIII = Coord_DIII, VarOut_VIII = Var_VIII)

     iSnapshot = iSnapshot + 1

     ! Write into new file if the snapshot is selected
     if(DoSave_I(iSnapshot))then
        call save_plot_file(NameFileOut, TypePosition, TypeFileOut, &
             StringHeader, nStep, Time, Param_I(1:nParam), NameVar, &
             nDimIn = nDim, CoordIn_DIII = Coord_DIII, &
             VarIn_VIII = Var_VIII)
        nSnapshotWritten = nSnapshotWritten + 1
        TypePosition = 'append'
     end if

     deallocate(Coord_DIII, Var_VIII)

     if(nSnapshotWritten == nSnapshot) EXIT 
     !when all selected snapshots are written  

  end do

  close(iUnit)

end program select_snapshot

!=========================================================================

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError

  stop

end subroutine CON_stop
