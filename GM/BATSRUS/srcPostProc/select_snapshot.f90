program IDL_Snapshot

  ! Read a *.outs IDL file of type ascii/real4/real8, select snapshots, 
  ! then write it into a new file, of any type of ascii/real4/real8 

  use ModIoUnit,   ONLY: io_unit_new
  use ModPlotFIle, ONLY: save_plot_file, read_plot_file

  implicit none

  character(len=100) :: NameFileIn
  character(len=100) :: NameFileOut
  character(len=10)  :: TypeFileIn
  character(len=10)  :: TypeFileOut
  character(len=100) :: StringSnapshot
  character(len=500) :: StringHeader
  character(len=500) :: NameVar
  integer            :: nStep, nDim, nParam, nVar
  real               :: Time
  real               :: Param_I(100)
  real,  allocatable :: Coord_DIII(:,:,:,:), Var_VIII(:,:,:,:)

  integer :: n_D(0:3), n1, n2, n3 
  integer :: i, i1, i2, n
  integer :: iSnapshot, iSnapshot1, iSnapshot2, DiSnapshot, nSnapshotWritten

  integer, parameter :: MaxSnapshot = 1000
  integer :: iSnapshot_I(MaxSnapshot), nSnapshot

  integer :: iUnit
  character(len=*), parameter:: NameSub = 'IDL_Snapshot'

  !-------------------------------------------------------------------------

  ! Name and type of the input file
  write(*,'(a)') 'Input filename:'
  read(*,'(a)') NameFileIn
  write(*,'(a)') 'Type of input file (ascii/real4/real8):'
  read(*,'(a)') TypeFileIn

  ! Input snapshots 
  write(*,*) '----------------------------------------------------------'  
  write(*,*) 'Example for snapshots selecting: 2; 3 ; 6:12:2; 20:25 ;'
  write(*,*) 'Please use ";" to indicate snapshots as the example shows.'  
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

  i = 1   ! index of the input string
  iSnapshot = 0   ! index of the iSnapshot_I array

  do while(i <= len_trim(StringSnapshot))
     if((StringSnapshot(i:i) >= '0' .and. StringSnapshot(i:i) <= '9') &
          .or. StringSnapshot(i:i) == ' ')then
        i = i + 1

        ! The cases of "x:y(:z)" and "x;" (as x:x:1)
     else if(StringSnapshot(i:i) == ':' .or. StringSnapshot(i:i) == ';' &
          .and. i /= 1)then
        iSnapshot = iSnapshot + 1
        ! Skip blanks before [:] or [;]
        i1 = 0
        do while(StringSnapshot((i-1):(i-1)) == ' ')
           i = i - 1
           i1 = i1 + 1
        end do     ! back i1 position till the number at i-1
        ! Get the beginning index of the snapshot
        i2 = i - 1
        do while(StringSnapshot(i2:i2) >= '0' .and. &
             StringSnapshot(i2:i2) <= '9')          
           i2 = i2 - 1
        end do
        read(StringSnapshot((i2+1):(i-1)), *) iSnapshot1  !beginning
        iSnapshot_I(iSnapshot) = iSnapshot1

        ! go forward until [;] or [:]
        i = i + i1
        if(StringSnapshot(i:i) == ';')then  !x; = x:x:1
           iSnapshot2 = iSnapshot1
           DiSnapshot = 1
        else                                !x:y(:z)
           ! skip [:]
           i = i + 1
           ! Skip blanks after [:]
           do while(StringSnapshot(i:i) == ' ')
              i = i + 1
           end do
           ! Get the ending index of snapshots
           i2 = i
           do while(StringSnapshot(i2:i2) >= '0' .and. &
                StringSnapshot(i2:i2) <= '9')          
              i2 = i2 + 1
           end do
           read(StringSnapshot(i:(i2-1)), *) iSnapshot2     !ending

           i = i2  ! to the current position
           ! Skip blanks before [:] or [;]
           do while(StringSnapshot(i:i) == ' ')
              i = i + 1
           end do

           if(StringSnapshot(i:i) == ':')then  !x:y:z
              i = i + 1
              ! Skip blanks after [:]
              do while(StringSnapshot(i:i) == ' ')
                 i = i + 1
              end do
              ! Get z
              i2 = i
              do while(StringSnapshot(i2:i2) >= '0' .and. &
                   StringSnapshot(i2:i2) <= '9')          
                 i2 = i2 + 1
              end do
              read(StringSnapshot(i:(i2-1)), *) DiSnapshot  !interval
              i = i2  
              ! Skip blanks after the number
              do while(StringSnapshot(i:i) == ' ')
                 i = i + 1
              end do

           else if(StringSnapshot(i:i) == ';')then  !x:y = x:y:1
              DiSnapshot = 1
           else
              call CON_stop('cannot identify n1:n2(:n3) at' &
                   //StringSnapshot(i:i))
           end if
        end if
        ! Get the index to fill in iSnapshot_I
        do while(iSnapshot_I(iSnapshot) < iSnapshot2)
           iSnapshot = iSnapshot + 1
           iSnapshot_I(iSnapshot) = iSnapshot_I(iSnapshot-1) + DiSnapshot
           if(iSnapshot_I(iSnapshot) > iSnapshot2)then
              iSnapshot = iSnapshot - 1
              EXIT    !when the calculated index larger than the ending
           end if
        end do
        ! Skip [;]
        i = i + 1 
     else
        call CON_stop('cannot identify the character ' &
             //StringSnapshot(i:i))
     end if
  end do

  ! Now we get the index of snapshots need to extract out. 
  nSnapshot = iSnapshot
  write(*,*)  nSnapshot, 'snapshots selected:', iSnapshot_I(1:nSnapshot)


  ! Open a I/O unit for writing the new file 
  iUnit = io_unit_new()
  n = 0  ! number of snapshots read from the *.outs file
  nSnapshotWritten = 0 ! number of snapshots have been written

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

     n = n + 1
     ! Write into new file if the snapshot is selected
     do iSnapshot = 1, nSnapshot
        if(iSnapshot_I(iSnapshot) == n) then
           call save_plot_file(NameFileOut, 'append', TypeFileOut, &
                StringHeader, nStep, Time, Param_I(1:nParam), NameVar, &
                nDimIn = nDim, CoordIn_DIII = Coord_DIII, VarIn_VIII = Var_VIII)
           nSnapshotWritten = nSnapshotWritten + 1
        end if
     end do

     deallocate(Coord_DIII, Var_VIII)

     if(nSnapshotWritten == nSnapshot) EXIT 
     !when all selected snapshots are written  

  end do

  close(iUnit)

end program IDL_Snapshot

!=========================================================================

subroutine CON_stop(StringError)

  implicit none

  character (len=*), intent(in) :: StringError

  write(*,'(a)') 'ERROR: '//StringError

  stop

end subroutine CON_stop
