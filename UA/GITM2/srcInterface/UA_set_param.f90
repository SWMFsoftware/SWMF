! Wrapper for GITM
!==============================================================================
subroutine UA_set_param(CompInfo, TypeAction)

  use ModInputs, only: cInputText
  use ModReadParam, only: read_text, n_line_read

  use ModTime, ONLY: StartTime, tSimulation, CurrentTime
  use ModInputs, only: iStartTime, IsFramework, iOutputUnit_, set_defaults, &
       nInputLines
  use ModTimeConvert, ONLY: time_real_to_int
  use CON_physics,    ONLY: get_time
  use ModIoUnit
  use ModProcUA
  use ModGITM, only: iCommGITM, nProcs, iProcGITM => iProc
  use ModPlanet, only: init_planet
  use CON_comp_info
  use ModUtilities, ONLY: check_dir

  implicit none

  character (len=*), parameter :: NameSub='UA_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: iError

  iError = 0

  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')

     call put(CompInfo,&
          Use=.true.,                                      &
          NameVersion='Global Iono-Thermo Model (Ridley)', &
          Version=2.0)

  case('MPI')

     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

     iCommGITM = iComm
     iProcGITM = iProc
     nProcs    = nProc

     if(iProc==0)then
        call check_dir("UA/DataIn")
        call check_dir("UA/data")
        call check_dir("UA/RestartOUT")
     end if

     IsFramework = .true.

     call init_planet
     call set_defaults

  case('READ')

     call read_text(cInputText)
     cInputText(n_line_read()+1) = "#END"
     nInputLines=n_line_read()+1

     call set_inputs

  case('CHECK')

     call check_param(iError)

  case('STDOUT')

!     iUnitStdOut=STDOUT_
!     if(nProc==1)then
!        StringPrefix='UA:'
!     else
!        write(StringPrefix,'(a,i3.3,a)')'UA',iProc,':'
!     end if

  case('FILEOUT')

     call get(CompInfo,iUnitOut=iOutputUnit_)
!     StringPrefix=''

  case('GRID')

     call UA_set_grid

  case default

     call CON_stop(NameSub//' UA_ERROR: invalid TypeAction='//TypeAction)

  end select

  if (iError /= 0) &
       call CON_stop(NameSub//' UA_ERROR in TypeAction='//TypeAction)

end subroutine UA_set_param

!=============================================================================

subroutine UA_set_grid

  ! Set the grid descriptor for UA
  ! Since UA has a static grid the descriptor has to be set once.
  ! There can be many couplers that attempt to set the descriptor,
  ! so we must check IsInitialized.
  use ModProcUA
  use CON_Coupler
  use CON_comp_info
  use ModNumConst
  use ModSizeGitm
  use ModSphereInterface, only: iStartBLK
  use ModInputs, only: nBlocksLat, nBlocksLon
  use ModUtilities, ONLY: check_allocate

  implicit none

  character (len=*), parameter :: NameSub='UA_set_grid'
  logical :: IsInitialized=.false.
  integer, allocatable :: iProc_A(:), iProcPE_A(:)

  integer :: iError, iBlock, iBlockPE

  real, allocatable :: CoLat_I(:), Lon_I(:), Alt_I(:)
  real, allocatable :: LatPE_I(:), LonPE_I(:)

  logical :: DoTest, DoTestMe, Done

  !------------------------------------------------------
  !    call CON_set_do_test(NameSub,DoTest, DoTestMe)
  !    if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
  if(IsInitialized) return

  IsInitialized=.true.

  if(iProc>=0)then
     allocate(CoLat_I(nBlocksLat*nLats), &
          Lon_I(nBlocksLon*nLons), &
          iProc_A(nBlocksLat*nBlocksLon), &
          LatPE_I(nBlocksLat*nBlocksLon*nLats), &
          LonPE_I(nBlocksLat*nBlocksLon*nLons), &
          iProcPE_A(nBlocksLat*nBlocksLon),&
          Alt_I(-1:nAlts+2),&
          stat = iError)

     if (iError /= 0) then
        write(*,*) NameSub, " Error in allocating variables"
        write(*,*) " Lat_I, Lon_I, iProc_A, LatPE_I, LonPE_I, iProcPE_A"
        call CON_stop(NameSub//' UA_ERROR')
     endif

     LatPE_I   = -1.0e32
     LonPE_I   = -1.0e32
     iProcPE_A = -1

     do iBlockPE = 1, nBlocks

        iBlock = iStartBLK + iBlockPE

!        LatPE_I((iBlock-1)*nLats+1:iBlock*nLats) = &
!             Latitude(1:nLats,iBlockPE)
!        LonPE_I((iBlock-1)*nLons+1:iBlock*nLons) = &
!             Longitude(1:nLons,iBlockPE)
        iProcPE_A(iBlock) = iProc

     enddo

!     call MPI_allreduce( LatPE_I, CoLat_I, nBlocksLon*nBlocksLat*nLats, &
!          MPI_REAL, MPI_MAX, iComm, iError)
!     ! Save into colatitudes instead of latitude
!     CoLat_I = cHalfPi - CoLat_I
!
!     call MPI_allreduce( LonPE_I, Lon_I, nBlocksLon*nBlocksLat*nLons, &
!          MPI_REAL, MPI_MAX, iComm, iError)

     call MPI_allreduce( iProcPE_A, iProc_A, nBlocksLon*nBlocksLat, &
          MPI_INTEGER, MPI_MAX, iComm, iError)
!     Alt_I=Altitude(:)
  else
     allocate( CoLat_I(1), Lon_I(1),iProc_A(1),Alt_I(1),stat=iError)
     call check_allocate(iError,NameSub)
  end if

  call set_grid_descriptor(                        &
       UA_,                                        &! component index
       nDim=3,                                     &! dimensionality
       nRootBlock_D=(/nBlocksLat,nBlocksLon,1/),     &! blocks
       nCell_D =(/nLats,nLons,nAlts/),             &! size of node based grid
       XyzMin_D=(/cHalf,cHalf,cHalf/),                   &! generalize coord
       XyzMax_D=(/nLats-cHalf,nLons-cHalf,nAlts-cHalf/), &! generalize coord
       TypeCoord='GEO',                            &! magnetic coordinates
!       Coord1_I= CoLat_I,                          &! colatitudes
!       Coord2_I= Lon_I,                            &! longitudes
!       Coord3_I= Alt_I,                            &! radial size in meters
       iProc_A = iProc_A,                          &! processor assigment
       IsPeriodic_D=(/.false.,.true.,.false./))     ! periodic in longitude

end subroutine UA_set_grid


!\
! -----------------------------------------------------------------------
!/

subroutine check_param(iError)

  implicit none

  integer, intent(out) :: iError

  iError = 0

end subroutine check_param
