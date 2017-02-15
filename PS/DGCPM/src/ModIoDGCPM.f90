!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModIoDGCPM
  !\
  ! Input/output variable definition module for the DGCPM program.
  ! Updated: Aron Dodger, January 2012 
  !/

  use ModSizeDGCPM
  use ModIoUnit, Only: iUnit => UNITTMP_, STDOUT_

  logical :: IsFramework

  ! Define a few time and geophysical index input variables
  ! Formerly: Common block PARAM2
  integer :: nstep,iwpi
  real    :: ut,r,time,tWriteOutput

  integer :: nst, nkp, nibc
  
  ! Define variables for setting Kp:
  character(len=5) :: NameSourceKp='const'
  real             :: kpConst=0.0, Kp=0.0
  
  ! Define some variables set with the input file
  ! Formerly: Common block PARAM3
  real :: TINJ
  character(len=5) ::  name

  ! Define convection input parameters
  ! Formerly: Common block PARAM5
  integer, parameter :: ilame=450, ippcm=5000
  real :: lamgam,lambe(ilame),tlame(ilame),ppc(ippcm),tppc(ippcm)

  integer, parameter        :: iCharLen_     = 100

  integer, parameter        :: nInputMaxLines = 10000
  integer                   :: nInputLines
  character (len=iCharLen_) :: cInputText(nInputMaxLines)

  character (len=iCharLen_) :: cInputFile = "input.dgcpm"
  character(len=*), parameter :: cOutputDir  = "PS/Output/"
  character(len=*), parameter :: cInputDir   = "PS/Input/"
  character(len=*), parameter :: cRestartIn  = "PS/restartIN/"
  character(len=*), parameter :: cRestartOut = "PS/restartOUT/"
  

  integer :: iUnitOut=STDOUT_
  integer, parameter :: lStringPrefix=6
  character (len=lStringPrefix) :: StringPrefix = '' 

  ! MLTSlice file parameters.
  logical :: DoMltSlice = .false.
  integer :: nMltSlice = 4
  real    :: DtMltSlice = 300.0
  integer, allocatable :: iUnitMlt(:)

  ! LSlice file parameters.
  integer :: iUnitSlice

contains

  !===========================================================================
  subroutine load_restart_file(NameFileIn)
    
    use ModIoUnit, ONLY: UnitTMP_  
    use ModSizeDGCPM, only: nthetacells, nphicells
    use ModMainDGCPM, only: vthetacells, vphicells, mgridden, &
         mgridx, mgridy, mgridoc, mgridn, mgridvol

    implicit none
    
    character(len=*), intent(IN) :: NameFileIn

    integer :: i, j, nThetaFile, nPhiFile
    real, allocatable, dimension(:)   :: vThetaFile, vPhiFile
    real, allocatable, dimension(:,:) :: mDenFile, mxFile, myFile, ocFile
    
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='load_restart_file'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Open restart file:
    if(DoTestMe) write(*,*) 'PS: Reading restart file ', NameFileIn
    open(file=NameFileIn, status='old', unit=UnitTMP_)

    ! Load the size of the grid within the file:
    read(UnitTMP_, *) nThetaFile, nPhiFile

    if(DoTestMe)then
       write(*,*)'PS: Restart file size (nTheta,nPhi)=', nThetaFile, nPhiFile
       write(*,*)'PS: Domain size (nTheta,nPhi)=',       nThetaCells,nPhiCells
    end if

    if ( (nThetaFile==nThetaCells).and.(nPhiFile==nPhiCells) ) then
       ! If restart grid same as current grid size,
       ! read directly into state variable arrays:
       if(DoTestMe) write(*,*)'PS: No restart file interpolation required'
       read(UnitTMP_,*) vthetacells
       read(UnitTMP_,*) vphicells
       read(UnitTMP_,*) mgridden
       read(UnitTMP_,*) mgridx
       read(UnitTMP_,*) mgridy
       read(UnitTMP_,*) mgridoc
    else
       ! If different grid sizes, buffer arrays are required.
       ! Allocate buffer arrays:
       if(DoTestMe) write(*,*)'PS: Restart file interpolation required'
       allocate(vThetaFile(nThetaFile), vPhiFile(nPhiFile))
       allocate(myFile(nThetaFile, nPhiFile), ocFile(nThetaFile, nPhiFile), &
            mDenFile(nThetaFile, nPhiFile), mxFile(nThetaFile, nPhiFile))

       ! Load data into buffer arrays:
       read(UnitTMP_,*) vThetaFile
       read(UnitTMP_,*) vPhiFile
       read(UnitTMP_,*) mDenFile
       read(UnitTMP_,*) mxFile
       read(UnitTMP_,*) myFile
       read(UnitTMP_,*) ocFile

       ! Interpolate from buffer to state arrays:
       call interpol2Dpolar(&
            vThetaFile, nThetaFile, vPhiFile, nPhiFile, mDenFile, &
            vthetacells, nthetacells, vphicells, nphicells, mgridden)
       call interpol2Dpolar(&
            vThetaFile, nThetaFile, vPhiFile, nPhiFile, mxFile, &
            vthetacells, nthetacells, vphicells, nphicells, mgridx)
       call interpol2Dpolar(&
            vThetaFile, nThetaFile, vPhiFile, nPhiFile, myFile, &
            vthetacells, nthetacells, vphicells, nphicells, mgridy)
       call interpol2Dpolar(&
            vThetaFile, nThetaFile, vPhiFile, nPhiFile, ocFile, &
            vthetacells, nthetacells, vphicells, nphicells, mgridoc)
     
       ! Deallocate buffer arrays:
       deallocate(myFile, ocFile, vThetaFile, vPhiFile, mDenFile, mxFile)
    end if

    ! Close restart file
    close(unit = UnitTMP_)

    ! Populate flux tube content (previously "Denton" subroutine):
    do i = 1, nthetacells
       do j = 1, nphicells
          if (mgridoc(i,j).gt.0.999) then
             mgridn(i,j) = mgridden(i,j) * mgridvol(i,j)
          endif
       enddo
    enddo
    
  end subroutine load_restart_file
  !===========================================================================
  subroutine write_prefix

    if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

  end subroutine write_prefix

end Module ModIoDGCPM
