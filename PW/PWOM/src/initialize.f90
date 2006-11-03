subroutine PW_initialize

  use ModIoUnit, ONLY: io_unit_new,UnitTmp_
  use ModPwom
  implicit none

  ! Temporary variables
  real:: ddt1, xxx
  integer:: ns
  !---------------------------------------------------------------------------
  !***************************************************************************
  !  Set the number of fieldlines that each processor solves for
  !***************************************************************************
  if (iproc < mod(MaxLine,nProc)) then
     nLine=int(ceiling(real(MaxLine)/real(nProc)))
  else
     nLine=int(floor(real(MaxLine)/real(nProc)))
  endif

  write(*,*)'!!! iProc, nProc, nLine=',iProc, nProc, nLine

  !**************************************************************************
  !  Define file names and unit numbers, and open for reading and writing.
  !***************************************************************************
  NameSourceGraphics = 'PW/plot_sources.out'
  NameCollision      = 'PW/plots_collision.out'
  NamePhiNorth       = 'PW/North.dat'
  NamePhiSouth       = 'PW/South.dat'

  do iLine=1,nLine
     if (iproc .lt. mod(MaxLine,nProc)) then
        iLineGlobal(iLine)=&
             iproc*ceiling(real(MaxLine)/real(nProc))+iLine
     else
        iLineGlobal(iLine)=&
             (mod(MaxLine,nProc))*ceiling(real(MaxLine)/real(nProc)) &
             + ((iproc)-mod(MaxLine,nProc))                        &
             *floor(real(MaxLine)/real(nProc))+iLine
     endif
     write(NameRestartIn(iLine),"(a,i4.4,a)") &
          'PW/restartIN/restart_iline',iLineGlobal(iLine),'.dat'
     write(NameRestart(iLine),"(a,i4.4,a)") &
          'PW/restartOUT/restart_iline',iLineGlobal(iLine),'.dat'

     write(NameGraphics(iLine),"(a,i4.4,a)") &
          'PW/plots/plots_iline',iLineGlobal(iLine),'.out'

     iUnitGraphics(iLine)  = io_unit_new()
     open(iUnitGraphics(iLine),FILE=NameGraphics(iLine))
  enddo

!******************************************************************************
!  Read the restart file
!******************************************************************************
  nDim= 390
  if(IsRestart)then

     do iLine=1,nLine
        OPEN(UNIT=UnitTmp_, FILE=NameRestartIn(iLine), STATUS='OLD')
        READ (UnitTmp_,*) TIME,DDT1,NS
        READ (UnitTmp_,*)    GeoMagLat(iLine),GeoMagLon(iLine)
        
        READ (UnitTmp_,*) &
             (XXX,uOxyg(i,iLine),pOxyg(i,iLine),dOxyg(i,iLine),TOxyg(i,iLine),&
             i=1,NDIM)
        READ (UnitTmp_,*) &
             (XXX,uHel(i,iLine),pHel(i,iLine),dHel(i,iLine),THel(i,iLine),    &
             i=1,NDIM)
        READ (UnitTmp_,*) &
             (XXX,uHyd(i,iLine),pHyd(i,iLine),dHyd(i,iLine),THyd(i,iLine),    &
             i=1,NDIM)
        READ (UnitTmp_,*) &
             (XXX,uElect(i,iLine),pElect(i,iLine),dElect(i,iLine),            &
             TElect(i,iLine),i=1,NDIM)
        CLOSE(UNIT=UnitTmp_)
     enddo
  else
     Time=0.0
  endif

  !****************************************************************************
  ! Use Get_GITM to bring in neutral atmosphere from GITM
  !****************************************************************************
  !call GetNeutralData

  !****************************************************************************
  !  Set parameters for reading in potential and time of simulation
  !****************************************************************************
  Dtheta  = 0.0242
  Dphi    = 0.0245

  Dt      =    50.0
  !  Dt      =    0.2
  !  Dt      =    Tmax
  !maxTime = 10000.0
  maxTime = Tmax
  Time    =     0.0

  !****************************************************************************
  ! Read information from IE file, and get the velocities
  !****************************************************************************

  call get_electrodynamicPW

  !****************************************************************************
  !  Move flux tube around
  !****************************************************************************

  !initialize field line locations

  call initial_line_location

end subroutine PW_initialize
