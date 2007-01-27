subroutine PW_initialize

  use ModMpi
  use ModIoUnit, ONLY: io_unit_new,UnitTmp_
  use ModPwom
  implicit none

  ! Temporary variables
  real:: ddt1, xxx
  integer:: ns, iPe, iError
  !---------------------------------------------------------------------------
  !***************************************************************************
  !  Set the number of fieldlines that each processor solves for
  !***************************************************************************
  if (iProc < mod(NTotalLine,nProc)) then
     nLine=int(ceiling(real(nTotalLine)/real(nProc)))
  else
     nLine=int(floor(real(nTotalLine)/real(nProc)))
  endif

  allocate(nLine_P(0:nProc-1), nLineBefore_P(0:nProc-1))
  call MPI_allgather(nLine,1,MPI_INTEGER, nLine_P,1,MPI_INTEGER,iComm,iError)
  nLineBefore_P(0) = 0
  do iPe = 1, nProc - 1
     nLineBefore_P(iPe) = sum(nLine_P(0:iPe-1))
  end do

  !**************************************************************************
  !  Define file names and unit numbers, and open for reading and writing.
  !***************************************************************************
  NameSourceGraphics = 'PW/plot_sources.out'
  NameCollision      = 'PW/plots_collision.out'
  NamePhiNorth       = 'PW/North.dat'
  NamePhiSouth       = 'PW/South.dat'

  do iLine=1,nLine
     if (iproc .lt. mod(nTotalLine,nProc)) then
        iLineGlobal(iLine)=&
             iproc*ceiling(real(nTotalLine)/real(nProc))+iLine
     else
        iLineGlobal(iLine)=&
             (mod(NTotalLine,nProc))*ceiling(real(nTotalLine)/real(nProc)) &
             + ((iproc)-mod(nTotalLine,nProc))                        &
             *floor(real(nTotalLine)/real(nProc))+iLine
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

  if(IsRestart)then

     do iLine=1,nLine
        OPEN(UNIT=UnitTmp_, FILE=NameRestartIn(iLine), STATUS='OLD')
        READ (UnitTmp_,*) TIME,DDT1,NS
        READ (UnitTmp_,*)    GeoMagLat_I(iLine),GeoMagLon_I(iLine)
        
        READ (UnitTmp_,*) &
             (XXX,uOxyg_CI(i,iLine),pOxyg_CI(i,iLine),dOxyg_CI(i,iLine),TOxyg(i,iLine),&
             i=1,nAlt)
        READ (UnitTmp_,*) &
             (XXX,uHel_CI(i,iLine),pHel_CI(i,iLine),dHel_CI(i,iLine),THel(i,iLine),    &
             i=1,nAlt)
        READ (UnitTmp_,*) &
             (XXX,uHyd_CI(i,iLine),pHyd_CI(i,iLine),dHyd_CI(i,iLine),THyd(i,iLine),    &
             i=1,nAlt)
        READ (UnitTmp_,*) &
             (XXX,uElect_CI(i,iLine),pElect_CI(i,iLine),dElect_CI(i,iLine),            &
             TElect(i,iLine),i=1,nAlt)
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
 
  TimeMax = Tmax
  Time    =     0.0

  !****************************************************************************
  ! Read information from IE file, and get the velocities
  !****************************************************************************
  
  if(.not.UseIe)then
     call PW_get_electrodynamics

     !initialize field line locations
     call initial_line_location
  end if

end subroutine PW_initialize
