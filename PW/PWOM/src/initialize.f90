subroutine PW_initialize

  use ModNumConst, ONLY: cDegToRad
  use ModMpi
  use ModIoUnit, ONLY: io_unit_new,UnitTmp_
  use ModPwom
  use ModCommonPlanet,ONLY: nIon,iRho_I,iU_I,iP_I,iT_I
  implicit none

  ! Temporary variables
  real:: ddt1, xxx
  integer:: ns, iPe, iError,iIon
  !---------------------------------------------------------------------------
  !***************************************************************************
  !  Set the number of fieldlines that each processor solves for
  !***************************************************************************
  if (iProc < mod(nTotalLine,nProc)) then
     nLine= (nTotalLine+nProc-1)/nProc
  else
     nLine= nTotalLine/nProc
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
     
     !Setup log files
     if (nLog == -1) then
        write(NameOutput(iLine),"(a,i4.4,a)") &
             'PW/log_iline',iLineGlobal(iLine),'.out'   
        iUnitOutput(iLine)  = io_unit_new()
        open(iUnitOutput(iLine),FILE=NameOutput(iLine))    
     elseif(nLog ==0) then
        !do nothing in this case
     elseif(nLog==iLineGlobal(iLine)) then
        write(NameOutput(iLine),"(a,i4.4,a)") &
             'PW/log_iline',iLineGlobal(iLine),'.out'      
        iUnitOutput(iLine)  = io_unit_new()
        open(iUnitOutput(iLine),FILE=NameOutput(iLine))     
     else
     end if

     
     
  enddo

!******************************************************************************
!  Read the restart file
!******************************************************************************

  if(IsRestart)then

     do iLine=1,nLine
        OPEN(UNIT=UnitTmp_, FILE=NameRestartIn(iLine), STATUS='OLD')
        READ (UnitTmp_,*) TIME,DDT1,nStep
        READ (UnitTmp_,*) GeoMagLat_I(iLine),GeoMagLon_I(iLine)

        ThetaLine_I (iLine) = (90.0-GeoMagLat_I(iLine)) * cDegToRad
        PhiLine_I   (iLine) = GeoMagLon_I(iLine)        * cDegToRad
        
        do iIon=1,nIon
           READ (UnitTmp_,*) &
                (XXX,State_CVI(i,iU_I(iIon),iLine),State_CVI(i,iP_I(iIon),iLine),&
                State_CVI(i,iRho_I(iIon),iLine),State_CVI(i,iT_I(iIon),iLine),&
                i=1,nAlt)
        enddo
        
        CLOSE(UNIT=UnitTmp_)
     enddo
  else
     do iLine = 1, nLine
        ThetaLine_I (iLine) = 10.0 * cDegToRad
        PhiLine_I   (iLine) = 0.0
     end do
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

  !****************************************************************************
  ! Read information from IE file, and get the velocities
  !****************************************************************************
  
  if(.not.UseIe)then
     call PW_get_electrodynamics

     !initialize field line locations
     call initial_line_location
  end if

end subroutine PW_initialize
