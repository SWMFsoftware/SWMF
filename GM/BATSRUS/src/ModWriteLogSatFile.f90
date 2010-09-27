!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_logfile(iSatIn,iFile)
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY  : State_VGB, tmp1_BLK
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModPhysics, ONLY  : Si2Io_V, Si2No_V, UnitT_
  use ModIO
  use ModIoUnit, ONLY   : io_unit_new
  use ModUtilities, ONLY: flush_unit, split_string
  use ModSatelliteFile, ONLY: Satellite_Name, Satellite_First_Write, &
       iUnitSat_I, Sat_Time, Satellite_Vars, DoTrackSatellite_I, xSatellite
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: iFile

  ! iSatIn = 0 -> write logfile 
  ! iSatIn >=1 -> write satellite output files (number = isat)
  integer, intent(in) :: iSatIn

  ! Logfile variables
  integer, parameter :: MaxLogVar=40
  integer, parameter :: MaxLogR = 10
  real :: LogVar_I(MaxLogVar), LogVarSum_I(MaxLogVar)
  character (len=10) :: NameLogVar_I(MaxLogVar)
  character (len=10) :: NameLogR_I(MaxLogR)
  real :: LogR_I(MaxLogR)
  integer :: nLogVar, nLogR, nLogTot, i,j,iSat
  integer :: nFluxVar       ! number of flux variables used
  integer :: iUnit          ! local unit number
  real :: Xyz_D(3)
  character (LEN=255) :: NameAll
  character (LEN=100) :: StringTime

  logical :: DoWritePosition
  logical :: DoTest,DoTestMe
  real :: pmin, pmax
  real, external :: maxval_loc_BLK, minval_loc_BLK
  integer :: loc(5)
  integer :: iTime_I(7) ! integer time: year,month,day,hour,minute,sec,msec
  integer :: iError
  !---------------------------------------------------------------------------
  call set_oktest('write_logfile',DoTest,DoTestMe)

  DoWritePosition = .false.

  if(iSatIn == 0 .and. ( index(test_string,'show_pmin')>0 .or. &
       index(test_string,'show_pmax')>0 ) ) then

     tmp1_BLK(1:nI,1:nJ,1:nK,1:nBlock) = State_VGB(P_,1:nI,1:nJ,1:nK,1:nBlock)

     if(index(test_string,'show_pmin')>0)then
        pMin = minval_loc_BLK(nProc,tmp1_BLK,loc)
        if(loc(5)==iProc)write(*,*)'pmin, loc, x, y, z=',pmin,loc, &
             x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     end if

     if(index(test_string,'show_pmax')>0)then
        pMax = maxval_loc_BLK(nProc,tmp1_BLK,loc)
        if(loc(5)==iProc)write(*,*)'pmax, loc, x, y, z=',pmax,loc, &
             x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     end if
  end if

  nFluxVar = 0
  nLogR    = 0
  LogR_I   = 0.0
  if (iSatIn == 0) then
     call split_string(log_vars,MaxLogVar,NameLogVar_I,nLogVar)
     do i=1,nLogVar
        if(index(NameLogVar_I(i),'flx')>0) nFluxVar = nFluxVar + 1
        if(index(NameLogVar_I(i),'pnt')>0 .or.  &
             index(NameLogVar_I(i),'test')>0)  DoWritePosition = .true.
     end do
     if (nFluxVar >0) then
        call split_string(log_R_str, MaxLogR, NameLogR_I,nLogR)
        do i=1,nLogR
           read(NameLogR_I(i),*) LogR_I(i)
        end do
     end if
     StringTime = log_time
  elseif (iSatIn>=1) then
     iSat = iSatIn
     if (.not. DoTrackSatellite_I(iSat)) return
     call split_string(satellite_vars(satellite_+iSat),MaxLogVar, &
          NameLogVar_I,nLogVar)
     StringTime = sat_time(satellite_+iSat)
     DoWritePosition = .true.
  end if

  ! check to make sure that the total number of Logfile variables is smaller
  ! than MaxLogVar.  If it is not write a warning message and truncate the
  ! list.
  if (nLogVar + nFluxVar*nLogR > MaxLogVar) then
     write(*,*)'Warning in write_logfile: Number of logfile variables exceeds '
     write(*,*)'the array dimensions.  Truncating list - recompile with larger'
     write(*,*)'array dimensions'
     if (nLogVar >= MaxLogVar) then
        nLogVar = MaxLogVar
        nLogR = 1
     else
        nLogR = (MaxLogVar-nLogVar)/nFluxVar + 1
     end if
  end if
  nLogTot = nLogVar + nFluxVar*(nLogR-1)

  ! load the time variables to the list of output
  NameAll = ''
  if(index(StringTime,'none')<=0) then
     if(index(StringTime,'step')>0) NameAll = 'it'
     if(index(StringTime,'date')>0) NameAll = trim(NameAll)// &
          ' year mo dy hr mn sc msc'
     if(index(StringTime,'time')>0) NameAll = trim(NameAll)//' t'
  end if

  if (DoWritePosition) then
     NameAll = trim(NameAll)//' X Y Z'
     if (iSatIn == 0) then
        Xyz_D = XyzTestCell_D
     elseif (iSatIn >= 1) then
        Xyz_D = XSatellite(iSat,1:3)
     end if
  end if

  do i=1,nLogVar
     if(index(NameLogVar_I(i),'flx')>0) then
        do j=1,nLogR
           NameAll = trim(NameAll)//' '//trim(NameLogVar_I(i))//'_R='// &
                trim(NameLogR_I(j))
        end do
     else
        NameAll = trim(NameAll)//' '//trim(NameLogVar_I(i))
     end if
  end do

  if(DoTestMe.and.n_step==1)then
     write(*,*)'nLogVar,nFluxVar,nLogR,nLogTot:',  &
          nLogVar,nFluxVar,nLogR,nLogTot
     write(*,*)'NameLogVar_I:',NameLogVar_I(1:nLogVar)
     write(*,*)'NameLogR_I:',NameLogR_I(1:nLogR)
     write(*,*)'LogR_I:',LogR_I(1:nLogR)
     write(*,*)'NameAll:',trim(NameAll)
  end if

  if(iProc==0) then
     if (iSatIn==0) then
        if(unit_log<0)then
           unit_log = io_unit_new()
           if(n_step < 1000000)then
              write(filename,'(a,i6.6,a)')&
                   trim(NamePlotDir)//'log_n',n_step,'.log'
           else
              write(filename,'(a,i8.8,a)')&
                   trim(NamePlotDir)//'log_n',n_step,'.log'
           end if
           open(unit_log,file=filename,status="replace")
           if (index(NameAll,'pnt')>0 .or. index(NameAll,'test')>0) then
   	      if (coord_test) then
   		 write(unit_log,'(a,3es13.5)')  &
                      'test point (X,Y,Z): ',Xtest,Ytest,Ztest
   	      else
   		 write(unit_log,'(a,5(i6))')  &
                      'test point(I,J,K,BLK,PROC): ', &
   		      Itest,Jtest,Ktest,BLKtest,PROCtest
   	      end if
           else
              write(unit_log,'(a)')'Volume averages, fluxes, etc'
           end if
           write(unit_log,'(a)') trim(NameAll)
        endif
        iUnit = unit_log
     elseif (iSatIn >= 1) then
        iUnit = iUnitSat_I(iSat)
        if (Satellite_first_write(iSat)) then
           write(iUnit,'(a)')  &
                'Satellite data for Satellite: '//trim(Satellite_name(isat))
           write(iUnit,'(a)')trim(NameAll)
           Satellite_first_write(iSat)=.false.
        end if
     end if
  endif

  call set_logvar(nLogVar,NameLogVar_I,nLogR,LogR_I,nLogTot,LogVar_I,iSatIn)

  if(nProc > 0)then
     ! Collect LogVar_I from all processors
     call MPI_reduce(LogVar_I, LogVarSum_I, nLogTot, MPI_REAL, MPI_SUM, 0, &
          iComm, iError)
     if(iProc == 0)LogVar_I = LogVarSum_I
  end if


  ! WRITE OUT THE LINE INTO THE LOGFILE OR THE SATELLITE FILE
  if(iProc==0) then

     if (plot_dimensional(iFile))  &
          call normalize_logvar(nLogVar,NameLogVar_I,nLogR,LogR_I, &
          nLogTot,LogVar_I)

     ! first output the appropriate time data
     if(index(StringTime,'none')>0) then
        ! do nothing
     else
        if(index(StringTime,'step')>0) &
             write(iUnit,'(i7)',ADVANCE='NO') n_step
        if(index(StringTime,'date')>0) then
           call get_date_time(iTime_I)
           write(iUnit,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') &
                iTime_I
        end if
        if(index(StringTime,'time')>0) then
           !note that Time_Simulation is in SI units.
           if(plot_dimensional(iFile)) then
              write(iUnit,'(es13.5)',ADVANCE='NO') &
                   Time_Simulation*Si2Io_V(UnitT_)
           else
              write(iUnit,'(es13.5)',ADVANCE='NO') &
                   Time_Simulation*Si2No_V(UnitT_)
           end if
        end if
     end if

     ! Now write the position of the test point or satellite if desired
     if(DoWritePosition) &
          write(iUnit,'(3es13.5)',ADVANCE='NO') Xyz_D

     ! Finally write out the data variables
     write(iUnit,'(100es13.5)') LogVar_I(1:nLogTot)

     call flush_unit(iUnit)
  end if

end subroutine write_logfile

!==============================================================================
subroutine set_logvar(nLogVar,NameLogVar_I,nLogR,LogR_I,nLogTot,LogVar_I,iSat)

  use ModProcMH
  use ModNumConst
  use ModMPI
  use ModMain, ONLY: n_step,Dt,Cfl,unusedBLK,nI,nJ,nK,nBlock,UseUserLogFiles,&
       iTest,jTest,kTest,ProcTest,BlkTest,optimize_message_pass,x_,y_,&
       UseRotatingFrame,UseB0
  use ModPhysics,    ONLY: rCurrents, inv_gm1, OMEGABody
  use ModVarIndexes
  use ModAdvance,    ONLY: tmp1_BLK, tmp2_BLK, &
       B0_DGB, State_VGB, Energy_GBI, DivB1_GB
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,x1,x2,y1,y2,z1,z2, &
       DomainVolume
  use ModRaytrace,   ONLY: ray  !^CFG  IF RAYTRACE
  use ModSatelliteFile, ONLY: get_satellite_ray !^CFG  IF RAYTRACE
  use ModSatelliteFile, ONLY: xSatellite
  use ModIO
  use ModMultiFluid, ONLY: iRho, iRhoUx, iRhoUy, iRhoUz, iP, iFluid

  implicit none

  integer, intent(in)            :: nLogVar, nLogR, nLogTot, iSat
  character (LEN=10), intent(in) :: NameLogVar_I(nLogVar)
  real, intent(in)               :: LogR_I(nLogR)
  real, intent(out)              :: LogVar_I(nLogTot)

  character (len=10) :: NameLogVar

  real :: StateIntegral_V(nVar)
  real :: SatRayVar_I(5), SatRayVarSum_I(5)

  integer :: iVar,iR,iVarTot, iBLK
  integer :: i,j,k
  integer :: iError
  real :: r

  logical :: DoTest,DoTestMe

  ! StateSat_V contains the weight (0), the state (1:nVar), 
  ! and the currents (nVar+1:nVar+3) at the position of the satellite
  ! B0Sat_D contains B0 at the satellite
  real :: StateSat_V(0:nVar+3), B0Sat_D(3)

  real, external :: integrate_BLK, maxval_BLK, minval_BLK
  real, external :: calc_sphere,integrate_circle

  !-------------------------------------------------------------------------
  call set_oktest('set_logvar',DoTest,DoTestMe)

  if(DoTestMe.and.n_step==1)then
     write(*,*)'nLogVar,nLogR,nLogTot:',nLogVar,nLogR,nLogTot
     write(*,*)'NameLogVar_I:',NameLogVar_I(1:nLogVar)
     write(*,*)'LogR_I:',LogR_I(1:nLogR)
  end if

  LogVar_I(1:nLogTot) = 0.0
  tmp1_BLK=1.00
  DomainVolume  =integrate_BLK(nProc,tmp1_BLK)

  ! Obtain data to calculate log variables
  if(iSat>=1)then
     ! Satellites need B0 and the state at the satellite position
     if(UseB0)then
        call get_b0(xSatellite(iSat,1),xSatellite(iSat,2),xSatellite(iSat,3),&
             B0Sat_D)
     else
        B0Sat_D=0.00
     end if
     call get_point_data(0.0,xSatellite(iSat,:),1,nBlock,1,nVar+3,StateSat_V)
     call collect_satellite_data(xSatellite(iSat,:),StateSat_V)
     if (UseRotatingFrame) then
        StateSat_V(rhoUx_)=StateSat_V(rhoUx_) &
                - StateSat_V(rho_)*OMEGAbody*xSatellite(iSat,y_)
        StateSat_V(rhoUy_)=StateSat_V(rhoUy_) &
                + StateSat_V(rho_)*OMEGAbody*xSatellite(iSat,x_)
     end if

     !^CFG IF RAYTRACE BEGIN
     !If any ray tracing satellite variables are present, collect ray data
     do iVar=1, nLogVar
        select case(NameLogVar_I(iVar))
        case('theta1','theta2','phi1','phi2','status')
           call get_satellite_ray(iSat, SatRayVar_I)
           call MPI_reduce(SatRayVar_I, SatRayVarSum_I, 6, MPI_REAL, MPI_SUM, &
                0, iComm, iError)
           EXIT
        end select
     enddo
     !^CFG END RAYTRACE

  else
     ! The logfile may need the integral of conservative variables
     ! Also extract the pressure into a variable for pmin and pmax
     do iVar = 1, nLogVar
        call normalize_name_log_var(NameLogVar_I(iVar), NameLogVar)
        select case(NameLogVar)
        case('rho','rhoux','rhouy','rhouz','bx','by','bz','p','pmin','pmax')
           call integrate_domain(StateIntegral_V, tmp2_BLK)
           EXIT
        end select
     end do
  end if

  iVarTot = 0
  do iVar=1,nLogVar

     iVarTot = iVarTot+1
     call normalize_name_log_var(NameLogVar_I(iVar), NameLogVar)

     ! If we are a satellite and not a logfile (iSat>=1) then we should
     ! do only satellite variables so append a 'sat' to the end of the 
     ! variables
     if (iSat >= 1) then
        call set_sat_var
     else
        call set_log_var
     end if
  end do

contains
  !============================================================================
  subroutine set_log_var

    use ModMain,      ONLY: x_, y_, z_
    use ModUser,      ONLY: user_get_log_var
    use ModUtilities, ONLY: lower_case
    use ModCurrent,   ONLY: get_current
    use ModWaves,     ONLY: UseWavePressure

    ! Local variables
    real :: Bx, By, Bz, RhoUx, RhoUy, RhoUz, bDotB, bDotU, qval, qval_all
    real :: Current_D(3)
    real :: FullB_DG(3,0:nI+1,0:nJ+1,0:nK+1)

    integer :: jVar
    character(len=10) :: NameVar, NameLogVarLower

    ! External function for ionosphere    !^CFG IF IONOSPHERE
    real, external :: logvar_ionosphere   !^CFG IF IONOSPHERE
    !------------------------------------------------------------------------

    select case(NameLogVar)

!!$! MHD variables averaged over the computational domain
    case('e')
       LogVar_I(iVarTot) = &
            integrate_BLK(1,Energy_GBI(:,:,:,:,iFluid))/DomainVolume
    case('pmin')
       ! Divide by nProc so that adding up the processors can work
       LogVar_I(iVarTot) = minval_BLK(nProc,tmp2_BLK)/nProc
    case('pmax')
       ! Divide by nProc so that adding up the processors can work
       LogVar_I(iVarTot) = maxval_BLK(nProc,tmp2_BLK)/nProc
    case('ux')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUx,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('uy')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUy,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('uz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUz,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('ekinx')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUx,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('ekiny')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUy,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('ekinz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(iRhoUz,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/DomainVolume
    case('ekin')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               (State_VGB(iRhoUx,1:nI,1:nJ,1:nK,iBLK)**2+&
               State_VGB(iRhoUy,1:nI,1:nJ,1:nK,iBLK)**2+&
               State_VGB(iRhoUz,1:nI,1:nJ,1:nK,iBLK)**2)&
               /State_VGB(iRho,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/DomainVolume

    case('jin','jout','jinmax','joutmax')

       if(index(optimize_message_pass,'opt')>0) &
            call stop_mpi('Spherical integral of J requires '// &
            'message passing edges and corners. Fix PARAM.in!')

       ! calculate the total/maximum current either into or out of the body
       do iBLK = 1, nBlock
          if(unusedBLK(iBLK))cycle

          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             ! Calculate radial current
             call get_current(i,j,k,iBlk,Current_D)
             tmp1_BLK(i,j,k,iBLK) = ( &
                  Current_D(x_)*x_BLK(i,j,k,iBLK)+ &
                  Current_D(y_)*y_BLK(i,j,k,iBLK)+ &
                  Current_D(z_)*z_BLK(i,j,k,iBLK) ) / r_BLK(i,j,k,iBLK)
          end do; end do; end do

          !now modify tmp1 according to the case we want
          select case(NameLogVar)
          case('jin')
             tmp1_BLK(:,:,:,iBLK) = min(tmp1_BLK(:,:,:,iBLK), 0.0)
          case('jout')
             tmp1_BLK(:,:,:,iBLK) = max(tmp1_BLK(:,:,:,iBLK), 0.0)
          end select
       end do
       select case(NameLogVar)
       case('jin','jout')
          LogVar_I(iVarTot) = calc_sphere('integrate',180,rCurrents,tmp1_BLK) &
               /(4.0*cPi*rCurrents**2)
       case('jinmax')
          qval = calc_sphere('minval',180,rCurrents,tmp1_BLK)
          if(nProc>1)then
             call MPI_allreduce(qval, qval_all, 1, MPI_REAL, MPI_MIN, &
             iComm, iError)
             ! Divide by nProc so that adding up the processors can work
             LogVar_I(iVarTot) = qval_all/nProc
          else
             LogVar_I(iVarTot) = qval/nProc
          endif
       case('joutmax')
          qval = calc_sphere('maxval',180,rCurrents,tmp1_BLK)
          if(nProc>1)then
             call MPI_allreduce(qval, qval_all, 1, MPI_REAL, MPI_MAX, &
             iComm, iError)
             ! Divide by nProc so that adding up the processors can work
             LogVar_I(iVarTot) = qval_all/nProc
          else
             LogVar_I(iVarTot) = qval/nProc
          endif
       end select

    case('dst')
       ! Calculate the Biot-Savart formula for the center of the Earth:
       ! B = (1/4 Pi)*integral( (curl B) x R / R**3 dV)
       ! where R is the vector FROM the CURRENT to the FIELD POINT!
       ! Since the field point is at the origin, R = (-x,-y,-z)
       ! Only the Z component is calculated here: (J x R)_z = -J_x*y + Jy*x

       ! Calculate 
       do iBLK = 1, nBlock
          if(unusedBLK(iBLK))cycle           
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if ( r_BLK(i,j,k,iBLK) < rCurrents .or. &
                  x_BLK(i+1,j,k,iBLK) > x2 .or.      &
                  x_BLK(i-1,j,k,iBLK) < x1 .or.      &
                  y_BLK(i,j+1,k,iBLK) > y2 .or.      &
                  y_BLK(i,j-1,k,iBLK) < y1 .or.      &
                  z_BLK(i,j,k+1,iBLK) > z2 .or.      &
                  z_BLK(i,j,k-1,iBLK) < z1 ) then
                tmp1_BLK(i,j,k,iBLK)=0.0
                CYCLE
             end if
             call get_current(i,j,k,iBlk,Current_D)
             tmp1_BLK(i,j,k,iBLK) = ( &
                  -Current_D(x_)*y_BLK(i,j,k,iBLK) &
                  +Current_D(y_)*x_BLK(i,j,k,iBLK) ) / r_BLK(i,j,k,iBLK)**3
          end do; end do; end do
       end do
       ! The /4pi is part of the Biot-Savart formula
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK) / (4*cPi)

    case('dstdivb')
       ! Calculate the contribution of Div B to the surface integral of DST
       ! Error = (1/4 Pi)*integral( (div B) R / R**3 dV)
       ! Since the field point is at the origin, R = (-x,-y,-z)
       ! Only the Z component is calculated here: z * div B/R**3

       ! Calculate 
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle           
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             if ( r_BLK(i,j,k,iBLK) < rCurrents .or. &
                  x_BLK(i+1,j,k,iBLK) > x2 .or.      &
                  x_BLK(i-1,j,k,iBLK) < x1 .or.      &
                  y_BLK(i,j+1,k,iBLK) > y2 .or.      &
                  y_BLK(i,j-1,k,iBLK) < y1 .or.      &
                  z_BLK(i,j,k+1,iBLK) > z2 .or.      &
                  z_BLK(i,j,k-1,iBLK) < z1 ) then
                tmp1_BLK(i,j,k,iBLK)=0.0
                CYCLE
             end if
             tmp1_BLK(i,j,k,iBLK) = &
                  z_BLK(i,j,k,iBLK) * DivB1_GB(i,j,k,iBLK) &
                  / r_BLK(i,j,k,iBLK)**3
          end do; end do; end do
       end do
       ! The 4*pi is part of the Biot-Savart formula
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK) / (4*cPi)

!!$! MHD variables at Itest, Jtest, Ktest, BLKtest, PROCtest
    case('rhopnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(iRho,iTest,jTest,kTest,BlkTest)
    case('rhouxpnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(iRhoUx,iTest,jTest,kTest,BlkTest)
    case('rhouypnt')     
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(iRhoUy,iTest,jTest,kTest,BlkTest)
    case('rhouzpnt')      
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(iRhoUz,iTest,jTest,kTest,BlkTest)
    case('b1xpnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(Bx_,iTest,jTest,kTest,BlkTest)
    case('b1ypnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(By_,iTest,jTest,kTest,BlkTest)
    case('b1zpnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(Bz_,iTest,jTest,kTest,BlkTest)
    case('bxpnt')
       if(iProc == ProcTest)then 
          if(UseB0)then
             LogVar_I(iVarTot) = &
                  B0_DGB(x_,iTest,jTest,kTest,BlkTest) + &
                  State_VGB(Bx_,iTest,jTest,kTest,BlkTest)
          else
             LogVar_I(iVarTot) = &
                  State_VGB(Bx_,iTest,jTest,kTest,BlkTest)
          end if
       end if
    case('bypnt')                      
       if(iProc == ProcTest) then
          if(UseB0)then
             LogVar_I(iVarTot) = &
                  B0_DGB(y_,iTest,jTest,kTest,BlkTest) + &
                  State_VGB(By_,iTest,jTest,kTest,BlkTest)
          else
             LogVar_I(iVarTot) = &
                  State_VGB(By_,iTest,jTest,kTest,BlkTest)
          end if
       end if
    case('bzpnt')                      
       if(iProc == ProcTest) then
          if(UseB0)then
             LogVar_I(iVarTot) = &
                  B0_DGB(z_,iTest,jTest,kTest,BlkTest) + &
                  State_VGB(Bz_,iTest,jTest,kTest,BlkTest)
          else
             LogVar_I(iVarTot) = &
                  State_VGB(Bz_,iTest,jTest,kTest,BlkTest)
          end if
       end if
    case('ppnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(iP,iTest,jTest,kTest,BlkTest)
    case('epnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = Energy_GBI(iTest,jTest,kTest,BlkTest,iFluid)
    case('uxpnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(iRhoUx,iTest,jTest,kTest,BlkTest) / &
            State_VGB(iRho,  iTest,jTest,kTest,BlkTest)
    case('uypnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(iRhoUy,iTest,jTest,kTest,BlkTest) / &
            State_VGB(iRho,  iTest,jTest,kTest,BlkTest)
    case('uzpnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(iRhoUz,iTest,jTest,kTest,BlkTest) / &
            State_VGB(iRho,  iTest,jTest,kTest,BlkTest)
!!$!RAYTRACE variables                                ^CFG  IF RAYTRACE BEGIN
       ! RAYTRACE variables averaged over volume
    case('theta1', 'theta2', 'phi1', 'phi2','status')
       select case(NameLogvar)
       case('theta1')
          i = 1; j = 1
       case('theta2')
          i = 1; j = 2
       case('phi1')
          i = 2; j = 1
       case('phi2')
          i = 2; j = 2
       case('status')
          i = 3; j = 1
       end select
       do iBlk = 1, nBlock
          if(unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBlk) = ray(i,j,1:nI,1:nJ,1:nK,iBlk)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/DomainVolume
       ! RAYTRACE variables at Itest, Jtest, Ktest, BLKtest, PROCtest
    case('theta1pnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(1,1,iTest,jTest,kTest,BlkTest)
    case('theta2pnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(1,2,iTest,jTest,kTest,BlkTest)
    case('phi1pnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(2,1,iTest,jTest,kTest,BlkTest)
    case('phi2pnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(2,2,iTest,jTest,kTest,BlkTest)
    case('statuspnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(3,1,iTest,jTest,kTest,BlkTest)
!!$!^CFG END RAYTRACE

!!$! Ionosphere values                                 ^CFG IF IONOSPHERE BEGIN
    case('cpcpn','cpcp_n','cpcp_north','cpcpnorth',&
         'cpcps','cpcp_s','cpcp_south','cpcpsouth')
       ! It is sufficient to calculate it on processor 0
       if(iProc == 0) &
            LogVar_I(iVarTot) = logvar_ionosphere(NameLogVar) 
       !                                               ^CFG END IONOSPHERE

!!$! Flux values integrated over the surface of a shpere
    case('aflx')
       ! just to check that the area is being computed correctly
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,1:nBlock) = 1.0
          LogVar_I(iVarTot) = calc_sphere('integrate',50, r, tmp1_BLK)
       end do
    case('rhoflx')
       ! rho*U_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK)) CYCLE
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                tmp1_BLK(i,j,k,iBLK) = &
               (State_VGB(iRhoUx,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
               +State_VGB(iRhoUy,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
               +State_VGB(iRhoUz,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
               )/R_BLK(i,j,k,iBLK)
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = calc_sphere('integrate',360, r, tmp1_BLK)
       end do
    case('dstflx')
       ! B1z
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,1:nBlock) = &
               State_VGB(Bz_,0:nI+1,0:nJ+1,0:nK+1,1:nBlock)
          LogVar_I(iVarTot) = calc_sphere('integrate',180, r, tmp1_BLK) / &
               (4*cPi*r**2)
       end do
    case('bflx')
       ! B_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK)) CYCLE
             FullB_DG=State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             if(UseB0)FullB_DG = FullB_DG &
                  +B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                tmp1_BLK(i,j,k,iBLK) = ( &
                     FullB_DG(x_,i,j,k)*x_BLK(i,j,k,iBLK) + &
                     FullB_DG(y_,i,j,k)*y_BLK(i,j,k,iBLK) + &
                     FullB_DG(z_,i,j,k)*z_BLK(i,j,k,iBLK) &
                     ) / R_BLK(i,j,k,iBLK)
             end do; end do; end do
          end do
 
          LogVar_I(iVarTot) = calc_sphere('integrate',360, r, tmp1_BLK)
       end do
    case('b2flx')
       ! B^2*u_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)           
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))cycle
             FullB_DG = State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             if(UseB0)FullB_DG = FullB_DG &
                  +B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBLK)           
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                tmp1_BLK(i,j,k,iBLK) = &
                     sum(FullB_DG(:,i,j,k)**2)
                tmp1_BLK(i,j,k,iBLK) = 0.5*tmp1_BLK(i,j,k,iBLK)* &
                     ( State_VGB(iRhoUx,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                     + State_VGB(iRhoUy,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                     + State_VGB(iRhoUz,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                     ) / (State_VGB(iRho,i,j,k,iBLK)*R_BLK(i,j,k,iBLK))
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = calc_sphere('integrate',360, r,tmp1_BLK)
       end do
    case('pvecflx')
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          r = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))CYCLE
             FullB_DG=State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             if(UseB0)FullB_DG = FullB_DG &
                  +B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                RhoUx = State_VGB(iRhoUx,i,j,k,iBLk)
                RhoUy = State_VGB(iRhoUy,i,j,k,iBLk)
                RhoUz = State_VGB(iRhoUz,i,j,k,iBLk)
                bDotb = sum(FullB_DG(:,i,j,k)**2)
                bDotU = FullB_DG(x_,i,j,k)*RhoUx +&
                     FullB_DG(y_,i,j,k)*RhoUy + FullB_DG(z_,i,j,k)*RhoUz
                tmp1_BLK(i,j,k,iBLk) = ( &
                  ( bDotb*rhoUx - bDotU*FullB_DG(x_,i,j,k))*X_BLK(i,j,k,iBLk)+&
                  ( bDotb*rhoUy - bDotU*FullB_DG(y_,i,j,k))*Y_BLK(i,j,k,iBLk)+&
                  ( bDotb*rhoUz - bDotU*FullB_DG(z_,i,j,k))*Z_BLK(i,j,k,iBLk) &
                  ) / (State_VGB(iRho,i,j,k,iBLk)*R_BLK(i,j,k,iBLk))
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = calc_sphere('integrate',360, r,tmp1_BLK)
       end do

       ! simple circular integrals
    case('e2dflx')
       !this is the azimuthal component of the electric field 
       !integrated around a circle
       ! Ex*y - Ey*x
       iVarTot = iVarTot-1
       do iR=1,nLogR
          iVarTot = iVarTot+1
          r = LogR_I(iR)
          do iBLK = 1,nBlock
             if(unusedBLK(iBLK))CYCLE
             FullB_DG = State_VGB(Bx_:Bz_,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             if(UseB0)FullB_DG = FullB_DG + B0_DGB(:,0:nI+1,0:nJ+1,0:nK+1,iBLK)
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                Bx = FullB_DG(x_,i,j,k)
                By = FullB_DG(y_,i,j,k)
                Bz = FullB_DG(z_,i,j,k)
                RhoUx = State_VGB(iRhoUx,i,j,k,iBLk)
                RhoUy = State_VGB(iRhoUy,i,j,k,iBLk)
                RhoUz = State_VGB(iRhoUz,i,j,k,iBLk)
                tmp1_BLK(i,j,k,iBLK) = &
                     ( (RhoUy*Bz-RhoUz*By)*y_BLK(i,j,k,iBLK) &
                     - (RhoUz*Bx-RhoUx*Bz)*x_BLK(i,j,k,iBLK) &
                     ) / (State_VGB(iRho,i,j,k,iBLK)*R_BLK(i,j,k,iBLK)) 
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = integrate_circle(r, 0.0, tmp1_BLK)
       end do

       ! OTHER VALUES
    case('dt')
       if(iProc == 0)LogVar_I(iVarTot) = dt

    case('cfl')
       if(iProc == 0)LogVar_I(iVarTot) = Cfl

    case('ew', 'erad')
       if(Ew_ == 1)then
          if(UseWavePressure)then
             LogVar_I(iVarTot) = &
                  sum(StateIntegral_V(WaveFirst_:WaveLast_))/DomainVolume
          else
             LogVar_I(iVarTot) = 0.0
          end if
       else
          LogVar_I(iVarTot) = StateIntegral_V(Ew_)/DomainVolume
       end if

    case default
       ! Check if the variable name is one of the state variables
       NameLogVarLower = NameLogVar_I(iVar)
       call lower_case(NameLogVarLower)

       do jVar = 1, nVar
          NameVar = NameVar_V(jVar)
          call lower_case(NameVar)
          if(NameVar /= NameLogVarLower) CYCLE
          LogVar_I(iVarTot) = StateIntegral_V(jVar)/DomainVolume
          RETURN
       end do
       if (UseUserLogFiles) then
          if (index(NameLogVar,'flx')>0)then
             iVarTot = iVarTot - 1
             do iR=1,nLogR
                iVarTot = iVarTot + 1
                r = LogR_I(iR)           
                call user_get_log_var(LogVar_I(iVarTot), NameLogVarLower, r)
             end do
          else
             call user_get_log_var(LogVar_I(iVarTot), NameLogVarLower)
          end if
       else
          if(iProc == 0)then
             LogVar_I(iVarTot) = -7777.
             call write_myname;
             write(*,*) 'WARNING in set_logvar: unknown logvarname=',NameLogVar
          end if
       end if
    end select
  end subroutine set_log_var

  !==========================================================================
  subroutine set_sat_var
    use ModProcMH
    use ModNumConst
    use ModVarIndexes
    use ModIO
    use ModUtilities, ONLY: lower_case
    use ModPhysics, ONLY:  AverageIonCharge, ElectronTemperatureRatio
    use ModMultiFluid, ONLY: UseMultiIon, IsMhd, iFluid, iRho, iP, iRhoIon_I, &
         MassIon_I
    implicit none

    integer :: jVar, jFluid
    character(len=10) :: NameVar, NameLogVarLower

    !-------------------------------------------------------------------------
    if (iProc/=0) RETURN

    select case(NameLogVar)
    case('rho')
       LogVar_I(iVarTot) = StateSat_V(IRho)
    case('rhoux')
       LogVar_I(iVarTot) = StateSat_V(iRhoUx)
    case('rhouy')
       LogVar_I(iVarTot) = StateSat_V(iRhoUy)
    case('rhouz')
       LogVar_I(iVarTot) = StateSat_V(iRhoUz)
    case('bx')
       LogVar_I(iVarTot) = StateSat_V(Bx_) + B0Sat_D(1)
    case('by')
       LogVar_I(iVarTot) = StateSat_V(By_) + B0Sat_D(2)
    case('bz')
       LogVar_I(iVarTot) = StateSat_V(Bz_) + B0Sat_D(3)
    case('e')
       LogVar_I(iVarTot) = inv_gm1*StateSat_V(iP) + 0.5*&
            sum(StateSat_V(iRhoUx:iRhoUz)**2)/StateSat_V(iRho)
       if(iFluid == 1 .and. IsMhd) &
            LogVar_I(iVarTot) = LogVar_I(iVarTot) &
            + 0.5*sum((StateSat_V(Bx_:Bz_) + B0Sat_D)**2)
    case('n','t','temp')
       ! Calculate the number density
       if(UseMultiSpecies)then
          LogVar_I(iVarTot)=0.0
          do jVar = SpeciesFirst_, SpeciesLast_
             LogVar_I(iVarTot) = LogVar_I(iVarTot) + &
                  StateSat_V(jVar)/MassSpecies_V(jVar)
          end do
       else if(UseMultiIon .and. IsMhd .and. iFluid==1)then
          LogVar_I(iVarTot) = sum(StateSat_V(iRhoIon_I)/MassIon_I)
       else
          LogVar_I(iVarTot) = StateSat_V(iRho)/MassFluid_I(iFluid)
       end if
       
       ! Calculate temperature from P = n*k*T + ne*k*Te = n*k*T*(1+ne/n*Te/T)
       if(NameLogVar /= 'n') LogVar_I(iVarTot) = &
            StateSat_V(iP) / LogVar_I(iVarTot) &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
    case('p')
       LogVar_I(iVarTot) = StateSat_V(iP)
    case('pperp')
       LogVar_I(iVarTot) = (3*StateSat_V(iP)-StateSat_V(Ppar_))/2.0
    case('ux')
       LogVar_I(iVarTot) = StateSat_V(iRhoUx)/StateSat_V(iRho)
    case('uy')
       LogVar_I(iVarTot) = StateSat_V(iRhoUy)/StateSat_V(iRho)
    case('uz')
       LogVar_I(iVarTot) = StateSat_V(iRhoUz)/StateSat_V(iRho)
    case('b0x')
       LogVar_I(iVarTot) = B0Sat_D(1)
    case('b0y')
       LogVar_I(iVarTot) = B0Sat_D(2)
    case('b0z')
       LogVar_I(iVarTot) = B0Sat_D(3)
    case('b1x')
       LogVar_I(iVarTot) = StateSat_V(Bx_)
    case('b1y')
       LogVar_I(iVarTot) = StateSat_V(By_)
    case('b1z')
       LogVar_I(iVarTot) = StateSat_V(Bz_)
    case('jx')
       LogVar_I(iVarTot) = StateSat_V(nVar+1)
    case('jy')
       LogVar_I(iVarTot) = StateSat_V(nVar+2)
    case('jz')
       LogVar_I(iVarTot) = StateSat_V(nVar+3)
    case('weight')
       ! Total weight of a point (normall 1.0)
       LogVar_I(iVarTot) = StateSat_V(0)
    case('order')
       ! Order of accuracy (2 if total weight is 1.0)
       if(abs(StateSat_V(0)-1)<1.e-5)then
          LogVar_I(iVarTot) = 2
       else
          LogVar_I(iVarTot) = 1
       end if
       
       !Raytracing footpoint values               !^CFG IF RAYTRACE BEGIN
    case('theta1')
       LogVar_I(iVarTot) = SatRayVarSum_I(1)
    case('phi1')
       LogVar_I(iVarTot) = SatRayVarSum_I(2)
    case('status')
       LogVar_I(iVarTot) = SatRayVarSum_I(3)
    case('theta2')
       LogVar_I(iVarTot) = SatRayVarSum_I(4)
    case('phi2')
       LogVar_I(iVarTot) = SatRayVarSum_I(5)      !^CFG END RAYTRACE

    case default
       ! Check if the variable name is one of the state variables
       NameLogVarLower = NameLogVar_I(iVar)
       call lower_case(NameLogVarLower)
       do jVar = 1, nVar
          NameVar = NameVar_V(jVar)
          call lower_case(NameVar)
          if(NameVar /= NameLogVarLower) CYCLE
          LogVar_I(iVarTot) = StateSat_V(jVar)
          RETURN
       end do
       LogVar_I(iVarTot) = -777.0
       if(iProc==0)write(*,*)'WARNING in var_sat: unknown variable ',&
            NameLogVar,' for iSat = ',iSat
    end select
  end subroutine set_sat_var
  !===========================================================================

end subroutine set_logvar

!=============================================================================
subroutine normalize_logvar(nLogVar,NameLogVar_I,nLogR,&
     LogR_I,nLogTot,LogVar_I)

  use ModPhysics
  use ModVarIndexes, ONLY: NameVar_V, UnitUser_V
  use ModUtilities, ONLY: lower_case
  implicit none

  integer, intent(in) :: nLogVar, nLogR, nLogTot
  character (LEN=10), intent(in) :: NameLogVar_I(nLogVar)
  real, intent(inout) :: LogVar_I(nLogTot)
  real, intent(in) :: LogR_I(nLogR)

  character (len=10) :: NameLogVar, NameVar
  integer :: iVar, iVarTot, jVar
  !-------------------------------------------------------------------------

  iVarTot = 0
  do iVar=1,nLogVar

     iVarTot = iVarTot+1
     NameLogVar = NameLogVar_I(iVar)
     call normalize_name_log_var(NameLogVar_I(iVar), NameLogVar)

     select case(NameLogVar)

!!$! BASIC MHD variables
     case('rho','rhopnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*No2Io_V(UnitRho_)
     case('rhoux','rhouy','rhouz', 'rhouxpnt','rhouypnt','rhouzpnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*No2Io_V(UnitRhoU_)
     case('bx','by','bz','bxpnt','bypnt','bzpnt','b1xpnt','b1ypnt','b1zpnt', &
          'b1x','b1y','b1z','b0x','b0y','b0z','dst','dstdivb')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*No2Io_V(UnitB_)
     case('e','epnt','ew','erad')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*No2Io_V(UnitEnergyDens_)
     case('p','ppnt','pmin','pmax','pperp')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*No2Io_V(UnitP_)
     case('ux','uy','uz','uxpnt','uypnt','uzpnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*No2Io_V(UnitU_)
     case('ekinx','ekiny','ekinz','ekin')
        LogVar_I(iVarTot)= LogVar_I(iVarTot) &
             *No2Io_V(UnitRho_)*No2Io_V(UnitU_)**2
     case('jx','jy','jz','jxpnt','jypnt','jzpnt',&
          'jin','jout','jinmax','joutmax')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*No2Io_V(UnitJ_)

     case('n')
        LogVar_I(iVarTot)=LogVar_I(iVarTot)*No2Io_V(UnitN_)
     case('t','temp')
        LogVar_I(iVarTot)=LogVar_I(iVarTot)*No2Io_V(UnitTemperature_)

!!$! Ionosphere values                                
        !^CFG IF IONOSPHERE BEGIN
     case('cpcpn','cpcp_n','cpcp_north','cpcpnorth',&
          'cpcps','cpcp_s','cpcp_south','cpcpsouth')
        ! User unit is kV = 1000 V
        LogVar_I(iVarTot) = LogVar_I(iVarTot) &
             *(No2Si_V(UnitElectric_)*No2Si_V(UnitX_))/1000.0       
        !^CFG END IONOSPHERE

!!$! Flux values
     case('aflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *No2Io_V(UnitX_)**2
        iVarTot = iVarTot+nLogR-1
     case('rhoflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(No2Si_V(UnitRho_)*No2Si_V(UnitU_)*No2Si_V(UnitX_)**2)
        iVarTot = iVarTot+nLogR-1
     case('dstflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *No2Io_V(UnitB_)
        iVarTot = iVarTot+nLogR-1
     case('bflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(No2Si_V(UnitB_)*No2Si_V(UnitX_)**2)
        iVarTot = iVarTot+nLogR-1
     case('b2flx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(No2Si_V(UnitPoynting_)*No2Si_V(UnitX_)**2)
        iVarTot = iVarTot+nLogR-1
     case('pvecflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(No2Si_V(UnitPoynting_)*No2Si_V(UnitX_)**2)
        iVarTot = iVarTot+nLogR-1
     case('e2dflx') 
        ! circular integral of the azimuthal component of the electric field
        ! on the surface of a sphere.
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *No2Si_V(UnitElectric_)*No2Si_V(UnitX_)
        iVarTot = iVarTot+nLogR-1

        ! OTHER VALUES
     case('dt')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*No2Io_V(UnitT_)

     case default
        do jVar = 1, nVar
           NameVar = NameVar_V(jVar)
           call lower_case(NameVar)
           if(NameLogVar /= NameVar) CYCLE
           LogVar_I(iVarTot)=LogVar_I(iVarTot)*UnitUser_V(jVar)
           EXIT
        end do
        ! no normalization
     end select
  end do ! iVar
end subroutine normalize_logvar

!==============================================================================

real function calc_sphere(TypeAction,nTheta,Radius,Array_GB)

  ! This function calculates the integral of the incomming variable Array_GB
  ! over the surface of a sphere centered at the origin radius Radius.  
  ! The resolution in the colatitude is determined by the nTheta parameter.

  use ModMain,           ONLY: nI,nJ,nK,nBLK,nBlock,unusedBLK,&
       optimize_message_pass
  use ModGeometry,       ONLY: x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK, &
       r_BLK, XyzStart_Blk
  use ModNumConst
  use ModCovariant, ONLY: TypeGeometry             
  use ModInterpolate, ONLY: trilinear
  implicit none

  ! Arguments

  character (len=*), intent(in) :: TypeAction
  integer, intent(in)           :: nTheta
  real, intent(in)              :: Radius 
  real, intent(in)              :: Array_GB(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK)

  ! Local variables
  real :: Result, Darea0, Darea, Average 

  ! Indices and coordinates
  integer :: iBlock,i,j,k,i1,i2
  integer :: MaxPhi, nPhi
  real    :: x, y, z, InvDxyz_D(3)
  real    :: Dr
  real    :: xMin, xMax, yMin, yMax, zMin, zMax, rMin, rMax
  real    :: dTheta, dPhi, Phi, Theta, SinTheta

  real :: Array_G(0:nI+1,0:nJ+1,0:nK+1)

  ! Store cartesian coordinates for sake of efficiency
  ! The x and y depend on iPhi,iTheta while z only depends on iTheta
  !  real, allocatable :: x_II(:,:), y_II(:,:), z_I(:), SinTheta_I(:)

  logical :: DoTest=.false.,DoTestMe=.false.
  !---------------------------------------------------------------------------

  call set_oktest('calc_sphere',DoTest,DoTestMe)
  call timing_start('calc_sphere_'//TypeAction)

  ! Initialize result
  select case(TypeAction)
  case('integrate')
     Result = 0.0
  case('maxval')
     Result = -Huge(0.0)
  case('minval')
     Result = Huge(0.0)
  case default
     call stop_mpi('ERROR in calc_sphere: Invalid action='//TypeAction)
  end select

  if(index(TypeGeometry,'spherical') > 0)then         
     ! For spherical geometry it is sufficient to 
     ! interpolate in the radial direction

     do iBlock = 1, nBlock
        if (unusedBLK(iBlock)) CYCLE
        rMin = cHalf*(R_BLK( 0, 1, 1,iBlock) + R_BLK( 1, 1, 1,iBlock))
        if(rMin > Radius) CYCLE
        rMax = cHalf*(R_BLK(NI, 1, 1,iBlock) + R_BLK(NI+1,1,1,iBlock))
        if(rMax <= Radius) CYCLE

        ! Set temporary array
        Array_G = Array_GB(0:nI+1,0:nJ+1,0:nK+1,iBlock)

        dTheta = dz_BLK(iBlock); dPhi=dy_BLK(iBlock)
        dArea0 = Radius**2 *dPhi *dTheta

        ! Find the radial index just after Radius
        i2=0
        do while ( Radius > R_BLK( i2, 1, 1, iBlock))
           i2 = i2+1
        end do
        i1=i2-1

        Dr = (Radius - R_BLK(i1, 1, 1, iBlock)) &
             /(R_BLK(i2, 1, 1, iBlock) - R_BLK(i1, 1, 1, iBlock))
        if(Dr<0.or.Dr>1)call stop_mpi('wrong index in calc_sphere')

        select case(TypeAction)
        case('integrate')
           ! Integrate in theta
           do k = 1, nK
              SinTheta = sqrt(x_BLK(1,1,k,iBlock)**2 + y_BLK(1,1,k,iBlock)**2)&
                   /r_BLK(1,1,k,iBlock)
              ! Integrate in Phi by adding up 1..nJ and interpolate in R
              Average = Dr * sum( Array_G(i2, 1:nJ, k)) + &
                   (1-Dr)  * sum( Array_G(i1, 1:nJ, k))
              Result = Result + dArea0*SinTheta * Average
           end do
        case('maxval')
           Average = maxval(   Dr * Array_G(i2, 1:nJ, 1:nK)  &
                +           (1-Dr)* Array_G(i1, 1:nJ, 1:nK))
           Result = max(Average, Result)
        case('minval')
           Average = minval(   Dr * Array_G(i2, 1:nJ, 1:nK)  &
                +           (1-Dr)* Array_G(i1, 1:nJ, 1:nK))
           Result = min(Average, Result)
        end select
     end do
  elseif(TypeGeometry /= 'cartesian')then
     call stop_mpi('ERROR in calc_sphere: Not implemented for geometry=' &
          //TypeGeometry)
  else                                                     
     ! Get the angular resolution from the input parameter nTheta
     MaxPhi = 2*nTheta
     dTheta = cPi/nTheta
     dArea0 = Radius**2 * 2 * sin(cHalf*dTheta)

     if (DoTestMe) write(*,*) 'nTheta,MaxPhi,dTheta[deg]:',nTheta,MaxPhi,&
          dTheta*cRadToDeg

     ! Calculate sin(theta) and x,y,z in advance
     !  allocate( x_II(MaxPhi,nTheta), y_II(MaxPhi, nTheta), z_I(nTheta), &
     !       SinTheta_I(nTheta) )
     !
     !  do i = 1, nTheta
     !     Theta         = (i - cHalf)*dTheta
     !     SinTheta      = sin(Theta)
     !     SinTheta_I(i) = SinTheta
     !     z_I(i)        = Radius*cos(Theta)
     !
     !     ! Number of Phi coordinates is proportional to 2*nTheta*SinTheta
     !     ! This keeps the area of the spherical cells roughly constant 
     !     ! and the shape roughly a square.
     !     ! Make sure that nPhi is a multiple of 4 to preserve symmetry
     !     nPhi  = min(MaxPhi, 4 * ceiling(cQuarter*MaxPhi*SinTheta))
     !     Dphi  = cTwoPi/nPhi
     !     do j = 1, nPhi
     !        Phi = j*dPhi
     !        ! Convert to Cartesian coordinates
     !        x_II(j,i) = Radius*SinTheta*cos(Phi)
     !        y_II(j,i) = Radius*SinTheta*sin(Phi)
     !     end do
     !  end do

     ! Sum all cells within range
     do iBlock = 1, nBlock

        if (unusedBLK(iBlock)) CYCLE

        ! get the max and min radial distance for this block so that 
        ! we can check whether or not this block contibutes to the sum.

        xMin = cHalf*(x_BLK( 0, 0, 0,iBlock) + x_BLK(   1,   1  , 1,iBlock))
        xMax = cHalf*(x_BLK(nI,nJ,nK,iBlock) + x_BLK(nI+1,nJ+1,nK+1,iBlock))
        yMin = cHalf*(y_BLK( 0, 0, 0,iBlock) + y_BLK(   1,   1,   1,iBlock))
        yMax = cHalf*(y_BLK(nI,nJ,nK,iBlock) + y_BLK(nI+1,nJ+1,nK+1,iBlock))
        zMin = cHalf*(z_BLK( 0, 0, 0,iBlock) + z_BLK(   1,   1,   1,iBlock))
        zMax = cHalf*(z_BLK(nI,nJ,nK,iBlock) + z_BLK(nI+1,nJ+1,nK+1,iBlock))

        if( minmod(xMin,xMax)**2+minmod(yMin,yMax)**2+minmod(zMin,zMax)**2 &
             > Radius**2) &
             CYCLE
        if( maxmod(xMin,xMax)**2+maxmod(yMin,yMax)**2+maxmod(zMin,zMax)**2 &
             < Radius**2) &
             CYCLE

        InvDxyz_D = 1.0/(/ Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock) /)

        ! Set temporary array
        Array_G = Array_GB(0:nI+1,0:nJ+1,0:nK+1,iBlock)

        ! Fill in edges and corners for the first layer so that bilinear 
        ! interpolation can be used without message passing these values

        if(index(optimize_message_pass,'opt')>0) call fill_edge_corner(Array_G)

        do i = 1, nTheta

           ! Check if z is inside the block
           Theta = dTheta*(i-cHalf)
           z     = Radius*cos(Theta)
           if(z <  zMin) CYCLE
           if(z >= zMax) CYCLE

           SinTheta = sin(Theta)

           ! Number of Phi coordinates is proportional to 2*nTheta*SinTheta
           nPhi = min(MaxPhi, 4 * ceiling(0.25*MaxPhi*SinTheta))
           dPhi = cTwoPi/nPhi

           ! Area of the surface element
           dArea = dArea0 * SinTheta * dPhi

           do j = 1, nPhi
              Phi = j*dPhi

              ! Check if x and y are inside the block
              x = Radius * SinTheta * cos(Phi)
              if(x <  xMin) CYCLE
              if(x >= xMax) CYCLE

              y = Radius * SinTheta * sin(Phi)
              if(y <  yMin) CYCLE
              if(y >= yMax) CYCLE

              ! Compute the interpolated values at the current location.
              ! Coordinates are normalized so that index=coordinate. 
              ! XyzStart corresponds to 1,1,1 so we have to add 1 to the index.

              Average = trilinear( Array_G, 0,nI+1, 0,nJ+1, 0,nK+1, &
                   1 + InvDxyz_D*((/ x, y, z /) - XyzStart_Blk(:,iBlock)) )

              select case(TypeAction)
              case('integrate')
                 Result = Result + dArea * Average
              case('maxval')
                 Result = max(Result, Average)
              case('minval')
                 Result = min(Result, Average)
              end select
              !if(iBlock==222.and.i==22.and.j==76)then
              !   write(*,*)'j, Result=',j, Result
              !   write(*,*)'x,y,z=',x,y,z
              !   write(*,*)'XyzStart_Blk=',XyzStart_Blk(:,iBlock)
              !   write(*,*)'InvDxyz_D=',InvDxyz_D
              !end if
           end do
        end do
     end do
  end if                                           
  ! deallocate(x_II, y_II, z_I, SinTheta_I)

  calc_sphere = Result
  call timing_stop('calc_sphere_'//TypeAction)

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

  real function maxmod(x,y)
    real, intent(in) :: x,y
    maxmod = max(abs(x),abs(y))
  end function maxmod

end function calc_sphere


!==============================================================================

real function integrate_circle(Radius,z,Array_GB)

  ! This function calculates the integral of the incomming variable Array_GB
  ! over a cirle parallel to the equitorial plane.  The radius of the circle
  ! for the z axis is defined by the radius Radius and the z position is
  ! is given by z.  

  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlock,unusedBLK,optimize_message_pass
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,Dx_BLK,Dy_BLK,Dz_BLK,XyzStart_Blk
  use ModNumConst
  use ModInterpolate, ONLY: trilinear
  implicit none

  ! Arguments

  real, intent(in) :: Array_GB(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK)
  real, intent(in) :: Radius, z 

  ! Local variables
  real :: Integral, Average 

  ! Indices and coordinates
  integer :: iBlock,i
  integer :: nPhi
  real :: xMin,xMax,yMin,yMax,zMin,zMax
  real :: x, y, InvDxyz_D(3)
  real :: dPhi,Phi
  real :: Array_G(0:nI+1,0:nJ+1,0:nK+1)

  logical :: DoTest,DoTestMe
  !---------------------------------------------------------------------------

  Integral = 0.0

  call set_oktest('integrate_circle',DoTest,DoTestMe)

  ! the angular resolution of the integral is hard coded
  nPhi = 720
  dPhi = cTwoPi/nPhi

  if (DoTestMe) write(*,*) 'nPhi,dPhi',nPhi,dPhi

  ! Sum all cells within range

  do iBlock = 1, nBlock

     if (unusedBLK(iBlock)) CYCLE
     ! get the max and min radial (cylindrical) distance for this block so 
     ! that we can check whether or not this block contibutes to the sum.

     zMin = 0.5*(z_BLK( 0, 0, 0,iBlock)+z_BLK(   1,   1,   1,iBlock))
     zMax = 0.5*(z_BLK(nI,nJ,nK,iBlock)+z_BLK(nI+1,nJ+1,nK+1,iBlock))

     if(z < zMin .or. z >= zMax ) CYCLE

     xMin = 0.5*(x_BLK( 0, 0, 0,iBlock)+x_BLK(   1,   1  , 1,iBlock))
     xMax = 0.5*(x_BLK(nI,nJ,nK,iBlock)+x_BLK(nI+1,nJ+1,nK+1,iBlock))
     yMin = 0.5*(y_BLK( 0, 0, 0,iBlock)+y_BLK(   1,   1,   1,iBlock))
     yMax = 0.5*(y_BLK(nI,nJ,nK,iBlock)+y_BLK(nI+1,nJ+1,nK+1,iBlock))

     if( minmod(xMin,xMax)**2 + minmod(yMin,yMax)**2 > Radius**2) CYCLE
     if( maxmod(xMin,xMax)**2 + maxmod(yMin,yMax)**2 < Radius**2) CYCLE

     Array_G = Array_GB(0:nI+1,0:nJ+1,0:nK+1,iBlock)

     if(index(optimize_message_pass,'opt')>0) call fill_edge_corner(Array_G)

     InvDxyz_D = 1.0/(/ Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock) /)

     do i = 1, nPhi
        Phi = (i-0.5)*dPhi

        ! get the xyz coordinates
        x = Radius*cos(Phi)
        if( x < xMin .or. x >= xMax) CYCLE

        y = Radius*sin(Phi)
        if( y < yMin .or. y >= yMax) CYCLE

        ! Compute the interpolated values at the current location.
        ! Coordinates are normalized so that index=coordinate. 
        ! XyzStart corresponds to 1,1,1 so we have to add 1 to the index.

        Average = trilinear( Array_GB(:,:,:,iBlock),0,nI+1,0,nJ+1,0,nK+1, &
             1 + InvDxyz_D*((/ x, y, z /) - XyzStart_Blk(:,iBlock)))

        Integral = Integral + Average
     end do
  end do

  ! Now multiply by the size of each interval in the integral.  
  ! Since they are all the same we can do this after the fact 
  ! and not inside the above loops

  integrate_circle = Integral*Radius*dPhi

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

  real function maxmod(x,y)
    real, intent(in) :: x,y
    maxmod = max(abs(x),abs(y))
  end function maxmod

end function integrate_circle

!==============================================================================

subroutine collect_satellite_data(Xyz_D, StateCurrent_V)

  use ModProcMH, ONLY: nProc, iProc, iComm
  use ModVarIndexes, ONLY : nVar
  use ModMpi
  implicit none

  !INPUT ARGUMENTS:
  real, intent(in) :: Xyz_D(3) ! The position of the interpolated state

  !INPUT/OUTPUT ARGUMENTS:
  ! On input StateCurrent_V contains the weight and the interpolated state 
  ! on a given PE.
  ! On output only PE 0 contains new data.
  ! In the first 0th element the total weight is returned.
  ! If weight is 1.0 then the returned state is second order accurate
  ! If weight is positive but not 1.0 then the returned state is first order
  ! If weight is 0.0 (or less?) then the point was not found and the 
  !     returned state is -777.0
  ! The rest of the elements contain the globally interpolated state.

  real, intent(inout) :: StateCurrent_V(0:nVar+3)

  !LOCAL VARIABLES:
  ! This is needed for MPI_reduce
  real :: StateCurrentAll_V(0:nVar+3)

  ! Temporary variables
  real    :: Weight
  integer :: iError
  !---------------------------------------------------------------------------
  ! Collect contributions from all the processors to PE 0
  if(nProc>1)then
     call MPI_reduce(StateCurrent_V, StateCurrentAll_V, nVar+4,&
          MPI_REAL, MPI_SUM, 0, iComm, iError)
     if(iProc==0)StateCurrent_V = StateCurrentAll_V
  end if

  ! Check total weight and divide by it if necessary
  if(iProc==0)then
     Weight = StateCurrent_V(0)
     if(Weight<=0.0)then
        write(*,*)'collect_satellite_data WARNING total weight =',&
             Weight,' at Xyz_D=',Xyz_D
        StateCurrent_V = -777.0
     elseif(abs(Weight - 1) > 1.e-5)then
        StateCurrent_V(1:nVar+3) = StateCurrent_V(1:nVar+3) / Weight
     end if
  end if

end subroutine collect_satellite_data

!==============================================================================

subroutine satellite_test

  use ModProcMH, ONLY: iProc
  use ModVarIndexes
  use ModMain,     ONLY: nBlock,xTest,yTest,zTest
  use ModAdvance,  ONLY: State_VGB
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  implicit none
  real :: State_V(0:nVar+3)

  State_VGB(Bx_,:,:,:,:) = y_BLK
  State_VGB(By_,:,:,:,:) = z_BLK
  State_VGB(Bz_,:,:,:,:) = x_BLK

  call get_point_data(0.0,(/xTest,yTest,zTest/),1,nBlock,1,nVar+3,State_V)
  call collect_satellite_data((/xTest,yTest,zTest/),State_V)

  if(iProc==0)then
     if(max(abs(State_V(Bx_)-yTest),abs(State_V(By_)-zTest),&
          abs(State_V(Bz_)-xTest)) > 1e-7)then
        write(*,*)'Satellite state=',State_V(1:nVar)
        write(*,*)'Weight, Difference=', State_V(0), &
             State_V(Bx_)-yTest,State_V(By_)-zTest,State_V(Bz_)-xTest
     else
        write(*,*)'Satellite state is correct'
     end if
     if(maxval(abs(State_V(nVar+1:nVar+3)+1)) > 1e-7)then
        write(*,*)'Satellite current (should be -1)=',State_V(nVar+1:nVar+3)
     else
        write(*,*)'Satellite current is correct'
     end if
  end if

end subroutine satellite_test

!=============================================================================

subroutine normalize_name_log_var(NameIn, NameOut)

  ! Normalize the logvar name to lower case and 
  ! replace alternative names with a unique name

  use ModUtilities, ONLY: lower_case
  use ModMultiFluid,ONLY: extract_fluid_name
  implicit none

  character(len=*), intent(in)  :: NameIn
  character(len=*), intent(out) :: NameOut
  integer :: l
  !---------------------------------------------------------------------------
  NameOut = NameIn

  ! Switch to all lower case
  call lower_case(NameOut)

  call extract_fluid_name(NameOut)

  ! Replace mx with rhoux, my with rhouy, mz with rhouz
  select case(NameOut(1:2))
  case('mx','my','mz')
     NameOut = 'rhou'//NameOut(2:len(NameOut)-3)
  end select

  ! Replace *test with *pnt
  l = len_trim(NameOut)
  if(l<4) RETURN
  if(NameOut(l-3:l) == 'test') NameOut(l-3:l)='pnt '

end subroutine normalize_name_log_var

!=============================================================================
subroutine integrate_domain(Sum_V, Pressure_GB)

  ! Since the State_VGB variable has the variables in the first index,
  ! it is not efficient to integrate the variables one by one.
  ! This subroutine integrates all the state variables at the same time
  ! for the processor only.

  use ModProcMH
  use ModAdvance,   ONLY: State_VGB, P_
  use ModMain,      ONLY: nI, nJ, nK, nBlock, MaxBlock, UnusedBLK
  use ModVarIndexes,ONLY: nVar
  use ModGeometry,  ONLY: cV_BLK, true_BLK, true_cell 
  use ModGeometry,  ONLY: UseCovariant                
  use ModGeometry,  ONLY: 
  use ModNumConst
  implicit none 

  ! Arguments
  real, intent(out) :: Sum_V(nVar)
  real, intent(out) :: Pressure_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)

  ! Local variables:
  integer :: iBlock, iVar

  logical :: DoTest, DoTestMe

  !---------------------------------------------------------------------------
  if(UseCovariant)then                            
     call integrate_domain_covar(Sum_V, Pressure_GB)
     return
  end if                                          
  call set_oktest('integrate_domain',DoTest, DoTestMe)
  call timing_start('int_domain')
  Sum_V = cZero
                                                     
  do iBlock = 1, nBlock
     if(UnusedBLK(iBlock)) CYCLE
     if(true_BLK(iBlock)) then
        do iVar=1,nVar
           Sum_V(iVar) = Sum_V(iVar) + &
                sum(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))*cV_BLK(iBlock)
        end do
     else
        do iVar=1,nVar
           Sum_V(iVar) = Sum_V(iVar) + sum(&
                State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))*cV_BLK(iBlock)
        end do
     end if
     Pressure_GB(1:nI,1:nJ,1:nK,iBlock) = State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)
  end do
  call timing_stop('int_domain')
end subroutine integrate_domain

