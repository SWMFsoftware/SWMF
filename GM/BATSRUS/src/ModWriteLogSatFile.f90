!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_logfile(iSatIn,iFile)
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY  : State_VGB, tmp1_BLK
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModPhysics, ONLY  : unitUSER_t,unitSI_t
  use ModIO
  use ModIoUnit, ONLY   : io_unit_new
  use ModUtilities, ONLY: flush_unit
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
  logical :: oktest,oktest_me
  real :: pmin, pmax
  real, external :: maxval_loc_BLK, minval_loc_BLK
  integer :: loc(5)
  integer :: iTime_I(7) ! integer time: year,month,day,hour,minute,sec,msec
  integer :: iError
  !---------------------------------------------------------------------------
  call set_oktest('write_logfile',oktest,oktest_me)

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
     call split_str(log_vars,MaxLogVar,NameLogVar_I,nLogVar)
     do i=1,nLogVar
        if(index(NameLogVar_I(i),'flx')>0) nFluxVar = nFluxVar + 1
        if(index(NameLogVar_I(i),'pnt')>0 .or.  &
             index(NameLogVar_I(i),'test')>0)  DoWritePosition = .true.
     end do
     if (nFluxVar >0) then
        call split_str(log_R_str, MaxLogR, NameLogR_I,nLogR)
        do i=1,nLogR
           read(NameLogR_I(i),*) LogR_I(i)
        end do
     end if
     StringTime = log_time
  elseif (iSatIn>=1) then
     iSat = iSatIn
     if (.not. DoTrackSatellite_I(iSat)) return
     call split_str(satellite_vars(satellite_+iSat),MaxLogVar, &
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
        Xyz_D = (/ xTest_mod, yTest_mod, zTest_mod /); 
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

  if(oktest_me.and.n_step==1)then
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
           write(filename,'(a,i6.6,a)')&
                trim(NamePlotDir)//'log_n',n_step,'.log'
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

  LogVar_I(1:nLogTot) = 0.0

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
           if(plot_dimensional(iFile)) then
              !note that the funny (unitUSER_t/unitSI_t) is because by default
              !Time_Simulation is in seconds(SI).  If the user wants output
              !in some other time units this is the conversion.  In most 
              !cases the two are the same the and (...)=1.
              write(iUnit,'(es13.5)',ADVANCE='NO') &
                   Time_Simulation*(unitUSER_t/unitSI_t)
           else
              write(iUnit,'(es13.5)',ADVANCE='NO') &
                   Time_Simulation/unitSI_t
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
  use ModMain, ONLY: n_step,dt,unusedBLK,nI,nJ,nK,nBlock,gcn,UseUserLogFiles,&
       iTest,jTest,kTest,ProcTest,BlkTest
  use ModPhysics,    ONLY: rCurrents, rBody
  use ModVarIndexes
  use ModAdvance,    ONLY: tmp1_BLK, tmp2_BLK, &
       B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, State_VGB, E_BLK, DivB1_GB
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
       x1,x2,y1,y2,z1,z2
  use ModRaytrace,   ONLY: ray  !^CFG  IF RAYTRACE
  use ModIO
  implicit none

  integer, intent(in) :: nLogVar, nLogR, nLogTot, iSat
  character (LEN=10), intent(in) :: NameLogVar_I(nLogVar)
  real, intent(out) :: LogVar_I(nLogTot)
  real, intent(in) :: LogR_I(nLogR)
  character (len=10) :: NameLogVar

  real :: volume
  real,dimension(nVar) :: StateIntegral_V
  real, external :: integrate_BLK, maxval_BLK, minval_BLK
  real, external :: integrate_sphere, integrate_flux_circ,test_cell_value

  integer :: iVar,iR,iVarTot,itmp,jtmp, iBLK
  integer :: i,j,k
  real :: R_log

  ! Logical 
  logical :: oktest,oktest_me

  ! B0, the state and the sum of weights at the position of the satellite
  real :: StateSat_V(0:nVar+3), B0Sat_D(3)

  !-------------------------------------------------------------------------
  call set_oktest('set_logvar',oktest,oktest_me)

  if(oktest_me.and.n_step==1)then
     write(*,*)'nLogVar,nLogR,nLogTot:',nLogVar,nLogR,nLogTot
     write(*,*)'NameLogVar_I:',NameLogVar_I(1:nLogVar)
     write(*,*)'LogR_I:',LogR_I(1:nLogR)
  end if

  tmp1_BLK=1.00
  volume  =integrate_BLK(nProc,tmp1_BLK)

  ! Obtain data to calculate log variables
  if(iSat>=1)then
     ! Satellites need B0 and the state at the satellite position
     call get_b0(xSatellite(iSat,1),xSatellite(iSat,2),xSatellite(iSat,3),&
          B0Sat_D)
     call get_point_data(0.0,xSatellite(iSat,:),1,nBlock,1,nVar+3,StateSat_V)
     call collect_satellite_data(xSatellite(iSat,:),StateSat_V)
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

    ! Local variables
    real :: Bx, By, Bz, RhoUx, RhoUy, RhoUz, bDotB, bDotU

    ! External function for ionosphere    !^CFG IF IONOSPHERE
    real, external :: logvar_ionosphere   !^CFG IF IONOSPHERE

    select case(NameLogVar)

!!$! MHD variables averaged over the computational domain
    case('rho')
       LogVar_I(iVarTot) = StateIntegral_V(rho_)/volume
    case('rhoux')
       LogVar_I(iVarTot) = StateIntegral_V(rhoUx_)/volume
    case('rhouy')
       LogVar_I(iVarTot) = StateIntegral_V(rhoUy_)/volume
    case('rhouz')
       LogVar_I(iVarTot) = StateIntegral_V(rhoUz_)/volume
    case('bx')
       LogVar_I(iVarTot) = StateIntegral_V(Bx_)/volume
    case('by')
       LogVar_I(iVarTot) = StateIntegral_V(By_)/volume
    case('bz')
       LogVar_I(iVarTot) = StateIntegral_V(Bz_)/volume
    case('p')
       LogVar_I(iVarTot) = StateIntegral_V(P_)/volume
    case('e')
       LogVar_I(iVarTot) = integrate_BLK(1,E_BLK)/volume
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
               State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/volume
    case('uy')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/volume
    case('uz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK) / &
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = integrate_BLK(1,tmp1_BLK)/volume
    case('ekinx')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/volume
    case('ekiny')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/volume
    case('ekinz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2/&
               State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/volume
    case('ekin')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               (State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2+&
               State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2+&
               State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2)&
               /State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)
       end do
       LogVar_I(iVarTot) = cHalf*integrate_BLK(1,tmp1_BLK)/volume

    case('dst')
       ! Calculate the Biot-Savart formula for the center of the Earth:
       ! B = (1/4 Pi)*integral( (curl B) x R / R**3 dV)
       ! where R is the vector FROM the CURRENT to the FIELD POINT!
       ! Since the field point is at the origin, R = (-x,-y,-z)
       ! Only the Z component is calculated here: (J x R)_z = -J_x*y + Jy*x

       ! Calculate 
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle           
          do k=1,nK; do j=1,nJ; do i=1,nI
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
             tmp1_BLK(i,j,k,iBLK) = (                             &
                  ((State_VGB(Bz_,i,j+1,k,iBLK)                   &
                  - State_VGB(Bz_,i,j-1,k,iBLK)) / dy_BLK(iBLK) - &
                  ( State_VGB(By_,i,j,k+1,iBLK)                   &
                  - State_VGB(By_,i,j,k-1,iBLK)) / dz_BLK(iBLK)   &
                  ) * y_BLK(i,j,k,iBLK)                           &
                  -                                               &
                  ((State_VGB(Bx_,i,j,k+1,iBLK)                   &
                  - State_VGB(Bx_,i,j,k-1,iBLK)) / dz_BLK(iBLK) - &
                  ( State_VGB(Bz_,i+1,j,k,iBLK)                   &
                  - State_VGB(Bz_,i-1,j,k,iBLK)) / dx_BLK(iBLK)   &
                  ) * x_BLK(i,j,k,iBLK)                           &
                  ) / r_BLK(i,j,k,iBLK)**3
          end do; end do; end do
       end do
       ! The negative sign is due to R = (-x,-y,-z), 
       ! the 0.5 is because the centered difference formulae above used 
       ! (1/Dx) instead of the correct 1/(2Dx), 
       ! the 4pi is part of the Biot-Savart formula
       LogVar_I(iVarTot) = -0.5 * integrate_BLK(1,tmp1_BLK) / (4*cPi)

    case('dstdivb')
       ! Calculate the contribution of Div B to the surface integral of DST
       ! Error = (1/4 Pi)*integral( (div B) R / R**3 dV)
       ! Since the field point is at the origin, R = (-x,-y,-z)
       ! Only the Z component is calculated here: z * div B/R**3

       ! Calculate 
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle           
          do k=1,nK; do j=1,nJ; do i=1,nI
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
            LogVar_I(iVarTot) = State_VGB(rho_,iTest,jTest,kTest,BlkTest)
    case('rhouxpnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(rhoUx_,iTest,jTest,kTest,BlkTest)
    case('rhouypnt')     
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(rhoUy_,iTest,jTest,kTest,BlkTest)
    case('rhouzpnt')      
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(rhoUz_,iTest,jTest,kTest,BlkTest)
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
       tmp1_BLK = B0xcell_BLK + State_VGB(Bx_,:,:,:,:) 
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            B0xcell_BLK(iTest,jTest,kTest,BlkTest) + &
            State_VGB(Bx_,iTest,jTest,kTest,BlkTest)
    case('bypnt')                      
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            B0ycell_BLK(iTest,jTest,kTest,BlkTest) + &
            State_VGB(By_,iTest,jTest,kTest,BlkTest)
    case('bzpnt')                      
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            B0zcell_BLK(iTest,jTest,kTest,BlkTest) + &
            State_VGB(Bz_,iTest,jTest,kTest,BlkTest)
    case('ppnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = State_VGB(P_,iTest,jTest,kTest,BlkTest)
    case('epnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = E_BLK(iTest,jTest,kTest,BlkTest)
    case('uxpnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(rhoUx_,iTest,jTest,kTest,BlkTest) / &
            State_VGB(rho_,  iTest,jTest,kTest,BlkTest)
    case('uypnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(rhoUy_,iTest,jTest,kTest,BlkTest) / &
            State_VGB(rho_,  iTest,jTest,kTest,BlkTest)
    case('uzpnt')
       if(iProc == ProcTest) LogVar_I(iVarTot) = &
            State_VGB(rhoUz_,iTest,jTest,kTest,BlkTest) / &
            State_VGB(rho_,  iTest,jTest,kTest,BlkTest)
!!$!RAYTRACE variables                                ^CFG  IF RAYTRACE BEGIN
       ! RAYTRACE variables at Itest, Jtest, Ktest, BLKtest, PROCtest
    case('theta1pnt')
       if(iProc == ProcTest) &
            LogVar_I(iVarTot) = ray(1,1,iTest,jTest,kTest,BlkTest)
    case('theta2pnt')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(1,2,1:nI,1:nJ,1:nK,:)
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
          R_log = LogR_I(iR)
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,1:nBlock) = 1.0
          LogVar_I(iVarTot) = integrate_sphere(360, R_log, tmp1_BLK)
       end do
    case('rhoflx')
       ! rho*U_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          R_log = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK)) CYCLE
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                tmp1_BLK(i,j,k,iBLK) = &
               (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
               +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
               +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
               )/R_BLK(i,j,k,iBLK)
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = integrate_sphere(360, R_log, tmp1_BLK)
       end do
    case('dstflx')
       ! B1z
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          R_log = LogR_I(iR)
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,1:nBlock) = &
               State_VGB(Bz_,0:nI+1,0:nJ+1,0:nK+1,1:nBlock)
          LogVar_I(iVarTot) = integrate_sphere(180, R_log, tmp1_BLK) / &
               (4*cPi*R_log**2)
       end do
    case('bflx')
       ! B_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          R_log = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK)) CYCLE
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                tmp1_BLK(i,j,k,iBLK) = ( &
                     (State_VGB(Bx_,i,j,k,iBLK) &
                     +B0xCell_BLK(i,j,k,iBLK))*x_BLK(i,j,k,iBLK) + &
                     (State_VGB(By_,i,j,k,iBLK) &
                     +B0yCell_BLK(i,j,k,iBLK))*y_BLK(i,j,k,iBLK) + &
                     (State_VGB(Bz_,i,j,k,iBLK) &
                     +B0zCell_BLK(i,j,k,iBLK))*z_BLK(i,j,k,iBLK) &
                     ) / R_BLK(i,j,k,iBLK)
             end do; end do; end do
          end do
 
          LogVar_I(iVarTot) = integrate_sphere(360, R_log, tmp1_BLK)
       end do
    case('b2flx')
       ! B^2*u_R
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          R_log = LogR_I(iR)           
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))cycle           
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = &
                  (State_VGB(Bx_,i,j,k,iBLK) + B0xCell_BLK(i,j,k,iBLK))**2 + &
                  (State_VGB(By_,i,j,k,iBLK) + B0yCell_BLK(i,j,k,iBLK))**2 + &
                  (State_VGB(Bz_,i,j,k,iBLK) + B0zCell_BLK(i,j,k,iBLK))**2 
             tmp1_BLK(i,j,k,iBLK) = 0.5*tmp1_BLK(i,j,k,iBLK)* &
                  ( State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  + State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  + State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  ) / (State_VGB(rho_,i,j,k,iBLK)*R_BLK(i,j,k,iBLK))
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = integrate_sphere(360, R_log,tmp1_BLK)
       end do
    case('pvecflx')
       iVarTot = iVarTot - 1
       do iR=1,nLogR
          iVarTot = iVarTot + 1
          R_log = LogR_I(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))CYCLE
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                Bx = State_VGB(Bx_,i,j,k,iBLk)+B0xCell_BLK(i,j,k,iBLk)
                By = State_VGB(By_,i,j,k,iBLk)+B0yCell_BLK(i,j,k,iBLk)
                Bz = State_VGB(Bz_,i,j,k,iBLk)+B0zCell_BLK(i,j,k,iBLk)
                RhoUx = State_VGB(rhoUx_,i,j,k,iBLk)
                RhoUy = State_VGB(rhoUy_,i,j,k,iBLk)
                RhoUz = State_VGB(rhoUz_,i,j,k,iBLk)
                bDotb = Bx**2 + By**2 + Bz**2
                bDotU = Bx*RhoUx + By*RhoUy + Bz*RhoUz
                tmp1_BLK(i,j,k,iBLk) = ( &
                  ( bDotb*rhoUx - bDotU*Bx )*X_BLK(i,j,k,iBLk) + &
                  ( bDotb*rhoUy - bDotU*By )*Y_BLK(i,j,k,iBLk) + &  
                  ( bDotb*rhoUz - bDotU*Bz )*Z_BLK(i,j,k,iBLk)   &  
                  ) / (State_VGB(rho_,i,j,k,iBLk)*R_BLK(i,j,k,iBLk))
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = integrate_sphere(360, R_log,tmp1_BLK)
       end do


       ! simple circular integrals
    case('e2dflx')
       !this is the azimuthal component of the electric field 
       !integrated around a circle
       ! Ex*y - Ey*x
       iVarTot = iVarTot-1
       do iR=1,nLogR
          iVarTot = iVarTot+1
          R_log = LogR_I(iR)
          do iBLK = 1,nBlock
             if(unusedBLK(iBLK))CYCLE
             do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
                Bx = State_VGB(Bx_,i,j,k,iBLk)+B0xCell_BLK(i,j,k,iBLk)
                By = State_VGB(By_,i,j,k,iBLk)+B0yCell_BLK(i,j,k,iBLk)
                Bz = State_VGB(Bz_,i,j,k,iBLk)+B0zCell_BLK(i,j,k,iBLk)
                RhoUx = State_VGB(rhoUx_,i,j,k,iBLk)
                RhoUy = State_VGB(rhoUy_,i,j,k,iBLk)
                RhoUz = State_VGB(rhoUz_,i,j,k,iBLk)
                tmp1_BLK(i,j,k,iBLK) = &
                     ( (RhoUy*Bz-RhoUz*By)*y_BLK(i,j,k,iBLK) &
                     - (RhoUz*Bx-RhoUx*Bz)*x_BLK(i,j,k,iBLK) &
                     ) / (State_VGB(rho_,i,j,k,iBLK)*R_BLK(i,j,k,iBLK)) 
             end do; end do; end do
          end do
          LogVar_I(iVarTot) = integrate_flux_circ(nProc,R_log,0.0,tmp1_BLK)
       end do

       ! OTHER VALUES
    case('dt')
       if(iProc == 0)LogVar_I(iVarTot) = dt

       ! DEFAULT
    case default
       if (UseUserLogFiles) then                 !^CFG IF USERFILES BEGIN
          call user_get_log_var( LogVar_I(iVarTot), NameLogVar )
       else                                      !^CFG END USERFILES
          if(iProc == 0)then
             LogVar_I(iVarTot) = -7777.
             call write_myname;
             write(*,*) 'WARNING in set_logvar: unknown logvarname=',NameLogVar
          end if
       end if                                    !^CFG IF USERFILES
    end select
  end subroutine set_log_var

  !==========================================================================
  subroutine set_sat_var
    use ModProcMH
    use ModNumConst
    use ModVarIndexes
    use ModIO
    implicit none
    !-------------------------------------------------------------------------
    if (iProc/=0) RETURN

    select case(NameLogVar)
    case('rho')
       LogVar_I(iVarTot) = StateSat_V(Rho_)
    case('rhoux')
       LogVar_I(iVarTot) = StateSat_V(RhoUx_)
    case('rhouy')
       LogVar_I(iVarTot) = StateSat_V(RhoUy_)
    case('rhouz')
       LogVar_I(iVarTot) = StateSat_V(RhoUz_)
    case('bx')
       LogVar_I(iVarTot) = StateSat_V(Bx_)+B0Sat_D(1)
    case('by')
       LogVar_I(iVarTot) = StateSat_V(By_)+B0Sat_D(2)
    case('bz')
       LogVar_I(iVarTot) = StateSat_V(Bz_)+B0Sat_D(3)
    case('e')
       LogVar_I(iVarTot) = 0.5*(&
            sum(StateSat_V(RhoUx_:RhoUz_)**2)/StateSat_V(rho_) + &
            sum((StateSat_V(Bx_:Bz_)+B0Sat_D)**2))
    case('p')
       LogVar_I(iVarTot) = StateSat_V(P_)
    case('ux')
       LogVar_I(iVarTot) = StateSat_V(RhoUx_)/StateSat_V(Rho_)
    case('uy')
       LogVar_I(iVarTot) = StateSat_V(RhoUy_)/StateSat_V(Rho_)
    case('uz')
       LogVar_I(iVarTot) = StateSat_V(RhoUz_)/StateSat_V(Rho_)
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
    case default
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
  implicit none

  integer, intent(in) :: nLogVar, nLogR, nLogTot
  character (LEN=10), intent(in) :: NameLogVar_I(nLogVar)
  real, intent(out) :: LogVar_I(nLogTot)
  real, intent(in) :: LogR_I(nLogR)
  character (len=10) :: NameLogVar

  integer :: iVar, iR, iVarTot
  !-------------------------------------------------------------------------

  iVarTot = 0
  do iVar=1,nLogVar

     iVarTot = iVarTot+1
     NameLogVar = NameLogVar_I(iVar)
     call normalize_name_log_var(NameLogVar_I(iVar), NameLogVar)

     select case(NameLogVar)

!!$! BASIC MHD variables
     case('rho','rhopnt') 
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_n
     case('rhoux','rhouy','rhouz', 'rhouxpnt','rhouypnt','rhouzpnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_rhoU
     case('bx','by','bz','bxpnt','bypnt','bzpnt','b1xpnt','b1ypnt','b1zpnt', &
          'dst','dstdivb')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_B
     case('e','epnt')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*unitUSER_energydens
     case('p','ppnt','pmin','pmax')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*unitUSER_p
     case('ux','uy','uz','uxpnt','uypnt','uzpnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_U
     case('ekinx','ekiny','ekinz','ekin')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_rho*unitUSER_U**2
     case('jx','jy','jz','jxpnt','jypnt','jzpnt')
        LogVar_I(iVarTot)= LogVar_I(iVarTot)*unitUSER_J

!!$! Ionosphere values                                !^CFG IF IONOSPHERE BEGIN
     case('cpcpn','cpcp_n','cpcp_north','cpcpnorth',&
          'cpcps','cpcp_s','cpcp_south','cpcpsouth')
        ! User unit is kV = 1000 V
        LogVar_I(iVarTot) = LogVar_I(iVarTot)/1000.0  !^CFG END IONOSPHERE

!!$! Flux values
     case('aflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitUSER_x**2)
        iVarTot = iVarTot+nLogR-1
     case('rhoflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitSI_n*unitSI_U*unitSI_x**2)
        iVarTot = iVarTot+nLogR-1
     case('dstflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *unitUSER_B
        iVarTot = iVarTot+nLogR-1
     case('bflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitSI_B*unitSI_x**2)
        iVarTot = iVarTot+nLogR-1
     case('b2flx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitSI_Poynting*unitSI_x**2)
        iVarTot = iVarTot+nLogR-1
     case('pvecflx')
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitSI_Poynting*unitSI_x**2)
        iVarTot = iVarTot+nLogR-1
     case('e2dflx') 
        ! circular integral if
        ! the azimuthal component of the electric field
        ! on the surface of a sphere
        LogVar_I(iVarTot:iVarTot+nLogR-1) = LogVar_I(iVarTot:iVarTot+nLogR-1)&
             *(unitSI_electric*unitSI_x)
        iVarTot = iVarTot+nLogR-1

        ! OTHER VALUES
     case('dt')
        LogVar_I(iVarTot) = LogVar_I(iVarTot)*unitUSER_t

     case default
        ! no normalization

     end select
  end do ! iVar
end subroutine normalize_logvar

!==============================================================================

real function integrate_sphere(nTheta,Radius,Array)

  ! This function calculates the integral of the incomming variable Array
  ! over the surface of a sphere centered at the origin radius Radius.  
  ! The resolution in the colatitude is determined by the nTheta parameter.

  use ModMain,           ONLY: nI,nJ,nK,nBLK,nBlock,unusedBLK
  use ModGeometry,       ONLY: x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModNumConst
  use ModCoordTransform, ONLY: sph_to_xyz
  implicit none

  ! Arguments

  integer, intent(in) :: nTheta
  real, intent(in)    :: Array(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK)
  real, intent(in)    :: Radius 

  ! Local variables
  real :: Integral, Darea0, Darea, Average 

  ! Indices and coordinates
  integer :: iBlock,i,j,i1,i2,j1,j2,k1,k2
  integer :: MaxPhi, nPhi
  real    :: x, y, z, DxInv, DyInv, DzInv, xNorm, yNorm, zNorm, Dx, Dy, Dz
  real    :: xMin, xMax, yMin, yMax, zMin, zMax
  real    :: dTheta, dPhi, Theta, Phi, SinTheta, CosTheta

  ! Store cartesian coordinates for sake of efficiency
  ! The x and y depend on iPhi,iTheta while z only depends on iTheta
  !  real, allocatable :: x_II(:,:), y_II(:,:), z_I(:), SinTheta_I(:)

  logical :: DoTest,DoTestMe
  !---------------------------------------------------------------------------

  call set_oktest('integrate_sphere',DoTest,DoTestMe)
  call timing_start('int_sphere')

  ! Get the angular resolution from the input parameter nTheta
  MaxPhi = 2*nTheta
  dTheta = cPi/nTheta
  dArea0 = Radius**2 * cTwo * sin(cHalf*dTheta)

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
  Integral = 0.0
  do iBlock = 1, nBlock

     if (unusedBLK(iBlock)) CYCLE

     ! get the max and min radial distance for this block so that we can check
     ! whether or not this block contibutes to the sum.

     xMin = cHalf*(x_BLK( 0, 0, 0,iBlock) + x_BLK(   1,   1  , 1,iBlock))
     xMax = cHalf*(x_BLK(nI,nJ,nK,iBlock) + x_BLK(nI+1,nJ+1,nK+1,iBlock))
     yMin = cHalf*(y_BLK( 0, 0, 0,iBlock) + y_BLK(   1,   1,   1,iBlock))
     yMax = cHalf*(y_BLK(nI,nJ,nK,iBlock) + y_BLK(nI+1,nJ+1,nK+1,iBlock))
     zMin = cHalf*(z_BLK( 0, 0, 0,iBlock) + z_BLK(   1,   1,   1,iBlock))
     zMax = cHalf*(z_BLK(nI,nJ,nK,iBlock) + z_BLK(nI+1,nJ+1,nK+1,iBlock))

     if( minmod(xMin,xMax)**2 + minmod(yMin,yMax)**2 + minmod(zMin,zMax)**2 &
          > Radius**2) &
          CYCLE
     if ( max(abs(xMin),abs(xMax))**2 + max(abs(yMin),abs(yMax))**2 + &
          max(abs(zMin),abs(zMax))**2 < Radius**2 ) &
          CYCLE

     DxInv = cOne/dx_BLK(iBlock)
     DyInv = cOne/dy_BLK(iBlock)
     DzInv = cOne/dz_BLK(iBlock)

     do i = 1, nTheta
        ! Check if z is inside the block
        Theta = dTheta*(i-cHalf)
        z     = Radius*cos(Theta)
        if(z <  zMin) CYCLE
        if(z >= zMax) CYCLE

        SinTheta = sin(Theta)

        ! Number of Phi coordinates is proportional to 2*nTheta*SinTheta
        nPhi = min(MaxPhi, 4 * ceiling(cQuarter*MaxPhi*SinTheta))
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

           ! Convert to normalized coordinates 
           ! (in cell centers the nomalized coordinate equals the index)
           xNorm = DxInv*(x - x_BLK(1,1,1,iBlock)) + cOne
           yNorm = DyInv*(y - y_BLK(1,1,1,iBlock)) + cOne
           zNorm = DzInv*(z - z_BLK(1,1,1,iBlock)) + cOne

           ! Determine cell indices corresponding to location 
           i1 = nint(xNorm)
           j1 = nint(yNorm)
           k1 = nint(zNorm)

           ! Distance relative to the cell centers
           Dx = xNorm-i1
           Dy = yNorm-j1
           Dz = zNorm-k1

           ! The indexes of the cells in the direction of the point
           i2 = i1 + sign(1.0, Dx)
           j2 = j1 + sign(1.0, Dy)
           k2 = k1 + sign(1.0, Dz)

           ! In the interpolation formula only the abs(distance) occurs
           Dx = abs(Dx)
           Dy = abs(Dy)
           Dz = abs(Dz)

           ! Second order interpolation in 3 directions
           Average = &
                Array(i1,j1,k1,iBlock)*(1-Dx-Dy-Dz) + &
                Array(i2,j1,k1,iBlock)*Dx             + &
                Array(i1,j2,k1,iBlock)*Dy             + &
                Array(i1,j1,k2,iBlock)*Dz

           Integral = Integral + dArea * Average

        end do
     end do
  end do
  ! deallocate(x_II, y_II, z_I, SinTheta_I)

  integrate_sphere = Integral
  call timing_stop('int_sphere')

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

end function integrate_sphere


!==============================================================================

real function integrate_flux_circ(qnum,qrad,qz,qa)

  ! This function calculates the integral of the incomming variable qa
  ! over a cirle parallel to the equitorial plane.  The radius of the circle
  ! for the z axis is defined by the radius qrad and the z position is
  ! is given by qz.  

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlock,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModNumConst
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: qnum
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
       intent(in) :: qa
  real, intent(in) :: qrad, qz 

  ! Local variables
  integer :: iError
  real :: integral_sum, GLOBAL_integral_sum, term_to_add 

  ! Indices and coordinates
  integer :: iBLK,i,j,i1,i2,j1,j2,k1,k2
  integer :: nphi_circ
  real :: xx1,xx2,yy1,yy2,zz1,zz2,minRblk, maxRblk
  real :: x,y,z,xx,yy,zz,dx1,dx2,dy1,dy2,dz1,dz2
  real :: dphi_circ,phi_circ

  logical :: oktest,oktest_me
  !---------------------------------------------------------------------------

  integral_sum = 0.0
  GLOBAL_integral_sum = 0.0

  call set_oktest('integrate_flux_circ',oktest,oktest_me)

  ! the angular resolution of the integral is hard coded
  ! in degrees
  dphi_circ = 0.5

  ! now comcircpute the number of points that correspond to the above
  ! intervals and then recompute the intervals in case an even multiple
  ! was not chosen.
  nphi_circ   = 360.0/dphi_circ
  dphi_circ   = 360.0/real(nphi_circ)
  !convert to radians
  dphi_circ   = dphi_circ*cPi/180.0

  if (oktest_me) write(*,*) 'nphi,dphi',nphi_circ,dphi_circ

  ! Sum all cells within range

  do iBLK = 1,nBlock

     if (unusedBLK(iBLK)) CYCLE
     ! get the max and min radial (cylindrical) distance for this block so 
     ! that we can check whether or not this block contibutes to the sum.
     xx1 = 0.50*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
     xx2 = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
     yy1 = 0.50*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
     yy2 = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
     zz1 = 0.50*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
     zz2 = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))
     minRblk = sqrt((min(abs(xx1),abs(xx2)))**2 + &
          (min(abs(yy1),abs(yy2)))**2)
     maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
          (max(abs(yy1),abs(yy2)))**2)
     if (minRblk > qrad .or. maxRblk < qrad ) CYCLE


     do i=1,nphi_circ
        phi_circ = (i-.5)*dphi_circ

        ! get the xyz coordinates
        x = qrad*cos(phi_circ)
        y = qrad*sin(phi_circ)
        z = qz

        ! check to see if this point is inside the block - if so add it to the sum
        if (x >= xx1 .and. x < xx2 .and. &
             y >= yy1 .and. y < yy2 .and. &
             z >= zz1 .and. z < zz2 ) then

           !compute the interpolated values at the current location.
           ! Convert to normalized coordinates (index and position are the same)
           xx=(x-x_BLK(1,1,1,iBLK))/dx_BLK(iBLK)+1.
           yy=(y-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1.
           zz=(z-z_BLK(1,1,1,iBLK))/dz_BLK(iBLK)+1.

           ! Determine cell indices corresponding to location 
           i1=floor(xx); i2=i1+1
           j1=floor(yy); j2=j1+1
           k1=floor(zz); k2=k1+1

           ! Distance relative to the cell centers
           dx1=xx-i1; dx2=1.-dx1
           dy1=yy-j1; dy2=1.-dy1
           dz1=zz-k1; dz2=1.-dz1

           ! Bilinear interpolation in 3D
           term_to_add = &
                dx1*(   dy1*(   dz1*qa(i2,j2,k2,iBLK)+&
                dz2*qa(i2,j2,k1,iBLK))+&
                dy2*(   dz1*qa(i2,j1,k2,iBLK)+&
                dz2*qa(i2,j1,k1,iBLK)))+&
                dx2*(   dy1*(   dz1*qa(i1,j2,k2,iBLK)+&
                dz2*qa(i1,j2,k1,iBLK))+&
                dy2*(   dz1*qa(i1,j1,k2,iBLK)+&
                dz2*qa(i1,j1,k1,iBLK)))

           integral_sum = integral_sum + term_to_add

        end if

     end do

  end do

  ! now multiply by the size of each interval in the integral.  Since they are all the
  ! same we can do this after the fact and not inside the above loops

  integral_sum = integral_sum*qrad*dphi_circ

  if(qnum>1)then
     call MPI_allreduce(integral_sum, GLOBAL_integral_sum, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
  else
     GLOBAL_integral_sum = integral_sum
  end if
  integrate_flux_circ = GLOBAL_integral_sum

end function integrate_flux_circ

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
  implicit none

  character(len=*), intent(in)  :: NameIn
  character(len=*), intent(out) :: NameOut
  integer :: l
  !---------------------------------------------------------------------------
  NameOut = NameIn

  ! Switch to all lower case
  call lower_case(NameOut)

  ! Replace mx with rhoux, my with rhouy, mz with rhouz
  select case(NameOut(1:2))
  case('mx','my','mz')
     NameOut = 'rhou'//NameOut(2:len(NameOut)-3)
  end select

  ! Replace *test with *pnt
  l = len_trim(NameOut)
  if(NameOut(l-3:l) == 'test') NameOut(l-3:l)='pnt '

end subroutine normalize_name_log_var

!^CFG IF CARTESIAN BEGIN
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
  use ModNumConst
  implicit none 

  ! Arguments
  real, intent(out) :: Sum_V(nVar)
  real, intent(out) :: Pressure_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)

  ! Local variables:
  integer :: iBlock, iVar

  logical :: DoTest, DoTestMe

  !---------------------------------------------------------------------------

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
!^CFG END CARTESIAN
