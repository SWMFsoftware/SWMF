!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_logfile(log_type,ifile)
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY  : State_VGB
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModPhysics, ONLY  : unitUSER_t,unitSI_t
  use ModIO
  use ModIoUnit, ONLY   : io_unit_new
  use ModUtilities, ONLY: flush_unit
  implicit none

  ! Arguments

  integer, intent(in) :: ifile

  ! log_type = 0 -> write logfile 
  ! log_type >=1 -> write satellite output files (number = isat)
  integer, intent(in) :: log_type

  ! Logfile variables
  integer, parameter :: nlogvarmax=40
  integer, parameter :: nlogRmax = 10
  real :: LogVar(nlogvarmax)
  character (len=10) :: logvarnames(nlogvarmax)
  character (len=10) :: logRnames(nlogRmax)
  real :: logRvalues(nlogRmax)
  integer :: nlogvar, nlogR, nlogTot, i,j,isat
  integer :: nfluxvar       ! number of flux variables used
  integer :: unit_outfile   ! local unit number
  real :: xLOC(3)
  character (LEN=255) :: allnames
  character (LEN=100) :: time_type

  ! Logical 
  logical :: fileopen, write_position
  logical :: oktest,oktest_me
  real :: pmin, pmax
  real, external :: maxval_loc_BLK, minval_loc_BLK
  integer :: loc(5)
  integer :: iTime_I(7) ! integer time: year,month,day,hour,minute,sec,msec

  !---------------------------------------------------------------------------
  call set_oktest('write_logfile',oktest,oktest_me)

  write_position = .false.

  if(log_type == 0) then
     if(index(test_string,'show_pmin')>0)then
        pmin = minval_loc_BLK(nProc,State_VGB(P_,:,:,:,:),loc)
        if(loc(5)==iProc)write(*,*)'pmin, loc, x, y, z=',pmin,loc, &
             x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     end if

     if(index(test_string,'show_pmax')>0)then
        pmax = maxval_loc_BLK(nProc,State_VGB(P_,:,:,:,:),loc)
        if(loc(5)==iProc)write(*,*)'pmax, loc, x, y, z=',pmax,loc, &
             x_BLK(loc(1),loc(2),loc(3),loc(4)),&
             y_BLK(loc(1),loc(2),loc(3),loc(4)),&
             z_BLK(loc(1),loc(2),loc(3),loc(4))
     end if
  end if

  if (log_type ==0) then
     call split_str(log_vars,nlogvarmax,logvarnames,nlogvar)
     nfluxvar = 0
     do i=1,nlogvar
        if(index(logvarnames(i),'flx')>0) nfluxvar = nfluxvar + 1
        if(index(logvarnames(i),'pnt')>0 .or.  &
             index(logvarnames(i),'test')>0)  write_position = .true.
     end do
     if (nfluxvar >0) then
        call split_str(log_R_str, nlogRmax,logRnames,nlogR)
        do i=1,nlogR
           read(logRnames(i),*) logRvalues(i)
        end do
     else
        nlogR = 0
        logRvalues = 0
     end if
     time_type = log_time
  elseif (log_type>=1) then
     isat = log_type
     if (.not. DoTrackSatellite_I(isat)) return
     call split_str(satellite_vars(satellite_+isat),nlogvarmax, &
          logvarnames,nlogvar)
     nlogR = 0
     logRvalues = 0
     nfluxvar = 0
     time_type = sat_time(satellite_+isat)
     write_position = .true.
  end if


  ! check to make sure that the total number of Logfile variables is smaller
  ! than nlogvarmax.  If it is not write a warning message and truncate the
  ! list.
  if (nlogvar + nfluxvar*nlogR > nlogvarmax) then
     write(*,*)'Warning in write_logfile: Number of logfile variables exceeds '
     write(*,*)'the array dimensions.  Truncating list - recompile with larger'
     write(*,*)'array dimensions'
     if (nlogvar >= nlogvarmax) then
        nlogvar = nlogvarmax
        nlogR = 1
     else
        nlogR = (nlogvarmax-nlogvar)/nfluxvar + 1
     end if
  end if
  nlogTot = nlogvar + nfluxvar*(nlogR-1)

  ! load the time variables to the list of output
  if(index(time_type,'none')>0) then
     allnames = ''
  else
     allnames = ''
     if(index(time_type,'step')>0) allnames = 'it'
     if(index(time_type,'date')>0) &
          write(allnames,'(a)') allnames(1:len_trim(allnames))// &
          ' year mo dy hr mn sc msc'
     if(index(time_type,'time')>0) &
          write(allnames,'(a)') allnames(1:len_trim(allnames))//' t'
  end if

  if (write_position) then
     write(allnames,'(a)') allnames(1:len_trim(allnames))//' X Y Z'
     if (log_type==0) then
        xLOC(1)=Xtest_mod; xLOC(2)=Ytest_mod; xLOC(3)=Ztest_mod; 
     elseif (log_type>=1) then
        xLOC(1)=XSatellite(isat,1)
        xLOC(2)=XSatellite(isat,2)
        xLOC(3)=Xsatellite(isat,3) 
     end if
  end if

  do i=1,nlogvar
     if(index(logvarnames(i),'flx')>0) then
        do j=1,nlogR
           allnames = trim(allnames)//' '//trim(logvarnames(i))//'_R='// &
                trim(logRnames(j))
        end do
     else
        allnames = trim(allnames)//' '//trim(logvarnames(i))
     end if
  end do

  if(oktest_me.and.n_step==1)then
     write(*,*)'nlogvar,nfluxvar,nlogR,nlogTot:',  &
          nlogvar,nfluxvar,nlogR,nlogTot
     write(*,*)'logvarnames:',logvarnames(1:nlogvar)
     write(*,*)'logRnames:',logRnames(1:nlogR)
     write(*,*)'logRvalues:',logRvalues(1:nlogR)
     write(*,*)'allnames:',trim(allnames)
  end if

  if(iProc==0) then
     if (log_type==0) then
        if(unit_log<0)then
           unit_log = io_unit_new()
           write(filename,'(a,i6.6,a)')&
                trim(NamePlotDir)//'log_n',n_step,'.log'
           open(unit_log,file=filename,status="replace")
           write(unit_log,'(a)')'Volume Averages, Fluxes, or point values'
           if (index(allnames,'pnt')>0 .or. index(allnames,'test')>0) then
   	      if (coord_test) then
   		 write(unit_log,'(a,3(1pe13.5),a)')  &
                      'Requested Test point (X,Y,Z): ', &
   		      Xtest,Ytest,Ztest,'--Using nearest cell center.'
   	      else
   		 write(unit_log,'(a,5(i6))')  &
                      'Requested Test point(I,J,K,BLK,PROC): ', &
   		      Itest,Jtest,Ktest,BLKtest,PROCtest
   	      end if
           end if
           write(unit_log,'(a)')allnames(1:len_trim(allnames))
        endif
        unit_outfile = unit_log
     elseif (log_type >= 1) then
        unit_outfile = iUnitSat_I(iSat)
        if (Satellite_first_write(isat)) then
           write(unit_outfile,'(a)')  &
                'Satellite data for Satellite: '//&
                Satellite_name(isat)(1:len_trim(Satellite_name(isat)))
           write(unit_outfile,'(a)')allnames(1:len_trim(allnames))
           Satellite_first_write(isat)=.false.
        end if
     end if
  endif

  call set_logvar(nlogvar,logvarnames,nlogR,logRvalues,nlogTot,LogVar,log_type)

  ! WRITE OUT THE LINE INTO THE LOGFILE OR THE SATELLITE FILE
  if(iProc==0) then

     if (plot_dimensional(ifile))  &
          call normalize_logvar(nlogvar,logvarnames,nlogR,logRvalues, &
          nlogTot,LogVar)

     ! first output the appropriate time data
     if(index(time_type,'none')>0) then
        ! do nothing
     else
        if(index(time_type,'step')>0) &
             write(unit_outfile,'(i7)',ADVANCE='NO') n_step
        if(index(time_type,'date')>0) then
           call get_date_time(iTime_I)
           write(unit_outfile,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') &
                iTime_I
        end if
        if(index(time_type,'time')>0) then
           if(plot_dimensional(ifile)) then
              !note that the funny (unitUSER_t/unitSI_t) is because by default
              !Time_Simulation is in seconds(SI).  If the user wants output
              !in some other time units this is the conversion.  In most 
              !cases the two are the same the and (...)=1.
              write(unit_outfile,'(1pe13.5)',ADVANCE='NO') &
                   Time_Simulation*(unitUSER_t/unitSI_t)
           else
              write(unit_outfile,'(1pe13.5)',ADVANCE='NO') &
                   Time_Simulation/unitSI_t
           end if
        end if
     end if

     ! Now write the position of the test point or satellite if desired
     if(write_position) &
          write(unit_outfile,'(3(1pe13.5))',ADVANCE='NO') (xLOC(i),i=1,3)

     ! Finally write out the data variables
     write(unit_outfile,'(100(1pe13.5))') LogVar(1:nlogTot)

     call flush_unit(unit_outfile)
  end if

end subroutine write_logfile


!==============================================================================
subroutine set_logvar(nlogvar,logvarnames,nlogR,logRvalues,nlogTot,LogVar,isat)

  use ModProcMH
  use ModNumConst
  use ModMain, ONLY: n_step,dt,unusedBLK,nI,nJ,nK,nBlock,gcn,UseUserLogFiles
  use ModPhysics,    ONLY: rCurrents, rBody
  use ModVarIndexes
  use ModAdvance,    ONLY: tmp1_BLK,tmp2_BLK,B0xCell_BLK,B0yCell_BLK, &
       B0zCell_BLK, State_VGB, E_BLK, DivB1_GB
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
       x1,x2,y1,y2,z1,z2
  use ModRaytrace,   ONLY: ray  !^CFG  IF RAYTRACE
  use ModIO
  implicit none

  integer, intent(in) :: Nlogvar, nlogR, nlogTot, isat
  character (LEN=10), intent(in) :: logvarnames(Nlogvar)
  real, intent(out) :: LogVar(nLogTot)
  real, intent(in) :: logRvalues(nlogR)
  character (len=10) :: NameLogVar

  real :: volume
  real,dimension(nVar) :: StateIntegral_V
  real, external :: integrate_BLK, maxval_BLK, minval_BLK
  real, external :: integrate_flux_sph, integrate_flux_circ,test_cell_value

  integer :: iVar,iR,iVarTot,itmp,jtmp, iBLK
  integer :: i,j,k
  real :: R_log

  real :: tmp2_G(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)

  ! Logical 
  logical :: oktest,oktest_me

  ! B0, the state and the sum of weights at the position of the satellite
  real :: StateSat_V(0:nVar+3), B0Sat_D(3)

  !-------------------------------------------------------------------------
  call set_oktest('set_logvar',oktest,oktest_me)

  if(oktest_me.and.n_step==1)then
     write(*,*)'nlogvar,nlogR,nlogTot:',nlogvar,nlogR,nlogTot
     write(*,*)'logvarnames:',logvarnames(1:nlogvar)
     write(*,*)'logRvalues:',logRvalues(1:nlogR)
  end if

  tmp1_BLK=1.00
  volume  =integrate_BLK(nProc,tmp1_BLK)

  ! Obtain data to calculate log variables
  if(iSat>=1)then
     ! Satellites need B0 and the state at the satellite position
     call get_b0(xSatellite(iSat,1),xSatellite(iSat,2),xSatellite(iSat,3),&
          B0Sat_D)
     call get_point_data(0.0,xSatellite(iSat,:),1,nVar+3,StateSat_V)
     call collect_satellite_data(xSatellite(iSat,:),StateSat_V)
  else
     ! The logfile usually needs the integral of conservative variables
     call integrate_cell_centered_vars(StateIntegral_V)
  end if

  iVarTot = 0
  do iVar=1,nLogVar

     iVarTot = iVarTot+1
     NameLogVar = logvarnames(iVar)

     ! If we are a satellite and not a logfile (isat>=1) then we should
     ! do only satellite variables so append a 'sat' to the end of the 
     ! variables
     if (isat >= 1) then
        call set_sat_var
     else
        call set_log_var
     end if
  end do

contains
  !============================================================================
  subroutine set_log_var

    ! External function for ionosphere    !^CFG IF IONOSPHERE
    real, external :: logvar_ionosphere   !^CFG IF IONOSPHERE

    select case(NameLogVar)

       ! BASIC MHD variables averaged over the computational domain
    case('rho')
       LogVar(iVarTot) = StateIntegral_V(rho_)/volume
    case('rhoUx','rhoux','mx')
       LogVar(iVarTot) = StateIntegral_V(rhoUx_)/volume
    case('rhoUy','rhouy','my')
       LogVar(iVarTot) = StateIntegral_V(rhoUy_)/volume
    case('rhoUz','rhouz','mz')
       LogVar(iVarTot) = StateIntegral_V(rhoUz_)/volume
    case('Bx','bx')
       LogVar(iVarTot) = StateIntegral_V(Bx_)/volume
    case('By','by')
       LogVar(iVarTot) = StateIntegral_V(By_)/volume
    case('Bz','bz')
       LogVar(iVarTot) = StateIntegral_V(Bz_)/volume
    case('E','e')
       LogVar(iVarTot) = integrate_BLK(nProc,E_BLK)/volume

       ! Maximum and minimum values of pressure 
       ! (tmp2_BLK=CellCenteredVar_VGB(P_,:,:,:,:))
    case('Pmin','pmin')
       LogVar(iVarTot) = minval_BLK(nProc,tmp2_BLK)
    case('Pmax','pmax')
       LogVar(iVarTot) = maxval_BLK(nProc,tmp2_BLK)

       ! Extra MHD variables averaged over the computational domain
    case('P','p')
       LogVar(iVarTot) = StateIntegral_V(P_)/volume
    case('Ux','ux')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUx_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = integrate_BLK(nProc,tmp1_BLK)/volume
    case('Uy','uy')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUy_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = integrate_BLK(nProc,tmp1_BLK)/volume
    case('Uz','uz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUz_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = integrate_BLK(nProc,tmp1_BLK)/volume
    case('Ekinx','ekinx')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUx_,:,:,:,iBLK)**2/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume
    case('Ekiny','ekiny')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUy_,:,:,:,iBLK)**2/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume
    case('Ekinz','ekinz')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUz_,:,:,:,iBLK)**2/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume
    case('Ekin','ekin')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               (State_VGB(rhoUx_,:,:,:,iBLK)**2+&
               State_VGB(rhoUy_,:,:,:,iBLK)**2+&
               State_VGB(rhoUz_,:,:,:,iBLK)**2)&
               /State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = cHalf*integrate_BLK(nProc,tmp1_BLK)/volume

    case('dst','Dst','DST')
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
       LogVar(iVarTot) = -0.5 * integrate_BLK(nProc,tmp1_BLK) / (4*cPi)

       ! Basic MHD variables at a single point given by Itest, Jtest, Ktest, 
       !   BLKtest, PROCtest

    case('dstdivb','DstDivb','DSTDIVB')
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
       ! The 4pi is part of the Biot-Savart formula
       LogVar(iVarTot) = integrate_BLK(nProc,tmp1_BLK) / (4*cPi)

       ! Basic MHD variables at a single point given by Itest, Jtest, Ktest, 
       !   BLKtest, PROCtest

    case('rhopnt','rhotest')
       LogVar(iVarTot) = test_cell_value(State_VGB(rho_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('rhoUxpnt','rhouxpnt','mxpnt','rhoUxtest','rhouxtest','mxtest')
       LogVar(iVarTot) = test_cell_value(State_VGB(rhoUx_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('rhoUypnt','rhouypnt','mypnt','rhoUytest','rhouytest','mytest')     
       LogVar(iVarTot) = test_cell_value(State_VGB(rhoUy_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('rhoUzpnt','rhouzpnt','mzpnt','rhoUztest','rhouztest','mztest')      
       LogVar(iVarTot) = test_cell_value(State_VGB(rhoUz_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Bxpnt','bxpnt','Bxtest','bxtest')                      
       tmp1_BLK = B0xcell_BLK + State_VGB(Bx_,:,:,:,:) 
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Bypnt','bypnt','Bytest','bytest')                      
       tmp1_BLK = B0ycell_BLK + State_VGB(By_,:,:,:,:) 
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Bzpnt','bzpnt','Bztest','bztest')                      
       tmp1_BLK = B0zcell_BLK + State_VGB(Bz_,:,:,:,:) 
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Epnt','epnt','Etest','etest')                        
       LogVar(iVarTot) = test_cell_value(E_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Ppnt','ppnt','Ptest','ptest', &
         'Pthpnt','pthpnt','Pthtest','pthtest')                        
       LogVar(iVarTot) = test_cell_value(State_VGB(P_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)

       ! EXTRA MHD variables at a single point given by Itest,Jtest,Ktest,
       !    BLKtest, PROCtest
    case('Uxpnt','uxpnt','Uxtest','uxtest')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = &
               State_VGB(rhoUx_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Uypnt','uypnt','Uytest','uytest')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUy_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Uzpnt','uzpnt','Uztest','uztest')
       do iBLK=1,nBlock
          if (unusedBLK(iBLK)) cycle
          tmp1_BLK(:,:,:,iBLK) = State_VGB(rhoUz_,:,:,:,iBLK)/&
               State_VGB(rho_,:,:,:,iBLK)
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('B1xpnt','b1xpnt','B1xtest','b1xtest')
       LogVar(iVarTot) = test_cell_value(&
            State_VGB(Bx_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('B1ypnt','b1ypnt','B1ytest','b1ytest')
       LogVar(iVarTot) = test_cell_value(&
            State_VGB(By_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('B1zpnt','b1zpnt','B1ztest','b1ztest')
       LogVar(iVarTot) = test_cell_value(&
            State_VGB(Bz_,:,:,:,:),-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Jxpnt','jxpnt','Jxtest','jxtest')
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle           
          !^CFG IF CARTESIAN BEGIN
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,iBLK)=0.5*(&
               (State_VGB(Bz_,0:nI+1,1:nJ+2,0:nK+1,iBLK)-&
               State_VGB(Bz_,0:nI+1,-1:nJ,0:nK+1,iBLK))&
               /dy_BLK(iBLK) &
               -(State_VGB(By_,0:nI+1,0:nJ+1,1:nK+2,iBLK)-&
               State_VGB(By_,0:nI+1,0:nJ+1,-1:nK,iBLK))&
               /dz_BLK(iBLK))
          !^CFG END CARTESIAN
          !call covar_curlb_plotvar(1,iBLK,tmp1_BLK(:,:,:,iBLK))!^CFG IF NOT CARTESIAN
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Jypnt','jypnt','Jytest','jytest')
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle
          !^CFG IF CARTESIAN BEGIN
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,iBLK)=0.5*(&
               (State_VGB(Bx_,0:nI+1,0:nJ+1,1:nK+2,iBLK)-&
               State_VGB(Bx_,0:nI+1,0:nJ+1,-1:nK,iBLK))&
               /dz_BLK(iBLK) &
               -(State_VGB(Bz_,1:nI+2,0:nJ+1,0:nK+1,iBLK)-&
               State_VGB(Bz_,-1:nI,0:nJ+1,0:nK+1,iBLK))&
               /dx_BLK(iBLK))
          !^CFG END CARTESIAN
          !call covar_curlb_plotvar(2,iBLK,tmp1_BLK(:,:,:,iBLK))!^CFG IF NOT CARTESIAN
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case('Jzpnt','jzpnt','Jztest','jztest')
       do iBLK=1,nBlock
          if(unusedBLK(iBLK))cycle
          !^CFG IF CARTESIAN BEGIN
          tmp1_BLK(0:nI+1,0:nJ+1,0:nK+1,iBLK)=0.5*(&
               (State_VGB(By_,1:nI+2,0:nJ+1,0:nK+1,iBLK)-&
               State_VGB(By_,-1:nI,0:nJ+1,0:nK+1,iBLK))&
               /dx_BLK(iBLK) &
               -(State_VGB(Bx_,0:nI+1,1:nJ+2,0:nK+1,iBLK)-&
               State_VGB(Bx_,0:nI+1,-1:nJ,0:nK+1,iBLK))&
               /dy_BLK(iBLK))
          !^CFG END CARTESIAN
          !call covar_curlb_plotvar(3,iBLK,tmp1_BLK(:,:,:,iBLK))!^CFG IF NOT CARTESIAN
       end do
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)

!!$!^CFG  IF RAYTRACE BEGIN
       ! BASIC RAYTRACE variables at a single point given by Itest,Jtest,Ktest
       !    BLKtest, PROCtest
       ! note that here we do not need the ghost cells so there is no
       ! interpolation as there is below for the satellite variables.
    case ('theta1pnt','theta1test')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(1,1,1:nI,1:nJ,1:nK,:)
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case ('theta2pnt','theta2test')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(1,2,1:nI,1:nJ,1:nK,:)
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case ('phi1pnt','phi1test')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(2,1,1:nI,1:nJ,1:nK,:)
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case ('phi2pnt','phi2test')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(2,2,1:nI,1:nJ,1:nK,:)
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
    case ('statuspnt','statustest')
       tmp1_BLK(1:nI,1:nJ,1:nK,:)=ray(3,1,1:nI,1:nJ,1:nK,:)
       LogVar(iVarTot) = test_cell_value(tmp1_BLK,-1,nI+2,-1,nJ+2,-1,nK+2)
!!$!^CFG END RAYTRACE

       ! Ionosphere values                            !^CFG IF IONOSPHERE BEGIN
    case('cpcpn','cpcp_n','cpcp_north','cpcpnorth',&
         'cpcps','cpcp_s','cpcp_south','cpcpsouth')
       LogVar(iVarTot) = logvar_ionosphere(NameLogVar) !^CFG END IONOSPHERE

       ! Flux values
    case('Aflx','aflx')   
       ! just to check that the area is being computed correctly
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)
          tmp1_BLK = 1.0
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log, tmp1_BLK)
       end do
    case('rhoflx')
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)
          tmp1_BLK = (State_VGB(rhoUx_,:,:,:,:)*x_BLK+&
               State_VGB(rhoUy_,:,:,:,:)*y_BLK+State_VGB(rhoUz_,:,:,:,:)*z_BLK)/R_BLK
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log, tmp1_BLK)
       end do
    case('dstflx')
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)
          tmp1_BLK = State_VGB(Bz_,:,:,:,:)
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log, tmp1_BLK) / &
               (4*cPi*R_log**2)
       end do
    case('Bflx','bflx')
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)
          tmp1_BLK = ( (State_VGB(Bx_,:,:,:,:)+B0xCell_BLK)*x_BLK + &
               (State_VGB(By_,:,:,:,:)+B0yCell_BLK)*y_BLK + &
               (State_VGB(Bz_,:,:,:,:)+B0zCell_BLK)*z_BLK )/R_BLK
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log, tmp1_BLK)
       end do
    case('B2flx','b2flx')
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)           
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))cycle           
             tmp1_BLK(:,:,:,iBLK) = (State_VGB(Bx_,:,:,:,iBLK)+&
                  B0xCell_BLK(:,:,:,iBLK))**2 + &
                  (State_VGB(By_,:,:,:,iBLK)+&
                  B0yCell_BLK(:,:,:,iBLK))**2 + &
                  (State_VGB(Bz_,:,:,:,iBLK)+&
                  B0zCell_BLK(:,:,:,iBLK))**2 
             tmp1_BLK(:,:,:,iBLK) = 0.5*tmp1_BLK(:,:,:,iBLK)* &
                  ( State_VGB(rhoUx_,:,:,:,iBLK)*x_BLK(:,:,:,iBLK) &
                  +State_VGB(rhoUy_,:,:,:,iBLK)*y_BLK(:,:,:,iBLK) &
                  +State_VGB(rhoUz_,:,:,:,iBLK)*z_BLK(:,:,:,iBLK))/ &
                  (State_VGB(rho_,:,:,:,iBLK)*R_BLK(:,:,:,iBLK))
          end do
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log,tmp1_BLK)
       end do
    case('pvecflx','Pvecflx','PVecflx','PVecFlx','pVecFlx','pvecFlx')
       iVarTot = iVarTot - 1
       do iR=1,nlogR
          iVarTot = iVarTot + 1
          R_log = logRvalues(iR)
          do iBLK=1,nBlock
             if(unusedBLK(iBLK))cycle           
             tmp1_BLK(:,:,:,iBLk) = (State_VGB(Bx_,:,:,:,iBLk)+B0xCell_BLK(:,:,:,iBLk))**2 + &
                  (State_VGB(By_,:,:,:,iBLk)+B0yCell_BLK(:,:,:,iBLk))**2 + &
                  (State_VGB(Bz_,:,:,:,iBLk)+B0zCell_BLK(:,:,:,iBLk))**2 
             tmp2_G = (State_VGB(Bx_,:,:,:,iBLk)+B0xCell_BLK(:,:,:,iBLk))* &
                  State_VGB(rhoUx_,:,:,:,iBLk) + &
                  (State_VGB(By_,:,:,:,iBLk)+B0yCell_BLK(:,:,:,iBLk))* &
                  State_VGB(rhoUy_,:,:,:,iBLk) + &
                  (State_VGB(Bz_,:,:,:,iBLk)+B0zCell_BLK(:,:,:,iBLk))* &
                  State_VGB(rhoUz_,:,:,:,iBLk) 
             tmp1_BLK(:,:,:,iBLk) = ( ( tmp1_BLK(:,:,:,iBLk)*State_VGB(rhoUx_,:,:,:,iBLk) &
                  -tmp2_G*(State_VGB(Bx_,:,:,:,iBLk)+ &
                  B0xCell_BLK(:,:,:,iBLk)))*X_BLK(:,:,:,iBLk) &
                  +( tmp1_BLK(:,:,:,iBLk)*State_VGB(rhoUy_,:,:,:,iBLk) &
                  -tmp2_G*(State_VGB(By_,:,:,:,iBLk)+ &
                  B0yCell_BLK(:,:,:,iBLk)))*Y_BLK(:,:,:,iBLk) &  
                  +( tmp1_BLK(:,:,:,iBLk)*State_VGB(rhoUz_,:,:,:,iBLk) &
                  -tmp2_G*(State_VGB(Bz_,:,:,:,iBLk)+ &
                  B0zCell_BLK(:,:,:,iBLk)))*Z_BLK(:,:,:,iBLk) )&   
                  /(State_VGB(rho_,:,:,:,iBLk)*R_BLK(:,:,:,iBLk))
          end do
          LogVar(iVarTot) = integrate_flux_sph(nProc, R_log,tmp1_BLK)
       end do


       ! simple circular integrals
    case('E2dflx','e2dflx','E2Dflx','E2DFlx')
       !this is the azimuthal component of the electric field 
       !integrated around a circle
       iVarTot = iVarTot-1
       do iR=1,nlogR
          iVarTot = iVarTot+1
          R_log = logRvalues(iR)
          do iBLK = 1,nBlock
             if(unusedBLK(iBLK))cycle
             tmp1_BLK(:,:,:,iBLK) = -(-( State_VGB(rhoUy_,:,:,:,iBLK)* &
                  (State_VGB(Bz_,:,:,:,iBLk)+B0zCell_BLK(:,:,:,iBLk)) &
                  -State_VGB(rhoUz_,:,:,:,iBLK)* &
                  (State_VGB(By_,:,:,:,iBLk)+B0yCell_BLK(:,:,:,iBLk)) )* &
                  Y_BLK(:,:,:,iBLK) &
                  +( State_VGB(rhoUz_,:,:,:,iBLK)* &
                  (State_VGB(Bx_,:,:,:,iBLk)+B0xCell_BLK(:,:,:,iBLk)) &
                  -State_VGB(rhoUx_,:,:,:,iBLK)* &
                  (State_VGB(Bz_,:,:,:,iBLk)+B0zCell_BLK(:,:,:,iBLk)) )* &
                  X_BLK(:,:,:,iBLK)  )/ &
                  (State_VGB(rho_,:,:,:,iBLK)*R_BLK(:,:,:,iBLK))
          end do
          LogVar(iVarTot) = integrate_flux_circ(nProc,R_log,0.0,tmp1_BLK)
       end do

       ! OTHER VALUES
    case('dt', 'DT')
       LogVar(iVarTot) = dt

       ! DEFAULT
    case default
       if (UseUserLogFiles) then                 !^CFG IF USERFILES BEGIN
          call user_get_log_var( LogVar(iVarTot), NameLogVar )
       else                                      !^CFG END USERFILES
          LogVar(iVarTot) = -7777.
          call write_myname;
          write(*,*) 'WARNING in set_logvar: unknown logvarname=',NameLogVar
       endif                                     !^CFG IF USERFILES
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

       ! Basic MHD variables at a single point along a satellite trajectory
    case('rho')
       LogVar(iVarTot) = StateSat_V(Rho_)
    case('rhoUx','rhoux','mx')
       LogVar(iVarTot) = StateSat_V(RhoUx_)
    case('rhoUy','rhouy','my')
       LogVar(iVarTot) = StateSat_V(RhoUy_)
    case('rhoUz','rhouz','mz')
       LogVar(iVarTot) = StateSat_V(RhoUz_)
    case('Bx','bx')
       LogVar(iVarTot) = StateSat_V(Bx_)+B0Sat_D(1)
    case('By','by')
       LogVar(iVarTot) = StateSat_V(By_)+B0Sat_D(2)
    case('Bz','bz')
       LogVar(iVarTot) = StateSat_V(Bz_)+B0Sat_D(3)
    case('E','e')
       LogVar(iVarTot) = 0.5*(&
            sum(StateSat_V(RhoUx_:RhoUz_)**2)/StateSat_V(rho_) + &
            sum((StateSat_V(Bx_:Bz_)+B0Sat_D)**2))
    case('P','p','Pth','pth')
       LogVar(iVarTot) = StateSat_V(P_)

       ! EXTRA MHD variables at a single point along a satellite trajectory,
    case('Ux','ux')
       LogVar(iVarTot) = StateSat_V(RhoUx_)/StateSat_V(Rho_)
    case('Uy','uy')
       LogVar(iVarTot) = StateSat_V(RhoUy_)/StateSat_V(Rho_)
    case('Uz','uz')
       LogVar(iVarTot) = StateSat_V(RhoUz_)/StateSat_V(Rho_)
    case('B1x','b1x')
       LogVar(iVarTot) = StateSat_V(Bx_)
    case('B1y','b1y')
       LogVar(iVarTot) = StateSat_V(By_)
    case('B1z','b1z')
       LogVar(iVarTot) = StateSat_V(Bz_)
    case('Jx','jx')
       LogVar(iVarTot) = StateSat_V(nVar+1)
    case('Jy','jy')
       LogVar(iVarTot) = StateSat_V(nVar+2)
    case('Jz','jz')
       LogVar(iVarTot) = StateSat_V(nVar+3)
    case('weight')
       LogVar(iVarTot) = StateSat_V(0)
    case('order')
       if(abs(StateSat_V(0)-1)<1.e-5)then
          LogVar(iVarTot) = 2
       else
          LogVar(iVarTot) = 1
       end if
    case default
       LogVar(iVarTot) = -777.0
       if(iProc==0)write(*,*)'WARNING in var_sat: unknown variable ',&
            NameLogVar,' for iSat = ',iSat
    end select
  end subroutine set_sat_var

end subroutine set_logvar

!==============================================================================
subroutine normalize_logvar(nlogvar,logvarnames,nlogR,logRvalues,nlogTot,LogVar)

  use ModPhysics
  implicit none

  integer, intent(in) :: Nlogvar, nlogR, nlogTot
  character (LEN=10), intent(in) :: logvarnames(Nlogvar)
  real, intent(out) :: LogVar(nLogTot)
  real, intent(in) :: logRvalues(nlogR)
  character (len=10) :: s

  integer :: iVar, iR, iVarTot


  !-------------------------------------------------------------------------

  iVarTot = 0
  do iVar=1,nLogVar

     iVarTot = iVarTot+1
     s=logvarnames(iVar)
     select case(s)

        ! BASIC MHD variables
     case('rho','rhopnt','rhotest') 
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_n
     case('rhoUx','rhoux','mx','rhoUy','rhouy','my','rhoUz','rhouz','mz', &
          'rhoUxpnt','rhouxpnt','mxpnt','rhoUypnt','rhouypnt','mypnt', &
          'rhoUzpnt','rhouzpnt','mzpnt', &
          'rhoUxtest','rhouxtest','mxtest','rhoUytest','rhouytest','mytest', &
          'rhoUztest','rhouztest','mztest')
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_rhoU
     case('Bx','bx','By','by','Bz','bz', &
          'Bxpnt','bxpnt','Bypnt','bypnt','Bzpnt','bzpnt', &
          'B1xpnt','b1xpnt','B1ypnt','b1ypnt','B1zpnt','b1zpnt', &
          'Bxtest','bxtest','Bytest','bytest','Bztest','bztest', &
          'B1xtest','b1xtest','B1ytest','b1ytest','B1ztest','b1ztest',&
          'dst','Dst','DST','dstdivb','DstDivb','DSTDIVB')
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_B
     case('E','e','Epnt','epnt','Etest','etest')
        LogVar(iVarTot) = LogVar(iVarTot)*unitUSER_energydens
     case('P','p','Pth','pth','Ppnt','ppnt','Pthpnt','pthpnt', &
          'Ptest','ptest','Pthtest','pthtest')
        LogVar(iVarTot) = LogVar(iVarTot)*unitUSER_p

        ! Extra MHD variables
     case('Ux','ux','Uy','uy','Uz','uz', &
          'Uxpnt','uxpnt','Uypnt','uypnt','Uzpnt','uzpnt', &
          'Uxtest','uxtest','Uytest','uytest','Uztest','uztest')
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_U
     case('Ekinx','ekinx','Ekiny','ekiny','Ekinz','ekinz','ekin')
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_rho*unitUSER_U**2
     case('Jx','jx','Jy','jy','Jz','jz', &
          'Jxpnt','jxpnt','Jypnt','jypnt','Jzpnt','jzpnt', &
          'Jxtest','jxtest','Jytest','jytest','Jztest','jztest')
        LogVar(iVarTot)= LogVar(iVarTot)*unitUSER_J

!!$!^CFG  IF RAYTRACE BEGIN
        ! Basic RAY TRACE variables
     case('theta1','theta2','phi1','phi2','status', &
          'theta1pnt','theta2pnt','phi1pnt','phi2pnt','statuspnt', &
          'theta1test','theta2test','phi1test','phi2test','statustest')
!!$!^CFG END RAYTRACE

        ! Maximum and minimum values
     case('Pmin','pmin','Pmax','pmax')
        LogVar(iVarTot) = LogVar(iVarTot)*unitUSER_p


        ! Ionosphere values                           !^CFG IF IONOSPHERE BEGIN
     case('cpcpn','cpcp_n','cpcp_north','cpcpnorth',&
          'cpcps','cpcp_s','cpcp_south','cpcpsouth')
        LogVar(iVarTot) = LogVar(iVarTot)/1000.0      !^CFG END IONOSPHERE

        ! Flux values
     case('Aflx','aflx')   
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitUSER_x**2)
        end do
     case('rhoflx')
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitSI_n*unitSI_U*unitSI_x**2)
        end do
     case('dstflx')
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*unitUSER_B
        end do
     case('Bflx','bflx')
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitSI_B*unitSI_x**2)
        end do
     case('B2flx','b2flx')
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitSI_Poynting*unitSI_x**2)
        end do
     case('pvecflx','Pvecflx','PVecflx','PVecFlx','pVecFlx','pvecFlx')
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitSI_Poynting*unitSI_x**2)
        end do

        ! simple circular integrals
     case('E2dflx','e2dflx','E2Dflx','E2DFlx')
        !this is the azimuthal component of the electric field integrated around a circle
        iVarTot = iVarTot - 1
        do iR=1,nlogR
           iVarTot = iVarTot + 1
           LogVar(iVarTot) = LogVar(iVarTot)*(unitSI_electric*unitSI_x)
        end do

        ! OTHER VALUES
     case('dt', 'DT')
        LogVar(iVarTot) = LogVar(iVarTot)*unitUSER_t


        ! DEFAULT FOR A BAD SELECTION

     case default
        ! no normalization

     end select
  end do ! iVar
end subroutine normalize_logvar


!==============================================================================

real function integrate_flux_sph(qnum,qrad,qa)

  ! This function calculates the integral of the incomming variable qa
  ! over the surface defined by the radius qrad.  This routine assumes that the
  ! value passed in qa is the radial component of the transport only.

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
  real, intent(in) :: qrad 

  ! Local variables
  integer :: iError
  real :: Flux_sum, GLOBAL_flux_sum, flux_to_add 

  ! Indices and coordinates
  integer :: iBLK,i,j,i1,i2,j1,j2,k1,k2
  integer :: ntheta_sph,nphi_sph
  real :: xx1,xx2,yy1,yy2,zz1,zz2,minRblk, maxRblk
  real :: x,y,z,xx,yy,zz,dx1,dx2,dy1,dy2,dz1,dz2
  real :: dtheta_sph, dphi_sph,theta_sph,phi_sph

  logical :: oktest,oktest_me
  !---------------------------------------------------------------------------

  Flux_sum = 0.0
  GLOBAL_Flux_sum = 0.0

  call set_oktest('integrate_flux_sph',oktest,oktest_me)

  ! the angular resolution of the integral is hard coded
  ! in degree
  dtheta_sph = 0.5
  dphi_sph   = 0.5

  ! now compute the number of points that correspond to the above
  ! intervals and then recompute the intervals in case an even multiple
  ! was not chosen.
  ntheta_sph = 180.0/dtheta_sph
  nphi_sph   = 360.0/dphi_sph
  dtheta_sph = 180.0/real(ntheta_sph)
  dphi_sph   = 360.0/real(nphi_sph)
  !convert to radians
  dtheta_sph = dtheta_sph*cPi/180.0
  dphi_sph   = dphi_sph*cPi/180.0

  if (oktest_me) write(*,*) 'ntheta,nphi,dtheta,dphi',ntheta_sph,nphi_sph,dtheta_sph,dphi_sph

  ! Sum all cells within range

  do iBLK = 1,nBlock

     if (unusedBLK(iBLK)) CYCLE
     ! get the max and min radial distance for this block so that we can check
     ! whether or not this block contibutes to the sum.
     xx1 = 0.50*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
     xx2 = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
     yy1 = 0.50*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
     yy2 = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
     zz1 = 0.50*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
     zz2 = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))
     minRblk = sqrt((min(abs(xx1),abs(xx2)))**2 + &
          (min(abs(yy1),abs(yy2)))**2 + &
          (min(abs(zz1),abs(zz2)))**2)
     maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
          (max(abs(yy1),abs(yy2)))**2 + &
          (max(abs(zz1),abs(zz2)))**2)
     if (minRblk > qrad .or. maxRblk < qrad ) CYCLE


     do i=1,ntheta_sph
        theta_sph = (i-.5)*dtheta_sph
        do j=1,nphi_sph
           phi_sph = j*dphi_sph

           ! get the xyz coordinates
           x = qrad*sin(theta_sph)*cos(phi_sph)
           y = qrad*sin(theta_sph)*sin(phi_sph)
           z = qrad*cos(theta_sph)

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
              flux_to_add = &
                   dx1*(   dy1*(   dz1*qa(i2,j2,k2,iBLK)   &
                   +               dz2*qa(i2,j2,k1,iBLK))  &
                   +       dy2*(   dz1*qa(i2,j1,k2,iBLK)   &
                   +               dz2*qa(i2,j1,k1,iBLK))) &
                   +dx2*(  dy1*(   dz1*qa(i1,j2,k2,iBLK)   &
                   +               dz2*qa(i1,j2,k1,iBLK))  &
                   +       dy2*(   dz1*qa(i1,j1,k2,iBLK)   &
                   +               dz2*qa(i1,j1,k1,iBLK)))

              flux_to_add = flux_to_add* &
                   2.0*(qrad**2)*dphi_sph*sin(dtheta_sph/2.0)*sin(theta_sph)
              Flux_sum = Flux_sum + flux_to_add

           end if

        end do
     end do

  end do

  if(qnum>1)then
     call MPI_allreduce(Flux_sum, GLOBAL_flux_sum, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
  else
     GLOBAL_flux_sum = Flux_sum
  end if
  integrate_flux_sph = GLOBAL_flux_sum

end function integrate_flux_sph


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
  use ModMain, ONLY: xTest,yTest,zTest
  use ModAdvance, ONLY: State_VGB
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  implicit none
  real :: State_V(0:nVar+3)

  State_VGB(Bx_,:,:,:,:) = y_BLK
  State_VGB(By_,:,:,:,:) = z_BLK
  State_VGB(Bz_,:,:,:,:) = x_BLK

  call get_point_data(0.0,(/xTest,yTest,zTest/),1,nVar+3,State_V)
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
