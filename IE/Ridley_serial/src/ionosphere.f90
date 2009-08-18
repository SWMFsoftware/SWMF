!^CFG COPYRIGHT UM
!
! This is the ionospheric solver for the BATS-R-US MHD code.
! Initially written by Clinton Groth as a multigrid solver.
! Modified by Aaron Ridley:
!   Initial Modifications:
!     - Changed to a GMRES solver, Sept. 99
!   - Added f107 dependence to the conductance
!   - Added a simple auroral oval
!

!====================================
!                                   |
!         Ionosphere Model          |
!                                   |
!====================================

subroutine ionosphere(iter, iAction)
  
  !\
  ! ionosphere driver (main or controlling) routine.
  !/
  
  use ModIonosphere
  use IE_ModMain
  use ModProcIE

  implicit none
  
  integer, intent(in) :: iter, iAction
  integer :: iModel
  real :: f107
  
  logical :: oktest, oktest_me
  
  call CON_set_do_test('ionosphere',oktest,oktest_me)
  if(oktest)write(*,*)'Ionosphere starting with action, me=',iAction, iProc

  iModel = conductance_model
  
  f107 = f107_flux
  
  select case (iAction)
     
  case (1)
     !\
     ! Create fine grids for north and south 
     ! hemisphere ionospheric solutions
     !/
     call ionosphere_fine_grid
     call ionosphere_init
     
  case(3)
     call ionosphere_read_restart_file(iter)
     
  case (5)
     !\
     ! Create a restart solution file containing the
     ! FAC driving the ionospheric solution.
     !/
     call ionosphere_write_restart_file(iter)

  case default
     write(*,*)'IE ionosphere: ERROR, iAction = ',iAction
     call CON_stop('IE ionosphere: ERROR invalid iAction')
     
  end select
  
end subroutine ionosphere

!*************************************************************************
!
! INITIALIZATION Routines
!
!*************************************************************************

subroutine ionosphere_fine_grid
  !\
  ! This routine sets the fine grid meshes for the
  ! northern and southern hemispheres.
  !/
  use ModProcIE
  use IE_ModIO, ONLY: iUnitOut, write_prefix
  use IE_ModMain, ONLY: conductance_model
  use ModIonosphere
  use CON_planet, ONLY: get_planet, lNamePlanet
  use ModNumConst, ONLY: cHalfPi, cTwoPi
  implicit none

  integer :: i,j
  real :: dTheta_l, dPsi_l
  character(len=*), parameter :: NameSub='ionosphere_fine_grid'
  character(len=lNamePlanet) :: NamePlanet
  !------------------------------------------------------------------------
  dTheta_l = cHalfPi/(IONO_nTheta-1)
  dPsi_l   = cTwoPi/(IONO_nPsi-1)

  call get_planet(NamePlanetOut = NamePlanet, RadiusPlanetOut = IONO_Radius)

  select case (NamePlanet)
  case ('EARTH')
     IONO_Height = 110000.00
  case ('SATURN')
     IONO_Height = 1000000.00
  case ('JUPITER')
     IONO_Height = 1000000.00
  case default            
     IONO_Height = IONO_Radius * 0.02  ! Most planets should work
  end select

  Radius = IONO_Radius + IONO_Height

  if(iProc==0)then
     call write_prefix; write(iUnitOut,*)
     call write_prefix; write(iUnitOut,*)'   Ionosphere Solution Parameters'
     call write_prefix; write(iUnitOut,*)'   ------------------------------'
     call write_prefix; write(iUnitOut,*)
     call write_prefix; write(iUnitOut,'(10X,A,I5,A,I5)') &
          'Ionospheric grids   : ',IONO_nTheta,' x ',IONO_nPsi
     call write_prefix; write(iUnitOut,'(10X,A,f7.1,A)') &
          'Height of ionosphere: ',IONO_Height/1000,' km'
     call write_prefix; write(iUnitOut,'(10X,A)',ADVANCE='NO') &
          'Conductance model   : '
     select case(conductance_model)
     case(0)
        write(iUnitOut,'(a)') "Constant Pedersen (no HALL)"
     case(1)
        write(iUnitOut,'(a)') "Constant Pedersen and Hall"
     case(2)
        write(iUnitOut,'(a)') "Solar EUV Only"
     case(3)
        write(iUnitOut,'(a)') "Supersmooth oval"
     case(4)
        write(iUnitOut,'(a)') "Restricted oval"
     case(5)
        write(iUnitOut,'(a)') "Realistic oval"
     case(6)
        write(iUnitOut,'(a)') "More Realistic oval"
     case default
        call CON_stop(NameSub//" IE_ERROR invalid conductance_model")
     end select
     call write_prefix; write(iUnitOut,*)
  end if

  do j = 1, IONO_nPsi
     IONO_NORTH_Theta(1,j) = IONO_Theta_0
     IONO_NORTH_Psi(1,j) = real(j-1)*dPsi_l
     IONO_NORTH_X(1,j) = sin(IONO_NORTH_Theta(1,j))* &
          cos(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Y(1,j) = sin(IONO_NORTH_Theta(1,j))* &
          sin(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Z(1,j) = cos(IONO_NORTH_Theta(1,j))

     do i = 2, IONO_nTheta
        IONO_NORTH_Theta(i,j) = cHalfPi*real(i-1)/real(IONO_nTheta-1)
        IONO_NORTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
        IONO_NORTH_X(i,j) = sin(IONO_NORTH_Theta(i,j))* &
             cos(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Y(i,j) = sin(IONO_NORTH_Theta(i,j))* &
             sin(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Z(i,j) = cos(IONO_NORTH_Theta(i,j))
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        IONO_SOUTH_Theta(i,j) = cPi - IONO_NORTH_Theta(IONO_nTheta-(i-1),j)
        IONO_SOUTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
        IONO_SOUTH_X(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
             cos(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Y(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
             sin(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Z(i,j) = cos(IONO_SOUTH_Theta(i,j))
     end do
  end do

  !
  ! dPsi is going to be defined as 2*dPsi for all of the points:
  !

  do j = 2, IONO_nPsi-1
     dPsi_North(j) = IONO_NORTH_Psi(1,j+1)-IONO_NORTH_Psi(1,j-1)
  enddo
  dPsi_North(1)    = IONO_NORTH_Psi(1,2) - &
       (IONO_NORTH_Psi(1,IONO_nPsi-1) - cTwoPi)
  dPsi_North(IONO_nPsi) = dPsi_North(1)

  !
  ! dTheta is going to be defined as 2*dTheta for all but the top and bottom
  !   points.  As these locations, we switch to 1st order
  !

  do i = 2, IONO_nTheta-1
     dTheta_North(i)        = IONO_NORTH_Theta(i+1,1) - IONO_NORTH_Theta(i-1,1)
  enddo
  dTheta_North(1)           = IONO_NORTH_Theta(2,1) - IONO_NORTH_Theta(1,1)
  dTheta_North(IONO_nTheta) =                                                &
       IONO_NORTH_Theta(IONO_nTheta,1)-IONO_NORTH_Theta(IONO_nTheta-1,1)

  !
  ! dPsi is going to be defined as 2*dPsi for all of the points:
  !

  do j = 2, IONO_nPsi-1
     dPsi_South(j)      = IONO_SOUTH_Psi(1,j+1)-IONO_SOUTH_Psi(1,j-1)
  enddo
  dPsi_South(1)         = IONO_SOUTH_Psi(1,2)-                               &
       (IONO_SOUTH_Psi(1,IONO_nPsi-1) - cTwoPi)
  dPsi_South(IONO_nPsi) = dPsi_South(1)

  !
  ! dTheta is going to be defined as 2*dTheta for all but the top and bottom
  !   points.  As these locations, we switch to 1st order
  !

  do i = 2, IONO_nTheta-1
     dTheta_South(i)        = IONO_SOUTH_Theta(i+1,1) - IONO_SOUTH_Theta(i-1,1)
  enddo
  dTheta_South(1)           = IONO_SOUTH_Theta(2,1) - IONO_SOUTH_Theta(1,1)
  dTheta_South(IONO_nTheta) = IONO_SOUTH_Theta(IONO_nTheta,1) -              &
       IONO_SOUTH_Theta(IONO_nTheta-1,1)

end subroutine ionosphere_fine_grid

!=============================================================================

subroutine ionosphere_init
  !\
  ! This routine initializes the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  implicit none

  IONO_NORTH_PHI = 0.00
  IONO_NORTH_JR = 0.00 
  IONO_NORTH_Jx = 0.00
  IONO_NORTH_Jy = 0.00
  IONO_NORTH_Jz = 0.00
  IONO_NORTH_Ex = 0.00
  IONO_NORTH_Ey = 0.00
  IONO_NORTH_Ez = 0.00
  IONO_NORTH_ETh = 0.00
  IONO_NORTH_EPs = 0.00
  IONO_NORTH_Ux = 0.00
  IONO_NORTH_Uy = 0.00
  IONO_NORTH_Uz = 0.00

  IONO_SOUTH_PHI = 0.00
  IONO_SOUTH_JR = 0.00
  IONO_SOUTH_Jx = 0.00
  IONO_SOUTH_Jy = 0.00
  IONO_SOUTH_Jz = 0.00
  IONO_SOUTH_Ex = 0.00
  IONO_SOUTH_Ey = 0.00
  IONO_SOUTH_Ez = 0.00
  IONO_SOUTH_ETh = 0.00
  IONO_SOUTH_EPs = 0.00
  IONO_SOUTH_Ux = 0.00
  IONO_SOUTH_Uy = 0.00
  IONO_SOUTH_Uz = 0.00  

  IONO_NORTH_TGCM_JR = 0.00                                !^CFG IF TIEGCM
  IONO_SOUTH_TGCM_JR = 0.00                                !^CFG IF TIEGCM

  SAVE_NORTH_SigmaH = 0.00
  SAVE_SOUTH_SigmaH = 0.00

  IONO_NORTH_EFlux = 0.00
  IONO_SOUTH_EFlux = 0.00

  IONO_NORTH_Ave_E = 0.00
  IONO_SOUTH_Ave_E = 0.00

end subroutine ionosphere_init

!*************************************************************************
!
! INPUT/OUTPUT Routines
!
!*************************************************************************

!^CFG  IF TIEGCM BEGIN
!-------------------------------------------------------------------------
! write_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine write_timegcm_file(iter, phi_north, phi_south,   &
     eflux_north, eflux_south, avee_north, avee_south)

  use ModIonosphere
  use IE_ModIo
  use IE_ModMain,ONLY: Time_Array
  implicit none

  integer, intent(in) :: iter
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::             &
       phi_north, phi_south, eflux_north, eflux_south,      &
       avee_north, avee_south
  real psi_offset

  integer :: i, j
  character (len=4), Parameter :: IO_ext=".dat"

  write(*,*) '=> Writing output datafiles for TIMEGCM.'

  open(unit=iUnit,file=trim(NameIonoDir)//"MHD_to_TIMEGCM"//IO_ext, &
       status="unknown")

!  Write Year, month, day, hour, minute, second

  write(iUnit,fmt="(6I5)") (Time_Array(i),i=1,6)

! Write number of thetas and psis

  write(iUnit,fmt="(2I10)") IONO_nTheta*2, IONO_nPsi

! Silliness to occur here:
!   This ionospheric code starts with noon at 0.0 longitude.  This is 
!   different than almost all other ionospheres, so we need to shift
!   the output by 180.0 degrees.
!   Also, we want 0.0 and 360.0 to both be printed - but the new ones,
!   so we have to ignore the old 0.0 and repeat the old 180.0.

  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = cRadToDeg*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(iUnit,fmt="(1P,5(E13.5))")  &
             cRadToDeg*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = cRadToDeg*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(iUnit,fmt="(1P,5(E13.5))")  &
             cRadToDeg*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
  end do
  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = cRadToDeg*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(iUnit,fmt="(1P,5(E13.5))")  &
             cRadToDeg*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = cRadToDeg*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(iUnit,fmt="(1P,5(E13.5))")  &
             cRadToDeg*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
  end do

  close(iUnit)
 
end subroutine write_timegcm_file
!^CFG END TIEGCM

!^CFG  IF TIEGCM BEGIN
!-------------------------------------------------------------------------
! read_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine read_timegcm_file

  use ModIonosphere
  use IE_ModIo
  implicit none

  integer :: i, j, iter, nTheta, nPsi
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=100) :: TIMEGCM_file
  real :: dum_lat, dum_lon, min_lat, max_lat, max_fac
  integer, dimension(1:6) :: time_array
  real :: tgcm_time

  min_lat = 60.0
  max_lat = 85.0
  max_fac = 0.0

  write(*,*) '=> Reading output file from the TIMEGCM.'

  TIMEGCM_file = trim(NameIonoDir)//"TIMEGCM_to_MHD.dat"

  open(unit=iUnit,file=TIMEGCM_file, status="old")

  read(iUnit,fmt="(6I5)") (time_array(i),i=1,6)

  read(iUnit,fmt="(2I10)") nTheta, nPsi

  do i = IONO_nTheta,1,-1
     do j = 1, IONO_nPsi
        read(iUnit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_SOUTH_TGCM_JR(i,j)

        if (dum_lat.lt.(90+min_lat)) then
           if (abs(iono_south_tgcm_jr(i,j)) > max_fac) &
                max_fac = abs(iono_south_tgcm_jr(i,j))
           IONO_SOUTH_TGCM_JR(i,j) = 0.0
        endif

        if (dum_lat.gt.(90+max_lat)) then
           IONO_SOUTH_TGCM_JR(i,j) = 0.0
        endif

     end do
  end do

  max_fac = 0.0
  do i = IONO_nTheta,1,-1
     do j = 1, IONO_nPsi
        read(iUnit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_NORTH_TGCM_JR(i,j)

        if (dum_lat.gt.(90-min_lat)) then
           if (abs(iono_north_tgcm_jr(i,j)) > max_fac) &
                max_fac = abs(iono_south_tgcm_jr(i,j))
           IONO_NORTH_TGCM_JR(i,j) = 0.0
        endif

        if (dum_lat.lt.(90-max_lat)) then
           IONO_NORTH_TGCM_JR(i,j) = 0.0
        endif
     end do
  end do

  close(iUnit)
 
end subroutine read_timegcm_file
!^CFG END TIEGCM

!-------------------------------------------------------------------------
! ionosphere_write_output
!
!
!
!-------------------------------------------------------------------------

subroutine IE_output

  use IE_ModIo
  use IE_ModMain
  use ModProcIE
  
  implicit none

  integer :: iFile
  logical :: DoSaveFile
  !--------------------------------------------------------------------------
  do iFile = 1, nFile
     DoSaveFile = .false.
     if(time_accurate .and. dt_output(ifile)>0.)then
        if(int(time_simulation/dt_output(ifile))>t_output_last(ifile))then
           t_output_last(ifile)=int(time_simulation/dt_output(ifile))
           DoSaveFile = .true.
        end if
     elseif(dn_output(ifile)>=0)then
        if(dn_output(ifile)==0)then
           DoSaveFile = .true.
        else if(nSolve>0.and.mod(nSolve,dn_output(ifile))==0)then
           DoSaveFile = .true.
        end if
     end if
     
     if(DoSaveFile)then
        if(iProc==0)      call ionosphere_write_output(iFile, 1)
        if(iProc==nProc-1)call ionosphere_write_output(iFile, 2)
     end if
  end do

end subroutine IE_output

subroutine ionosphere_write_output(iFile, iBlock)
  !\
  ! This routine writes out the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres to
  ! a data file.
  !/
  use ModIonosphere
  use IE_ModIo
  use IE_ModMain,ONLY: Time_Array, NSOLVE, Time_Accurate, Time_Simulation,&
       ThetaTilt
  use ModNumConst,ONLY: cRadToDeg

  implicit none

  integer, intent(in) :: ifile
  integer, intent(in) :: iBlock ! 1 for north 2 for south

  integer, parameter :: tec_type = 1
  integer, parameter :: idl_type = 2
  integer, parameter :: all_vars = 1
  integer, parameter :: min_vars = 2
  integer, parameter :: uam_vars = 3                     !^CFG  IF TIEGCM
  integer, parameter :: aur_vars = 4

  integer :: output_type, variables, year

  integer :: i, j
  character(len=*), parameter :: NameSub='IE_write_output'
  character(len=4) :: IO_ext
  character (len=23) :: textNandT
  !------------------------------------------------------------------------

  variables = min_vars
  output_type = idl_type

  select case(plot_vars(ifile))
  case('minimum')
     variables = min_vars
  case('maximum')
     variables = all_vars
  case('uam')                                   !^CFG  IF TIEGCM
     variables = uam_vars                       !^CFG  IF TIEGCM
  case('aur')
     variables = aur_vars
  end select

  select case(plot_form(ifile))
  case('idl')
     output_type = idl_type
  case('tec')
     output_type = tec_type
  case default
     call CON_stop(NameSub//' IE_ERROR invalid plot_form='//plot_form(ifile))
  end select

  if (output_type == idl_type) then
     IO_ext = ".idl"
  elseif (output_type == tec_type) then
     IO_ext = ".tec"
  endif

  if (.not.time_accurate) then

     write(NameFile,'(a,i6.6,a,i1,a)')trim(NameIonoDir)//"in",nSolve,&
          "_b",iBlock,IO_ext

     write(textNandT,'(a,i7.7)') "N=",nSolve

  else
     write(NameFile,'(a,3i2.2,"_",3i2.2,"_",i3.3,"_b",i1,a)') &
          trim(NameIonoDir)//"it",&
          mod(time_array(1),100),time_array(2:7),&
          iBlock,IO_ext

     write(textNandT,'(a,i7.7,a,i4.4,a,i2.2,a,i2.2)') &
          "N=",nSolve," T=", &
          int( Time_Simulation/3600.),":",&
          int((Time_Simulation-(3600.*int(Time_Simulation/3600.)))/60.),":",&
          int( Time_Simulation-(  60.*int(Time_Simulation/  60.)))

  end if

  open(unit=iUnit,file=NameFile,status="unknown")

  select case(iBlock)
  case(1) ! North writes header and data

     if (output_type == idl_type) then

        write(iUnit, '(a)')  'TITLE'
        write(iUnit, '(a,i4.4,"-",5(i2.2,"-"),i3.3,a,2f7.2,a)') &
             '  "BATSRUS: Ionospheric Potential Solution,  ', &
             Time_Array(1:7),', BTiltDeg=',ThetaTilt*180.0/cPi,0.,'"'

        write(iUnit, *) ' '

        write(iUnit, '(a)')  'NUMERICAL VALUES'
        select case(variables)
        case(min_vars)
           write(iUnit, '(I5,a)')  6, ' nvars'
        case(all_vars)
           write(iUnit, '(I5,a)') 22, ' nvars'
        case(uam_vars)                             !^CFG  IF TIEGCM
           write(iUnit, '(I5,a)')  9, ' nvars'     !^CFG  IF TIEGCM
        case(aur_vars)
           write(iUnit, '(I5,a)') 13, ' nvars'
        end select
        write(iUnit, '(I5,a)') IONO_nTheta, ' nTheta'
        write(iUnit, '(I5,a)')   IONO_nPsi, ' nPhi'

        write(iUnit, *) ' '

        write(iUnit, '(a)') 'VARIABLE LIST'

        select case(variables)
        case(min_vars)
           write(iUnit, '(I5,a)')  1, ' Theta [deg]'
           write(iUnit, '(I5,a)')  2, ' Psi [deg]'
           write(iUnit, '(I5,a)')  3, ' SigmaH [mhos]'
           write(iUnit, '(I5,a)')  4, ' SigmaP [mhos]'
           write(iUnit, '(I5,a)')  5, ' Jr [mA/m^2]'
           write(iUnit, '(I5,a)')  6, ' Phi [kV]'

        case(all_vars)
           write(iUnit, '(I5,a)')  1, ' X [Re]'
           write(iUnit, '(I5,a)')  2, ' Y [Re]'
           write(iUnit, '(I5,a)')  3, ' Z [Re]'
           write(iUnit, '(I5,a)')  4, ' Theta [deg]'
           write(iUnit, '(I5,a)')  5, ' Psi [deg]'
           write(iUnit, '(I5,a)')  6, ' SigmaH [mhos]'
           write(iUnit, '(I5,a)')  7, ' SigmaP [mhos]'
           write(iUnit, '(I5,a)')  8, ' E-Flux [W/m2]'
           write(iUnit, '(I5,a)')  9, ' Ave-E [eV]'
           write(iUnit, '(I5,a)') 10, ' Jr [mA/m^2]'
           write(iUnit, '(I5,a)') 11, ' Phi [kV]'
           write(iUnit, '(I5,a)') 12, ' Ex [mV/m]'
           write(iUnit, '(I5,a)') 13, ' Ey [mV/m]'
           write(iUnit, '(I5,a)') 14, ' Ez [mV/m]'
           write(iUnit, '(I5,a)') 15, ' Jx [microA/m2]'
           write(iUnit, '(I5,a)') 16, ' Jy [microA/m2]'
           write(iUnit, '(I5,a)') 17, ' Jz [microA/m2]'
           write(iUnit, '(I5,a)') 18, ' Ux [km/s]'
           write(iUnit, '(I5,a)') 19, ' Uy [km/s]'
           write(iUnit, '(I5,a)') 20, ' Uz [km/s]'
! add joule heating output --Yiqun Sep 2008
           write(iUnit, '(I5,a)') 21, ' JouleHeat [mW/m2]'
           write(iUnit, '(I5,a)') 22, ' IonNumFlux [/cm2/s]'

        case(uam_vars)                                  !^CFG  IF TIEGCM BEGIN
           write(iUnit, '(I5,a)')  1, ' Theta [deg]'
           write(iUnit, '(I5,a)')  2, ' Psi [deg]'
           write(iUnit, '(I5,a)')  3, ' SigmaH [mhos]'
           write(iUnit, '(I5,a)')  4, ' SigmaP [mhos]'
           write(iUnit, '(I5,a)')  5, ' Jr [mA/m^2]'
           write(iUnit, '(I5,a)')  6, ' Jr(NW) [mA/m^2]'
           write(iUnit, '(I5,a)')  7, ' E-Flux [W/m2]'
           write(iUnit, '(I5,a)')  8, ' Ave-E [eV]'
           write(iUnit, '(I5,a)')  9, ' Phi [kV]'        !^CFG END TIEGCM

        case(aur_vars)
           write(iUnit, '(I5,a)')  1, ' Theta [deg]'
           write(iUnit, '(I5,a)')  2, ' Psi [deg]'
           write(iUnit, '(I5,a)')  3, ' SigmaH [mhos]'
           write(iUnit, '(I5,a)')  4, ' SigmaP [mhos]'
           write(iUnit, '(I5,a)')  5, ' Jr [mA/m^2]'
           write(iUnit, '(I5,a)')  6, ' Phi [kV]'
           write(iUnit, '(I5,a)')  7, ' E-Flux [W/m2]'
           write(iUnit, '(I5,a)')  8, ' Ave-E [eV]'
           write(iUnit, '(I5,a)')  9, ' RT 1/B [1/T]'
           write(iUnit, '(I5,a)') 10, ' RT Rho [amu/cm^3]'
           write(iUnit, '(I5,a)') 11, ' RT P [Pa]'
! add joule heating output --Yiqun Sep 2008
           write(iUnit, '(I5,a)') 12, ' JouleHeat [mW/m2]'
           write(iUnit, '(I5,a)') 13, ' IonNumFlux [/cm2/s]'

        end select

        write(iUnit, *) ' '

        write(iUnit, '(a)') 'TIME'
        write(iUnit, '(I5,a)')  Time_Array(1), ' Year'
        write(iUnit, '(I5,a)')  Time_Array(2), ' Month'
        write(iUnit, '(I5,a)')  Time_Array(3), ' Day'
        write(iUnit, '(I5,a)')  Time_Array(4), ' Hour'
        write(iUnit, '(I5,a)')  Time_Array(5), ' Minute'
        write(iUnit, '(I5,a)')  Time_Array(6), ' Second'
        write(iUnit, '(I5,a)')  Time_Array(7), ' Millisecond'

        write(iUnit, *) ' '

        write(iUnit, '(a)') 'SIMULATION'
        write(iUnit, '(I13,a)')     nSolve,          ' nSolve'
        write(iUnit, '(1pe13.5,a)') time_simulation, ' Time_Simulation'

        write(iUnit, *) ' '

        write(iUnit, '(a)') 'DIPOLE TILT'
        write(iUnit, '(F7.2,a)')  ThetaTilt*180.0/cPi, ' ThetaTiltDeg'
        write(iUnit, '(F7.2,a)')           0.0, ' PhiTiltDeg'

        write(iUnit, *) ' '
        write(iUnit, '(a)') 'BEGIN NORTHERN HEMISPHERE'

     elseif (output_type == tec_type) then

        write(iUnit, '(a,i4.4,"-",5(i2.2,"-"),i3.3,a,2f7.2,a)') &
             'TITLE="BATSRUS: Ionospheric Potential Solution,  ', &
             Time_Array(1:7),', BTiltDeg=',ThetaTilt*180.0/cPi,0.,'"'

        if (variables == min_vars) then
           write(iUnit, *)  'VARIABLES= "Theta [deg]","Psi [deg]"'
           write(iUnit, *)  ' "SigmaH [S]","SigmaP [S]"'
           write(iUnit, *)  ' "JR [`mA/m^2]","PHI [kV]"'

        elseif (variables == all_vars) then
           write(iUnit, *)  'VARIABLES= "X [R]","Y [R]","Z [R]"'
           write(iUnit, *)  ' "Theta [deg]","Psi [deg]"'
           write(iUnit, *)  ' "SigmaH [S]","SigmaP [S]"'
           write(iUnit, *)  ' "E-Flux [W/m2]"'
           write(iUnit, *)  ' "Ave-E [eV]"'
           write(iUnit, *)  ' "JR [`mA/m^2]","PHI [kV]"'
           write(iUnit, *)  ' "Ex [mV/m]","Ey [mV/m]","Ez [mV/m]"'
           write(iUnit, *)  ' "Jx [`mA/m^2]","Jy [`mA/m^2]","Jz [`mA/m^2]"'
           write(iUnit, *)  ' "Ux [km/s]","Uy [km/s]","Uz [km/s]"'
! add joule heating output --Yiqun Sep 2008
           write(iUnit, *)  ' "JouleHeat [mW/m2]"'
           write(iUnit, *)  ' "IonNumFlux [/cm2/s]"'

        elseif (variables == uam_vars) then             !^CFG  IF TIEGCM BEGIN
           write(iUnit, *)  'VARIABLES= "Theta [deg]","Psi [deg]"'
           write(iUnit, *)  ' "SigmaH [S]","SigmaP [S]"'
           write(iUnit, *)  ' "JR [`mA/m^2]","JR (NW) [`mA/m^2]",'
           write(iUnit, *)  ' "PHI [kV]"'               !^CFG END TIEGCM

        elseif (variables == aur_vars) then
           write(iUnit, *)  'VARIABLES= "Theta [deg]","Psi [deg]"'
           write(iUnit, *)  ' "SigmaH [S]","SigmaP [S]"'
           write(iUnit, *)  ' "JR [`mA/m^2]","PHI [kV]"'
           write(iUnit, *)  ' "E-Flux [W/m2]"'
           write(iUnit, *)  ' "Ave-E [eV]"'
           write(iUnit, *)  ' "RT 1/B [1/T]","RT Rho [amu/cm^3]","RT P [Pa]"'
! add joule heating output --Yiqun Sep 2008
           write(iUnit, *)  ' "JouleHeat [mW/m2]"'
           write(iUnit, *)  ' "IonNumFlux [/cm2/s]"'

        endif

        write(iUnit,'(a)') 'ZONE T="IonN '//textNandT//'"'
        write(iUnit, *)  'I= ',IONO_nTheta,' J= ',IONO_nPsi,' F=POINT'

     endif

     if (variables == min_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(6(E13.5))")  &
                   cRadToDeg*IONO_NORTH_Theta(i,j), &
                   cRadToDeg*IONO_NORTH_Psi(i,j), &
                   IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                   1.0e06*IONO_NORTH_JR(i,j),   &
                   1.0e-03*IONO_NORTH_PHI(i,j)
           end do
        end do

     elseif (variables == all_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(22(E13.5))")  &
                   IONO_NORTH_X(i,j),IONO_NORTH_Y(i,j),IONO_NORTH_Z(i,j), &
                   cRadToDeg*IONO_NORTH_Theta(i,j), &
                   cRadToDeg*IONO_NORTH_Psi(i,j), &
                   IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                   IONO_NORTH_EFlux(i,j), &
                   IONO_NORTH_Ave_E(i,j), &
                   1.0e06*IONO_NORTH_JR(i,j),1.0e-03*IONO_NORTH_PHI(i,j), &
                   1.0e03*IONO_NORTH_Ex(i,j),1.0e03*IONO_NORTH_Ey(i,j), &
                   1.0e03*IONO_NORTH_Ez(i,j), &
                   1.0e06*IONO_NORTH_Jx(i,j),1.0e06*IONO_NORTH_Jy(i,j), &
                   1.0e06*IONO_NORTH_Jz(i,j), &
                   1.0e-03*IONO_NORTH_Ux(i,j),1.0e-03*IONO_NORTH_Uy(i,j), &
                   1.0e-03*IONO_NORTH_Uz(i,j), &
! add joule heating Yiqun Sep 2008
                   1.0e03*IONO_NORTH_Joule(i,j), &
                   1.0e-04*IONO_NORTH_IonNumFlux(i,j)
           end do
        end do

     elseif (variables == uam_vars) then              !^CFG  IF TIEGCM BEGIN
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(6(E13.5))")  &
                   cRadToDeg*IONO_NORTH_Theta(i,j), &
                   cRadToDeg*IONO_NORTH_Psi(i,j), &
                   IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                   1.0e06*IONO_NORTH_JR(i,j), &
                   1.0e06*IONO_NORTH_TGCM_JR(i,j), &
                   IONO_NORTH_EFlux(i,j), &
                   IONO_NORTH_Ave_E(i,j), &
                   1.0e-03*IONO_NORTH_PHI(i,j)
           end do
        end do                                        !^CFG END TIEGCM

     elseif (variables == aur_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(13(E13.5))")  &
                   cRadToDeg*IONO_NORTH_Theta(i,j), &
                   cRadToDeg*IONO_NORTH_Psi(i,j), &
                   IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
                   1.0e06*IONO_NORTH_JR(i,j),   &
                   1.0e-03*IONO_NORTH_PHI(i,j), &
                   IONO_NORTH_EFlux(i,j), &
                   IONO_NORTH_Ave_E(i,j), &
                   IONO_NORTH_invB(i,j),IONO_NORTH_rho(i,j),IONO_NORTH_p(i,j), &
! add joule heating Yiqun Sep 2008
                   1.0e03*IONO_NORTH_Joule(i,j), &
                   1.0e-04*IONO_NORTH_IonNumFlux(i,j)
           end do
        end do
     endif

  case(2) ! South writes data only

     if (output_type == idl_type) then
        write(iUnit, *) ' '
        write(iUnit, '(a)') 'BEGIN SOUTHERN HEMISPHERE'
     elseif (output_type == tec_type) then

        write(iUnit,'(a)') 'ZONE T="IonS '//textNandT//'"'
        write(iUnit, *) 'I= ',IONO_nTheta,' J= ',IONO_nPsi,' F=POINT'
     endif

     if (variables == min_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(6(E13.5))")  &
                   cRadToDeg*IONO_SOUTH_Theta(i,j), &
                   cRadToDeg*IONO_SOUTH_Psi(i,j), &
                   IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                   1.0e06*IONO_SOUTH_JR(i,j),   &
                   1.0e-03*IONO_SOUTH_PHI(i,j)
           end do
        end do

     elseif (variables == all_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(22(E13.5))")  &
                   IONO_SOUTH_X(i,j),IONO_SOUTH_Y(i,j),IONO_SOUTH_Z(i,j), &
                   cRadToDeg*IONO_SOUTH_Theta(i,j), &
                   cRadToDeg*IONO_SOUTH_Psi(i,j), &
                   IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                   IONO_SOUTH_EFlux(i,j), &
                   IONO_SOUTH_Ave_E(i,j), &
                   1.0e06*IONO_SOUTH_JR(i,j),1.0e-03*IONO_SOUTH_PHI(i,j), &
                   1.0e03*IONO_SOUTH_Ex(i,j),1.0e03*IONO_SOUTH_Ey(i,j), &
                   1.0e03*IONO_SOUTH_Ez(i,j), &
                   1.0e06*IONO_SOUTH_Jx(i,j),1.0e06*IONO_SOUTH_Jy(i,j), &
                   1.0e06*IONO_SOUTH_Jz(i,j), &
                   1.0e-03*IONO_SOUTH_Ux(i,j),1.0e-03*IONO_SOUTH_Uy(i,j), &
                   1.0e-03*IONO_SOUTH_Uz(i,j), &
! add joule heating Yiqun Sep 2008
                   1.0e03*IONO_SOUTH_Joule(i,j), &
                   1.0e-04*IONO_SOUTH_IonNumFlux(i,j)
           end do
        end do

     elseif (variables == uam_vars) then        !^CFG  IF TIEGCM BEGIN
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(6(E13.5))")  &
                   cRadToDeg*IONO_SOUTH_Theta(i,j), &
                   cRadToDeg*IONO_SOUTH_Psi(i,j), &
                   IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                   1.0e06*IONO_SOUTH_JR(i,j), &
                   1.0e06*IONO_SOUTH_TGCM_JR(i,j),  &
                   IONO_SOUTH_EFlux(i,j), &
                   IONO_SOUTH_Ave_E(i,j), &
                   1.0e-03*IONO_SOUTH_PHI(i,j)
           end do
        end do                                           !^CFG END TIEGCM

     elseif (variables == aur_vars) then
        do j = 1, IONO_nPsi
           do i = 1, IONO_nTheta
              write(iUnit,fmt="(13(E13.5))")  &
                   cRadToDeg*IONO_SOUTH_Theta(i,j), &
                   cRadToDeg*IONO_SOUTH_Psi(i,j), &
                   IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
                   1.0e06*IONO_SOUTH_JR(i,j),   &
                   1.0e-03*IONO_SOUTH_PHI(i,j), &
                   IONO_SOUTH_EFlux(i,j), &
                   IONO_SOUTH_Ave_E(i,j), &
                   IONO_SOUTH_invB(i,j),IONO_SOUTH_rho(i,j),IONO_SOUTH_p(i,j), &
! add joule heating Yiqun Sep 2008
                   1.0e03*IONO_SOUTH_Joule(i,j), &
                   1.0e-04*IONO_SOUTH_IonNumFlux(i,j)
           end do
        end do
     endif
  case default
     call CON_stop(NameSub//' invalid iBlock value')
  end select
  close(iUnit)

end subroutine ionosphere_write_output

!=============================================================================

subroutine IE_save_logfile

  use ModIonosphere
  use IE_ModIo
  use IE_ModMain
  use ModProcIE
  use ModMpi
  
  implicit none
  integer :: iError
  logical :: IsFirstTime = .true.
  !--------------------------------------------------------------------------
  if (.not.DoSaveLogfile) return

  if(nProc>1) call MPI_bcast(cpcp_south,1,MPI_REAL,1,iComm,iError)

  if(iProc/=0) return

  write(NameFile,'(a)') trim(NameIonoDir)//"logIE.log"
  if(IsFirstTime) then
    IsFirstTime = .false.
    open(unit=iUnit,file=NameFile,status="replace")
    write(iUnit,fmt="(a)")  'Ridley Ionosphere Model, [deg] and [kV]'
    write(iUnit,fmt="(a)") &
         'nsolve t yy mm dd hh mm ss ms ttilt ptilt cpcpn cpcps'
  else
    open(unit=iUnit,file=NameFile,status="old",position="append")
  end if
  write(iUnit,fmt="(i8,es13.5,i5,5i3,i4,f8.2,2es13.5)") &
       nSolve, Time_Simulation, Time_Array(1:7), &
       ThetaTilt*cRadToDeg, 0.0, cpcp_north, cpcp_south
  close(iUnit)

end subroutine IE_save_logfile

!-------------------------------------------------------------------------
! ionosphere_read_restart_file
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_read_restart_file(iter)
  !\
  ! This routine creates an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  use IE_ModIo
  implicit none

  integer, intent(in) :: iter

  integer :: i, j, i2, j2, nhemi, ierror
  real :: xx,yy
  character (len=4), Parameter :: IO_ext=".rst"

  write(*,*) '=> Reading restart file for ionosphere.'

  write(NameFile,'(a,i6.6,a)')trim(NameIonoDir)//"ionosphere_n",iter,IO_ext

  open(unit=iUnit,file=NameFile,status="old",iostat=ierror)

  if (ierror==0) then

     read(iUnit,*) IONO_Radius, IONO_Height, &
          IONO_Bdp, IONO_Radius_Mag_Boundary, &
          IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(iUnit,*) i2,j2,xx,yy,IONO_NORTH_JR(i,j),                  &
                SAVE_NORTH_SigmaH(i,j),SAVE_NORTH_SigmaP(i,j)
        end do
     end do

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(iUnit,*) i2,j2,xx,yy,IONO_SOUTH_JR(i,j),                  &
                SAVE_SOUTH_SigmaH(i,j),SAVE_SOUTH_SigmaP(i,j) 
        end do
     end do

     close(iUnit)

  else

     write(6,*) "==================================================="
     write(6,*) "= Error in finding ionosphere restart file.       ="
     write(6,*) "==================================================="

  endif

end subroutine ionosphere_read_restart_file

!-------------------------------------------------------------------------
! ionosphere_write_restart_file
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_write_restart_file(iter)
  !\
  ! This routine reads in an ionospheric 
  ! restart solution file.
  !/
  use ModIonosphere
  use IE_ModIo
  implicit none

  integer, intent(in) :: iter

  integer :: i, j
  character (len=4), Parameter :: IO_ext=".rst"

  write(*,*) '=> Writing restart file for ionosphere.'

  write(NameFile,'(a,i6.6,a)')trim(NameIonoDir)//"ionosphere_n",iter,IO_ext

  open(unit=iUnit,file=NameFile,status="unknown")

  write(iUnit,*) IONO_Radius, IONO_Height,                                    &
                 IONO_Bdp, IONO_Radius_Mag_Boundary,                          &
                 IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(iUnit,*) i,j,                                                   &
                       IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j),             &
                       IONO_NORTH_JR(i,j),SAVE_NORTH_SigmaH(i,j),             &
                       SAVE_NORTH_SigmaP(i,j)
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(iUnit,*) i,j,                                                   &
                       IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j),             &
                       IONO_SOUTH_JR(i,j),SAVE_SOUTH_SigmaH(i,j),             &
                       SAVE_SOUTH_SigmaP(i,j)
     end do
  end do

  close(iUnit)

end subroutine ionosphere_write_restart_file

 
subroutine iono_getpot(isize,jsize,MHD_lat,MHD_lon,MHD_pot,MHD_Jr)
  use ModIonosphere
  implicit none

  integer, intent(in) :: isize,jsize
  real, intent(in), dimension(isize) :: MHD_lat
  real, intent(in), dimension(jsize) :: MHD_lon
  real, intent(out), dimension(isize,jsize) :: MHD_pot,MHD_Jr

  integer :: i,j, ii,jj
  integer                  :: inter_lat, inter_lon
  real                     :: colat, dx_lat, dx_lon, mlt
  logical                  :: done
  logical                  :: DoTest

  do i=1,isize
     do j=1,jsize

        DoTest=.false.!!! i==10.and.j==5 !!!

        colat = (90.0 - MHD_lat(i))*(cPi/180.)
        mlt   = (       MHD_lon(j))*(cPi/180.)

        if(DoTest)write(*,*)'iono_getpot: i,j,colat,mlt=',i,j,colat,mlt

        do inter_lat=IONO_nTheta,1,-1
           if (iono_north_theta(inter_lat,1) <= colat) EXIT
        enddo

        if(inter_lat == 0)then
           write(*,*)'WARNING: IM colat is below theta(1)=',&
                colat,iono_north_theta(1,1)
           inter_lat = 1
           dx_lat = 0.0
        elseif(inter_lat == IONO_nTheta) then
           write(*,*)'WARNING: IM colat is above theta(n)=',&
                colat,iono_north_theta(IONO_nTheta,1)
           inter_lat = IONO_nTheta - 1
           dx_lat = 1.0
        else
           dx_lat = (colat - iono_north_theta(inter_lat,1))/           &
                (iono_north_theta(inter_lat+1,1)-iono_north_theta(inter_lat,1))
        endif

        if(DoTest)write(*,*)'inter_lat,theta1,theta2,dx_lat=',&
             inter_lat,iono_north_theta(inter_lat:inter_lat+1,1),dx_lat

        do inter_lon=IONO_nPsi,1,-1
           if (iono_north_psi(1,inter_lon) <= mlt) EXIT
        end do
        
        if(inter_lon==0)then
           write(*,*)'WARNING: IM longitude is below psi(1)=',&
                mlt,iono_north_psi(1,1)
           inter_lon=1
           dx_lon=0.0
        elseif(inter_lon==IONO_nPsi)then
           write(*,*)'WARNING: IM longitude is above psi(n)=',&
                mlt,iono_north_psi(1,IONO_nPsi)
           inter_lon=IONO_nPsi-1
           dx_lon=1.0
        else
           dx_lon = (mlt - iono_north_psi(1,inter_lon)) /                  &
                (iono_north_psi(1,inter_lon+1)-iono_north_psi(1,inter_lon))
        end if

        if(DoTest)write(*,*)'inter_lon,psi1,psi2,dx_lon=',&
             inter_lon,iono_north_psi(1,inter_lon:inter_lon+1),dx_lon

        MHD_pot(i,j) = &
             IONO_NORTH_PHI(inter_lat+1,inter_lon+1)*(    dx_lat)*(    dx_lon) + &
             IONO_NORTH_PHI(inter_lat  ,inter_lon+1)*(1.0-dx_lat)*(    dx_lon) + &
             IONO_NORTH_PHI(inter_lat+1,inter_lon  )*(    dx_lat)*(1.0-dx_lon) + &
             IONO_NORTH_PHI(inter_lat  ,inter_lon  )*(1.0-dx_lat)*(1.0-dx_lon)

        MHD_Jr(i,j) = &
             IONO_NORTH_JR(inter_lat+1,inter_lon+1)*(    dx_lat)*(    dx_lon) + &
             IONO_NORTH_JR(inter_lat  ,inter_lon+1)*(1.0-dx_lat)*(    dx_lon) + &
             IONO_NORTH_JR(inter_lat+1,inter_lon  )*(    dx_lat)*(1.0-dx_lon) + &
             IONO_NORTH_JR(inter_lat  ,inter_lon  )*(1.0-dx_lat)*(1.0-dx_lon)

     end do
  end do

end subroutine iono_getpot
