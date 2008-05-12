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
  use ModMain
  implicit none

  include "../TIEGCM/mhd_f90.h"

  integer, intent(in) :: iter,iAction

  integer :: iModel, input_stat, iunit_mhd,i,j
  real :: f107

  real :: Radius
  real :: New_seconds

  if (ionosphere_type.lt.0) then
     iModel = abs(ionosphere_type)
  else
     iModel = ionosphere_type
  endif

  f107 = f107_flux

  select case (iAction)

     !\
     ! Create fine grids for north and south 
     ! hemisphere ionospheric solutions
     !/

     case (1)
        IONO_PI = 2.00*asin(1.00)
        call ionosphere_fine_grid
        call ionosphere_init

        IONO_Old_Seconds = 3600.00

        if (restart.and.time_accurate) then
           New_Seconds = Start_Time_Array(6)+Start_Time_Array(7)/1000.0
           IONO_Old_Seconds = New_Seconds
        endif

     case(2)

        if (IONO_NORTH_nMagBndPts >= 32 .and. &
            IONO_SOUTH_nMagBndPts >= 32) call ionosphere_fac

     case(3)

        call ionosphere_read_restart_file(iter)

     !\
     ! Write ionospheric solution for north and south 
     ! hemispheres to output data file.
     !/
     case (4)
        call ionosphere_write_output(iter)

     !\
     ! Create a restart solution file containing the
     ! FAC driving the ionospheric solution.
     !/
     case (5)
        call ionosphere_write_restart_file(iter)

     !\
     ! Perform ionophere calculations for north and south 
     ! hemispheres using magnetospheric FAC solution.
     ! A multigrid technique is used to solve the elliptic
     ! equations for the ionospheric potential on each of
     ! the hemispheres.
     !/
     case (6)

!        if (IONO_NORTH_nMagBndPts >= 32 .and. &
!            IONO_SOUTH_nMagBndPts >= 32) then

           Radius = IONO_Radius + IONO_Height

           write(6,*) "=> Calling Ionosphere at Real World Time :"
           write(6,"(7I5)") Start_Time_Array

           New_Seconds = Start_Time_Array(6)+Start_Time_Array(7)/1000.0

           input_stat = 0
           if (New_Seconds <= IONO_Old_Seconds) then
              if (TIMEGCM_Coupled) then
                 write(6,*) "=> Coupling with the TIMEGCM at this time"
                 write(6,*) "=> New Seconds : ",New_Seconds
                 write(6,*) "=> Old Seconds : ",IONO_Old_Seconds
                 input_stat = 0
                 iunit_mhd = 23
                 open(unit = iunit_mhd, file="./ionosphere/MHD_to_TIMEGCM.dat", &
                      status="OLD",iostat=input_stat)
                 if (input_stat.ne.0) then
                    write(6,*) "Skipping TIEGCM this time, since there is no"
                    write(6,*) "MHD potential file. Will run it next time"
                 else 
                    close(iunit_mhd)
                    call tgcm(first_time,mhd_flag)
                 endif
              endif
              if (RCM_Coupled) then
                 write(6,*) "=> Coupling with the RCM at this time"
                 !              call rcm
              endif

           endif

           if (input_stat.eq.0) IONO_Old_Seconds = New_Seconds

           if (TIMEGCM_Coupled) then
              call read_timegcm_file
           endif
           if (RCM_Coupled) then
              call read_rcm_file
           endif

           call FACs_to_fluxes

           open(unit=13,file='facs.txt',status='unknown')
           write(13,fmt="(2I10)") IONO_nTheta/2, IONO_nPsi
           do i = 1, IONO_nTheta/2
              do j = 1, IONO_nPsi
                 write(13,fmt="(1P,6(E13.5))")  &
                      (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
                      (360.00/(2.00*IONO_PI))*IONO_NORTH_PSI(i,j), &
                      iono_north_jr(i,j),iono_north_tgcm_jr(i,j), &
                      iono_north_eflux(i,j), iono_north_ave_e(i,j)
              end do
           enddo
           close(13)

           if (TIMEGCM_Coupled)                                               &
                IONO_NORTH_JR = IONO_NORTH_JR - IONO_NORTH_TGCM_JR

           if (RCM_Coupled)                                                   &
                IONO_NORTH_JR = IONO_NORTH_JR + IONO_NORTH_RCM_JR

           IONO_NORTH_PHI = IONO_NORTH_JR

           call ionosphere_conductance(IONO_NORTH_Sigma0,               &
                                    IONO_NORTH_SigmaH, IONO_NORTH_SigmaP,    &
                                    IONO_NORTH_SigmaThTh,                    &
                                    IONO_NORTH_SigmaThPs,                    &
                                    IONO_NORTH_SigmaPsPs,                    &
                                    IONO_NORTH_dSigmaThTh_dTheta,            &
                                    IONO_NORTH_dSigmaThPs_dTheta,            &
                                    IONO_NORTH_dSigmaPsPs_dTheta,            &
                                    IONO_NORTH_dSigmaThTh_dPsi,              &
                                    IONO_NORTH_dSigmaThPs_dPsi,              &
                                    IONO_NORTH_dSigmaPsPs_dPsi,              &
                                    IONO_NORTH_EFlux, IONO_NORTH_Ave_E,      &
                                    IONO_NORTH_Theta, IONO_NORTH_Psi,        &
                                    IONO_nTheta, IONO_nPsi,                  &
                                    dTheta_North, dPsi_North,                &
                                    iModel, f107)

           call ionosphere_solver(IONO_NORTH_PHI,                               &
                               IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs,   &
                               IONO_NORTH_SigmaPsPs,                         &
                               IONO_NORTH_dSigmaThTh_dTheta,                 &
                               IONO_NORTH_dSigmaThPs_dTheta, &
                               IONO_NORTH_dSigmaPsPs_dTheta, &
                               IONO_NORTH_dSigmaThTh_dPsi,  &
                               IONO_NORTH_dSigmaThPs_dPsi, &
                               IONO_NORTH_dSigmaPsPs_dPsi, &
                               IONO_NORTH_Theta, IONO_NORTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi,               &
                               dTheta_North, dPsi_North,                     &
                               4, 0)

           call ionosphere_currents(IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz, &
                                 IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
                                 IONO_NORTH_ETh,IONO_NORTH_EPs, &
                                 IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
                                 IONO_NORTH_PHI, &
                                 IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
                                 IONO_NORTH_SigmaPsPs, &
                                 IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                                 IONO_NORTH_Theta, IONO_NORTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi,              &
                                 dTheta_North, dPsi_North)

           call ionosphere_magBCs(IONO_NORTH_PHI_BC, IONO_NORTH_ETh_BC, &
                               IONO_NORTH_EPs_BC, &
                               IONO_NORTH_UR_BC, IONO_NORTH_UTh_BC, &
                               IONO_NORTH_UPs_BC, &
                               IONO_Radius_Mag_Boundary, IONO_NORTH_PHI, &
                               IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
                               IONO_NORTH_Theta, IONO_NORTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi,               &
                               dTheta_North, dPsi_North)

           if (TIMEGCM_Coupled)                                                 &
                IONO_SOUTH_JR = IONO_SOUTH_JR - IONO_SOUTH_TGCM_JR

           if (RCM_Coupled)                                                     &
                IONO_SOUTH_JR = IONO_SOUTH_JR + IONO_SOUTH_RCM_JR

           IONO_SOUTH_PHI = IONO_SOUTH_JR

           call ionosphere_conductance(IONO_SOUTH_Sigma0, IONO_SOUTH_SigmaH, &
                                    IONO_SOUTH_SigmaP, &
                                    IONO_SOUTH_SigmaThTh, &
                                    IONO_SOUTH_SigmaThPs, &
                                    IONO_SOUTH_SigmaPsPs, &
                                    IONO_SOUTH_dSigmaThTh_dTheta, &
                                    IONO_SOUTH_dSigmaThPs_dTheta, &
                                    IONO_SOUTH_dSigmaPsPs_dTheta, &
                                    IONO_SOUTH_dSigmaThTh_dPsi, &
                                    IONO_SOUTH_dSigmaThPs_dPsi, &
                                    IONO_SOUTH_dSigmaPsPs_dPsi, &
                                    IONO_SOUTH_EFlux, IONO_SOUTH_Ave_E,  &
                                    IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                                    IONO_nTheta, IONO_nPsi,                  &
                                    dTheta_South, dPsi_South,                &
                                    iModel, f107)

           call ionosphere_solver(IONO_SOUTH_PHI, &
                               IONO_SOUTH_SigmaThTh, &
                               IONO_SOUTH_SigmaThPs, &
                               IONO_SOUTH_SigmaPsPs, &
                               IONO_SOUTH_dSigmaThTh_dTheta, &
                               IONO_SOUTH_dSigmaThPs_dTheta, &
                               IONO_SOUTH_dSigmaPsPs_dTheta, &
                               IONO_SOUTH_dSigmaThTh_dPsi, &
                               IONO_SOUTH_dSigmaThPs_dPsi, &
                               IONO_SOUTH_dSigmaPsPs_dPsi, &
                               IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi,               &
                               dTheta_South, dPsi_South,                     &
                               4, 1)

           call ionosphere_currents(IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz, &
                                 IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
                                 IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
                                 IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
                                 IONO_SOUTH_PHI, &
                                 IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
                                 IONO_SOUTH_SigmaPsPs, &
                                 IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
                                 IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                                 Radius, IONO_nTheta, IONO_nPsi,              &
                                 dTheta_South, dPsi_South)

           call ionosphere_magBCs(IONO_SOUTH_PHI_BC, IONO_SOUTH_ETh_BC, &
                               IONO_SOUTH_EPs_BC, &
                               IONO_SOUTH_UR_BC, IONO_SOUTH_UTh_BC, &
                               IONO_SOUTH_UPs_BC, &
                               IONO_Radius_Mag_Boundary, IONO_SOUTH_PHI, &
                               IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
                               IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
                               Radius, IONO_nTheta, IONO_nPsi,               &
                               dTheta_South, dPsi_South)

           if (TIMEGCM_Coupled)                                                 &
                call write_timegcm_file(iter, IONO_NORTH_PHI, IONO_SOUTH_PHI,   &
                                     IONO_NORTH_EFlux, IONO_SOUTH_EFlux,     &
                                     IONO_NORTH_Ave_E, IONO_SOUTH_Ave_E)


           if (ionosphere_type.lt.0) new_conductance = .false.

!        endif

  end select  

end subroutine ionosphere

!*************************************************************************
!
! INITIALIZATION Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! ionosphere_fine_grid
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_fine_grid
  !\
  ! This routine sets the fine grid meshes for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  implicit none

  integer :: i,j
  real :: dTheta_l, dPsi_l
        
  dTheta_l = 0.50*IONO_PI/real(IONO_nTheta-1)
  dPsi_l = dTheta_l

  do j = 1, IONO_nPsi
     IONO_NORTH_Theta(1,j) = IONO_Theta_0
     IONO_NORTH_Psi(1,j) = real(j-1)*dPsi_l
     IONO_NORTH_X(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         cos(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Y(1,j) = sin(IONO_NORTH_Theta(1,j))* &
                         sin(IONO_NORTH_Psi(1,j))
     IONO_NORTH_Z(1,j) = cos(IONO_NORTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_NORTH_Theta(i,j) = IONO_PI*(0.5*real(i-1)/real(IONO_nTheta))
        IONO_NORTH_Psi(i,j) = IONO_NORTH_Psi(1,j)
        IONO_NORTH_X(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            cos(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Y(i,j) = sin(IONO_NORTH_Theta(i,j))* &
                            sin(IONO_NORTH_Psi(i,j))
        IONO_NORTH_Z(i,j) = cos(IONO_NORTH_Theta(i,j))
     end do  

     IONO_NORTH_Theta(IONO_nTheta,j) = 0.50*IONO_PI
     IONO_NORTH_Psi(IONO_nTheta,j) = IONO_NORTH_Psi(1,j)
     IONO_NORTH_X(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Y(IONO_nTheta,j) = sin(IONO_NORTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_NORTH_Psi(IONO_nTheta,j))
     IONO_NORTH_Z(IONO_nTheta,j) = cos(IONO_NORTH_Theta(IONO_nTheta,j))
  end do

  do j = 1, IONO_nPsi
     IONO_SOUTH_Theta(1,j) = IONO_PI-0.50*IONO_PI
     IONO_SOUTH_Psi(1,j) = real(j-1)*dPsi_l
     IONO_SOUTH_X(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         cos(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Y(1,j) = sin(IONO_SOUTH_Theta(1,j))* &
                         sin(IONO_SOUTH_Psi(1,j))
     IONO_SOUTH_Z(1,j) = cos(IONO_SOUTH_Theta(1,j))

     do i = 1, IONO_nTheta
        IONO_SOUTH_Theta(i,j) = IONO_PI*(0.5+0.5*real(i-1)/real(IONO_nTheta))
        IONO_SOUTH_Psi(i,j) = IONO_SOUTH_Psi(1,j)
        IONO_SOUTH_X(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            cos(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Y(i,j) = sin(IONO_SOUTH_Theta(i,j))* &
                            sin(IONO_SOUTH_Psi(i,j))
        IONO_SOUTH_Z(i,j) = cos(IONO_SOUTH_Theta(i,j))
     end do  

     IONO_SOUTH_Theta(IONO_nTheta,j) = IONO_PI-IONO_Theta_0
     IONO_SOUTH_Psi(IONO_nTheta,j) = IONO_SOUTH_Psi(1,j)
     IONO_SOUTH_X(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   cos(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Y(IONO_nTheta,j) = sin(IONO_SOUTH_Theta(IONO_nTheta,j))* &
                                   sin(IONO_SOUTH_Psi(IONO_nTheta,j))
     IONO_SOUTH_Z(IONO_nTheta,j) = cos(IONO_SOUTH_Theta(IONO_nTheta,j))
  end do

!
! dPsi is going to be defined as 2*dPsi for all of the points:
!

  do j = 2, IONO_nPsi-1
     dPsi_North(j) = IONO_NORTH_Psi(1,j+1)-IONO_NORTH_Psi(1,j-1)
  enddo
  dPsi_North(1)    = IONO_NORTH_Psi(1,2) -                                    &
       (IONO_NORTH_Psi(1,IONO_nPsi-1) - 2.0*IONO_PI)
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
       (IONO_SOUTH_Psi(1,IONO_nPsi-1)-2.0*IONO_PI)
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

!-------------------------------------------------------------------------
! ionosphere_init
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_init
  !\
  ! This routine initializes the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres.
  !/
  use ModIonosphere
  implicit none
  include "../TIEGCM/mhd_f90.h"

  mhd_flag = .true.
  first_time = .true.

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

  IONO_NORTH_PHI_BC = 0.00
  IONO_NORTH_ETh_BC = 0.00
  IONO_NORTH_EPs_BC = 0.00
  IONO_NORTH_UR_BC = 0.00 
  IONO_NORTH_UTh_BC = 0.00 
  IONO_NORTH_UPs_BC = 0.00

  IONO_SOUTH_PHI_BC = 0.00
  IONO_SOUTH_ETh_BC = 0.00
  IONO_SOUTH_EPs_BC = 0.00
  IONO_SOUTH_UR_BC = 0.00
  IONO_SOUTH_UTh_BC = 0.00 
  IONO_SOUTH_UPs_BC = 0.00

  IONO_NORTH_TGCM_JR = 0.00
  IONO_SOUTH_TGCM_JR = 0.00

  new_conductance = .true.

end subroutine ionosphere_init

!*************************************************************************
!
! INPUT/OUTPUT Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! write_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine write_timegcm_file(iter, phi_north, phi_south,   &
     eflux_north, eflux_south, avee_north, avee_south)

  use ModIonosphere
  use ModMain
  implicit none

  integer, intent(in) :: iter
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::             &
       phi_north, phi_south, eflux_north, eflux_south,      &
       avee_north, avee_south
  real psi_offset

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".dat"

  write(*,*) '=> Writing output datafiles for TIMEGCM.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"MHD_to_TIMEGCM"//IO_ext, &
       status="unknown")

!  Write Year, month, day, hour, minute, second

  write(Iunit,fmt="(6I5)") (Start_Time_Array(i),i=1,6)

! Write number of thetas and psis

  write(Iunit,fmt="(2I10)") IONO_nTheta*2, IONO_nPsi

! Silliness to occur here:
!   This ionospheric code starts with noon at 0.0 longitude.  This is 
!   different than almost all other ionospheres, so we need to shift
!   the output by 180.0 degrees.
!   Also, we want 0.0 and 360.0 to both be printed - but the new ones,
!   so we have to ignore the old 0.0 and repeat the old 180.0.

  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_south(i,j),eflux_south(i,j), &
             avee_south(i,j)
     end do
  end do
  do i = IONO_nTheta, 1, -1
     do j = IONO_nPsi/2+1, IONO_nPsi
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
     do j = 2, IONO_nPsi/2+1
        psi_offset = (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j) - 180.0
        if (psi_offset < 0) psi_offset = psi_offset + 360.0
        write(Iunit,fmt="(1P,5(E13.5))")  &
             (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
             psi_offset, &
             1.0e-03*phi_north(i,j),eflux_north(i,j), &
             avee_north(i,j)
     end do
  end do

  close(UNIT=Iunit)
 
end subroutine write_timegcm_file

!-------------------------------------------------------------------------
! read_timegcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine read_timegcm_file

  use ModIonosphere
  implicit none

  integer :: Iunit, i, j, iter, nTheta, nPsi
  character (len=13), Parameter :: iono_dir="./"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=100) :: TIMEGCM_file
  real :: dum_lat, dum_lon, min_lat, max_lat, max_fac
  integer, dimension(1:6) :: time_array
  real :: tgcm_time
  integer :: ierror

  ierror = 0

  min_lat = 55.0
  max_lat = 85.0
  max_fac = 0.0

  write(*,*) '=> Reading output file from the TIMEGCM.'

  TIMEGCM_file = "./ionosphere/TIMEGCM_to_MHD.dat"

  Iunit = 24
  open(unit=Iunit,file=TIMEGCM_file, status="old",iostat=ierror)

  if (ierror.ne.0) then
     write(6,*) "Could not open TIMEGCM output file. Skipping"
  else

     read(Iunit,fmt="(6I5)") (time_array(i),i=1,6)

     read(Iunit,fmt="(2I10)") nTheta, nPsi

     do i = IONO_nTheta,1,-1
        do j = 1, IONO_nPsi
           read(Iunit,fmt="(3(E13.5))")  &
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

     if (min_lat > 0.0) &
          write(6,*) "=> Maximum abs(FAC) below ",min_lat," (SH) is ",max_fac

     max_fac = 0.0
     do i = IONO_nTheta,1,-1
        do j = 1, IONO_nPsi
           read(Iunit,fmt="(3(E13.5))")  &
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

     if (min_lat > 0.0) &
          write(6,*) "=> Maximum abs(FAC) below ",min_lat," (NH) is ",max_fac

     close(UNIT=Iunit)
 
   endif

end subroutine read_timegcm_file

!-------------------------------------------------------------------------
! read_rcm_file
!
!
!
!-------------------------------------------------------------------------

subroutine read_rcm_file

  use ModIonosphere
  implicit none

  integer :: Iunit, i, j, iter, nTheta, nPsi
  character (len=13), Parameter :: iono_dir="../d.timegcm/"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=100) :: TIMEGCM_file
  real :: dum_lat, dum_lon, min_lat, max_lat
  integer, dimension(1:6) :: time_array
  real :: tgcm_time

  min_lat = 5.0
  max_lat = 180.0 - min_lat

  write(*,*) '=> Reading output file from the TIMEGCM.'

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"TIMEGCM_to_MHD_filemanager"//IO_ext, &
       status="old")
  read(Iunit,fmt="(A100)") TIMEGCM_file
  close(UNIT=Iunit)

  open(unit=Iunit,file=iono_dir//TIMEGCM_file, status="old")

  read(Iunit,fmt="(6I5)") (time_array(i),i=1,6)

  read(Iunit,fmt="(2I10)") nTheta, nPsi

  do i = 1, IONO_nTheta
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_NORTH_TGCM_JR(i,j)
        if (dum_lat.lt.min_lat) IONO_NORTH_TGCM_JR(i,j) = 0.0
     end do
  end do

  do i = 1, IONO_nTheta
     do j = 1, IONO_nPsi
        read(Iunit,fmt="(3(E13.5))")  &
             dum_lat, dum_lon, IONO_SOUTH_TGCM_JR(i,j)
        if (dum_lat.gt.max_lat) IONO_SOUTH_TGCM_JR(i,j) = 0.0
     end do
  end do

  close(UNIT=Iunit)
 
end subroutine read_rcm_file

!-------------------------------------------------------------------------
! ionosphere_write_output
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_write_output(iter)
  !\
  ! This routine writes out the fine grid 
  ! ionospheric solutions for the
  ! northern and southern hemispheres to
  ! a data file.
  !/
  use ModIonosphere
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".dat"
  character (len=7) :: NUMiter

  write(*,*) '=> Writing output datafiles for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_output: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//"_"//NUMiter//IO_ext, &
       status="unknown")

  write(Iunit, *)  'TITLE="BATSRUS: Ionospheric Potential Solution"'
  write(Iunit, *)  'VARIABLES= "X [R]","Y [R]","Z [R]"'
  write(Iunit, *)  ' "Theta [deg]","Psi [deg]"'
  write(Iunit, *)  ' "SigmaH [S]","SigmaP [S]"'
  write(Iunit, *)  ' "JR [`mA/m^2]","PHI [kV]"'
  write(Iunit, *)  ' "Ex [mV/m]","Ey [mV/m]","Ez [mV/m]"'
  write(Iunit, *)  ' "Jx [`mA/m^2]","Jy [`mA/m^2]","Jz [`mA/m^2]"'
  write(Iunit, *)  ' "Ux [km/s]","Uy [km/s]","Uz [km/s]"'
  write(Iunit, *)  'ZONE T="Northern Hemisphere" ','I= ',IONO_nTheta, &
                   ' J= ',IONO_nPsi,' F=POINT'
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,fmt="(18(E13.5))")  &
              IONO_NORTH_X(i,j),IONO_NORTH_Y(i,j),IONO_NORTH_Z(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_NORTH_Theta(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_NORTH_Psi(i,j), &
              IONO_NORTH_SigmaH(i,j),IONO_NORTH_SigmaP(i,j), &
              1.0e06*IONO_NORTH_JR(i,j),1.0e-03*IONO_NORTH_PHI(i,j), &
              1.0e03*IONO_NORTH_Ex(i,j),1.0e03*IONO_NORTH_Ey(i,j), &
              1.0e03*IONO_NORTH_Ez(i,j), &
              1.0e06*IONO_NORTH_Jx(i,j),1.0e06*IONO_NORTH_Jy(i,j), &
              1.0e06*IONO_NORTH_Jz(i,j), &
              1.0e-03*IONO_NORTH_Ux(i,j),1.0e-03*IONO_NORTH_Uy(i,j), &
              1.0e-03*IONO_NORTH_Uz(i,j)
     end do
  end do
  write(Iunit, *)  'ZONE T="Southern Hemisphere" ','I= ',IONO_nTheta, &
                   ' J= ',IONO_nPsi,' F=POINT'
  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,fmt="(18(E13.5))")  &
              IONO_SOUTH_X(i,j),IONO_SOUTH_Y(i,j),IONO_SOUTH_Z(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_SOUTH_Theta(i,j), &
              (360.00/(2.00*IONO_PI))*IONO_SOUTH_Psi(i,j), &
              IONO_SOUTH_SigmaH(i,j),IONO_SOUTH_SigmaP(i,j), &
              1.0e06*IONO_SOUTH_JR(i,j),1.0e-03*IONO_SOUTH_PHI(i,j), &
              1.0e03*IONO_SOUTH_Ex(i,j),1.0e03*IONO_SOUTH_Ey(i,j), &
              1.0e03*IONO_SOUTH_Ez(i,j), &
              1.0e06*IONO_SOUTH_Jx(i,j),1.0e06*IONO_SOUTH_Jy(i,j), &
              1.0e06*IONO_SOUTH_Jz(i,j), &
              1.0e-03*IONO_SOUTH_Ux(i,j),1.0e-03*IONO_SOUTH_Uy(i,j), &
              1.0e-03*IONO_SOUTH_Uz(i,j)
     end do
  end do

  close(UNIT=Iunit)

end subroutine ionosphere_write_output

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
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j, i2, j2, nhemi, ierror
  real :: xx,yy
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"
  character (len=7) :: NUMiter

  write(*,*) '=> Reading restart file for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_restart: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere"//"_"//NUMiter//IO_ext, &
       status="old",iostat=ierror)

  if (ierror.eq.0) then

     read(Iunit,*) IONO_Radius, IONO_Height, &
          IONO_Ref_Density, IONO_Ref_SoundSpeed, &
          IONO_Bdp, IONO_Radius_Mag_Boundary, &
          IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(Iunit,*) i2,j2,xx,yy,IONO_NORTH_JR(i,j),                  &
                SAVE_NORTH_SigmaH(i,j),SAVE_NORTH_SigmaP(i,j)
        end do
     end do

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta
           read(Iunit,*) i2,j2,xx,yy,IONO_SOUTH_JR(i,j),                  &
                SAVE_SOUTH_SigmaH(i,j),SAVE_SOUTH_SigmaP(i,j) 
        end do
     end do

     close(UNIT=Iunit)

     skip_ionosphere_restart = .false.
     write(6,*) skip_ionosphere_restart

  else

     write(6,*) "==================================================="
     write(6,*) "= Error in finding ionosphere restart file.       ="
     write(6,*) "==================================================="

     skip_ionosphere_restart = .true.

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
  implicit none

  integer, intent(in) :: iter

  integer :: Iunit, i, j
  character (len=13), Parameter :: iono_dir="./ionosphere/"
  character (len=4), Parameter :: IO_ext=".rst"
  character (len=7) :: NUMiter

  write(*,*) '=> Writing restart file for ionosphere.'

  if (iter < 10) then
     write (NUMiter,'(a6,i1)') 'n00000',iter
  elseif (iter < 100) then
     write (NUMiter,'(a5,i2)') 'n0000',iter
  elseif (iter < 1000) then
     write (NUMiter,'(a4,i3)') 'n000',iter
  elseif (iter < 10000) then
     write (NUMiter,'(a3,i4)') 'n00',iter
  elseif (iter < 100000) then
     write (NUMiter,'(a2,i5)') 'n0',iter
  elseif (iter < 1000000) then
     write (NUMiter,'(a1,i6)') 'n',iter
  else
     write(*,*) "ionosphere_restart: iter = ",iter, &
                " Error, too many iterations for filename"
     stop
  end if

  Iunit = 23
  open(unit=Iunit,file=iono_dir//"ionosphere_"//NUMiter//IO_ext,              &
       status="unknown")

  write(Iunit,*) IONO_Radius, IONO_Height,                                    &
                 IONO_Ref_Density, IONO_Ref_SoundSpeed,                       &
                 IONO_Bdp, IONO_Radius_Mag_Boundary,                          &
                 IONO_NORTH_Theta_Max, IONO_SOUTH_Theta_Min

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j,                                                   &
                       IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j),             &
                       IONO_NORTH_JR(i,j),SAVE_NORTH_SigmaH(i,j),             &
                       SAVE_NORTH_SigmaP(i,j)
     end do
  end do

  do j = 1, IONO_nPsi
     do i = 1, IONO_nTheta
        write(Iunit,*) i,j,                                                   &
                       IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j),             &
                       IONO_SOUTH_JR(i,j),SAVE_SOUTH_SigmaH(i,j),             &
                       SAVE_SOUTH_SigmaP(i,j)
     end do
  end do

  close(UNIT=Iunit)

end subroutine ionosphere_write_restart_file

 
