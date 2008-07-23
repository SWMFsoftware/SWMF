subroutine epencalc(utime,f107,bc_choice,IMF_By, IMF_Bz, SW_V)

  use ModIonoHeidi

!  use ModIonosphere
!  use ModMain
!  use ModIO
  implicit none

  real, intent(in)    :: utime,f107

  integer, intent(in) :: bc_choice
  real, intent(in)    :: IMF_By, IMF_Bz, SW_V

  integer :: iError

  logical :: debug

  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  t, m

  real*8 :: TimeOfDay

  character (len=100), dimension(100) :: Lines
  character (len=100) :: cAMIEFileNorth
  character (len=100) :: cAMIEFileSouth
  logical :: IsFirstTime = .true.

  debug = .false.
  
  plot_vars(1) = 'minimum'
  plot_form(1) = 'idl'

  IONO_NORTH_JR = 0.0
  IONO_SOUTH_JR = 0.0

!  if (debug) write(*,*) "increment_real_world_time"
!  call increment_real_world_time(utime)

!
! Get FAC Data
!

  if (debug) write(*,*) "read_ring_current"

  call read_ring_current

!!!  IONO_NORTH_JR = IONO_NORTH_JR + IONO_NORTH_RCM_JR
!!!  IONO_SOUTH_JR = IONO_SOUTH_JR + IONO_SOUTH_RCM_JR
!!!
!!!  ! If bc_choice < 0  then use weimer/AMIE over the whole thing
!!!  ! If bc_choice > 0 then use weimer/AMIE as a high lat boundary
!!!  ! If bc_choice == 0 then use no high latitude boundary
!!!
!!!  if (bc_choice /= 0) then
!!!
!!!     IONO_NORTH_PHI = 0.0
!!!     IONO_SOUTH_PHI = 0.0
!!!
!!!     if (IsFirstTime) then
!!!
!!!        if (abs(bc_choice) == 1) then
!!!
!!!           Lines(1) = "#BACKGROUND"
!!!           Lines(2) = "../IM/HEIDI/src/"
!!!           Lines(3) = "weimer96"
!!!           Lines(4) = "ihp"
!!!           Lines(5) = "idontknow"
!!!           Lines(6) = ""
!!!           Lines(7) = "#DEBUG"
!!!           Lines(8) = "2"
!!!           Lines(9) = "0"
!!!           Lines(10) = ""
!!!           Lines(11) = "#END"
!!!
!!!        endif
!!!
!!!        if (abs(bc_choice) == 2) then
!!!
!!!           cAMIEFileNorth = "b20020418.wmsf"
!!!           cAMIEFileSouth = "b20020418.wmsf"
!!!
!!!           Lines(1) = "#AMIEFILES"
!!!           Lines(2) = cAMIEFileNorth
!!!           Lines(3) = cAMIEFileSouth
!!!           Lines(4) = ""
!!!           Lines(5) = ""
!!!           Lines(6) = ""
!!!           Lines(7) = "#DEBUG"
!!!           Lines(8) = "2"
!!!           Lines(9) = "0"
!!!           Lines(10) = ""
!!!           Lines(11) = "#END"
!!!
!!!        endif
!!!
!!!        call EIE_set_inputs(Lines)
!!!        
!!!        call EIE_Initialize(iError)
!!!        
!!!        if (iError /=0) then
!!!           write(*,*) "Error in IE_Initialize, called from get_potential"
!!!           return
!!!        endif
!!!
!!!        IsFirstTime = .false.
!!!
!!!     endif
!!!
!!!     write(*,*) "time : ", Real_Time_of_Year, utime
!!!
!!!     call IO_SetTime(Real_Time_of_Year)
!!!
!!!
!!!     TimeOfDay = mod(utime,24.0*3600.0)
!!!     write(*,*) "Time to AMIE : ",TimeOfDay
!!!     if (abs(bc_choice) == 2) call get_AMIE_values(TimeOfDay)
!!!
!!!     call IO_SetIMFBz(imf_bz)
!!!     call IO_SetIMFBy(imf_by)
!!!     call IO_SetSWV  (abs(sw_v)/1000.0)
!!!     
!!!     ! I think that the IE library assumes that all arrays are 
!!!     ! nLong, nLat.  The ionosphere is nLat, nLong, so we need to
!!!     ! reverse these below...
!!!
!!!     call IO_SetnMLTs(IONO_nTheta)
!!!     call IO_SetnLats(IONO_nPsi)
!!!
!!!     t = 90.0 - IONO_NORTH_Theta*180.0/IONO_PI
!!!     m = mod(IONO_NORTH_Psi*12.0/IONO_PI + 12.0, 24.0)
!!!
!!!     call IO_SetGrid(m, t, iError)
!!!     call IO_GetPotential(IONO_NORTH_PHI, iError)
!!!
!!!write(*,*) minval(iono_north_phi),maxval(iono_north_phi)
!!!
!!!     t = IONO_SOUTH_Theta*180.0/IONO_PI - 90.0
!!!     m = mod(IONO_SOUTH_Psi*12.0/IONO_PI + 12.0, 24.0)
!!!
!!!     call IO_SetGrid(m, t, iError)
!!!     call IO_GetPotential(IONO_SOUTH_PHI, iError)
!!!
!!!     if (bc_choice < 0) then
!!!
!!!        call write_ring_current
!!!        return
!!!
!!!     endif
!!!
!!!  else
!!!
!!!     IONO_NORTH_PHI = 0.0
!!!     IONO_SOUTH_PHI = 0.0
!!!
!!!  endif
!!!
!!!  ionosphere_type = 3
!!!  f107_flux = f107
!!!
!!!  PolarCapPedConductance = 0.5
!!!  StarLightPedConductance = 1.0  ! Aaron's default value = 0.5
!!!
!!!  if (debug) write(*,*) "facs_to_fluxes"
!!!
!!!  call FACs_to_fluxes(ionosphere_type)
!!!
!!!  if (debug) write(*,*) "ionosphere_conductance"
!!!
!!!  call ionosphere_conductance(IONO_NORTH_Sigma0,               &
!!!       IONO_NORTH_SigmaH, IONO_NORTH_SigmaP,    &
!!!       IONO_NORTH_SigmaThTh,                    &
!!!       IONO_NORTH_SigmaThPs,                    &
!!!       IONO_NORTH_SigmaPsPs,                    &
!!!       IONO_NORTH_dSigmaThTh_dTheta,            &
!!!       IONO_NORTH_dSigmaThPs_dTheta,            &
!!!       IONO_NORTH_dSigmaPsPs_dTheta,            &
!!!       IONO_NORTH_dSigmaThTh_dPsi,              &
!!!       IONO_NORTH_dSigmaThPs_dPsi,              &
!!!       IONO_NORTH_dSigmaPsPs_dPsi,              &
!!!       IONO_NORTH_EFlux, IONO_NORTH_Ave_E,      &
!!!       IONO_NORTH_Theta, IONO_NORTH_Psi,        &
!!!       IONO_nTheta, IONO_nPsi,                  &
!!!       dTheta_North, dPsi_North,                &
!!!       ionosphere_type, f107_flux)
!!!        
!!!  if (debug) write(*,*) "ionosphere_solver"
!!!
!!!  if (.not.debug) then
!!!
!!!     call ionosphere_solver(IONO_NORTH_PHI,  IONO_NORTH_JR, &
!!!          IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs,   &
!!!          IONO_NORTH_SigmaPsPs,                         &
!!!          IONO_NORTH_dSigmaThTh_dTheta,                 &
!!!          IONO_NORTH_dSigmaThPs_dTheta, &
!!!          IONO_NORTH_dSigmaPsPs_dTheta, &
!!!          IONO_NORTH_dSigmaThTh_dPsi,  &
!!!          IONO_NORTH_dSigmaThPs_dPsi, &
!!!          IONO_NORTH_dSigmaPsPs_dPsi, &
!!!          IONO_NORTH_Theta, IONO_NORTH_Psi, &
!!!          IONO_Radius, IONO_nTheta, IONO_nPsi,               &
!!!          dTheta_North, dPsi_North,                     &
!!!          4, 0)
!!!
!!!!     call ionosphere_currents(IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz,&
!!!!          IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
!!!!          IONO_NORTH_ETh,IONO_NORTH_EPs, &
!!!!          IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
!!!!          IONO_NORTH_PHI, &
!!!!          IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
!!!!          IONO_NORTH_SigmaPsPs, &
!!!!          IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
!!!!          IONO_NORTH_Theta, IONO_NORTH_Psi, &
!!!!          IONO_Radius, IONO_nTheta, IONO_nPsi,              &
!!!!          dTheta_North, dPsi_North)
!!!
!!!  endif
!!!
!!!  if (debug) write(*,*) "ionosphere_conducatance"
!!!
!!!  call ionosphere_conductance(IONO_SOUTH_Sigma0,               &
!!!       IONO_SOUTH_SigmaH, &
!!!       IONO_SOUTH_SigmaP, &
!!!       IONO_SOUTH_SigmaThTh, &
!!!       IONO_SOUTH_SigmaThPs, &
!!!       IONO_SOUTH_SigmaPsPs, &
!!!       IONO_SOUTH_dSigmaThTh_dTheta, &
!!!       IONO_SOUTH_dSigmaThPs_dTheta, &
!!!       IONO_SOUTH_dSigmaPsPs_dTheta, &
!!!       IONO_SOUTH_dSigmaThTh_dPsi, &
!!!       IONO_SOUTH_dSigmaThPs_dPsi, &
!!!       IONO_SOUTH_dSigmaPsPs_dPsi, &
!!!       IONO_SOUTH_EFlux, IONO_SOUTH_Ave_E,  &
!!!       IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
!!!       IONO_nTheta, IONO_nPsi,                  &
!!!       dTheta_South, dPsi_South,                &
!!!       ionosphere_type, f107_flux)
!!!        
!!!  if (debug) write(*,*) "ionosphere_solver"
!!!
!!!  if (.not.debug) then
!!!
!!!     call ionosphere_solver(IONO_SOUTH_PHI, IONO_South_JR, &
!!!          IONO_SOUTH_SigmaThTh, &
!!!          IONO_SOUTH_SigmaThPs, &
!!!          IONO_SOUTH_SigmaPsPs, &
!!!          IONO_SOUTH_dSigmaThTh_dTheta, &
!!!          IONO_SOUTH_dSigmaThPs_dTheta, &
!!!          IONO_SOUTH_dSigmaPsPs_dTheta, &
!!!          IONO_SOUTH_dSigmaThTh_dPsi, &
!!!          IONO_SOUTH_dSigmaThPs_dPsi, &
!!!          IONO_SOUTH_dSigmaPsPs_dPsi, &
!!!          IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
!!!          IONO_Radius, IONO_nTheta, IONO_nPsi,               &
!!!          dTheta_South, dPsi_South,                     &
!!!          4, 1)
!!!
!!!!        call ionosphere_currents(IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz,&
!!!!             IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
!!!!             IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
!!!!             IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
!!!!             IONO_SOUTH_PHI, &
!!!!             IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
!!!!             IONO_SOUTH_SigmaPsPs, &
!!!!             IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
!!!!             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
!!!!             IONO_Radius, IONO_nTheta, IONO_nPsi,              &
!!!!             dTheta_South, dPsi_South)
!!!
!!!  endif

  if (debug) write(*,*) "write_ring_current"

  call write_ring_current

end subroutine epencalc
