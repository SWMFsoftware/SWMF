subroutine epencalc(SimulationTime,f107,bc_choice,IMF_By, IMF_Bz, SW_V)

  use ModIonoHeidi
  use ModKind,        ONLY: Real8_
  use ModConst,       ONLY: cPi, cRadToDeg
  use ModTimeConvert, ONLY : time_int_to_real
  use ModHeidiIO,     ONLY :TimeArray

  implicit none

  real, intent(in)    :: SimulationTime,f107
  real, intent(in)    :: IMF_By, IMF_Bz, SW_V
  integer, intent(in) :: bc_choice
  integer :: iError
  logical :: debug
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  t, m
  character (len=100), dimension(100) :: Lines
  character (len=100) :: cAMIEFileNorth
  character (len=100) :: cAMIEFileSouth
  logical :: IsFirstTime = .true.
  real(Real8_) :: RealTime,TimeOfDay
  !-----------------------------------------------------------
  debug = .false.

  plot_vars(1) = 'minimum'
  plot_form(1) = 'idl'

  IONO_NORTH_JR = 0.0
  IONO_SOUTH_JR = 0.0

  if (debug) write(*,*) "increment_real_world_time"
  call time_int_to_real(TimeArray,RealTime)
  !\
  ! Get FAC Data
  !/
  if (debug) write(*,*) "read_ring_current"

  call read_ring_current

  IONO_NORTH_JR = IONO_NORTH_JR + IONO_NORTH_RCM_JR
  IONO_SOUTH_JR = IONO_SOUTH_JR + IONO_SOUTH_RCM_JR
  !\
  ! If bc_choice < 0  then use weimer/AMIE over the whole thing
  ! If bc_choice > 0 then use weimer/AMIE as a high lat boundary
  ! If bc_choice == 0 then use no high latitude boundary
  !/
  if (bc_choice /= 0) then

     IONO_NORTH_PHI = 0.0
     IONO_SOUTH_PHI = 0.0

     if (IsFirstTime) then

        if (abs(bc_choice) == 1) then
           Lines(1) = "#BACKGROUND"
           Lines(2) = "EIE/"
           Lines(3) = "weimer96"
           Lines(4) = "ihp"
           Lines(5) = "idontknow"
           Lines(6) = ""
           Lines(7) = "#DEBUG"
           Lines(8) = "2"
           Lines(9) = "0"
           Lines(10) = ""
           Lines(11) = "#END"
        endif

        if (abs(bc_choice) == 2) then
           cAMIEFileNorth = "b20020418.wmsf"
           cAMIEFileSouth = "b20020418.wmsf"
           Lines(1) = "#AMIEFILES"
           Lines(2) = cAMIEFileNorth
           Lines(3) = cAMIEFileSouth
           Lines(4) = ""
           Lines(5) = ""
           Lines(6) = ""
           Lines(7) = "#DEBUG"
           Lines(8) = "2"
           Lines(9) = "0"
           Lines(10) = ""
           Lines(11) = "#END"
        endif

        call EIE_set_inputs(Lines)
        call EIE_Initialize(iError)

        if (iError /=0) then
           write(*,*) "Error in EIE_Initialize, called from get_potential"
           return
        endif

        IsFirstTime = .false.

     endif

     RealTime = RealTime + SimulationTime

     write(*,*) "time: RealTime, SimulationTime", RealTime, SimulationTime

     call IO_SetTime(RealTime)

     TimeOfDay = mod(SimulationTime,24.0*3600.0)
     write(*,*) "Time to AMIE : ",TimeOfDay
     if (abs(bc_choice) == 2) call get_AMIE_values(TimeOfDay)

     call IO_SetIMFBz(imf_bz)
     call IO_SetIMFBy(imf_by)
     call IO_SetSWV  (abs(sw_v)/1000.0)

     ! I think that the IE library assumes that all arrays are 
     ! nLong, nLat.  The ionosphere is nLat, nLong, so we need to
     ! reverse these below...

     call IO_SetnMLTs(IONO_nTheta)
     call IO_SetnLats(IONO_nPsi)

     t = 90.0 - IONO_NORTH_Theta*cRadToDeg
     m = mod(IONO_NORTH_Psi*12.0/cPi + 12.0, 24.0)   

     call IO_SetGrid(m, t, iError)
     call IO_GetPotential(IONO_NORTH_PHI, iError)

     write(*,*) minval(iono_north_phi),maxval(iono_north_phi)

     t = IONO_SOUTH_Theta*cRadToDeg - 90.0
     m = mod(IONO_SOUTH_Psi*12.0/cPi + 12.0, 24.0)

     call IO_SetGrid(m, t, iError)
     call IO_GetPotential(IONO_SOUTH_PHI, iError)

     if (bc_choice < 0) then
        call write_ring_current
        return
     endif

  else

     IONO_NORTH_PHI = 0.0
     IONO_SOUTH_PHI = 0.0

  endif

  if (debug) write(*,*) "write_ring_current"

  call write_ring_current

end subroutine epencalc
