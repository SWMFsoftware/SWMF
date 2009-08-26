module ModWaves
  implicit none
  !
  ! Intended to simulate waves and wave turbulences. The wave propagate
  ! with respect to the background with some speed. The only implemented
  ! case is the propagation with Alfven speed, \pm V_A=\sqrt{B^2/rho}

  logical:: UseAlfvenSpeed = .false.

  integer :: AlfvenSpeedPlusFirst_ = 2
  integer :: AlfvenSpeedPlusLast_  = 1

  integer :: AlfvenSpeedMinusFirst_ = 2
  integer :: AlfvenSpeedMinusLast_  = 1
  real :: AlfvenSpeed  !Auxiliary variable
  
  real :: FreqMin = 0.0
  real :: FreqMax = 0.0
  real :: FreqInertialRange = 0.0

  !To swotch this option, set UseAlfvenSpeed = .true. in the user routines 
  !and modify accordingly the named indexes to set the state variables for
  !which the advection with \pm V_A should be applied

  !The background affects the wave propagation and, accordingly, the wave
  !should affect the background motion. The only case is implemented:
  !wave contribute to the total pressure, the contribution being equal
  !to (gamma_W-1)\sum{state vector components assigned to the waves}
  !This assumes a proper normalization of the wave state variables
  !(they should have a sense of the wave energy, should be multiplied if
  !needed by \omega, \Delta \log \omega etc (Rona's {\cal I} is an example).

  !Again, the wave pressure logicals and named indexes should be modified
  !only via the user routines.
  
  logical:: UseWavePressure = .false.

  real :: GammaWave = 1.50
 
 
  real :: WaveEnergy = 0.0 !Auxiliary variable
  real :: DivU = 0.0       !Auxiliary variable

  character(len=*),dimension(0:99), parameter :: NameNumber_I=(/ &
       '00','01','02','03','04','05','06','07','08','09', &
       '10','11','12','13','14','15','16','17','18','19', &
       '20','21','22','23','24','25','26','27','28','29', &
       '30','31','32','33','34','35','36','37','38','39', &
       '40','41','42','43','44','45','46','47','48','49', &
       '50','51','52','53','54','55','56','57','58','59', &
       '60','61','62','63','64','65','66','67','68','69', &
       '70','71','72','73','74','75','76','77','78','79', &
       '80','81','82','83','84','85','86','87','88','89', &
       '90','91','92','93','94','95','96','97','98','99'     /)

contains
  !============================================================================
  subroutine read_alfven_speed
    use ModReadParam,  ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseAlfvenSpeed',UseAlfvenSpeed)
    if(UseAlfvenSpeed)then
       call read_var(' AlfvenSpeedPlusFirst' , AlfvenSpeedPlusFirst_ )
       call read_var(' AlfvenSpeedPlusLast'  , AlfvenSpeedPlusLast_  )
       call read_var(' AlfvenSpeedMinusFirst', AlfvenSpeedMinusFirst_)
       call read_var(' AlfvenSpeedMinusLast' , AlfvenSpeedMinusLast_ )
    end if
  end subroutine read_alfven_speed
  !============================================================================
  subroutine read_wave_pressure
    use ModReadParam,  ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseWavePressure',UseWavePressure)
  end subroutine read_wave_pressure
  !============================================================================
  subroutine read_frequency
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('FreqMin',FreqMin)
    call read_var('FreqMax',FreqMax)
    call read_var('FreqInertialRange',FreqInertialRange)
  end subroutine read_freq_grid
end module ModWaves
