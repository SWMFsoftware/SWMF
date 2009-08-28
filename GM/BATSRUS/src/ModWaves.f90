module ModWaves
  use ModVarIndexes, ONLY: nWave, WaveFirst_, WaveLast_
  implicit none
  !
  ! Intended to simulate waves and wave turbulences. The wave propagate
  ! with respect to the background with some speed. The only implemented
  ! case is the propagation with Alfven speed, \pm V_A=\sqrt{B^2/rho}

  logical:: UseAlfvenSpeed = .false.

  integer, private, parameter :: nWaveHalf = max(nWave/2, 1)

  integer, parameter :: AlfvenSpeedPlusFirst_ = WaveFirst_
  integer, parameter :: AlfvenSpeedPlusLast_  = WaveFirst_ - 1 + nWaveHalf

  integer, parameter :: AlfvenSpeedMinusFirst_ = WaveLast_ + 1 - nWaveHalf
  integer, parameter :: AlfvenSpeedMinusLast_  = WaveLast_

  real :: AlfvenSpeed  !Auxiliary variable
  
  real :: FreqMinSI = -1.0
  real :: FreqMaxSI = -1.0

  !To switch this option, set UseAlfvenSpeed = .true. in the user routines 
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

 
contains
  !============================================================================
  subroutine read_alfven_speed
    use ModReadParam,  ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseAlfvenSpeed',UseAlfvenSpeed)
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
    call read_var('FreqMinSI',FreqMinSI)
    call read_var('FreqMaxSI',FreqMaxSI)
  end subroutine read_frequency
end module ModWaves
