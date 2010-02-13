module ModWaves
  use ModVarIndexes, ONLY: nWave, WaveFirst_, WaveLast_
  use ModSize,       ONLY: nI, nJ, nK
  implicit none
  !
  ! Intended to simulate waves and wave turbulences. The wave propagate
  ! with respect to the background with some speed. The only implemented
  ! case is the propagation with Alfven speed, \pm V_A=\sqrt{B^2/rho}

  logical:: UseAlfvenWaves = .false.

  integer, private, parameter :: nWaveHalf = max(nWave/2, 1)

  integer, parameter :: AlfvenWavePlusFirst_ = WaveFirst_
  integer, parameter :: AlfvenWavePlusLast_  = WaveFirst_ - 1 + nWaveHalf

  integer, parameter :: AlfvenWaveMinusFirst_ = WaveLast_ + 1 - nWaveHalf
  integer, parameter :: AlfvenWaveMinusLast_  = WaveLast_


  real :: AlfvenSpeed  !Auxiliary variable
  
  real :: FreqMinSI = -1.0
  real :: FreqMaxSI = -1.0
  
  !Centered frequency, for each bin
  !The logarithmic grid is set by default. Can be reset in use_user_perturbation
  real, dimension(WaveFirst_:WaveLast_):: FrequencySi_W = -1.0

  
  real :: DeltaLogFrequency = 0.0

  !To switch this option, set UseAlfvenWaves = .true.  
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

  !For multi-component wave pressure the limitatation of different
  !frequency groups does not ensure the proper limitation of the total 
  !pressure therefore, it is reasonble to collect the total wave pressure
  !and limit it accordingly 
  logical:: UseWavePressureLtd = .false.
  real :: GammaWave = 1.50
 
  !Spectral functions: all functions are normalized by unity

  real,dimension(WaveFirst_:WaveLast_):: Spectrum_W  = 1.0/nWave , &
                          SpectrumPlus_W  = 2.0/nWave , &
                          SpectrumMinus_W = 2.0/nWave

  character(LEN=10)::  NameSpectralFunction = 'uniform'

  !Parameters for different kinds of the specrral function
  !For Planckian: not implemented
  !real:: TRadSpectrum, EnergyMin

  !For power law: I\propto 1/f^{PowerIndex}, f> FreqStartSi
  real:: PowerIndex = 5.0/3.0, FreqStartSi = -1.0

  real :: WaveEnergy = 0.0 !Auxiliary variable
  real :: DivU_C(nI,nJ,nK) = 0.0       !Auxiliary variable

 
contains
  subroutine read_spectrum
    use ModReadParam,  ONLY: read_var
    
  end subroutine read_spectrum
  !============================================================================
  subroutine check_waves
    integer:: iWave
    !--------------------------------------------------------------------------
    if(UseAlfvenWaves .and. 2*nWaveHalf /= nWave)call stop_mpi(&
         'nWave should be positive even when Alfven waves are used')

    if(UseAlfvenWaves          .and. &
         FreqMinSI > 0.0       .and. &
         FreqMaxSi > FreqMinSi .and. &
         nWaveHalf > 1             ) &
         DeltaLogFrequency = log( FreqMaxSI / FreqMinSI ) / nWaveHalf

    if( (.not.UseAlfvenWaves)  .and. &
         FreqMinSI > 0.0       .and. &
         FreqMaxSi > FreqMinSi .and. &
         nWave > 1           ) &
         DeltaLogFrequency = log( FreqMaxSI / FreqMinSI ) / nWave

    if(FreqMinSi > 0           .and. &
       FreqMaxSi > FreqMinSi   .and. &
       DeltaLogFrequency > 0.0)then
       
       !set bin-centered frequencies

       if(UseAlfvenWaves)then

          FrequencySi_W(AlfvenWavePlusFirst_) = &
               FreqMinSi * exp( 0.50 * DeltaLogFrequency)

          do iWave = AlfvenWavePlusFirst_+ 1, AlfvenWavePlusLast_
             FrequencySi_W(iWave) = FrequencySi_W(iWave - 1) * &
                  exp( DeltaLogFrequency )
          end do

          FrequencySi_W(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) = &
               FrequencySi_W(AlfvenWavePlusFirst_:AlfvenWavePlusLast_)
       else

          FrequencySi_W(WaveFirst_) = FreqMinSi * exp( 0.50 * DeltaLogFrequency)

          do iWave = WaveFirst_ + 1, WaveLast_
             FrequencySi_W(iWave) = FrequencySi_W(iWave - 1) * &
                  exp( DeltaLogFrequency )
          end do
       end if
    end if
    !Set spectral functions

    SpectrumMinus_W = 0.0; SpectrumPlus_W = 0.0

    select case(NameSpectralFunction)
       
    case('uniform')
       Spectrum_W = 1.0/nWave
       
       if(UseAlfvenWaves)then        
          SpectrumPlus_W(AlfvenWavePlusFirst_:AlfvenWavePlusLast_) &
               = 1.0/nWaveHalf
          SpectrumMinus_W(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) &
               = 1.0/nWaveHalf
       end if
    case('powerlaw')
       if(any(FrequencySi_W <= 0.0))call stop_mpi(&
            'Power law spectrum cannot be set for not-positive frequencies')
       if(all(FrequencySi_W < FreqStartSi))call stop_mpi(&
            'All bin-centered frequencies are below the start frequency for spectrum')

       where(FrequencySi_W >= FreqStartSi)&
            Spectrum_W = 1.0 / FrequencySi_W**(PowerIndex - 1.0)

       !Normalize by unity:
       Spectrum_W = Spectrum_W/sum(Spectrum_W)

       if(UseAlfvenWaves)then

          where(FrequencySi_W(&
               AlfvenWavePlusFirst_:AlfvenWavePlusLast_) >= FreqStartSi)&
               SpectrumPlus_W(&
               AlfvenWavePlusFirst_:AlfvenWavePlusLast_) = 1.0 / &
               FrequencySi_W(&
               AlfvenWavePlusFirst_:AlfvenWavePlusLast_)**(PowerIndex - 1.0)

          where(FrequencySi_W(&
               AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) >= FreqStartSi)&
               SpectrumMinus_W(&
               AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) = 1.0 / &
               FrequencySi_W(&
               AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_)**(PowerIndex - 1.0)
          
          !Normalize by unity:
          SpectrumPlus_W  = SpectrumPlus_W  /sum(SpectrumPlus_W)
          SpectrumMinus_W = SpectrumMinus_W /sum(SpectrumMinus_W)
       end if
    case default
       call stop_mpi(&
            'Spectral function '//NameSpectralFunction//' is not implemented')
    end select
  end subroutine check_waves
  !============================================================================
  subroutine read_wave_pressure
    use ModReadParam,  ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseWavePressure'   ,UseWavePressure)
    call read_var('UseWavePressureLtd',UseWavePressureLtd)
  end subroutine read_wave_pressure
  !============================================================================
  subroutine read_frequency
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('FreqMinSI',FreqMinSI)
    call read_var('FreqMaxSI',FreqMaxSI)
  end subroutine read_frequency
  !============================================================================
  subroutine update_wave_group_advection(iBlock)
    use ModAdvance,           ONLY: State_VGB, time_blk
    use ModGeometry,          ONLY: true_cell
    use ModLinearAdvection,   ONLY: advance_lin_advection_plus, &
                                    advance_lin_advection_minus
    use ModMain,              ONLY: CFL
    Use ModFaceValue,         ONLY:BetaLimiter
   
    integer,intent(in)    ::iBlock

    !\
    !Soution vector, to account for zero boundary condition
    !a single layer of the ghost cels is added
    !/
    real:: F_I(0:nWave+1), F2_I( 0: nWaveHalf+1)
    

    !\
    !Auxiiary vector of CFL numbers:
    !/
    real:: CFL_I(1:nWave), CFL2_I(1:nWaveHalf)
    
    !\
    !Loop variables:
    !/
    integer :: i,j,k
    !---------------------------
    if(DeltaLogFrequency<= 0.0                   &
         .or. (UseAlfvenWaves .and. nWaveHalf==1) &
         .or.((.not.UseAlfvenWaves ).and. nWave==1))return
    if(UseAlfvenWaves)then

       do k= 1, nK; do j= 1,nJ; do i= 1,nI
          if(.not. true_cell(i,j,k,iBlock)) CYCLE
          CFL2_I = abs( DivU_C(i,j,k) ) * (GammaWave - 1.0)/&
               DeltaLogFrequency * CFL * time_blk(i,j,k, iBlock)
          ! Boundary conditions
          F2_I(0) = 0.0 ; F2_I(nWaveHalf+1) = 0.0
          if(DivU_C(i,j,k)>0.0)then
         
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenWavePlusFirst_:AlfvenWavePlusLast_, i,j,k, iBlock)
             F2_I(nWaveHalf+1)=F2_I(nWaveHalf)
             call advance_lin_advection_minus( CFL2_I, nWaveHalf, 1, 1, F2_I, &
                  BetaLimiter, UseConservativeBC= .true.) 
             State_VGB(AlfvenWavePlusFirst_:AlfvenWavePlusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)

             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_, i,j,k, iBlock)
             F2_I(nWaveHalf+1) = F2_I(nWaveHalf)
             call advance_lin_advection_minus( CFL2_I, nWaveHalf, 1, 1, F2_I, &
                  BetaLimiter, UseConservativeBC= .true.) 
             State_VGB(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)
          else
         
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenWavePlusFirst_:AlfvenWavePlusLast_, i,j,k, iBlock)
             F2_I(0) = F2_I(1) 
             call advance_lin_advection_plus( CFL2_I, nWaveHalf, 1, 1, F2_I, &
                  BetaLimiter, UseConservativeBC= .true.) 
             State_VGB(AlfvenWavePlusFirst_:AlfvenWavePlusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)
             
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_, i,j,k, iBlock)
             F2_I(0) = F2_I(1)
             call advance_lin_advection_plus( CFL2_I, nWaveHalf, 1, 1, F2_I, &
                  BetaLimiter, UseConservativeBC= .true.)
             State_VGB(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)
          end if
          
          
       end do; end do; end do
    else
       !Boundary conditions at very low and very high frequency:
       F_I(0) = 0.0; F_I(nWave+1) = 0.0
       do k= 1, nK; do j= 1,nJ; do i= 1,nI
          
          
          CFL_I = abs( DivU_C(i,j,k) ) * (GammaWave - 1.0)/&
               DeltaLogFrequency * CFL * time_blk(i,j,k, iBlock)
          
          if(DivU_C(i,j,k)>0.0)then

             F_I( 1:nWave) = &
                  State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock)
             
             call advance_lin_advection_minus( CFL_I, nWave, 1, 1, F_I, UseConservativeBC= .true.) 
             State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock) = &
                  F_I( 1:nWave)

  
          else
             F_I( 1:nWave) = &
                  State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock)
             
             call advance_lin_advection_plus( CFL_I, nWave, 1, 1, F_I, UseConservativeBC= .true.)
             State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock) = &
                  F_I( 1:nWave)

          end if
       end do; end do; end do
    end if
       
   
  end subroutine update_wave_group_advection
end module ModWaves
