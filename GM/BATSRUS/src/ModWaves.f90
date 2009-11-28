module ModWaves
  use ModVarIndexes, ONLY: nWave, WaveFirst_, WaveLast_
  use ModSize,       ONLY: nI, nJ, nK
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

  
  real :: DeltaLogFrequency = 0.0

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
  real :: DivU_C(nI,nJ,nK) = 0.0       !Auxiliary variable

 
contains
  !============================================================================
  subroutine read_alfven_speed
    use ModReadParam,  ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseAlfvenSpeed',UseAlfvenSpeed)
    if(UseAlfvenSpeed.and. FreqMinSI > 0.0.and. nWaveHalf > 1) &
         DeltaLogFrequency = log(FreqMaxSI / FreqMinSI)/nWaveHalf
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
    if(UseAlfvenSpeed.and. FreqMinSI > 0.0.and. nWaveHalf > 1) &
         DeltaLogFrequency = log( FreqMaxSI / FreqMinSI ) / nWaveHalf
    if( (.not.UseAlfvenSpeed).and. FreqMinSI > 0.0.and. nWave > 1) &
         DeltaLogFrequency = log( FreqMaxSI / FreqMinSI ) / nWave
  end subroutine read_frequency
  !============================================================================
  subroutine update_wave_group_advection(iBlock)
    use ModAdvance, ONLY: State_VGB, time_blk
    use ModLinearAdvection, ONLY: advance_lin_advection_plus, &
         advance_lin_advection_minus
    use ModMain, ONLY: CFL
    Use ModFaceValue,ONLY:BetaLimiter
    integer,intent(in)::iBlock

    !\
    !Soution vector, to account for zero boundary condition
    !a singe ayer of the ghost cels is added
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
         .or. (UseAlfvenSpeed .and. nWaveHalf==1) &
         .or.((.not.UseAlfvenSpeed ).and. nWave==1))return
    if(UseAlfvenSpeed)then

       !Boundary conditions at very low and very high frequency:
       F2_I(0) = 0.0; F2_I(nWaveHalf+1) = 0.0

       do k= 1, nK; do j= 1,nJ; do i= 1,nI
          CFL2_I = abs( DivU_C(i,j,k) ) * (GammaWave - 1.0)/&
               DeltaLogFrequency * CFL * time_blk(i,j,k, iBlock)

          if(DivU_C(i,j,k)>0.0)then
             F2_I(nWaveHalf+1) = F2_I(nWaveHalf)
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenSpeedPlusFirst_:AlfvenSpeedPlusLast_, i,j,k, iBlock)
             
             call advance_lin_advection_minus( CFL2_I, nWaveHalf, 1, 1, F2_I,BetaLimiter) 
             State_VGB(AlfvenSpeedPlusFirst_:AlfvenSpeedPlusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)

             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenSpeedMinusFirst_:AlfvenSpeedMinusLast_, i,j,k, iBlock)
             
             call advance_lin_advection_minus( CFL2_I, nWaveHalf, 1, 1, F2_I,BetaLimiter) 
             State_VGB(AlfvenSpeedMinusFirst_:AlfvenSpeedMinusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)
          else
             F2_I(0) = F2_I(1)
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenSpeedPlusFirst_:AlfvenSpeedPlusLast_, i,j,k, iBlock)
             
             call advance_lin_advection_plus( CFL2_I, nWaveHalf, 1, 1, F2_I,BetaLimiter) 
             State_VGB(AlfvenSpeedPlusFirst_:AlfvenSpeedPlusLast_, i,j,k, iBlock) = &
                  F2_I( 1:nWaveHalf)
             
             F2_I( 1:nWaveHalf) = &
                  State_VGB(AlfvenSpeedMinusFirst_:AlfvenSpeedMinusLast_, i,j,k, iBlock)
             
             call advance_lin_advection_plus( CFL2_I, nWaveHalf, 1, 1, F2_I,BetaLimiter) 
             State_VGB(AlfvenSpeedMinusFirst_:AlfvenSpeedMinusLast_, i,j,k, iBlock) = &
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
             
             call advance_lin_advection_minus( CFL_I, nWave, 1, 1, F_I) 
             State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock) = &
                  F_I( 1:nWave)

  
          else
             F_I( 1:nWave) = &
                  State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock)
             
             call advance_lin_advection_plus( CFL_I, nWave, 1, 1, F_I) 
             State_VGB(WaveFirst_:WaveLast_, i,j,k, iBlock) = &
                  F_I( 1:nWave)

          end if
       end do; end do; end do
    end if
       
   
  end subroutine update_wave_group_advection
end module ModWaves
