!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModDistribution
  ! The module contains the velocity (momentum) distribution function
  ! and methods for initializing it as well as the offset routine.
  use ModNumConst,ONLY: cTiny
  use ModConst,   ONLY: cLightSpeed, energy_in
  use SP_ModSize, ONLY: nParticleMax, nP=>nMomentum
  use SP_ModUnit, ONLY: UnitEnergy=>UnitParticleEnergy, &
       kinetic_energy_to_momentum, momentum_to_energy
  use SP_ModGrid, ONLY: nBlock, nParticle_B
  implicit none
  SAVE
  PRIVATE ! except
  !Public members:
  public:: init              !Initialize Distribution_IIB
  public:: offset
  public:: get_integral_flux
  public:: nP,  TotalEnergyInj, MomentumInj, MomentumMax,  DLogP, &
       EnergyInj, EnergyMax
  !\
  !!!!!!!!!!!!!!!Grid along the nomentum axis              !!!!!!
  ! Injection and maximal energy in the simulation
  ! To be read from the PARAM.in file: KINETIC energies
  real:: EnergyInj=10.0, EnergyMax=1.0E+07
  ! Injection and max momentum in the simulation
  real:: MomentumInj, MomentumMax
  ! Total energy, including the rest mass energy. Velocity in terms
  ! of momentum and total energy equals:
  ! velocity =  momentum*cLightSpeed**2/TotalEnergy
  real:: TotalEnergyInj
  ! Size of a  log-momentum mesh. For momentum we use both the 
  ! denotaion, P, and a word, momentum - whichever is more covenient
  real:: DLogP        !log(MomentumMax/MomentumInj)/nP
  !/
  !\
  ! Functions to convert the grid index to momentum and energy
  !/
  !\
  ! Momentum and energy at the grid points
  real, public:: Momentum_I(0:nP+1)
  real, public:: Energy_I(0:nP+1)
  ! Total energy, including the rest mass energy
  real, public:: TotalEnergy_I(0:nP+1)
  !/
  !\
  !!!!!!!!!!!!!!Grid in the momentum space                 !!!!!!
  !iP     0     1                         nP   nP+1
  !       |     |    ....                 |     | 
  !P      P_inj P_inj*exp(\Delta(log P))  P_Max P_Max*exp(\Delta(log P))
  ! This is because we put two boundary conditions: the background value at
  ! the right one and the physical condition at the left one, for the 
  ! velocity distribution function
  !/
  !\
  ! Velosity Distribution Function (VDF) 
  ! Number of points along the momentum axis is set in ModSize
  ! 1st index - log(momentum)
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, public, allocatable:: Distribution_IIB(:,:,:)
  !/
  logical :: DoInit = .true.
contains
  !=================================================================
  subroutine init
    use SP_ModUnit,   ONLY: momentum_to_kinetic_energy
    use ModUtilities, ONLY: check_allocate
    ! set the initial distribution on all lines
    integer:: iBlock, iParticle, iP, iError
    !----------------------------------------------------------
    if(.not.DoInit)RETURN
    DoInit = .false.
    ! convert energies to momenta
    MomentumInj  = kinetic_energy_to_momentum(EnergyInj*UnitEnergy)
    MomentumMax  = kinetic_energy_to_momentum(EnergyMax*UnitEnergy)
    ! total injection energy (including the rest mass energy
    TotalEnergyInj = momentum_to_energy(MomentumInj)
    DLogP = log(MomentumMax/MomentumInj)/nP
    !\
    ! Functions to convert the grid index to momentum and energy
    !/
    do iP = 0, nP +1
       Momentum_I(iP) = MomentumInj*exp(iP*DLogP)
       Energy_I(iP) = momentum_to_kinetic_energy(Momentum_I(iP))
       TotalEnergy_I(iP) = momentum_to_energy(Momentum_I(iP))
    end do
    !\
    ! Distribution function
    !/    
    allocate(Distribution_IIB(&
         0:nP+1,1:nParticleMax,nBlock), stat=iError)
    call check_allocate(iError, 'Distribution_IIB')
    ! initialization depends on momentum, however, this corresponds
    ! to a constant differential flux (intensity), thus ensuring
    ! uniform backgound while visualizing this quantity
    do iBlock = 1, nBlock
       do iParticle = 1, nParticleMax
          ! Overall density of the fast particles is of the order 
          ! of 10^-6 m^-3. Integral flux is less than 100 per
          ! (m^2 ster s). Differential background flux is constant.
          do iP = 0, nP +1
             Distribution_IIB(iP,iParticle,iBlock) = 0.10*&
                  cTiny/(MomentumMax*Momentum_I(iP)**2)
          end do
       end do
    end do
  end subroutine init
  !============================================================
  subroutine offset(iBlock, iOffset, AlphaIn)
    use SP_ModGrid, ONLY: &
         BOld_, RhoOld_, nParticle_B, ShockOld_, iShock_IB,  State_VIB
    ! shift in the data arrays is required if the grid point(s) is  
    ! appended or removed at the foot point of the magnetic field line
    !SHIFTED ARE:  State_VIB((/RhoOld_,BOld_/),:,:), Distribution_IIB,
    !ShockOld, nParticle_B
    integer, intent(in)        :: iBlock
    integer, intent(in)        :: iOffset
    real, optional, intent(in) :: AlphaIn
    real :: Alpha
    character(len=*), parameter :: NameSub = "SP: offset"
    !------------
    Alpha = 0
    if(present(AlphaIn))Alpha=AlphaIn
    if(iOffset==0)then
       RETURN
    elseif(iOffset==1)then
       State_VIB((/RhoOld_,BOld_/),2:nParticle_B(iBlock)+1,iBlock) &
            = State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock),iBlock)
       Distribution_IIB(:,2:nParticle_B( iBlock)+iOffset, iBlock)&
            = Distribution_IIB(:,1:nParticle_B( iBlock), iBlock)
       State_VIB((/RhoOld_, BOld_/), 1, iBlock) = &
            (Alpha + 1)*State_VIB((/RhoOld_, BOld_/), 2, iBlock) &
            -Alpha     * State_VIB((/RhoOld_, BOld_/), 3, iBlock)
       Distribution_IIB(:,1,iBlock) = Distribution_IIB(:,2,iBlock) + &
            Alpha*(Distribution_IIB(:,2,iBlock) - Distribution_IIB(:,3,iBlock))
    elseif(iOffset < 0)then
       State_VIB((/RhoOld_,BOld_/),1:nParticle_B(iBlock)+iOffset,iBlock) &
            =  State_VIB((/RhoOld_,BOld_/),1-iOffset:nParticle_B(iBlock),&
            iBlock)
       Distribution_IIB(:,1:nParticle_B( iBlock)+iOffset, iBlock)&
            = Distribution_IIB(:,1-iOffset:nParticle_B( iBlock), iBlock)
    else
       call CON_stop('No algorithm for iOffset >1 in '//NameSub)
    end if
    if(iShock_IB(ShockOld_, iBlock)>1)&
         iShock_IB(ShockOld_, iBlock) = &
         max(iShock_IB(ShockOld_, iBlock) + iOffset, 1)
    nParticle_B(iBlock) = nParticle_B( iBlock) + iOffset
  end subroutine offset
  !===========================================================================
  subroutine get_integral_flux
    use SP_ModGrid, ONLY: EFlux_, Flux0_, Flux1_, Flux2_, Flux3_, Flux4_,&
         Flux5_, Flux6_, FluxMax_, Flux_VIB 
    use ModConst, ONLY: energy_in
    ! compute the total (simulated) integral flux of particles as well as
    ! particle flux in the 6 GOES channels; also compute total energy flux
    !------------------------------------------------------------------------
    integer:: iBlock, iParticle, iP, iFlux ! loop variables
    real   :: EFlux ! the value of energy flux
    real   :: EChannel_I(6) ! energy limits of GOES channels
    real   :: dFlux, dFlux1 ! increments
    real   :: Flux, Flux_I(6) ! the value of particle flux
    real   :: Norm  ! normalization factor
    !-------------------------------------------------------------------------
    ! energy limits of GOES channels
    EChannel_I = (/5,10,30,50,60,100/) * energy_in('MeV')
    do iBlock = 1, nBlock
       do iParticle = 1, nParticle_B( iBlock)
          !\
          ! Integration loop with midpoint rule
          !/
          ! reset values
          EFlux = 0.0
          Flux_I= 0.0
          Norm  = 0.0
          Flux  = 0.0
          do iP = 1, nP - 1
             ! the flux increment from iP
             dFlux = 0.5 * &
                  (Energy_I(iP+1) - Energy_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle,iBlock)*&
                  Momentum_I(iP)**2 &
                  +&
                  Distribution_IIB(iP+1,iParticle,iBlock)*&
                  Momentum_I(iP+1)**2)

             ! increase the total flux
             Flux = Flux + dFlux

             ! increase FOES channels' fluxes
             do iFlux = 1, 6
                ! check whether reached the channel's cut-off level
                if(Energy_I(iP+1) < EChannel_I(iFlux))&
                     CYCLE

                if(Energy_I(iP+1) >= EChannel_I(iFlux))then
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux
                else
                   ! channel cutoff level is often in the middle of a bin;
                   ! compute partial flux increment
                   dFlux1 =&
                        ((-0.50*(Energy_I(iP) + EChannel_I(iFlux)) + &
                        Energy_I(iP+1) )*&
                        Distribution_IIB(iP,iParticle,iBlock)*&
                        Momentum_I(iP)**2  &
                        -0.50*(Energy_I(iP)-EChannel_I(iFlux))*&
                        Distribution_IIB(iP+1,iParticle,iBlock)*&
                        Momentum_I(iP+1)**2)*&
                        (Energy_I(iP)-EChannel_I(iFlux))/&
                        (Energy_I(iP+1)-Energy_I(iP))
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux1
                end if
             end do

             ! increase total energy flux
             EFlux = EFlux + 0.5 * &
                  (Energy_I(iP+1) - Energy_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle,iBlock)*&
                  Energy_I(iP) * &
                  Momentum_I(iP)**2 &
                  +&
                  Distribution_IIB(iP+1,iParticle,iBlock)*&
                  Energy_I(iP+1) * &
                  Momentum_I(iP+1)**2)

             ! normalization factor
             Norm = Norm + 0.5 * &
                  (Momentum_I(iP+1) - Momentum_I(iP)) * (&
                  Distribution_IIB(iP,  iParticle, iBlock) * &
                  Momentum_I(iP)**2 &
                  + &
                  Distribution_IIB(iP+1,iParticle, iBlock) * &
                  Momentum_I(iP+1)**2)
          end do

          ! store the results
          Flux_VIB(Flux0_,        iParticle, iBlock) = Flux        / Norm
          Flux_VIB(Flux1_:Flux6_, iParticle, iBlock) = Flux_I(1:6) / Norm
          Flux_VIB(EFlux_,        iParticle, iBlock) = EFlux       / Norm
       end do
    end do
  end subroutine get_integral_flux
end module SP_ModDistribution
