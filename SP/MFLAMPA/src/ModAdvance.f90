module SP_ModAdvance

  ! The module contains methods for advancing the solution in time

  use ModNumConst, ONLY: cTiny

  use ModConst, ONLY: cPi, cMu, cProtonMass, cGyroradius, cLightSpeed, RSun

  use ModCoordTransform, ONLY: rlonlat_to_xyz

  use SP_ModSize, ONLY: iParticleMin, iParticleMax, nMomentumBin, nDim

  use SP_ModGrid, ONLY: &
       X_,Y_,Z_, D_,S_,Rho_,RhoOld_, Ux_,Uy_,Uz_,U_, Bx_,By_,Bz_,B_,BOld_, T_,&
       Begin_, End_, Shock_, ShockOld_, EFlux_,&
       Flux0_, Flux1_, Flux2_, Flux3_, Flux4_, Flux5_, Flux6_, &
       nBlock, &
       State_VIB, Distribution_IIB, iGridLocal_IB, &
       MomentumScale_I, LogMomentumScale_I, EnergyScale_I, LogEnergyScale_I,&
       DMomentumOverDEnergy_I, &
       distance_to_next

  use SP_ModDiffusion, ONLY: advance_diffusion

  use SP_ModLogAdvection, ONLY: advance_log_advection

  implicit none

  SAVE
  
  private ! except

  public:: TimeGlobal, iIterGlobal, &
       set_injection_param, init_advance_const, advance
  public:: DoTraceShock, UseDiffusion

  real, external:: &
       kinetic_energy_to_momentum, momentum_to_kinetic_energy, &
       momentum_to_energy, energy_in


  !\
  ! Global interation and time
  !-----------------------------
  real   :: TimeGlobal  = -1.0
  integer:: iIterGlobal = -1
  !-----------------------------
  ! units of energy
  character(len=*), parameter:: NameEUnit = 'kev'
  real:: UnitEnergy
  ! simulated particles
  character(len=*), parameter:: NameParticle = 'proton'
  !-----------------------------
  ! Injection and max energy in the simulation
  real:: EnergyInj=10.0, EnergyMax=1.0E+07
  real:: TotalEnergyInj
  !-----------------------------
  ! Injection and max momentum in the simulation
  real:: MomentumInj, MomentumMax
  !-----------------------------
  ! Size of a momentum bin on a log-scale
  real:: DLogMomentum
  !-----------------------------
  ! Injection efficiency
  real:: CInj = 1.0
  !-----------------------------
  ! Spectral index in the BC
  real:: SpectralIndex = 5.0
  !-----------------------------
  ! limitation of CFL number
  real:: CFLMax = 0.9
  !-----------------------------
  ! variables shared between subroutines in this module
  integer:: iBegin, iEnd, iBlock, iShock
  real, dimension(iParticleMin:iParticleMax):: Rho_I, RhoOld_I, U_I, T_I
  real, dimension(iParticleMin:iParticleMax):: DLogRho_I
  real, dimension(iParticleMin:iParticleMax):: Radius_I, B_I,   BOld_I
  real, dimension(iParticleMin:iParticleMax):: FermiFirst_I
  !-----------------------------
  ! df/dt = DOuter * d(DInner * df/dx)/dx
  real, dimension(iParticleMin:iParticleMax):: DOuter_I, DInner_I, DInnerInj_I
  !-----------------------------
  ! level of turbulence
  real:: BOverDeltaB2 = 1.0
  !-----------------------------
  real:: MachAlfven
  !-----------------------------
  integer:: nWidth = 50
  !-----------------------------
  logical:: UseRealDiffusionUpstream = .true.
  logical:: DoTraceShock = .true., UseDiffusion = .true.
  !/


contains
  
  subroutine set_injection_param
    use ModReadParam, ONLY: read_var
    !---------------------------------------------
    call read_var('EnergyInj',    EnergyInj)
    call read_var('EnergyMax',    EnergyMax)
    call read_var('SpectralIndex',SpectralIndex)
    call read_var('Efficiency',   CInj)
  end subroutine set_injection_param

  !============================================================================

  subroutine init_advance_const
    ! compute all needed constants
    integer:: iMomentumBin, iBlock, iParticle
    !---------------------------------------------
    ! account for units of energy
    UnitEnergy = energy_in(NameEUnit)
    EnergyInj = EnergyInj * UnitEnergy
    EnergyMax = EnergyMax * UnitEnergy
    ! convert energies to momenta
    MomentumInj  = kinetic_energy_to_momentum(EnergyInj, NameParticle)
    MomentumMax  = kinetic_energy_to_momentum(EnergyMax, NameParticle)
    ! also compute total injection energy
    TotalEnergyInj = momentum_to_energy(MomentumInj, NameParticle)
    DLogMomentum = log(MomentumMax/MomentumInj) / nMomentumBin
    do iMomentumBin = 1, nMomentumBin
       LogMomentumScale_I(iMomentumBin) = &
            log(MomentumInj) + (iMomentumBin-1) * DLogMomentum
       MomentumScale_I(iMomentumBin) = exp(LogMomentumScale_I(iMomentumBin))
       EnergyScale_I(iMomentumBin) = momentum_to_kinetic_energy(&
            MomentumScale_I(iMomentumBin), NameParticle)
       LogEnergyScale_I(iMomentumBin) = log(EnergyScale_I(iMomentumBin))
       DMomentumOverDEnergy_I(iMomentumBin) = &
            momentum_to_energy(MomentumScale_I(iMomentumBin), NameParticle) / &
            (MomentumScale_I(iMomentumBin) * cLightSpeed**2)
    end do
  end subroutine init_advance_const

  !============================================================================
  
  subroutine set_initial_condition
    ! set the initial distribution on all lines
    integer:: iBlock, iParticle, iMomentumBin
    !--------------------------------------------------------------------------
    do iBlock = 1, nBlock
       do iParticle = iGridLocal_IB(Begin_,iBlock), iGridLocal_IB(End_,iBlock)
          do iMomentumBin = 1, nMomentumBin
             Distribution_IIB(iMomentumBin,iParticle,iBlock) = &
                  cTiny / kinetic_energy_to_momentum(EnergyMax,NameParticle)/&
               (MomentumScale_I(iMomentumBin))**2
          end do
       end do
    end do
  end subroutine set_initial_condition

  !============================================================================

  subroutine get_shock_location
    ! find location of a shock wave on every field line
    !--------------------------------------------------------------------------
    ! loop variable
    integer:: iSearchMin, iSearchMax
    integer:: iShockCandidate
    !--------------------------------------------------------------------------
    if(.not.DoTraceShock)then
       iGridLocal_IB(Shock_, iBlock) = iGridLocal_IB(Begin_, iBlock)
       RETURN
    end if

    ! shock front is assumed to be location of max gradient log(Rho1/Rho2);
    ! shock never moves back
    iSearchMin= max(iGridLocal_IB(ShockOld_, iBlock), &
         iGridLocal_IB(Begin_, iBlock) + nWidth )
    iSearchMax= iGridLocal_IB(End_, iBlock) - nWidth - 1
    iShockCandidate = iSearchMin - 1 + maxloc(&
         DLogRho_I(iSearchMin:iSearchMax),&
         1, MASK = Radius_I(iSearchMin:iSearchMax) > 1.2)

    if(DLogRho_I(iShockCandidate) > 0.0)&
         iGridLocal_IB(Shock_, iBlock) = iShockCandidate
  end subroutine get_shock_location

  !===========================================================================

  subroutine get_integral_flux
    ! compute the total (simulated) integral flux of particles as well as
    ! particle flux in the 6 GOES channels; also compute total energy flux
    !------------------------------------------------------------------------
    integer:: iBlock, iParticle, iBin, iFlux ! loop variables
    real   :: EFlux ! the value of energy flux
    real   :: EChannel_I(6) ! energy limits of GOES channels
    real   :: dFlux, dFlux1 ! increments
    real   :: Flux, Flux_I(6) ! the value of particle flux
    real   :: Norm  ! normalization factor
    !-------------------------------------------------------------------------
    ! energy limits of GOES channels
    EChannel_I = (/5,10,30,50,60,100/) * energy_in('MeV')
    ! reset values
    EFlux = 0.0
    Flux_I= 0.0
    Norm  = 0.0
    do iBlock = 1, nBlock
       do iParticle = &
            iGridLocal_IB(Begin_, iBlock), &
            iGridLocal_IB(End_, iBlock)
          !\
          ! Integration loop with midpoint rule
          !/
          do iBin = 1, nMomentumBin - 1
             ! the flux increment from iBin
             dFlux = 0.5 * &
                  (EnergyScale_I(iBin+1) - EnergyScale_I(iBin)) * (&
                  Distribution_IIB(iBin,  iParticle,iBlock)*&
                  MomentumScale_I(iBin)**2 &
                  +&
                  Distribution_IIB(iBin+1,iParticle,iBlock)*&
                  MomentumScale_I(iBin+1)**2)

             ! increase the total flux
             Flux = Flux + dFlux

             ! increase FOES channels' fluxes
             do iFlux = 1, 6
                ! check whether reached the channel's cut-off level
                if(EnergyScale_I(iBin+1) < EChannel_I(iFlux))&
                     CYCLE

                if(EnergyScale_I(iBin+1) >= EChannel_I(iFlux))then
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux
                else
                   ! channel cutoff level is often in the middle of a bin;
                   ! compute partial flux increment
                   dFlux1 =&
                        ((-0.50*(EnergyScale_I(iBin) + EChannel_I(iFlux)) + &
                        EnergyScale_I(iBin+1) )*&
                        Distribution_IIB(iBin,iParticle,iBlock)*&
                        MomentumScale_I(iBin)**2  &
                        -0.50*(EnergyScale_I(iBin)-EChannel_I(iFlux))*&
                        Distribution_IIB(iBin+1,iParticle,iBlock)*&
                        MomentumScale_I(iBin+1)**2)*&
                        (EnergyScale_I(iBin)-EChannel_I(iFlux))/&
                        (EnergyScale_I(iBin+1)-EnergyScale_I(iBin))
                   Flux_I(iFlux) = Flux_I(iFlux) + dFlux1
                end if
             end do

             ! increase total energy flux
             EFlux = EFlux + 0.5 * &
                  (EnergyScale_I(iBin+1) - EnergyScale_I(iBin)) * (&
                  Distribution_IIB(iBin,  iParticle,iBlock)*&
                  EnergyScale_I(iBin) * &
                  MomentumScale_I(iBin)**2 &
                  +&
                  Distribution_IIB(iBin+1,iParticle,iBlock)*&
                  EnergyScale_I(iBin+1) * &
                  MomentumScale_I(iBin+1)**2)

             ! normalization factor
             Norm = Norm + 0.5 * &
                  (MomentumScale_I(iBin+1) - MomentumScale_I(iBin)) * (&
                  Distribution_IIB(iBin,  iParticle, iBlock) * &
                  MomentumScale_I(iBin)**2 &
                  + &
                  Distribution_IIB(iBin+1,iParticle, iBlock) * &
                  MomentumScale_I(iBin+1)**2)
          end do

          ! store the results
          State_VIB(Flux0_,        iParticle, iBlock) = Flux        / Norm
          State_VIB(Flux1_:Flux6_, iParticle, iBlock) = Flux_I(1:6) / Norm
          State_VIB(EFlux_,        iParticle, iBlock) = EFlux       / Norm
       end do
    end do
  end subroutine get_integral_flux

  !============================================================================

  function mach_alfven() result(MachAlfven)
    ! alfvenic mach number for the current line
    real:: MachAlfven

    real:: SpeedAlfven, SpeedUpstream
    !--------------------------------------------
    ! speed upstream is relative to the shock:
    ! \rho_u * (U_u - V_{shock}) = \rho_d * (U_d - V_{shock})
    SpeedUpstream = Rho_I(iShock+1-nWidth)*&
         (U_I(  iShock+1-nWidth) - U_I(  iShock+nWidth))/ &
         (Rho_I(iShock+1-nWidth) - Rho_I(iShock+nWidth))
    SpeedAlfven = &
         B_I(iShock+nWidth) / sqrt(cMu*cProtonMass*Rho_I(iShock+nWidth))
    MachAlfven = SpeedUpstream / SpeedAlfven
  end function mach_alfven

  !===========================================================================

  subroutine steepen_shock
    ! change the density profile near the shock front so it becomes steeper
    ! for the current line
    integer:: iParticle ! loop variable
    real   :: DLogRhoBackground, DLogRhoExcess, Misc, Length
    !--------------------------------------------------------------------------
    ! compute the background value of DLogRho as average in the upstream 
    DLogRhoBackground = 0.0

    ! find the excess of DLogRho within the shock compared to background
    ! averaged over length
    DLogRhoExcess = 0.0
    Length = 0.0
    do iParticle = iShock - nWidth, iShock + nWidth - 1
       Misc = 0.5 * (DLogRho_I(iParticle) + DLogRho_I(iParticle+1)) - &
            DLogRhoBackground
       if(Misc > 0.0)then
          DLogRhoExcess = DLogRhoExcess + &
               Misc * State_VIB(D_, iParticle, iBlock)
          Length = Length + State_VIB(D_, iParticle, iBlock)
       end if
    end do

    ! check for zero excess
    if(DLogRhoExcess == 0.0)&
         RETURN
    ! otherwise, get the averaged value
    DLogRhoExcess = DLogRhoExcess / Length

    ! apply the result within the shock width
    DLogRho_I(iShock-nWidth:iShock+nWidth) = min(&
         DLogRhoBackground, &
         DLogRho_I(iShock-nWidth:iShock+nWidth))
    DLogRho_I(iShock) = DLogRhoExcess + DLogRhoBackground
    
    ! also, sharpen the magnitude of the magnetic field
    ! post shock part
    B_I(  iShock+1-nWidth:iShock+1) = maxval(B_I(  iShock+1-nWidth:iShock+1))
    ! pre shock part
    B_I(  iShock+1:iShock+nWidth) = minval(B_I(  iShock+1:iShock+nWidth))
  end subroutine steepen_shock

  !===========================================================================

  subroutine set_advection_boundary_condition
    ! set boundary conditions on each particle on the current line
    !----------------------------------------
    ! loop variable
    integer:: iParticle
    real:: Density, Momentum
    !----------------------------------------
    do iParticle = iBegin, iEnd
       ! default injection distribution, see Sokolov et al., 2004, eq (3)
       Density = Rho_I(iParticle)
       Momentum= kinetic_energy_to_momentum(&
            State_VIB(T_,iParticle,iBlock)*UnitEnergy,NameParticle)
       Distribution_IIB(1,iParticle,iBlock) = &
            0.25/cPi/(SpectralIndex-3) * CInj * Density / &
            Momentum**3 * (Momentum/MomentumInj)**SpectralIndex
    end do
  end subroutine set_advection_boundary_condition

  !===========================================================================

  subroutine set_diffusion
    ! set diffusion coefficient for the current line
    !--------------------------------------------------------------------------
    DOuter_I(iBegin:iEnd) = B_I(iBegin:iEnd)
    if(.not.UseRealDiffusionUpstream)then
       ! Sokolov et al., 2004: eq (4), 
       ! note: P = Energy * Vel / C**2
       DInnerInj_I(iBegin:iEnd) = &
            BOverDeltaB2*&
            cGyroRadius*(MomentumInj*cLightSpeed)**2/&
            (B_I(iBegin:iEnd)**2*TotalEnergyInj) / RSun**2
    else
       ! diffusion is different up- and down-stream
       ! Sokolov et al. 2004, paragraphs before and after eq (4)
       where(Radius_I(iBegin:iEnd) > 0.9 * Radius_I(iShock))
          ! upstream:
          DInnerInj_I(iBegin:iEnd) = &
               0.2/3.0 * Radius_I(iBegin:iEnd) / RSun * &
               (MomentumInj*cLightSpeed**2)/(B_I(iBegin:iEnd)*TotalEnergyInj)*&
               (MomentumInj*cLightSpeed/energy_in('GeV'))**(1.0/3)
       elsewhere
          ! downstream
          DInnerInj_I(iBegin:iEnd)=&
               cGyroRadius*(MomentumInj*cLightSpeed)**2 / RSun**2/&
               (B_I(iBegin:iEnd)**2 * TotalEnergyInj)/&
               (10.0*CInj*MachAlfven) / &
               min(1.0, 1.0/0.9 * Radius_I(iBegin:iEnd)/Radius_I(iShock))
       end where
    end if

    ! set the boundary condition for diffusion
    Distribution_IIB(2:nMomentumBin, iBegin, iBlock) = &
         Distribution_IIB(1, iBegin, iBlock) * &
         (MomentumScale_I(1) / MomentumScale_I(2:nMomentumBin))**SpectralIndex
  end subroutine set_diffusion

  !===========================================================================

  subroutine set_fermi_first_order
    ! first order Fermi acceleration for the current line
    !--------------------------------------------------------------------------
    FermiFirst_I(iBegin:iEnd) = DLogRho_I(iBegin:iEnd) / (3*DLogMomentum)
 end subroutine set_fermi_first_order

  !===========================================================================

  subroutine advance(TimeLimit)
    ! advance the solution in time
    real, intent(in):: TimeLimit
    integer:: iShockOld, iMomentumBin
    integer:: nProgress, iProgress
    integer:: iStep, nStep
    integer:: iParticle
    real   :: Alpha
    real:: Momentum, DiffCoeffMin =1.0E+04 /RSun
    real:: DtFull, DtProgress, Dt
    logical, save:: IsFirstCall = .true.

    character(len=*), parameter:: NameSub = 'SP:advance'
    !--------------------------------------------------------------------------
    if(IsFirstCall)then
       IsFirstCall = .false.
       call set_initial_condition
    end if

    ! the full time step
    DtFull = TimeLimit - TimeGlobal
    ! go line by line and advance the solution
    do iBlock = 1, nBlock

       ! the active particles on the line
       iBegin = iGridLocal_IB(Begin_, iBlock)
       iEnd   = iGridLocal_IB(End_, iBlock)
       ! various data along the line
       Radius_I(iBegin:iEnd) = &
            sqrt(sum(State_VIB(X_:Z_, iBegin:iEnd, iBlock)**2, 1))
       Rho_I(    iBegin:iEnd) = State_VIB(Rho_,   iBegin:iEnd,iBlock)
       U_I(      iBegin:iEnd) = State_VIB(U_,     iBegin:iEnd,iBlock)
       T_I(      iBegin:iEnd) = State_VIB(T_,     iBegin:iEnd,iBlock)
       BOld_I(   iBegin:iEnd) = State_VIB(BOld_,  iBegin:iEnd,iBlock)
       RhoOld_I( iBegin:iEnd) = State_VIB(RhoOld_,iBegin:iEnd,iBlock)
       DLogRho_I(iBegin:iEnd) = log(Rho_I(iBegin:iEnd)/RhoOld_I(iBegin:iEnd))

       ! identify shock in the data
       call get_shock_location

       ! find how far shock has travelled on this line: nProgress
       iShock    = iGridLocal_IB(Shock_,   iBlock)
       iShockOld = iGridLocal_IB(ShockOld_,iBlock)
       nProgress = MAX(1, iShock - iShockOld)
       iShockOld = MIN(iShockOld, iShock-1)
       iGridLocal_IB(ShockOld_,iBlock) = iShock

       ! each particles shock has crossed should be
       ! processed separately => reduce the time step
       DtProgress = DtFull / nProgress

       ! go over each crossed particle
       do iProgress = 1, nProgress
          ! account for change in the background up to the current moment
          Alpha = real(iProgress) / real(nProgress)
          Rho_I(iBegin:iEnd) = State_VIB(RhoOld_,iBegin:iEnd,iBlock) +Alpha *&
               (State_VIB(Rho_,  iBegin:iEnd,iBlock) - &
               State_VIB(RhoOld_,iBegin:iEnd,iBlock))

          DLogRho_I(iBegin:iEnd)=log(Rho_I(iBegin:iEnd)/RhoOld_I(iBegin:iEnd))

          B_I(iBegin:iEnd) = State_VIB(BOld_,iBegin:iEnd,iBlock) + Alpha*&
               (State_VIB(B_,  iBegin:iEnd,iBlock) - &
               State_VIB(BOld_,iBegin:iEnd,iBlock))
          iShock = iShockOld + iProgress

          ! find the shock alfven mach number, also steepen the shock
          if(iShock < iEnd - nWidth .and. iShock > nWidth)then
             MachAlfven = mach_alfven()
             call steepen_shock
          else
             MachAlfven = 1.0
          end if

          ! 1st order Fermi acceleration is responsible for advection 
          ! in momentum space
          call set_fermi_first_order

          RhoOld_I(iBegin:iEnd) = Rho_I(iBegin:iEnd)
          BOld_I(  iBegin:iEnd) = B_I(  iBegin:iEnd)

          Dt = DtProgress
          nStep = 1
          if(nStep>1)then
             Dt = Dt / nStep
             FermiFirst_I(iBegin:iEnd) = FermiFirst_I(iBegin:iEnd) / nStep
          end if

          ! compute diffusion along the field line
          call set_diffusion
          do iStep = 1, nStep
             ! update bc for advection
             call set_advection_boundary_condition

             ! advection in the momentum space
             do iParticle = iBegin+1, iEnd
                call advance_log_advection(&
                     FermiFirst_I(iParticle),nMomentumBin-2,1,1,&
                     Distribution_IIB(1:nMomentumBin,iParticle,iBlock), &
                     .true.,DLogMomentum)
             end do

             ! diffusion along the field line
             if(.not.UseDiffusion) CYCLE
             do iMomentumBin = 1, nMomentumBin
                Momentum = exp((iMomentumBin-1) * DLogMomentum)
                DInner_I(iBegin:iEnd) =&
                     DInnerInj_I(iBegin:iEnd) * Momentum**2 * &
                     TotalEnergyInj/&
                     momentum_to_energy(Momentum*MomentumInj,NameParticle)
                if(UseRealDiffusionUpstream)then
                   where(Radius_I(iBegin:iEnd) > 0.9*Radius_I(iShock))
                      ! upstream:
                      DInner_I(iBegin:iEnd) = &
                           DInner_I(iBegin:iEnd) / Momentum**(2.0/3)
                   end where
                end if
                DInner_I(iBegin:iEnd) = max(DInner_I(iBegin:iEnd),&
                     DiffCoeffMin/DOuter_I(iBegin:iEnd))
                call advance_diffusion(Dt, iEnd-iBegin+1,&
                     State_VIB(D_,iBegin:iEnd,iBlock), &
                     Distribution_IIB(iMomentumBin,iBegin:iEnd,iBlock),&
                     DOuter_I(iBegin:iEnd), DInner_I(iBegin:iEnd))
             end do

          end do
          
       end do
    end do
    
    call get_integral_flux
  end subroutine advance
  
end module SP_ModAdvance
