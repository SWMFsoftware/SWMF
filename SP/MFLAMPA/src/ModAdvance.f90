!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_ModAdvance
  ! The module contains methods for advancing the solution in time
  use ModNumConst,ONLY: cPi
  use ModConst,   ONLY: cLightSpeed, energy_in
  use SP_ModSize, ONLY: nParticleMax
  use SP_ModDistribution, ONLY: nP, Distribution_IIB, Momentum_I,&
        TotalEnergy_I, MomentumMax, MomentumInj,&
        DLogP, EnergyInjIo, EnergyMaxIo
  use SP_ModGrid, ONLY: State_VIB, iShock_IB,   R_,   &
       Shock_, ShockOld_, DLogRho_, nBlock, nParticle_B 
  implicit none
  SAVE
  PRIVATE ! except
  !Public members:
  public:: init        !Initialize reused variables 
  public:: read_param  !read injection parameters
  public:: advance     !Advance solution Distribution_IIB
  !\
  !||||||||||||||Boundary condition at the injection energy!!!!!!
  ! Injection efficiency and assumed spectral index with the energy
  ! range k_BT_i< Energy < EnergyInjection, to be read from PARAM.in 
  real:: CInj = 1.0, SpectralIndex = 5.0
  !/
  !\
  ! Total energy, including the rest mass energy. Velocity in terms
  ! of momentum and total energy equals:
  ! velocity =  momentum*cLightSpeed**2/TotalEnergy
  real:: TotalEnergyInj
  !/
  !\
  !!!!!!!!!!!!!!!!!!!!!!!!!Local parameters!!!!!!!!!!!!!!!
  real:: CFL=0.9        !Controls the maximum allowed time step
  ! level of turbulence: Ratio of regular to irregual magnetic field,
  ! squared. 
  real, parameter    :: BOverDeltaB2 = 1.0
  !\
  integer, public, parameter :: nWidth = 50
  !/
  !\
  logical:: UseRealDiffusionUpstream = .true.
  logical, public:: DoTraceShock = .true., UseDiffusion = .true.
  !/
  logical :: DoInit = .true.
contains
  subroutine init
    !----------------------------------------------------------
    if(.not.DoInit)RETURN
    DoInit = .false.
    ! total injection energy (including the rest mass energy
    TotalEnergyInj =  TotalEnergy_I(0) 
  end subroutine init
  !===========================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in  
    character(len=*), parameter :: NameSub='SP:read_param_adv'
    !---------------------------------------------
    select case(NameCommand)
    case('#INJECTION')
       call read_var('EnergyInjIo',    EnergyInjIo)
       call read_var('EnergyMaxIo',    EnergyMaxIo)
       call read_var('SpectralIndex',SpectralIndex)
       call read_var('Efficiency',   CInj)
    case('#CFL')
       call read_var('Cfl',CFL)
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================
  subroutine advance(TimeLimit)
    ! advance the solution of the diffusive kinetic equation:               
    !            f_t+[(1/3)*(d(ln rho)/dt]*f_{ln p}=B*d/ds[D/B*df/ds]
    ! with accounting for diffusion and Fermi acceleration 
    ! from SPTime to TimeLimit
    ! Prototype: FLAMPA/src/SP_main, case("RUN"), Roussev&Sokolov2008
    ! Version: Borovikov&Sokolov, Dec.19 2017, distinctions:
    ! (1) no turbulence (2) new shock finder moved to SP_ModMain,
    ! and (3) new steepen_shock
    use SP_ModTime,         ONLY: SPTime
    use SP_ModDiffusion,    ONLY: advance_diffusion
    use SP_ModLogAdvection, ONLY: advance_log_advection
    use ModConst,           ONLY: cMu, cProtonMass, cGyroradius, RSun
    use SP_ModGrid,         ONLY: D_, Rho_, RhoOld_, B_, BOld_, U_, T_
    real, intent(in):: TimeLimit
    !\
    ! Loop variables
    integer  :: iP, iParticle, iBlock
    !/
    !\
    !For a given block: nParticle_B, iShock_IB:
    integer  :: iEnd, iShock, iShockOld
    !/
    !\
    ! Upper limit and variable for the Loop which makes
    ! a time step so short that the shock wave passes
    ! a single grid interval per a progress step:
    integer  :: iProgress, nProgress
    ! Upper limit and variable for the loop which makes a 
    ! time step so short that the CFL for the (explicit) 
    ! advection is less that CFL declared abobe. 
    integer  :: iStep, nStep
    !/
    !\
    ! coefficient to interpolate  "old" and "new"
    real     :: Alpha
    !/
    !\
    ! Lower imit to floor the spatial diffusion coefficient For a
    ! given spatial and temporal resolution, the value of the  
    ! diffusion coefficient should be artificially increased to get       
    ! the diffusion length to be larger than the shock front width, 
    ! which even for the steepened shock is as wide as a mesh size
    ! of the Largangian grid, State_VIB(D_,:,:)*RSun. In this way, 
    ! the Lagrangian grid resolution is sufficient to resolve a 
    ! precursor in the upstream distribution of 
    ! DiffCoeffMin=0 (default value), we do NOT enhance   
    ! the diffusion coefficient! Physically, DiffCoeffMin 
    ! should be given by the product of shock wave speed  
    ! and local grid spacing. 
    real, parameter::  DiffCoeffMin =1.0E+04 /RSun
    ! Full difference between DataInputTime and SPTime
    real      :: DtFull
    ! Time step in the PROGRESS Loop, DtFull/nProgress
    real      :: DtProgress 
    ! Time step in the STEP Loop, DtProgress/nStep
    real      :: Dt
    real      :: MachAlfven
    !\
    ! Local arrays to store the state vector
    ! Heliocentric distance, in R_S
    real      :: Radius_I(1:nParticleMax)
    ! Plasma speed, [m/s]
    real      :: U_I(     1:nParticleMax)
    ! Plasma density, [a.m.u./m3]
    real, dimension(1:nParticleMax):: Rho_I, RhoOld_I
    ! Interplanetary magnetic field, [T]
    real, dimension(1:nParticleMax):: B_I, BOld_I
    ! Lagrangian derivatives 
    real, dimension(1:nParticleMax):: DLogRho_I, FermiFirst_I
    !\
    ! Coefficients in the diffusion operator
    ! df/dt = DOuter * d(DInner * df/dx)/dx
    real, dimension(1:nParticleMax):: &
         DOuter_I, DInner_I, DInnerInj_I
    !/
    character(len=*), parameter:: NameSub = 'SP:advance'
    !---------------------------------------------------------------
    ! the full time interval
    DtFull = TimeLimit - SPTime
    ! go line by line and advance the solution
    BLOCK:do iBlock = 1, nBlock

       ! the active particles on the line
       iEnd   = nParticle_B( iBlock)
       ! various data along the line
       Radius_I( 1:iEnd) = State_VIB(R_,      1:iEnd,iBlock)
       U_I(      1:iEnd) = State_VIB(U_,      1:iEnd,iBlock)
       BOld_I(   1:iEnd) = State_VIB(BOld_,   1:iEnd,iBlock)
       RhoOld_I( 1:iEnd) = State_VIB(RhoOld_, 1:iEnd,iBlock)

       ! find how far shock has travelled on this line: nProgress
       iShock    = iShock_IB(Shock_,   iBlock)
       iShockOld = iShock_IB(ShockOld_,iBlock)
       if(DoTraceShock)then
          nProgress = MAX(1, iShock - iShockOld)
          iShockOld = MIN(iShockOld, iShock-1)
       else
          nProgress = 1
          iShockOld = 0
       end if

       ! each particles shock has crossed should be
       ! processed separately => reduce the time step
       DtProgress = DtFull/nProgress

       ! go over each crossed particle
       PROGRESS:do iProgress = 1, nProgress
          ! account for change in the background up to the current moment
          Alpha = real(iProgress) / real(nProgress)
          Rho_I(1:iEnd) = State_VIB(RhoOld_,1:iEnd,iBlock) +Alpha *&
               (State_VIB(Rho_,  1:iEnd,iBlock) - &
               State_VIB(RhoOld_,1:iEnd,iBlock))

          DLogRho_I(1:iEnd)=log(Rho_I(1:iEnd)/RhoOld_I(1:iEnd))

          B_I(1:iEnd) = State_VIB(BOld_,1:iEnd,iBlock) + Alpha*&
               (State_VIB(B_,  1:iEnd,iBlock) - &
               State_VIB(BOld_,1:iEnd,iBlock))
          iShock = iShockOld + iProgress

          ! find the shock alfven mach number, also steepen the shock
          if(iShock < iEnd - nWidth.and.iShock > nWidth.and.DoTraceShock)then
             MachAlfven = mach_alfven()
             call steepen_shock
          else
             MachAlfven = 1.0
          end if

          ! 1st order Fermi acceleration is responsible for advection 
          ! in momentum space
          ! first order Fermi acceleration for the current line
          !--------------------------------------------------------
          FermiFirst_I(1:iEnd) = DLogRho_I(1:iEnd) / (3*DLogP)
          !\
          ! Check if the number of time steps is greater than 1:
          !/
          nStep = 1+int(maxval(abs(FermiFirst_I(2:iEnd)))/CFL)
          RhoOld_I(1:iEnd) = Rho_I(1:iEnd)
          BOld_I(  1:iEnd) = B_I(  1:iEnd)

          Dt = DtProgress
          if(nStep>1)then 
             Dt = Dt / nStep
             FermiFirst_I(1:iEnd) = FermiFirst_I(1:iEnd) / nStep
          end if
          ! compute diffusion along the field line
          ! we calculate: "Outer diffusion"=B_I and
          ! "Inner diffusion at the injection energy (iP=0)
          call set_diffusion
          STEP:do iStep = 1, nStep !Currently nStep = 1
             ! update bc for advection
             call set_advection_bc

             ! advection in the momentum space
             do iParticle = 2, iEnd
                call advance_log_advection(&
                     FermiFirst_I(iParticle), nP, 1, 1,        &
                     Distribution_IIB(0:nP+1,iParticle,iBlock),&
                     .false.)
             end do
             if(.not.UseDiffusion) CYCLE STEP 
             ! diffusion along the field line
             MOMENTUM:do iP = 1, nP
                ! For each momentum account for dependence
                ! of the diffusion coefficient on momentum
                ! D\propto r_L*v\propto Momentum**2/TotalEnergy
                DInner_I(1:iEnd) = &
                     DInnerInj_I(1:iEnd)*Momentum_I(iP)**2*    &
                     TotalEnergyInj/(TotalEnergy_I(iP)*MomentumInj**2)
                if(UseRealDiffusionUpstream)then
                   where(Radius_I(1:iEnd) > 0.9*Radius_I(iShock))
                      ! upstream: reset the diffusion coefficient to
                      ! (1/6)(0.4AU)*(R/1AU)*v*(pc/1GeV)^(1/3) 
                      DInner_I(1:iEnd) = DInner_I(1:iEnd)/&
                           (Momentum_I(iP)/MomentumInj)**(2.0/3)
                   end where
                end if
                DInner_I(1:iEnd) = max(DInner_I(1:iEnd),&
                     DiffCoeffMin/DOuter_I(1:iEnd))
                call advance_diffusion(Dt, iEnd,&
                     State_VIB(D_,1:iEnd,iBlock), &
                     Distribution_IIB(iP,1:iEnd,iBlock),&
                     DOuter_I(1:iEnd), DInner_I(1:iEnd))
             end do MOMENTUM
          end do STEP
       end do PROGRESS
    end do BLOCK
  contains
    function mach_alfven() result(MachAlfven)
      ! alfvenic mach number for the current line
      real:: MachAlfven
      
      real:: SpeedAlfvenUpstream, SpeedUpstream
      !--------------------------------------------
      ! speed upstream is relative to the shock:
      ! \rho_u * (U_u - V_{shock}) = \rho_d * (U_d - V_{shock})
      SpeedUpstream = Rho_I(iShock+1-nWidth)*&
           (U_I(  iShock + 1 - nWidth) - U_I(  iShock + nWidth))/ &
           (Rho_I(iShock + 1 - nWidth) - Rho_I(iShock + nWidth))
      SpeedAlfvenUpstream = B_I(iShock + nWidth)/ &
           sqrt(cMu*cProtonMass*Rho_I(iShock + nWidth))
      MachAlfven = SpeedUpstream / SpeedAlfvenUpstream
    end function mach_alfven
    !=======================
    subroutine steepen_shock
      ! change the density profile near the shock front so it becomes steeper
      ! for the current line
      integer:: iParticle ! loop variable
      real   :: DLogRhoExcessIntegral, DLogRhoExcess
      real, parameter:: DLogRhoBackground = 0.0
      !---------------------------------------------------------------------
      ! find the excess of DLogRho within the shock compared to background
      ! averaged over length
      DLogRhoExcessIntegral = 0.0
      do iParticle = iShock - nWidth, iShock + nWidth - 1
         DLogRhoExcess = 0.5*(DLogRho_I(iParticle) + DLogRho_I(iParticle+1)) &
              - DLogRhoBackground !D log(rho)/Dt*\Delta t = -\div U*\Delta t
         if(DLogRhoExcess>0)then
            !This is a jump in velocity accross the shock wave * \Delta t
            DLogRhoExcessIntegral = DLogRhoExcessIntegral + &
                 DLogRhoExcess*State_VIB(D_, iParticle, iBlock)
         end if
      end do

      ! check for zero excess
      if(DLogRhoExcessIntegral == 0.0)RETURN
      ! nullify excess  within the smoothed shock 
      DLogRho_I(iShock-nWidth:iShock+nWidth) = min(&
           DLogRhoBackground, &
           DLogRho_I(iShock-nWidth:iShock+nWidth))
      ! ...and concetrate it at the shock front, applying the whole jump
      ! in the velocity at a single grid point 
      DLogRho_I(iShock) = DLogRhoBackground + DLogRhoExcessIntegral/&
           State_VIB(D_, iParticle, iBlock)
      ! also, sharpen the magnetic field magnitude
      ! post shock part
      B_I(iShock+1-nWidth:iShock+1) = maxval(B_I(iShock+1-nWidth:iShock+1))
      ! pre shock part
      B_I(iShock+1:iShock+nWidth  ) = minval(B_I(iShock+1:iShock+nWidth))
    end subroutine steepen_shock
    !=============================================================
    subroutine set_diffusion
      ! set diffusion coefficient for the current line
      !-----------------------------------------------------------
      DOuter_I(1:iEnd) = B_I(1:iEnd)
      !DInner = DiffusionCoeff/B_)
      !\
      ! Compute the diffusion coefficient at the
      ! injection energy:
      !/
      if(.not.UseRealDiffusionUpstream)then
         ! Sokolov et al., 2004: eq (4), 
         ! note: Momentum = TotalEnergy * Vel / C**2
         ! Gyroradius = cGyroRadius * momentum / |B|
         ! DInner = (B/\delta B)**2*Gyroradius*Vel/|B| 
         DInnerInj_I(1:iEnd) = BOverDeltaB2*&
              cGyroRadius*(MomentumInj*cLightSpeed)**2/&
              (B_I(1:iEnd)**2*TotalEnergyInj)/RSun**2
      else
         ! diffusion is different up- and down-stream
         ! Sokolov et al. 2004, paragraphs before and after eq (4)
         where(Radius_I(1:iEnd) > 0.9 * Radius_I(iShock))
            ! upstream: reset the diffusion coefficient to 
            ! (1/6)(0.4AU)*(R/1AU)*v*(pc/1GeV)^(1/3) 
            DInnerInj_I(1:iEnd) = &
                 0.20/3.0 * Radius_I(1:iEnd)/RSun*&
                 (MomentumInj*cLightSpeed**2)/(B_I(1:iEnd)*TotalEnergyInj)*&
                 (MomentumInj*cLightSpeed/energy_in('GeV'))**(1.0/3)
         elsewhere
            ! downstream we use the estimate for (\delta B/B) 
            ! (see detail in Sokolov, 2004):
            ! (\delta B/B)**2 = (\delta B/B)**2_SW*(R/R_SW)
            ! where SW means "shock wave". For the SW we use an
            ! estimate from  from Lee (1983): 
            ! (\delta B/B)**2_SW = const(=10 below)*CInj*MachAlfven
            DInnerInj_I(1:iEnd)=&
                 cGyroRadius*(MomentumInj*cLightSpeed)**2/RSun**2/&
                 (B_I(1:iEnd)**2*TotalEnergyInj)/&
                 (10.0*CInj*MachAlfven)/&
                 min(1.0, 1.0/0.9 * Radius_I(1:iEnd)/Radius_I(iShock))
         end where
      end if

      ! set the left boundary condition for diffusion
      Distribution_IIB(1:nP+1, 1, iBlock) = &
           Distribution_IIB(0, 1, iBlock) * &
           (Momentum_I(0)/Momentum_I(1:nP+1))**SpectralIndex
    end subroutine set_diffusion
    !=========================================================================
    subroutine set_advection_bc
      use SP_ModUnit, ONLY: UnitParticleEnergy, kinetic_energy_to_momentum  
      ! set boundary conditions on grid point on the current line
      ! LOCAL VARIABLES:
      integer:: iParticle   ! loop variable
      real   :: MomentumTi  !Momentum for the thermal energy k_BTi
      !----------------------------------------
      do iParticle = 1, iEnd
         ! injection(Ti, Rho), see Sokolov et al., 2004, eq (3)
         MomentumTi = kinetic_energy_to_momentum(&
              State_VIB(T_,iParticle,iBlock)*UnitParticleEnergy)
         Distribution_IIB(0,iParticle,iBlock) = &
              0.25/cPi/(SpectralIndex-3)*CInj*Rho_I(iParticle)/ &
              MomentumTi**3 * (MomentumTi/MomentumInj)**SpectralIndex
      end do
    end subroutine set_advection_bc
  end subroutine advance
end module SP_ModAdvance
