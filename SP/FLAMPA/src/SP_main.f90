module SP_ModMain
  use ModConst
  use ModUtilities,ONLY:check_allocate
  implicit none
  integer:: nP=200      !Number of grids along the (ln p)-coordinate.        !
  integer:: nX=100      !Number of points along the spatial coordinate.      !
  real:: CFL=0.9        !Maximum allowed time step.                          !
  real:: PInjection,EInjection= 1.0e+00,EnergyMax=1.0e+07,&
       PMaxLog,DeltaLnP !Injection momentum, injection energy, maximum range !
                        !with respect to (ln p)-coordinate, and step size    !
                        !with respect to (ln p)-coordinate.                  !
  real:: BOverDeltaB2=cOne
                        !Ratio of regular to irregual magnetic field,        !
                        !squared.                                            !
  real:: DiffCoeffMin=cZero
                        !For a given spatial and temporal resolution, the    !
                        !value of the diffusion coefficient should be        !
                        !artificially enhanced in order to get the diffusion !
                        !length larger than the shock front width - governed !
                        !by artificially enhanced numerical dissipation. By  !
                        !doing this, the Fermi acceleration rate decreases   !
                        !accordingly for lower energy particles. By setting  !
                        !DiffCoeffMin=0 (default value), we do NOT enhance   !
                        !the diffusion coefficient! Physically, DiffCoeffMin !
                        !should be given by the product of shock wave speed  !
                        !and local grid spacing.                             !
  logical:: UseRealDiffusionUpstream=.false.
                        !Default valuse is .false. With this logical set to  !
                        !.true., the non-turbulent diffusion is used upstream!
  real:: DsResolution=cHalf/10.0
  real,allocatable,save,dimension(:):: DInner_I,DOuter_I,DInnerInj_I
                        !Diffusion coefficient at momentum p, Laplacian      !
                        !multiplier, and diffusion coefficient at PInjection.!
  real,allocatable,save,dimension(:):: Rho_I,RhoOld_I,RhoSmooth_I,&          !
                                       RhoSmoothOld_I,B_I,FermiA_I           !
                        !Mass density, old value of mass density,            !
                        !and Fermi acceleration rate.                        !
  real,allocatable,save,dimension(:):: T_I
                        !Temperature in units of NameEnergyUnit [KeV]        !
  real,allocatable,save,dimension(:):: U_I
                        !Velocity, in m/s                                    !
  real,allocatable,save,dimension(:,:):: X_DI,State_VI,F_II
                        !Three Lagrangian coordinates, MHD state and         !
                        !full distribution function, f(x,p).                 !
  integer:: iDataSet=0 !Index in file name from IH_ data                    !
  real:: DataInputTime  !Time of simulation from IH_ data                    !
  real:: SP_Time=cZero,SP_Dt
                        !Starting time of simulation and size of time step.  !
  real:: FInjection=0.0120
                        !Value of the distribution function at PInjection,   !
                        !for test simulation, or constant value of f_0 in    !
                        !the distribution for suprathermal particles.        !
  real:: SuprathIndex=4.0+cOne
                        !Spectral index of the suprathermal particles        !
  real,parameter:: CInjection=0.0280
                        !Injection efficiency                                !
  real:: AlfvenMach=cOne!Alfven Mach number of the shock wave                !
  real,parameter:: FTolerance=cTiny**2
                        !Default value for the distribution function.        !
  character(LEN=10):: NameEnergyUnit='KeV',NameParticle='Proton'
                        !Energy unit to be used and type of particles.       !
  integer:: iLogFile    !File unit number for log file.                      !
  logical:: DoLogFile=.true. 
                        !If set to *true* this logical switch enables to     !
                        !save a log file that contains integral fluxes at    !
                        !three energy bands: >1MeV, >10MeV, and >100MeV.     !
  integer:: SP_iPlot=0  !Initialize plot index for SP data.                  !
  character(LEN=50):: SP_TypePlot='cdf def ind tec '
                        !Type of output SP data to be saved in a file.       !
  character(LEN=10):: SP_DirOut='./SP/IO2/'
                        !Output directory for SP plot and log files.         !
  real:: SP_TimePlot=6.0e+01
                        !Time interval to save plot file.                    ! 
  logical:: DoTest=.false.
                        !If set to *true* a test case will be run only.      !
                        !All the complicated physics about the diffusion     !
                        !coefficient will be ignored.                        !
  logical:: UseRelativistic=.true.
                        !Logical switch that determines whether a            !
                        !relativistic or a non-relativistic case will be     !
                        !treated. The default value is *true*.               !
  real,external:: momentum_to_energy,momentum_to_kinetic_energy
  real,external:: energy_to_momentum,kinetic_energy_to_momentum
  real,external:: energy_in
                        !These are functions used to transform momentum into !
                        !(kinetic)energy, and vice versa. The last function  !
                        !is used to set the energy unit for the simulations. !
  integer,save::iShock,iShockOld=1
  integer:: iStdOut=6
  character(LEN=4):: prefix='SP: '
  logical:: DoWriteAll
Contains
  subroutine SP_allocate
    !------------------------------------------------------------------------!
    integer:: iError,iX,iLnP
    character(LEN=*),parameter:: NameSub='SP_allocate'
    !------------------------------------------------------------------------!
    allocate(DInner_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'DInner_I')
    DInner_I = cOne
    allocate(DInnerInj_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'DInnerInj_I')
    DInnerInj_I = cOne
    allocate(DOuter_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'DOuter_I')
    DOuter_I = cOne
    allocate(FermiA_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'FermiA_I')
    FermiA_I = cZero
    allocate(Rho_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'Rho_I')
    Rho_I = cOne
    allocate(RhoOld_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'RhoOld_I')
    RhoOld_I = cOne
    allocate(RhoSmooth_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'RhoSmooth_I')
    RhoSmooth_I = cOne
    allocate(RhoSmoothOld_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'RhoSmoothOld_I')
    RhoSmoothOld_I = cOne
    allocate(T_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'T_I')
    T_I = cOne
    allocate(B_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'B_I')
    B_I = cOne
    allocate(U_I(1:nX),stat=iError)
    call check_allocate(iError,NameSub//'U_I')
    U_I = cZero
    allocate(X_DI(3,1:nX),stat=iError)
    call check_allocate(iError,NameSub//'X_DI')
    allocate(State_VI(8,1:nX),stat=iError)
    call check_allocate(iError,NameSub//'State_VI')
    do iX=1,nX
       X_DI(1  ,iX) = real(iX)
       X_DI(2:3,iX) = cZero
    end do
    allocate(F_II(0:nP+1,1:nX),stat=iError)
    call check_allocate(iError,NameSub//'F_II')
    if(DoTest)then
       F_II(1:nP+1,:) = FTolerance
    else
       do iLnP=1,nP+1
          F_II(iLnP,:) = cTiny/10.0/kinetic_energy_to_momentum(&
               EnergyMax*energy_in(NameEnergyUnit),NameParticle)/&
               (PInjection*exp(iLnP*DeltaLnP))**2
       end do
       !Overall density of the fast particles is of the order of 10^-6 m^-3
       !Integral flux is less than 100 m^-2 ster^-1 s
       !Differential background flux is constant.
       end if
    F_II(0     ,:) = FInjection
  end subroutine SP_allocate
  !------------------------------------ DONE --------------------------------!
  !==========================================================================!
  integer function i_shock(RhoIn_I)
    real,dimension(nX),intent(in)::RhoIn_I
    !\
    ! Find position of shock wave and steepen shock front::
    ! The new position of the shock wave cannot be smaller
    ! than the old position, iShockOld(>=1)!!!
    !/
    i_shock=iShockOld-1+maxloc(&
         log(RhoIn_I(iShockOld:nX-1)/RhoIn_I(iShockOld+1:nX))/      &
         sqrt((X_DI(1,iShockOld+1:nX)-X_DI(1,iShockOld:nX-1))**2+   &
              (X_DI(2,iShockOld+1:nX)-X_DI(2,iShockOld:nX-1))**2+   &
              (X_DI(3,iShockOld+1:nX)-X_DI(3,iShockOld:nX-1))**2),1,&
          MASK=X_DI(1,iShockOld+1:nX)**2+X_DI(2,iShockOld+1:nX)**2+ &
               X_DI(3,iShockOld+1:nX)**2>12.0*RSun**2)
    !   call write_shock(i_shock,RhoIn_I)
  end function i_shock
  subroutine write_shock(i_shock,RhoIn_I)
    integer,intent(in)::i_shock
    real,dimension(nX),intent(in)::RhoIn_I
    !\
    ! Find position of shock wave and steepen shock front::
    !/
    write(iStdOut,*)prefix//'Shock wave front at i_shock = ',i_shock,&
         ', coords = ',X_DI(:,i_shock),', r = ',sqrt(sum(X_DI(:,i_shock)**2))/&
         Rsun,' Rsun',', density = ',RhoIn_I(i_shock),', compr. ratio = ',&
         maxval(RhoIn_I(i_shock+1-nint(cOne/DsResolution):iShock))/&
         minval(RhoIn_I(i_shock+1:i_shock+nint(cOne/DsResolution))),&
         ', Mach number = ',AlfvenMach
  end subroutine write_shock

!============================================================================!
!============================================================================!
subroutine SP_diffusive_shock(&
     TypeAction,              &                  
     TimeToFinish,            &                
     nToFinish)
  implicit none
  !--------------------------------------------------------------------------!
  integer:: nStep=1,iStep,iX,iLnP !Number of time steps, current index of    !
                                  !time step, current index of x-coordinate, !
                                  !and current index of (ln p)-coordinate.   !
  integer:: iFile                 !File unit to be assigned.                 !
  real:: Energy,Momentum          !Energy and momentum of particles.         !
  integer::nProgress,iProgress    !To split time interval between snapshots  !
  real::DtProgress                !in order for the shock wave to propagate  !
                                  !not more than one grid cell per time step !  
  save
  character(LEN=*),intent(in):: TypeAction
  real,intent(in),optional   :: TimeToFinish  
  integer,intent(in),optional:: nToFinish
  !--------------------------------------------------------------------------!
  write(iStdOut,*)prefix//TypeAction
  select case(TypeAction)
  !---------------------------------- INIT ----------------------------------!
  case("INIT")
    
     write(iStdOut,*)prefix,'EInjection = ',EInjection,' '//NameEnergyUnit
     EInjection = EInjection*&
         energy_in(NameEnergyUnit)!Convert injection energy in MKS units.    !
     write(iStdOut,*)prefix,'EInjection = ',EInjection,' J'
     PInjection = &
         kinetic_energy_to_momentum(EInjection,&
         NameParticle)            !Compute injection momentum in MKS units.  !
     write(iStdOut,*)prefix,'PInjection = ',PInjection,' Kg*m/s'
     Energy = EnergyMax*&
         energy_in(NameEnergyUnit)!Convert maximum energy in MKS units.      !
     write(iStdOut,*)prefix,'EnergyMax = ',EnergyMax,' '//NameEnergyUnit,&
          ' = ',Energy, ' J'
     Momentum = &
         kinetic_energy_to_momentum(Energy,&
         NameParticle)            !Compute maximum momentum in MKS units.    !
     write(iStdOut,*)prefix,'PMax = ',Momentum,'Kg*m/s'
     PMaxLog = &
         log(Momentum/PInjection) !Compute the maximum range with respect to !
                                  !the (ln p)-coordinate.                    !
     DeltaLnP = PMaxLog/nP        !Step-size with respect to the             !
                                  !(ln p)-coordinate.                        !
     
     write(iStdOut,*)prefix,' Log(Momentum) interval', PMaxLog,&
          ' is divided by ',nP,' equaly spaced intervals, DeltaLnP = ',&
          DeltaLnP
     if (DoLogFile) call write_logfile_SP('OPEN')
  !----------------------------------- RUN ----------------------------------!
  case("RUN")
     !-----------------------------------------------------------------------!
     ! Advance the solution of the diffusive kinetic equation:               !
     !            f_t+[(1/3)*(d(ln rho)/dt]*f_{ln p}=B*d/ds[D/B*df/ds]       !
     ! from strat time, SP_Time, to end time, TimeToFinish.                  !
     !-----------------------------------------------------------------------!
     nProgress = max(iShock-iShockOld,1)
     iShockOld = min(iShockOld,iShock-1)
     DtProgress = (TimeToFinish-SP_Time)/real(nProgress)
     do iProgress=1,nProgress
        Rho_I(1:nX)=RhoSmoothOld_I(1:nX)+&
             (RhoSmooth_I(1:nX)-RhoSmoothOld_I(1:nX))*&
             real(iProgress)/real(nProgress)
        !\
        ! Steepen shock wave front +/- (cOne/DsResolution) points from
        ! position of shock front::
        !/
        iShock=iShockOld+iProgress
        if(iShock>nint(cOne/DsResolution).and.&
             iShock<nX-nint(cOne/DsResolution))then
           !\
           ! Obtain the Alfven Mach number of the shock wave::
           !/
           AlfvenMach=Rho_I(iShock+1-nint(cOne/DsResolution))*&
                (U_I(iShock+1-nint(cOne/DsResolution))-  &
                 U_I(iShock+nint(cOne/DsResolution)))/   &
                (Rho_I(iShock+1-nint(cOne/DsResolution))-&
                 Rho_I(iShock+nint(cOne/DsResolution)))/ &
                (B_I(iShock+nint(cOne/DsResolution))/    &
                 sqrt(cMu*cProtonMass*Rho_I(iShock+nint(cOne/DsResolution))))
           call write_shock(iShock,Rho_I(1:nX))
           Rho_I(iShock+1-nint(cOne/DsResolution):iShock)=&
                maxval(Rho_I(iShock+1-nint(cOne/DsResolution):iShock))
           Rho_I(iShock+1:iShock+nint(cOne/DsResolution))=&
                minval(Rho_I(iShock+1:iShock+nint(cOne/DsResolution)))
        end if
        !\
        ! Obtain the Fermi acceleration rate divided by
        ! the time step, and then divided by DeltaLnP: 
        !/
        FermiA_I(1:nX) = log(Rho_I(1:nX)/RhoOld_I(1:nX))*((1.0/3)/DeltaLnP)
        RhoOld_I = Rho_I

        !\
        ! Check if the number of time steps is greater than 1:
        !/
        nStep = 1+int(maxval(abs(FermiA_I(1:nX)))/CFL)
        SP_Dt = DtProgress/real(nStep)
        write(iStdOut,*)prefix,' Time step is set to ',SP_Dt,' s'
        !\
        ! If nStep>1, compute the Fermi acceleration rate
        ! divided by the time step, and then divided by DeltaLnP:
        !/
        if (nStep>1) FermiA_I = FermiA_I/real(nStep)
        !\
        ! Approximate the diffusion coefficient and define
        ! the Laplacian multiplier:
        !/
        DOuter_I(1:nX) = B_I(1:nX)
        !\
        ! Compute the diffusion coefficient at the
        ! injection energy:
        !/
        if (DoTest) then
           DInnerInj_I(1:nX) = BOverDeltaB2   !Value for test purposes only.
        else
           Energy = momentum_to_energy(PInjection,NameParticle)
           DInnerInj_I(1:nX) = BOverDeltaB2*& !Physically meaningful value.
                cGyroRadius*(PInjection*cLightSpeed)**2/&
                (B_I(1:nX)**2*Energy)
           !           write(iStdOut,*)prefix//&
           !                'Min and Max value for diffusuon coefficient at injection energy = ',&
           !                minval(DInnerInj_I(1:nX)*DOuter_I(1:nX)),&
           !                maxval(DInnerInj_I(1:nX)*DOuter_I(1:nX))
           if(UseRealDiffusionUpstream)then
              ! reset the diffusion coefficient upstream to be equal to
              !(1/6)(0.4AU)*(R/1AU)*v*(pc/1GeV)^(1/3) 
              where(X_DI(1,1:nX)**2+&
                    X_DI(2,1:nX)**2+&
                    X_DI(3,1:nX)**2>&
                   (X_DI(1,iShock)**2+&
                    X_DI(2,iShock)**2+&
                    X_DI(3,iShock)**2)*1.2)
              DInnerInj_I(1:nX)=&
                   sqrt(X_DI(1,1:nX)**2+&
                        X_DI(2,1:nX)**2+&
                        X_DI(3,1:nX)**2)*&
                   6.6667e-2*&
                   (PInjection*cLightSpeed**2)/&
                   (B_I(1:nX)*Energy)*&
                   (PInjection*cLightSpeed/energy_in('GeV'))**(1.0/3)
              elsewhere
                 DInnerInj_I(1:nX)=DInnerInj_I(1:nX)/min(cOne,&
                      sqrt((X_DI(1,1:nX)**2+&
                            X_DI(2,1:nX)**2+&
                            X_DI(3,1:nX)**2)/&
                           (X_DI(1,iShock)**2+&
                            X_DI(2,iShock)**2+&
                            X_DI(3,iShock)**2))/0.9)/&
                           (BOverDeltaB2*10.0*CInjection*AlfvenMach)
              end where
           end if
        end if
        do iStep=1,nStep
           do iX=1,nX
              !\
              ! Apply the boundary condition for F_II at injection momentum::
              ! We assume distribution function for suprathermal particles
              ! given by: 
              ! f(p_thermal<p<p_injection,x)=CInjection*(1/4pi/(SuprathIndex-3))*&
              !                   f_injection*(Rho(x)/(2TM_i)^(3/2))*&
              !                   ((2TM_i)^(1/2)/p_injection)^SuprathIndex
              ! where f_injection is a dimensionless constant of the order of
              ! unity, equal to 1 by default, and the default spectral index
              ! for suprathermal particles is SuprathIndex=5.
              !/
              if (.not.DoTest) F_II(0,iX)=FInjection*CInjection*&
                   0.250/cPi/(SuprathIndex-3.0)*Rho_I(iX)/&
                   kinetic_energy_to_momentum(&
                   T_I(iX)*energy_in(NameEnergyUnit),NameParticle)**3*&
                   (kinetic_energy_to_momentum(&
                   T_I(iX)*energy_in(NameEnergyUnit),NameParticle)    &
                   /PInjection)**SuprathIndex
              !\
              ! Advance the advection part of the DKE:
              !/
              call advance_advection(FermiA_I(iX),nP,F_II(:,iX))
           end do
           !           write(iStdOut,*)prefix,' '
           !           write(iStdOut,*)prefix,'Min/Max of F_II at PInjection = ',&
           !           minval(F_II(0,1:nX)),maxval(F_II(0,1:nX)), 's^3 kg^-^3 m^-^6'
           !           write(iStdOut,*)prefix,' '
           
           do iLnP=1,nP
              Momentum = exp(real(iLnP)*DeltaLnP) !This gives P/PInjection.
              if (.not.UseRelativistic) then
                 DInner_I = DInnerInj_I*Momentum**2
              else
                 DInner_I = DInnerInj_I*Momentum**2 *&
                      momentum_to_energy(cOne*    PInjection,NameParticle)/&
                      momentum_to_energy(Momentum*PInjection,NameParticle)
              end if
              if(UseRealDiffusionUpstream)&
              ! reset the diffusion coefficient upstream to be equal to
              !(1/6)(0.4AU)*(R/1AU)*v*(pc/1GeV)^(1/3) 
              where(X_DI(1,iShock+1:nX)**2+&
                    X_DI(2,iShock+1:nX)**2+&
                    X_DI(3,iShock+1:nX)**2>&
                   (X_DI(1,iShock)**2+&
                    X_DI(2,iShock)**2+&
                    X_DI(3,iShock)**2)*1.2)&
                    DInner_I(iShock+1:nX)=&
                    DInner_I(iShock+1:nX)*Momentum**((1.0/3)-cOne)
              !\
              ! For a given spatial and temporal resolution, the value of
              ! the diffusion coefficient should be artificially enhanced
              ! in order to get the diffusion length larger than the shock
              ! front width - governed by artificially enhanced numerical
              ! dissipation. By doing this, the Fermi acceleration rate
              ! decreases accordingly for lower particle energies.
              !/
              call enhance_diffusion(DInner_I,DOuter_I)
              !\
              ! Advance the diffusion part of the DKE:
              !/
              call advance_diffusion(SP_Dt,nX,X_DI(:,1:nX),F_II(iLnP,1:nX),&
                   DOuter_I(1:nX),DInner_I(1:nX))
           end do
           !\
           ! Update simulated time for SP:
           !/
           SP_Time = SP_Time+SP_Dt
           if (DoLogFile) call write_logfile_SP('WRITE')
           call write_plotfile_SP(&
                int(SP_Time/SP_TimePlot)/=SP_iPlot,SP_TypePlot)
        end do  !iStep loop
     end do     !iProgress loop
  !--------------------------------- FINALIZE -------------------------------!
  case("FINALIZE")
     if (DoLogFile) call write_logfile_SP('CLOSE')
     call write_plotfile_SP(.true.,SP_TypePlot)
  case default
     call CON_stop('Unknown action in SP_difussive_shock'//TypeAction)
  end select
  !------------------------------------ DONE --------------------------------!
end subroutine SP_diffusive_shock
end module SP_ModMain
