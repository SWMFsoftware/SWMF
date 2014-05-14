!#NOTPUBLIC  email:xzjia@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModSize
  use ModVarIndexes, ONLY: nVar
  use ModAdvance,    ONLY: Pe_, UseElectronPressure

  use ModUserEmpty,								&
       IMPLEMENTED1 => user_read_inputs,		&
       IMPLEMENTED2 => user_init_session,		&
       IMPLEMENTED3 => user_set_ICs,			&
       IMPLEMENTED4 => user_calc_sources,		&
       IMPLEMENTED5 => user_set_boundary_cells, 	&
       IMPLEMENTED6 => user_set_face_boundary, 		&
       IMPLEMENTED7 => user_set_plot_var

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.2
  character (len=*), parameter :: NameUserModule = 'Europa Single-Fluid MHD Pe, Xianzhe Jia, June 2013'

  logical :: UseResistivePlanet = .false. 

  real, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1) :: NnNeutral, &
       UnxNeutral, UnyNeutral, UnzNeutral, TnNeutral, fioniz_asymfac

  real, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1) :: RhoSdot, &
       CXSdot, CXrate, AlphaSdot, Le
  real :: NeutralMass, Nn1, H1, Nn2, p2, Nr0, rcos, Tmin, mi_mean
  real :: fioniz, CX_sigma, alpha
  real :: RhoSdot_Norm, CXSdot_Norm, AlphaSdot_Norm
  real :: fioniz_Norm, CX_sigma_Norm, alpha_Norm, n_Norm

  real :: SW_Pi, SW_Pe, BodyN_upstream, BodyT_upstream, BodyN_downstream, &
          BodyT_downstream

  character*30 :: type_innerbcs, innerbcs_rho, innerbcs_p, innerbcs_u &
                  , innerbcs_b, Ionization_model


contains

  !===========================================================================

  subroutine user_read_inputs

    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam
    use ModPhysics
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character(len=100) :: NameCommand
    !-------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input Europa starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#EUROPA")
          call read_var('Nn1' , Nn1)             !! Neutral surface density for the exponential falling-off component [1/cm^3]
          call read_var('H1'  , H1)              !! Neutral scale height for the exponential falling-off component  [km]
          call read_var('Nn2', Nn2)              !! Neutral surface density for the power-law falling-off component [1/cm^3]
          call read_var('p2', p2)                !! Index for the power-law falling-off compoent
          call read_var('Nr0', Nr0)              !! Radius of the neutral atmos. base [planet radius]
          call read_var('rcos', rcos)            !! max fraction in cosine relative to uniform distr. [0..infinity]
          call read_var('mi_mean' , mi_mean)     !! mean ion mass [amu]
          call read_var('Tmin', Tmin)            !! Minimum ion temperature (Europa's nightside surface temperature)
          call read_var('fioniz',fioniz)		 !! Ionization frequency (1/s)
          call read_var('CX_sigma',CX_sigma)	 !! Charge-exchange cross-section (cm^2)
          call read_var('alpha',alpha)			 !! recombination rate (cm^3/s)

          H1=H1*1E3                              !! conversion to SI  
          Nn1=Nn1*1E6                            !! conversion to SI
          Nn2=Nn2*1E6                            !! conversion to SI

       case('#IONIZATION')
          call read_var('Ionization_model' , Ionization_model)             !! whether or not use asymmetric IONIZATION

       case('#INNERBCS')
          call read_var('type_innerbcs',type_innerbcs)
          if(type_innerbcs == 'custom') then
            call read_var('innerbcs_rho',innerbcs_rho)
            call read_var('innerbcs_p',innerbcs_p)
            call read_var('innerbcs_u',innerbcs_u)
            call read_var('innerbcs_b',innerbcs_b)
          end if

       case('#BODYPARAM')
          call read_var('BodyN_upstream',BodyN_upstream)        !! Body density (in /cc) on upstream side
          call read_var('BodyT_upstream',BodyT_upstream)            !! Body temperature (in K) on upstream side
          call read_var('BodyN_downstream',BodyN_downstream)    !! Body density (in /cc) on downstream side
          call read_var('BodyT_downstream',BodyT_downstream)        !! Body temperature (in K) on downstream side

          ! BodyN_upstream = BodyN_upstream*1E6           !! conversion to SI (1/m^3)
          ! BodyN_downstream = BodyN_downstream*1E6           !! conversion to SI (1/m^3)

       case('#USERINPUTEND')
          EXIT

       case default
          if(iProc==0) then
             write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) ' *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !===========================================================================

  subroutine user_init_session

    use CON_planet,     ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,    ONLY: cPi
    use ModPhysics,     ONLY: Si2No_V,No2Si_V,UnitRho_, &
         UnitP_, UnitX_, UnitN_, UnitRhoU_, UnitT_
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModProcMH,      ONLY: iProc
    use ModResistivity, ONLY: Si2NoEta
    use ModGeometry,    ONLY: TypeGeometry

    use ModBlockData,   ONLY: MaxBlockData
    use ModSize,        ONLY: nIJK
    use ModVarIndexes,  ONLY: nVar


    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A22,E10.3)"
    !-------------------------------------------------------------------

    MaxBlockData = nVar*nIJK

    !    call set_oktest('user_init_session',DoTest,DoTestMe)

    fioniz_Norm=1/Si2No_V(UnitT_)                           	!! conversion from SI to unitless
    alpha_Norm=1E-6*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)    	!! conversion from SI to unitless
    CX_sigma_Norm=1E-4*Si2No_V(UnitX_)**2      			!! conversion from SI to unitless
    n_Norm=Si2No_V(UnitN_)                             		!! conversion from SI to unitless

    if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,*) 'Conducting Planet (eta =0)'
       write(*,*) ''
       write(*,*) ''
       write(*,*)'Units'
       write(*,*)'No2SI_V(UnitX_)     =',No2SI_V(UnitX_)
       write(*,*)'No2SI_V(UnitRho_)   =',No2SI_V(UnitRho_)
       write(*,*)'No2SI_V(UnitRhoU_)  =',No2SI_V(UnitRhoU_)
       write(*,*)'No2SI_V(UnitP_)     =',No2SI_V(UnitP_)
       write(*,*)'No2SI_V(UnitT_)     =',No2SI_V(UnitT_)
       write(*,*)'fioniz_Norm   =',fioniz_Norm
       write(*,*)'alpha_Norm  =',alpha_Norm
       write(*,*)'CX_sigma_Norm     =',CX_sigma_Norm
       write(*,*)'n_Norm    =',n_Norm
    end if
  end subroutine user_init_session

  !===========================================================================

  subroutine user_set_ICs(iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: Unused_B, nBlockMax
    use ModGeometry,   ONLY: R_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG, nIJK
    use ModVarIndexes
    ! use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, P_, Pe_, Bx_, By_, Bz_
    use ModNumConst,   ONLY: cPi
    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModPhysics


    !    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP, iEnergy, &
    !         iRho, iRhoUx, iRhoUy, iRhoUz

    integer, intent(in) :: iBlock

    integer :: i,j,k
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !-------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

        State_VGB(Rho_,i,j,k,iBlock)      = SW_n*mi_mean
        State_VGB(RhoUx_,i,j,k,iBlock)    = SW_n*mi_mean*SW_Ux
        State_VGB(RhoUy_,i,j,k,iBlock)    = SW_n*mi_mean*SW_Uy
        State_VGB(RhoUz_,i,j,k,iBlock)    = SW_n*mi_mean*SW_Uz
	
        if(UseElectronPressure) then
            if(ElectronPressureRatio.le.0.) call stop_mpi('ERROR: Electron Pressure Ratio > 0 for init!')
            State_VGB(P_,i,j,k,iBlock)      = SW_n*SW_T_dim*Io2No_V(UnitTemperature_)
            State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
        else
            State_VGB(P_,i,j,k,iBlock)      = SW_n*SW_T_dim*Io2No_V(UnitTemperature_) &
                                              *(1.+ElectronPressureRatio)
        end if

    end do; end do; end do

  end subroutine user_set_ics

  !========================================================================
  subroutine user_neutral_atmosphere(iBlock)
    use ModBlockData,  ONLY: get_block_data, set_block_data, put_block_data, &
         use_block_data
    use ModPhysics
    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number, Body1 
    use ModProcMH,     ONLY: iProc 
    use ModGeometry,   ONLY: R_BLK, Xyz_DGB
    !    use ModMultiFluid, ONLY: MassIon_I
    use ModIO,         ONLY: iUnitOut

    integer,intent(in) :: iBlock
    logical :: DoTest, DoTestMe=.true., init=.true.
    real :: theta, cos_theta, Dir_length, term
    real, dimension(3) :: Dir_I
    integer :: i, j, k

    ! Dayside neutral distribution centered at vector direction
    Dir_I(1) =  1.0 !! -SW_Ux ! for direction of undisturbed plasma inflow
    Dir_I(2) =  0.0 !! -SW_Uz
    Dir_I(3) =  0.0 !! -SW_Uz

    !----------------------------------------------------------------------

    !    if(iProc==PROCtest .and. iBlock == BlkTest) then
    !       call set_oktest('user_neutral_atmosphere',DoTest,DoTestMe)
    !    else
    !       DoTest=.false.; DoTestMe=.false.
    !    end if

    if (.not.use_block_data(iBlock)) then
       ! calculate and print out total mass loading (integral from rE to infinity)
       !       if(iProc==0.and.init) then
       !          term = 2*nH2*H2**2*rPlanetSi+nH1*H1*rPlanetSi**2+2*nH2*H2**3+& 
       !               2*nH1*H1**2*rPlanetSi+2*nH1*H1**3+nH2*H2*rPlanetSi**2 
       !          write(iUnitOut,*)'Total Mass Loading = ',&
       !               term*cPi*(4+rcos)*(vO2toOp*MassIon_I(Op_)+vO2toO2p*MassIon_I(O2p_))&
       !               *cProtonMass,' [kg/s]'
       !          init=.false.
       !       end if

       Dir_length = sqrt(dot_product(Dir_I,Dir_I))

       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          ! Initialize neutral atmosphere
          NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)  = 0.
          UnxNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnyNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnzNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          TnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)  = 0.

          ! Initialize the aymmtric ionization multiplication factor
          fioniz_asymfac(i-MinI+1,j-MinJ+1,k-MinK+1) = 1.0

          ! angle of cell position relative to ram direction
          cos_theta=(Dir_I(1)*Xyz_DGB(x_,i,j,k,iBlock)+Dir_I(2)*Xyz_DGB(y_,i,j,k,iBlock)+&
               Dir_I(3)*Xyz_DGB(z_,i,j,k,iBlock))/(R_BLK(i,j,k,iBlock)*Dir_length)

          ! two symmetric distributions w/ different scale heights
          NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = &
               Nn1*exp(-(R_BLK(i,j,k,iBlock)-Nr0)/(H1/rPlanetSi))+ & ! H1 scale length contribution
               Nn2*(Nr0/R_BLK(i,j,k,iBlock))**p2               ! p2 power-law contribution

          !! ??? test w/o body on
          !if (NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) > nH1+nH2) then
          !   NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = nH1+nH2
          !end if

          ! relative increase in relation to symmetric part, cosine distribution
          if(cos_theta>=0.) then
             NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)*(1.+rcos*cos_theta) 
          end if

          ! Calculate asymmetric ionization rate
          ! Power-law decay with radial distance
          fioniz_asymfac(i-MinI+1,j-MinJ+1,k-MinK+1) = 1-0.9*(Nr0/R_BLK(i,j,k,iBlock))**p2  ! p2 power-law contribution

          ! relative increase in relation to symmetric part, cosine distribution
          if(cos_theta<=0.) then
              fioniz_asymfac(i-MinI+1,j-MinJ+1,k-MinK+1) = fioniz_asymfac(i-MinI+1,j-MinJ+1,k-MinK+1)*(1.-rcos*cos_theta) 
          end if

          UnxNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.   !! set bulk flow speed of neutral gas to 0 m/s
          UnyNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnzNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          TnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)  = 600.  !! Kliore et al., Science, 277, 355-358, 1997
       end do;  end do;  end do

       NeutralMass = mi_mean*cProtonMass

       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, NnNeutral)
       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnxNeutral)
       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnyNeutral)
       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnzNeutral)
       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, TnNeutral)
       call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, fioniz_asymfac)

       call set_block_data(iBlock); !! has to be set in case data is accessed before first iteration is finished

       if(DoTestMe) then
          write(*,*)'user_neutral_atmosphere:'
          theta=acos((Dir_I(1)*Xyz_DGB(x_,iTest,jTest,kTest,BlkTest)+Dir_I(2)*Xyz_DGB(y_,iTest,jTest,kTest,BlkTest)+&
               Dir_I(3)*Xyz_DGB(z_,iTest,jTest,kTest,BlkTest))/(R_BLK(iTest,jTest,kTest,BlkTest)*Dir_length))
          write(*,*)'x      = ',Xyz_DGB(x_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'y      = ',Xyz_DGB(y_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'z      = ',Xyz_DGB(z_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'r      = ',R_BLK(iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'theta  = ',theta," [radians]"
          write(*,*)'rcos   = ',rcos," [ ]"
          write(*,*)'Nn1    = ',Nn1," [m^-3]"
          write(*,*)'H1     = ',H1," [m]"
          write(*,*)'Nn2    = ',Nn2," [m^-3]"
          write(*,*)'p2     = ',p2," [R]"

          write(*,*)'Neutral parameters : '
          write(*,*)'N_n    = ',NnNeutral(iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m^-3]"
          write(*,*)'M_n    = ',NeutralMass," [kg]"
          write(*,*)'unx    = ',UnxNeutral(iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
          write(*,*)'uny    = ',UnyNeutral(iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
          write(*,*)'unz    = ',UnzNeutral(iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
          write(*,*)'Tn     = ',TnNeutral(iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [K]"
          write(*,*)''
       end if

    else
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, NnNeutral)
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnxNeutral)
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnyNeutral)
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, UnzNeutral)
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, TnNeutral)
       call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, fioniz_asymfac)
    end if

  end subroutine user_neutral_atmosphere

  !========================================================================

  subroutine user_calc_sources(iBlock)

    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number, Dt_BLK
    use ModAdvance,    ONLY: State_VGB, Source_VC, Rho_, RhoUx_, RhoUy_, RhoUz_, &
         Bx_,By_,Bz_, P_, Energy_
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cElectronCharge, cProtonMass
    use ModGeometry,   ONLY: Rmin_BLK, r_BLK, Xyz_DGB
    use ModCurrent,    ONLY: get_current
    use ModProcMH,     ONLY: iProc
    use ModPhysics
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit, IsPointImplSource
    use ModBlockData,  ONLY: get_block_data, set_block_data, put_block_data, &
                             use_block_data
    integer, intent(in) :: iBlock


    real, dimension(1:3,1:nI,1:nJ,1:nK) :: Current_DC, uIon, uElec
    real, dimension(1:nI,1:nJ,1:nK) :: uIon_sq, uElec_sq, nElec, nIon		
    real, dimension(1:nI,1:nJ,1:nK) :: SRho, SRhoUx, SRhoUy, SRhoUz, SBx, SBy, SBz, &
         SP, SPe, SE 
    real, dimension(1:nI,1:nJ,1:nK) :: Ti, Te

    integer :: i,j,k

    !----------------------------------------------------------------------

!    if(iBlock == BlkTest) then
!       call set_oktest('user_calc_sources',DoTest,DoTestMe)
!    else
!       DoTest=.false.; DoTestMe=.false.
!    end if


    call user_neutral_atmosphere(iBlock)

    !! Set the source arrays for this block to zero
    SRho   = 0.0
    SRhoUx = 0.0
    SRhoUy = 0.0
    SRhoUz = 0.0
    SBx    = 0.0
    SBy    = 0.0
    SBz    = 0.0
    SP     = 0.0
    SPe	   = 0.0
    SE     = 0.0

    ! Initialize the mass-loading terms
    RhoSdot = 0.
    CXrate = 0.
    CXSdot = 0.
    Le = 0.
    AlphaSdot = 0.

    ! Get the electron and ion number densities in SI units ( n_e = n_i )
    do k=1,nK; do j=1,nJ; do i=1,nI
       nIon(i,j,k) = State_VGB(Rho_,i,j,k,iBlock)/mi_mean*No2SI_V(UnitN_)
       nElec(i,j,k) = nIon(i,j,k)
    end do; end do; end do

    ! Get the ion fluid velocity in SI unit
    uIon(1,1:nI,1:nJ,1:nK)=State_VGB(RhoUx_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon(2,1:nI,1:nJ,1:nK)=State_VGB(RhoUy_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon(3,1:nI,1:nJ,1:nK)=State_VGB(RhoUz_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)

    ! Get the ion fluid speed and electron velocty in SI unit	 
    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Calculate ion fluid speed
       uIon_sq(i,j,k) = sum(uIon(1:3,i,j,k)**2)

       ! Calculate electron velocity "uElec" from the Hall velocity -J/(e*n) [m/s]
       call get_current(i,j,k,iBlock,Current_DC(:,i,j,k))
       uElec(1:3,i,j,k) = uIon(1:3,i,j,k)-Current_DC(1:3,i,j,k)/(nElec(i,j,k)*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_)
       uElec_sq(i,j,k) = sum(uElec(1:3,i,j,k)**2)
    end do; end do; end do

    ! Get the ion and electron temperature
    ! Ion temperature is calculated from ion pressure
    ! Electron temperature is calculated from electron pressure
    if (UseElectronPressure) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti(i,j,k) = State_VGB(P_,i,j,k,iBlock)*NO2SI_V(UnitP_)/&
               (cBoltzmann*nIon(i,j,k))
          Te(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/&
               (cBoltzmann*nElec(i,j,k))
       end do; end do; end do
    else
       ! Electron temperature calculated from pressure assuming Te_C=Ti_IC*ElectronTemperatureRatio:
       ! p=nkT with n_e=n_i*Z_i (quasi-neutrality), n=n_e+n_i and p=p_e+p_i=p_i*(1+ElectronPressureRatio)
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti(i,j,k) = State_VGB(P_,i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*nIon(i,j,k))
          Te(i,j,k) = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
               NO2SI_V(UnitP_)/(cBoltzmann*nElec(i,j,k))
       end do; end do; end do
    end if

    ! Calculate the source and loss terms associated with ionizaiton, charge-exchange and recombination
    do k=1,nK; do j=1,nJ; do i=1,nI
       
       RhoSdot(i,j,k) = NnNeutral(i,j,k)*mi_mean*cProtonMass*Si2No_V(UnitRho_) &
                        *fioniz*fioniz_Norm      !! Ionization in dimensionless unit

       if(Ionization_model == 'asymioniz') then
           RhoSdot(i,j,k) = RhoSdot(i,j,k)*fioniz_asymfac(i,j,k)
       end if

       !CXrate(i,j,k) =	CX_sigma*CX_sigma_Norm*NnNeutral(i,j,k)*n_Norm*uIon_sq(i,j,k)**0.5 &
       !  		*Si2No_V(UnitU_)	!! Charge-exchange in dimensionless unit
        
       CXrate(i,j,k) =	CX_sigma*1E-4*NnNeutral(i,j,k)*uIon_sq(i,j,k)**0.5 &
        /Si2No_V(UnitT_)	!! Charge-exchange in dimensionless unit
       CXSdot(i,j,k) = CXrate(i,j,k)*State_VGB(Rho_,i,j,k,iBlocK) 	!! Charge-exchange rate in dimensionless unit
       Le(i,j,k) = alpha*alpha_Norm*nElec(i,j,k)*Si2No_V(UnitN_) 	!! Recombination in dimensionless unit
       AlphaSdot(i,j,k) = Le(i,j,k)*State_VGB(Rho_,i,j,k,iBlocK) 	!! Recombination rate in dimensionless unit

       ! Source and Loss terms in dimensionless unit	
       SRho(i,j,k) = RhoSdot(i,j,k) - AlphaSdot(i,j,k) !! source due to ionization and loss due to recombination

       SRhoUx(i,j,k) = - State_VGB(RhoUx_,i,j,k,iBlocK)*( &
            CXrate(i,j,k) + Le(i,j,k))   !! Momentum loss due to charge exchange & recombination

       SRhoUy(i,j,k) = - State_VGB(RhoUy_,i,j,k,iBlocK)*( &
            CXrate(i,j,k) + Le(i,j,k))   !! Momentum loss due to charge exchange & recombination

       SRhoUz(i,j,k) = - State_VGB(RhoUz_,i,j,k,iBlocK)*( &
            CXrate(i,j,k) + Le(i,j,k))   !! Momentum loss due to charge exchange & recombination

       if(UseElectronPressure) then
          ! Modified pressure terms to take the addition and loss of thermal energy
          SP(i,j,k) = -1*(CXrate(i,j,k) + Le(i,j,k))*State_VGB(P_,i,j,k,iBlocK) &
              + (RhoSdot(i,j,k) + CXSdot(i,j,k))*uIon_sq(i,j,k)*Si2No_V(UnitU_)**2/3

          ! SP(i,j,k) = -1*(CXrate(i,j,k) + Le(i,j,k))* &		!! Ion pressure change due to charge exchange & recombination
          !     State_VGB(P_,i,j,k,iBlocK)
          
          SPe(i,j,k) = -1*Le(i,j,k)*State_VGB(Pe_,i,j,k,iBlocK)		!! Electron pressure change due to recombination
          SE(i,j,k) = -1*(CXrate(i,j,k) + Le(i,j,k))* &		!! Ion pressure change due to charge exchange & recombination
               (0.5*State_VGB(Rho_,i,j,k,iBlocK)*uIon_sq(i,j,k)*Si2No_V(UnitU_)**2 + &
               1.5*State_VGB(P_,i,j,k,iBlocK) + 1.5*State_VGB(Pe_,i,j,k,iBlocK))
       else
          SP(i,j,k) = -1*(CXrate(i,j,k) + Le(i,j,k))*State_VGB(P_,i,j,k,iBlocK) & 
                      *(1.+ElectronPressureRatio) + (RhoSdot(i,j,k) +  & 
                      CXSdot(i,j,k))*uIon_sq(i,j,k)*Si2No_V(UnitU_)**2/3
          ! SP(i,j,k) = -1.5*(CXrate(i,j,k) + Le(i,j,k))* &		!! Ion pressure change due to charge exchange & recombination
          !     State_VGB(P_,i,j,k,iBlocK)*(1.+ElectronPressureRatio)
          SE(i,j,k) = -1*(CXrate(i,j,k) + Le(i,j,k))* &		!! Ion pressure change due to charge exchange & recombination
               (0.5*State_VGB(Rho_,i,j,k,iBlocK)*uIon_sq(i,j,k)*Si2No_V(UnitU_)**2 + &
               1.5*State_VGB(P_,i,j,k,iBlocK)*(1.+ElectronPressureRatio))
       end if

       Source_VC(Rho_   ,i,j,k) = SRho(i,j,k)     + Source_VC(Rho_   ,i,j,k)
       Source_VC(RhoUx_ ,i,j,k) = SRhoUx(i,j,k)   + Source_VC(RhoUx_ ,i,j,k)
       Source_VC(RhoUy_ ,i,j,k) = SRhoUy(i,j,k)   + Source_VC(RhoUy_ ,i,j,k)
       Source_VC(RhoUz_ ,i,j,k) = SRhoUz(i,j,k)   + Source_VC(RhoUz_ ,i,j,k)
       Source_VC(Bx_    ,i,j,k) = SBx(i,j,k)      + Source_VC(Bx_    ,i,j,k)
       Source_VC(By_    ,i,j,k) = SBy(i,j,k)      + Source_VC(By_    ,i,j,k)
       Source_VC(Bz_    ,i,j,k) = SBz(i,j,k)      + Source_VC(Bz_    ,i,j,k)
       Source_VC(Energy_,i,j,k) = SE(i,j,k)       + Source_VC(Energy_,i,j,k)

       if(UseElectronPressure) then
          Source_VC(P_     ,i,j,k) = SP(i,j,k)   + Source_VC(P_     ,i,j,k)
          Source_VC(Pe_    ,i,j,k) = SPe(i,j,k)  + Source_VC(Pe_    ,i,j,k)
       else
          Source_VC(P_     ,i,j,k) = SP(i,j,k)   + Source_VC(P_     ,i,j,k)
       end if

    end do;  end do;  end do

 !   if (.not.use_block_data(iBlock)) then
     !   call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, RhoSdot)
     !   call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, CXSdot)
     !   call put_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, AlphaSdot)
 !   end if

  end subroutine user_calc_sources

  !========================================================================
  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, x1, x2

    implicit none

    integer, intent(in):: iBlock

    character (len=*), parameter :: Name='user_set_boundary_cells'

    !--------------------------------------------------------------------------
    ! For inflow in positive x direction
    ! IsBoundaryCell_GI(:,:,:,ExtraBc_) = Xyz_DGB(x_,:,:,:,iBlock) < x1
    ! For inflow in negative x direction
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = Xyz_DGB(x_,:,:,:,iBlock) > x2

  end subroutine user_set_boundary_cells

  !========================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModVarIndexes
    use ModPhysics,    ONLY: SW_Ux, SW_Uy, SW_Uz, SW_rho, SW_p, SW_T_dim, &
                             BodyNDim_I, BodyRho_I, BodyP_I, FaceState_VI, &
                             Si2No_V, Io2No_V, UnitN_, UnitTemperature_

    use ModFaceBoundary, ONLY: FaceCoords_D, VarsTrueFace_V, iBoundary, &
                               B0Face_D
    use ModB0, ONLY: B0_DX, B0_DY, B0_DZ
    use ModNumConst
    use ModConst,      ONLY: cBoltzmann, cProtonMass 

    real, intent(out):: VarsGhostFace_V(nVar)

    real:: UdotR, BdotR, URefl_D(1:3), Brefl_D(1:3), Borig_D(1:3)
    real:: FaceState_V(nVar)
    real:: XFace,YFace,ZFace, rFace, rFace2, ram_angle

    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)
    rFace2 = XFace**2 + YFace**2 + ZFace**2
    rFace  = sqrt(rFace2)

    ram_angle = acos((-SW_Ux*XFace-SW_Uy*YFace-SW_Uz*ZFace)/rFace/&
    (SW_Ux**2+SW_Uy**2+SW_Uz**2)**0.5)

    !--------------------------------------------------------------------------
    !UdotR = dot_product(VarsTrueFace_V(Ux_:Uz_),FaceCoords_D)* &
    ! 2.0/dot_product(FaceCoords_D,FaceCoords_D)
    !URefl_D = FaceCoords_D*UdotR

    UDotR = sum(VarsTrueFace_V(Ux_:Uz_)*FaceCoords_D)/rFace2

    !B0Face_D = 0.0
    
    Borig_D = VarsTrueFace_V(Bx_:Bz_) ! Reflect B1
    ! Borig_D = Borig_D + B0Face_D  ! Reflect Full B
    BDotR = sum(Borig_D*FaceCoords_D)/rFace2
    Brefl_D = 2*FaceCoords_D*BDotR

    ! Default fixed/initial state for this boundary
    FaceState_V = FaceState_VI(:,iBoundary)

    select case (type_innerbcs)
    case('float') ! Floating boundary
        VarsGhostFace_V = VarsTrueFace_V
        VarsGhostFace_V(Bx_:Bz_)= VarsTrueFace_V(Bx_:Bz_) + B0Face_D
    case('fixed') ! Infinitly conducting body
        VarsGhostFace_V = FaceState_V
        VarsGhostFace_V(Ux_:Uz_)= -VarsTrueFace_V(Ux_:Uz_) ! - 2*UDotR*FaceCoords_D
        VarsGhostFace_V(Bx_:Bz_)= VarsTrueFace_V(Bx_:Bz_) - 2*BDotR*FaceCoords_D
    case('reflect') ! Infinitly conducting body
        VarsGhostFace_V = VarsTrueFace_V
        ! VarsGhostFace_V(Rho_)  = BodyRho_I(1)
        ! VarsGhostFace_V(P_) = BodyP_I(1)
        VarsGhostFace_V(Ux_:Uz_)= -VarsTrueFace_V(Ux_:Uz_)
        VarsGhostFace_V(Bx_:Bz_)= VarsTrueFace_V(Bx_:Bz_) - 2*BDotR*FaceCoords_D
    case('reflectB') ! Reflect B but use float boundary conditions for other parameters        
        ! On upstream, reflect B but use float boundary conditions for other parameters
        VarsGhostFace_V = VarsTrueFace_V
        VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_) - 2*BDotR*FaceCoords_D

        ! On downstream,  set radial velocity to be zero and fix the density and pressure	
        if(ram_angle > cPi/2) then
            VarsGhostFace_V(Ux_:Uz_) = VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D
            VarsGhostFace_V(Rho_)  = BodyRho_I(1)
            VarsGhostFace_V(P_) = BodyP_I(1)
        end if
    case('test') ! Test boundary conditions
        
        VarsGhostFace_V = VarsTrueFace_V    !Floating everything first
        ! VarsGhostFace_V(RhoUx_:RhoUz_) = 0.0

        VarsGhostFace_V(Rho_)  = BodyRho_I(1) ! Fixed density and pressure to body values
        VarsGhostFace_V(P_) = BodyP_I(1)
        ! VarsGhostFace_V(Ux_:Uz_)= VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D ! Reflect azimuthal c  omponents of the velocity
        VarsGhostFace_V(Ux_:Uz_)= -1*VarsTrueFace_V(Ux_:Uz_) ! Reflect all components of the velocity
        ! VarsGhostFace_V(Bx_:Bz_) =  VarsTrueFace_V(Bx_:Bz_) - BRefl_D ! Reflect normal componnet of B
    case('asymfixed') ! Fixed boundary conditions with different density and pressure on upstream/downstream
        VarsGhostFace_V = VarsTrueFace_V    !Floating everything first
        VarsGhostFace_V(Ux_:Uz_)= -1*VarsTrueFace_V(Ux_:Uz_) ! Reflect all components of the velocity
        ! VarsGhostFace_V(Bx_:Bz_)= VarsTrueFace_V(Bx_:Bz_) - 2*BDotR*FaceCoords_D

        ! VarsGhostFace_V(Rho_)  = BodyN_upstream*Io2No_V(UnitN_)*mi_mean
        ! VarsGhostFace_V(P_) = BodyN_upstream*Io2No_V(UnitN_)*BodyT_upstream &
        !                      *Io2No_V(UnitTemperature_)

        ! For downsteam hemisphere, fix the density and pressure	
        !if(ram_angle > cPi/2) then
        if(UDotR > 0) then
            ! VarsGhostFace_V(Ux_:Uz_) = VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D
            VarsGhostFace_V(Rho_)  = BodyN_downstream*Io2No_V(UnitN_)*mi_mean
            VarsGhostFace_V(P_) = BodyN_downstream*Io2No_V(UnitN_)*BodyT_downstream &
                                  *Io2No_V(UnitTemperature_)
        end if
    case('custom') ! user customized inner boundary conditions
        VarsGhostFace_V = VarsTrueFace_V    !Floating everything first

        if(innerbcs_rho == 'fixed') then
            VarsGhostFace_V(Rho_)  = BodyN_upstream*Io2No_V(UnitN_)*mi_mean
        end if

        if(innerbcs_p == 'fixed') then
            VarsGhostFace_V(P_)  = BodyN_upstream*Io2No_V(UnitN_)*BodyT_upstream &
                                     *Io2No_V(UnitTemperature_)
        end if

        if(innerbcs_u == 'zero') then
            VarsGhostFace_V(Ux_:Uz_)= 0.0
        else if(innerbcs_u == 'reflectall') then
            VarsGhostFace_V(Ux_:Uz_)= -1*VarsTrueFace_V(Ux_:Uz_) ! Reflect all components of the velocity
        else if(innerbcs_u == 'reflect') then
            VarsGhostFace_V(Ux_:Uz_)= -1*VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D ! Reflect azimuthal components of the velocity
        end if

        if(innerbcs_b == 'reflect') then
            VarsGhostFace_V(Bx_:Bz_)= VarsTrueFace_V(Bx_:Bz_) - 2*BDotR*FaceCoords_D
        else if(innerbcs_b == 'zero') then
            VarsGhostFace_V(Bx_:Bz_)= 0.0
        end if

        ! For downsteam hemisphere, always fix the density and pressure	
        !if(ram_angle > cPi/2) then
        if(UDotR > 0) then
            ! VarsGhostFace_V(Ux_:Uz_) = VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D
            VarsGhostFace_V(Rho_)  = BodyN_downstream*Io2No_V(UnitN_)*mi_mean
            VarsGhostFace_V(P_) = BodyN_downstream*Io2No_V(UnitN_)*BodyT_downstream &
            *Io2No_V(UnitTemperature_)
        end if

    case('zeroB')
        VarsGhostFace_V(Ux_:Uz_)= VarsTrueFace_V(Ux_:Uz_) - 2*UDotR*FaceCoords_D
        VarsGhostFace_V(Bx_:Bz_)= 0.0
    case('default')
        write(*,*)'unknown type of user inner bcs'
    end select

  end subroutine user_set_face_boundary

 !========================================================================
 subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
    PlotVar_G, PlotVarBody, UsePlotVarBody, &
    NameTecVar, NameTecUnit, NameIdlUnit, IsFound) 


    use ModAdvance,    ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitU_, UnitT_, &
                             UnitRho_, ElectronCharge, ElectronPressureRatio
    use ModVarIndexes, ONLY: Rho_, P_, Pe_
    use ModConst,      ONLY: cBoltzmann, cProtonMass
    use ModCurrent,    ONLY: get_current
!    use ModMultiFluid, ONLY: MassIon_I
    use ModMain,       ONLY: Dt_BLK
    use ModBlockData,  ONLY: get_block_data, set_block_data, put_block_data, &
                             use_block_data
    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k
    real :: nElec, nIon, uIon_sq
    integer, save :: iBlockLast = -1
    real, dimension(3) :: Current_I
    real, dimension(3) :: uIon

    !--------------------------------------------------------------------------

    if(iBlock /= iBlockLast) then
       iBlockLast = iBlock
       call user_neutral_atmosphere(iBlock)
       ! call user_calc_sources(iBlock)
    end if

    IsFound = .true.

    select case(NameVar)
    case('nn')
        NameIdlUnit = '1/cm^3'
        NameTecUnit = '[1/cm^3]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        ! do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
            PlotVar_G(i,j,k) = NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)/1E6
        end do; end do; end do

    case('unx')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[cm/s]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            PlotVar_G(i,j,k) = UnxNeutral(i,j,k)/1E3 !! x direction
        end do; end do; end do

    case('uny')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[cm/s]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            PlotVar_G(i,j,k) = UnyNeutral(i,j,k)/1E3 !! y direction
        end do; end do; end do

    case('unz')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[km/s]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            PlotVar_G(i,j,k) = UnzNeutral(i,j,k)/1E3 !! z direction
        end do; end do; end do

    case('tn')
        NameIdlUnit = 'K'
        NameTecUnit = '[K]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            PlotVar_G(i,j,k) = TnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)
        end do; end do; end do

    case('te')
        NameIdlUnit = 'K'
        NameTecUnit = '[K]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            nIon = State_VGB(Rho_,i,j,k,iBlock)/mi_mean*No2SI_V(UnitN_)
            nElec = nIon
            if(UseElectronPressure)then
                PlotVar_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)/&
                (cBoltzmann*nElec)
            else
                PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)*ElectronPressureRatio/&
                (1.+ElectronPressureRatio)/(cBoltzmann*nElec)
            end if
        end do; end do; end do

    case('ti')
        NameIdlUnit = 'K'
        NameTecUnit = '[K]'
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
            if(UseElectronPressure) then
                PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
                (cBoltzmann*State_VGB(Rho_,i,j,k,iBlock)/mi_mean*NO2SI_V(UnitN_))
            else
                PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(1.+ElectronPressureRatio)/ &
                (cBoltzmann*State_VGB(Rho_,i,j,k,iBlock)/mi_mean*NO2SI_V(UnitN_))
            end if
        end do; end do; end do

    case('uex')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[km/s]'
        do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
            call get_current(i,j,k,iBlock,Current_I)
            ! Get the ion velocity and density in SI unit
            uIon(1)=State_VGB(RhoUx_,i,j,k,iBlock) / &
                    State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitU_)
            nIon = State_VGB(Rho_,i,j,k,iBlock)/mi_mean*No2SI_V(UnitN_)
            nElec = nIon
            PlotVar_G(i,j,k) = (uIon(1)-Current_I(1)/(nElec*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_))/1E3
        end do; end do; end do

    case('uey')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[km/s]'
        do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
            call get_current(i,j,k,iBlock,Current_I)
            ! Get the ion velocity and density in SI unit
            uIon(2)=State_VGB(RhoUy_,i,j,k,iBlock) / &
                    State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitU_)
            nIon = State_VGB(Rho_,i,j,k,iBlock)/mi_mean*No2SI_V(UnitN_)
            nElec = nIon
            PlotVar_G(i,j,k) = (uIon(2)-Current_I(2)/(nElec*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_))/1E3
        end do; end do; end do

    case('uez')
        NameIdlUnit = 'km/s'
        NameTecUnit = '[km/s]'
        do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
            call get_current(i,j,k,iBlock,Current_I)
            ! Get the ion velocity and density in SI unit
            uIon(3)=State_VGB(RhoUz_,i,j,k,iBlock) / &
                    State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitU_)
            nIon = State_VGB(Rho_,i,j,k,iBlock)/mi_mean*No2SI_V(UnitN_)
            nElec = nIon
            PlotVar_G(i,j,k) = (uIon(3)-Current_I(3)/(nElec*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_))/1E3
        end do; end do; end do

    case('rhosi')
        NameIdlUnit = 'amu/cm^3s'
        NameTecUnit = '[amu/cm^3 s]'

        ! call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, RhoSdot)

        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        ! do k=1,nK; do j=1,nJ; do i=1,nI            
           ! PlotVar_G(i,j,k)=RhoSdot(i,j,k)/(n_Norm*mi_mean*fioniz*fioniz_Norm)/1E6 
           ! PlotVar_G(i,j,k)=RhoSdot(i,j,k)*No2SI_V(UnitRho_)/(No2SI_V(UnitT_))/ &
           !                 (cProtonMass*1E6) 
           PlotVar_G(i,j,k)=NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1)*mi_mean* &
                            cProtonMass*Si2No_V(UnitRho_)*fioniz*fioniz_Norm* &
                            No2SI_V(UnitRho_)/(No2SI_V(UnitT_))/(cProtonMass*1E6)

           if(Ionization_model == 'asymioniz') then
                PlotVar_G(i,j,k) = PlotVar_G(i,j,k) &
                                   *fioniz_asymfac(i-MinI+1,j-MinJ+1,k-MinK+1)
           end if
        end do; end do; end do

    case('qsd')
        NameIdlUnit = 'amu/cm^3s'
        NameTecUnit = '[amu/cm^3 s]'
        ! call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, CXSdot)

        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        ! do k=1,nK; do j=1,nJ; do i=1,nI            
           ! PlotVar_G(i,j,k)=CXSdot(i,j,k)/(n_Norm*mi_mean*fioniz*fioniz_Norm)/1E6 
            uIon(1)=State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
            uIon(2)=State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
            uIon(3)=State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
            uIon_sq = sum(uIon(1:3)**2)*No2SI_V(UnitU_)**2  !! convert to SI unit

            PlotVar_G(i,j,k)= CX_sigma*1E-4*NnNeutral(i-MinI+1,j-MinJ+1,k-MinK+1) &
                              *uIon_sq**0.5*State_VGB(Rho_,i,j,k,iBlocK) &
                              *No2SI_V(UnitRho_)/(cProtonMass*1E6)
            
            !PlotVar_G(i,j,k)= CX_sigma*CX_sigma_Norm*NnNeutral(i-MinI+1,j-MinJ+1, &
            !            k-MinK+1)*n_Norm*uIon_sq**0.5* &
            !            State_VGB(Rho_,i,j,k,iBlocK)*No2SI_V(UnitRho_)/ &
            !            (No2SI_V(UnitT_))/(cProtonMass*1E6) 
        end do; end do; end do

    case('alpha')
        NameIdlUnit = 'amu/cm^3s'
        NameTecUnit = '[amu/cm^3 s]'
        ! call get_block_data(iBlock, MaxI-MinI+1, MaxJ-MinJ+1, MaxK-MinK+1, AlphaSdot)
        do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
        ! do k=1,nK; do j=1,nJ; do i=1,nI        
            ! PlotVar_G(i,j,k)=AlphaSdot(i,j,k)/(n_Norm*mi_mean*fioniz*fioniz_Norm)/1E6

            nIon = State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)/(cProtonMass*mi_mean)  !! in SI unit
            nElec = nIon  !! in SI unit
            PlotVar_G(i,j,k)=alpha*alpha_Norm*nElec*Si2No_V(UnitN_)* &
                             State_VGB(Rho_,i,j,k,iBlocK)*No2SI_V(UnitRho_)/ &
                             (No2SI_V(UnitT_))/(cProtonMass*1E6) 
        end do; end do; end do
    case default
        IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

 end subroutine user_set_plot_var

end module ModUser
