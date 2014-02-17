!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:rubinmar@umich.edu  expires:12/31/2099
!========================================================================

module ModUser

  use ModMultiFluid
  use ModSize,       ONLY: nI, nJ, nK, nBlk
  use ModVarIndexes, ONLY: nVar
  use ModAdvance,    ONLY: Pe_, UseElectronPressure
  use ModUserEmpty,                               &
       IMPLEMENTED1  => user_read_inputs,         &
       IMPLEMENTED2  => user_calc_sources,        &
       IMPLEMENTED3  => user_update_states,       &
       IMPLEMENTED4  => user_set_face_boundary,   &
       IMPLEMENTED5  => user_set_resistivity,     &
       IMPLEMENTED6  => user_material_properties, &
       IMPLEMENTED7  => user_init_point_implicit, &
       IMPLEMENTED8  => user_init_session,        &
       IMPLEMENTED9  => user_set_plot_var,        &
       IMPLEMENTED10 => user_set_ICs,             &
       IMPLEMENTED11 => user_get_log_var,         & 
       IMPLEMENTED12 => user_set_boundary_cells,  &
       IMPLEMENTED13 => user_set_cell_boundary

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.1
  character (len=*), parameter :: NameUserModule = &
       'Rubin, 3-fluid Comet MHD module, Dec 2012'

  integer, parameter, public :: nNeutral = 2 !! number of neutral species
  !! Neutral species names
  character (len=6), parameter, public :: NameNeutral_I(nNeutral) = &
       (/ 'H2O  ', 'H    ' /)
  integer, parameter :: H2O_  =  1
  integer, parameter :: H_    =  2
  !! Ion species names
  integer, parameter :: SW_   =  1
  integer, parameter :: Hp_   =  2
  integer, parameter :: H2Op_ =  3

  real, dimension(nNeutral,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1) :: NnNeutral_IG, &
       UnxNeutral_IG, UnyNeutral_IG, UnzNeutral_IG, TnNeutral_IG
  real, dimension(1:nNeutral) :: NeutralMass_I(nNeutral)
  integer :: iNeutralBlockLast = -1
  real :: Qprod, Tmin, rHelio, vHI
  real, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: ne20eV_GB = 0.

  !! Plotting array to be used for testing
  ! real, dimension(4,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: TestArray = 0.


  !! To plot the cells where the heat flux is capped by the limiter
  ! real, public, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: FluxLimited_GB = 0.
  !! Add the following lines to get_impl_heat_cond_state 
  ! use ModUser, ONLY: FluxLimited_GB
  !! and after "The threshold heat flux limiter model"
  ! FluxLimited_GB(i,j,k,iBlock) = 0.
  ! if(HeatCoef*GradTe > FreeStreamFlux) &
  !    FluxLimited_GB(i,j,k,iBlock) = 1.
  !! before limiter is applied


contains

  !========================================================================

  subroutine user_read_inputs

    use ModMain,      ONLY: lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var, read_line, read_command
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand

    !-----------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input Comet starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#COMET")
          call read_var('Qprod', Qprod)             !! Neutral gas production rate [1/s]
          call read_var('rHelio', rHelio)           !! Heliocentric distance [AU]
          call read_var('vHI', vHI)                 !! Ionization frequency for cometary heavy ions (1/lifetime of cometary heavy neutrals)
          call read_var('Tmin', Tmin)               !! Minimum ion temperature (enforced in update states)
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input Comet ends'
          endif
          EXIT
       case default
          if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !========================================================================

  subroutine user_neutral_atmosphere(iBlock)
    use ModBlockData,  ONLY: get_block_data, set_block_data, put_block_data, &
         use_block_data
    use ModPhysics,    ONLY: cPi, rPlanetSi, cProtonMass, No2SI_V, &
         cRadToDeg, UnitX_
    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number 
    use ModProcMH,     ONLY: iProc 
    use ModGeometry,   ONLY: R_BLK, Xyz_DGB,x1,y1,y2,z1,z2
    use ModMultiFluid, ONLY: MassIon_I
    use ModIO,         ONLY: iUnitOut
    use ModConst,      ONLY: mSun, cGravitation
    use ModLookupTable,ONLY: i_lookup_table, interpolate_lookup_table

    integer,intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true., init=.true.
    integer :: i, j, k, iNeutral
    ! real :: RadAccl, uH
    real :: xUp
    real, dimension(nNeutral,nNeutral) :: DissocRate_II
    real, dimension(nNeutral) :: DestructRate_I, uNeutr_I

    ! Variables for neutral tables
    logical:: DoCheckTable = .true.
    integer, save:: iTableNeutral_I(nNeutral)
    real:: r, Phi, x, y, z, Yz, Neutral_V(6)
    integer, parameter:: Nn_=1, Uxn_ = 2, Uyn_ = 3, Tn_=4

    character(len=*), parameter:: NameSub = 'user_neutral_atmosphere'
    !----------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest) then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !! Neutral dissociation/destruction rates
    DissocRate_II = 0.0 ; DestructRate_I = 0.0 
    DestructRate_I(H2O_) = vHI  !! from PARAM.in file

    iNeutralBlockLast = iBlock

    if(DoCheckTable)then
       ! Check for neutral data tables
       do iNeutral = 1, nNeutral
          iTableNeutral_I(iNeutral) = i_lookup_table(NameNeutral_I(iNeutral))
       end do

       if(iProc==0)write(*,*) NameSub,' iTableNeutral_I=', iTableNeutral_I

       DoCheckTable = .false.
    end if

    if (.not.use_block_data(iBlock)) then

       ! Fill analytic solution if any neutral data table is missing
       if(any(iTableNeutral_I < 0))then

          !if(iProc==0.and.init) then
          !   !! Init, placeholder
          !   init=.false.
          !end if

          uNeutr_I(H2O_) = 1000. ! [m/s]
          uNeutr_I(H_)   = 8000. ! [m/s]

          ! Average hydrogen atom velocity (18 km/s)  
          ! uH = 18.E3
          ! Hydrogen radiation pressure acceleration (beta = 1.0)
          ! RadAccl = 1.0*cGravitation*mSun/(rHelio**2*cAU**2)

          xUp = 4.0 ! upstream limit for hydrogen distribution in Mkm

          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

             ! Haser model for H2O
             NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1) = Qprod/(4.*cPi*(R_BLK(i,j,k,iBlock)*rPlanetSI)**2&
                  *uNeutr_I(H2O_))*exp(-DestructRate_I(H2O_)/uNeutr_I(H2O_)*R_BLK(i,j,k,iBlock)*rPlanetSI)

             ! ! Limits hydrogen upstream distance due to radiation pressure
             ! ! (parabola x<-0.5*RadAccl/uH^2*(y^2+z^2)+0.5*uH^2/RadAccl)
             ! if(Xyz_DGB(x_,i,j,k,iBlock) < -0.5*RadAccl/(uH**2)*rPlanetSI*&
             !      (Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)+0.5*uH**2/RadAccl/rPlanetSI) then
             !    ! H (after Combi 1996 with [cm^-3] -> [m^-3] and r[m] -> r[cm])
             !    NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1) = Qprod*1.205E-11*1e6*(1E2*R_BLK(i,j,k,iBlock)*rPlanetSI)**(-1.6103)
             ! else
             !    NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
             ! end if

             ! Limits hydrogen upstream distance by parabola to avoid interaction with the side boundary
             ! (parabola x<(x_min-x_upstream)/min_width^2*(y^2+z^2)+x_upstream
             if(Xyz_DGB(x_,i,j,k,iBlock) < (x1-xUp)/min(abs(y1),abs(y2),abs(z1),abs(z2))**2.*&
                  (Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)+xUp) then
                ! H (after Combi 1996 with [cm^-3] -> [m^-3] and r[m] -> r[cm])
                NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1) = Qprod*1.205E-11*1e6*(1E2*R_BLK(i,j,k,iBlock)*rPlanetSI)**(-1.6103)
             else
                NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
             end if

             ! H from Haser Model (assuming vn_H2O = vn_H)
             !   DissocRate_II(H2O_,H_) = 1.64E-5/(rHelio**2) ! H2O + hv -> OH + H (Huebner et al. 1992)
             !NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)    = Qprod*DissocRate_II(H2O_,H_)/ &
             !     (DestructRate_I(H2O_)-DissocRate_II(H2O_,H_))* &
             !     (exp(-DissocRate_II(H2O_,H_)/uNeutr_I(H2O_)*R_BLK(i,j,k,iBlock)*rPlanetSI)-&
             !     exp(-DestructRate_I(H2O_)/uNeutr_I(H2O_)*R_BLK(i,j,k,iBlock)*rPlanetSI))/&
             !     (4.*cPi*(R_BLK(i,j,k,iBlock)*rPlanetSI)**2*uNeutr_I(H2O_))     

             UnxNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1) = uNeutr_I(H2O_)*Xyz_DGB(x_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)
             UnyNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1) = uNeutr_I(H2O_)*Xyz_DGB(y_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)
             UnzNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1) = uNeutr_I(H2O_)*Xyz_DGB(z_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)

             UnxNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)   = uNeutr_I(H_)*Xyz_DGB(x_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)
             UnyNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)   = uNeutr_I(H_)*Xyz_DGB(y_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)
             UnzNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)   = uNeutr_I(H_)*Xyz_DGB(z_,i,j,k,iBlock)/ &
                  R_BLK(i,j,k,iBlock)


             TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)  =  100.   !! estimate
             TnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)    =  1000.  !! estimate
          end do;  end do;  end do
       end if

       ! Fill neutral state by interpolating lookup table
       do iNeutral = 1, nNeutral

          if(iTableNeutral_I(iNeutral) <= 0) CYCLE

          do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
             x = Xyz_DGB(1,i,j,k,iBlock)
             y = Xyz_DGB(2,i,j,k,iBlock)
             z = Xyz_DGB(3,i,j,k,iBlock)

             Yz = sqrt(y**2 + z**2)

             ! Radial distance in meters for lookup
             r  = sqrt(x**2 + y**2 + z**2)*No2Si_V(UnitX_)

             ! Angle relative to x axis in degrees for lookup
             Phi  = atan2(Yz, x)*cRadToDeg

             call interpolate_lookup_table( &
                  iTableNeutral_I(iNeutral), r, Phi, Neutral_V)

             NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1) = Neutral_V(Nn_)
             UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)= Neutral_V(Uxn_)
             UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)= Neutral_V(Uyn_)*y/Yz
             UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)= Neutral_V(Uyn_)*z/Yz
             TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1) = Neutral_V(Tn_)
          end do; end do; end do
       end do


       NeutralMass_I(H2O_) = 18.*cProtonMass
       NeutralMass_I(H_)   =  1.*cProtonMass

       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, NnNeutral_IG)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnxNeutral_IG)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnyNeutral_IG)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnzNeutral_IG)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, TnNeutral_IG)

       call set_block_data(iBlock); !! has to be set in case data is accessed before first iteration is finished

       if(DoTestMe) then
          write(*,*)'user_neutral_atmosphere:'
          write(*,*)'x      = ',Xyz_DGB(x_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'y      = ',Xyz_DGB(y_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'z      = ',Xyz_DGB(z_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'r      = ',R_BLK(iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'Qprod  = ',Qprod," [1/s]"
          do iNeutral=1,nNeutral
             write(*,*)'Neutral species # ',iNeutral,': ', NameNeutral_I(iNeutral)
             write(*,*)'n_n    = ',NnNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m^-3]"
             write(*,*)'m_n    = ',NeutralMass_I(iNeutral)," [kg]"
             write(*,*)'unx    = ',UnxNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'uny    = ',UnyNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'unz    = ',UnzNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'Tn     = ',TnNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [K]"
          end do
          write(*,*)''
       end if
    else
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, NnNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnxNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnyNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnzNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, TnNeutral_IG)
    end if

  end subroutine user_neutral_atmosphere

  !========================================================================

  subroutine calc_electron_collision_rates(Te,nElec,i,j,k,iBlock,fen_I,fei_I)

    ! calculate all collision rates for electrons (fen, fei)
    ! (used for sources & resistivity)

    use ModAdvance,    ONLY: State_VGB, Rho_
    use ModPhysics,    ONLY: No2SI_V, UnitN_
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I

    integer,intent(in) :: i,j,k,iBlock   
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(out)   :: fen_I(nNeutral)
    real,intent(out)   :: fei_I(nIonFluid)

    real :: sqrtTe
    !----------------------------------------------------------------------

    !! electron-neutral and electron-ion collision rates
    !! provide all rates in SI units

    !! reduced temperature ~ Te 
    sqrtTe = sqrt(Te)

    !! initialize all collision rates with zero
    fei_I = 0. ; fen_I = 0.

    !! Electron - neutral collision rates
    !! e - H2O,  Itikawa, Planet. Space Sci., 1971 and Itikawa, Phys. Fluids 1983
    fen_I(H2O_) = 2.745E-5*NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)/1e6*Te**(-0.62)      !! rate in [1/s]
    !! e - H, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fen_I(H_) = 4.5E-9*NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)/1e6*(1.-1.35E-4*Te)*sqrtTe !! rate in [1/s]
    if (fen_I(H_) < 0.) then
       fen_I(H_) = 0.
    end if

    !!Electron - ion collision rates
    !! e - H2Op, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fei_I(H2Op_) = 54.5*ChargeIon_I(H2Op_)**2*State_VGB(H2OpRho_,i,j,k,iBlock)/MassIon_I(H2Op_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)      !! rate in [1/s]
    !! e - Hp, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fei_I(Hp_) = 54.5*ChargeIon_I(Hp_)**2*State_VGB(HpRho_,i,j,k,iBlock)/MassIon_I(Hp_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)      !! rate in [1/s]
    fei_I(SW_) = 54.5*ChargeIon_I(SW_)**2*State_VGB(SwRho_,i,j,k,iBlock)/MassIon_I(SW_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)      !! rate in [1/s]

  end subroutine calc_electron_collision_rates

  !========================================================================
  subroutine user_calc_rates(Ti_I,Te,i,j,k,iBlock,nElec,nIon_I,fin_II,fii_II,fie_I,alpha_I,kin_IIII,v_II,&
       ve_II,uElec_D,uIon_DI,Qexc_II,Qion_II)

    ! calculate all rates not involving electron collisions

    use ModPhysics,  ONLY: SI2No_V, UnitN_, rPlanetSI, rBody
    use ModConst,    ONLY: cElectronCharge, cBoltzmann, cElectronMass, cProtonMass
    use ModMain,     ONLY: Body1, iTest, jTest, kTest, BlkTest
    use ModNumConst, ONLY: cPi
    use ModGeometry, ONLY: R_BLK, Xyz_DGB

    integer,intent(in) :: i,j,k,iBlock
    real,intent(in)    :: Ti_I(nIonFluid)
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(in)    :: nIon_I(nIonFluid)
    real,intent(in)    :: uElec_D(3)
    real,intent(in)    :: uIon_DI(3,nIonFluid)
    real,intent(out)   :: fin_II(nIonFluid,nNeutral)
    real,intent(out)   :: fii_II(nIonFluid,nIonFluid)
    real,intent(out)   :: fie_I(nIonFluid)
    real,intent(out)   :: alpha_I(nIonFluid)
    real,intent(out)   :: kin_IIII(nIonFluid,nNeutral,nNeutral,nIonFluid)
    real,intent(out)   :: v_II(nNeutral,nIonFluid)
    real,intent(out)   :: ve_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qexc_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qion_II(nNeutral,nIonFluid)


    real :: Tr, ueBulk2, ueTherm2, Ee, A(7), Tred, Mred
    real :: Dist, NCol, sigma, J3, uNeutr, log10Te, sqrtTe
    real,dimension(nNeutral,nIonFluid) :: sigma_e
    integer :: n
    real, save :: ElImpRate_I(nNeutral,61)!, ElCrossSect_I(61)


    real, save :: sigmaeh2o = 4.53E-21  !! Ionization cross section for 20 eV electrons [m^2]
    real, save :: ve = 2.65E6           !! Speed of 20 eV electrons [m/s]

    ! H2O and H electron impact rate depending on electron temperature (Cravens et al. 1987)
    ElImpRate_I(H2O_,1:61) = (/ 0.00E+00, 1.14E-16, 2.03E-16, 3.04E-16, 4.37E-16, 6.34E-16, 9.07E-16, &
         1.28E-15, 1.79E-15, 2.34E-15, 3.15E-15, 4.35E-15, 5.54E-15, 6.90E-15, 8.47E-15, 1.05E-14, 1.25E-14, &
         1.51E-14, 1.80E-14, 2.09E-14, 2.41E-14, 2.74E-14, 3.09E-14, 3.46E-14, 3.90E-14, 4.34E-14, 4.80E-14, &
         5.24E-14, 5.70E-14, 6.15E-14, 6.63E-14, 7.15E-14, 7.61E-14, 8.03E-14, 8.47E-14, 8.90E-14, 9.30E-14, &
         9.72E-14, 1.01E-13, 1.05E-13, 1.08E-13, 1.11E-13, 1.14E-13, 1.16E-13, 1.18E-13, 1.21E-13, 1.23E-13, &
         1.24E-13, 1.25E-13, 1.25E-13, 1.25E-13, 1.23E-13, 1.23E-13, 1.23E-13, 1.21E-13, 1.20E-13, 1.17E-13, &
         1.15E-13, 1.12E-13, 1.10E-13, 1.07E-13 /)
    ElImpRate_I(H_,1:61) = (/ 0.00E+00, 1.74E-16, 2.89E-16, 4.14E-16, 6.01E-16, 8.18E-16, 1.02E-15, 1.32E-15, &
         1.72E-15, 2.20E-15, 2.84E-15, 3.48E-15, 4.20E-15, 5.12E-15, 5.97E-15, 6.83E-15, 7.82E-15, 8.95E-15, &
         1.01E-14, 1.12E-14, 1.24E-14, 1.36E-14, 1.48E-14, 1.61E-14, 1.75E-14, 1.86E-14, 1.98E-14, 2.12E-14, &
         2.24E-14, 2.34E-14, 2.48E-14, 2.56E-14, 2.64E-14, 2.73E-14, 2.81E-14, 2.90E-14, 2.95E-14, 2.98E-14, &
         3.01E-14, 3.06E-14, 3.06E-14, 3.04E-14, 3.01E-14, 2.95E-14, 2.90E-14, 2.87E-14, 2.77E-14, 2.72E-14, &
         2.63E-14, 2.55E-14, 2.48E-14, 2.37E-14, 2.28E-14, 2.17E-14, 2.07E-14, 1.99E-14, 1.89E-14, 1.80E-14, &
         1.71E-14, 1.64E-14, 1.56E-14 /)

    ! ! Hydrogen-electron impact ionization cross section depending on electron energy
    ! ElCrossSect_I(1:61) = (/ 0.000E-00, 0.114E-20, 0.189E-20, 0.257E-20, 0.320E-20, 0.377E-20, &
    !      0.429E-20, 0.473E-20, 0.510E-20, 0.542E-20, 0.568E-20, 0.587E-20, 0.600E-20, 0.609E-20, 0.613E-20, &
    !      0.613E-20, 0.608E-20, 0.600E-20, 0.588E-20, 0.575E-20, 0.559E-20, 0.541E-20, 0.521E-20, 0.501E-20, &
    !      0.480E-20, 0.457E-20, 0.435E-20, 0.414E-20, 0.392E-20, 0.369E-20, 0.348E-20, 0.328E-20, 0.307E-20, &
    !      0.288E-20, 0.270E-20, 0.252E-20, 0.235E-20, 0.219E-20, 0.204E-20, 0.190E-20, 0.176E-20, 0.163E-20, &
    !      0.152E-20, 0.141E-20, 0.130E-20, 0.120E-20, 0.111E-20, 0.103E-20, 0.095E-20, 0.088E-20, 0.081E-20, &
    !      0.075E-20, 0.069E-20, 0.063E-20, 0.059E-20, 0.054E-20, 0.050E-20, 0.045E-20, 0.042E-20, 0.039E-20, &
    !      0.036E-20 /)
    !----------------------------------------------------------------------

    !! provide all rates in SI units

    !! initialize all collision rates with zero
    fin_II = 0. ; fii_II = 0. ; fie_I = 0. ; kin_IIII = 0. ; alpha_I = 0. ; v_II = 0. ; ve_II = 0.
    sigma_e = 0. ; Qexc_II = 0. ; Qion_II = 0.

    !! Ionization rates (add opacity correction/shadowing when needed)
    v_II(H2O_,Hp_)   = (1.30E-8+3.26E-8)/(rHelio**2) !! Photoionization rate producing Hp from H2O & OH [1/s] (Huebner et al. 1992)
    v_II(H_,Hp_)     = 7.30E-8/(rHelio**2)           !! Photoionization rate producing Hp from H [1/s] (Huebner et al. 1992)
    v_II(H2O_,H2Op_) = vHI                           !! from PARAM.in
    v_II(H_,SW_)     = v_II(H_,Hp_)*1e-9             !! set production of solar wind protons to a small but non-zero value

    !! Electron excess energies from ionization (increases electron pressure)
    Qexc_II(H2O_,Hp_)   = 4.0054E-18 ! 25.0 eV, Huebner 1992
    Qexc_II(H_,Hp_)     = 5.6076E-19 !  3.5 eV, Huebner 1992
    Qexc_II(H2O_,H2Op_) = 1.9226E-18 ! 12.0 eV, Huebner 1992

    !! Ionization potential for electron impact ionization (needs to be delivered by the ionizing electron)
    Qion_II(H2O_,Hp_)   = 3.70e-18 ! 23.1 eV (estimate) 5.11 eV bond strength OH-H and 4.40 eV O-H and 13.6 eV H-> H+ (Wikipedia)
    Qion_II(H2O_,H2Op_) = 2.02e-18 ! 12.6 eV Joshipura et al. (2007)
    Qion_II(H_,Hp_)     = 2.18e-18 ! 13.6 eV

    ! ! UV opacity
    J3 = 4.5E14/(rHelio**2) ! J3 = 4.5E14 [m^-2*s^-1]: lambda < 984A solar flux @ 1 AU, Marconi, 1982)
    sigma = (v_II(H2O_,Hp_)+v_II(H_,Hp_)+v_II(H2O_,H2Op_))/J3
    ! Alternative:
    ! Cross section of J3 flux for ionization (lambda < 984A) [m^2], Marconi, 1982)
    ! sigma_13=2.4E-18 cm^2: H2O + hv -> OH + H
    ! sigma_23=1.0E-18 cm^2: H2O + hv -> H2 + O(1D)
    ! sigma_33=8.4E-18 cm^2: H2O + hv -> H2Op + e
    ! sigma = 1.18E-21 ! sum(sigma_i3)

    ! Dist distance from sun-comet line, only neutral H2O considered
    Dist = sqrt(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)*&
         rPlanetSI+0.1
    if (Dist.le.rBody*rPlanetSI .and. Xyz_DGB(x_,i,j,k,iBlock).le.0. .and. Body1) then
       v_II = v_II*1e-9 ! Inside the body's shadow
    else
       ! N total number of water-type molecules in upstream column
       uNeutr = sqrt(UnxNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2+UnyNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2+&
            UnzNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2)
       NCol = Qprod/uNeutr/Dist/4.*(-atan(Xyz_DGB(x_,i,j,k,iBlock)*rPlanetSI/Dist)/cPi+0.5)
       v_II = v_II*exp(-sigma*NCol) + v_II*1e-9
       !if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
       !   write(*,*)'sigma      = ',sigma
       !   write(*,*)'Dist       = ',Dist
       !   write(*,*)'NCol       = ',NCol
       !   write(*,*)'Correction = ',exp(-sigma*NCol)
       !end if
    end if

    log10Te = log10(Te)
    ! H2O electron impact ionization cross section after Cravens et al. 1987: H2O & e -> H2Op + 2e
    n = int((log10Te-4.45053516)/0.0415 + 1.0)
    if (n > 60) n = 60
    if (n < 1) n = 1
    ve_II(H2O_,H2Op_) = max(nElec*((log10Te-((n-1.0)*0.0415+4.45053516))/0.0415*(ElImpRate_I(H2O_,n+1)-ElImpRate_I(H2O_,n))+&
         ElImpRate_I(H2O_,n)),0.0) ! linear interpolation
    v_II(H2O_,H2Op_) = v_II(H2O_,H2Op_) + ve_II(H2O_,H2Op_) ! v_II is the total ionization rate, photons and electrons!
    ! TestArray(1,i,j,k,iBlock) = ve_II(H2O_,H2Op_)/v_II(H2O_,H2Op_)

    ! Hydrogen electron impact ionization cross section after Cravens et al. 1987: H & e -> Hp + 2e
    n = int((log10Te-4.527)/0.0406 + 1.0)
    if (n > 60) n = 60
    if (n < 1) n = 1
    ve_II(H_,Hp_) = max(nElec*((log10Te-((n-1.0)*0.0406+4.527))/0.0406*(ElImpRate_I(H_,n+1)-ElImpRate_I(H_,n))+&
         ElImpRate_I(H_,n)),0.0) ! linear interpolation
    v_II(H_,Hp_) = v_II(H_,Hp_) + ve_II(H_,Hp_) ! v_II is the total ionization rate, photons and electrons!
    ! TestArray(1,i,j,k,iBlock) = ve_II(H_,Hp_)/v_II(H_,Hp_)

    ! Number of energetic electrons (number of equivalent 20 eV electrons)
    ne20eV_GB(i,j,k,iBlock) = ve_II(H2O_,H2Op_)/(sigmaeh2o*ve)



    ! ! Evaluate the number of energetic electrons
    ! ! ! !!Hydrogen electron impact ionization rate: H & e -> Hp + 2e
    ! ueBulk2  = sum(uElec_D(:)**2)                                ! Electron bulk speed qubed in [m^2/s^2]
    ! ueTherm2 = 3.*cBoltzmann*Te/(cElectronMass)                  ! Electron thermal speed qubed in [m^2/s^2]
    ! Ee = 0.5*cElectronMass*(ueBulk2+ueTherm2)/cElectronCharge    ! Electron energy in [eV]
    ! !! BEQ electron impact cross section from Y.-K. Kim and M.E. Rudd, Phys. Rev. A 50, 3954 (1994)
    ! n = int((log10(Ee)-1.15286905)/0.0408 + 1.0)
    ! if (n > 61) n = 61
    ! if (n < 1) n = 1
    ! sigma_e(H_,Hp_) = ElCrossSect_I(n)             ! electron impact cross section in [m^2]
    ! if (n > 1) then
    !    ne20eV_GB(i,j,k,iBlock) = ve_II(H_,Hp_)/(sigma_e(H_,Hp_)*sqrt(ueBulk2+ueTherm2))
    ! else
    !    ne20eV_GB(i,j,k,iBlock) = 0.
    ! end if
    ! ! !! Gombosi book f12=n2*sigma12*sqrt(v1**2+v2**2)           
    ! ! ve_II(H_,Hp_) = nElec*sigma_e(H_,Hp_)*sqrt(ueBulk2+ueTherm2) ! neutral velocity neglected
    ! ! v_II(H_,Hp_) = v_II(H_,Hp_) + ve_II(H_,Hp_)                  ! v_II is the total ionization rate, photons and electrons!
    ! ! write(*,*)'Ee = ', Ee, ' Te = ', Te,' ElCrossSect_I(n) = ', ElCrossSect_I(n),' n = ', n
    ! ! STOP



    !! ********** Ion-neutral collision/charge exchange rates ********** 
    !! Example(s)
    ! resonant H+ & O -> O+ & H  subtracts H+ and adds O+
    ! kin_IIII(Hp_,O_,Op_,H_) = 6.61E-11/1E6*sqrt(Ti_I(Hp_))*(1.0-0.047*log10(Ti_I(Hp)))**2    !! rate in [m^3/s]
    ! resonant O+ & H -> H+ & O  subtracts O+ and adds H+
    ! kin_IIII(Op_,H_,Hp_,O_) = 4.63E-12/1E6*sqrt(TnNeutral(H_,i-MinI+1,j-MinJ+1,k-MinK+1)+TOp_/16.)  !! rate in [m^3/s]


    !! H2Op & H2O -> H2Op & H2O    ! non-resonant
    !! fin_II(H2Op_,H2O_) = 0.
    !! H2Op & H2O -> H2O & H2Op    ! resonant
    kin_IIII(H2Op_,H2O_,H2O_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996)

    !! H2Op & H -> H2Op & H    ! non-resonant
    !!fin_II(H2Op_,H_) = 0.
    !! H2Op & H -> H2O & Hp    ! resonant
    kin_IIII(H2Op_,H_,H2O_,Hp_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated

    !! Hp & H -> H & Hp  ! resonant, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tr=(Ti_I(Hp_)+TnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1))/2.
    kin_IIII(Hp_,H_,H_,Hp_)  = 2.65E-10/1E6*sqrt(Tr)*(1.0-0.083*log10(Tr))**2  !! rate in [m^3/s]
    !! SWp & H -> H & Hp  ! resonant, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tr=(Ti_I(SW_)+TnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1))/2.
    kin_IIII(SW_,H_,H_,Hp_) = 2.65E-10/1E6*sqrt(Tr)*(1.0-0.083*log10(Tr))**2  !! rate in [m^3/s]
    !! Hp & H -> H & Hp  ! non-resonant
    !fin_II(Hp_,H_) = 0.

    !! Hp & H2O -> H & H2Op    ! resonant, Benna et al., Plan. Sp. Sci., (2007)
    !!    uHpBulk = sqrt(sum(State_VGB(HpRhoUx_:HpRhoUz_,i,j,k,iBlock)**2)) / &
    !!         State_VGB(HpRho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    !!    uHpTherm = sqrt(3.*cBoltzmann*Tion/(MassIon(Hp_)*cProton))
    !! Benna et al., Plan. Sp. Sci., (2007)
    !! vt = sqrt(v*v + 2.11E8*t_p) ??? chemistry code
    !!???kin_IIII(Hp_,H2O_,H_,H2Op_) = 1E-4*2.1E-15*(uHpBulk + uHpTherm) < -calc??? Benna et al., Plan. Sp. Sci., (2007)

    !! Hp & H2O -> H & H2Op    ! resonant, estimate to get the same drag on SW-protons
    kin_IIII(Hp_,H2O_,H_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated
    kin_IIII(SW_,H2O_,H_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated

    !! ********** Ion-ion collision rates ********** 
    ! SWp - SWp is left zero because the they do not result in a change in the source terms
    ! Hp - Hp is left zero because the they do not result in a change in the source terms
    ! H2Op - H2Op is left zero because the they do not result in a change in the source terms

    ! H2Op - Hp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000    
    Tred = (MassIon_I(H2Op_)*Ti_I(Hp_)+MassIon_I(Hp_)*Ti_I(H2Op_))/(MassIon_I(H2Op_)+MassIon_I(Hp_)) ! reduced temp
    Mred = MassIon_I(H2Op_)*MassIon_I(Hp_)/(MassIon_I(H2Op_)+MassIon_I(Hp_)) ! reduced mass
    fii_II(H2Op_,Hp_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(Hp_)**2/MassIon_I(H2Op_)*&
         sqrt(Mred)*1e-6*nIon_I(Hp_)/(Tred*sqrt(Tred))
    ! Hp - H2Op, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(Hp_,H2Op_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(Hp_)**2/MassIon_I(Hp_)*&
         sqrt(Mred)*1e-6*nIon_I(H2Op_)/(Tred*sqrt(Tred))

    ! H2Op - SWp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tred = (MassIon_I(H2Op_)*Ti_I(SW_)+MassIon_I(SW_)*Ti_I(H2Op_))/(MassIon_I(H2Op_)+MassIon_I(SW_)) ! reduced temp
    Mred = MassIon_I(H2Op_)*MassIon_I(SW_)/(MassIon_I(H2Op_)+MassIon_I(SW_)) ! reduced mass
    fii_II(H2Op_,SW_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SW_)**2/MassIon_I(H2Op_)*&
         sqrt(Mred)*1e-6*nIon_I(SW_)/(Tred*sqrt(Tred))
    ! SWp - H2Op, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(SW_,H2Op_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SW_)**2/MassIon_I(SW_)*&
         sqrt(Mred)*1e-6*nIon_I(H2Op_)/(Tred*sqrt(Tred))

    ! SWp - Hp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tred = (MassIon_I(SW_)*Ti_I(Hp_)+MassIon_I(Hp_)*Ti_I(SW_))/(MassIon_I(SW_)+MassIon_I(Hp_)) ! reduced temp
    Mred = MassIon_I(SW_)*MassIon_I(Hp_)/(MassIon_I(SW_)+MassIon_I(Hp_)) ! reduced mass
    fii_II(SW_,Hp_) = 1.27*ChargeIon_I(SW_)**2*ChargeIon_I(Hp_)**2/MassIon_I(SW_)*&
         sqrt(Mred)*1e-6*nIon_I(Hp_)/(Tred*sqrt(Tred))
    ! Hp - SWp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(Hp_,SW_) = 1.27*ChargeIon_I(SW_)**2*ChargeIon_I(Hp_)**2/MassIon_I(Hp_)*&
         sqrt(Mred)*1e-6*nIon_I(SW_)/(Tred*sqrt(Tred))

    !! Ion - electron collision rates, reduced mass=~me and reduced temperature=~Te
    sqrtTe = sqrt(Te)
    !! H2Op - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(H2Op_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(H2Op_)*ChargeIon_I(H2Op_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]
    !! Hp - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(Hp_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(Hp_)*ChargeIon_I(Hp_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]
    fie_I(SW_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(SW_)*ChargeIon_I(SW_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]

    !! ********** Ion-electron recombination rates ********** 

    !! Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    if (Te < 800.) then 
       alpha_I(H2Op_) = 1E-6*1.57E-5*Te**(-0.569) !! rate in [m^3/s]
    elseif (Te<4000) then
       alpha_I(H2Op_) = 1E-6*4.73E-5*Te**(-0.74)  !! rate in [m^3/s]
    else
       alpha_I(H2Op_) = 1E-6*1.03E-3*Te**(-1.111) !! rate in [m^3/s]
    end if

    !     if (Te < 200.) then 
    !        alpha_I(H2Op_) = 1E-6*7E-7*sqrt(300./Te) !! rate in [m^3/s]
    !     else
    !        alpha_I(H2Op_) = 2.342*1E-6*7E-7*Te**(0.2553-0.1633*log10(Te)) !! rate in [m^3/s]
    !     end if

    alpha_I(Hp_)   = 1E-6*4.8E-12*(250/Te)**0.7  !! Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    alpha_I(SW_)  = alpha_I(Hp_)
    !alpha_I(Hp_)   = 1E-6*3.5E-12*(Te/300)**(-0.7)  !! Schmidt et al., Comput. Phys. Commun. (1988)

  end subroutine user_calc_rates

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
    use ModPhysics,    ONLY: SW_Ux, SW_Uy, SW_Uz, UnitN_, UnitRho_, UnitU_, UnitP_, UnitT_, UnitB_, &
         ElectronPressureRatio, ElectronCharge, Si2No_V, No2Si_V, UnitEnergyDens_, UnitJ_, UnitRhoU_, UnitX_
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit, IsPointImplSource

    integer, intent(in) :: iBlock

    real, dimension(1:nI,1:nJ,1:nK) :: nElec_C, Te_C, SBx_C, SBy_C, SBz_C, SPe_C
    real, dimension(4,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoTerm_IIC
    real, dimension(5,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoUxTerm_IIC, SRhoUyTerm_IIC, SRhoUzTerm_IIC
    real, dimension(8,1:nIonFluid,1:nI,1:nJ,1:nK) :: SPTerm_IIC
    real, dimension(8,1:nI,1:nJ,1:nK) :: SPeTerm_IC

    real, dimension(1:3,1:nI,1:nJ,1:nK) :: Current_DC, uIonMean_DC, uElec_DC
    real, dimension(1:3,1:nIonFluid,1:nI,1:nJ,1:nK) :: uIon_DIC
    real, dimension(1:nIonFluid,1:nNeutral,1:nI,1:nJ,1:nK) :: fin_IIC, uIonNeu2_IIC
    real, dimension(1:nNeutral,1:nI,1:nJ,1:nK) :: fen_IC, uNeuElec2_IC
    real, dimension(1:nIonFluid,1:nIonFluid,1:nI,1:nJ,1:nK) :: fii_IIC, uIonIon2_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: Ti_IC, uIonElec2_IC, fei_IC, fie_IC, &
         nIon_IC, SRho_IC, SRhoUx_IC, SRhoUy_IC, SRhoUz_IC, SP_IC
    real, dimension(1:nNeutral,1:nIonFluid) :: Qexc_II, Qion_II
    real, dimension(1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: v_IIC, ve_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: alpha_IC
    real, dimension(1:nIonFluid) :: fiiTot_I, finTot_I, vAdd_I, kinAdd_I, kinSub_I
    real, dimension(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: kin_IIIIC

    logical :: DoTest, DoTestMe=.true.
    real :: theta, fenTot, feiTot,logTe
    integer :: i,j,k,iNeutral,jNeutral,iIonFluid,jIonFluid,iTerm,iDim

    real :: a, b, d, l, f=5./8., MaxHeat
    !----------------------------------------------------------------------

    ! Do not evaluate any source terms explicitly when running pointimplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Evaluate source terms explicitly even when running pointimplicit
    !if(UsePointImplicit .and. IsPointImplSource) RETURN

    !! Limit region for evaluation for source term evaluation
    !! if(RMin_BLK(iBlock) > 2.) RETURN

    if(iBlock == BlkTest.and.iProc==ProcTest) then
       call set_oktest('user_calc_sources',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    !! Set the source arrays for this block to zero
    SRho_IC        = 0.
    SRhoTerm_IIC   = 0.
    SRhoUx_IC      = 0.
    SRhoUxTerm_IIC = 0.
    SRhoUy_IC      = 0.
    SRhoUyTerm_IIC = 0.
    SRhoUz_IC      = 0.
    SRhoUzTerm_IIC = 0.
    SBx_C          = 0.
    SBy_C          = 0.
    SBz_C          = 0.
    SP_IC          = 0.
    SPTerm_IIC     = 0.
    SPe_C          = 0.
    SPeTerm_IC     = 0.


    ! nElec_C is the electron/ion density in SI units ( n_e=sum(n_i*Zi) )
    do k=1,nK; do j=1,nJ; do i=1,nI
       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    !! ion velocity components in SI

    uIon_DIC(1,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUxIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(2,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUyIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(3,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUzIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIonMean_DC(1:3,1:nI,1:nJ,1:nK) = 0.
    do iIonFluid=1,nIonFluid
       uIonMean_DC(1,1:nI,1:nJ,1:nK) = uIonMean_DC(1,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(2,1:nI,1:nJ,1:nK) = uIonMean_DC(2,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(3,1:nI,1:nJ,1:nK) = uIonMean_DC(3,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid) 
    end do

    !! (u_i-u_n)^2 in SI
    do iIonFluid=1,nIonFluid
       do iNeutral=1,nNeutral
          uIonNeu2_IIC(iIonFluid,iNeutral,1:nI,1:nJ,1:nK) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnxNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnyNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnzNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2
       end do
    end do

    !! (u_i1-u_i2)^2 in SI
    do iIonFluid=1,nIonFluid
       do jIonFluid=1,nIonFluid
          uIonIon2_IIC(iIonFluid,jIonFluid,:,:,:) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(1,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(2,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(3,jIonFluid,1:nI,1:nJ,1:nK))**2
       end do
    end do

    if (UseElectronPressure) then
       ! Electron temperature calculated from electron pressure
       ! Ion temperature is calculated from ion pressure
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/&
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
               nElec_C(i,j,k))
       end do; end do; end do
    else
       ! Electron temperature calculated from pressure assuming Te_C=Ti_IC*ElectronTemperatureRatio:
       ! p=nkT with n_e=n_i*Z_i (quasi-neutrality), n=n_e+n_i and p=p_e+p_i=p_i*(1+ElectronPressureRatio)
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
               NO2SI_V(UnitP_)/(cBoltzmann*nElec_C(i,j,k))
       end do; end do; end do
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! No need to evaluate source terms for cells inside the body
       ! if((Body1).and.(R_BLK(i,j,k,iBlock)<rBody)) CYCLE

       call get_current(i,j,k,iBlock,Current_DC(:,i,j,k))

       ! calculate uElec_DC from Hall velocity -J/(e*n) [m/s]
       uElec_DC(1:3,i,j,k) = uIonMean_DC(1:3,i,j,k)-Current_DC(1:3,i,j,k)/(nElec_C(i,j,k)*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_)

       call calc_electron_collision_rates(Te_C(i,j,k),nElec_C(i,j,k),i,j,k,iBlock,fen_IC(1:nNeutral,i,j,k), &
            fei_IC(1:nIonFluid,i,j,k))
       call user_calc_rates(Ti_IC(1:nIonFluid,i,j,k),Te_C(i,j,k),i,j,k,iBlock,nElec_C(i,j,k),nIon_IC(1:nIonFluid,i,j,k),&
            fin_IIC(1:nIonFluid,1:nNeutral,i,j,k),fii_IIC(1:nIonFluid,1:nIonFluid,i,j,k),fie_IC(1:nIonFluid,i,j,k),&
            alpha_IC(1:nIonFluid,i,j,k),kin_IIIIC(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,i,j,k),&
            v_IIC(1:nNeutral,1:nIonFluid,i,j,k),ve_IIC(1:nNeutral,1:nIonFluid,i,j,k),uElec_DC(1:3,i,j,k),&
            uIon_DIC(1:3,1:nIonFluid,i,j,k),Qexc_II(1:nNeutral,1:nIonFluid),Qion_II(1:nNeutral,1:nIonFluid))

       !! Zeroth moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinAdd_I = 0. ; kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do

       !! Sources divideded into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       SRhoTerm_IIC(1,1:nIonFluid,i,j,k) = vAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I              !! newly ionized neutrals
       SRhoTerm_IIC(2,1:nIonFluid,i,j,k) = kinAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I            !! mass added through ion-neutral charge exchange
       SRhoTerm_IIC(3,1:nIonFluid,i,j,k) = -kinSub_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I           !! mass removed through ion-neutral charge exchange
       SRhoTerm_IIC(4,1:nIonFluid,i,j,k) = -alpha_IC(1:nIonFluid,i,j,k)*(nElec_C(i,j,k)* &                            !! loss due to recombination
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)/Si2No_V(UnitT_))*MassIon_I


       !! First moment, x component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(1,iIonFluid,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUxTerm_IIC(1,iIonFluid,i,j,k) = SRhoUxTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(1,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(1,i,j,k)
       SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(1,1:nIonFluid,i,j,k)-uIonMean_DC(1,i,j,k))*Si2No_V(UnitU_)
       SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, y component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(2,iIonFluid,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do


       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*&
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUyTerm_IIC(1,iIonFluid,i,j,k) = SRhoUyTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(2,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(2,i,j,k)
       SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(2,1:nIonFluid,i,j,k)-uIonMean_DC(2,i,j,k))*Si2No_V(UnitU_)
       SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, z component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(3,iIonFluid,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                 !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                      !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUzTerm_IIC(1,iIonFluid,i,j,k) = SRhoUzTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(3,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &            !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(3,i,j,k)
       SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(3,1:nIonFluid,i,j,k)-uIonMean_DC(3,i,j,k))*Si2No_V(UnitU_)
       SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &   !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &   !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       ! (u_n-u_e)^2 difference in neutral and electron speeds qubed [m^2/s^2]
       do iNeutral=1,nNeutral
          uNeuElec2_IC(iNeutral,i,j,k) = ((UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(1,i,j,k))**2 &
               +(UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(2,i,j,k))**2 &
               +(UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(3,i,j,k))**2)
       end do

       ! (u_i-u_e)^2 difference in ion and electron speeds qubed [m^2/s^2]
       do iIonFluid=1,nIonFluid
          uIonElec2_IC(iIonFluid,i,j,k) = (uIon_DIC(1,iIonFluid,i,j,k)-uElec_DC(1,i,j,k))**2+&
               (uIon_DIC(2,iIonFluid,i,j,k)-uElec_DC(2,i,j,k))**2+&
               (uIon_DIC(3,iIonFluid,i,j,k)-uElec_DC(3,i,j,k))**2
       end do

       !! Second moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + &!!nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)*&
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do

       SPTerm_IIC(1,1:nIonFluid,i,j,k) = -(kinSub_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! lost ions through charge exchange and recombination
            alpha_IC(1:nIonFluid,i,j,k)*nElec_C(i,j,k)/Si2No_V(UnitT_))*State_VGB(iPIon_I,i,j,k,iBlock)

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               cBoltzmann*(Ti_IC(iIonFluid,i,j,k)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(2,1:nIonFluid,i,j,k) = 2.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)        

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+NeutralMass_I(iNeutral)/cProtonMass)*&
               cBoltzmann*(TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(3,1:nIonFluid,i,j,k) = 2.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(4,1:nIonFluid,i,j,k) = 2.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*nIon_IC(1:nIonFluid,i,j,k)*&
            cBoltzmann*(Te_C(i,j,k)-Ti_IC(1:nIonFluid,i,j,k))*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(5,1:nIonFluid,i,j,k) = 2./3.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*cElectronMass*&             !! ion-electron collisional exchange (due to Hall velocity)
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*uIonElec2_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitU_)**2

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*MassIon_I(iIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               uIonIon2_IIC(1:nIonFluid,iIonFluid,i,j,k)
       end do

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*NeutralMass_I(iNeutral)/(MassIon_I(1:nIonFluid)+NeutralMass_I(iNeutral)/cProtonMass)*&
               uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)/cProtonMass
       end do

       SPTerm_IIC(6,1:nIonFluid,i,j,k) = 2./3.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2
       SPTerm_IIC(7,1:nIonFluid,i,j,k) = 2./3.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2


       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)
       end do
       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m*s^2)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)*&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(jIonFluid,iNeutral,i,j,k)
                end do
             end do
          end do
       end do

       SPTerm_IIC(8,1:nIonFluid,i,j,k) = 1./3.*(vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*&
            MassIon_I(1:nIonFluid)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2

       if (UseElectronPressure) then
          SPeTerm_IC(1,i,j,k) = -sum(alpha_IC(1:nIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k))/ &                           !! lost electrons through recombination
               Si2No_V(UnitT_)*State_VGB(Pe_,i,j,k,iBlock)

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uNeuElec2_IC(iNeutral,i,j,k)
          end do
          SPeTerm_IC(2,i,j,k) = 1./3.*cElectronMass*sum(vAdd_I)*Si2No_V(UnitRho_)/Si2No_V(UnitT_)*Si2No_V(UnitU_)**2      !! new electrons through photoionized neutrals

          feiTot = 0.
          do iIonFluid=1,nIonFluid
             feiTot = feiTot+fei_IC(iIonFluid,i,j,k)/MassIon_I(iIonFluid)*&
                  (Ti_IC(iIonFluid,i,j,k)-Te_C(i,j,k))
          end do
          SPeTerm_IC(3,i,j,k) = 2.*cElectronMass*Si2No_V(UnitRho_)/Si2No_V(UnitN_)*&                                      !! ion-electron collisional exchange (thermal motion)
               nElec_C(i,j,k)*cBoltzmann*Si2No_V(UnitEnergyDens_)*feiTot/Si2No_V(UnitT_)

          fenTot = 0.
          do iNeutral=1,nNeutral
             fenTot = fenTot+fen_IC(iNeutral,i,j,k)/NeutralMass_I(iNeutral)*&
                  (TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Te_C(i,j,k))
          end do
          SPeTerm_IC(4,i,j,k) = 2.*cElectronMass*nElec_C(i,j,k)*cBoltzmann*&                                              !! electron-neutral collisional exchange (thermal motion)
               Si2No_V(UnitEnergyDens_)*fenTot/Si2No_V(UnitT_)

          SPeTerm_IC(5,i,j,k) = 2./3.*sum(fei_IC(1:nIonFluid,i,j,k)*uIonElec2_IC(1:nIonFluid,i,j,k))/ &                   !! ion-electron collisional exchange (due to Hall velocity)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          SPeTerm_IC(6,i,j,k) = 2./3.*sum(fen_IC(1:nNeutral,i,j,k)*uNeuElec2_IC(1:nNeutral,i,j,k))/&                      !! electron-neutral collisional exchange (bulk motion)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          ! ! Deriving the heating according to Gan and Cravens 1990
          ! ! Parabola x=a*d^2+b
          ! d = sqrt(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)*No2SI_V(UnitX_) ! distance d=sqrt(y^2+z^2) from Sun-comet axis [m]
          ! b = f*d + Xyz_DGB(x_,i,j,k,iBlock)*No2SI_V(UnitX_) ! parabola's vertex point (b, 0.0, 0.0) with f=5/8 estimate for parobola shape
          ! a = (Xyz_DGB(x_,i,j,k,iBlock)*No2SI_V(UnitX_)-b)/(d**2+1e-9)
          ! MaxHeat = max(-7.326e-20*b+1.744e-12, 0.0) ! Maximum heating at the parabola's vertex point (b, 0.0, 0.0)
          ! SPeTerm_IC(7,i,j,k) = MaxHeat 
          ! l = 0.0 ! path length along the parabola from vertex point
          ! if (a < 0.) then
          !    l = (2.*d*sqrt(1.+4.*a**2*d**2)*a+log(2.*a*d+sqrt(4.*a**2*d**2+1.)))/(4.*a)
          !    if (b > 0.) then
          !       SPeTerm_IC(7,i,j,k) = min(MaxHeat, 50.29*l**(-1.981))
          !    else
          !       SPeTerm_IC(7,i,j,k) = 0.0
          !    end if
          ! end if
          ! ! if ((Xyz_DGB(x_,i,j,k,iBlock)<0.0).and.&
          ! !     (Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2 < 0.005**2)) then
          ! !   SPeTerm_IC(7,i,j,k) = 0.
          ! !end if
          ! SPeTerm_IC(7,i,j,k) = SPeTerm_IC(7,i,j,k)*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
          ! ! TestArray(1,i,j,k,iBlock) = l
          ! ! TestArray(2,i,j,k,iBlock) = MaxHeat
          ! ! TestArray(3,i,j,k,iBlock) = b

          ! if (Xyz_DGB(x_,i,j,k,iBlock) < -160.*(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)+0.004) then
          !    SPeTerm_IC(7,i,j,k) = 0.0
          !    ! TestArray(1,i,j,k,iBlock) = TestArray(1,i,j,k,iBlock)/3.                                                                                                                                          
          ! end if

          ! TestArray(1,i,j,k,iBlock) = SPeTerm_IC(7,i,j,k)

          ! ! if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
          ! !    write(*,*)'x         = ',Xyz_DGB(x_,i,j,k,iBlock)," [rPlanet]"
          ! !    write(*,*)'y         = ',Xyz_DGB(y_,i,j,k,iBlock)," [rPlanet]"
          ! !    write(*,*)'z         = ',Xyz_DGB(z_,i,j,k,iBlock)," [rPlanet]"
          ! !    write(*,*)'l = ', l
          ! !    write(*,*)'MaxHeat = ', MaxHeat
          ! !    write(*,*)'b = ', b
          ! !    write(*,*)'SPeTerm = ', SPeTerm_IC(7,i,j,k)*No2SI_V(UnitEnergyDens_)/No2SI_V(UnitT_)
          ! ! end if

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+((v_IIC(iNeutral,1:nIonFluid,i,j,k)-&
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k))*Qexc_II(iNeutral,1:nIonFluid)- &
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k)*Qion_II(iNeutral,1:nIonFluid))*ChargeIon_I(1:nIonFluid)* &          
                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
          end do
          SPeTerm_IC(7,i,j,k) = 2./3.*sum(vAdd_I)*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                                ! heating of electrons due to ionization excess energy

          SPeTerm_IC(7,i,j,k) = SPeTerm_IC(7,i,j,k)*1./max(1.,((-300.*R_BLK(i,j,k,iBlock)+7.65)*vHI/2e-6))                ! fit to heating ratio of Gan & Cravens 1990

          if (Xyz_DGB(x_,i,j,k,iBlock) < -160.*(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)+0.004) then
             SPeTerm_IC(7,i,j,k) = 0.0
             ! TestArray(1,i,j,k,iBlock) = TestArray(1,i,j,k,iBlock)/3.                                                                                                                                          
          end if

          logTe = log(Te_C(i,j,k))
          SPeTerm_IC(8,i,j,k) = exp(-188.4701+33.2547*logTe-2.0792*logTe**2+0.0425*logTe**3)                              !! electron cooling due to collisions w/ water vapor
          if(Te_C(i,j,k)<1.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)) then
             SPeTerm_IC(8,i,j,k)=4.5e-9/(0.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))* &
                  (Te_C(i,j,k)-TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))
          else
             SPeTerm_IC(8,i,j,k)=SPeTerm_IC(8,i,j,k)+4.5e-9
          end if

          SPeTerm_IC(8,i,j,k) = -2./3.*NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*nElec_C(i,j,k)* &
               SPeTerm_IC(8,i,j,k)/1e6*1.60217733e-19*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                           

       end if

       !! sum up individual terms
       do iTerm=1,4
          SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(4,1:nIonFluid,i,j,k)
       do iTerm=1,5
          SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k)
       do iTerm=1,8
          SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(6,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(7,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(8,1:nIonFluid,i,j,k)
       if(UseElectronPressure) then
          SPe_C(i,j,k) = sum(SPeTerm_IC(1:8,i,j,k))
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(1,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(2,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(3,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(4,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(5,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(6,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(7,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(8,i,j,k)
       end if

       Source_VC(iRhoIon_I   ,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)    + Source_VC(iRhoIon_I   ,i,j,k)
       Source_VC(iRhoUxIon_I ,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUxIon_I ,i,j,k)
       Source_VC(iRhoUyIon_I ,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUyIon_I ,i,j,k)
       Source_VC(iRhoUzIon_I ,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUzIon_I ,i,j,k)
       Source_VC(iPIon_I     ,i,j,k) = SP_IC(1:nIonFluid,i,j,k)      + Source_VC(iPIon_I     ,i,j,k)

       Source_VC(Rho_   ,i,j,k) = sum(SRho_IC(1:nIonFluid,i,j,k))    + Source_VC(Rho_   ,i,j,k)
       Source_VC(rhoUx_ ,i,j,k) = sum(SRhoUx_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUx_ ,i,j,k)
       Source_VC(rhoUy_ ,i,j,k) = sum(SRhoUy_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUy_ ,i,j,k)
       Source_VC(rhoUz_ ,i,j,k) = sum(SRhoUz_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUz_ ,i,j,k)
       Source_VC(Bx_    ,i,j,k) = SBx_C(i,j,k)                       + Source_VC(Bx_    ,i,j,k)
       Source_VC(By_    ,i,j,k) = SBy_C(i,j,k)                       + Source_VC(By_    ,i,j,k)
       Source_VC(Bz_    ,i,j,k) = SBz_C(i,j,k)                       + Source_VC(Bz_    ,i,j,k)
       if(UseElectronPressure) then
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))   + Source_VC(P_     ,i,j,k)
          Source_VC(Pe_    ,i,j,k) = SPe_C(i,j,k)                    + Source_VC(Pe_    ,i,j,k)
       else
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))*(1.+ElectronPressureRatio) + &
               Source_VC(P_     ,i,j,k)
       end if

    end do;  end do;  end do

    if(DoTestMe) then
       write(*,*)'user_calc_sources:'
       write(*,*)'Inputs: '
       i=iTest ; j=jTest ; k=kTest
       theta=acos((-SW_Ux*Xyz_DGB(x_,i,j,k,iBlock)-SW_Uy*Xyz_DGB(y_,i,j,k,iBlock)&
            -SW_Uz*Xyz_DGB(z_,i,j,k,iBlock))/R_BLK(i,j,k,iBlock)/&
            (SW_Ux**2+SW_Uy**2+SW_Uz**2)**0.5)
123    format (A13,ES25.16,A15,A3,F7.2,A3)
       write(*,123)'x         = ',Xyz_DGB(x_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'y         = ',Xyz_DGB(y_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'z         = ',Xyz_DGB(z_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'r         = ',R_BLK(i,j,k,iBlock)," [rPlanet]"
       write(*,123)'SW_Ux     = ',SW_Ux*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uy     = ',SW_Uy*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uz     = ',SW_Uz*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'Tmin      = ',Tmin," [K]"
       write(*,*)''
       write(*,*)'Neutrals:'
       do iNeutral=1,nNeutral
          write(*,124)'Neutral species #',iNeutral,': ', NameNeutral_I(iNeutral)," (",&
               NeutralMass_I(iNeutral)/cProtonMass," amu)"
          write(*,123)'n_n       = ',NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m^-3]"
          write(*,123)'m_n       = ',NeutralMass_I(iNeutral)," [kg]"
          write(*,123)'unx       = ',UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'uny       = ',UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'unz       = ',UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'Tn        = ',TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [K]"
       end do
       write(*,*)''
       write(*,*)'Total plasma phase (e- and i+):'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'uRhox     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'Ptot      = ',(State_VGB(P_,i,j,k,iBlock)+State_VGB(Pe_,i,j,k,iBlock))*&
               No2SI_V(UnitP_)," [kg/(m*s^2)]"
       else
          write(*,123)'Ptot      = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [kg/(m*s^2)]"
       end if
       write(*,123)'Bx        = ',State_VGB(Bx_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'By        = ',State_VGB(By_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'Bz        = ',State_VGB(Bz_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'uMeanx    = ',uIonMean_DC(1,i,j,k)," [m/s]"
       write(*,123)'uMeany    = ',uIonMean_DC(2,i,j,k)," [m/s]"
       write(*,123)'uMeanz    = ',uIonMean_DC(3,i,j,k)," [m/s]"       
       write(*,123)'jx        = ',Current_DC(1,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jy        = ',Current_DC(2,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jz        = ',Current_DC(3,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,*)''
       write(*,123)'SRho      = ',sum(SRho_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
            100.*sum(SRho_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(Rho_,i,j,k,iBlock)),"%)"
       if (State_VGB(RhoUx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUx_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUy_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUy_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUy_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUz_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(Bx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBx_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(By_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBy_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(By_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(Bz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBz_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if(UseElectronPressure) then
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(P_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)*&
               (1+ElectronPressureRatio)," [Pa/s]"," (",100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/ &
               (State_VGB(P_,i,j,k,iBlock))*(1+ElectronPressureRatio),"%)"
       end if
       write(*,*)''
       write(*,123)'dt        = ',Dt_BLK(iBlock)*No2SI_V(UnitT_)," [s]"
       write(*,*)''
       write(*,*)'Individual ion fluids:'
       do iIonFluid=1,nIonFluid
          write(*,124)'Ion species     #',iIonFluid,': ',NameFluid_I(iIonFluid+1)," (",&
               MassIon_I(iIonFluid)," amu/",ChargeIon_I(iIonFluid)," e)"
124       format (A17,I2,A3,A7,A3,F5.2,A5,F5.1,A3)
          write(*,123)'Ux        = ',uIon_DIC(1,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uy        = ',uIon_DIC(2,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uz        = ',uIon_DIC(3,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'ni        = ',nIon_IC(iIonFluid,i,j,k)," [m^-3]"
          write(*,123)'Ti        = ',Ti_IC(iIonFluid,i,j,k)," [K]"
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'SRho      = ',SRho_IC(iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRho_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT1   = ',SRhoTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT2   = ',SRhoTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT3   = ',SRhoTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT4   = ',SRhoTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          if (State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUx_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUy_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUz_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          write(*,123)'SP        = ',SP_IC(iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SP_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT1     = ',SPTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT2     = ',SPTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT3     = ',SPTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT4     = ',SPTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT5     = ',SPTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT6     = ',SPTerm_IIC(6,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(6,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT7     = ',SPTerm_IIC(7,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(7,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT8     = ',SPTerm_IIC(8,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(8,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
       end do
       write(*,*)''
       write(*,*)'Electrons:'
       write(*,123)'n_e       = ',nElec_C(i,j,k)," [m^-3]"
       if (UseElectronPressure) then
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if
       write(*,123)'Uex       = ',uElec_DC(1,i,j,k)," [m/s]"
       write(*,123)'Uey       = ',uElec_DC(2,i,j,k)," [m/s]"
       write(*,123)'Uez       = ',uElec_DC(3,i,j,k)," [m/s]"
       write(*,123)'Te        = ',Te_C(i,j,k)," [K]"
       if(UseElectronPressure) then
          if (State_VGB(Pe_,i,j,k,iBlock).gt.0.) then
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPe_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(1,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(2,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(3,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(4,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(5,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(6,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(7,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT8    = ',SPeTerm_IC(8,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(8,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
          else
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT8    = ',SPeTerm_IC(8,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
          end if
       end if
       write(*,*)''
       write(*,*)'Ion-electron combinations:'
       do iIonFluid=1,nIonFluid
          write(*,*)NameFluid_I(iIonFluid+1), '&  e'
          write(*,123)'fei       = ',fei_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'fie       = ',fie_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'|u_ime|   = ',sqrt(uIonElec2_IC(iIonFluid,i,j,k))," [m/s]"
          write(*,123)'alpha     = ',alpha_IC(iIonFluid,i,j,k)," [m^3/s]"
       end do
       write(*,*)''
       write(*,*)'Ion-ion combinations:'
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             if (iIonFluid.ne.jIonFluid) then
                write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameFluid_I(jIonFluid+1)
                write(*,123)' fii      = ',fii_IIC(iIonFluid,jIonFluid,i,j,k)," [1/s]"
                write(*,123)' |u_imi|  = ',sqrt(uIonIon2_IIC(iIonFluid,jIonFluid,i,j,k))," [m/s]"
             end if
          end do
       end do
       write(*,*)''
       write(*,*)'Ion-neutral combinations:'
       do iIonFluid=1,nIonFluid
          do iNeutral=1,nNeutral
             write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameNeutral_I(iNeutral)
             write(*,123)' v_phio   = ',v_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' v_eio    = ',ve_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' fin      = ',fin_IIC(iIonFluid,iNeutral,i,j,k)," [1/s]"
             write(*,123)' |u_imn|  = ',sqrt(uIonNeu2_IIC(iIonFluid,iNeutral,i,j,k))," [m/s]"
             write(*,*)' kin (Ion & Neutral-> Neutral & Ion):'
             do jIonFluid=1,nIonFluid
                do jNeutral=1,nNeutral
                   write(*,*)' ',NameFluid_I(iIonFluid+1),'&  ',NameNeutral_I(iNeutral),'->  ', &
                        NameNeutral_I(jNeutral),'&  ',NameFluid_I(jIonFluid+1),'=',&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)," [m^3/s]"  
                end do
             end do
          end do
       end do
       write(*,*)''
       write(*,*)'Electron-neutral combinations:'
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral)
          write(*,123)'fen      = ',fen_IC(iNeutral,i,j,k)," [1/s]"
          write(*,123)'|u_nme|  = ',sqrt(uNeuElec2_IC(iNeutral,i,j,k))," [m/s]"
       end do
       write(*,*)''
    end if

  end subroutine user_calc_sources

  !========================================================================

  subroutine user_update_states(iStage,iBlock)
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: SW_N, LowDensityRatio, cBoltzmann, ElectronPressureRatio, Si2No_V, &
         No2Si_V, UnitN_, UnitP_!, UnitB_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: r_BLK

    integer,intent(in) :: iStage, iBlock
    integer :: i,j,k,iIonFluid

    real, dimension(1:nI,1:nJ,1:nK) :: nElec_C
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) ::nIon_IC

    !----------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    ! Enforce minimum temperature (pressure), Tmin, if temperatures Ti_IC or Te_C are below

    do k=1,nK; do j=1,nJ; do i=1,nI

       do iIonFluid=1,nIonFluid
          ! set minimum mass density (and in these locations Ti = Tmin and vi=vbulkplasma)
          if(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) < SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2) then
             State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) = SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2
             State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)/ &
                  MassIon_I(iIonFluid)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       end do

       ! ! fix solar wind inside cavity to minimum value
       ! if(sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)*No2Si_V(UnitB_)**2 < 1e-9**2) then
       ! if(R_BLK(i,j,k,iBlock)<2.e-5) then
       !    State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) = SW_n*MassIon_I(SW_)*LowDensityRatio**2
       !    State_VGB(iRhoUxIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iRhoUyIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iRhoUzIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iPIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock)/ &
       !         MassIon_I(SW_)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)          
       ! end if

       State_VGB(Rho_,i,j,k,iBlock) = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))

       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))

       do iIonFluid=1,nIonFluid
          ! set minimum pressure
          if(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*NO2SI_V(UnitP_) < &
               nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin) then
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = &
                  nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       end do

       if(UseElectronPressure) then
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          if (State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_) < nElec_C(i,j,k)*cBoltzmann*Tmin) then
             State_VGB(Pe_,i,j,k,iBlock) = nElec_C(i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       else
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
       end if


    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !========================================================================

  subroutine derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I, nElec, EtaSi)
    use ModResistivity,  ONLY: Eta0SI
    use ModAdvance,      ONLY: State_VGB
    use ModVarIndexes,   ONLY: Rho_, Bx_, Bz_
    use ModPhysics,      ONLY: No2Si_V, UnitP_, UnitN_, UnitTemperature_, ElectronPressureRatio, &
         UnitB_, UnitRho_
    use ModConst,        ONLY: cBoltzmann, cElectronMass, cElectronCharge, &
         cMu, cProtonMass, cPi, cEps!, cLightSpeed, cTwoPi
    use ModB0,           ONLY: B0_DGB
    use ModMain,         ONLY: UseB0, iTest, jTest, kTest, BlkTest, ProcTest
    use ModProcMH,       ONLY: iProc

    integer, intent(in)  :: iBlock, i, j, k
    real,    intent(in)  :: TeSI
    real,    intent(in)  :: nIon_I(1:nIonFluid)
    real,    intent(in)  :: nElec
    real,    intent(out) :: EtaSi

    real :: EtaSiColl!, EtaSiSpitzer, lnL
    !    real, save :: SpitzerCoef, EtaPerpSpitzerSi
    logical, save :: FirstCall = .true.
    logical :: DoTest, DoTestMe=.true.
    real :: eeSigma!, B0_D(3)
    real, dimension(nIonFluid) :: fei_I, eiSigma_I
    real, dimension(nNeutral)  :: fen_I, enSigma_I
    integer :: iIonFluid, iNeutral

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('derive_cell_diffusivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! Spitzer formulation from Stoecker "Taschenbuch der Physik", Verlag "Harri Deutsch"
    ! lnL = log(1e7*TeSI**1.5/sqrt(nElec))
    ! EtaSiSpitzer = cElectronCharge**2*lnL/(32.*sqrt(2*cPi/cElectronMass*(cBoltzmann*TeSI)**3)*cEps**2)/cMu


    !! Collisional type resisitivity/diffusivity
    call calc_electron_collision_rates(TeSI,nElec,i,j,k,iBlock,fen_I(1:nNeutral),fei_I(1:nIonFluid))
    eiSigma_I(1:nIonFluid) = cElectronCharge**2*nElec/((fei_I(1:nIonFluid)+1E-20)*cElectronMass) 
    enSigma_I(1:nNeutral) = cElectronCharge**2*nElec/((fen_I(1:nNeutral)+1E-20)*cElectronMass)
    !! Eta_G is calculated from both conductivities using Kirchhoff's rule:
    !! 1/sigma_tot = 1/eiSigma_I+1/enSigma_I
    !! The resulting conductivity is close to Spitzer conductivity far from the comet and
    !! decreases due to abundant electron-neutral collisions close to the nucleus
    !! EtaSiColl = 1/(sigma_tot*mu_0) magnetic diffusivity [m^2/s]
    EtaSiColl = (sum(1/eiSigma_I(1:nIonFluid))+sum(1/enSigma_I(1:nNeutral)))/cMu
    !! Total diffusivity [m^2/s]
    EtaSi = Eta0SI + EtaSiColl

    ! TestArray(1,i,j,k,iBlock) = EtaSiColl
    ! TestArray(2,i,j,k,iBlock) = EtaSiSpitzer
    ! TestArray(3,i,j,k,iBlock) = EtaSiSpitzer/EtaSiColl

    if(DoTestMe) then
       write(*,*)'derive_cell_diffusivity:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       do iIonFluid=1,nIonFluid
          write(*,*)'e & ',NameFluid_I(iIonFluid+1),':'
          write(*,*)'s_ei  = ',eiSigma_I(iIonFluid)," [1/(Ohm*m)]"
       end do
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral),':'
          write(*,*)'s_en  = ',enSigma_I(iNeutral)," [1/(Ohm*m)]"
       end do
       write(*,*)'e & e:'
       write(*,*)'s_ee  = ',eeSigma," [1/(Ohm*m)]"
       write(*,*)''
       write(*,*)'Eta0   = ',Eta0Si," [m^2/s]"
       write(*,*)'Eta_en = ',sum(1/enSigma_I(1:nNeutral))/cMu," [m^2/s]"
       write(*,*)'Eta_ei = ',sum(1/eiSigma_I(1:nIonFluid))/cMu," [m^2/s]"
       write(*,*)'Eta_ee = ',1/eeSigma/cMu," [m^2/s]"
       write(*,*)'Eta_eX = ',EtaSiColl," [m^2/s]"
       write(*,*)'EtaTot = ',EtaSi," [m^2/s]"
       write(*,*)''
    end if

  end subroutine derive_cell_diffusivity

  !========================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModPhysics,     ONLY: No2Io_V, Io2No_V, No2Si_V, Si2No_V, &
         UnitN_, UnitTemperature_, UnitX_, UnitT_, UnitP_, ElectronPressureRatio
    use ModProcMH,      ONLY: iProc
    use ModMain,        ONLY: ProcTest, BlkTest, iTest, jTest, kTest, nBlockMax
    use ModAdvance,     ONLY: State_VGB
    use ModGeometry,    ONLY: Rmin_BLK, R_BLK
    use ModVarIndexes,  ONLY: Rho_, Pe_, P_
    use ModConst,       ONLY: cMu, cBoltzmann, cElectronMass, cElectronCharge
    use ModMultiFluid,  ONLY: MassIon_I
    use ModResistivity, ONLY: Eta0

    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) 

    integer :: i, j, k
    logical :: DoTest, DoTestMe=.true.
    real, dimension(1:nNeutral,MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: enSigma_IG
    real, dimension(1:nIonFluid,MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: eiSigma_IG, nIon_IG
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: Te_G, nElec_G
    real, dimension(1:nNeutral) :: fen_I
    real, dimension(1:nIonFluid) :: fei_I
    integer :: iIonFluid, iNeutral
    real :: EtaSi

    !---------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest) then
       call set_oktest('user_set_resistivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    ! nElec_G is the electron/ion density in SI units (n_e=n_itot)
    nElec_G = 0.
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       nIon_IG(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
       nElec_G(i,j,k) = sum(nIon_IG(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    if (UseElectronPressure) then
       Te_G = State_VGB(Pe_,:,:,:,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
            nElec_G)
    else
       Te_G(:,:,:) = State_VGB(P_,:,:,:,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec_G(:,:,:))
    end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

       call derive_cell_diffusivity(iBlock, i, j, k, Te_G(i,j,k), nIon_IG(1:nIonFluid,i,j,k), nElec_G(i,j,k), EtaSi)
       Eta_G(i,j,k) = EtaSi*SI2No_V(UnitX_)**2/SI2No_V(UnitT_)

    end do; end do; end do

    if(DoTestMe) then
       write(*,*)'user_set_resistivity:'
       write(*,*)'Te    = ',Te_G(iTest,jTest,kTest)," [K]"
       write(*,*)'n_e   = ',nElec_G(iTest,jTest,kTest)," [m^-3]"
       write(*,*)'Eta   = ',Eta_G(iTest,jTest,kTest)*No2SI_V(UnitX_)**2/No2SI_V(UnitT_)," [m^2/s]"
    end if

  end subroutine user_set_resistivity

  !========================================================================

  subroutine user_material_properties(State_V, i,j,k,iBlock,iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut,   &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W, &
       EntropyOut)

    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitX_, ElectronPressureRatio, inv_gm1
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, ExtraEInt_
    use ModConst,        ONLY: cElectronCharge, cBoltzmann, cMu, cElectronMass
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: ProcTest, iTest, jTest, kTest, BlkTest
    use ModResistivity,  ONLY: Eta0SI
    use ModProcMH,       ONLY: iProc
    use ModGeometry,     ONLY: Xyz_DGB
    use ModNumConst,     ONLY: cTiny

    !------------------------------------------------------------------------
    ! The State_V vector is in normalized units
    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in) :: i, j, k, iBlock, iDir
    real, optional, intent(in)  :: EinternalIn                  ! [J/m^3]
    real, optional, intent(in)  :: TeIn                         ! [K]
    real, optional, intent(out) :: NatomicOut                   ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut          ! dimensionless
    real, optional, intent(out) :: EinternalOut                 ! [J/m^3]
    real, optional, intent(out) :: TeOut                        ! [K]
    real, optional, intent(out) :: PressureOut                  ! [Pa]   
    real, optional, intent(out) :: CvOut                        ! [J/(K*m^3)]  
    real, optional, intent(out) :: GammaOut                     ! dimensionless
    real, optional, intent(out) :: HeatCondOut                  ! [W/(m*K)]   
    real, optional, intent(out) :: IonHeatCondOut               ! [W/(m*K)]
    real, optional, intent(out) :: TeTiRelaxOut                 ! [1/s]  
    real, optional, intent(out) :: OpacityPlanckOut_W(nWave)    ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave) ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)           ! [J/m^3] 
    real, optional, intent(out) :: EntropyOut

    real, save :: KappaCoeffSI = (cBoltzmann/cElectronCharge)**2/cMu
    real :: nElec, EtaSI, TeSI
    real, dimension(nIonFluid) :: nIon_I
    integer :: iIonFluid
    logical :: DoTest, DoTestMe=.true.

    real :: xmin, xmax, HeatCondFactor, widthmax, widthmin, xMaxyz, xMinyz

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('user_material_properties',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
    nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
    if (UseElectronPressure) then
       TeSI = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    else
       TeSI = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    end if
    if(present(CvOut)) CvOut = cBoltzmann*nElec*inv_gm1
    if(present(TeOut)) TeOut = TeSI
    if(present(AverageIonChargeOut).or.present(NatomicOut)) AverageIonChargeOut = nElec/sum(nIon_I(1:nIonFluid))
    if(present(NatomicOut)) NatomicOut = nElec/AverageIonChargeOut

    if(present(HeatCondOut)) then
       !!write(*,*)'iBlock = ',iBlock,'iNeutralBlockLast = ',iNeutralBlockLast

       xmin =  75000e3*Si2No_V(UnitX_) ! cometopause 75'000 km
       xmax = 100000e3*Si2No_V(UnitX_) ! cometopause max
       widthmin = 2.*xmin              ! flaring ratio 2 (Cravens 1989)
       widthmax = 2.*xmax              ! flaring ratio 2 (Cravens 1989)

       xMaxyz = -(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)/widthmax**2*xmax+xmax

       if(Xyz_DGB(x_,i,j,k,iBlock) > xMaxyz) then
          HeatCondFactor = 0.0
          HeatCondOut = 0.0
       else
          if (iBlock.ne.iNeutralBlockLast) then
             call user_neutral_atmosphere(iBlock)
          end if

          call derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I(1:nIonFluid), nElec, EtaSi)

          xMinyz = -(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)/widthmin**2*xmin+xmin
          if (Xyz_DGB(x_,i,j,k,iBlock) < xMinyz) then
             HeatCondFactor = 1.0
          else
             HeatCondFactor=(xMaxyz-Xyz_DGB(x_,i,j,k,iBlock))/(xMaxyz-xMinyz)
          end if
          HeatCondOut = TeSI/EtaSI*KappaCoeffSI*HeatCondFactor
       end if

       ! TestArray(1,i,j,k,iBlock) = nElec*sqrt(cBoltzmann*TeSi/cElectronMass)*cBoltzmann*TeSI/HeatCondOut
       ! TestArray(1,i,j,k,iBlock) = HeatCondOut
       ! TestArray(1,i,j,k,iBlock) = HeatCondFactor

    end if

    if(DoTestMe) then
       write(*,*)'user_material_properties:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       if(present(CvOut)) write(*,*)'Cv     = ',CvOut,' [J/(K*m^3)]'
       if(present(HeatCondOut)) then
          write(*,*)'Eta    = ',EtaSI," [m^2/s]"
          write(*,*)'Kappa  = ',HeatCondOut," [W/(m*K)]"
          write(*,*)'Ffree  = ',nElec*sqrt(cBoltzmann*TeSi/cElectronMass)*cBoltzmann*TeSI," [W/m^2]"
       end if
       write(*,*)''
    end if
  end subroutine user_material_properties

  !========================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !----------------------------------------------------------------------

    !! Source terms are evaluated explicitly!
    !RETURN

    ! All ion momenta are implicit
    if(UseElectronPressure)then
       allocate(iVarPointImpl_I(5*nIonFluid + 1))
       iVarPointImpl_I(5*nIonFluid + 1) = Pe_
    else
       allocate(iVarPointImpl_I(5*nIonFluid))
    end if

    do iFluid = 1, nIonFluid
       iVarPointImpl_I(5*iFluid-4) = iRhoIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-3) = iRhoUxIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-2) = iRhoUyIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-1) = iRhoUzIon_I(iFluid)
       iVarPointImpl_I(5*iFluid)   = iPIon_I(iFluid)
    end do

    IsPointImplMatrixSet = .false.
    !IsAsymmetric= .false.

  end subroutine user_init_point_implicit

  !========================================================================

  subroutine user_init_session
    use ModPhysics,    ONLY: ElectronPressureRatio

    !----------------------------------------------------------------------
    if (ElectronPressureRatio.le.0.) call stop_mpi('ERROR: Electron Pressure Ratio > 0 for init!')

  end subroutine user_init_session

  !========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional,&
       PlotVar_G, PlotVarBody, UsePlotVarBody,&
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitU_, UnitT_, &
         ElectronCharge, ElectronPressureRatio
    use ModVarIndexes, ONLY: Rho_, P_, Pe_
    use ModConst,      ONLY: cBoltzmann
    use ModCurrent,    ONLY: get_current
    use ModMultiFluid, ONLY: MassIon_I
    use ModMain,       ONLY: Dt_BLK

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k, iIonFluid
    real :: nElec
    real, dimension(3)           :: Current_I, uIonMean_I
    real, dimension(nIonFluid)   :: nIon_I
    real, dimension(3,nIonFluid) :: uIon_I


    !--------------------------------------------------------------------------

    IsFound = .true.

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    select case(NameVar)
    case('nn1')
       NameIdlUnit = '1/cm^3'
       NameTecUnit = '[1/cm^3]'
       !do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = NnNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6
       end do; end do; end do

    case('nn2')
       NameIdlUnit = '1/cm^3'
       NameTecUnit = '[1/cm^3]'
       !do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = NnNeutral_IG(2,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6
       end do; end do; end do

    case('unx1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[cm/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnxNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! x direction
       end do; end do; end do

    case('uny1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[cm/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnyNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! y direction
       end do; end do; end do

    case('unz1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnzNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! z direction
       end do; end do; end do

    case('tn1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = TnNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do; end do; end do

    case('te')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          if(UseElectronPressure)then
             PlotVar_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)/&
                  (cBoltzmann*nElec)
          else
             PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)*ElectronPressureRatio/&
                  (1.+ElectronPressureRatio)/(cBoltzmann*nElec)
          end if
       end do; end do; end do

    case('ti1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(iPIon_I(1),i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I(1),i,j,k,iBlock)/MassIon_I(1)*NO2SI_V(UnitN_))
       end do; end do; end do

    case('uex')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(1,1:nIonFluid) = State_VGB(iRhoUxIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(1) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(1) = uIonMean_I(1)+nIon_I(iIonFluid)* &
                  uIon_I(1,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(1)-Current_I(1)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uey')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(2,1:nIonFluid) = State_VGB(iRhoUyIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(2) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(2) = uIonMean_I(2)+nIon_I(iIonFluid)* &
                  uIon_I(2,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(2)-Current_I(2)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uez')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(3,1:nIonFluid) = State_VGB(iRhoUzIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(3) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(3) = uIonMean_I(3)+nIon_I(iIonFluid)* &
                  uIon_I(3,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(3)-Current_I(3)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('dt')
       NameIdlUnit = 's'
       NameTecUnit = '[s]'
       PlotVar_G(:,:,:) = Dt_BLK(iBlock)*No2SI_V(UnitT_)

    case('ns')
       NameIdlUnit = ' '   
       NameTecUnit = '[ ]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = ne20eV_GB(i,j,k,iBlock)
       end do; end do; end do
       ! case('testarray1')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(1,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray2')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(2,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray3')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(3,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray4')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(4,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('fluxlim')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = FluxLimited_GB(i,j,k,iBlock)
       !    end do; end do; end do
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var



  !========================================================================

  subroutine user_preset_conditions(i,j,k,iBlock)
    ! This is applied as initial conditions and in the upstream boundary for the semi-implicit heat conduction
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: SW_N, SW_Ux, SW_Uy, SW_Uz, SW_T_Dim, LowDensityRatio, ElectronPressureRatio, Io2No_V, &
         UnitTemperature_

    integer,intent(in) :: i, j, k, iBlock

    State_VGB(SwRho_,i,j,k,iBlock)     = SW_n*MassIon_I(SW_)
    State_VGB(SwRhoUx_,i,j,k,iBlock)   = SW_n*MassIon_I(SW_)*SW_Ux
    State_VGB(SwRhoUy_,i,j,k,iBlock)   = SW_n*MassIon_I(SW_)*SW_Uy
    State_VGB(SwRhoUz_,i,j,k,iBlock)   = SW_n*MassIon_I(SW_)*SW_Uz
    State_VGB(SwP_,i,j,k,iBlock)       = SW_n*SW_T_dim*Io2No_V(UnitTemperature_)

    State_VGB(HpRho_,i,j,k,iBlock)     = SW_n*LowDensityRatio*MassIon_I(Hp_)
    State_VGB(HpRhoUx_,i,j,k,iBlock)   = SW_n*LowDensityRatio*MassIon_I(Hp_)*SW_Ux
    State_VGB(HpRhoUy_,i,j,k,iBlock)   = SW_n*LowDensityRatio*MassIon_I(Hp_)*SW_Uy
    State_VGB(HpRhoUz_,i,j,k,iBlock)   = SW_n*LowDensityRatio*MassIon_I(Hp_)*SW_Uz
    State_VGB(HpP_,i,j,k,iBlock)       = SW_n*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_)

    State_VGB(H2OpRho_,i,j,k,iBlock)   = SW_n*LowDensityRatio*MassIon_I(H2Op_)
    State_VGB(H2OpRhoUx_,i,j,k,iBlock) = SW_n*LowDensityRatio*MassIon_I(H2Op_)*SW_Ux
    State_VGB(H2OpRhoUy_,i,j,k,iBlock) = SW_n*LowDensityRatio*MassIon_I(H2Op_)*SW_Uy
    State_VGB(H2OpRhoUz_,i,j,k,iBlock) = SW_n*LowDensityRatio*MassIon_I(H2Op_)*SW_Uz
    State_VGB(H2OpP_,i,j,k,iBlock)     = SW_n*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_)          

    State_VGB(Rho_,i,j,k,iBlock)       = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
    State_VGB(RhoUx_,i,j,k,iBlock)     = sum(State_VGB(iRhoUxIon_I,i,j,k,iBlock))
    State_VGB(RhoUy_,i,j,k,iBlock)     = sum(State_VGB(iRhoUyIon_I,i,j,k,iBlock))
    State_VGB(RhoUz_,i,j,k,iBlock)     = sum(State_VGB(iRhoUzIon_I,i,j,k,iBlock))

    if(UseElectronPressure) then
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))
       State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
    else
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))* &
            (1.+ElectronPressureRatio)
    end if

  end subroutine user_preset_conditions

  !========================================================================

  subroutine user_set_ICs(iBlock)
    use ModIO,       ONLY: restart
    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest, Body1_, Body1
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: ElectronPressureRatio, No2Si_V, UnitRho_, UnitRhoU_, UnitP_, rBody
    use ModConst,    ONLY: cBoltzmann
    use ModGeometry, ONLY: R_BLK

    integer, intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true.
    integer :: i, j, k, iIonFluid
    ! !-------------------------------------------------------------------------
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_ICs', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !if (iBlock.ne.iNeutralBlockLast) then
    !   call user_neutral_atmosphere(iBlock)
    !end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       if ((Body1) .and. (R_BLK(i,j,k,iBlock) < Rbody)) then
          ! State_VGB(:,i,j,k,iBlock) = CellState_VI(:,Body1_)
          ! State_VGB(iUx_I,i,j,k,iBlock) = 0.
          ! State_VGB(iUy_I,i,j,k,iBlock) = 0.
          ! State_VGB(iUz_I,i,j,k,iBlock) = 0.
          ! State_VGB(iPIon_I,i,j,k,iBlock) = BodyNDim_I*cBoltzmann*BodyTDim_I*1e6*SI2No_V(UnitP_)
          ! if (UseElectronPressure) then
          !    State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          !    State_VGB(Pe_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock)*&
          !         ElectronPressureRatio
          ! else
          !    State_VGB(P_,i,j,k,iBlock) = &
          !         sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
          ! end if
       else

          call user_preset_conditions(i,j,k,iBlock)

       end if
    end do; end do ; end do


    if(DoTestMe) then
       i=iTest ; j=jTest ; k=kTest
123    format (A13,ES25.16,A15)
       do iIonFluid=1,nIonFluid
          write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       end do
       write(*,*)''
       write(*,*)'Total:'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'rhoUx     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)/(1.+ElectronPressureRatio)* &
               No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if

    end if

  end subroutine user_set_ICs

  !========================================================================

  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,  ONLY: State_VGB
    use ModImplicit, ONLY: StateSemi_VGB, iTeImpl
    use ModSize,     ONLY: nI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModPhysics,  ONLY: Si2No_V, UnitTemperature_, SW_n, SW_Ux, SW_Uy, SW_Uz, &
         LowDensityRatio, ElectronPressureRatio, SW_T_dim, Io2No_V

    integer,          intent(in)  :: iBlock, iSide
    character(len=*),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer:: i, j, k
    real :: TeSi

    character(len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !-------------------------------------------------------------------
    IsFound = .true.

    if(TypeBc == 'usersemi')then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI

          call user_preset_conditions(i,j,k,iBlock)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, TeOut=TeSi)
          StateSemi_VGB(iTeImpl,i,j,k,iBlock) = TeSi*Si2No_V(UnitTemperature_)


       end do; end do; end do

       RETURN
    elseif(TypeBc == 'usersemilinear')then
       RETURN
    end if


  end subroutine user_set_cell_boundary

  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,       ONLY: Dt_BLK, BLKtest, ProcTest
    use ModPhysics,    ONLY: No2SI_V, UnitT_
    use ModProcMH,     ONLY: iProc

    real, intent(out)             :: VarValue
    character (len=*), intent(in) :: TypeVar
    real, intent(in), optional    :: Radius

    character (len=*), parameter  :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------------
    select case(TypeVar)
    case('dtpnt')
       ! local time step
       VarValue = 0.
       if (iProc == ProcTest) VarValue = Dt_BLK(BLKtest)*No2SI_V(UnitT_)

    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  !============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModSize,         ONLY: x_
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_
    use ModFaceBoundary, ONLY: TimeBc, iFace, jFace, kFace, FaceCoords_D, &
         iBoundary, VarsTrueFace_V
    use ModSolarwind,    ONLY: get_solar_wind_point
    !    use ModB0,           ONLY: B0_DX
    use ModMain,         ONLY: body1_
    use ModPhysics,      ONLY: LowDensityRatio, SW_Ux, SW_Uy, SW_Uz, SW_n, SW_T_dim, &
         ElectronPressureRatio, UnitTemperature_, Io2No_V, SW_Bx, SW_By, SW_Bz, &
         NO2SI_V, UnitP_, UnitRho_, UnitRhoU_, UnitU_, UnitB_, UnitN_, Io2SI_V, &
         BodyRho_I, BodyP_I

    logical :: FirstCall = .true., DoTest, DoTestMe=.true.
    integer :: iIonFluid
    real    :: UdotR(nIonFluid), URefl_D(1:3,nIonFluid)
    !real    :: BdotR, BRefl_D(1:3)

    real, intent(out):: VarsGhostFace_V(nVar)

    !------------------------------------------------------------------------

    if(DoTestMe.and.FirstCall) then
       call set_oktest('user_set_face_boundary',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !! Outer boundaries
    if(iBoundary >= 0) then

       call get_solar_wind_point(TimeBc, FaceCoords_D(x_),VarsGhostFace_V)

       ! Solar wind protons     
       VarsGhostFace_V(SwRho_)     = SW_n*MassIon_I(SW_)*Io2No_V(UnitRho_)
       VarsGhostFace_V(SwRhoUx_)   = SW_Ux
       VarsGhostFace_V(SwRhoUy_)   = SW_Uy
       VarsGhostFace_V(SwRhoUz_)   = SW_Uz
       VarsGhostFace_V(SwP_)       = SW_n*Io2No_V(UnitN_)*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)
       ! Cometary protons
       VarsGhostFace_V(HpRho_)     = SW_n*MassIon_I(Hp_)*LowDensityRatio*Io2No_V(UnitRho_)
       VarsGhostFace_V(HpRhoUx_)   = SW_Ux
       VarsGhostFace_V(HpRhoUy_)   = SW_Uy
       VarsGhostFace_V(HpRhoUz_)   = SW_Uz
       VarsGhostFace_V(HpP_)       = SW_n*Io2No_V(UnitN_)*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)

       ! Cometary heavies
       VarsGhostFace_V(H2OpRho_)   = SW_n*MassIon_I(H2Op_)*LowDensityRatio*Io2No_V(UnitRho_)
       VarsGhostFace_V(H2OpRhoUx_) = SW_Ux
       VarsGhostFace_V(H2OpRhoUy_) = SW_Uy
       VarsGhostFace_V(H2OpRhoUz_) = SW_Uz
       VarsGhostFace_V(H2OpP_)     = SW_n*Io2No_V(UnitN_)*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)

       VarsGhostFace_V(Rho_)       = sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(RhoUx_)     = SW_Ux
       VarsGhostFace_V(RhoUy_)     = SW_Uy
       VarsGhostFace_V(RhoUz_)     = SW_Uz

       ! VarsGhostFace_V(Bx_:Bz_)    = VarsGhostFace_V(Bx_:Bz_) - B0_DX(:,iFace, jFace, kFace)
       ! VarsGhostFace_V(Bx_)        = SW_Bx - B0_DX(1,iFace, jFace, kFace)
       ! VarsGhostFace_V(By_)        = SW_By - B0_DX(2,iFace, jFace, kFace)
       ! VarsGhostFace_V(Bz_)        = SW_Bz - B0_DX(3,iFace, jFace, kFace)

       if(UseElectronPressure) then
          VarsGhostFace_V(P_)      = sum(VarsGhostFace_V(iPIon_I))
          VarsGhostFace_V(Pe_)     = VarsGhostFace_V(P_)*ElectronPressureRatio
       else
          VarsGhostFace_V(P_)      = sum(VarsGhostFace_V(iPIon_I))*(1.+ElectronPressureRatio)
       end if

       !! Body boundaries
    else if (iBoundary <= body1_) then

       !! Projection length of U_ and B_ on the local surface radius vector
       !! in units of the surface radius vector [rBody]

       do iIonFluid=1,nIonFluid
          UdotR(iIonFluid) = dot_product(VarsTrueFace_V(iRhoUx_I(iIonFluid):iRhoUz_I(iIonFluid)), &
               FaceCoords_D)/dot_product(FaceCoords_D,FaceCoords_D)
          !! Projection vectors
          URefl_D(1:3,iIonFluid) = UdotR(iIonFluid)*FaceCoords_D(1:3)      
       end do
       !BdotR = dot_product(VarsTrueFace_V(Bx_:Bz_),FaceCoords_D)/ &
       !     dot_product(FaceCoords_D,FaceCoords_D)


       !! Projection vectors
       !BRefl_D = BdotR*FaceCoords_D

       !! Floating boundary conditions allowing inflow
       VarsGhostFace_V = VarsTrueFace_V

       !! Bz component propagated through moon, Bx and By didn't
       !  VarsGhostFace_V(Bx_:By_) = 0.0
       !  VarsGhostFace_V(Bz_)     = SW_Bz

       !! set outward flux body value (Comet's surface not considered as plasma source)
       !! leave inward flux untouched
       do iIonFluid=1,nIonFluid
          if (UdotR(iIonFluid) > 0.0) then
             VarsGhostFace_V(iUx_I(iIonFluid):iUz_I(iIonFluid)) = 0.0
             VarsGhostFace_V(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid)) = 0.0
             VarsGhostFace_V(iRhoIon_I(iIonFluid)) = BodyRho_I(iIonFluid)
             VarsGhostFace_V(iPIon_I(iIonFluid)) = BodyP_I(iIonFluid)
          endif
          VarsGhostFace_V(Rho_)   = sum(VarsGhostFace_V(iRhoIon_I))
          VarsGhostFace_V(RhoUx_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUxIon_I))/ &
               sum(VarsGhostFace_V(iRhoIon_I))
          VarsGhostFace_V(RhoUy_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUyIon_I))/ &
               sum(VarsGhostFace_V(iRhoIon_I))
          VarsGhostFace_V(RhoUz_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUzIon_I))/ &
               sum(VarsGhostFace_V(iRhoIon_I))
          VarsGhostFace_V(P_)     = sum(VarsGhostFace_V(iPIon_I))
       end do

    end if

    if(DoTestMe) then
       !FirstCall = .false.
       write(*,*)'Boundary No =',iBoundary
       write(*,*)'VarsGhostFaces'
       do iIonFluid=1,nIonFluid
          write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
          write(*,*)'VarsGhostFaceRho   = ',VarsGhostFace_V(iRhoIon_I(iIonFluid))*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,*)'VarsGhostFaceUx    = ',VarsGhostFace_V(iRhoUxIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceUy    = ',VarsGhostFace_V(iRhoUyIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceUz    = ',VarsGhostFace_V(iRhoUzIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceP     = ',VarsGhostFace_V(iPIon_I(iIonFluid))*No2SI_V(UnitP_)," [Pa]"
       end do
       write(*,*)''
       write(*,*)'Total ion fluid:'
       write(*,*)'VarsGhostFaceRho   = ',sum(VarsGhostFace_V(iRhoIon_I)*No2SI_V(UnitRho_))," [kg/m^3]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUxIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUyIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUzIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceP   = ',VarsGhostFace_V(P_)*No2SI_V(UnitP_)," [Pa]"
       if (UseElectronPressure) then
          write(*,*) ''
          write(*,*)'VarsGhostFacePe  = ',VarsGhostFace_V(Pe_)*No2SI_V(UnitP_)," [Pa]"
       end if
       write(*,*)''
       write(*,*)'Magnetic field:'
       write(*,*)'VarsGhostFaceBx    = ',VarsGhostFace_V(Bx_)*No2SI_V(UnitB_)," [T]"
       write(*,*)'VarsGhostFaceBy    = ',VarsGhostFace_V(By_)*No2SI_V(UnitB_)," [T]"
       write(*,*)'VarsGhostFaceBz    = ',VarsGhostFace_V(Bz_)*No2SI_V(UnitB_)," [T]"
       write(*,*)''
    end if

  end subroutine user_set_face_boundary

  !============================================================================

  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, x2

    implicit none

    integer, intent(in):: iBlock

    character (len=*), parameter :: Name='user_set_boundary_cells'

    !--------------------------------------------------------------------------
    ! For inflow in negative x direction
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = Xyz_DGB(x_,:,:,:,iBlock) > x2

  end subroutine user_set_boundary_cells

end module ModUser
