!#NOTPUBLIC  email:rubinmar@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!========================================================================

module ModUser

  use ModSize
  use ModVarIndexes, ONLY: nVar
  use ModAdvance,    ONLY: Pe_, UseElectronPressure
  use ModUserEmpty,                              &
       IMPLEMENTED1  => user_read_inputs,         &
       IMPLEMENTED2  => user_calc_sources,        &
       IMPLEMENTED3  => user_update_states,       &
       IMPLEMENTED4  => user_set_face_boundary,   &
       IMPLEMENTED5  => user_set_resistivity,     &
       IMPLEMENTED6  => user_material_properties, &
       IMPLEMENTED7  => user_init_point_implicit, &
       IMPLEMENTED8  => user_init_session,        &
       IMPLEMENTED9  => user_set_plot_var,        &
       IMPLEMANTED10 => user_set_ICs,             &
       IMPLEMANTED11 => user_get_log_var,         & 
       IMPLEMENTED12 => user_set_boundary_cells,  &
       IMPLEMENTED13 => user_set_cell_boundary


  use ModMultiFluid

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.1
  character (len=*), parameter :: NameUserModule = &
       'Rubin, 3-fluid Comet MHD module, Dec 2011'

  integer, parameter, public :: nNeutral = 2 !! number of neutral species
  !! Neutral species names
  character (len=6), parameter, public :: NameNeutral_I(nNeutral) = &
       (/ 'H2O  ', 'H    ' /)
  integer, parameter, public :: H2O_  =  1
  integer, parameter, public :: H_    =  2
  !! Ion species names
  integer, parameter, public :: SWp_  =  1
  integer, parameter, public :: Hp_   =  2
  integer, parameter, public :: H2Op_ =  3

  real, dimension(nNeutral,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1) :: NnNeutral_IG, &
       UnxNeutral_IG, UnyNeutral_IG, UnzNeutral_IG, TnNeutral_IG
  real, dimension(1:nNeutral) :: NeutralMass_I(nNeutral)
  integer :: iNeutralBlockLast = -1
  real :: Qprod, Tmin, rHelio
  !! ??? testing
  ! real, dimension(4,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: TestArray = 0.

contains

  !========================================================================
  subroutine user_read_inputs

    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
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
          call read_var('Qprod' , Qprod)            !! Neutral gas production rate [1/s]
          call read_var('rHelio', rHelio)           !! Heliocentric distance [AU]
          call read_var('Tmin' , Tmin)              !! Minimum ion temperature (enforced in update states)
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
    use ModPhysics
    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number 
    use ModProcMH,     ONLY: iProc 
    use ModGeometry,   ONLY: R_BLK, Xyz_DGB
    use ModMultiFluid, ONLY: MassIon_I
    use ModIO,         ONLY: iUnitOut

    integer,intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true., init=.true.
    integer :: i, j, k, iNeutral
    real, dimension(nNeutral,nNeutral) :: DissocRate_II
    real, dimension(nNeutral) :: DestructRate_I, uNeutr_I

    !! Neutral dissociation/destruction rates
    DissocRate_II = 0.0 ; DestructRate_I = 0.0 
    DestructRate_I(H2O_) = 1.00E-6 !! Gombosi et al., J. Geophys. Res., (1996) 

    !----------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest) then
       call set_oktest('user_neutral_atmosphere',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    iNeutralBlockLast = iBlock

    if (.not.use_block_data(iBlock)) then

       if(iProc==0.and.init) then
          !! Init, placeholder
          init=.false.
       end if

       uNeutr_I(H2O_) = 800.  ! [m/s]
       uNeutr_I(H_)   = 8000. ! [m/s]

       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

          ! Haser model for H2O
          NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1) = Qprod/(4.*cPi*(R_BLK(i,j,k,iBlock)*rPlanetSI)**2&
               *uNeutr_I(H2O_))*exp(-DestructRate_I(H2O_)/uNeutr_I(H2O_)*R_BLK(i,j,k,iBlock)*rPlanetSI)

          ! H (after Combi 1996 with [cm^-3] -> [m^-3] and r[m] -> r[cm])
          NnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1) = Qprod*1.205E-11*1e6*(1E2*R_BLK(i,j,k,iBlock)*rPlanetSI)**(-1.6103)

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


          TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)  =  25.    !! estimate
          TnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1)    =  1000.  !! estimate
       end do;  end do;  end do

       NeutralMass_I(H2O_) = 17.*cProtonMass
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
  subroutine calc_electron_collision_rates(Te,i,j,k,iBlock,fen_I,fei_I)

    ! calculate all collision rates involving electrons 
    ! (used for sources & resistivity)

    use ModAdvance,    ONLY: State_VGB, Rho_
    use ModPhysics,    ONLY: No2SI_V, UnitN_
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I

    integer,intent(in) :: i,j,k,iBlock   
    real,intent(in)    :: Te
    real,intent(out)   :: fen_I(nNeutral)
    real,intent(out)   :: fei_I(nIonFluid)
    real :: sqrtTe
    !----------------------------------------------------------------------

    !! electron-neutral and electron-ion collision rates
    !! provide all rates in SI units

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
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)
    !! e - Hp, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fei_I(Hp_) = 54.5*ChargeIon_I(Hp_)**2*State_VGB(HpRho_,i,j,k,iBlock)/MassIon_I(Hp_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)
    fei_I(SWp_) = fei_I(Hp_)

  end subroutine calc_electron_collision_rates

  !========================================================================
  subroutine user_calc_rates(Ti_I,Te,i,j,k,iBlock,nElec,nIon_I,fin_II,fii_II,alpha_I,kin_IIII,v_II,uElec_D,uIon_DI,Qexc_II)

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
    real,intent(out)   :: alpha_I(nIonFluid)
    real,intent(out)   :: kin_IIII(nIonFluid,nNeutral,nNeutral,nIonFluid)
    real,intent(out)   :: v_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qexc_II(nNeutral,nIonFluid)


    real :: Tr, ueBulk2, ueTherm2, Ee, A(7), uiBulk2(nIonFluid), uiTherm2(nIonFluid), Tred, Mred
    real :: Dist, NCol, sigma, J3, uNeutr
    real,dimension(nNeutral,nIonFluid) :: sigma_e
    integer :: n

    !----------------------------------------------------------------------

    !! provide all rates in SI units

    !! initialize all collision rates with zero
    fin_II = 0. ; fii_II = 0. ; kin_IIII = 0. ; alpha_I = 0. ; v_II = 0.; sigma_e = 0. ; Qexc_II = 0.

    !! Ionization rates (add opacity correction/shadowing when needed)
    v_II(H2O_,Hp_)   = (1.30E-8+3.26E-8)/(rHelio**2) !! Ionization rate producing Hp from H2O & OH [1/s] (Huebner et al. 1992 & estimate)
    v_II(H_,Hp_)     = 7.30E-8/(rHelio**2) !! Ionization rate producing Hp from H [1/s] (Huebner et al. 1992)
    v_II(H2O_,H2Op_) = 1.00E-6 !! Gombosi et al., J. Geophys. Res., (1996) 
    v_II(H_,SWp_)    = v_II(H_,Hp_)*1e-9 ! set production of solar wind protons to a small but non-zero value

    !! Electron excess energies from ionization (increases electron pressure)
    Qexc_II(H2O_,Hp_)   = 4.0054E-18 ! 25.0 eV, Huebner 1992
    Qexc_II(H_,Hp_)     = 5.6076E-19 !  3.5 eV, Huebner 1992
    Qexc_II(H2O_,H2Op_) = 1.9226E-18 ! 12.0 eV, Huebner 1992
    
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

    !! electron impact ionization rate
    ueBulk2  = sum(uElec_D(:)**2)                             ! Electron bulk speed qubed in [m^2/s^2]
    ueTherm2 = 3.*cBoltzmann*Te/(cElectronMass)               ! Electron thermal speed qubed in [m^2/s^2]
    Ee = 0.5*cElectronMass*(ueBulk2+ueTherm2)/cElectronCharge ! Electron energy in [eV]
    !! electron impact cross section fitted according to 
    ! Talukder et al., Empirical model for electron impact ionization cross sections of neutral atoms, The European Physical Journal D
    if (Ee > 13.) then
       A = (/ 13.3569, 2457.69, -1910.67, -5405.82, 17242.7, -27453.3, 15429.2 /) ! fit parameters
       sigma_e(H2O_,H2Op_) = A(2)*log(Ee/A(1))
       do n=1,5
          sigma_e(H2O_,H2Op_) = sigma_e(H2O_,H2Op_)+A(n+2)*(1-A(1)/Ee)**n
       end do
       sigma_e(H2O_,H2Op_) = 1e-20*1/(Ee*A(1))*sigma_e(H2O_,H2Op_)                ! electron impact cross section in [m^2]
       if (sigma_e(H2O_,H2Op_).gt.0.) then
          uiBulk2(H2Op_)  = sum(uIon_DI(:,H2Op_)**2)                                 ! Water bulk speed qubed in [m^2/s^2]
          uiTherm2(H2Op_) = 3.*cBoltzmann*Ti_I(H2Op_)/(MassIon_I(H2Op_)*cProtonMass) ! Water thermal speed qubed in [m^2/s^2]
          !! Gombosi book f12=n2*sigma12*sqrt(v1**2+v2**2)                           ! add electron impact ionization rate [1/s]
          v_II(H2O_,H2Op_) = v_II(H2O_,H2Op_) + nElec*sigma_e(H2O_,H2Op_)*sqrt(ueBulk2+ueTherm2+uiBulk2(H2Op_)+uiTherm2(H2Op_))
       end if
    end if

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
    Tr=(Ti_I(SWp_)+TnNeutral_IG(H_,i-MinI+1,j-MinJ+1,k-MinK+1))/2.
    kin_IIII(SWp_,H_,H_,Hp_) = 2.65E-10/1E6*sqrt(Tr)*(1.0-0.083*log10(Tr))**2  !! rate in [m^3/s]
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
    kin_IIII(Hp_,H2O_,H_,H2Op_)  = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated
    kin_IIII(SWp_,H2O_,H_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated

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
    Tred = (MassIon_I(H2Op_)*Ti_I(SWp_)+MassIon_I(SWp_)*Ti_I(H2Op_))/(MassIon_I(H2Op_)+MassIon_I(SWp_)) ! reduced temp
    Mred = MassIon_I(H2Op_)*MassIon_I(SWp_)/(MassIon_I(H2Op_)+MassIon_I(SWp_)) ! reduced mass
    fii_II(H2Op_,SWp_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SWp_)**2/MassIon_I(H2Op_)*&
         sqrt(Mred)*1e-6*nIon_I(SWp_)/(Tred*sqrt(Tred))
    ! SWp - H2Op, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(SWp_,H2Op_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SWp_)**2/MassIon_I(SWp_)*&
         sqrt(Mred)*1e-6*nIon_I(H2Op_)/(Tred*sqrt(Tred))

    ! SWp - Hp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tred = (MassIon_I(SWp_)*Ti_I(Hp_)+MassIon_I(Hp_)*Ti_I(SWp_))/(MassIon_I(SWp_)+MassIon_I(Hp_)) ! reduced temp
    Mred = MassIon_I(SWp_)*MassIon_I(Hp_)/(MassIon_I(SWp_)+MassIon_I(Hp_)) ! reduced mass
    fii_II(SWp_,Hp_) = 1.27*ChargeIon_I(SWp_)**2*ChargeIon_I(Hp_)**2/MassIon_I(SWp_)*&
         sqrt(Mred)*1e-6*nIon_I(Hp_)/(Tred*sqrt(Tred))
    ! Hp - SWp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(Hp_,SWp_) = 1.27*ChargeIon_I(SWp_)**2*ChargeIon_I(Hp_)**2/MassIon_I(Hp_)*&
         sqrt(Mred)*1e-6*nIon_I(SWp_)/(Tred*sqrt(Tred))

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
    alpha_I(SWp_)  = alpha_I(Hp_)
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
    use ModPhysics
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
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: Ti_IC, uIonElec2_IC, fei_IC, nIon_IC, SRho_IC, &
         SRhoUx_IC, SRhoUy_IC, SRhoUz_IC, SP_IC
    real, dimension(1:nNeutral,1:nIonFluid) :: Qexc_II
    real, dimension(1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: v_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: alpha_IC
    real, dimension(1:nIonFluid) :: fiiTot_I, finTot_I, vAdd_I, kinAdd_I, kinSub_I
    real, dimension(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: kin_IIIIC

    logical :: DoTest, DoTestMe=.true.
    real :: theta, fenTot, feiTot,logTe
    integer :: i,j,k,iNeutral,jNeutral,iIonFluid,jIonFluid,iTerm,iDim


    !----------------------------------------------------------------------

    ! Do not evaluate any source terms explicitly when running pointimplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Evaluate source terms explicitly even when running pointimplicit
    !if(UsePointImplicit .and. IsPointImplSource) RETURN

    !! Limit region for evaluation for source term evaluation
    !! if(RMin_BLK(iBlock) > 2.) RETURN

    if(iBlock == BlkTest) then
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


          !???
          ! do iIonFluid=1,nIonFluid
          !    if ((State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock).le.0.).or.(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock).le.0.)) then
          !       write(*,*)'Fluid = ',NameFluid_I(iIonFluid+1),'Pi = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock),&
          !            'Rhoi = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
          !       write(*,*)'x = ',Xyz_DGB(x_,i,j,k,iBlock),'y = ',Xyz_DGB(y_,i,j,k,iBlock),'z = ',Xyz_DGB(z_,i,j,k,iBlock),&
          !            'r   = ',r_BLK(i,j,k,iBlock)
          !       write(*,*)'i = ',i,'j = ',j,'k = ',k,'iBlock = ',iBlock
          !    end if
          ! end do
          ! if (State_VGB(Pe_,i,j,k,iBlock).le.0.) then
          !       write(*,*)'Pe = ',State_VGB(Pe_,i,j,k,iBlock)
          !       write(*,*)'x = ',Xyz_DGB(x_,i,j,k,iBlock),'y = ',Xyz_DGB(y_,i,j,k,iBlock),'z = ',Xyz_DGB(z_,i,j,k,iBlock),&
          !            'r   = ',r_BLK(i,j,k,iBlock)
          !       write(*,*)'i = ',i,'j = ',j,'k = ',k,'iBlock = ',iBlock
          ! end if


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
       ! No need to evaluate source terms for cells between conducting core and surface
       !!if((R_BLK(i,j,k,iBlock) < PlanetRadius).and.(UseResistivePlanet)) CYCLE

       call get_current(i,j,k,iBlock,Current_DC(:,i,j,k))

       ! calculate uElec_DC from Hall velocity -J/(e*n) [m/s]
       uElec_DC(1:3,i,j,k) = uIonMean_DC(1:3,i,j,k)-Current_DC(1:3,i,j,k)/(nElec_C(i,j,k)*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_)

       call calc_electron_collision_rates(Te_C(i,j,k),i,j,k,iBlock,fen_IC(1:nNeutral,i,j,k),fei_IC(1:nIonFluid,i,j,k))
       call user_calc_rates(Ti_IC(1:nIonFluid,i,j,k),Te_C(i,j,k),i,j,k,iBlock,nElec_C(i,j,k),nIon_IC(1:nIonFluid,i,j,k),&
            fin_IIC(1:nIonFluid,1:nNeutral,i,j,k),fii_IIC(1:nIonFluid,1:nIonFluid,i,j,k),alpha_IC(1:nIonFluid,i,j,k),&
            kin_IIIIC(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,i,j,k),v_IIC(1:nNeutral,1:nIonFluid,i,j,k),uElec_DC(1:3,i,j,k),&
            uIon_DIC(1:3,1:nIonFluid,i,j,k),Qexc_II(1:nNeutral,1:nIonFluid))

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
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+MassIon_I(iIonFluid)*fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&   !! ion-ion collisions
               (uIon_DIC(1,iIonFluid,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+NeutralMass_I(iNeutral)*fin_IIC(1:nIonFluid,iNeutral,i,j,k)*& !! ion-neutral collisions
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do

       !Current_DC(1,i,j,k) = 0.0001*Si2No_V(UnitJ_)
       !Current_DC(2,i,j,k) = 1*Si2No_V(UnitJ_)
       !Current_DC(3,i,j,k) = 1*Si2No_V(UnitJ_)


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
       SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(1,i,j,k)
       SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(1,1:nIonFluid,i,j,k)-uIonMean_DC(1,i,j,k))*Si2No_V(UnitU_)
       SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass

       !if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
       !   write(*,*)'un-ui     = ',(UnxNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)- uIon_DIC(1,1,i,j,k))
       !   write(*,*)'finTot = ',finTot_I(1)
       !   write(*,*)'ni    = ',nIon_IC(1,i,j,k)
       !end if
       SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)

       !! First moment, y component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+MassIon_I(iIonFluid)*fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&   !! ion-ion collisions
               (uIon_DIC(2,iIonFluid,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+NeutralMass_I(iNeutral)*fin_IIC(1:nIonFluid,iNeutral,i,j,k)*& !! ion-neutral collisions
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
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                               !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUyTerm_IIC(1,iIonFluid,i,j,k) = SRhoUyTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(2,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(2,i,j,k)
       SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(2,1:nIonFluid,i,j,k)-uIonMean_DC(2,i,j,k))*Si2No_V(UnitU_)
       SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass
       SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)

       !! First moment, z component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+MassIon_I(iIonFluid)*fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&   !! ion-ion collisions
               (uIon_DIC(3,iIonFluid,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+NeutralMass_I(iNeutral)*fin_IIC(1:nIonFluid,iNeutral,i,j,k)*& !! ion-neutral collisions
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
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                              !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUzTerm_IIC(1,iIonFluid,i,j,k) = SRhoUzTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(3,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &            !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(3,i,j,k)
       SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k) = -fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(3,1:nIonFluid,i,j,k)-uIonMean_DC(3,i,j,k))*Si2No_V(UnitU_)
       SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &   !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass
       SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &   !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)


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

       !alpha_IC(1:nIonFluid,i,j,k) = 0.

       !! Second moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       !!kinAdd_I = 0. ; 
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

       ! if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
       !    write(*,*)'Ti-Ti     = ',(Ti_IC(2,i,j,k)-Ti_IC(1,i,j,k))
       !    write(*,*)'fiiTot = ',fiiTot_I(1)
       !    write(*,*)'ni    = ',nIon_IC(1,i,j,k)
       ! end if

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+NeutralMass_I(iNeutral)/cProtonMass)*&
               cBoltzmann*(TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(3,1:nIonFluid,i,j,k) = 2.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)

       !if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
       !   write(*,*)'un-ui     = ',(UnxNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)- uIon_DIC(1,1,i,j,k))
       !   write(*,*)'finTot = ',finTot_I(1)
       !   write(*,*)'ni    = ',nIon_IC(1,i,j,k)
       !end if

       SPTerm_IIC(4,1:nIonFluid,i,j,k) = 2.*fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*nIon_IC(1:nIonFluid,i,j,k)*&
            cBoltzmann*(Te_C(i,j,k)-Ti_IC(1:nIonFluid,i,j,k))*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(5,1:nIonFluid,i,j,k) = 2./3.*fei_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*cElectronMass*&             !! ion-electron collisional exchange (due to Hall velocity)
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


       !kinAdd_I = 0.
       !vAdd_I = 0.
       ! if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
       !    write(*,*)'un-ui     = ',(UnxNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)- uIon_DIC(1,1,i,j,k))
       !    write(*,*)'kinAdd_I  = ',kinAdd_I(2)
       ! end if


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

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid)* &
                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*Qexc_II(iNeutral,1:nIonFluid)
          end do
          SPeTerm_IC(7,i,j,k) = 2./3.*sum(vAdd_I)*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                                !! heating of electrons due to ionization excess energy.

          ! TestArray(1,i,j,k,iBlock) = Te_C(i,j,k)
          ! TestArray(2,i,j,k,iBlock) = nElec_C(i,j,k)
          logTe = log(Te_C(i,j,k))
          SPeTerm_IC(8,i,j,k) = exp(-188.4701+33.2547*logTe-2.0792*logTe**2+0.0425*logTe**3)                              !! electron cooling due to collisions w/ water vapor
          if(Te_C(i,j,k)<1.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)) then
             SPeTerm_IC(8,i,j,k)=4.5e-9/(0.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))* &
                  (Te_C(i,j,k)-TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))
          else
             SPeTerm_IC(8,i,j,k)=SPeTerm_IC(8,i,j,k)+4.5e-9
          end if
          ! TestArray(3,i,j,k,iBlock) = SPeTerm_IC(8,i,j,k)
          SPeTerm_IC(8,i,j,k) = -2./3.*NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*nElec_C(i,j,k)* &
               SPeTerm_IC(8,i,j,k)/1e6*1.60217733e-19*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                           
          ! TestArray(4,i,j,k,iBlock) = SPeTerm_IC(8,i,j,k)
          ! SPeTerm_IC(8,i,j,k) = 0.


          !! Test capping heating in ramp up phase ???
          ! if (SPeTerm_IC(7,i,j,k) > -SPeTerm_IC(8,i,j,k)) SPeTerm_IC(7,i,j,k) =- SPeTerm_IC(8,i,j,k)
          ! if (SPeTerm_IC(7,i,j,k) < 0.) SPeTerm_IC(7,i,j,k) = SPeTerm_IC(8,i,j,k)

          !if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest) then
          ! write(*,*)'vAdd_I = ',sum(vAdd_I)
          ! write(*,*)'Si2No_V(UnitEnergyDens_) = ',Si2No_V(UnitEnergyDens_)
          ! write(*,*)'Si2No_V(UnitT_) = ',Si2No_V(UnitT_)
          ! write(*,*)'NnNeutral_I = ', NnNeutral_IG(:,i-MinI+1,j-MinJ+1,k-MinK+1)
          ! write(*,*)'Sp3  = ',SPTerm_IIC(3,:,i,j,k)
          ! write(*,*)'Sp7  = ',SPTerm_IIC(7,:,i,j,k)
          ! write(*,*)'Spe7 = ',SPeTerm_IC(7,i,j,k)
          ! write(*,*)'Spe8 = ',SPeTerm_IC(8,i,j,k)
          !  write(*,*)'Test = ',TestArray(i,j,k,iBlock)
          !end if

       end if

       !! sum up individual terms
       do iTerm=1,4
          SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(4,1:nIonFluid,i,j,k)
       ! do iTerm=1,5
       !    SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       !    SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       !    SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       ! end do
       SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k)
       SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k)
       SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k)
       SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k)
       SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k)
       SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k)
       SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k)
       SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k)
       SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k)
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
124       format (A17,I2,A3,A7,A3,F5.1,A5,F5.1,A3)
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
             write(*,*)'lr  = ',exp(-188.4701+33.2547*log(Te_C(i,j,k))-2.0792*(log(Te_C(i,j,k)))**2+0.0425*(log(Te_C(i,j,k)))**3)
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
             write(*,123)' v_io     = ',v_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
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
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB
    use ModPhysics
    use ModEnergy
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
          ! set minimum mass density
          if(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) < SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2) then
             State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) = SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2
             State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
          end if
       end do
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
       call calc_electron_collision_rates(Te_G(i,j,k),i,j,k,iBlock,fen_I(1:nNeutral),fei_I(1:nIonFluid))


       ! classical conductivities due to ion-electron and electron-neutral collisions (1E-20 to avoid div by 0)
       eiSigma_IG(1:nIonFluid,i,j,k) = cElectronCharge**2*nElec_G(i,j,k)/((fei_I(1:nIonFluid)+1E-20)*cElectronMass) 
       enSigma_IG(1:nNeutral,i,j,k) = cElectronCharge**2*nElec_G(i,j,k)/((fen_I(1:nNeutral)+1E-20)*cElectronMass)

       !! Eta_G is calculated from both conductivities using Kirchhoff's rule:
       !! 1/sigma_tot = 1/eiSigma_IG+1/enSigma_IG
       !! The resulting conductivity is close to Spitzer conductivity far from the comet and
       !! decreases due to abundant electron-neutral collisions close to the nucleus

       !! Eta_G = 1/(sigma_tot*mu_0) magnetic diffusivity
       Eta_G(i,j,k) = Eta0 + (sum(1/eiSigma_IG(1:nIonFluid,i,j,k))+sum(1/enSigma_IG(1:nNeutral,i,j,k)))*&
            SI2No_V(UnitX_)**2/SI2No_V(UnitT_)/cMu

    end do; end do; end do

    if(DoTestMe) then
       write(*,*)'Te    = ',Te_G(iTest,jTest,kTest)," [K]"
       write(*,*)'n_e   = ',nElec_G(iTest,jTest,kTest)," [m^-3]"
       do iIonFluid=1,nIonFluid
          write(*,*)'e & ',NameFluid_I(iIonFluid+1),':'
          write(*,*)'s_ei  = ',eiSigma_IG(iIonFluid,iTest,jTest,kTest)," [1/(Ohm*m)]"
       end do
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral),':'
          write(*,*)'s_en  = ',enSigma_IG(iNeutral,iTest,jTest,kTest)," [1/(Ohm*m)]"
       end do
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

    use ModPhysics,     ONLY: No2Si_V, UnitP_, UnitN_, ElectronPressureRatio, inv_gm1
    use ModVarIndexes,  ONLY: nVar, Rho_, p_, ExtraEInt_
    use ModConst,       ONLY: cElectronCharge, cBoltzmann, cMu, cElectronMass
    use ModAdvance,     ONLY: State_VGB
    use ModMain,        ONLY: iTest, jTest, kTest, BlkTest
    use ModResistivity, ONLY: Eta0SI

    !------------------------------------------------------------------------
    ! The State_V vector is in normalized units
    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in) :: i, j, k, iBlock, iDir
    real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
    real, optional, intent(in)  :: TeIn                    ! [K]
    real, optional, intent(out) :: NatomicOut              ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut     ! dimensionless
    real, optional, intent(out) :: EinternalOut            ! [J/m^3]
    real, optional, intent(out) :: TeOut                   ! [K]
    real, optional, intent(out) :: PressureOut             ! [Pa]   
    real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]  
    real, optional, intent(out) :: GammaOut
    real, optional, intent(out) :: HeatCondOut             ! [Jm^2/(Ks)]   
    real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]  
    real, optional, intent(out) :: OpacityPlanckOut_W(nWave)      ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave)   ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3] 
    real, optional, intent(out) :: EntropyOut
    
    real, save :: KappaCoeffSI = (cBoltzmann/cElectronCharge)**2/cMu
    real :: nElec, EtaSI, TeSI!, HeatCond
    real, dimension(nIonFluid) :: nIon_I, fei_I, eiSigma_I
    real, dimension(nNeutral) :: fen_I, enSigma_I
    integer :: iIonFluid, iNeutral
    logical :: DoTest, DoTestMe=.true.

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest) then
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
 

    if(present(HeatCondOut)) then
       !!write(*,*)'iBlock = ',iBlock,'iNeutralBlockLast = ',iNeutralBlockLast
       if (iBlock.ne.iNeutralBlockLast) then
          call user_neutral_atmosphere(iBlock)
       end if

       call calc_electron_collision_rates(TeSI,i,j,k,iBlock,fen_I(1:nNeutral),fei_I(1:nIonFluid))
       eiSigma_I(1:nIonFluid) = cElectronCharge**2*nElec/((fei_I(1:nIonFluid)+1E-20)*cElectronMass) 
       enSigma_I(1:nNeutral) = cElectronCharge**2*nElec/((fen_I(1:nNeutral)+1E-20)*cElectronMass)
    
       !! EtaSI is calculated from both conductivities using Kirchhoff's rule:
       !! 1/sigma_tot = 1/eiSigma_IG+1/enSigma_IG
       !! The resulting conductivity is close to Spitzer conductivity far from the comet and
       !! decreases due to abundant electron-neutral collisions close to the nucleus
       
       !! EtaSI = 1/(sigma_tot*mu_0) magnetic diffusivity [m^2/s]
       EtaSI = Eta0SI + (sum(1/eiSigma_I(1:nIonFluid))+sum(1/enSigma_I(1:nNeutral)))/cMu
       HeatCondOut = TeSI/EtaSI*KappaCoeffSI
    end if

    if(DoTestMe) then
       write(*,*)'user_material_properties:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       if(present(CvOut)) write(*,*)'Cv     = ',CvOut,' [J/(K*m^3)]'
       if(present(HeatCondOut)) then
          do iIonFluid=1,nIonFluid
             write(*,*)'e & ',NameFluid_I(iIonFluid+1),':'
             write(*,*)'s_ei  = ',eiSigma_I(iIonFluid)," [1/(Ohm*m)]"
          end do
          do iNeutral=1,nNeutral
             write(*,*)'e & ',NameNeutral_I(iNeutral),':'
             write(*,*)'s_en  = ',enSigma_I(iNeutral)," [1/(Ohm*m)]"
          end do
          write(*,*)'Eta    = ',EtaSI," [m^2/s]"
          write(*,*)'Kappa  = ',HeatCondOut," [W/(m*K)]"
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
          PlotVar_G(i,j,k) = uIonMean_I(1)-Current_I(1)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_)/1E3
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
          PlotVar_G(i,j,k) = uIonMean_I(2)-Current_I(2)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_)/1E3
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
          PlotVar_G(i,j,k) = uIonMean_I(3)-Current_I(3)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_)/1E3
       end do; end do; end do

    case('dt')
       NameIdlUnit = 's'
       NameTecUnit = '[s]'
       PlotVar_G(:,:,:) = Dt_BLK(iBlock)*No2SI_V(UnitT_)
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
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var

  !========================================================================

  subroutine user_set_ICs(iBlock)
    use ModIO,       ONLY: restart
    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest, Body1_, Body1
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics
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

    !call user_neutral_atmosphere(iBlock)

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
          ! State_VGB(:,i,j,k,iBlock) = CellState_VI(:,1)

          State_VGB(SWpRho_,i,j,k,iBlock)    = SW_n*MassIon_I(SWp_)
          State_VGB(SWpRhoUx_,i,j,k,iBlock)  = SW_n*MassIon_I(SWp_)*SW_Ux
          State_VGB(SWpRhoUy_,i,j,k,iBlock)  = SW_n*MassIon_I(SWp_)*SW_Uy
          State_VGB(SWpRhoUz_,i,j,k,iBlock)  = SW_n*MassIon_I(SWp_)*SW_Uz
          State_VGB(SWpP_,i,j,k,iBlock)      = SW_n*SW_T_dim*Io2No_V(UnitTemperature_)

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
             if(ElectronPressureRatio.le.0.) call stop_mpi('ERROR: Electron Pressure Ratio > 0 for init!')
             State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))
             State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
          else
             State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))* &
                  (1.+ElectronPressureRatio)
          end if

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
    use ModPhysics,  ONLY: Si2No_V, UnitTemperature_

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

    use ModMain,       ONLY: Dt_BLK, BLKtest
    use ModPhysics,    ONLY: No2SI_V, UnitT_

    real, intent(out)             :: VarValue
    character (len=*), intent(in) :: TypeVar
    real, intent(in), optional    :: Radius

    character (len=*), parameter  :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------------
    select case(TypeVar)
    case('dtpnt')
       VarValue = Dt_BLK(BLKtest)*No2SI_V(UnitT_)
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
         BodyRho_I, BodyP_I, SW_Bz

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
       VarsGhostFace_V(SWpRho_)    = SW_n*MassIon_I(SWp_)*Io2No_V(UnitRho_)
       VarsGhostFace_V(SWpRhoUx_)  = SW_Ux
       VarsGhostFace_V(SWpRhoUy_)  = SW_Uy
       VarsGhostFace_V(SWpRhoUz_)  = SW_Uz
       VarsGhostFace_V(SWpP_)      = SW_n*Io2No_V(UnitN_)*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)
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

       ! ???
       !CellState_VI(:,iBoundary)=FaceState_VI(:,iBoundary)

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

end module ModUser
