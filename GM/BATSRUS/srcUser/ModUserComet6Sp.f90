!^CFG COPYRIGHT UM
!this file contains the ModUser for general Cometary MHD
!MassSpecies is defined at [SpeciesFirst:SpeciesLast]
! Last 01.15 committed in
! This 01.16 change Te calculation, add ini rates for solar max.
! This ModUser.f90 requires copying ModEquationsMHDComet.f90 into ModEquation.f90.
! jpattern = 0, 10, 11, 12, sets the electron impact ionization term, and
! the way electron temperature is treated for dissociative recombination processes[Gombosi96].
!==============================================================================

module ModUser
  ! This is the user module for Comets
  use ModSize
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_init_point_implicit,        &
       IMPLEMENTED5 => user_calc_sources,               &
       IMPLEMENTED6 => user_update_states

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real, parameter :: VersionUserModule = 01.13
  character (len=*), parameter :: NameUserModule = &
       '6 Species Cometary MHD Module, Yingdong Jia 2008'

  real ::  Qprod=7.e-9, ionization_rate=1.e-6, Unr, Unr_km=1.0, kin, kin_cc=1.7e-9
  real ::  Tion, mbar, lambda
  integer :: jpattern
  logical :: DoTest, DoTestMe

  integer ::  nSpecies=1
  integer, parameter :: MaxSpecies=6, MaxNuSpecies=6,  &
       MaxIni=11, MaxCx=13
  integer ::  nNuSpecies=6	!nNuSpecies<6 doesn't work for this version yet

  logical :: DoInitialize = .True.
  real, public, dimension(0:nI+1, 0:nJ+1, 0:nK+1, nBLK, MaxNuSpecies) :: &
       NNeu_BLK = 1.
  real, public, dimension(0:nI+1, 0:nJ+1, 0:nK+1, nBLK, 3           ) :: &
       UNeu_BLK = 0.

  integer, parameter :: H2Op_ =1, Hp_ =2, H3Op_ =3, OHp_ =4, Op_ =5, COp_ =6

  integer, parameter :: H2O_=1, H_=2, OH_=3, O_=4, CO_=5, H2_=6

  real, dimension(MaxNuSpecies)::  NuMassSp_I=(/18.,1.,17.,16.,28.,2./)
  !at least for single neutral fluid version this neutral mass is useless.

  integer, parameter :: &	!ionization reaction number
       H2O_H2Op_ = 1 ,&
       H2O_Hp_   = 2 ,&
       H_Hp_     = 3 ,&
       OH_Hp_    = 4 ,&
       H2O_OHp_  = 5 ,&
       OH_OHp_   = 6 ,&
       H2O_Op_   = 7 ,&
       OH_Op_    = 8 ,&
       O_Op_     = 9 ,&
       CO_Op_    = 10,&
       CO_COp_   = 11

  !The ionization rates that put down here are midified in user_init by
  ! a factor (ionization_rate) that reflecting the solar distance.
  real, dimension(MaxIni) :: IoniRate_I=(/5.4e-7, .2e-7, &
       1.2e-7, .6e-7, .9e-7, 3.94e-7, 9.5e-9, 5.3e-8, 3.45e-7, &
       .4e-7, 6.25e-7/)  !s^(-1), solar minimum conditions at 1AU

  !  real, dimension(MaxIni) :: IoniRate_I=(/5.4e-7, .2e-7, &
  !	1.2e-7, .6e-7, .9e-7, 3.94e-7, 9.5e-9, 5.3e-8, 3.45e-7, &
  !	.4e-7, 6.25e-7/)  !s^(-1), solar maximum conditions at 1AU

  integer, parameter :: &	!CX reaction number
       H2Op_H2O__H3Op_OH_ = 1 ,&
       H2Op_CO__COp_H2O_  = 2 ,&
       H2Op_H2__H3Op_H_   = 3 ,&
       Hp_H2O__H2Op_H_    = 4 ,&
       Hp_OH__OHp_H_      = 5 ,&
       Hp_O__Op_H_        = 6 ,&
       OHp_H2O__H3Op_O_   = 7 ,&
       OHp_H2O__H2Op_OH_  = 8 ,&
       OHp_CO__COp_OH_    = 9 ,&
       Op_H2O__H2Op_O_    = 10,&
       Op_OH__OHp_O_      = 11,&
       COp_H2O__H2Op_CO_  = 12,&
       COp_OH__OHp_CO_    = 13

  real, dimension(MaxCx) :: CxRate_I = (/2.e-15, 4.3e-16, 0.8e-15, 8.2e-15, &
       4.4e-15, 0.38e-15, 1.3e-15, 1.6e-15, 0.8e-15, 3.2e-15, 3.2e-15, &
       1.6e-15, 0.3e-15/)  !m^3 s^(-1), already in SI unit.

contains

  !========================================================================
  !  SUBROUTINE user_read_inputs
  !========================================================================
  subroutine user_read_inputs

    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i, j, k, ireadvar
    character (len=100) :: NameCommand, line
    !-------------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input COMET starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#COMETPARAM")
          call read_var('Qprod' ,Qprod)
          call read_var('Unr_km' ,Unr_km)
          call read_var('mbar',mbar)
          call read_var('ionization_rate',ionization_rate)
          call read_var('kin_cc',kin_cc)
          call read_var('Tion', Tion)
          call read_var('jpattern' ,jpattern)
          DoInitialize = .True.

       case("#MultiSP")		!ini rates hardwired in parameters.
          !         call read_var('UseMultiSpecies',UseMultiSpecies)
          !nNuSpecies<6 doesn't work for this version yet
          !         if(UseMultiSpecies)then
          call read_var('nSpecies'  , nSpecies  )
          call read_var('nNuSpecies', nNuSpecies)
          !	 endif

       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input COMET ends'
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

    !    call stop_user(Name)
  end subroutine user_read_inputs

  !========================================================================
  subroutine user_init_session		!called before reading the Param.in
    !========================================================================
    use ModMain
    use ModPhysics
    use ModVarIndexes

    integer::iBoundary
    !--------------------------------------------------------------------------
    !For Outer Boundaries
    do iBoundary=East_,Top_
       FaceState_VI(SpeciesFirst_:SpeciesLast_,iBoundary) = cTiny8
       FaceState_VI(rhoHp_,iBoundary) = SW_rho
    end do
    FaceState_VI(SpeciesFirst_:SpeciesLast_,body1_) = BodyRho_I(1)/100.
    FaceState_VI(rhoH2Op_,body1_) = BodyRho_I(H2Op_)
    FaceState_VI(rhoHp_  ,body1_) = cTiny

    CellState_VI(SpeciesFirst_:SpeciesLast_,:) = &
         FaceState_VI(SpeciesFirst_:SpeciesLast_,:)

    kin=kin_cc*1E-6
    Unr=Unr_km*1E3
    lambda = Unr/ionization_rate	!in meters 
    ! should be lambda water(5.4e-7), but this is the same with 1sp model.

    !    if(UseMultiSpecies) mbar = NuMassSp_I(H2O_)

    !convert to real solar distance value
    IoniRate_I = IoniRate_I*ionization_rate*1e6

  end subroutine user_init_session

  !========================================================================
  subroutine user_set_ICs
    !========================================================================

    use ModMain, ONLY: globalBLK
    use ModPhysics, ONLY: CellState_VI
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes, ONLY: nVar
    use ModIO,ONLY: restart
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModNumConst
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitU_

    integer :: i,j,k, iVar
    !--------------------------------------------------------------------------

    if(restart)return	!jet and Bx parameters has to be set in 1st session.

    if(.not.restart)then
       do iVar=1,nVar
          State_VGB(iVar,:,:,:,globalBLK)   = CellState_VI(iVar,1)
       end do
    end if

    !define neutral speed array, necessary for special neutral distribution types only
    UNeu_BLK(:,:,:,globalBLK,1) = Unr* &
         x_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)/ &
         R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)
    UNeu_BLK(:,:,:,globalBLK,2) = Unr* &
         y_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)/ &
         R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)
    UNeu_BLK(:,:,:,globalBLK,3) = Unr* &
         z_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)/ &
         R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)
    UNeu_BLK(:,:,:,globalBLK,:) = &
         UNeu_BLK(:,:,:,globalBLK,:)/No2Si_V(UnitU_)

    !define neutral density array, necessary for special neutral distribution types only
    NNeu_BLK(:,:,:,globalBLK,H2O_) = Qprod * &	!!!Number Density!!!
         exp(-R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)*No2Si_V(UnitX_)/lambda)/ &
         (4*cPi*Unr*R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)**2*No2Si_V(UnitX_)**2)
    !minor species using estimation only [combi96,Haberli96], despite of chemistry/scale hight.
    NNeu_BLK(:,:,:,globalBLK,OH_ ) = NNeu_BLK(:,:,:,globalBLK,H2O_)/30.
    NNeu_BLK(:,:,:,globalBLK,O_  ) = NNeu_BLK(:,:,:,globalBLK,H2O_)/100.
    NNeu_BLK(:,:,:,globalBLK,CO_ ) = NNeu_BLK(:,:,:,globalBLK,H2O_)/6.
    NNeu_BLK(:,:,:,globalBLK,H2_ ) = NNeu_BLK(:,:,:,globalBLK,H2O_)/500.
    if ( jpattern == 22 ) then		!uniform estimated H density
       NNeu_BLK(:,:,:,globalBLK,H_) = Qprod * &
            exp(-R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)* &
            No2Si_V(UnitX_)/lambda/12.)/ &
            (4*cPi*Unr*R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)**2* &
            No2Si_V(UnitX_)**2 * 6)
    else		!Haser model for H density, should be used for other species later, esp. O and OH.
       NNeu_BLK(:,:,:,globalBLK,H_) = Qprod * 1.e7/(1.e5-1.e7)* &
            ( exp(-R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)* &
            No2Si_V(UnitX_)/lambda*10.) - &
            exp(-R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)* &
            No2Si_V(UnitX_)/lambda/12.) )/ &
            (4*cPi*Unr*R_BLK(0:nI+1,0:nJ+1,0:nK+1,globalBLK)**2* &
            No2Si_V(UnitX_)**2 * 6)
    endif
    if(nNuSpecies < 2) NNeu_BLK(:,:,:,globalBLK,2:MaxNuSpecies) = cTiny

  end subroutine user_set_ICs

  !============================================================================

  subroutine user_init_point_implicit

    use ModVarIndexes
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet

    ! Allocate and set iVarPointImpl_I
    allocate(iVarPointImpl_I(4+nSpecies))

    iVarPointImpl_I = (/RhoUx_, RhoUy_, RhoUz_, p_, &
         rhoH2Op_, rhoHp_, rhoH3Op_, rhoOHp_, rhoOp_, rhoCOp_/)

    ! Note that energy is not an independent variable for the 
    ! point implicit scheme. The pressure is an independent variable,
    ! and in this example there is no implicit pressure source term.

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.
    ! Initialization for comet implicit sources: false => numerical ptimplicit.
    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

  !========================================================================
  !  SUBROUTINE user_calc_sources
  !========================================================================
  subroutine user_calc_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModPointImplicit, ONLY: UsePointImplicit, IsPointImplSource
    use ModVarIndexes, ONLY: UseMultiSpecies

    if(.not.UsePointImplicit)then
       ! Add all source terms if we do not use the point implicit scheme
       call user_expl_source
       if(.not. UseMultiSpecies)then
	  call user_impl_source_sgl
       else
	  call user_impl_source_multi
       endif
    elseif(IsPointImplSource)then
       ! Add implicit sources only
       if(.not. UseMultiSpecies)then
	  call user_impl_source_sgl
       else
	  call user_impl_source_multi
       endif
    else
       ! Add explicit sources only
       call user_expl_source
    end if

  end subroutine user_calc_sources
  !==========================================================================
  subroutine user_expl_source

    ! Here come the explicit source terms
    ! The energy source is only needed in the explicit source term.

  end subroutine user_expl_source
  !==========================================================================
  subroutine user_impl_source_sgl
  end subroutine user_impl_source_sgl

  !========================================================================
  !  SUBROUTINE user_impl_source_multi
  !========================================================================
  subroutine user_impl_source_multi

    ! This is a test and example for using point implicit source terms
    ! Note that the energy is a dependent variable in the
    ! point implicit scheme, so there is no energy source here.
    ! The pressure source is zero.

    use ModPointImplicit, ONLY: &
         UsePointImplicit, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC

    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, n_step, BlkTest, &
         iTest, jTest, kTest, PROCTest
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModVarIndexes
    use ModPhysics
    use ModProcMH

    integer :: i, j, k, m
    real :: eta, usqr, Unsqr, totalNumRho, Te, alphaTe, fi, &
         totalSourceRho, totalLossRho, totalNumLoss, Rkm, logR
    real :: Source_H3Op, Loss_H3Op, NumLoss_H3Op, collisn	!H3Op special rxn
    real :: Un(3), U(3)
    real, dimension(MaxSpecies) :: Source_I, Loss_I, AlphaN_I, ni_I
    real, dimension(MaxNuSpecies) :: nn_I
    !-------------------------------------------------------------------------

    do k=1,nK ;   do j=1,nJ ;    do i=1,nI

       State_VGB(rho_,i,j,k,globalBLK) = &	!essential for pt-implicit
            sum( State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,globalBLK) )

       nn_I(1:MaxSpecies) = NNeu_BLK(i,j,k,globalBLK,1:MaxSpecies)*No2Si_V(UnitT_)
       ni_I(1:MaxSpecies) = State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,globalBLK) &
            /MassSpecies_V(SpeciesFirst_:SpeciesLast_)
       ! use of nSpecies is an example of how to change this into
       ! nSpecies version instead of 6 species version
       totalNumRho = sum(ni_I(1:nSpecies))

       !The following is for the special impact ionization and recombination Te [gombosi96].
       Rkm  = R_BLK(i,j,k,globalBLK)*No2Si_V(UnitX_)/1E3	!km unit
       logR = log10(Rkm)
       if (jpattern == 10 .or. jpattern == 11 ) then    ! jpattern= 0, regular comet[Gombosi97]
          if (rkm >= 5000. .and. rkm < 10000.) then     ! jpattern=10, Comet Halley[Gombosi96]
	     fi = 1.0+0.77*log(rkm/5000.)               ! jpattern=11, regular comet with electron impact ionization
          elseif (rkm >= 10000. .and. rkm < 50000.) then
             fi = 1.5-0.31067*log(rkm/10000.)           ! jpattern=12, regluar case with electron temperature
          else                                          !      profile measureed by Giotto
	     fi = cOne
          endif
       else
	  fi = cOne
       endif	! jpattern
       if (jpattern == 10 .or. jpattern == 12 ) then
	  if ( Rkm <= 1584.893 ) then
	     Te = 100.
          else if ( Rkm <= 6918.310 ) then
	     Te = 10.0**( 1.143*logR-1.667 )
          else if ( Rkm <= 1.e4 ) then
	     Te = 10.0**( 10.965*logR-39.3725 )
          else if ( Rkm <= 1.e5 ) then
	     Te = 10.0**( .5135*logR+2.4325 )
          else
	     Te = 1.e5
          endif
       else
	  Te = State_VGB(p_,i,j,k,globalBLK)*No2Si_V(UnitP_) / &
               ( 2*No2Si_V(UnitN_)*cBoltzmann*totalNumRho )
       endif	!jpattern

       If(Te < 200.) then
	  alphaTe  = 7.E-7*sqrt(300./Te)
       else
	  alphaTe  = 2.342*7.E-7*Te**(0.2553-0.1633*log10(Te))
       endif

       !define recombination terms ion by ion
       AlphaN_I(H2Op_) = alphaTe
       AlphaN_I(Hp_  ) = cZero
       AlphaN_I(H3Op_) = alphaTe
       AlphaN_I(Op_  ) = cZero
       AlphaN_I(OHp_ ) = 3.8E-8*sqrt(300./Te)
       AlphaN_I(COp_ ) = 1.E-7*(300./Te)**0.46
       ! normalize
       AlphaN_I = AlphaN_I*totalNumRho/1E6*No2Si_V(UnitN_)*No2Si_V(UnitT_)	

!!!!calculate source terms case by case---!!!
       !calculate ionization source terms
       Source_I(H2Op_) = IoniRate_I(H2O_H2Op_)*fi*nn_I(H2O_)
       Source_I(Hp_  ) = IoniRate_I(H2O_Hp_  )*nn_I(H2O_) + &
            IoniRate_I(H_Hp_  )*nn_I(H_  ) + &
            IoniRate_I(OH_Hp_ )*nn_I(OH_ )
       Source_I(H3Op_) = cZero
       Source_I(OHp_ ) = IoniRate_I(H2O_OHp_ )*nn_I(H2O_) + &
            IoniRate_I(OH_OHp_)*nn_I(OH_ )
       Source_I(Op_  ) = IoniRate_I(H2O_Op_  )*nn_I(H2O_) + &
            IoniRate_I(OH_Op_ )*nn_I(OH_ ) + &
            IoniRate_I(O_Op_  )*nn_I(O_  ) + &
            IoniRate_I(CO_Op_ )*nn_I(CO_ )
       Source_I(COp_ ) = IoniRate_I(CO_COp_  )*nn_I(CO_)
       Source_I = Source_I / No2Si_V(UnitN_)	!normalize

       !add CX source terms
       Source_I(H2Op_) = Source_I(H2Op_) + &
            ( ni_I(Hp_  )*CxRate_I(Hp_H2O__H2Op_H_   ) + &
            ni_I(OHp_ )*CxRate_I(OHp_H2O__H2Op_OH_ ) + &
            ni_I(Op_  )*CxRate_I(Op_H2O__H2Op_O_   ) + &
            ni_I(COp_ )*CxRate_I(COp_H2O__H2Op_CO_ ) )*nn_I(H2O_)
       Source_I(Hp_  ) = Source_I(Hp_  ) + &	!Temperarily using H2O mirror for H mirror
            ( ni_I(Hp_  )*CxRate_I(H2Op_H2O__H3Op_OH_) )*nn_I(H_  )
       Source_I(H3Op_) = Source_I(H3Op_) + &
            ( ni_I(H2Op_)*CxRate_I(H2Op_H2O__H3Op_OH_) + &
            ni_I(OHp_ )*CxRate_I(OHp_H2O__H3Op_O_  ) )*nn_I(H2O_) + &
            ni_I(H2Op_  )*CxRate_I(H2Op_H2__H3Op_H_  )  *nn_I(H2_ )
       Source_I(OHp_ ) = Source_I(OHp_ ) + &
            ( ni_I(Hp_  )*CxRate_I(Hp_OH__OHp_H_     ) + &
            ni_I(Op_  )*CxRate_I(Op_OH__OHp_O_     ) + &
            ni_I(COp_ )*CxRate_I(COp_OH__OHp_CO_   ) )*nn_I(OH_ )
       Source_I(Op_  ) = Source_I(Op_  ) + &
            ( ni_I(Hp_  )*CxRate_I(Hp_O__Op_H_       ) )*nn_I(O_  )
       Source_I(COp_ ) = Source_I(COp_ ) + &
            ( ni_I(H2Op_)*CxRate_I(H2Op_CO__COp_H2O_ ) + &
            ni_I(OHp_ )*CxRate_I(OHp_CO__COp_OH_   ) )*nn_I(CO_ )

       !define special charge exchange loss term
       Source_H3Op = MassSpecies_V(RhoH3Op_) * &
            ( ni_I(H2Op_)*CxRate_I(H2Op_H2O__H3Op_OH_)  *nn_I(H2O_) + &
            ni_I(H2Op_)*CxRate_I(H2Op_H2__H3Op_H_  )  *nn_I(H2_ ) )


       Source_I = Source_I * MassSpecies_V	!*m_s

!!!!calculate loss terms case by case---!!!!
       !calculate ionization source terms
       Loss_I(H2Op_) = CxRate_I(H2Op_H2O__H3Op_OH_)*nn_I(H2O_) + &
            CxRate_I(H2Op_CO__COp_H2O_)*nn_I(CO_ ) + &
            CxRate_I(H2Op_H2__H3Op_H_ )*nn_I(H2_ )
       Loss_I(Hp_  ) = CxRate_I(Hp_H2O__H2Op_H_   )*nn_I(H2O_) + &
            CxRate_I(H2Op_H2O__H3Op_OH_)*nn_I(H_ ) + &	!using H2O mirror for H mirror
            CxRate_I(Hp_OH__OHp_H_    )*nn_I(OH_ ) + &
            CxRate_I(Hp_O__Op_H_      )*nn_I(O_  )
       Loss_I(H3Op_) = cZero
       Loss_I(OHp_ ) = CxRate_I(OHp_H2O__H3Op_O_  )*nn_I(H2O_) + &
            CxRate_I(OHp_H2O__H2Op_OH_)*nn_I(H2O_) + &
            CxRate_I(OHp_CO__COp_OH_  )*nn_I(CO_ )
       Loss_I(Op_  ) = CxRate_I(Op_H2O__H2Op_O_   )*nn_I(H2O_) + &
            CxRate_I(Op_OH__OHp_O_    )*nn_I(OH_ )
       Loss_I(COp_ ) = CxRate_I(COp_H2O__H2Op_CO_ )*nn_I(H2O_) + &
            CxRate_I(COp_OH__OHp_CO_  )*nn_I(OH_ )

       Loss_I = ( Loss_I + AlphaN_I ) * &		!*ms*ns
            State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k, globalBLK)

       !define special charge exchange loss term
       Loss_H3Op = ( CxRate_I(H2Op_H2O__H3Op_OH_)*nn_I(H2O_) + &
            CxRate_I(H2Op_H2__H3Op_H_ )*nn_I(H2_ ) ) * &
            State_VGB(rhoH2Op_,i,j,k, globalBLK)
       NumLoss_H3Op = ( CxRate_I(H2Op_H2O__H3Op_OH_)*nn_I(H2O_) + &
            CxRate_I(H2Op_H2__H3Op_H_ )*nn_I(H2_ ) ) * ni_I(H2Op_)

       !calculate Source_VC
       totalSourceRho  = sum( Source_I )
       totalLossRho    = sum( Loss_I   )
       totalNumLoss    = sum( Loss_I / MassSpecies_V )
       collisn         = State_VGB(rho_,i,j,k, globalBLK)*kin*nn_I(H2O_)

       Un   = UNeu_BLK(i,j,k,globalBLK,1:3)
       U    = State_VGB(rhoUx_:rhoUz_,i,j,k,globalBLK)/ &
            State_VGB(rho_,i,j,k,globalBLK)
       unsqr= dot_product( Un, Un )
       usqr = dot_product( U , U  )

       Source_VC(SpeciesFirst_:SpeciesLast_,i,j,k) = &
            Source_VC(SpeciesFirst_:SpeciesLast_,i,j,k) + &
            Source_I - Loss_I

       do m = 1,3
	  Source_VC(rhoU_+m,i,j,k) = Source_VC(rhoU_+m,i,j,k) + &
               (totalSourceRho-Source_H3Op+collisn)*Un(m) - &
               (totalLossRho-Loss_H3Op+collisn)*U(m)
       enddo	!added H3O corrections and elastic collisions [ThesisJia07]

       Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + &
            (1.0/3.0)*(totalSourceRho-Source_H3Op+collisn)* &
            dot_product( (Un-U) , (Un-U) ) - &
            cHalf  *State_VGB(p_,i,j,k,globalBLK)/totalNumRho* &
            ( (totalNumLoss-NumLoss_H3Op+collisn/mbar) + sum(AlphaN_I/MassSpecies_V* &
            State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k, globalBLK)) )
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + &
            cHalf*( (totalSourceRho-Source_H3Op+collisn)*Unsqr - &
            (totalLossRho-Loss_H3Op+collisn)*usqr - &
            inv_gm1*State_VGB(p_,i,j,k,globalBLK)/totalNumRho* &
            ((totalNumLoss-NumLoss_H3Op+collisn/mbar) + sum(AlphaN_I/MassSpecies_V* &
            State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k, globalBLK))) )

       ! source energy for explicit, p for implicit
    end do; enddo; enddo

  end subroutine user_impl_source_multi


  !========================================================================
  !  SUBROUTINE user_update_states(iStage,iBlock)
  !========================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain
    use ModPhysics, ONLY: cBoltzmann, No2Si_V, UnitN_, UnitP_
    use CON_planet,  ONLY: NamePlanet
    use ModEnergy
    use ModGeometry,ONLY: x_BLK,y_BLK,z_BLK
    use ModProcMH

    integer,intent(in):: iStage,iBlock
    integer:: i,j,k
    real :: Pthmin_VC

    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update check of temperature::
    !/

    if (NamePlanet == 'HALLEY') then
       ! now check to see if the temperature is less than some
       ! prescribed minimum. If it is set it to the minimum value
       do k=1,nK ;   do j=1,nJ ;    do i=1,nI
          if ( UseMultiSpecies ) then
             Pthmin_VC = sum( &
                  State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,iBlock)/ &
                  MassSpecies_V(SpeciesFirst_:SpeciesLast_) )* &
                  No2Si_V(UnitN_)*cBoltzmann*Tion/No2Si_V(UnitP_)
          else
             Pthmin_VC = State_VGB(rho_,i,j,k,iBlock)/ &
                  mbar*No2Si_V(UnitN_)*cBoltzmann*Tion/No2Si_V(UnitP_)
          endif
          if( State_VGB(p_,i,j,k,iBlock) < Pthmin_VC ) &
               State_VGB(p_,i,j,k,iBlock) = Pthmin_VC
       enddo; enddo; enddo
    end if

    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/

  end subroutine user_update_states

end module ModUser
