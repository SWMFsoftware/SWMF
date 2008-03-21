!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the user module for Venus

  use ModSize
  use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_,p_,Bx_, By_, Bz_,&
       MassSpecies_V,SpeciesFirst_,SpeciesLast_  
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_boundary_cells,         &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_calc_sources,               &
       IMPLEMENTED7 => user_init_point_implicit,        &
       IMPLEMENTED10=> user_get_log_var,                &
       IMPLEMENTED11=> user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Venus 4 species MHD code, Yingjuan Ma'

  ! Radius within which the point implicit scheme should be used
  real :: rPointImplicit = 2.0

  ! Venus stuff
  logical ::  UseMultiSpecies=.true.
  integer, parameter :: MaxSpecies=4, MaxNuSpecies=3,  &
       MaxReactions=10
  integer :: nSpecies=4, nNuSpecies=3, &
       nReactions=10
  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK,MaxNuSpecies) :: &
       nDenNuSpecies_CBI    !number density of neutral Species

  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK) :: &
       Productrate_CB    !production rate according to optical depth

  real, dimension(MaxReactions) :: ReactionRate_I
  real, dimension(MaxReactions,MaxSpecies):: CoeffSpecies_II, &
       dSdRho_II !, dLdRho_II
  real, dimension(MaxSpecies)::LossSpecies_I, &
       SiSpecies_I,  LiSpecies_I,  PhoIon_I, Recb_I
  !        dStndRho_I,  dLtdRho_I,  dLtndNumRho_I, &
  real:: totalNumRho, totalLossRho,totalSourceRho, totalLossNumRho, &
       totalSourceNumRho, totalLossx, totalLossNumx

  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK) :: &
       MaxSiSpecies_CB,  MaxLiSpecies_CB
  common /TimeBlock/ MaxSiSpecies_CB,  MaxLiSpecies_CB

  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK) :: &
       MaxSLSpecies_CB

  !the reactions considered:(p means ion, em means electron)
  !the prefered order of a reaction is ions, Nus, hv and electrons
  integer, parameter :: &!reaction number
       CO2_hv__CO2p_em_=1 ,&!CO2+hv-->CO2p+em  
       O_hv__Op_em_=2     ,&!O+hv-->Op+em       
       CO2p_O__O2p_CO_=3  ,&!CO2p+O-->O2p+CO   
       Op_CO2__O2p_CO_=4  ,&!Op+CO2-->O2p+CO
       CO2p_O__Op_CO2_=5  ,&!CO2p+O-->Op+CO2
       O2p_em__O_O_=6     ,&!O2p+em-->O+O 
       CO2p_em__CO_O_=7   ,&!CO2p+em-->CO+O
       Hp_O__Op_H_=8      ,&!Hp+O-->Op+H
       Op_H__Hp_O_=9      ,&!Op+H-->Hp+O   
       H_hv__Hp_em_=10      !H+hv-->Hp+em

  real, dimension(MaxReactions) :: Rate_I
  real, dimension(MaxReactions) :: &  !for solar maximum condition, venus
       Ratedim_I=(/3.270e-6,1.224e-6, 1.64e-10, 1.1e-9, &
       9.60e-11, 7.38e-8, 3.1e-7, 5.084e-10, 6.4e-10, 5.58e-8 /)  !cm^3 s^(-1)

!CO2: 1.695e-6/0.72^2=3.26968  ; 6.346/0.72^2 = 12.24
  integer, parameter :: &! order of ion species
       Hp_  =1, &
       O2p_ =2, &
       Op_  =3, &
       CO2p_=4

  character (len=10), dimension(MaxSpecies):: &
       ion_name_I=(/'Hp  ', 'O2p ', 'Op  ','CO2p'/)

  real, dimension(MaxSpecies)::  &
       MassSpecies_I=(/1., 32., 16., 44. /)  !atm

  !  MassSpecies_I(Hp_)=1	 !atm
  !  MassSpecies_I(CO2p_)=44     !atm
  !  MassSpecies_I(O2p_)=32      !atm
  !  MassSpecies_I(Op_)=16       !atm

  integer, parameter :: & ! order of Neutral species
       CO2_=1 ,&
       O_  =2 ,&   
       H_  =3

  real, dimension(MaxNuSpecies)::CrossSection_I,&
       CrossSectionDim_I=(/2.6e-17,1.5e-17,0.0/)
  real:: Productrate0,Optdep
  real, dimension(MaxNuSpecies)::  NuMassSpecies_I=(/44,16,1/)
  !  NuMassSpecies_I(CO2_)=44	!atm
  !  NuMassSpecies_I(O_)=16	!atm

  real, dimension(MaxNuSpecies):: HNuSpecies_I,&
       HNuSpecies_dim_I=(/6.7e3,18.4e3, 1000.e3/)
  !HNuSpecies_dim_I(CO2_)=6.7e3   !m
  !HNuSpecies_dim_I(O_)=18.4e3    !m

  real, dimension(MaxNuSpecies):: BodynDenNuSpecies_I,&
       BodynDenNuSpDim_I=(/5.0e10,1.0e10, 0.0/)
  !BodynDenNuSpecies_dim_I(CO2_)=5e10 !cm^(-3)
  !BodynDenNuSpecies_dim_I(O_)=1e10  !cm^(-3)

  real, dimension(MaxSpecies):: BodyRhoSpecies_I
  integer, parameter :: & ! other numbers
       em_=-1 ,&
       hv_=-2   

  real :: TNu_body_dim = 1000.0, TNu_body, Tnu, Tnu_dim ! neutral temperature 
  real :: Ti_body_dim=1000.0, Ti_body  !ion temperature at the body
  real :: Tp_body_dim=2000.0, Tp_body  !dimensionless temperature 
  !of new created ions / plasma (Tp_body=2.0*Ti_body)

!  real :: XiT0, XiTx !dimensionless temperature of new created ions
!  real :: Ti_body_dim=300.0  !ion temperature at the body
!  real :: Ti_body_dim=1000.0  !ion temperature at the body
!  real :: Tnu_body_dim = 1000.0, Tnu_body, Tnu, Tnu_dim ! neutral temperature 
  real :: T300_dim = 300.0, T300 , Ti_dim =1000., Tp, Te0_dim=3000., Te0
  real,  dimension(1:nI,1:nJ,1:nK,nBLK) :: nu_BLK
  real :: nu0_dim=1.0e-9,nu0  !value from Tanaka, 1997, JGR

  !\
  ! The following are needed in user_sources::
  !/
  real, dimension(1:nI,1:nJ,1:nK):: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SE
  real, dimension(MaxSpecies,1:nI,1:nJ,1:nK) :: &
       SrhoSpecies

  character*30 :: SolarCondition, type_innerbcs='reflect'

contains
  !============================================================================

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case('#POINTIMPLICITREGION')
          call read_var('rPointImplicit',rPointImplicit)

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  !============================================================================
  subroutine user_init_point_implicit
    
    use ModVarIndexes
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet

    ! Allocate and set iVarPointImpl_I
    allocate(iVarPointImpl_I(8))

    iVarPointImpl_I = (/RhoHp_, RhoO2p_, RhoOp_, RhoCO2p_, &
         RhoUx_, RhoUy_, RhoUz_, P_/)

    ! Note that energy is not an independent variable for the 
    ! point implicit scheme. The pressure is an independent variable,
    ! and in this example there is no implicit pressure source term.

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.
    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

  !============================================================================
  
  subroutine user_calc_sources
    use ModAdvance,  ONLY: Source_VC,Energy_
    use ModNumConst, ONLY: cZero
    use ModVarIndexes, ONLY: rhoUx_, rhoUy_, rhoUz_
    use ModMain, ONLY: iTest, jTest, kTest, ProcTest, BlkTest, &
         GLOBALBLK
    use ModProcMH,   ONLY: iProc
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit, &
         IsPointImplSource
    use ModPhysics, ONLY: Rbody
    use ModGeometry,ONLY: R_BLK
    
    logical :: oktest,oktest_me
    !------------------------------------------------------------------------  
    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest('user_calc_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if
    if(UsePointImplicit)&
         UsePointImplicit_B(globalBLK) = &
         R_BLK(1,1,1,globalBLK) <= rPointImplicit &
         .and. R_BLK(nI,1,1,globalBLK) > rBody

    if(.not.(UsePointImplicit .and. UsePointImplicit_B(globalBLK)) )then
       ! Add all source terms if we do not use the point implicit 
       ! scheme for the Block
       call user_expl_source
       call user_impl_source

    elseif(IsPointImplSource)then
       ! Add implicit sources only
       call user_impl_source
    else
       ! Add explicit sources only
       call user_expl_source
    end if

    if(oktest_me)then
       write(*,*)'After Source(rho, rhoSp)=', &
            Source_VC(rho_:8,iTest,jTest,kTest)
       write(*,*)'Source(rhoU)=', Source_VC(5:8,iTest,jTest,kTest)
       write(*,*)'Source(B)=', Source_VC(9:11,iTest,jTest,kTest)
       write(*,*)'Source(p,E)', Source_VC(P_:P_+1,iTest,jTest,kTest)
    end if

  end subroutine user_calc_sources
  !=========================================================================

  subroutine user_impl_source
    use ModPointImplicit, ONLY: UsePointImplicit_B, &
         iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK
    use ModPhysics, ONLY: inv_gm1
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModGeometry,ONLY: r_BLK
    use ModNumConst,ONLY: cZero
    use ModVarIndexes,ONLY: Rho_, RhoHp_, RhoO2p_, RhoOp_, RhoCO2p_, &
         RhoUx_, RhoUy_, RhoUz_, P_, Energy_, Bx_, By_, Bz_
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest
    use ModProcMH,   ONLY: iProc
    !    use ModAdvance,  ONLY: Source_VC,Energy_
    !    use ModNumConst, ONLY: cZero
    logical :: oktest,oktest_me
    integer :: iBlock, i, j, k
    real    :: Coef
    !--------------------------------------------------------------------
    
    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest('user_imp_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if
    
    iBlock = GlobalBlk
     
    Srho   = cZero
    SrhoSpecies=cZero
    SrhoUx = cZero
    SrhoUy = cZero
    SrhoUz = cZero
    SBx    = cZero
    SBy    = cZero
    SBz    = cZero
    SP     = cZero
    SE     = cZero

    if(oktest_me)then
       !   write(*,*)'before Source(rhoU)=', Source_VC(6:8,itest,jtest,ktest)
       write(*,*)'Source(p,E)', Source_VC(P_:P_+1,iTest,jTest,kTest)
    end if

    call user_sources
    Source_VC(rho_       ,:,:,:) = Srho+Source_VC(rho_,:,:,:)
    Source_VC(rho_+1:rho_+MaxSpecies,:,:,:) = &
         SrhoSpecies+Source_VC(rho_+1:rho_+MaxSpecies,:,:,:)
    Source_VC(rhoUx_     ,:,:,:) = SrhoUx+Source_VC(rhoUx_,:,:,:)
    Source_VC(rhoUy_     ,:,:,:) = SrhoUy+Source_VC(rhoUy_,:,:,:)
    Source_VC(rhoUz_     ,:,:,:) = SrhoUz+Source_VC(rhoUz_,:,:,:)
    Source_VC(Bx_        ,:,:,:) = SBx+Source_VC(Bx_,:,:,:)
    Source_VC(By_        ,:,:,:) = SBy+Source_VC(By_,:,:,:)
    Source_VC(Bz_        ,:,:,:) = SBz+Source_VC(Bz_,:,:,:)
    Source_VC(P_     ,:,:,:) = SP+Source_VC(P_,:,:,:)
    Source_VC(Energy_,:,:,:) = SE+Source_VC(Energy_,:,:,:)

  end subroutine user_impl_source

  !===========================================================================

  subroutine user_expl_source
    !    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK
    !    use ModPointImplicit,ONLY: UsePointImplicit, UsePointImplicit_B

    !---------------------------------------------------------------------
    ! Here come the explicit source terms

  end subroutine user_expl_source

  !========================================================================

  subroutine user_sources
    use ModMain, ONLY: PROCTEST,GLOBALBLK,BLKTEST, iTest,jTest,kTest 
    use ModAdvance,  ONLY: State_VGB,VdtFace_x,VdtFace_y,VdtFace_z
    use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_,p_,Bx_, By_, Bz_
    use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK,R_BLK,&
         vInv_CB
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: Rbody, inv_gm1, gm1
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit

    ! Variables required by this user subroutine
    integer:: i,j,k,iSpecies, iBlock
    real :: inv_rho, inv_rho2, uu2, cosSZA, Productrate
    real :: alt, Te_dim = 300.0
    real :: totalPSNumRho=0.0,totalRLNumRhox=0.0, temps
    logical:: oktest,oktest_me
    real :: SourceLossMax, vdtmin
    !
    !\
    ! Variable meanings:
    !   Srho: Source terms for the continuity equation
    !   SE,SP: Source terms for the energy (conservative) and presure
    !          (primative) equations
    !   SrhoUx,SrhoUy,SrhoUz:  Source terms for the momentum equation
    !   SBx,SBy,SBz:  Souce terms for the magnetic field equations 
    !/
    !
    !--------------------------------------------------------------------------

    iBlock = globalBlk

    if (iProc==PROCtest.and.globalBLK==BLKtest) then
       call set_oktest('user_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    do k = 1, nK ;   do j = 1, nJ ;  do i = 1, nI
       inv_rho = 1.00/State_VGB(rho_,i,j,k,globalBLK)
       inv_rho2 = inv_rho*inv_rho
       uu2 =(State_VGB(Ux_,i,j,k,globalBLK)*State_VGB(Ux_,i,j,k,globalBLK)  &
            +State_VGB(Uy_,i,j,k,globalBLK)*State_VGB(Uy_,i,j,k,globalBLK)  &
            +State_VGB(Uz_,i,j,k,globalBLK)*State_VGB(Uz_,i,j,k,globalBLK)) &
            *inv_rho2
       
       SrhoUx(i,j,k) = SrhoUx(i,j,k) &
            -nu_BLK(i,j,k,globalBLK)*State_VGB(Ux_,i,j,k,globalBLK)
       SrhoUy(i,j,k) = SrhoUy(i,j,k)  &
            -nu_BLK(i,j,k,globalBLK)*State_VGB(Uy_,i,j,k,globalBLK)
       SrhoUz(i,j,k) = SrhoUz(i,j,k)  &
            -nu_BLK(i,j,k,globalBLK)*State_VGB(Uz_,i,j,k,globalBLK)
       
       if (R_BLK(i,j,k,globalBLK) > Rbody &
            .and.R_BLK(i,j,k,globalBLK) < 2.0 ) then
          
          ReactionRate_I=0.0
          CoeffSpecies_II(:,:)=0.0
          PhoIon_I(:)=0.0
          Recb_I(:)=0.0
          LossSpecies_I=0.0
          !totalNumRho=0.0
          SiSpecies_I(:)=0.0
          LiSpecies_I(:)=0.0
          
          totalLossRho=0.0
          totalLossNumRho=0.0
          totalSourceNumRho=0.0
          totalLossx=0.0
          totalLossNumx=0.0
          totalPSNumRho=0.0
          totalRLNumRhox=0.0
          

          totalNumRho=sum(State_VGB(rho_+1:rho_+nSpecies,i,j,k,globalBLK) &
               /MassSpecies_I(1:nSpecies))
          MaxSLSpecies_CB(i,j,k,globalBLK)=1.0e-3
         
          Productrate= Productrate_CB(i,j,k,globalBLK)
          ReactionRate_I(CO2_hv__CO2p_em_)= &
               Rate_I(CO2_hv__CO2p_em_)&
               *nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)
          PhoIon_I(CO2p_)=ReactionRate_I(CO2_hv__CO2p_em_) &
               *Productrate
          ReactionRate_I(O_hv__Op_em_)= &
               Rate_I(O_hv__Op_em_)&
               *nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
          PhoIon_I(Op_)=ReactionRate_I(O_hv__Op_em_) &
               *Productrate

!hydrogen density is neglected in the case
!!!              ReactionRate_I(H_hv__Hp_em_)= &
!!!                   Rate_I(H_hv__Hp_em_)&
!!!                   *nDenNuSpecies_CBI(i,j,k,globalBLK,H_)
!!!              PhoIon_I(Hp_)=ReactionRate_I(H_hv__Hp_em_) &
!!!                   *Productrate

!!!              Alt = (R_BLK(i,j,k,globalBLK)-1.0)*6052.0
!!!              if (Alt < 200.0 )then
!!!                 Te_dim = 300.0 + (Alt - 140.0)*3.7e3/60.0
!!!              else if( Alt < 800.0)then
!!!                 Te_dim = 4.0e3 + (Alt - 200.0)*5.0           
!!!              else
!!!                 Te_dim =7.0e3
!!!              end if
          Te_dim = 1000.0
          Ti_dim= Te_dim
          
          !charge exchange
          ReactionRate_I(CO2p_O__O2p_CO_)= &
               Rate_I(CO2p_O__O2p_CO_)&
               * nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
          CoeffSpecies_II(O2p_,CO2p_)=ReactionRate_I(CO2p_O__O2p_CO_)
          
          ReactionRate_I(Op_CO2__O2p_CO_)= &
               Rate_I(Op_CO2__O2p_CO_)&
               * nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)&
               *exp(log(Tnu_body_dim/Ti_dim)*0.39)
          CoeffSpecies_II(O2p_, Op_)=ReactionRate_I(Op_CO2__O2p_CO_)
          
          ReactionRate_I(CO2p_O__Op_CO2_)= &
               Rate_I(CO2p_O__Op_CO2_)&
               * nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
          CoeffSpecies_II(Op_,CO2p_)=ReactionRate_I(CO2p_O__Op_CO2_)
          
!!!              ReactionRate_I(Hp_O__Op_H_)= &
!!!                   Rate_I(Hp_O__Op_H_)* nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
!!!              CoeffSpecies_II(Op_,Hp_)=ReactionRate_I(Hp_O__Op_H_)
!!!
!!!              ReactionRate_I(Op_H__Hp_O_)= &
!!!                   Rate_I(Op_H__Hp_O_)* nDenNuSpecies_CBI(i,j,k,globalBLK,H_)
!!!              CoeffSpecies_II(Hp_,Op_)=ReactionRate_I(Op_H__Hp_O_)


          ! Recombination
          
          !              ReactionRate_I(O2p_em__O_O_)=Rate_I(O2p_em__O_O_)
          !              Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)
          
          !              ReactionRate_I(CO2p_em__CO_O_)=Rate_I(CO2p_em__CO_O_)
          !              Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)
          ! Recombination
          Tp=State_VGB(P_,i,j,k,globalBLK)/totalNumRho
          ReactionRate_I(O2p_em__O_O_)=Rate_I(O2p_em__O_O_)
          !          Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)*&
          !exp(log(TNu_body_dim/Te_dim)*0.56)
          Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)*&
               exp(log(TNu_body*2.0/Tp)*0.56)
          
          ReactionRate_I(CO2p_em__CO_O_)=Rate_I(CO2p_em__CO_O_)
          !          Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)*&
          !               sqrt(TNu_body_dim/Te_dim)
          Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)*&
               sqrt(TNu_body*2.0/Tp)
          
          
          !end if  !(x>0.0)
          
          do iSpecies=1, nSpecies
             LossSpecies_I=LossSpecies_I &
                  +CoeffSpecies_II(iSpecies, :)
             !                dStndRho_I=dStndRho_I  &
             !                     +CoeffSpecies_II(iSpecies, :)/MassSpecies_I(:)
             dSdRho_II(1:nSpecies, iSpecies)= &
                  CoeffSpecies_II(1:nSpecies, iSpecies)&
                  *MassSpecies_I(1:nSpecies)/MassSpecies_I(iSpecies)
             
          enddo
          
!!!              do iSpecies=1, nSpecies
!!!                 dLdRho_II(1:nSpecies, iSpecies)=Recb_I(1:nSpecies)&
!!!                      *rhoSpecies_GBI(i,j,k,globalBLK,1:nSpecies) &
!!!                      /MassSpecies_I(iSpecies)
!!!                 dLdRho_II(iSpecies, iSpecies)=  &
!!!                      dLdRho_II(iSpecies, iSpecies) &
!!!                      +LossSpecies_I(iSpecies)  &
!!!                      +Recb_I(iSpecies)*totalNumRho
!!!              enddo
!!!              !              dSLdRho_II=dSdRho_II-dLdRho_II
!!!
!!!              do iSpecies=1, nSpecies
!!!                 dLtdRho_I(:)=dLtdRho_I(:) +dLdRho_II(iSpecies,:)
!!!                 dLtndNumRho_I(:)=dLtndNumRho_I(:) &
!!!                      +dLdRho_II(iSpecies,:)*MassSpecies_I(:)/MassSpecies_I(iSpecies)
!!!              enddo              

          SiSpecies_I(:)=PhoIon_I(:)*MassSpecies_I(:)
          
          do iSpecies=1, nSpecies
             SiSpecies_I(1:nSpecies)=&
                  SiSpecies_I(1:nSpecies)  &
                  +dSdRho_II(1:nSpecies, iSpecies) &
                  *State_VGB(rho_+iSpecies, i,j,k, globalBLK)
             LiSpecies_I(iSpecies)= &
                  LiSpecies_I(iSpecies)  &
                  +(LossSpecies_I(iSpecies) +Recb_I(iSpecies)*totalNumRho)&
                  *State_VGB(rho_+iSpecies, i,j,k, globalBLK)
          enddo
          
          totalSourceRho=sum(SiSpecies_I(1:nSpecies))    
          totalLossRho=sum(LiSpecies_I(1:nSpecies))    
          !sum of the (Loss term) of all ion species
          totalLossNumRho=sum(LiSpecies_I(1:nSpecies)&
               /MassSpecies_I(1:nSpecies))   
          !sum of the (loss term/atom mass) of all ..
          totalSourceNumRho=sum(SiSpecies_I(1:nSpecies)&
               /MassSpecies_I(1:nSpecies))
          ! sum of the (Source term/atom mass) of all..
          totalLossx=totalLossRho*inv_rho
          totalLossNumx=totalLossNumRho/totalNumRho
          totalPSNumRho=sum(PhoIon_I(:)) 
          ! sum of the photonionziation source/atom mass) of all..
          totalRLNumRhox=sum(Recb_I(:) &
               *State_VGB(rho_+1:rho_+nSpecies, i,j,k, globalBLK)/MassSpecies_I(:))
          !sum of the (loss term/atom mass) due to recombination
          
          
          
          MaxSLSpecies_CB(i,j,k,globalBLK)=maxval(abs(SiSpecies_I(1:nSpecies)+&
               LiSpecies_I(1:nSpecies) ) /&
               (State_VGB(rho_+1:rho_+MaxSpecies, i,j,k, globalBLK)+1e-20))&
               /vInv_CB(i,j,k,globalBLK)
          
          if(.not.UsePointImplicit_B(iBlock) )then
             !sum of the (loss term/atom mass) due to recombination             
             SourceLossMax = 10.0*maxval(abs(SiSpecies_I(1:nSpecies)-&
                  LiSpecies_I(1:nSpecies) ) /&
                  (State_VGB(rho_+1:rho_+nSpecies, i,j,k, iBlock)+1e-20))&
                  /vInv_CB(i,j,k,iBlock)
             vdtmin=min(VdtFace_x(i,j,k),VdtFace_y(i,j,k),VdtFace_z(i,j,k))
             if(SourceLossMax > Vdtmin) then
                !UsePointImplicit_B(iBlock)=.true.                              
                !write(*,*)'should use Point-implicit or increase its region'   
                VdtFace_x(i,j,k) = max (SourceLossMax, VdtFace_x(i,j,k) )
                VdtFace_y(i,j,k) = max (SourceLossMax, VdtFace_y(i,j,k) )
                VdtFace_z(i,j,k) = max (SourceLossMax, VdtFace_z(i,j,k) )
             end if
          end if

          SrhoSpecies(1:nSpecies,i,j,k)=SrhoSpecies(1:nSpecies,i,j,k)&
               +SiSpecies_I(1:nSpecies) &
               -LiSpecies_I(1:nSpecies)
          
          Srho(i,j,k)=Srho(i,j,k)&
               +sum(SiSpecies_I(1:MaxSpecies))&
               -sum(LiSpecies_I(1:MaxSpecies))
          
          SrhoUx(i,j,k) = SrhoUx(i,j,k) &
               -State_VGB(Ux_,i,j,k,globalBLK)*totalLossx  
          
          SrhoUy(i,j,k) = SrhoUy(i,j,k)  &
               -State_VGB(Uy_,i,j,k,globalBLK)*totalLossx 
          
          SrhoUz(i,j,k) = SrhoUz(i,j,k)  &
               -State_VGB(Uz_,i,j,k,globalBLK)*totalLossx 
          
 !----- pressure and energy source terms
          temps = (totalSourceNumRho*Tp_body + totalPSNumRho*Te0) &
               -(totalLossNumx+totalRLNumRhox)*State_VGB(P_,i,j,k,globalBLK)
               

          SE(i,j,k) = SE(i,j,k) + (inv_gm1*temps-0.50*uu2*(totalLossRho)     ) 
          SP(i,j,k) = SP(i,j,k) + (temps        +0.50*uu2*(totalSourceRho)*gm1)

          SE(i,j,k) = SE(i,j,k)  &
               -0.5*State_VGB(rho_,i,j,k,globalBLK)*uu2*&
               nu_BLK(i,j,k,globalBLK) 

          SP(i,j,k) = SP(i,j,k)  &
               +gm1*0.5*State_VGB(rho_,i,j,k,globalBLK)*uu2*&
               nu_BLK(i,j,k,globalBLK) 

          if((State_VGB(P_,i,j,k,globalBLK)/totalNumRho).gt.(2.0*Tnu_body))then
             SE(i,j,k) = SE(i,j,k)  &
                  -nu_BLK(i,j,k,globalBLK)*totalNumRho*inv_gm1&
                  *(State_VGB(P_,i,j,k,globalBLK)/totalNumRho/2.0-Tnu_body)
             SP(i,j,k) = SP(i,j,k)  &
                  -nu_BLK(i,j,k,globalBLK)*totalNumRho &
                  *(State_VGB(P_,i,j,k,globalBLK)/totalNumRho/2.0-Tnu_body)
          end if
       else         



       endif !R_BLK(i,j,k,globalBLK) >= Rbody?
       
    end do; end do; end do     ! end of the i,j,k loop

  end subroutine user_sources

  !============================================================================
  subroutine user_init_session
    use ModMain
    use ModPhysics
    use ModVarIndexes

    integer::iBoundary
    !--------------------------------------------------------------------------
    !For Outer Boundaries
    do iBoundary=East_,Top_
       FaceState_VI(rhoHp_,iBoundary)    = SW_rho
       FaceState_VI(rhoO2p_,iBoundary)   = cTiny8
       FaceState_VI(rhoOp_,iBoundary)    = cTiny8
       FaceState_VI(rhoCO2p_,iBoundary)  = cTiny8     
    end do
    call set_multiSp_ICs  
    !    Rbody = 1.0 + 140.0e3/Rvenus
    BodyRho_I(1) = sum(BodyRhoSpecies_I(1:MaxSpecies))
    BodyP_I(1)   = SW_p*sum(BodyRhoSpecies_I(1:MaxSpecies)&
         /MassSpecies_I(1:MaxSpecies))*Ti_body_dim/SW_T_dim

    FaceState_VI(rho_,body1_)=BodyRho_I(1)
    FaceState_VI(rhoHp_:rhoCO2p_,body1_) = BodyRhoSpecies_I
    FaceState_VI(P_,body1_)=BodyP_I(1)
    CellState_VI(:,body1_:Top_)=FaceState_VI(:,body1_:Top_)
    do iBoundary=body1_,Top_  
       CellState_VI(rhoUx_:rhoUz_,iBoundary) = &
            FaceState_VI(Ux_:Uz_,iBoundary)*FaceState_VI(rho_,iBoundary)
    end do

    UnitUser_V(rhoHp_:rhoCO2p_)   = No2Io_V(UnitRho_)/MassSpecies_V

  end subroutine user_init_session

  !========================================================================

  subroutine user_set_ICs
    use ModProcMH, ONLY : iProc
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
    use ModIO, ONLY : restart
    use ModPhysics

    real :: Rmax, SinSlope, CosSlope,CosSZA
    real :: B4, dB4dx, zeta4, q4, epsi4, plobe, &
         XFace, YFace, ZFace
    integer :: i,j,k
    integer:: iBoundary
    logical::okTestMe=.false., okTest=.false.
    !-------------------------------------------------------------------------
    if(globalBLK==BLKtest .and. iProc==PROCtest)then
       call set_oktest('user_set_ics',oktest,oktestme)
    else
       oktest=.false.; oktestme=.false.
    endif

    if(okTestMe)then
       write(*,*)'BodynDenNuSpecies_I(:)=',&
            BodynDenNuSpecies_I(:)
       WRITE(*,*)''
       write(*,*)'HNuSpecies_I(1:nNuSpecies)=',HNuSpecies_I(:)
       WRITE(*,*)''
       write(*,*)'Rbody=', Rbody
       write(*,*)''
    end if

    !calculate neutral density
    do k=1,nK; do j=1,nJ; do i=1,nI
       if(R_BLK(i,j,k,globalBLK)<= Rbody)then
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=&
               BodynDenNuSpecies_I(:)
       else if(R_BLK(i,j,k,globalBLK)< 2.0) then
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=&
               BodynDenNuSpecies_I(:)* & 
               exp(-(R_BLK(i,j,k,globalBLK)-Rbody)&
               /HNuSpecies_I(:))
       else
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=0.0
       end if
    end do; end do; end do
    call neutral_density_averages
    do k=1,nK; do j=1,nJ; do i=1,nI
       nu_BLK(i,j,k,globalBLK)=&
            sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:nNuSPecies))*nu0
    end do; end do; end do

    ! calculate optical depth and producation rate
    do k=1,nK; do j=1,nJ; do i=1,nI
       cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
            x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)&
            +5.0e-4

       Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:MaxNuSpecies)*&
            CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
            6.0e-3)/cosSZA
       if( Optdep<11.5 .and. x_BLK(i,j,k,globalBLK) > 0.0) then 
          Productrate_CB(i,j,k,globalBLK) = max(exp(-Optdep), 1.0e-5)
       else
          Productrate_CB(i,j,k,globalBLK) = 1.0e-5
       end if

    end do; end do; end do

    if(.not.restart)then

       do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          if (R_BLK(i,j,k,globalBLK)< Rbody) then
             cosSZA=(0.5+sign(0.5,x_BLK(i,j,k,globalBLK)))*&
                  x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
                  1.0e-3
             State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
             !           State_VGB(rhoOp_,i,j,k,globalBLK)= 0.0
             !           State_VGB(rhoO2p_,i,j,k,globalBLK)= 0.0
             !           State_VGB(rhoCO2p_,i,j,k,globalBLK)= 0.0

             State_VGB(rhoOp_,i,j,k,globalBLK)= &
                  CellState_VI(rhoOp_,body1_)*cosSZA
             State_VGB(rhoO2p_,i,j,k,globalBLK)= &
                  CellState_VI(rhoOp_,body1_)*sqrt(cosSZA)
             State_VGB(rhoCO2p_,i,j,k,globalBLK)= &
                  CellState_VI(rhoOp_,body1_)*cosSZA
             State_VGB(rho_,i,j,k,globalBLK)  = &
                  sum( State_VGB(rho_+1:rho_+MaxSpecies,i,j,k,globalBLK))
             State_VGB(P_,i,j,k,globalBLK) = max(SW_p, &
                  sum(State_VGB(rho_+1:rho_+MaxSpecies,i,j,k,globalBLK) &
                  /MassSpecies_I(1:MaxSpecies))*Tp_body )

          else
             State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)
          end if
       end do;end do; end do;

       if(1==2.and.globalBLK==43)&
            write(*,*)'state_VGB(body1_)=',&
            CellState_VI(:,body1_),'(1)=',CellState_VI(:,1)

       do k=1,nK; do j=1,nJ; do i=1,nI
          
          if (true_cell(i,j,k,globalBLK).and. &
               R_BLK(i,j,k,globalBLK)<1.2*Rbody) then

             cosSZA=(0.5+sign(0.5,x_BLK(i,j,k,globalBLK)))*&
                  x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
                  1.0e-3

             State_VGB(rhoCO2p_,i,j,k,globalBLK)= Rate_I(CO2_hv__CO2p_em_)*&
                  cosSZA &
                  *nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)/&
                  nDenNuSpecies_CBI(i,j,k,globalBLK,O_)/&
                  (Rate_I(CO2p_O__O2p_CO_)+Rate_I(CO2p_O__Op_CO2_))
             State_VGB(rhoOp_,i,j,k,globalBLK)= (Rate_I(O_hv__Op_em_)*&
                  cosSZA+&
                  Rate_I(CO2p_O__Op_CO2_)*&
                  State_VGB(rhoCO2p_,i,j,k,globalBLK))&
                  *nDenNuSpecies_CBI(i,j,k,globalBLK,O_)&
                  /(nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)+3.0e5)&
                  /Rate_I(Op_CO2__O2p_CO_)
             State_VGB(rhoO2p_,i,j,k,globalBLK)= &
                  SQRT((nDenNuSpecies_CBI(i,j,k,globalBLK,O_)*&
                  State_VGB(rhoCO2p_,i,j,k,globalBLK)*&
                  Rate_I(CO2p_O__O2p_CO_)+&
                  nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)*&
                  State_VGB(rhoOp_,i,j,k,globalBLK)*&
                  Rate_I(Op_CO2__O2p_CO_)+1e-10)/Rate_I(O2p_em__O_O_))

             State_VGB(rhoO2p_:rhoCO2p_,i,j,k,globalBLK)=&
                  State_VGB(rhoO2p_:rhoCO2p_,i,j,k,globalBLK)*&
                  MassSpecies_I(O2p_:CO2p_)

          end if !(true_cell?)

       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI

          if(.not.true_cell(i,j,k,globalBLK))CYCLE 
          State_VGB(rho_,i,j,k,globalBLK)   =&
               sum(State_VGB(rho_+1:rho_+MaxSpecies,i,j,k,globalBLK))
          State_VGB(P_,i,j,k,globalBLK)= &
               max(SW_p, sum(State_VGB(rho_+1:rho_+MaxSpecies,i,j,k,globalBLK)&
               /MassSpecies_I(1:MaxSpecies))*Tp_body)
       end do; end do; end do

    end if
    time_BLK(:,:,:,globalBLK) = 0.00

    if(okTestMe)then
       write(*,*)'initial set up'
       write(*,*)'Rate_I=',Rate_I
       write(*,*)''
       write(*,*)'rhoSpecies_GBI(testcell,1:nSpecies)=',&
            State_VGB(rho_+1:rho_+MaxSpecies,itest,jtest,ktest,BLKtest)
       write(*,*)'p_BLK(itest,jtest,ktest,BLKtest)=',&
            State_VGB(P_,itest,jtest,ktest,BLKtest) 
       write(*,*)''
       write(*,*)'rhoSpecies_GBI(testcell+1,1:nSpecies)=',&
            State_VGB(rho_+1:rho_+MaxSpecies,itest+1,jtest,ktest,BLKtest)
       write(*,*)'p_BLK(itest+1,jtest,ktest,BLKtest)=',&
            State_VGB(P_,itest+1,jtest,ktest,BLKtest) 
       write(*,*)''
    end if


  end subroutine user_set_ICs

  !========================================================================

  subroutine set_multiSp_ICs

    ! Calculate the scale height of ion and neutal species and 
    ! intial boundary value of ion species

    use ModMain
    use ModConst
    use ModIO
    use ModPhysics

    real :: Productrate
    logical::oktest=.true., oktestme=.false.
    !---------------------------------------------------------------
    if(oktestme)then
       write(*,*)'in set_multisp_ICs, No2Io_V(UnitN_),t=',&
            No2Io_V(UnitN_),No2Io_V(UnitT_)
       write(*,*)'No2Si_V(UnitX_), temperature=',&
            No2Si_V(UnitX_), No2Si_V(UnitTemperature_)
       write(*,*)'kTp=',SW_p*Tp_body_dim/SW_T_dim, &
            Tp_body_dim/No2Si_V(UnitTemperature_)
       write(*,*)'BodynDenNuSpecies_dim_I(:)',&
            BodynDenNuSpdim_I(:)
    end if

    nu0=nu0_dim*No2Io_V(UnitN_)*No2Io_V(UnitT_)
    Tp_body = Tp_body_dim*Si2No_V(UnitTemperature_)

    Tnu_body = Tnu_body_dim *Si2No_V(UnitTemperature_)
    T300 = T300_dim*Si2No_V(UnitTemperature_)
    Te0 = Te0_dim*Si2No_V(UnitTemperature_)
    BodynDenNuSpecies_I(1:nNuSpecies)=&
         BodynDenNuSpDim_I(1:nNuSpecies)/No2Io_V(UnitN_)
    HNuSpecies_I(1:nNuSpecies)=&
         HNuSpecies_dim_I(1:nNuSpecies)*Si2No_V(UnitX_)

    ! normlize the reaction rate
    Rate_I(CO2_hv__CO2p_em_)= &
         Ratedim_I(CO2_hv__CO2p_em_)*No2Io_V(UnitT_)
    Rate_I(O_hv__Op_em_)=  &
         Ratedim_I(O_hv__Op_em_)*No2Io_V(UnitT_)
    Rate_I(CO2p_O__O2p_CO_)=  &
         Ratedim_I(CO2p_O__O2p_CO_)  &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(Op_CO2__O2p_CO_)=  &
         Ratedim_I(Op_CO2__O2p_CO_)*exp(log(8.0/3.0*T300/Tnu_body)*0.39) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_O__Op_CO2_)=  &
         Ratedim_I(CO2p_O__Op_CO2_) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)

    Rate_I(O2p_em__O_O_)=  &
         Ratedim_I(O2p_em__O_O_)*exp(log(4.0*T300/Tnu_body)*0.56)&
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_em__CO_O_)=  &
         Ratedim_I(CO2p_em__CO_O_)*sqrt(T300/Tnu_body)&
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)

    Rate_I(H_hv__Hp_em_)=  &
         Ratedim_I(H_hv__Hp_em_)*No2Io_V(UnitT_)
    Rate_I(Hp_O__Op_H_)=  &
         Ratedim_I(Hp_O__Op_H_) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)    
    Rate_I(Op_H__Hp_O_)=  &
         Ratedim_I(Op_H__Hp_O_) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)


    ReactionRate_I(CO2_hv__CO2p_em_)= &
         Rate_I(CO2_hv__CO2p_em_)*BodynDenNuSpecies_I(CO2_)
    PhoIon_I(CO2p_)=ReactionRate_I(CO2_hv__CO2p_em_) 

    ReactionRate_I(O_hv__Op_em_)= &
         Rate_I(O_hv__Op_em_)*BodynDenNuSpecies_I(O_)
    PhoIon_I(Op_)=ReactionRate_I(O_hv__Op_em_) 

    !charge exchange
    ReactionRate_I(CO2p_O__O2p_CO_)= &
         Rate_I(CO2p_O__O2p_CO_)* BodynDenNuSpecies_I(O_)
    CoeffSpecies_II(O2p_,CO2p_)=ReactionRate_I(CO2p_O__O2p_CO_)

    ReactionRate_I(Op_CO2__O2p_CO_)= &
         Rate_I(Op_CO2__O2p_CO_)* BodynDenNuSpecies_I(CO2_)
    CoeffSpecies_II(O2p_, Op_)=ReactionRate_I(Op_CO2__O2p_CO_)

    ReactionRate_I(CO2p_O__Op_CO2_)= &
         Rate_I(CO2p_O__Op_CO2_)* BodynDenNuSpecies_I(O_)
    CoeffSpecies_II(Op_,CO2p_)=ReactionRate_I(CO2p_O__Op_CO2_)

    !ion density at the body
    CrossSection_I=CrossSectionDim_I*No2Io_V(unitN_)*No2Si_V(unitX_)*1.0e2
    Productrate =1.0
    BodyRhoSpecies_I(Hp_)=SW_rho*0.1

    BodyRhoSpecies_I(CO2p_)= Rate_I(CO2_hv__CO2p_em_)*Productrate*&
         BodynDenNuSpecies_I(CO2_)/BodynDenNuSpecies_I(O_)/&
         (Rate_I(CO2p_O__O2p_CO_)+Rate_I(CO2p_O__Op_CO2_))
    BodyRhoSpecies_I(Op_)= (Rate_I(O_hv__Op_em_)*Productrate+&
         Rate_I(CO2p_O__Op_CO2_)*BodyRhoSpecies_I(CO2p_))&
         *BodynDenNuSpecies_I(O_)/(BodynDenNuSpecies_I(CO2_)+3.0e5)/&
         Rate_I(Op_CO2__O2p_CO_)
    BodyRhoSpecies_I(O2p_)= SQRT((BodynDenNuSpecies_I(O_)*&
         BodyRhoSpecies_I(CO2p_)*Rate_I(CO2p_O__O2p_CO_)+ &
         BodynDenNuSpecies_I(CO2_)*BodyRhoSpecies_I(Op_)*&
         Rate_I(Op_CO2__O2p_CO_))/Rate_I(O2p_em__O_O_))
    BodyRhoSpecies_I(:)=BodyRhoSpecies_I(:)*&
         MassSpecies_I(:)

    if(oktest)then
       write(*,*)' set parameters of Mars: BodyRhoSpecies_I(i)=',&
            BodyRhoSpecies_I(1:nSpecies)
       write(*,*)'neutral density=', &
            BodynDenNuSpecies_I(:)
       write(*,*)'nu0=',nu0
       write(*,*)'Rate_I=', Rate_I
       write(*,*)'Rate_dim_I=', Ratedim_I       
!       call stop_mpi('end')  
    end if

  end subroutine set_multiSp_ICs

  !========================================================================
  subroutine user_set_boundary_cells(iBLK)
    use ModGeometry
    use ModMain	

    integer,intent(in)::iBLK
    !-----------------------------------------------------------------------
    !  SHOULD define IsBoundaryCell_GI(:,:,:,ExtraBc_) using
    !  a boundary condition for iBLK block
    !  EXAMPLE: OUTER SPHERICAL BOUNDARY of radius of 100.
    !  IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<100.
    if (index(TypeGeometry,'spherical')>0)then
       if(XyzStart_BLK(Theta_,iBLK)<dz_BLK(iBLK))then
          !	IsBoundaryCell_GI(:,:,1-gcn:0,ExtraBc_)=.true.
          !	IsBoundaryCell_GI(1:nI,1:nJ,1-gcn:0,ExtraBc_)=.false.

          !	IsBoundaryCell_GI(:,:,1-gcn:0,ExtraBc_)=.true.
          IsBoundaryCell_GI(nI+1:nI+gcn,:,1-gcn:0,ExtraBc_)=.true.
          IsBoundaryCell_GI(1-gcn:0,:,1-gcn:0,ExtraBc_)=.true.	
       elseif(XyzStart_BLK(Theta_,iBLK)+nK*dz_BLK(iBLK)>cPi)then
          !        IsBoundaryCell_GI(:,:,nK+1:nK+gcn,ExtraBc_)=.true.
          !        IsBoundaryCell_GI(1:nI,1:nJ,nK+1:nK+gcn,ExtraBc_)=.false.

          !        IsBoundaryCell_GI(:,:,nK+1:nK+gcn,ExtraBc_)=.true.
          IsBoundaryCell_GI(nI+1:nI+gcn,:,nK+1:nK+gcn,ExtraBc_)=.true.
          IsBoundaryCell_GI(1-gcn:0,:,nK+1:nK+gcn,ExtraBc_)=.true.
       end if
    end if
  end subroutine user_set_boundary_cells

  !========================================================================

  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: West_, North_, Top_
    use ModVarIndexes, ONLY: nVar, RhoOp_, RhoO2p_, RhoCO2p_, RhoHp_
    use ModPhysics,    ONLY: SW_rho, SW_p, SW_T_dim
    use ModFaceBc,     ONLY: FaceCoords_D, VarsTrueFace_V

    real, intent(out):: VarsGhostFace_V(nVar)

    real:: XFace,YFace,ZFace, rFace, rFace2
    ! real:: v_phi(3) 
    real:: cosSZA 
    real:: uDotR, bDotR
    !--------------------------------------------------------------------------

    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)

    rFace2 = XFace**2 + YFace**2 + ZFace**2
    rFace  = sqrt(rFace2)

    !Apply boundary conditions
    cosSZA = (0.5+sign(0.5,XFace)) * XFace/max(rFace,1.0e-3) + 1.0e-3

    VarsGhostFace_V(rhoOp_) =  BodyRhoSpecies_I(Op_) *cosSZA

    VarsGhostFace_V(rhoO2p_) = BodyRhoSpecies_I(O2p_)*sqrt(cosSZA)

    VarsGhostFace_V(rhoCO2p_)=BodyRhoSpecies_I(CO2p_)*cosSZA

    VarsGhostFace_V(rhoHp_)=SW_rho*0.3

    VarsGhostFace_V(rho_) = sum(VarsGhostFace_V(rho_+1:rho_+MaxSpecies))
    VarsGhostFace_V(P_)=sum(VarsGhostFace_V(rho_+1:rho_+MaxSpecies)&
         /MassSpecies_I)*SW_p*Ti_body_dim/SW_T_dim 

    ! Reflective in radial direction
    uDotR = sum(VarsTrueFace_V(Ux_:Uz_)*FaceCoords_D)/rFace2
    bDotR = sum(VarsTrueFace_V(Bx_:Bz_)*FaceCoords_D)/rFace2

    select case (type_innerbcs) 
    case('float')
       VarsGhostFace_V(Ux_:Uz_) = &
            VarsTrueFace_V(Ux_:Uz_)
       VarsGhostFace_V(Bx_:Bz_) = &
            VarsTrueFace_V(Bx_:Bz_)
    case('reflect')   
       VarsGhostFace_V(Ux_:Uz_) = &
            VarsTrueFace_V(Ux_:Uz_) - 2*uDotR*FaceCoords_D
       VarsGhostFace_V(Bx_:Bz_) = &
            VarsTrueFace_V(Bx_:Bz_) - 2*bDotR*FaceCoords_D
    case('zeroB')
       VarsGhostFace_V(Ux_:Uz_) = &
            VarsTrueFace_V(Ux_:Uz_) - 2*uDotR*FaceCoords_D
       VarsGhostFace_V(Bx_:Bz_) = 0.0       
    case('default')
       write(*,*)'unknown type of user inner bcs'
    end select

    ! Apply corotation?
    !if (UseRotatingBcHere) then
    !   call calc_corotation_velocities(FaceCoords_D, v_phi)
    !   VarsGhostFace_V(Ux_:Uz_) = VarsGhostFace_V(Ux_:Uz_) + 2*v_phi
    !end if

  end subroutine user_face_bcs

  !====================================================================

  subroutine user_set_plot_var1(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: RhoHp_, RhoCO2p_, RhoO2p_, RhoOp_ 
    use ModPhysics,    ONLY: No2Io_V, UnitN_, NameTecUnit_V, NameIdlUnit_V
    use ModAdvance,    ONLY: State_VGB, Rho_

    ! Returns dimensional number density (/cc)

    integer,          intent(in) :: iBlock
    character(len=*), intent(in) :: NameVar
    logical,          intent(in) :: IsDimensional
    real,             intent(out):: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out):: PlotVarBody
    logical,          intent(out):: UsePlotVarBody
    character(len=*), intent(out):: NameTecVar
    character(len=*), intent(out):: NameTecUnit
    character(len=*), intent(out):: NameIdlUnit
    logical,          intent(out):: IsFound

    integer :: iVar
    real :: Coeff
    !--------------------------------------------------------------------------

    IsFound = .true.
    select case(NameVar)
    case('nhp')
       iVar = Hp_
       NameTecVar = 'H^+'
    case('nco2p')
       iVar = CO2p_
       NameTecVar = 'CO_2^+'
    case('no2p')
       iVar = O2p_
       NameTecVar = 'O_2^+'
    case('nop')
       iVar = Op_
       NameTecVar = 'O^+'
    case default
       IsFound = .false.  
    end select
    NameTecUnit = NameTecUnit_V(UnitN_)
    NameIdlUnit = NameIdlUnit_V(UnitN_)

    Coeff          = No2Io_V(UnitN_)/MassSpecies_I(iVar)
    PlotVar_G      = Coeff*State_VGB(Rho_+iVar,:,:,:,iBlock)
    PlotVarBody    = BodyRhoSpecies_I(iVar)
    UsePlotVarBody = .True.

  end subroutine user_set_plot_var1

  !====================================================================
  subroutine neutral_density_averages
    use ModMain
    use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK, true_cell,vInv_CB, R_BLK
    use ModNumConst
    use ModPhysics, ONLY : cTolerance
    use ModCovariant, ONLY : FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB
    integer :: i,j,k,iNu
    real:: FaceArea_DS(3,east_:top_),VInv

    real ::  density_IS(6,nNuSpecies),x,y,z,R0, factor
    ! real :: neutral_density
    ! Note: using true_cell to replace an Rbody test does not apply here
    !----------------------------------------------------------------

    do k=1,nK; do j=1,nJ;  do i=1,nI  
       VInv=vInv_CB(i,j,k,globalBLK)

       if(.not.true_cell(i,j,k,globalBLK))cycle
       !-------------------East----------------------------------
       x = 0.5*(x_BLK(i-1,j,k,globalBLK) + x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i-1,j,k,globalBLK) + y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i-1,j,k,globalBLK) + z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,East_)= FaceAreaI_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,East_)*x+ &
            FaceArea_DS(2,East_)*y+ &
            FaceArea_DS(3,East_)*z)/R0
       do iNu = 1, nNuSpecies 
          density_IS(East_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------West----------------------------------
       x = 0.5*(x_BLK(i+1,j,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i+1,j,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i+1,j,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,West_)= FaceAreaI_DFB(:,i+1,j,k,globalBLK)
       factor = (FaceArea_DS(1,West_)*x+ &
            FaceArea_DS(2,West_)*y+ &
            FaceArea_DS(3,West_)*z)/R0     
       do iNu = 1, nNuSpecies 
          density_IS(West_,iNu) =-neutral_density(R0,iNu)*factor
       end do

       !-------------------South----------------------------------
       x = 0.5*(x_BLK(i,j-1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i,j-1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i,j-1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,South_)=FaceAreaJ_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,South_)*x+ &
            FaceArea_DS(2,South_)*y+ &
            FaceArea_DS(3,South_)*z)/R0  
       do iNu = 1, nNuSpecies 
          density_IS(South_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------North----------------------------------
       x = 0.5*(x_BLK(i,j+1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i,j+1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i,j+1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,North_)=FaceAreaJ_DFB(:,i,j+1,k,globalBLK)
       factor = (FaceArea_DS(1,North_)*x+ &
            FaceArea_DS(2,North_)*y+ &
            FaceArea_DS(3,North_)*z)/R0     
       do iNu = 1, nNuSpecies 
          density_IS(North_,iNu) = -neutral_density(R0,iNu)*factor
       end do

       !-------------------Bot----------------------------------
       x = 0.5*(x_BLK(i,j,k-1,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i,j,k-1,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i,j,k-1,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,Bot_)= FaceAreaK_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,Bot_)*x+ &
            FaceArea_DS(2,Bot_)*y+ &
            FaceArea_DS(3,Bot_)*z)/R0
       do iNu = 1, nNuSpecies 
          density_IS(Bot_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------Top----------------------------------
       x = 0.5*(x_BLK(i,j,k+1,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = 0.5*(y_BLK(i,j,k+1,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = 0.5*(z_BLK(i,j,k+1,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,Top_)= FaceAreaK_DFB(:,i,j,k+1,globalBLK)
       factor = (FaceArea_DS(1,Top_)*x+ &
            FaceArea_DS(2,Top_)*y+ &
            FaceArea_DS(3,Top_)*z)/R0 
       do iNu = 1, nNuSpecies 
          density_IS(Top_,iNu) = -neutral_density(R0,iNu)*factor
       end do

       !-------------------SUM----------------------------------
       do iNu = 1, nNuSpecies 
          nDenNuSpecies_CBI(i,j,k,globalBLK,iNu)=VInv* &
               sum(density_IS(:,iNu))&
               *HNuSpecies_I(iNu)*BodynDenNuSpecies_I(iNu)
          if(nDenNuSpecies_CBI(i,j,k,globalBLK,iNu)<0)then
             write(*,*)'wrong sign, i,j,k,golablBLK, iNu',&
                  i,j,k,globalBLK,iNu, R_BLK(i,j,k,globalBLK)
          end if
       end do

    end do; end do ;end do 

  end subroutine neutral_density_averages

  !==============================================================================
  real function neutral_density(R0,iNu)
    !  use ModUser, ONLY : BodynDenNuSpecies_I, HNuSpecies_I
    use ModPhysics, ONLY : Rbody

    real, intent(in) :: R0
    integer, intent(in) :: iNu

    !-----------------------------------------------------------------------
    neutral_density = 0.0
    if( R0 >= 0.9*Rbody .and. R0< 3.0*Rbody ) &
         neutral_density= exp(-(R0-Rbody)/HNuSpecies_I(iNu))

  end function neutral_density
  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,&
         dx_BLK,dy_BLK,dz_BLK
    use ModMain,       ONLY: unusedBLK
    use ModVarIndexes, ONLY: Rho_, rhoHp_, rhoO2p_, RhoOp_,RhoCO2p_,&
         rhoUx_,rhoUy_,rhoUz_
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK
    use ModPhysics,ONLY: No2Si_V, UnitN_, UnitX_, UnitU_

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    real, external :: calc_sphere
    real ::mass
    integer:: i,j,k,iBLK
    character (len=*), parameter :: Name='user_get_log_var'
    logical:: oktest=.false.,oktest_me
    !-------------------------------------------------------------------
    call set_oktest('user_get_log_var',oktest,oktest_me)
    if(oktest)write(*,*)'in user_get_log_var: TypeVar=',TypeVar
    select case(TypeVar)
    case('hpflx')
       mass=1.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(rhoHp_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
          end do; end do; end do          
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)

    case('opflx')
       mass=16.
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(rhoOp_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
          end do; end do; end do          
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)    
    case('o2pflx')
       mass=32.
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(rhoO2p_,i,j,k,iBLK)*&
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
          end do; end do; end do          
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
    
    case('co2pflx')
       mass=44.
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(rhoCO2p_,i,j,k,iBLK)*&
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
          end do; end do; end do          
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
    
    case default
       call stop_mpi('wrong logvarname')
    end select
    !change to user value from normalized flux
    !    write(*,*)'varvalue, unitSI_n, unitSI_x, unitSI_U, mass, unitSI_t=',&
    !         varvalue, unitSI_n, unitSI_x, unitSI_U, mass, unitSI_t
    VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass

  end subroutine user_get_log_var

!=====================================================================
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)

    use ModPhysics, ONLY: Rbody

    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'

    !-------------------------------------------------------------------
    !    select case (InitialRefineType)
    !    case ('Venus3Dbodyfocus')
    ! Refine, focusing on body
    found=.true.
    refineBlock=.false.
    if (maxR > Rbody.and.(lev <= 1 .or. minR < 1.5*Rbody))&
         refineBlock = .true.

    !    case default

    !   end select

  end subroutine user_specify_initial_refinement

end module ModUser
