!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the user module for Mars 

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
       IMPLEMENTED9 => user_get_b0,                     &
       IMPLEMENTED10 => user_get_log_var,               &
       IMPLEMENTED11 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: NameUserModule = &
       'Mars 4 species MHD code, Yingjuan Ma'

  character (len=10) :: SolarCond='solarmax  '

  ! Radius within which the point implicit scheme should be used
  real :: rPointImplicit = 2.0

  ! Mars stuff
  logical ::  UseMultiSpecies=.true.
  integer, parameter :: MaxSpecies=4, MaxNuSpecies=8,  &
       MaxReactions=10
  integer :: nSpecies=4, nNuSpecies=3, &
       nReactions=10
  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK,MaxNuSpecies) :: &
       nDenNuSpecies_CBI    !number density of neutral Species
  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK) :: &
       TempNuSpecies_CBI    !tempature of neutral Species
  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK) :: &
       Productrate_CB    !production rate according to optical depth
  real,  dimension(1:nI, 1:nJ, 1:nK, nBLK,2) :: &
       Ionizationrate_CBI !Ionization rate !!according to TGCM

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
       O_hv__Op_em_    =2 ,&!O+hv-->Op+em       
       CO2p_O__O2p_CO_ =3 ,&!CO2p+O-->O2p+CO   
       Op_CO2__O2p_CO_ =4 ,&!Op+CO2-->O2p+CO
       CO2p_O__Op_CO2_ =5 ,&!CO2p+O-->Op+CO2
       O2p_em__O_O_    =6 ,&!O2p+em-->O+O 
       CO2p_em__CO_O_  =7 ,&!CO2p+em-->CO+O
       Hp_O__Op_H_     =8 ,&!Hp+O-->Op+H
       Op_H__Hp_O_     =9 ,&!Op+H-->Hp+O   
       H_hv__Hp_em_    =10  !H+hv-->Hp+em

  real, dimension(MaxReactions) :: Rate_I
  real, dimension(MaxReactions) :: &
       Ratedim_I=(/ 2.47e-7, 8.89e-8, 1.64e-10, 1.1e-9, &
     9.60e-11, 7.38e-8, 3.1e-7, 5.084e-10, 6.4e-10, 5.58e-8 /)  !cm^3 s^(-1)

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
       O_=2   ,&    
       H_=3, &
       Oh_=4   ,&
       Ohx_=5 , &
       Hx_=6, &
       Ox_=7 ,&
       CO2x_=8   

  real, dimension(MaxNuSpecies)::CrossSection_I,&
       CrossSectionDim_I=(/2.6e-17,1.5e-17,0.0,1.5e-17,&
       1.5e-17,0.0,1.5e-17,2.6e-17/)

  real:: Productrate0,Optdep
  real, dimension(MaxNuSpecies)::NuMassSpecies_I=(/44,16,1,16,16,1,16, 44/)
  !  NuMassSpecies_I(CO2_)=44	!atm
  !  NuMassSpecies_I(O_)=16	!atm

  real, dimension(MaxNuSpecies):: HNuSpecies_I=1.0,HNuSpeciesDim_I=1.0
  !HNuSpecies_dim_I(CO2_)=6.7e3   !m
  !HNuSpecies_dim_I(O_)=18.4e3    !m

  real, dimension(MaxNuSpecies):: BodynDenNuSpecies_I,&
       BodynDenNuSpDim_I=(/1.1593e12, 3.2278e9, 1.1307e7, 1.951e4, &
       1.5248e3, 9.4936e5, 5.2695e8, 2.2258e11/)

  real, dimension(MaxSpecies):: BodyRhoSpecies_I
  integer, parameter :: & ! other numbers
       em_=-1 ,&
       hv_=-2   

  real :: TNu_body_dim = 300.0, kTn, Tnu, Tnu_dim ! neutral temperature 
  real :: Ti_body_dim=300.0, kTi0  !ion temperature at the body
  real :: Tp_body_dim=600.0, kTp0  !dimensionless temperature 
  !of new created ions / plasma (Tp_body=2.0*Ti_body)

  real :: Te_new_dim=1000., KTe0

  real :: T300_dim = 300.0, T300 , Ti_dim =300.
  real,  dimension(1:nI,1:nJ,1:nK,nBLK) :: nu_BLK,nu1_BLK
  real :: nu0_dim=4.0e-10,nu0

  ! coefficient of Mars magnetic field
  real, dimension(0:61,0:61) :: cmars, dmars
  integer :: NNm
  real :: mars_eps=1e-4
  logical :: UseMarsB0 = .false.
  real :: rot = 1.0, thetilt = 0.0
  logical :: UseHotO = .false.
  logical :: UseTempCont=.false.
  logical :: UseSolarMax=.false.
  logical :: UseImpactIon=.false.
  real, dimension(32,MaxSpecies)::Impact_ION,Impact_ION_dim=0.0 
  real, dimension(32):: Temp_dim
  logical :: UseChargeEx=.true.

  integer,parameter::NLong=73, NLat=36, MaxAlt=21
  real :: Long_I(NLong), Lat_I(NLat), Alt_I(MaxAlt)
  real :: Temp(NLong, NLat, MaxAlt)
  real :: Den_CO2(NLong, NLat, MaxAlt)!,Den_CO2_dim(NLong, NLat, NAlt)
  real :: Den_O(NLong, NLat, MaxAlt)!,Den_O_dim(NLong, NLat, NAlt)
  real :: ICO2p(NLong, NLat, MaxAlt)!,ICO2p_dim(NLong, NLat, NAlt)
  real :: IOp(NLong, NLat, MaxAlt)!,IOp_dim(NLong, NLat, NAlt)
  logical ::UseMarsAtm=.false.
  integer :: NAlt=21


  !\
  ! The following are needed in user_sources::
  !/
  real, dimension(1:nI,1:nJ,1:nK):: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SE
  real, dimension(MaxSpecies,1:nI,1:nJ,1:nK) :: &
       SrhoSpecies
  logical:: UseOldEnergy=.true.

contains
  !=============================================================================

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam

    character (len=100) :: NameCommand
    integer:: i, j, k, n, m
    character (len=60):: TGCMFilename  
    character (len=100) :: line
    !------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case('#USEOLDENERGY')
          call read_var('UseOldEnergy',UseOldEnergy)
          if(.not.UseOldEnergy)then
             call read_var('Te_new_dim',Te_new_dim)
             !change temperature from ev to k
             Te_new_dim = Te_new_dim * 11610.0  
          end if

       case("#UseMarsB0")  !if or not include crustal magnetic field of Mars
          call read_var('UseMarsB0',UseMarsB0)
          if(UseMarsB0) then
             call read_var('NNm', NNm)
             call read_var('rot',rot)
             call read_var('thetilt', thetilt)
             rot= rot/180.0*3.141592653589793238462643383279
             thetilt= thetilt/180.0*3.141592653589793238462643383279
             cmars = 0.0
             dmars = 0.0
             open(15,file='marsmgsp.txt')
             do i=0,NNm 
                read(15,*)n,(cmars(n-1,m),m=0,n-1),(dmars(n-1,m),m=0,n-1)
             end do
             close(15)
          endif

       case("#SOLARCON") !solar cycle condition
          call read_var('SolarCon',SolarCond)

       case("#UseSolarMax") !solar cycle condition
          call read_var('UseSolarMax',UseSolarMax)
       
       case("#UseHotO")  !adding hot Oxygen or not
          call read_var('UseHotO',UseHotO)
       
       case("#UseTempCont") !add hoc term of the energy source
          call read_var('UseTempCont',UseTempCont)          

       case('#REACTIONS')
          call read_var('UseImpactIon',UseImpactIon)
          call read_var('UseChargeEx',UseChargeEx)
          open(15,file='read_in.dat')
          do i=1,32
             read(15,*)Temp_dim(i),Impact_ION_dim(i,Op_),Impact_ION_dim(i,Hp_)
          end do
          close(15)
          
       case("#UseMarsAtm")
          call read_var('UseMarsAtm',UseMarsAtm)
          if(UseMarsAtm)then
             call read_var('TGCMFilename',TGCMFilename)
             call read_var('NAlt', Nalt)
             open(15,file=TGCMFilename,status="old")
             read(15,*)line
             write(*,*)line, Nalt
             do k = 1, NAlt
                do j=1, NLat
                   do i=1, NLong
                      read(15,*)Long_I(i),Lat_I(j),Alt_I(k),Temp(i,j,k),Den_CO2(i,j,k),&
                           Den_O(i,j,k),ICO2p(i,j,k),IOp(i,j,k)
                   end do
                end do
             end do
             close(15)
             write(*,*)Long_I(Nlong),Lat_I(NLat),Alt_I(Nalt)
             write(*,*)Long_I(1),Lat_I(1),Alt_I(1)
             write(*,*)'Den_O(i,j,k),ICO2p(i,j,k),IOp(i,j,k)=',&
                  Den_O(Nlong,Nlat,Nalt),ICO2p(Nlong,Nlat,Nalt),&
                  IOp(Nlong,Nlat,Nalt)
             
          end if
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
  !========================================================================
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
  !  SUBROUTINE USER_SOURCES
  !========================================================================
  !\
  ! This subroutine is used to calculate sources for the MHD equations.  The
  ! routine is called for each block separately so that the user would typically
  ! need only to code the source term calculation for a single block (in other
  ! words inside the the k,j,i loop below).  As with all user subroutines, the
  ! variables declared in ModUser are available here.  Again, as with other
  ! user subroutines DO NOT MODIFY ANY GLOBAL VARIABLE DEFINED IN THE MODULES
  ! INCLUDED IN THIS SUBROUTINE UNLESS SPECIFIED!!
  !
  ! The user should load the global variables:
  !      Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SE,SP,SEw
  !
  ! Note that SE (energy) and SP (pressure) must both be loaded if the code is 
  ! going to use both the primitive and the conservative MHD equation advance  
  ! (see the USER MANUAL and the DESIGN document).  If using only primitive SP 
  ! must be loaded.  If using only conservative SE must be loaded.  The safe
  ! approach is to load both.
  !/
  subroutine user_sources
    use ModMain, ONLY: PROCTEST,GLOBALBLK,BLKTEST, iTest,jTest,kTest 
    use ModAdvance,  ONLY: State_VGB,VdtFace_x,VdtFace_y,VdtFace_z
    use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_,p_,Bx_, By_, Bz_
    use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK,R_BLK,&
         vInv_CB
    use ModConst,    ONLY: cZero,cHalf,cOne,cTwo,cTolerance
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: Rbody, inv_gm1, gm1
!    use ModBlockData,ONLY: use_block_data, put_block_data, get_block_data
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit

    ! Variables required by this user subroutine
    integer:: i,j,k,iSpecies, iBlock
    real :: inv_rho, inv_rho2, uu2,Productrate,kTi,kTe
    real :: alt, Te_dim = 300.0, temp
    real :: totalPSNumRho=0.0,totalRLNumRhox=0.0, temps
    logical:: oktest,oktest_me
    real :: SourceLossMax, vdtmin

    !
    !---------------------------------------------------------------------------
    !\
    ! Variable meanings:
    !   Srho: Source terms for the continuity equation
    !   SE,SP: Source terms for the energy (conservative) and presure
    !          (primative) equations
    !   SrhoUx,SrhoUy,SrhoUz:  Source terms for the momentum equation
    !   SBx,SBy,SBz:  Souce terms for the magnetic field equations 
    !/
    !---------------------------------------------------------------------------
    !
    iBlock = globalBlk

    if (iProc==PROCtest.and.iBlock==BLKtest) then
       call set_oktest('user_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    if (R_BLK(1,1,1,iBlock) > 3.0*Rbody) RETURN
    
    do k = 1, nK ;   do j = 1, nJ ;  do i = 1, nI

       if (R_BLK(i,j,k,iBlock) < Rbody) CYCLE

       inv_rho = 1.00/State_VGB(rho_,i,j,k,iBlock)
       inv_rho2 = inv_rho*inv_rho
       uu2 =(State_VGB(Ux_,i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)  &
            +State_VGB(Uy_,i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)  &
            +State_VGB(Uz_,i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)) &
            *inv_rho2

       ReactionRate_I=0.0
       CoeffSpecies_II(:,:)=0.0
       PhoIon_I(:)=0.0
       Recb_I(:)=0.0
       LossSpecies_I=0.0
       SiSpecies_I(:)=0.0
       LiSpecies_I(:)=0.0
       
       totalSourceRho=0.0
       totalLossRho=0.0
       totalLossNumRho=0.0
       totalSourceNumRho=0.0
       totalLossx=0.0
       totalLossNumx=0.0
       totalPSNumRho=0.0
       totalRLNumRhox=0.0
       
          
       totalNumRho=sum(State_VGB(rho_+1:rho_+nSpecies,i,j,k,iBlock) &
            /MassSpecies_I(1:nSpecies))
       MaxSLSpecies_CB(i,j,k,iBlock)=1.0e-3
       
       Productrate= Productrate_CB(i,j,k,iBlock)
          
       !             ReactionRate_I(CO2_hv__CO2p_em_)= &
       !                  Rate_I(CO2_hv__CO2p_em_)&
       !                  *nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)
       !             PhoIon_I(CO2p_)=ReactionRate_I(CO2_hv__CO2p_em_) &
       !                  *Productrate
       !             ReactionRate_I(O_hv__Op_em_)= &
       !                  Rate_I(O_hv__Op_em_)&
       !                  *nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       !             PhoIon_I(Op_)=ReactionRate_I(O_hv__Op_em_) &
       !                  *Productrate
       
       ReactionRate_I(H_hv__Hp_em_)= &
            Rate_I(H_hv__Hp_em_)&
            *nDenNuSpecies_CBI(i,j,k,iBlock,H_)
       PhoIon_I(Hp_)=ReactionRate_I(H_hv__Hp_em_) &
            *Productrate
       
       PhoIon_I(CO2p_)=Ionizationrate_CBI(i,j,k,iBlock,CO2_)
       PhoIon_I(Op_)=Ionizationrate_CBI(i,j,k,iBlock,O_)
       
       
!!!              Alt = (R_BLK(i,j,k,iBlock)-1.0)*6052.0
!!!              if (Alt < 200.0 )then
!!!                 Te_dim = 300.0 + (Alt - 140.0)*3.7e3/60.0
!!!              else if( Alt < 800.0)then
!!!                 Te_dim = 4.0e3 + (Alt - 200.0)*5.0           
!!!              else
!!!                 Te_dim =7.0e3
!!!              end if
!!!          Te_dim = 300.0
!!!          Ti_dim= Te_dim

       !charge exchange
       ReactionRate_I(CO2p_O__O2p_CO_)= &
            Rate_I(CO2p_O__O2p_CO_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       CoeffSpecies_II(O2p_,CO2p_)=ReactionRate_I(CO2p_O__O2p_CO_)
       
       ReactionRate_I(Op_CO2__O2p_CO_)= &
            Rate_I(Op_CO2__O2p_CO_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)&
            *exp(log(Tnu_body_dim/Ti_dim)*0.39)
       CoeffSpecies_II(O2p_, Op_)=ReactionRate_I(Op_CO2__O2p_CO_)
       
       ReactionRate_I(CO2p_O__Op_CO2_)= &
            Rate_I(CO2p_O__Op_CO2_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       CoeffSpecies_II(Op_,CO2p_)=ReactionRate_I(CO2p_O__Op_CO2_)
       
!!!              ReactionRate_I(Hp_O__Op_H_)= &
!!!                   Rate_I(Hp_O__Op_H_)* nDenNuSpecies_CBI(i,j,k,iBlock,O_)
!!!              CoeffSpecies_II(Op_,Hp_)=ReactionRate_I(Hp_O__Op_H_)
!!!
!!!              ReactionRate_I(Op_H__Hp_O_)= &
!!!                   Rate_I(Op_H__Hp_O_)* nDenNuSpecies_CBI(i,j,k,iBlock,H_)
!!!              CoeffSpecies_II(Hp_,Op_)=ReactionRate_I(Op_H__Hp_O_)
          
          
       ! Recombination
       
       !              ReactionRate_I(O2p_em__O_O_)=Rate_I(O2p_em__O_O_)
       !              Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)
          
       !              ReactionRate_I(CO2p_em__CO_O_)=Rate_I(CO2p_em__CO_O_)
       !              Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)
       ! Recombination
       
       kTi=State_VGB(P_,i,j,k,iBlock)/totalNumRho/2.0
       kTe=kTi

       ReactionRate_I(O2p_em__O_O_)=Rate_I(O2p_em__O_O_)
       !          Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)*exp(log(TNu_body_dim/Te_dim)*0.56)
       Recb_I(O2p_)=ReactionRate_I(O2p_em__O_O_)*exp(log(kTn/kTi)*0.56)
       
       ReactionRate_I(CO2p_em__CO_O_)=Rate_I(CO2p_em__CO_O_)
       !          Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)*&
       !               sqrt(TNu_body_dim/Te_dim)
       Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)*&
            sqrt(kTn/kTi)
       
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
!!!                      *rhoSpecies_GBI(i,j,k,iBlock,1:nSpecies) &
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
               *State_VGB(rho_+iSpecies, i,j,k, iBlock)
          LiSpecies_I(iSpecies)= &
               LiSpecies_I(iSpecies)  &
               +(LossSpecies_I(iSpecies) +Recb_I(iSpecies)*totalNumRho)&
               *State_VGB(rho_+iSpecies, i,j,k, iBlock)
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
            *State_VGB(rho_+1:rho_+nSpecies, i,j,k, iBlock)/MassSpecies_I(:))
       !sum of the (loss term/atom mass) due to recombination
       
       
       
       MaxSLSpecies_CB(i,j,k,iBlock)=maxval(abs(SiSpecies_I(1:nSpecies)+&
            LiSpecies_I(1:nSpecies) ) /&
            (State_VGB(rho_+1:rho_+MaxSpecies, i,j,k, iBlock)+1e-20))&
            /vInv_CB(i,j,k,iBlock)

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
       
       !          VdtFace_x(i,j,k) = max (MaxSLSpecies_CB(i,j,k,iBlock),&
       !               VdtFace_x(i,j,k) )
       !          VdtFace_y(i,j,k) = max (MaxSLSpecies_CB(i,j,k,iBlock),&
       !               VdtFace_y(i,j,k) )
       !          VdtFace_z(i,j,k) = max (MaxSLSpecies_CB(i,j,k,iBlock),&
       !               VdtFace_z(i,j,k) )

       SrhoSpecies(1:nSpecies,i,j,k)=SrhoSpecies(1:nSpecies,i,j,k)&
            +SiSpecies_I(1:nSpecies) &
            -LiSpecies_I(1:nSpecies)

       Srho(i,j,k)=Srho(i,j,k)&
            +sum(SiSpecies_I(1:MaxSpecies))&
            -sum(LiSpecies_I(1:MaxSpecies))

       SrhoUx(i,j,k) = SrhoUx(i,j,k) &
            -State_VGB(Ux_,i,j,k,iBlock)*totalLossx  

       SrhoUy(i,j,k) = SrhoUy(i,j,k)  &
            -State_VGB(Uy_,i,j,k,iBlock)*totalLossx 

       SrhoUz(i,j,k) = SrhoUz(i,j,k)  &
            -State_VGB(Uz_,i,j,k,iBlock)*totalLossx 

       SrhoUx(i,j,k) = SrhoUx(i,j,k) &
            -nu_BLK(i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)
       SrhoUy(i,j,k) = SrhoUy(i,j,k)  &
            -nu_BLK(i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)
       SrhoUz(i,j,k) = SrhoUz(i,j,k)  &
            -nu_BLK(i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)

          
!----- pressure and energy source terms
       if(UseOldEnergy)then 
          temp = (totalSourceNumRho*kTp0 + totalPSNumRho*T300*20.) &
               -(totalLossNumx+totalRLNumRhox)*State_VGB(P_,i,j,k,iBlock)
               

          SE(i,j,k) = SE(i,j,k) + (inv_gm1*temp-0.50*uu2*(totalLossRho)     ) 
          SP(i,j,k) = SP(i,j,k) + (temp        +0.50*uu2*(totalSourceRho)*gm1)

          SE(i,j,k) = SE(i,j,k)  &
               -0.5*State_VGB(rho_,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock) 

          SP(i,j,k) = SP(i,j,k)  &
               +gm1*0.5*State_VGB(rho_,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock) 

          if(kTi > kTn)then
             SE(i,j,k) = SE(i,j,k)  &
                  -nu_BLK(i,j,k,iBlock)*totalNumRho*inv_gm1&
                  *(kTi-kTn)
             SP(i,j,k) = SP(i,j,k)  &
                  -nu_BLK(i,j,k,iBlock)*totalNumRho &
                  *(kTi-kTn)
          end if
       else
          temps = totalSourceNumRho*kTn            &
               + totalNumRho*(kTn-KTi)*nu_BLK(i,j,k,iBlock) &
               + totalPSNumRho*kTe0                &
               - totalLossNumRho*kTi               &
               - totalRLNumRhox*totalNumRho*KTe
          
!!!          if(UseTempControl.and.kTi > kT300)&
!!!               temps = temps+totalNumRho*(kT1000-KTi)*Nu_C(i,j,k)*5.0
          
          SE(i,j,k) = SE(i,j,k)  &
               -0.5*State_VGB(rho_,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock)  &
               -0.50*uu2*(totalLossRho) &
               +inv_gm1*temps
             
          SP(i,j,k) = SP(i,j,k)  &
               +0.5*gm1*State_VGB(rho_,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock)  &
               +0.50*(gm1)*uu2*(totalSourceRho) &
               +temps
          
       end if
                
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
    !    Rbody = 1.0 + 140.0e3/Mars
    BodyRho_I(1) = sum(BodyRhoSpecies_I(1:MaxSpecies))
    BodyP_I(1)   = sum(BodyRhoSpecies_I(1:MaxSpecies)&
         /MassSpecies_I(1:MaxSpecies))*kTp0

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

  !===========================================================================

  subroutine user_set_ICs
    use ModProcMH, ONLY : iProc
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
    use ModIO, ONLY : restart
    use ModPhysics
    use ModNumConst

    real :: Rmax, SinSlope, CosSlope,CosSZA
    real :: B4, dB4dx, zeta4, q4, epsi4, plobe, &
         XFace, YFace, ZFace
    integer :: i,j,k
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
       else if(R_BLK(i,j,k,globalBLK)< 3.0) then
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=&
               BodynDenNuSpecies_I(:)* & 
               exp(-(R_BLK(i,j,k,globalBLK)-Rbody)&
               /HNuSpecies_I(:))
       else
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=0.0
       end if
    end do; end do; end do

!    if(okTestMe)then
!       write(*,*)'nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,1:nNuSPecies)=',&
!            nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,1:nNuSPecies) 
!       WRITE(*,*)''
!       write(*,*)'nu(testcell)=', nu_BLK(itest,jtest,ktest,BLKtest)
!       WRITE(*,*)''
!    end if



!    call neutral_density_averages  !calculate averaged neutral density

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

!    if(okTestMe)then
!       write(*,*)'nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,:)=',&
!            nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,:) 
!       WRITE(*,*)''
!       write(*,*)'Productrate_CB(testcell)=',&
!            Productrate_CB(itest,jtest,ktest,BLKtest)
!       write(*,*)''  
!    end if

    do k=1,nK; do j=1,nJ; do i=1,nI
       if(UseHotO) then
          nu_BLK(i,j,k,globalBLK)=&
               sum(nDenNuSpecies_CBI(i,j,k,globalBLK,:))*nu0

          nDenNuSpecies_CBI(i,j,k,globalBLK,O_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,Ox_)
          
          nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,CO2x_)
          
          nDenNuSpecies_CBI(i,j,k,globalBLK,O_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,Oh_)+&
               nDenNuSpecies_CBI(i,j,k,globalBLK,Ohx_)
          
          nDenNuSpecies_CBI(i,j,k,globalBLK,H_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,H_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,Hx_)

       else
          nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,CO2x_)
          
          nDenNuSpecies_CBI(i,j,k,globalBLK,O_)= &
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_)+ &
               nDenNuSpecies_CBI(i,j,k,globalBLK,Ox_)
          
          nu_BLK(i,j,k,globalBLK)=(nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)+&
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_))*nu0
          
          nDenNuSpecies_CBI(i,j,k,globalBLK,H_)= 1.0e-5
          
       end if

    end do; end do; end do 


    if(UseMarsAtm)then
       if(maxval(R_BLK(:,:,:,globalBLK))<3.0*Rbody) call Mars_input
 
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(UseHotO) then
             nDenNuSpecies_CBI(i,j,k,globalBLK,Oh_)= &
                  nDenNuSpecies_CBI(i,j,k,globalBLK,Oh_)+&
                  nDenNuSpecies_CBI(i,j,k,globalBLK,Ohx_)
             
             nDenNuSpecies_CBI(i,j,k,globalBLK,O_)= &
                  nDenNuSpecies_CBI(i,j,k,globalBLK,O_)+ &
                  nDenNuSpecies_CBI(i,j,k,globalBLK,Oh_)
             
             nu_BLK(i,j,k,globalBLK)=(nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)+&
                  nDenNuSpecies_CBI(i,j,k,globalBLK,O_)+&
                  nDenNuSpecies_CBI(i,j,k,globalBLK,H_) )*nu0
          else
              
             nu_BLK(i,j,k,globalBLK)=(nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)+&
                  nDenNuSpecies_CBI(i,j,k,globalBLK,O_))*nu0
             
             nDenNuSpecies_CBI(i,j,k,globalBLK,H_)= 1.0e-5
             
          end if
          
          Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=&
               Ionizationrate_CBI(i,j,k,globalBLK,CO2_)*&
               nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)
          Ionizationrate_CBI(i,j,k,globalBLK,O_)=&
               Ionizationrate_CBI(i,j,k,globalBLK,O_)*&
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
          
       end do; end do; end do 
    else
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ionizationrate_CBI(i,j,k,globalBLK,O_)= &
               Rate_I(O_hv__Op_em_)&
               *nDenNuSpecies_CBI(i,j,k,globalBLK,O_)&
               *Productrate_CB(i,j,k,globalBLK)

          Ionizationrate_CBI(i,j,k,globalBLK,CO2_)= &
               Rate_I(CO2_hv__CO2p_em_)&
               *nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)&
               *Productrate_CB(i,j,k,globalBLK)
       end do;end do; end do
    end if
    nu1_BLK(:,:,:,globalBLK)=nu_BLK(:,:,:,globalBLK)
    
    if(okTestMe)then
       write(*,*)'usehoto=',UseHotO
       write(*,*)'nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,:)=',&
            nDenNuSpecies_CBI(itest,jtest,ktest,BLKtest,:) 
       WRITE(*,*)''
       write(*,*)'nu(testcell)=', nu_BLK(itest,jtest,ktest,BLKtest)
       WRITE(*,*)''
       write(*,*)'Ionizationrate_CBI(testcell,CO2_)=',&
            Ionizationrate_CBI(itest,jtest,ktest,BLKtest,CO2_)
       write(*,*)'Ionizationrate_CBI(testcell,O_)=',&
            Ionizationrate_CBI(itest,jtest,ktest,BLKtest,O_)
       write(*,*)''
    end if

    do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
       if (R_BLK(i,j,k,globalBLK)< Rbody) then
          cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
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
          State_VGB(P_,i,j,k,globalBLK) = &
               max(SW_p, sum(State_VGB(rho_+1:rho_+MaxSpecies,i,j,k,globalBLK)&
               /MassSpecies_I(1:MaxSpecies))*kTp0 )
          
       else
          State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)
          State_VGB(Ux_:bz_,i,j,k,globalBLK)   =0.0          
       end if
    end do;end do; end do;
    
    if(OkTestMe)&
         write(*,*)'state_VGB(body1_)=',&
         CellState_VI(:,body1_),'cell_state_VI(:,1)=',CellState_VI(:,1)
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       
       if (true_cell(i,j,k,globalBLK).and. &
            R_BLK(i,j,k,globalBLK)<1.5*Rbody) then
          
          cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
               x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
               1.0e-3
          
          State_VGB(rhoCO2p_,i,j,k,globalBLK)= &
               Ionizationrate_CBI(i,j,k,globalBLK,CO2_) &
               /nDenNuSpecies_CBI(i,j,k,globalBLK,O_)   &
               /(Rate_I(CO2p_O__O2p_CO_)+Rate_I(CO2p_O__Op_CO2_))
          
          State_VGB(rhoOp_,i,j,k,globalBLK)= &
               (Ionizationrate_CBI(i,j,k,globalBLK,O_) &
               +Rate_I(CO2p_O__Op_CO2_)                &
               *State_VGB(rhoCO2p_,i,j,k,globalBLK)    &
               *nDenNuSpecies_CBI(i,j,k,globalBLK,O_)) &
               /(nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)+4.0e6)&
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
            /MassSpecies_I(1:MaxSpecies))*kTp0)
    end do; end do; end do
    
    
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !THIS SUBROUTINE calculate the scale height of ion and neutal species and 
  !intial boundary value of ion species
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_multiSp_ICs
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
       write(*,*)'kTp=',SW_p*Tnu_body_dim*2.0/SW_T_dim, &
            2.0*Tnu_body_dim/No2Si_V(UnitTemperature_)
       write(*,*)'BodynDenNuSpecies_dim_I(:)',&
            BodynDenNuSpdim_I(:)
    end if

    select case(SolarCond)

    case('solarmax')
       Tnu_body_dim = 134.0      ! neutral temperature 
       BodynDenNuSpDim_I(CO2_)= 4.435e12
       BodynDenNuSpDim_I(O_)= 8.0283e9
       BodynDenNuSpDim_I(H_)= 1.8374e6
       BodynDenNuSpDim_I(Oh_)= 6.3119e4
       BodynDenNuSpDim_I(Ohx_)= 3.9646e3
       BodynDenNuSpDim_I(Hx_)= 7.3638e4
       BodynDenNuSpDim_I(Ox_)= 5.1736e8
       BodynDenNuSpDim_I(CO2x_)= 8.0807e10
      
       HNuSpeciesDim_I(O_)=13.34 !scale height in KM
       HNuSpeciesDim_I(Ox_)=50.025 
       HNuSpeciesDim_I(Oh_)=290.5
       HNuSpeciesDim_I(Ohx_)=2436.6

       HNuSpeciesDim_I(CO2_)=6.5631
       HNuSpeciesDim_I(CO2x_)=17.064

       HNuSpeciesDim_I(H_)=13.133
       HNuSpeciesDim_I(Hx_)=610.0
       
       RateDim_I(CO2_hv__CO2p_em_)=7.3e-7
       RateDim_I(O_hv__Op_em_) = 2.734e-7
       RateDim_I(H_hv__Hp_em_) = 8.59e-8

    case('solarmin')  ! for solar min condition

       Tnu_body_dim = 117.0      ! neutral temperature 

       BodynDenNuSpDim_I(CO2_)= 1.1593e12
       BodynDenNuSpDim_I(O_)= 3.2278e9
       BodynDenNuSpDim_I(H_)= 1.1307e7
       BodynDenNuSpDim_I(Oh_)= 1.951e4
       BodynDenNuSpDim_I(Ohx_)= 1.5248e3
       BodynDenNuSpDim_I(Hx_)= 9.4936e5
       BodynDenNuSpDim_I(Ox_)= 5.2695e8
       BodynDenNuSpDim_I(CO2x_)= 2.2258e11

       HNuSpeciesDim_I(O_)=9.486  !scale height in km
       HNuSpeciesDim_I(Ox_)=30.45  
       HNuSpeciesDim_I(Oh_)=290.5
       HNuSpeciesDim_I(Ohx_)=2436.6

       HNuSpeciesDim_I(CO2_)=5.2667
       HNuSpeciesDim_I(CO2x_)=10.533

       HNuSpeciesDim_I(H_)=13.133
       HNuSpeciesDim_I(Hx_)=586.6

       RateDim_I(CO2_hv__CO2p_em_)=2.47e-7
       RateDim_I(O_hv__Op_em_) = 8.89e-8
       RateDim_I(H_hv__Hp_em_) = 5.58e-8

    case ('earlymars1')
       Tnu_body_dim = 180.0      ! neutral temperature 
       !increase to 1000k to exobase

       RateDim_I(CO2_hv__CO2p_em_)=6.6e-7*6.0/2.25
       RateDim_I(O_hv__Op_em_) = 2.0e-7*6.0/2.25
       RateDim_I(H_hv__Hp_em_) = 8.59e-8*6.0/2.25       

       BodynDenNuSpDim_I(CO2_)= 5.1e11
       HNuSpeciesDim_I(CO2_)=11.0  !scale height in km
       BodynDenNuSpDim_I(O_)= 3.8e10
       HNuSpeciesDim_I(O_)=2.3e5

       BodynDenNuSpDim_I(CO2_)= 5.1e11
       BodynDenNuSpDim_I(O_)= 3.8e10
       BodynDenNuSpDim_I(H_)= 2.3e5
       BodynDenNuSpDim_I(Hx_)= 5.0e4
       BodynDenNuSpDim_I(Ox_)= 2.2e9
       BodynDenNuSpDim_I(CO2x_)= 6.8e9

       HNuSpeciesDim_I(O_)=24.5  !scale height in km
       HNuSpeciesDim_I(Ox_)=100.
       HNuSpeciesDim_I(CO2_)=11.0
       HNuSpeciesDim_I(CO2x_)=32.0
       HNuSpeciesDim_I(H_)=24.0
       HNuSpeciesDim_I(Hx_)=650.

    case ('earlymars2')
       Tnu_body_dim = 180.0      ! neutral temperature 
       !increase to 1000k to exobase

       RateDim_I(CO2_hv__CO2p_em_)=6.6e-7*6.0/2.25
       RateDim_I(O_hv__Op_em_) = 2.0e-7*6.0/2.25
       RateDim_I(H_hv__Hp_em_) = 8.59e-8*6.0/2.25       

       BodynDenNuSpDim_I(CO2_)= 0.0
       HNuSpeciesDim_I(CO2_)=0.0  !scale height in km
       BodynDenNuSpDim_I(O_)= 0.0
       HNuSpeciesDim_I(O_)=0.0

    case default
       call stop_mpi('unknow solar condition',SolarCond)
    end select

    kTn = TNu_body_dim*Si2No_V(UnitTemperature_)
    kTi0 = kTn
    kTp0 = 2.0*kTn

    kTe0=max(Te_new_dim, Tnu_body_dim)*Si2No_V(UnitTemperature_)   

    T300 = T300_dim*Si2No_V(UnitTemperature_)

    if(oktest)then
       write(*,*)'Tnu_body=',kTn, TNu_body_dim
       write(*,*)'T300=', T300, T300_dim
       write(*,*)'Tp_body=', kTp0
    end if

    nu0=nu0_dim*No2Io_V(UnitN_)*No2Io_V(UnitT_)
    BodynDenNuSpecies_I(:)=&
         BodynDenNuSpDim_I(:)*Io2No_V(UnitN_)
    HNuSpecies_I(:)=&
         HNuSpeciesDim_I(:)*1.0e3*Si2No_V(UnitX_)

    ! normlize the reaction rate
    Rate_I(CO2_hv__CO2p_em_)= &
         Ratedim_I(CO2_hv__CO2p_em_)*No2Io_V(UnitT_)
    Rate_I(O_hv__Op_em_)=  &
         Ratedim_I(O_hv__Op_em_)*No2Io_V(UnitT_)
    Rate_I(CO2p_O__O2p_CO_)=  &
         Ratedim_I(CO2p_O__O2p_CO_)  &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(Op_CO2__O2p_CO_)=  &
         Ratedim_I(Op_CO2__O2p_CO_)*exp(log(8.0/3.0*T300/kTn)*0.39) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_O__Op_CO2_)=  &
         Ratedim_I(CO2p_O__Op_CO2_) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)

    Rate_I(O2p_em__O_O_)=  &
         Ratedim_I(O2p_em__O_O_)*exp(log(4.0*T300/kTn)*0.56)&
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_em__CO_O_)=  &
         Ratedim_I(CO2p_em__CO_O_)*sqrt(T300/kTn)&
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


    CrossSection_I=CrossSectiondim_I*No2Io_V(unitN_)*No2Si_V(unitX_)*1.0e2

    Optdep =  sum(BodynDenNuSpecies_I*CrossSection_I*HNuSpecies_I)
    Productrate0 = max(exp(-Optdep), 1.0e-5)

    if(oktest)then
       write(*,*)'=======in set_multisp=============='
       write(*,*)'BodynDenNuSpecies_I=',BodynDenNuSpecies_I
       write(*,*)'HNuSpecies_I=',HNuSpecies_I
       write(*,*)'solar min, Procductrate=', productrate0, Optdep
       write(*,*)'CrossSection_dim_I*unitUSER_n*unitSI_x=',CrossSectiondim_I,&
            No2Io_V(unitN_),No2Si_V(unitX_)  
       write(*,*)''
    end if

    !ion density at the body
    BodyRhoSpecies_I(Hp_)=SW_rho*0.3

    BodyRhoSpecies_I(CO2p_)= Rate_I(CO2_hv__CO2p_em_)*Productrate0*&
         BodynDenNuSpecies_I(CO2_)/BodynDenNuSpecies_I(O_)/&
         (Rate_I(CO2p_O__O2p_CO_)+Rate_I(CO2p_O__Op_CO2_))
    BodyRhoSpecies_I(Op_)= (Rate_I(O_hv__Op_em_)*Productrate0+&
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
    use ModNumConst	

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

    use ModSize,       ONLY: nDim,West_,North_,Top_	
    use ModMain,       ONLY: UseRotatingBc
    use ModVarIndexes, ONLY: nVar, RhoOp_, RhoO2p_, RhoCO2p_, RhoHp_
    use ModPhysics,    ONLY: SW_rho, SW_p, SW_T_dim
    use ModFaceBc,     ONLY: FaceCoords_D, VarsTrueFace_V

    real, intent(out):: VarsGhostFace_V(nVar)

    real:: XFace,YFace,ZFace,rFace,rFace2
    real:: v_phi(3)
    real:: cosSZA 
    real:: uDotR, bDotR
    !--------------------------------------------------------------------------

    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)

    rFace2 = XFace**2 + YFace**2 + ZFace**2
    rFace  = sqrt(rFace2)

    !Apply boundary conditions
    cosSZA=(0.5+sign(0.5,XFace)) * XFace/max(RFace,1.0e-3) + 1.0e-3

    VarsGhostFace_V(rhoOp_)  = BodyRhoSpecies_I(Op_) *cosSZA

    VarsGhostFace_V(rhoO2p_) = BodyRhoSpecies_I(O2p_)*sqrt(cosSZA)

    VarsGhostFace_V(rhoCO2p_)= BodyRhoSpecies_I(CO2p_)*cosSZA

    VarsGhostFace_V(rhoHp_)  = SW_rho*0.3


    VarsGhostFace_V(rho_) = sum(VarsGhostFace_V(rho_+1:rho_+MaxSpecies))
    VarsGhostFace_V(P_)=sum(VarsGhostFace_V(rho_+1:rho_+MaxSpecies)&
         /MassSpecies_I)*kTp0

    ! Reflective in radial direction
    uDotR = sum(VarsTrueFace_V(Ux_:Uz_)*FaceCoords_D)/rFace2
    ! bDotR = sum(VarsTrueFace_V(Bx_:Bz_)*FaceCoords_D)/rFace2

    VarsGhostFace_V(Ux_:Uz_) = VarsTrueFace_V(Ux_:Uz_) - 2*uDotR*FaceCoords_D
    ! VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_) - 2*bDotR*FaceCoords_D
    VarsGhostFace_V(Bx_:Bz_) = 0.0

    ! Apply corotation?
    if (UseRotatingBc) then
       call calc_corotation_velocities(FaceCoords_D, v_phi)
       VarsGhostFace_V(Ux_:Uz_) = VarsGhostFace_V(Ux_:Uz_) + 2*v_phi
    end if

  end subroutine user_face_bcs

  !====================================================================
  !                     neutral_density_averages
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
    !  real :: neutral_density
    !true_cell note: using true_cell to replace an Rbody test does not apply here
    !----------------------------------------------------------------

    do k=1,nK; do j=1,nJ;  do i=1,nI  
       VInv=vInv_CB(i,j,k,globalBLK)

       if(.not.true_cell(i,j,k,globalBLK))cycle
       !-------------------East----------------------------------
       x = cHalf*(x_BLK(i-1,j,k,globalBLK) + x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i-1,j,k,globalBLK) + y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i-1,j,k,globalBLK) + z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,East_)= FaceAreaI_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,East_)*x+ &
            FaceArea_DS(2,East_)*y+ &
            FaceArea_DS(3,East_)*z)/R0
       do iNu = 1, nNuSpecies 
          density_IS(East_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------West----------------------------------
       x = cHalf*(x_BLK(i+1,j,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i+1,j,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i+1,j,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,West_)= FaceAreaI_DFB(:,i+1,j,k,globalBLK)
       factor = (FaceArea_DS(1,West_)*x+ &
            FaceArea_DS(2,West_)*y+ &
            FaceArea_DS(3,West_)*z)/R0     
       do iNu = 1, nNuSpecies 
          density_IS(West_,iNu) =-neutral_density(R0,iNu)*factor
       end do

       !-------------------South----------------------------------
       x = cHalf*(x_BLK(i,j-1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i,j-1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i,j-1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,South_)=FaceAreaJ_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,South_)*x+ &
            FaceArea_DS(2,South_)*y+ &
            FaceArea_DS(3,South_)*z)/R0  
       do iNu = 1, nNuSpecies 
          density_IS(South_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------North----------------------------------
       x = cHalf*(x_BLK(i,j+1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i,j+1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i,j+1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,North_)=FaceAreaJ_DFB(:,i,j+1,k,globalBLK)
       factor = (FaceArea_DS(1,North_)*x+ &
            FaceArea_DS(2,North_)*y+ &
            FaceArea_DS(3,North_)*z)/R0     
       do iNu = 1, nNuSpecies 
          density_IS(North_,iNu) = -neutral_density(R0,iNu)*factor
       end do

       !-------------------Bot----------------------------------
       x = cHalf*(x_BLK(i,j,k-1,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i,j,k-1,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i,j,k-1,globalBLK)+z_BLK(i,j,k,globalBLK))
       R0 = sqrt(x*x + y*y + z*z+cTolerance**2)
       FaceArea_DS(:,Bot_)= FaceAreaK_DFB(:,i,j,k,globalBLK)
       factor = (FaceArea_DS(1,Bot_)*x+ &
            FaceArea_DS(2,Bot_)*y+ &
            FaceArea_DS(3,Bot_)*z)/R0
       do iNu = 1, nNuSpecies 
          density_IS(Bot_,iNu) = neutral_density(R0,iNu)*factor
       end do

       !-------------------Top----------------------------------
       x = cHalf*(x_BLK(i,j,k+1,globalBLK)+x_BLK(i,j,k,globalBLK))
       y = cHalf*(y_BLK(i,j,k+1,globalBLK)+y_BLK(i,j,k,globalBLK))
       z = cHalf*(z_BLK(i,j,k+1,globalBLK)+z_BLK(i,j,k,globalBLK))
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
    use ModPhysics, ONLY :Rbody,cZero

    real, intent(in) :: R0
    integer, intent(in) :: iNu

    !-----------------------------------------------------------------------
    neutral_density = cZero
    if( R0 >= 0.9*Rbody .and. R0< 3.0*Rbody ) &
         neutral_density= exp(-(R0-Rbody)/HNuSpecies_I(iNu))

  end function neutral_density

  !===================================================================
  subroutine user_get_b0(X1,Y1,Z1,B1)
    use ModMain
    use ModPhysics 
    use ModNumConst
    implicit none
  
    real, intent(in) :: X1,Y1,Z1
    real, intent(out), dimension(3) :: B1
    
    real :: R0, theta, phi, rr, X0, Y0, Z0
    real, dimension(3) :: bb, B0
    real :: sint, sinp, cost, cosp, uB

    !-------------------------------------------------------------------------
    call timing_start('user_get_b0')

    X0 = X1*cos(thetilt)-Z1*sin(thetilt)
    Y0 = Y1
    Z0 = X1*sin(thetilt)+Z1*cos(thetilt)
    
    R0 = sqrt(X0*X0 + Y0*Y0 + Z0*Z0)
    rr = max(R0, 1.00E-6)
    if(abs(X0).lt.1e-6) then
       if(Y0.lt.0) then
          phi=-cPi/2.
       else
          phi=cPi/2.
       endif
    else 
       if(X0.gt.0) then		
          phi=atan(Y0/X0)
       else 
          phi=cPi+atan(Y0/X0)
       endif
    endif
   
    !rot=cPi	
   	
    theta=acos(Z0/rr)

    call MarsB0(R0,theta, phi+rot, bb)

    sint=sin(theta)
    cost=cos(theta)
    sinp=sin(phi)
    cosp=cos(phi)

    B0(1) = bb(1)*sint*cosp+bb(2)*cost*cosp-bb(3)*sinp
    B0(2) = bb(1)*sint*sinp+bb(2)*cost*sinp+bb(3)*cosp 
    B0(3) = bb(1)*cost-bb(2)*sint
    
    B1(1) = B0(1)*cos(thetilt)+B0(3)*sin(thetilt)
    B1(2) = B0(2)
    B1(3) = -B0(1)*sin(thetilt)+B0(3)*cos(thetilt)


    !unit of magnetic field is uB=1.677600/0.263661
    !write(*,*)'unit=',No2Io_V(UnitB_)
    uB=No2Io_V(UnitB_)
    !   write(*,*)'UB=',uB
    !uB=6.3627
    B1(1)=B1(1)/uB
    B1(2)=B1(2)/uB
    B1(3)=B1(3)/uB

    call timing_stop('user_get_b0')
  end subroutine user_get_b0
  !===========================================================================
  subroutine MarsB0(r,theta, phi, bb)
    implicit none  
    
    integer, parameter:: nMax=62
    real, intent(in) :: r, theta, phi
    real, dimension(1:3),intent(out) :: bb
    !real :: Rlgndr, dRlgndr
    integer :: NN, n, m, im
    real :: dRnm, signsx, Rmm
    real :: xtcos,xtsin,xtabs, xx
    real, dimension(0:nMax-1) :: xpcos, xpsin
    real :: a, arr, arrn, arrm, somx2, fact, temp
    real,dimension(0:nMax,0:nMax) :: Rnm
    real,dimension(0:nMax), save  :: Factor1_I, Factor2_I, Factor3_I
    real,dimension(0:nMax,0:nMax), save :: Factor1_II, Factor2_II, Factor3_II
    logical :: DoSetFactor = .true.

    !-------------------------------------------------------------------------
    if(DoSetFactor)then
       DoSetFactor = .false.
       do m = 0, nMax
          Factor1_I(m) = sqrt((2.*m+2.)*(2.*m+1.))
          Factor2_I(m) = sqrt(4.*m+6.)
          Factor3_I(m) = sqrt(2.*m+5.)
          do n = m, nMax
             if(n>m+2)then
                temp= sqrt((n-1.-m)*(n+m+1.)/(2.*n+1.)/(n-m-2.))
                Factor1_II(n,m) = sqrt((2.*n-1.)/(n-m-2.))/temp
                Factor2_II(n,m) = sqrt((n+m)/(2.*n-3.))/temp
             end if
             Factor3_II(n,m) = sqrt((n-m)*(n+m+1.))
          end do
       end do
    end if

    a=1.035336
    arr=a/r
    
    !NNm=8
       
    mars_eps=1e-3
       
    if(r.lt.1.0) then 
       NN=0       
    else 
       NN=NNm-1
    endif
	        
    xtcos=cos(theta)
    xtsin=sin(theta)
                            
    do im=0,NN
       xpcos(im)=cos(im*phi)
       xpsin(im)=sin(im*phi)
    end do
    
    bb(1)=0.0
    bb(2)=0.0
    bb(3)=0.0		
	
    !	    somx2=sqrt((1.-xtcos)*(1.+xtcos))
    somx2=abs(xtsin)
    signsx=sign(1., xtsin)
    
    fact=1.
    Rmm=1.
    Rnm(0,0)=sqrt(2.)
    Rnm(1,0)=xtcos*sqrt(3.)*Rnm(0,0)
    do n=2, NN
       Rnm(n, 0)=(xtcos*sqrt((2.*n-1.)*(2.*n+1.))*Rnm(n-1,0)-&
            (n-1)*sqrt((2.*n+1.)/(2.*n-3.))*Rnm(n-2, 0))/n
       
    enddo !n
	    
    arrm=1.0

    call timing_start('crustal')

    do m=0, NN

       Rmm=Rmm*fact*somx2/Factor1_I(m)
       
       Rnm(m+1,m+1)=Rmm*Factor2_I(m) 
       Rnm(m, m+1)=0
       
       fact=fact+2.
	       
       arrm=arr*arrm
       arrn=arrm
	       
       do n=m,NN
          arrn=arr*arrn
          !write(*,*) 'arrn=', arrn, ' n=', n
          if(n> (m+2)) then
             Rnm(n,m+1) = xtcos*Factor1_II(n,m)*Rnm(n-1,m+1)-&
                  Factor2_II(n,m)*Rnm(n-2,m+1)
             
          else if(n > (m+1)) then
             Rnm(n,m+1)=xtcos*Factor3_I(m)*Rnm(m+1,m+1)
          endif

          dRnm=m*xtcos*Rnm(n,m)/xtsin-Rnm(n, m+1)*signsx* Factor3_II(n,m)
          
          bb(1)=bb(1)+(n+1)*arrn*Rnm(n,m)*(cmars(n,m)*xpcos(m)&
               +dmars(n,m)*xpsin(m))
          bb(2)=bb(2)-arrn*dRnm*(cmars(n,m)*&
               xpcos(m)+dmars(n,m)*xpsin(m))
          if(xtsin <= 1e-6) then
             bb(3)=0.
          else	
             bb(3)=bb(3)-arrn*Rnm(n,m)*m/xtsin*(-cmars(n,m)*xpsin(m)&
                  +dmars(n,m)*xpcos(m))
          endif
       end do !n
    end do !m

    call timing_stop('crustal')

  end subroutine MarsB0

!=====================================================================
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

  !============================================================================
  subroutine Mars_Input
    use ModMain
    use ModPhysics
    use ModConst
    use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
         XyzStart_BLK,TypeGeometry
    implicit none
    
    real, parameter :: TINY=1.0E-12 
    real :: hh, theta, phi, dR, dtheta, dphi, dH, Hscale, HCO2, HO, grav
    real:: tempICO2p, tempIOp
    real:: xLat, xLong,xAlt
    integer :: i,j,k,n, m
    integer:: iAlt, jLong, kLat, ip1,jp1,kp1
    logical:: oktest, oktestme=.true.
    !------ Interpolation/Expolation for Tn,nCO2,nO,PCO2p,POp ----- 
    
    dR=dx_BLK(globalBLK)
    dPhi=dy_BLK(globalBLK)
    dTheta=dz_BLK(globalBLK)

    select case(TypeGeometry)                                   
    case('cartesian')                                           
       call stop_mpi('Unknown geometry type = '//TypeGeometry)
       
    case('spherical','spherical_lnr')
       ! at least part of the block is outside the body 
       if (R_BLK(nI,1,1,globalBLK) >= Rbody) then  
          
          do k=1,nK
             Theta = (k-1)*dTheta  + xyzStart_BLK(Theta_,globalBLK)      
             Theta =  180*(0.5-Theta/cPi)
             kLat=int((theta+87.5)/5.0+1.0)
             kp1=min(kLat+1, NLat)
             kLat = max(kLat,1)
             
             do j=1,nJ  
                Phi = (j-1)*dPhi  + xyzStart_BLK(Phi_,globalBLK)
                if(phi>cPi)then 
                   phi=phi-2*cPi
                end if
                Phi = 180*(Phi/cPi) 
                jLong=int((phi+180)/5.0+1.0)                 
                jp1=min(jLong+1,NLong)
                jLong=max(jLong,1)
                
                do i=nI,1,-1                    
                   hh = (R_BLK(i,j,k,globalBLK)-1.00)*3396.00
                   !                 write(*,*)'hh=', hh, i,j,k,globalBLK
                   xLong=0.2*(Phi-Long_I(jLong))
                   xLat=0.2*(Theta-Lat_I(kLat))
                   if(hh.le.100.0)then  !inside the body
                      tempNuSpecies_CBI(i,j,k,globalBLK)= &
                           tempNuSpecies_CBI(i+1,j,k,globalBLK)
                      nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)=&
                           nDenNuSpecies_CBI(i+1,j,k,globalBLK,CO2_)
                      nDenNuSpecies_CBI(i,j,k,globalBLK,O_)= &
                           nDenNuSpecies_CBI(i+1,j,k,globalBLK,O_)
                      
                      !                    tempICO2p=max(tempICO2p,TINY)
                      !                    tempIOP=max(tempIOp,TINY)
                      Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=&
                           Ionizationrate_CBI(i+1,j,k,globalBLK,CO2_)
                      Ionizationrate_CBI(i,j,k,globalBLK,O_)=&
                           Ionizationrate_CBI(i+1,j,k,globalBLK,O_)
                   elseif(hh.le.Alt_I(NAlt))then
                      iAlt=int((hh -100.0)/10.0+1.0)
                      ip1=min(iAlt+1,NAlt)
                      if(iAlt.lt.1)then 
                         write(*,*)'wrong ialt',iAlt
                      end if
                      xalt=0.1*(hh-Alt_I(iAlt))
                      !interpolate
                      tempNuSpecies_CBI(i,j,k,globalBLK)=          &
                           ((Temp(jLong,kLat,iAlt)*(1-xLong)       &
                           + xLong*Temp(jp1, kLat, ialt))*(1-xLat) &
                           +(Temp(jLong,kp1,iAlt)*(1-xLong)        &
                           + xLong*Temp(jp1, kp1, ialt))*xLat)*(1-xAlt)&
                           +((Temp(jLong,kLat,ip1)*(1-xLong)       &
                           + xLong*Temp(jp1, kLat, ip1))*(1-xLat)   &
                           +(Temp(jLong,kp1,ip1)*(1-xLong)         &
                           + xLong*Temp(jp1, kp1, ip1))*xLat)*xAlt
                    
                      nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)=&         
                           ((Den_CO2(jLong,kLat,iAlt)*(1-xLong)+xLong*Den_CO2(jp1, kLat, ialt))*(1-xLat)+&
                           (Den_CO2(jLong,kp1,iAlt)*(1-xLong)+xLong*Den_CO2(jp1, kp1, ialt))*xLat)*(1-xAlt)+&
                           ((Den_CO2(jLong,kLat,ip1)*(1-xLong)+xLong*Den_CO2(jp1, kLat, ip1))*(1-xLat)+&
                           (Den_CO2(jLong,kp1,ip1)*(1-xLong)+xLong*Den_CO2(jp1, kp1, ip1))*xLat)*xAlt
                       
                      nDenNuSpecies_CBI(i,j,k,globalBLK,O_)=&
                           ((Den_O(jLong,kLat,iAlt)*(1-xLong)+xLong*Den_O(jp1, kLat, ialt))*(1-xLat)+&
                           (Den_O(jLong,kp1,iAlt)*(1-xLong)+xLong*Den_O(jp1, kp1, ialt))*xLat)*(1-xAlt)+&
                           ((Den_O(jLong,kLat,ip1)*(1-xLong)+xLong*Den_O(jp1, kLat, ip1))*(1-xLat)+&
                           (Den_O(jLong,kp1,ip1)*(1-xLong)+xLong*Den_O(jp1, kp1, ip1))*xLat)*xAlt
                    
                      tempICO2p=&
                           ((ICO2p(jLong,kLat,iAlt)*(1-xLong)+xLong*ICO2p(jp1, kLat, ialt))*(1-xLat)+&
                           (ICO2p(jLong,kp1,iAlt)*(1-xLong)+xLong*ICO2p(jp1, kp1, ialt))*xLat)*(1-xAlt)+&
                           ((ICO2p(jLong,kLat,ip1)*(1-xLong)+xLong*ICO2p(jp1, kLat, ip1))*(1-xLat)+&
                           (ICO2p(jLong,kp1,ip1)*(1-xLong)+xLong*ICO2p(jp1, kp1, ip1))*xLat)*xAlt
                      
                      tempIOP=&
                           ((IOp(jLong,kLat,iAlt)*(1-xLong)+xLong*IOp(jp1, kLat, ialt))*(1-xLat)+&
                           (IOp(jLong,kp1,iAlt)*(1-xLong)+xLong*IOp(jp1, kp1, ialt))*xLat)*(1-xAlt)+&
                           ((IOp(jLong,kLat,ip1)*(1-xLong)+xLong*IOp(jp1, kLat, ip1))*(1-xLat)+&
                           (IOp(jLong,kp1,ip1)*(1-xLong)+xLong*IOp(jp1, kp1, ip1))*xLat)*xAlt
                      
                      tempICO2p=max(tempICO2p,TINY)
                      tempIOP=max(tempIOp,TINY)
                      !                   Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=tempICO2p*nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)
                      !                   Ionizationrate_CBI(i,j,k,globalBLK,O_)=tempIOP*nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
                      Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=tempICO2p
                      Ionizationrate_CBI(i,j,k,globalBLK,O_)=tempIOP
                   else  !hh.gt.Alt_I(NAlt)
                      
                      dH= hh - Alt_I(NAlt)
                      tempNuSpecies_CBI(i,j,k,globalBLK)= &
                           (Temp(jLong,kLat,NAlt)*(1-xLong)+xLong*Temp(jp1, kLat,NAlt))*(1-xLat)+&
                           (Temp(jLong,kp1,NAlt)*(1-xLong)+xLong*Temp(jp1,kp1,NAlt))*xLat
                      
                      tempICO2p=&
                           (ICO2p(jLong,kLat,NAlt)*(1-xLong)+xLong*ICO2p(jp1, kLat, NAlt))*(1-xLat)+&
                           (ICO2p(jLong,kp1,NAlt)*(1-xLong)+xLong*ICO2p(jp1, kp1, NAlt))*xLat
                      
                      tempIOP=&
                           (IOp(jLong,kLat,NAlt)*(1-xLong)+xLong*IOp(jp1, kLat, NAlt))*(1-xLat)+&
                           (IOp(jLong,kp1,NAlt)*(1-xLong)+xLong*IOp(jp1, kp1, NAlt))*xLat
                      
                      ! grav=3.72/R_BLK(i,j,k,globalBLK)/R_BLK(i,j,k,globalBLK)
                      grav=3.72/(1.0+300.0/3396.0)/(1.0+300.0/3396.0)
                      
                      Hscale=cBoltzmann*&
                           tempNuSpecies_CBI(i,j,k,globalBLK)/grav/cProtonMass!in m unit
                      
                      HCO2= Hscale/NuMassSpecies_I(CO2_)/1.0e3
                      HO= Hscale/NuMassSpecies_I(O_)/1.0e3
                      
                      nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)=&
                           ((Den_CO2(jLong,kLat,NAlt)*(1-xLong)+xLong*Den_CO2(jp1, kLat,Nalt))*(1-xLat)+&
                           (Den_CO2(jLong,kp1,NAlt)*(1-xLong)+xLong*Den_CO2(jp1,kp1,Nalt))*xLat)&
                           *exp(-dH/HCO2)
                      nDenNuSpecies_CBI(i,j,k,globalBLK,O_)=&
                           ((Den_O(jLong,kLat,NAlt)*(1-xLong)+xLong*Den_O(jp1, kLat,Nalt))*(1-xLat)+&
                           (Den_O(jLong,kp1,NAlt)*(1-xLong)+xLong*Den_O(jp1,kp1,Nalt))*xLat)&
                           *exp(-dH/HO)
                      
                      tempICO2p=max(tempICO2p,TINY)
                      tempIOP=max(tempIOp,TINY)
                      !                    Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=tempICO2p*nDenNuSpecies_CBI(i,j,k,globalBLK,CO2_)
                      !                    Ionizationrate_CBI(i,j,k,globalBLK,O_)=tempIOP*nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
                      Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=tempICO2p
                      Ionizationrate_CBI(i,j,k,globalBLK,O_)=tempIOP

                   end if !hh.lt.or.gt.300km
                end do
             end do
          end do
       end if
    case default
       
       call stop_mpi('Unknown geometry type = '//TypeGeometry)
       
    end select
    if(oktestme)then
       write(*,*)'Mars input end', &
            dR,dPhi,dTheta, globalBLK, &
            maxval(nDenNuSpecies_CBI(nI,:,:,globalBLK,CO2_)),&
            minval(nDenNuSpecies_CBI(nI,:,:,globalBLK,CO2_)),&
            maxval(R_BLK(nI,:,:,globalBLK)),&
            minval(R_BLK(1,:,:,globalBLK))
       write(*,*)'Mars input end2',&
            globalBLK, maxval(nDenNuSpecies_CBI(nI,:,:,globalBLK,O_)),&
            minval(nDenNuSpecies_CBI(nI,:,:,globalBLK,O_)),&
            maxval(Ionizationrate_CBI(nI,:,:,globalBLK,CO2_)),&
            minval(Ionizationrate_CBI(nI,:,:,globalBLK,O_)),&
            maxval(R_BLK(nI,:,:,globalBLK)),&
            minval(R_BLK(1,:,:,globalBLK))
    end if
  end subroutine Mars_input
                   
  !============================================================================

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
