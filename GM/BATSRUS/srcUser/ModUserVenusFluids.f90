module ModUser
  ! This is the user module for Venus 
  !WARNING: It is still a work in progresss
  use ModSize
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics,                    &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_init_point_implicit,        &
       IMPLEMENTED4 => user_init_session,               &
       IMPLEMENTED5 => user_read_inputs,                &
       IMPLEMENTED6 => user_face_bcs,                   &
       IMPLEMENTED7 => user_set_plot_var,               &
       IMPLEMENTED8 => user_get_log_var,                &
       IMPLEMENTED9 => user_set_boundary_cells         

  use ModMultiFluid

  use ModNumConst, ONLY: cSqrtHalf, cDegToRad

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: NameUserModule = &
       'Venus 5 Fluids MHD code, Dalal Najib'

  character (len=10) :: SolarCond='solarmax  '

  real    :: CollisionCoefDim = 1.0, CollisionCoef
  ! Radius within which the point implicit scheme should be used
  real :: rPointImplicit = 4.0

  !Venus stuff etc

  !logical ::  UseMultiSpecies=.true.
 logical :: UseChapman= .false.
  integer, parameter :: MaxSpecies=nIonFluid, MaxNuSpecies=5,  &
       MaxReactions=10, nNuSpecies=3
  integer :: nSpecies=4,&
       nReactions=10

  real, allocatable:: nDenNuSpecies_CBI(:,:,:,:,:), &
       TempNuSpecies_CBI(:,:,:,:),  Productrate_CB(:,:,:,:),  &
       Ionizationrate_CBI (:,:,:,:,:), MaxSiSpecies_CB(:,:,:,:),&
       MaxLiSpecies_CB(:,:,:,:),MaxSLSpecies_CB(:,:,:,:),&
       nu_BLK(:,:,:,:),nu1_BLK(:,:,:,:)

  real, dimension(MaxReactions) :: ReactionRate_I
  real, dimension(MaxReactions,nIonFluid):: CoeffSpecies_II, &
       dSdRho_II !, dLdRho_II
  real, dimension(nIonFluid)::LossSpecies_I, &
       SiSpecies_I,  LiSpecies_I,  PhoIon_I, Recb_I, ImpIon_I


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

 !      CO2_hv__Op_CO_em_=11 !CO2+hv-->Op+CO_em

  real, dimension(MaxReactions) :: Rate_I
  real, dimension(MaxReactions) :: &
       Ratedim_I=(/3.270e-6,1.224e-6, 1.64e-10, 1.1e-9, &
       9.60e-11, 7.38e-8, 3.1e-7, 5.084e-10, 6.4e-10, 0.0 /)  !cm^3 s^(-1)

       !Ratedim_I=(/ 2.47e-7, 8.89e-8, 1.64e-10, 1.1e-9, &
       !9.60e-11, 7.38e-8, 3.1e-7, 5.084e-10, 6.4e-10, 5.58e-8, 0.0 /)  !cm^3 s^(-1)

  integer, parameter :: &! order of ion species
       Hp_  =1, &
       O2p_ =2, &
       Op_  =3, &
       CO2p_=4

  character (len=10), dimension(nIonFluid):: &
       ion_name_I=(/'Hp  ', 'O2p ', 'Op  ','CO2p'/)

 real, dimension(nNuSpecies)::  &
        MassNeutral_I=(/44., 16., 1. /)  !atm

  ! real, dimension(MaxSpecies)::  &
  !      MassSpecies_I=(/1., 32., 16., 44. /)  !atm

  !  MassSpecies_I(Hp_)=1	 !atm
  !  MassSpecies_I(CO2p_)=44     !atm
  !  MassSpecies_I(O2p_)=32      !atm
  !  MassSpecies_I(Op_)=16       !atm

  integer, parameter :: & ! order of Neutral species
       CO2_=1 ,&
       O_=2   ,&    
       H_=3, &
       Oh_=4   ,&
       Ohx_=5 

!       Hx_=6, &
!       Ox_=7 ,&
!       CO2x_=8,&
!       Oh2x_ =9

  ! Given in cm^2
  real, dimension(MaxNuSpecies)::CrossSection_I,&
       CrossSectionDim_I=(/2.6e-17,1.5e-17,0.0,1.5e-17,&
       1.5e-17/)

 real, dimension(nIonFluid,nNuSpecies)::ReducedMassNeutral_II

 real, dimension(nIonFluid,nNuSpecies)::ReducedMassNeutral2_II


 real, dimension(nIonFluid,nIonFluid)::ReducedMassIon_II


  real:: Productrate0,Optdep
  real, dimension(MaxNuSpecies)::NuMassSpecies_I=(/44,16,1,16,16/)
  !  NuMassSpecies_I(CO2_)=44	!atm
  !  NuMassSpecies_I(O_)=16	!atm

  real, dimension(MaxNuSpecies):: HNuSpecies_I=1.0,HNuSpeciesDim_I=1.0
  !HNuSpecies_dim_I(CO2_)=6.7e3   !m
  !HNuSpecies_dim_I(O_)=18.4e3    !m


  real, dimension(MaxNuSpecies):: BodynDenNuSpecies_I,&
       BodynDenNuSpDim_I=(/2.5e12,7.0e10, 0.0, 1.951e4, &
       1.5248e3/) !density at 100km



  real:: Altitude0=100.0e3 !altitude correspondint to neutral density
  real, dimension(MaxSpecies):: BodyRhoSpecies_I
  integer, parameter :: & ! other numbers
       em_=-1 ,&
       hv_=-2   

 real :: TNu_body_dim, TNu_body, Tnu, Tnu_dim  !dimensionless temperature of neutral, &
                           !new created ions, plasma at the body
  real :: Te_new_dim=3000.0, KTe0 !temperature of new created e
 real :: Ti_body_dim=300.0, Ti_body  !ion temperature at the body

 real :: T300_dim = 300.0, T300, Ti_dim=300.0 
!  real,  dimension(1:nI,1:nJ,1:nK,nBLK) :: nu_BLK
  real :: nu0_dim=1.0e-9,nu0  !value from Terada 2009 Tanaka, 1997, JGR

  real:: IonNeuCoeff_II(1:nIonFluid,1:nNuSpecies)
  real::IonIonCoeff_II(1:nIonFluid,1:nIonFluid)

  real, parameter:: IonNeuCoeffDim_II(1:nIonFluid,1:nNuSpecies)= reshape( (/ &
       41.4, 5.63, 8.95, 4., &
       30., 2.31, 4., 1.76, &
       4., 0.65, 1., 0.47 /), (/4, 3/) )

   real, parameter:: IonIonCoeffDim_II(1:nIonFluid,1:nIonFluid)=reshape( (/ &
       0.90, 0.039, 0.077, 0.029, &
       1.25, 0.16, 0.26, 0.12, &
       1.23, 0.13, 0.22, 0.10,&
       1.26, 0.17, 0.30, 0.14 /), (/4, 4/) )



  logical :: UseTempCont=.false.
  logical :: UseImpactIon=.false.
  real, dimension(32,nIonFluid)::Impact_ION,Impact_ION_dim=0.0 
  real, dimension(32):: Temp_dim
  logical :: UseChargeEx=.true.
  !logical ::InitSession= .false.

  integer,parameter::NLong=73, NLat=36, MaxAlt=21
  real :: Long_I(NLong), Lat_I(NLat), Alt_I(MaxAlt), Alt0
  real :: Temp(NLong, NLat, MaxAlt)
  real :: Den_CO2(NLong, NLat, MaxAlt)!,Den_CO2_dim(NLong, NLat, NAlt)
  real :: Den_O(NLong, NLat, MaxAlt)!,Den_O_dim(NLong, NLat, NAlt)
  real :: ICO2p(NLong, NLat, MaxAlt)!,ICO2p_dim(NLong, NLat, NAlt)
  real :: IOp(NLong, NLat, MaxAlt)!,IOp_dim(NLong, NLat, NAlt)
  
  integer :: NAlt=21


  !end Venus stuff

  logical:: UseOldEnergy=.false.

contains

  !============================================================================
  subroutine user_calc_sources

    use ModPointImplicit, ONLY: UsePointImplicit_B,UsePointImplicit,&
         IsPointImplSource, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, iNewGrid, iNewDecomposition, &
         PROCTEST,GLOBALBLK,BLKTEST, iTest,jTest,kTest
    use ModPhysics, ONLY: inv_gm1,Rbody,gm1,UnitTemperature_,Si2No_V, No2Io_V, No2Si_V, LowDensityRatio,&
         UnitT_,UnitN_,ElectronPressureRatio
    use ModAdvance, ONLY: State_VGB, Source_VC,VdtFace_x,&
         VdtFace_y,VdtFace_z
    use ModGeometry,ONLY: r_BLK,x_BLK,y_BLK,z_BLK,R_BLK,&
         vInv_CB
    use ModNumConst,ONLY: cZero,cHalf,cOne,cTolerance
!!$    use ModVarIndexes,ONLY: Rho_, HpRho_, O2pRho_, OpRho_, CO2pRho_, &
!!$         RhoUx_, RhoUy_, RhoUz_, HpP_,O2pP_,OpP_,CO2pP_, P_, Energy_, Bx_, By_, Bz_
    use ModVarIndexes
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest
    use ModProcMH,   ONLY: iProc
    !    use ModAdvance,  ONLY: Source_VC,Energy_
    !    use ModNumConst, ONLY: cZero
    integer :: iBlock, i, j, k,iSpecies,iFluid
    real    :: Coef,SourcesLossMax,vdtmin
    real :: inv_rho,inv_rho2,uu2,Productrate,kTi,kTe
    real    :: NumDens, InvNumDens
    real, dimension(nIonFluid) :: NumDens_I, InvRho_I,InvRho2_I,uu2_I,&
         Ux_I, Uy_I, Uz_I, Temp_I,temps_I,&
         LossNumRho_I, SourceNumRho_I, Lossx_I,LossNumx_I,&
         RLNumRhox_I, tmp_I, col_ii_I
    real :: alt, Te_dim = 300.0, tmp, num=4.0e-11, cosSZA
    real::  totalLossRho,totalSourceRho, totalLossNumRho, &
         totalSourceNumRho, totalLossx, totalLossNumx, SourceLossMax
    real :: totalPSNumRho=0.0,totalRLNumRhox=0.0, temps, testvar
    real :: X1, X2, X3, X4, Alt0
    real :: A0=-1143.33, A1=323.767, A2=-35.0431, A3=1.69073, A4=-0.0306575
    real :: B0=-1233.29, B1=347.764, B2=-37.4128, B3=1.79337, B4=-0.0322777
    character (len=*), parameter :: NameSub = 'user_calc_sources'
    logical :: DoTest, DoTestMe, DoTestCell

    integer :: iLastGrid=-100, iLastDecomposition=-100
    !--------------------------------------------------------------------

    iBlock = GlobalBlk

    !write(*,*)'nDenNuSpecies_CBI(1, 1, 1, iBlock, MaxNuSpecies=',nDenNuSpecies_CBI(1, 1, 1, iBlock, 1)
    if( nDenNuSpecies_CBI(1, 1, 1, iBlock, 1) < 0.0 & 
         .or. iLastGrid /= iNewGrid &
         .or. iLastDecomposition /= iNewDecomposition)then
      call set_neutral_density(iBlock)
      

       if(iBlock == nBLK)then
          iLastGrid          = iNewGrid
          iLastDecomposition = iNewDecomposition
       end if

    end if

    ! Do not provide explicit source term when point-implicit scheme is used
    ! IsPointImplSource is true only when called from ModPointImplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Check if inside rPointImplicit. 
    ! Set UsePointImplicit_B=F so the point implicit
    ! evaluation is not done at all.
!!!    UsePointImplicit_B(iBlock) = 

    if(r_BLK(1,1,1,iBlock) > rPointImplicit) RETURN

!!!!    if(.not.UsePointImplicit_B(iBlock)) RETURN

    if(iProc==PROCtest.and.iBlock==BLKtest)then
       call set_oktest(NameSub,DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !    if(DoTestMe)then
    !   write(*,*)'before Source(rhoU)=', Source_VC(6:8,itest,jtest,ktest)
    !      write(*,*)NameSub,' Source(p,E)', Source_VC(iP,iTest,jTest,kTest)
    !   end if

    !chemistry etc

!    do k=1,nK; do j=1,nJ; do i=1,nI
!       cosSZA=(0.5+sign(0.5,x_BLK(i,j,k,iBlock)))*&
!            x_BLK(i,j,k,iBlock)/max(R_BLK(i,j,k,iBlock),1.0e-3)&
!            +5.0e-4
 !      Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,iBlock,1:MaxNuSpecies)*&
 !           CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
 !           6.0e-3)/cosSZA
 !      if( Optdep<11.5 .and. x_BLK(i,j,k,iBlock) > 0.0) then
 !         Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-5)
 !      else
 !         Productrate_CB(i,j,k,iBlock) = 1.0e-5
 !      end if
 !   end do; end do; end do


      do k=1,nK; do j=1,nJ; do i=1,nI
             cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
                  x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-4)&
                  +1.0e-4
             Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:MaxNuSpecies)*&
                  CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
                  6.0e-3)/cosSZA
             if( Optdep<11.5 .and. x_BLK(i,j,k,globalBLK) > 0.0) then
                Productrate_CB(i,j,k,globalBLK) = max(exp(-Optdep), 1.0e-30)
             else
                Productrate_CB(i,j,k,globalBLK) = 1.0e-30
             end if

          end do; end do; end do


    do k = 1, nK ;   do j = 1, nJ ;  do i = 1, nI

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest


       NumDens_I  = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I
       NumDens    = sum(NumDens_I)
       InvNumDens = 1.0/NumDens

       Temp_I     = State_VGB(iPIon_I,i,j,k,iBlock)/NumDens_I
       !Temp_I=TNu_body
       !State_VGB(iPIon_I,i,j,k,iBlock)=2*Temp_I*NumDens_I

       !if(DoTestCell)then
       !write(*,*)'iblock=',iblock
       !write(*,*)' in user_calc_source,State_VGB(iRhoIon_I,i,j,k,iBlock)=',State_VGB(iRhoIon_I,iTest,jTest,kTest,iBlock)
       !write(*,*)' in user_calc_source,State_VGB(iPIon_I,i,j,k,iBlock)=',State_VGB(iPIon_I,iTest,jTest,kTest,iBlock)
       !end if


       InvRho_I = 1.0/State_VGB(iRho_I(2:nFluid),i,j,k,iBlock)

       if (R_BLK(i,j,k,iBlock) < Rbody) CYCLE

       InvRho_I = 1.0/State_VGB(iRhoIon_I,i,j,k,iBlock)
       inv_rho = 1.0/sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
       inv_rho2 = inv_rho**2
       uu2 = sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2) * inv_rho2

       do iFluid = 1, nIonFluid
          uu2_I(iFluid) = &
               ( State_VGB(iRhoUxIon_I(iFluid),i,j,k,iBlock)**2  &
               + State_VGB(iRhoUyIon_I(iFluid),i,j,k,iBlock)**2  &
               + State_VGB(iRhoUzIon_I(iFluid),i,j,k,iBlock)**2) &
               * InvRho_I(iFluid)**2
       end do

       ReactionRate_I=0.0
       CoeffSpecies_II(:,:)=0.0
       PhoIon_I = 0.0
       Recb_I   = 0.0
       LossSpecies_I=0.0
       SiSpecies_I=0.0
       LiSpecies_I=0.0

       totalSourceRho=0.0
       totalLossRho=0.0
       totalLossNumRho=0.0
       totalSourceNumRho=0.0
       totalLossx=0.0
       totalLossNumx=0.0
       totalPSNumRho=0.0
       totalRLNumRhox=0.0

       LossNumRho_I=0.0
       SourceNumRho_I=0.0
       Lossx_I=0.0
       LossNumx_I=0.0
       RLNumRhox_I=0.0


       MaxSLSpecies_CB(i,j,k,iBlock)=1.0e-3

       ! calculate optical depth and producation rate


       Productrate= Productrate_CB(i,j,k,iBlock)

       if(DoTestCell)then
          write(*,*)'nDenNuSpecies_CBI=',nDenNuSpecies_CBI(iTest,jTest,kTest,iBlock,:)
          write(*,*)'Productrate_CBI=',Productrate_CB(iTest,jTest,kTest,iBlock)
          write(*,*)'Ionizationrate_CBI=',Ionizationrate_CBI(iTest,jTest,kTest,iBlock,:)
          write(*,*)' Rate_I(H_hv__Hp_em_)', Rate_I(H_hv__Hp_em_)
          write(*,*)'ReactionRate_I=',ReactionRate_I
          write(*,*)'MassIon_I=',MassIon_I
       end if

       ReactionRate_I(H_hv__Hp_em_)= &
            Rate_I(H_hv__Hp_em_)*nDenNuSpecies_CBI(i,j,k,iBlock,H_)
       PhoIon_I(Hp_)=ReactionRate_I(H_hv__Hp_em_)*Productrate
       PhoIon_I(CO2p_)=Ionizationrate_CBI(i,j,k,iBlock,CO2_)
       PhoIon_I(Op_)=Ionizationrate_CBI(i,j,k,iBlock,O_)

       !IMPACT Ionization
       ImpIon_I=0.0
       if(UseImpactIon)then
          X1=log(Te_dim)
          X2=X1*X1
          X3=X2*X1
          X4=X2*X2
          ImpIOn_I(Hp_)=exp(A0+A1*X1+A2*X2+A3*X3+A4*X4)&
               *No2Io_V(UnitT_)*No2Io_V(UnitN_)
          ImpIOn_I(Op_)=exp(B0+B1*X1+B2*X2+B3*X3+B4*X4)&
               *No2Io_V(UnitT_)*No2Io_V(UnitN_)
          ImpIon_I(Hp_)=ImpIon_I(Hp_)*&
               nDenNuSpecies_CBI(i,j,k,globalBLK,H_)
          ImpIon_I(Op_)=ImpIon_I(Op_)*&
               nDenNuSpecies_CBI(i,j,k,globalBLK,O_)
          ImpIon_I=ImpIon_I*NumDens


       end if

       State_VGB(P_,i,j,k,iBlock)=sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1+ElectronPressureRatio)
       
       kTi=State_VGB(P_,i,j,k,iBlock)/NumDens/(1+ElectronPressureRatio)
       kTe=kTi

       Ti_dim=kTi*No2Si_V(UnitTemperature_)

       !charge exchange

       ReactionRate_I(CO2p_O__O2p_CO_)= &
            Rate_I(CO2p_O__O2p_CO_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       CoeffSpecies_II(O2p_,CO2p_)=ReactionRate_I(CO2p_O__O2p_CO_)

       ReactionRate_I(Op_CO2__O2p_CO_)= &
            Rate_I(Op_CO2__O2p_CO_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)&
            *exp(log(Tnu_body_dim/(max(Temp_I(3)/Si2No_V(UnitTemperature_),T300_dim)))*0.39)
       !write(*,*)'temp_I(3)dim=',(Temp_I(3)/Si2No_V(UnitTemperature_))
       CoeffSpecies_II(O2p_, Op_)=ReactionRate_I(Op_CO2__O2p_CO_)

       ReactionRate_I(CO2p_O__Op_CO2_)= &
            Rate_I(CO2p_O__Op_CO2_)&
            * nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       CoeffSpecies_II(Op_,CO2p_)=ReactionRate_I(CO2p_O__Op_CO2_)

       ReactionRate_I(Hp_O__Op_H_)= &
       Rate_I(Hp_O__Op_H_)* nDenNuSpecies_CBI(i,j,k,iBlock,O_)
       CoeffSpecies_II(Op_,Hp_)=ReactionRate_I(Hp_O__Op_H_)

       ReactionRate_I(Op_H__Hp_O_)= &
       Rate_I(Op_H__Hp_O_)* nDenNuSpecies_CBI(i,j,k,iBlock,H_)
       CoeffSpecies_II(Hp_,Op_)=ReactionRate_I(Op_H__Hp_O_)


       ! Recombination
       kTi=State_VGB(P_,i,j,k,iBlock)/NumDens/(1+ElectronPressureRatio)
       kTe=kTi

       if(kTi <= 0.0)then
          write(*,*)'i,j,k,iBlock=',i,j,k,iBlock
          write(*,*)'xyz=',x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock)
          write(*,*)'NumDens_I=',NumDens_I
          write(*,*)'NumDens=',NumDens
          write(*,*)'State_VGB(iPIon_I,i,j,k,iBlock)=',State_VGB(iPIon_I,i,j,k,iBlock)
          write(*,*)'kTi,p=',kTi,State_VGB(P_,i,j,k,iBlock)
          call stop_mpi('DEBUG')
       end if

       ReactionRate_I(O2p_em__O_O_)=Rate_I(O2p_em__O_O_)       
       Recb_I(O2p_) = ReactionRate_I(O2p_em__O_O_)*(TNu_body/(max(Temp_I(2),T300)))**0.56

       !            ReactionRate_I(O2p_em__O_O_)*exp(log(TNu_body/Temp_I(2))*0.56)

       ReactionRate_I(CO2p_em__CO_O_)=Rate_I(CO2p_em__CO_O_)       
       Recb_I(CO2p_)=ReactionRate_I(CO2p_em__CO_O_)*&
            sqrt(TNu_body/kTi)

       !end if  !(x>0.0)

       LossSpecies_I = sum( CoeffSpecies_II, DIM=1 )

       do iFluid=1, nIonFluid
          !LossSpecies_I = LossSpecies_I + CoeffSpecies_II(iFluid, :)

          dSdRho_II(1:nIonFluid, iFluid)= &
               CoeffSpecies_II(1:nIonFluid,iFluid)*MassIon_I/MassIon_I(iFluid)

       enddo

       SiSpecies_I =(PhoIon_I+ImpIon_I)*MassIon_I

       do iFluid=1, nIonFluid
          SiSpecies_I(1:nIonFluid)=&
               SiSpecies_I(1:nIonFluid)  &
               +dSdRho_II(1:nIonFluid, iFluid) &
               *State_VGB(iRhoIon_I(iFluid), i,j,k, iBlock)


          LiSpecies_I(iFluid)= &
               LiSpecies_I(iFluid)  &
               +(LossSpecies_I(iFluid) +Recb_I(iFluid)*NumDens)&
               *State_VGB(iRhoIon_I(iFluid), i,j,k, iBlock)
       enddo


       if(DoTestCell)then
          write(*,*)NameSub,' State_VGB(iRhoIon_I)',State_VGB(iRhoIon_I,iTest,jTest,kTest,BLKTest)
          do iFluid=1,nIonFluid
              write(*,*)'iFluid=',iFluid
              write(*,*)NameSub,' CoeffSpecies_II(iFluid,:)',CoeffSpecies_II(iFluid,:)
          end do 
           write(*,*)'SiSpecies_I=',SiSpecies_I
           write(*,*)'LiSpecies_I=',LiSpecies_I
           write(*,*)'Recb_I=',Recb_I
        end if



       !SiSpecies_I(3)=SiSpecies_I(2)
       !LiSpecies_I(3)=SiSpecies_I(2)


        !totalIMPNumRho=sum(ImpIon_I(:))

       totalSourceRho=sum(SiSpecies_I(1:nIonFluid))    
       totalLossRho=sum(LiSpecies_I(1:nIonFluid))    
       !sum of the (Loss term) of all ion species
       totalLossNumRho = sum(LiSpecies_I/MassIon_I)
       LossNumRho_I = LiSpecies_I/MassIon_I
       !sum of the (loss term/atom mass) of all ..
       totalSourceNumRho = sum(SiSpecies_I/MassIon_I)
       SourceNumRho_I    = SiSpecies_I/MassIon_I
       ! sum of the (Source term/atom mass) of all..
       totalLossx=totalLossRho*inv_rho
       Lossx_I = LiSpecies_I*invRho_I
       totalLossNumx = totalLossNumRho/NumDens
       LossNumx_I = LossNumRho_I/NumDens_I
       totalPSNumRho = sum(PhoIon_I)
       ! sum of the photonionziation source/atom mass) of all..
       totalRLNumRhox = sum(Recb_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
       RLNumRhox_I    =     Recb_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I
       !sum of the (loss term/atom mass) due to recombination    

       MaxSLSpecies_CB(i,j,k,iBlock)=maxval(abs(SiSpecies_I(1:nIonFluid)+&
            LiSpecies_I(1:nIonFluid) ) /&
            (State_VGB(iRhoIon_I(1:nIonFluid), i,j,k, iBlock)+1e-20))&
            /vInv_CB(i,j,k,iBlock)

       if(.not.UsePointImplicit_B(iBlock) )then
           !sum of the (loss term/atom mass) due to recombination                                                  \
                                                                                                                    
           SourceLossMax = 10.0*maxval(abs(SiSpecies_I(1:nSpecies)-&
                LiSpecies_I(1:nSpecies) ) /&
                (State_VGB(iRhoIon_I(1:nIonFluid), i,j,k, iBlock)+1e-20))&
                /vInv_CB(i,j,k,iBlock)
           vdtmin=min(VdtFace_x(i,j,k),VdtFace_y(i,j,k),VdtFace_z(i,j,k))



            if(SourceLossMax > Vdtmin) then
              !UsePointImplicit_B(iBlock)=.true.                                                                   \
              VdtFace_x(i,j,k) = max (SourceLossMax, VdtFace_x(i,j,k) )
              VdtFace_y(i,j,k) = max (SourceLossMax, VdtFace_y(i,j,k) )
              VdtFace_z(i,j,k) = max (SourceLossMax, VdtFace_z(i,j,k) )
           end if

        end if


       Source_VC(iRhoIon_I,i,j,k) =Source_VC(iRhoIon_I,i,j,k) &
            +SiSpecies_I &
            -LiSpecies_I


       Source_VC(rho_     ,i,j,k)=Source_VC(rho_     ,i,j,k)&
            +sum(SiSpecies_I(1:nIonFluid))&
            -sum(LiSpecies_I(1:nIonFluid))

       Source_VC(rhoUx_     ,i,j,k)= Source_VC(rhoUx_     ,i,j,k) &
            -State_VGB(Ux_,i,j,k,iBlock)*totalLossx&
            -nu_BLK(i,j,k,iBlock)*State_VGB(Ux_,i,j,k,iBlock)

       Source_VC(iRhoUxIon_I,i,j,k)=Source_VC(iRhoUxIon_I,i,j,k)&
            -State_VGB(iRhoUxIon_I,i,j,k,iBlock)*Lossx_I&
            -nu_BLK(i,j,k,iBlock)*State_VGB(iRhoUxIon_I,i,j,k,iBlock) 

       Source_VC(rhoUy_     ,i,j,k) = Source_VC(rhoUy_     ,i,j,k)  &
            -State_VGB(Uy_,i,j,k,iBlock)*totalLossx&
            -nu_BLK(i,j,k,iBlock)*State_VGB(Uy_,i,j,k,iBlock)

       Source_VC(iRhoUyIon_I,i,j,k)= Source_VC(iRhoUyIon_I,i,j,k)&
            -State_VGB(iRhoUyIon_I,i,j,k,iBlock)*Lossx_I&
            -nu_BLK(i,j,k,iBlock)*State_VGB(iRhoUyIon_I,i,j,k,iBlock)

       Source_VC(rhoUz_     ,i,j,k)= Source_VC(rhoUz_     ,i,j,k)  &
            -State_VGB(Uz_,i,j,k,iBlock)*totalLossx&
            -nu_BLK(i,j,k,iBlock)*State_VGB(Uz_,i,j,k,iBlock)

       Source_VC(iRhoUzIon_I,i,j,k)= Source_VC(iRhoUzIon_I,i,j,k)&
            -State_VGB(iRhoUzIon_I,i,j,k,iBlock)*Lossx_I&
            -nu_BLK(i,j,k,iBlock)*State_VGB(iRhoUzIon_I,i,j,k,iBlock)
       !----- pressure and energy source terms

      
          tmp = totalSourceNumRho*TNu_body            &
               + NumDens*(TNu_body-KTi)*nu_BLK(i,j,k,iBlock) &
               + totalPSNumRho*kTe0                &
               - totalLossNumRho*kTi               &
               - totalRLNumRhox*NumDens*KTe        
          tmp_I = SourceNumRho_I*TNu_body            &
               + NumDens_I*(TNu_body-Temp_I)*nu_BLK(i,j,k,iBlock) &
               + PhoIon_I*kTe0                &
               - LossNumRho_I*Temp_I               &
               - RLNumRhox_I*NumDens_I*KTe         

          if(DoTestCell)then
             write(*,*)NameSub,'SourceNumRho_I=',SourceNumRho_I
             write(*,*)NameSub,'NumDens_I=',NumDens_I
             write(*,*)NameSub,'PhoIon_I=',PhoIon_I
             write(*,*)NameSub,'LossNumRho_I=',LossNumRho_I
             write(*,*)NameSub,'RLNumRhox_I=',RLNumRhox_I
          end if


          do iFluid=1,nIonFluid
             col_ii_I(iFluid)=sum(2*State_VGB(iRhoIon_I(iFluid),i,j,k,iBlock)*State_VGB(iRhoIon_I,i,j,k,iBlock)&
                  /(MassIon_I)*(Temp_I-Temp_I(iFluid))*&
               IonIonCoeff_II(iFluid,:)*ReducedMassIon_II(:,iFluid)&
               /(Temp_I*No2Si_V(UnitTemperature_))**(3/2))

          end do



          tmp_I=tmp_I+ col_ii_I

          tmp= tmp+sum(col_ii_I)


          Source_VC(Energy_ ,i,j,k) =  Source_VC(Energy_ ,i,j,k)&
               -0.5*State_VGB(rho_ ,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock)-0.50*uu2*(totalLossRho) +inv_gm1*tmp

          do iFluid=1,nIonFluid
             Source_VC(Energy_ + iFluid,i,j,k) =  Source_VC(Energy_ + iFluid,i,j,k)&
                  -0.5*State_VGB(iRhoIon_I(iFluid),i,j,k,iBlock)*uu2_I(iFluid)*&
                  nu_BLK(i,j,k,iBlock)-0.50*uu2_I(iFluid)*(LiSpecies_I(iFluid)) &
                  +inv_gm1*tmp_I(iFluid)
          end do


          Source_VC(P_     ,i,j,k)  = Source_VC(P_     ,i,j,k)   &
               +0.5*gm1*State_VGB(rho_ ,i,j,k,iBlock)*uu2*&
               nu_BLK(i,j,k,iBlock)  &
               +0.50*(gm1)*uu2*(totalSourceRho) &
               +tmp



          Source_VC(iPIon_I,i,j,k)  = Source_VC(iPIon_I,i,j,k)   &
               + 0.5*gm1*State_VGB(iRhoIon_I,i,j,k,iBlock)*uu2_I*&
               nu_BLK(i,j,k,iBlock)  &
               + 0.5*gm1*uu2_I*SiSpecies_I &
               + tmp_I
       
    end do; end do; end do     ! end of the i,j,k loop

    !end of chemistry

  end subroutine user_calc_sources
  !=============================================================================
  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !------------------------------------------------------------------------

    ! All ion momenta are implicit
    allocate(iVarPointImpl_I(5*nIonFluid))

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
  ! ===========================================================================

  subroutine user_init_session
    use ModMain
    use ModPhysics
    use ModVarIndexes

    integer::iBoundary
    !--------------------------------------------------------------------------
    call set_multiSp_ICs  
    !    Rbody = 1.0 + 140.0e3/Venus
    BodyRho_I(1) = sum(BodyRhoSpecies_I(1:MaxSpecies))
    !BodyRho_I(2:nFluid)=BodyRhoSpecies_I(1:MaxSpecies)
    BodyP_I(1)   =sum(BodyRhoSpecies_I(1:MaxSpecies)&
         /MassIon_I)*Ti_body*(1+ElectronPressureRatio)
   
    BodyP_I(2:nFluid)=Ti_body*(BodyRhoSpecies_I(1:MaxSpecies)&
         /MassIon_I)

    FaceState_VI(rho_,body1_)=BodyRho_I(1)
    FaceState_VI(iRhoIon_I,body1_) = BodyRhoSpecies_I
    FaceState_VI(P_,body1_)=BodyP_I(1)

    CellState_VI(:,body1_:Top_)=FaceState_VI(:,body1_:Top_)
    do iBoundary=body1_,Top_  
       CellState_VI(rhoUx_:rhoUz_,iBoundary) = &
            FaceState_VI(Ux_:Uz_,iBoundary)*FaceState_VI(rho_,iBoundary)       
    end do

    if(.not.allocated(nDenNuSpecies_CBI))then
       allocate(nDenNuSpecies_CBI(nI, nJ, nK, nBLK, MaxNuSpecies))
       nDenNuSpecies_CBI = -1.0

       allocate( TempNuSpecies_CBI(1:nI, 1:nJ, 1:nK, nBLK))
       allocate(Productrate_CB(1:nI, 1:nJ, 1:nK, nBLK))
       allocate(Ionizationrate_CBI(1:nI, 1:nJ, 1:nK, nBLK,2))
       allocate(MaxSiSpecies_CB(1:nI, 1:nJ, 1:nK, nBLK))
       allocate(MaxLiSpecies_CB(1:nI, 1:nJ, 1:nK, nBLK))
       allocate(MaxSLSpecies_CB(1:nI, 1:nJ, 1:nK, nBLK))
       allocate(nu_BLK(1:nI,1:nJ,1:nK,nBLK))
       allocate(nu1_BLK(1:nI,1:nJ,1:nK,nBLK))
    end if

  end subroutine user_init_session
  !============================================================================
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

       case("#SOLARCON") !solar cycle condition
          call read_var('SolarCon',SolarCond)

       !case("#UseHotO")  !adding hot Oxygen or not
          !call read_var('UseHotO',UseHotO)
          !call read_var('CoeffHotO',CoeffHotO)

    !   case("#UseTempCont") !add hoc term of the energy source
    !      call read_var('UseTempCont',UseTempCont)          

       case('#REACTIONS')
          call read_var('UseImpactIon',UseImpactIon)
          call read_var('UseChargeEx',UseChargeEx)
          ! open(15,file='read_in.dat')
          ! do i=1,32
          !    read(15,*)Temp_dim(i),Impact_ION_dim(i,Op_),Impact_ION_dim(i,Hp_)
          ! end do
          ! close(15)


          case('#IonNeuCollision')
          call read_var('nu0_dim',nu0_dim)


          !in ModUserVenus single fluid
      !    case('#REACTION')
      !    call read_var('number',number)
      !    Ratedim_I(number)=0.0  

      !case('#INNERBCS')
      !    call read_var('type_innerbcs',type_innerbcs)

          case('#USECHAPMAN')
          call read_var('UseChapman',UseChapman)

       case('#POINTIMPLICITREGION')
          call read_var('rPointImplicit',rPointImplicit)
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  !============================================================================
  subroutine set_multiSp_ICs
    use ModMain
    use ModConst
    use ModIO
    use ModPhysics

    real :: Productrate
    logical::oktest=.false., oktestme=.false.
    real::alt0
    !---------------------------------------------------------------
  
  if(oktestme)then
       write(*,*)'in set_multisp_ICs, No2Io_V(UnitN_),t=',&
            No2Io_V(UnitN_),No2Io_V(UnitT_)
       write(*,*)'No2Si_V(UnitX_), temperature=',&
            No2Si_V(UnitX_), No2Si_V(UnitTemperature_)
       write(*,*)'BodynDenNuSpecies_dim_I(:)',&
            BodynDenNuSpdim_I(:)
    end if
    select case(SolarCond)

    case('solarmax')

       
      ! Tnu_body_dim = 134.0      ! neutral temperature
       Tnu_body_dim =300

       RateDim_I(CO2_hv__CO2p_em_)=1.695e-6/0.723/0.723/2.5 !scale to Venus                                                                    
       RateDim_I(O_hv__Op_em_) = 6.346e-7/0.723/0.723
       BodynDenNuSpDim_I(CO2_)= 1.0e15
       BodynDenNuSpDim_I(O_  )= 2.0e11
       BodynDenNuSpDim_I(Oh_)= 6.3119e4
       BodynDenNuSpDim_I(Ohx_)= 3.9646e3
       HNuSpeciesDim_I(CO2_)=5.5e3
       HNuSpeciesDim_I(O_  )=17.0e3
       HNuSpeciesDim_I(Oh_)=290.5e3
       HNuSpeciesDim_I(Ohx_)=2436.6e3

        RateDim_I(CO2_hv__CO2p_em_)=1.30e-6
       RateDim_I(O_hv__Op_em_) = 1.21e-6

case('solarmin')

!   Tnu_body_dim = 117.0      ! neutral temperature

   Tnu_body_dim =300

       RateDim_I(CO2_hv__CO2p_em_)=6.696e-7/0.723/0.723/3.0 !scale to Venus                                                                                            
       RateDim_I(O_hv__Op_em_) = 2.44e-7/0.723/0.723
       BodynDenNuSpDim_I(CO2_)=1.0e15 !neutral density in cm^-3                                                                                                        
       BodynDenNuSpDim_I(O_  )=1.3e11
       HNuSpeciesDim_I(CO2_)=5.07e3  !scale height in meter                                                                                                            
       HNuSpeciesDim_I(O_  )=15.5e3
       
       RateDim_I(CO2_hv__CO2p_em_)=4.27e-7
       RateDim_I(O_hv__Op_em_) = 4.67e-7



    case default
       call stop_mpi('unknow solar condition='//SolarCond)
    end select



    T300 = T300_dim*Si2No_V(UnitTemperature_)

    Ti_body_dim = Tnu_body_dim  !ion temperature at the body                                                                              


    TNu_body= TNu_body_dim*Si2No_V(UnitTemperature_)
    Ti_body = Ti_body_dim*Si2No_V(UnitTemperature_)
    T300 = T300_dim*Si2No_V(UnitTemperature_)

    kTe0=max(Te_new_dim, Tnu_body_dim)*Si2No_V(UnitTemperature_)



    nu0=nu0_dim*No2Io_V(UnitN_)*No2Io_V(UnitT_)




    BodynDenNuSpecies_I(1:nNuSpecies)=&
         BodynDenNuSpDim_I(1:nNuSpecies)/No2Io_V(UnitN_)
 
    HNuSpecies_I(1:nNuSpecies)=&
         HNuSpeciesDim_I(1:nNuSpecies)*Si2No_V(UnitX_)

   
    alt0= Rbody-Altitude0*Si2No_V(UnitX_)-1.0

    !write(*,*)'!!!Altitude0*Si2No_V(UnitX_=',Altitude0*Si2No_V(UnitX_)
    !write(*,*)'Rbody -1=',Rbody-1
    !write(*,*)'alt0=',alt0

   ! BodynDenNuSpecies_I(1:nNuSpecies)=&
   !      BodynDenNuSpecies_I(1:nNuSpecies)* &
   !      exp(-alt0/HNuSpecies_I(1:nNuSpecies))





    ! normlize the reaction rate                                                                                                                                       
    Rate_I(CO2_hv__CO2p_em_)= &
         Ratedim_I(CO2_hv__CO2p_em_)*No2Io_V(UnitT_)
    Rate_I(O_hv__Op_em_)=  &
         Ratedim_I(O_hv__Op_em_)*No2Io_V(UnitT_)
    Rate_I(CO2p_O__O2p_CO_)=  &
         Ratedim_I(CO2p_O__O2p_CO_)  &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(Op_CO2__O2p_CO_)=  &
         Ratedim_I(Op_CO2__O2p_CO_)*exp(log((8.0/3.0)*T300/TNu_body)*0.39) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_O__Op_CO2_)=  &
         Ratedim_I(CO2p_O__Op_CO2_) &
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)

    Rate_I(O2p_em__O_O_)=  &
         Ratedim_I(O2p_em__O_O_)*exp(log(4.0*T300/TNu_body)*0.56)&
         *No2Io_V(UnitT_)*No2Io_V(UnitN_)
    Rate_I(CO2p_em__CO_O_)=  &
         Ratedim_I(CO2p_em__CO_O_)*sqrt(T300/TNu_body)&
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

    write(*,*)'!!!! 1 CrossSectiondim_I=', CrossSection_I
    

     CrossSection_I = CrossSectiondim_I*1.0e-4 / &
         ( Si2No_V(UnitX_)*Si2No_V(UnitN_) )

     write(*,*)'!!!! 2 CrossSectiondim_I=', CrossSection_I


    IonNeuCoeff_II= IonNeuCoeffDim_II*No2Io_V(UnitN_)*No2Io_V(UnitT_)*1.e-10

    IonIonCoeff_II=IonIonCoeffDim_II*No2Io_V(UnitN_)*No2Io_V(UnitT_)


    Optdep =  sum(BodynDenNuSpecies_I*CrossSection_I*HNuSpecies_I)
    Productrate0 = max(exp(-Optdep), 1.0e-30)
    Productrate = Productrate0

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
    BodyRhoSpecies_I=BodyRhoSpecies_I*&
         MassIon_I

    do iFluid=1,nIonFluid
       ReducedMassNeutral2_II(iFluid,:)=1/&
            (MassIon_I(iFluid) + MassNeutral_I)
       ReducedMassNeutral_II(iFluid,:)=MassNeutral_I/&
            (MassIon_I(iFluid) + MassNeutral_I)
       ReducedMassIon_II(iFluid,:)=1/&
            (MassIon_I(iFluid) + MassIon_I)
    end do



!if(DoTestMe)then
  !     write(*,*)'crosssection= ',CrossSection_I, 'optdep=', Optdep
  !     write(*,*)'producationrate0=',Productrate0
  !     write(*,*)'hnuspecies_I=',HNuSpecies_I(1:nNuSpecies)
       write(*,*)' set parameters of Venus: BodyRhoSpecies_I(i)=',&
            BodyRhoSpecies_I(1:nSpecies)
       write(*,*)'BodynDenNuSpecies_I=',BodynDenNuSpecies_I
       write(*,*)'Rate_I=',Rate_I
  !     write(*,*)'neutral density=', &
  !          BodynDenNuSpecies_I(:)
  !     write(*,*)'nu0=',nu0
  !     write(*,*)'Rate_I=', Rate_I
  !     write(*,*)'Rate_dim_I=', Ratedim_I
!       call stop_mpi('end')                                                                                                                                           
 !   end if

  end subroutine set_multiSp_ICs
  !============================================================================

  ! This subroutine allows the user to apply initial conditions to the domain
  ! which are problem specific and cannot be created using the predefined
  ! options in BATSRUS.
  ! The variables specific to the problem are loaded from ModUser

  subroutine user_set_ICs

    use ModProcMH, ONLY : iProc
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
    use ModIO, ONLY : restart
    use ModPhysics
    use ModNumConst

    real :: Rmax, SinSlope, CosSlope, tempo
    real :: B4, dB4dx, zeta4, q4, epsi4, plobe, &
         XFace, YFace, ZFace
    real :: temp1,temp2,temp3
    integer :: i,j,k,q
    real, dimension(nIonFluid)::Temp_I
    integer:: iBoundary
   real :: CosSZA,Xp, chap_y, chap, sinSZA, iBlock
    character (len=*), parameter :: NameSub = 'user_set_ics'
    logical:: DoTest, DoTestMe, DoTestCell
    !-------------------------------------------------------------------------
    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    call set_neutral_density(globalBlk)

    if(DoTestMe)then
       write(*,*)'in set_ics'
       write(*,*)'BodynDenNuSpecies_I=',&
            BodynDenNuSpecies_I
       WRITE(*,*)''
       write(*,*)'HNuSpecies_I(1:nNuSpecies)=',HNuSpecies_I
       WRITE(*,*)''
       write(*,*)'Rbody=', Rbody
       write(*,*)''
    end if


    ! calculate optical depth and producation rate                                                                                                             
  !  do k=1,nK; do j=1,nJ; do i=1,nI
      ! cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
      !      x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)&
      !      +5.0e-4

     !  Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:MaxNuSpecies)*&
     !       CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
     !       6.0e-3)/cosSZA
     !  if( Optdep<11.5 .and. x_BLK(i,j,k,globalBLK) > 0.0) then
     !     Productrate_CB(i,j,k,globalBLK) = max(exp(-Optdep), 1.0e-5)
     !  else
     !     Productrate_CB(i,j,k,globalBLK) = 1.0e-5
     !  end if

!    end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          iBlock=globalBLK
           !  cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
           !       x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-4)&
           !       +1.0e-4
           !  Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:MaxNuSpecies)*&
           !       CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
           !       6.0e-3)/cosSZA
           !  if( Optdep<11.5 .and. x_BLK(i,j,k,globalBLK) > 0.0) then
           !     Productrate_CB(i,j,k,globalBLK) = max(exp(-Optdep), 1.0e-30)
           !  else
           !     Productrate_CB(i,j,k,globalBLK) = 1.0e-30
           !  end if

           Optdep=HNuSpecies_I(CO2_)*nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)*&
            CrossSection_I(CO2_)
       cosSZA = X_BLK(i,j,k,iBlock) /R_BLK(i,j,k,iBlock)
       if(UseChapman)then            
          if(Optdep > 13.8) then 
             Productrate_CB(i,j,k,iBlock) = 1.0e-6*&
                  max(cosSZA, 1.0e-6)                
          else
             Xp=R_BLK(i,j,k,iBlock)/HNuSpecies_I(1)
             chap_y= sqrt(0.5*Xp)*abs(cosSZA)
             if(cosSZA .gt.0.0)then   !SZA<90 deg (equation 13)
                if(chap_y<8.0)then
                   chap = sqrt(3.1415926/2.0*Xp)*&
                        (1.0606963+0.5564383*chap_y)/&
                        (1.0619896+1.7245609*chap_y+chap_y*chap_y)
                elseif(chap_y<100.0)then
                   chap = sqrt(3.1415926/2.0*Xp)*&
                        0.56498823/(0.6651874+chap_y)
                else
                   chap=0.0
                end if
             elseif(cosSZA > -0.5) then
                ! 120> SZA > 90 deg (equation 15) Smith and Smith, 1972

                sinSZA = sqrt(1.0 - cosSZA**2)
                if(chap_y<8.0)then
                   chap =sqrt(2.0*3.1415926*Xp)* &
                        (sqrt(sinSZA)*exp(Xp*(1.0-sinSZA)) &
                        -0.5*(1.0606963+0.5564383*chap_y)/ &
                        (1.0619896+1.7245609*chap_y+chap_y*chap_y))
                elseif(chap_y<100.0)then
                   chap =sqrt(2.0*3.1415926*Xp)* &
                        (sqrt(sinSZA)*exp(Xp*(1.0-sinSZA)) &
                        -0.5*0.56498823/(0.6651874+chap_y))
                else
                   chap =0.0
                end if
             else
                chap =1.0e10
             end if
             
             Optdep=Optdep*chap      
             Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-6)
          end if
       else  !old way of optical depth calculation, 
          !production rate drops to 1.0e-5 for SZA > 90
          
          cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,iBlock)))*&
               x_BLK(i,j,k,iBlock)/max(R_BLK(i,j,k,iBlock),1.0e-5)&
               +1.0e-5
          Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,iBlock,1:MaxNuSpecies)*&
               CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
               6.0e-3)/cosSZA
          if( Optdep<11.5 .and. x_BLK(i,j,k,iBlock) > 0.0) then 
             Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-5)
          else
             Productrate_CB(i,j,k,iBlock) = 1.0e-30
          end if
       end if
          end do; end do; end do

    nu1_BLK(:,:,:,globalBLK)=nu_BLK(:,:,:,globalBLK)
    do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
       if (R_BLK(i,j,k,globalBLK)< Rbody) then
          cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
               x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
               1.0e-3
          State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
          State_VGB(OpRho_,i,j,k,globalBLK)= &
               CellState_VI(OpRho_,body1_)*cosSZA
          State_VGB(O2pRho_,i,j,k,globalBLK)= &
               CellState_VI(OpRho_,body1_)*sqrt(cosSZA)
          State_VGB(CO2pRho_,i,j,k,globalBLK)= &
               CellState_VI(OpRho_,body1_)*cosSZA
          State_VGB(rho_,i,j,k,globalBLK)  = &
               sum( State_VGB(iRho_I(2:nFluid),i,j,k,globalBLK))

          State_VGB(iPIon_I,i,j,k,globalBLK) = Ti_body*State_VGB(iRhoIon_I,i,j,k,globalBLK) &
               /MassIon_I 
          State_VGB(P_,i,j,k,globalBLK) = &
               max(SW_p,sum(State_VGB(iPIon_I,i,j,k,globalBLK))*(1+ElectronPressureRatio))
       else

          State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)

          !write(*,*)'CellState_VI(HpRho_,1_)=',CellState_VI(HpRho_,1)
          !write(*,*)'State_VGB(HpRho_,i,j,k,globalBLK)=',State_VGB(HpRho_,i,j,k,globalBLK)
!!$           write(*,*)'State_VGB(rho_,i,j,k,globalBLK)=',State_VGB(rho_,i,j,k,globalBLK)
          State_VGB(Ux_:bz_,i,j,k,globalBLK)   =0.0
       end if
       

    end do;end do; end do;

    !if(DoTestMe)then
       !write(*,*)'!!!!in set_ics 0 '
 !      write(*,*)'State_VGB(iRhoIon_I)=',State_VGB(iRhoIon_I,iTest,jTest,kTest,BLKtest)
       !write(*,*)'nDenNuSpecies_CBI(i,j,k,globalBLK,:)=',nDenNuSpecies_CBI(iTest,jTest,kTest,globalBLK,:)
       !write(*,*)'Rate_I=',Rate_I
       !write(*,*)'Ionizationrate_CBI(i,j,k,globalBLK,CO2_)=',Ionizationrate_CBI(iTest,jTest,kTest,globalBLK,:)
     !         write(*,*)'State_VGB(iRhoIon_I)=',State_VGB(iRhoIon_I,iTest,jTest,kTest,BLKtest)
      ! write(*,*)'State_VGB(iPIon_I(1:nIonFluid),i,j,k,globalBLK)==',State_VGB(iPIon_I(1:nIonFluid),iTest,jTest,kTest,BLKtest)
      ! write(*,*)'sum(State_VGB(iRhoIon_I))=',sum(State_VGB(iRhoIon_I,iTest,jTest,kTest,BLKtest))
      ! write(*,*)'State_VGB(Rho_)=',State_VGB(Rho_,iTest,jTest,kTest,BLKtest)
   ! end if


!!$    if(DoTestMe)&
!!$         write(*,*)'state_VGB(body1_)=',&
!!$         CellState_VI(:,body1_),'cell_state_VI(:,1)=',CellState_VI(:,1)

    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(iRhoUx_I,i,j,k,globalBLK) = 0.0
       State_VGB(iRhoUy_I,i,j,k,globalBLK) = 0.0
       State_VGB(iRhoUz_I,i,j,k,globalBLK) = 0.0



       if (.not. (true_cell(i,j,k,globalBLK).and. &
            R_BLK(i,j,k,globalBLK)<1.5*Rbody) ) CYCLE


       cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
            x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
            1.0e-3

       State_VGB(CO2pRho_,i,j,k,globalBLK)= &
            Ionizationrate_CBI(i,j,k,globalBLK,CO2_) &
            /nDenNuSpecies_CBI(i,j,k,globalBLK,O_)   &
            /(Rate_I(CO2p_O__O2p_CO_)+Rate_I(CO2p_O__Op_CO2_))

       !write(*,*)'State_VGB(CO2pRho_,i,j,k,globalBLK)',State_VGB(CO2pRho_,i,j,k,globalBLK)

       State_VGB(OpRho_,i,j,k,globalBLK)= &
            (Ionizationrate_CBI(i,j,k,globalBLK,O_) &
            +Rate_I(CO2p_O__Op_CO2_)                &
            *State_VGB(CO2pRho_,i,j,k,globalBLK)    &
            *nDenNuSpecies_CBI(i,j,k,globalBLK,O_)) &
            /(nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)+4.0e6)&
            /Rate_I(Op_CO2__O2p_CO_)

       !write(*,*)'State_VGB(OpRho_,i,j,k,globalBLK)',State_VGB(OpRho_,i,j,k,globalBLK)

       State_VGB(O2pRho_,i,j,k,globalBLK)= &
            SQRT((nDenNuSpecies_CBI(i,j,k,globalBLK,O_)*&
            State_VGB(CO2pRho_,i,j,k,globalBLK)*&
            Rate_I(CO2p_O__O2p_CO_)+&
            nDenNuSpecies_CBI(i,j,k,globalBLK, CO2_)*&
            State_VGB(OpRho_,i,j,k,globalBLK)*&
            Rate_I(Op_CO2__O2p_CO_)+1e-10)/Rate_I(O2p_em__O_O_))

       !write(*,*)'State_VGB(O2pRho_,i,j,k,globalBLK)',State_VGB(O2pRho_,i,j,k,globalBLK)

       ! Convert to mass densities
       State_VGB(iRho_I(2:nFluid),i,j,k,globalBLK)=&
            State_VGB(iRho_I(2:nFluid),i,j,k,globalBLK)*&
            MassIon_I

    end do; end do; end do

    !


    do k=1,nK; do j=1,nJ; do i=1,nI

       if(.not.true_cell(i,j,k,globalBLK))CYCLE

       !IC for velocity
       State_VGB(iRhoUx_I,i,j,k,globalBLK) = 0.0
       State_VGB(iRhoUy_I,i,j,k,globalBLK) = 0.0
       State_VGB(iRhoUz_I,i,j,k,globalBLK) = 0.0


       State_VGB(rho_,i,j,k,globalBLK)   =&
            sum(State_VGB(iRhoIon_I,i,j,k,globalBLK))

       do q=1,nSpecies
          !write(*,*)'I got to the low density ratio'
          if(State_VGB(iRhoIon_I(q),i,j,k,globalBLK) < &
               LowDensityRatio* State_VGB(Rho_,i,j,k,globalBLK))then
             State_VGB(iRhoIon_I(q),i,j,k,globalBLK)= LowDensityRatio*&
                  State_VGB(Rho_,i,j,k,globalBLK)
          end if
       end do

       State_VGB(rho_,i,j,k,globalBLK)   =&
            sum(State_VGB(iRhoIon_I,i,j,k,globalBLK))

       State_VGB(P_,i,j,k,globalBLK) = &
            max(SW_p,(sum(State_VGB(iRhoIon_I,i,j,k,globalBLK)/(MassIon_I)))*Ti_body*(1+ElectronPressureRatio))

       State_VGB(iPIon_I,i,j,k,globalBLK)=State_VGB(P_,i,j,k,globalBLK)&
            /(sum(State_VGB(iRhoIon_I,i,j,k,globalBLK)/(MassIon_I)))*&
            State_VGB(iRhoIon_I,i,j,k,globalBLK)/MassIon_I/(1+ElectronPressureRatio)

       Temp_I=State_VGB(iPIon_I,i,j,k,globalBLK)/&
            (State_VGB(iRhoIon_I,i,j,k,globalBLK)/MassIon_I)

    end do; end do; end do

  end subroutine user_set_ICs

  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,       ONLY: nDim,West_,North_,Top_	
    use ModMain,       ONLY: UseRotatingBc, iTest, jTest, kTest, ProcTest, &
         BlkTest, GLOBALBLK, ExtraBc_, Body1_
    use ModProcMH,   ONLY: iProc
    use ModVarIndexes, ONLY: nVar, OpRho_, O2pRho_, CO2pRho_, HpRho_,HpP_,O2pP_,OpP_,CO2pP_,iRhoUx_I,iRhoUy_I,iRhoUz_I
    use ModPhysics,    ONLY: SW_rho, SW_p, SW_T_dim,ElectronPressureRatio,FaceState_VI
    use ModFaceBc,     ONLY: FaceCoords_D, VarsTrueFace_V, iFace,jFace,kFace, iBoundary
    

    real, intent(out):: VarsGhostFace_V(nVar)

    real:: XFace,YFace,ZFace,rFace,rFace2
    real:: v_phi(3)
    real:: cosSZA 
    real:: uDotR_I(nFluid), bDotR
    integer:: i,j,k
    character (len=*), parameter :: NameSub = 'user_face_bcs'
    logical:: DoTest, DoTestMe, DoTestCell
    !-------------------------------------------------------------------------

    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if(iBoundary == ExtraBc_)then
       VarsGhostFace_V = FaceState_VI(:,east_)
       RETURN
    elseif(iBoundary /= Body1_)then       
       call stop_mpi(NameSub//' invalid iBoundary value')
    end if

    DoTestCell = DoTestMe &
         .and. iFace==iTest .and. jFace==jTest .and. kFace==kTest

    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)

    rFace2 = XFace**2 + YFace**2 + ZFace**2
    rFace  = sqrt(rFace2)

    !Apply boundary conditions
    cosSZA=(0.5+sign(0.5,XFace)) * XFace/max(RFace,1.0e-3) + 1.0e-3

    VarsGhostFace_V(OpRho_)  = BodyRhoSpecies_I(Op_) *cosSZA

    VarsGhostFace_V(O2pRho_) = BodyRhoSpecies_I(O2p_)*sqrt(cosSZA)

    VarsGhostFace_V(CO2pRho_)= BodyRhoSpecies_I(CO2p_)*cosSZA

    VarsGhostFace_V(HpRho_)  = SW_rho*0.3
    VarsGhostFace_V(rho_) = sum(VarsGhostFace_V(iRhoIon_I))    

    VarsGhostFace_V(iPIon_I)=Ti_body*VarsGhostFace_V(iRhoIon_I)/MassIon_I

    VarsGhostFace_V(P_)=sum(VarsGhostFace_V(iPIon_I))*(1+ElectronPressureRatio)


    ! Reflective in radial direction
    uDotR_I = (VarsTrueFace_V(iRhoUx_I)*FaceCoords_D(1)+ &
         VarsTrueFace_V(iRhoUy_I)*FaceCoords_D(2)+ &
         VarsTrueFace_V(iRhoUz_I)*FaceCoords_D(3))/ rFace2

    ! bDotR = sum(VarsTrueFace_V(Bx_:Bz_)*FaceCoords_D)/rFace2

    VarsGhostFace_V(iRhoUx_I) = VarsTrueFace_V(iRhoUx_I) - 2*uDotR_I*FaceCoords_D(1)
    VarsGhostFace_V(iRhoUy_I) = VarsTrueFace_V(iRhoUy_I) - 2*uDotR_I*FaceCoords_D(2)
    VarsGhostFace_V(iRhoUz_I) = VarsTrueFace_V(iRhoUz_I) - 2*uDotR_I*FaceCoords_D(3)

    ! VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_) - 2*bDotR*FaceCoords_D
    VarsGhostFace_V(Bx_:Bz_) = 0.0
   
    !let ues try float like Yingjuan
    ! VarsGhostFace_V(iRhoUx_I) = VarsTrueFace_V(iRhoUx_I) !-2*uDotR_I*FaceCoords_D(1)
    !VarsGhostFace_V(iRhoUy_I) = VarsTrueFace_V(iRhoUy_I) !- 2*uDotR_I*FaceCoords_D(2)
    !VarsGhostFace_V(iRhoUz_I) = VarsTrueFace_V(iRhoUz_I) !- 2*uDotR_I*FaceCoords_D(3)
    !VarsGhostFace_V(Bx_:Bz_) = &
    !VarsTrueFace_V(Bx_:Bz_)

    
    if(DoTestCell) then
       write(*,*)'VarsGhostFace_V(iRhoIon_I)=',VarsGhostFace_V(iRhoIon_I)
       write(*,*)'VarsGhostFace_V(iRhoUx_I(1))=',VarsGhostFace_V(iRhoUx_I(1))
       write(*,*)'VarsGhostFace_V(iRhoUy_I(1))=',VarsGhostFace_V(iRhoUy_I(1))
       write(*,*)'VarsGhostFace_V(iRhoUz_I(1))=',VarsGhostFace_V(iRhoUz_I(1))
       write(*,*)'VarsGhostFace_V(P_)=',VarsGhostFace_V(P_)
       write(*,*)'VarsGhostFace_V(iPIon_I)=',VarsGhostFace_V(iPIon_I)
       write(*,*)'VarsGhostFace_V(Bx_ : Bz_)=',VarsGhostFace_V(Bx_: Bz_)
    end if

    ! Apply corotation?
    if (UseRotatingBc) then
       call calc_corotation_velocities(FaceCoords_D, v_phi)
       VarsGhostFace_V(iRhoUx_I) = VarsGhostFace_V(iRhoUx_I) + 2*v_phi(1)
       VarsGhostFace_V(iRhoUy_I) = VarsGhostFace_V(iRhoUy_I) + 2*v_phi(2)
       VarsGhostFace_V(iRhoUz_I) = VarsGhostFace_V(iRhoUz_I) + 2*v_phi(3)
    end if

  end subroutine user_face_bcs

  !========================================================================

  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, x_Blk, x2
    use ModBoundaryCells, ONLY: SaveBoundaryCells

    implicit none

    integer, intent(in):: iBlock

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = x_Blk(:,:,:,iBlock) > x2

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

  !========================================================================
 subroutine user_get_log_var(VarValue, TypeVar, Radius)
     use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK,&
          dx_BLK,dy_BLK,dz_BLK
     use ModMain
     use ModVarIndexes
     use ModAdvance,    ONLY: State_VGB,tmp1_BLK
     use ModPhysics,ONLY: No2Si_V, UnitN_, UnitX_, UnitU_

     real, intent(out)            :: VarValue
     character (len=*), intent(in):: TypeVar
     real, intent(in), optional :: Radius

     real, external :: calc_sphere
     real ::mass,value1,value2
     integer:: i,j,k,iBLK
     character (len=*), parameter :: Name='user_get_log_var'
     logical:: oktest=.false.,oktest_me

!---------------------------------------------------------------------------
     write(*,*)'user_get_log_var called'
     call set_oktest('user_get_log_var',oktest,oktest_me)
     if(oktest)write(*,*)'in user_get_log_var: TypeVar=',TypeVar

     select case(TypeVar)
        
case('hplflx')
       mass=1.0
       !write(*,*)'hola 2'
       do iBLK=1,nBLK
         ! write(*,*)'iBlock=',iBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             !write(*,*)'i,j,k=',i,j,k
             tmp1_BLK(i,j,k,iBLK) = State_VGB(HpRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)>0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


      case('hprflx')
       mass=1.0
       !write(*,*)'hola 2'
       do iBLK=1,nBLK
         ! write(*,*)'iBlock=',iBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             !write(*,*)'i,j,k=',i,j,k
             tmp1_BLK(i,j,k,iBLK) = State_VGB(HpRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


       case('hpflx')
       mass=1.0
       do iBLK=1,nBLK
         ! write(*,*)'iBlock=',iBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             !write(*,*)'i,j,k=',i,j,k
             tmp1_BLK(i,j,k,iBLK) =&
                  (State_VGB(HpRhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(HpRhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(HpRhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)
             !if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass

      case('oplflx')
       mass=16.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(OpRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)>0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


    case('oprflx')
       mass=16.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(OpRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass

      case('opflx')
       mass=16.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) =  &
                  (State_VGB(OpRhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(OpRhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(OpRhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)
             !if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
      VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass




      case('o2prflx')
       mass=32.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(O2pRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
       VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


       case('o2plflx')
       mass=32.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(O2pRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)>0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
       VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass

    case('o2pflx')
       mass=32.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) =  &
                  (State_VGB(O2pRhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(O2pRhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(O2pRhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)
             !if(x_BLK(i,j,k,iBLK)>0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
       VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
       VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass

        case('co2plflx')
       mass=44.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(CO2pRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)>0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
   VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


case('co2prflx')
       mass=44.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) = State_VGB(CO2pRho_,i,j,k,iBLK)* &
                  (State_VGB(rhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(rhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
             if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
             !tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
   VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


case('co2pflx')
       mass=44.0
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
             tmp1_BLK(i,j,k,iBLK) =  &
                  (State_VGB(CO2pRhoUx_,i,j,k,iBLK)*x_BLK(i,j,k,iBLK) &
                  +State_VGB(CO2pRhoUy_,i,j,k,iBLK)*y_BLK(i,j,k,iBLK) &
                  +State_VGB(CO2pRhoUz_,i,j,k,iBLK)*z_BLK(i,j,k,iBLK) &
                  )/R_BLK(i,j,k,iBLK)
            !if(x_BLK(i,j,k,iBLK)<0.0)tmp1_BLK(i,j,k,iBLK)=0.0
            ! tmp1_BLK(i,j,k,iBLK) = max(0.0, tmp1_BLK(i,j,k,iBLK))
          end do; end do; end do
       end do
      VarValue = calc_sphere('integrate', 360, Radius, tmp1_BLK)
   VarValue=VarValue*No2Si_V(UnitN_)*No2Si_V(UnitX_)**2*No2Si_V(UnitU_)/mass


   case default
       call stop_mpi('wrong logvarname')
    end select

   
 end subroutine user_get_log_var

 !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModPhysics, ONLY: rBody, No2Io_V, UnitRho_, BodyRho_I,UnitT_,UnitN_
    use ModMain, ONLY: Body1_
    use ModAdvance, ONLY: State_VGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, r_BLK, IsBoundaryBlock_IB
    use ModMain, ONLY: iTest, jTest, kTest, ProcTest, BlkTest, &
         GLOBALBLK
    use ModProcMH,   ONLY: iProc
    use ModSize, ONLY: nI, nJ, nK
    use ModMultiFluid

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

    character (len=*), parameter :: NameSub = 'user_set_plot_var'

    integer :: iVar, i, j, k, iIon
    logical :: oktest,oktest_me

    !-------------------------------------------------------------------
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_plot_var',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    IsFound=.true.

    select case(NameVar)
       !  case('hp')
       !     iVar=HpRho_
       !     iIon=1
       !  case('o2p')
       !     iVar=O2pRho_
       !     iIon=2
       !  case('op')
       !     iVar=OpRho_
       !     iIon=3
       !  case('co2p')
       !     iVar=CO2pRho_
       !     iIon=4

    case('co2')
       iVar=CO2_
    case('o')
       iVar=O_
!    case('h')
!       iVar=H_
!    case('oh')
!       iVar=Oh_
!    case('ohx')
!       iVar=Ohx_
!    case('hx')
!       iVar=Hx_
!    case('ox')
!       iVar=Ox_
!    case('co2x')
!       iVar=CO2x_

    case('rateo') 
       iVar=O_
    case('rateco2')
       iVar=CO2_
    !case('product') 
     !  PlotVar_G(1:nI,1:nJ,1:nK)=Productrate_CB(:,:,:,iBlock)

    case default
       IsFound= .false.
       call stop_mpi(NameSub//': unimplemented variable='//NameVar)
    end select
    !NameTecUnit = '[amu/cm3]'
    !NameIdlUnit = 'amu/cm3'

    NameTecUnit = '[cm-3.s-1]'
    NameIdlUnit = 'cm-3.s-1'
    !PlotVar_G(1:nI,1:nJ,1:nK)=Productrate_CB(:,:,:,iBlock)
    !PlotVar_G(1:nI,1:nJ,1:nK)=Ionizationrate_CBI(:,:,:,iBlock,iVar)/No2Io_V(UnitT_)*No2Io_V(UnitN_)
      PlotVar_G(1:nI,1:nJ,1:nK)= nDenNuSpecies_CBI(:,:,:,iBlock,iVar)*No2Io_V(UnitN_)
    !  PlotVar_G   = State_VGB(iVar,:,:,:,iBlock)/MassIon_I(iIon)


   ! PlotVarBody = BodyRho_I(iIon+1)

    !if(IsDimensional) PlotVar_G = PlotVar_G*No2Io_V(UnitRho_)

  end subroutine user_set_plot_var
  !============================================================================
  subroutine set_neutral_density(iBlock)

    use ModProcMH, ONLY : iProc
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
    use ModIO, ONLY : restart
    use ModPhysics
    use ModNumConst

    integer, intent(in):: iBlock

    real :: Rmax,CosSZA,Xp, chap_y, chap, sinSZA
    integer:: i, j, k, q
    character (len=*), parameter :: NameSub = 'set_neutral_density'
    logical:: DoTest, DoTestMe, DoTestCell

    !--------------------------------------------------------------------------

    if(iProc==PROCtest.and. iBlock==BLKtest)then
       call set_oktest(NameSub,DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if


    !if(DoTestCell)write(*,*)'set_neutral_density called,iBlock', iBlock


    !calculate neutral 

    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell= DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       
       !if(DoTestCell)then
       !   write(*,*)'iProc, iBlock, iTest, jTest, kTest=',iProc, iBlock, i, j, k 
       !   write(*,*)'I am inside testcell' 
       !   write(*,*)'beginning of set_neutral_density'
       !   write(*,*)'nDenNuSpecies_CBI=',nDenNuSpecies_CBI(iTest,jTest,kTest,iBlock,:)
       !   write(*,*)'Productrate_CBI=',Productrate_CB(iTest,jTest,kTest,iBlock)
       !   write(*,*)'Ionizationrate_CBI=',Ionizationrate_CBI(iTest,jTest,kTest,iBlock,:)
       !end if

       !write(*,*)' BodynDenNuSpecies_I=', BodynDenNuSpecies_I
       if(R_BLK(i,j,k,iBlock)<= Rbody)then
          nDenNuSpecies_CBI(i,j,k,iBlock,:)=&
               BodynDenNuSpecies_I

       else if(R_BLK(i,j,k,iBlock)< 5.0) then         
          nDenNuSpecies_CBI(i,j,k,iBlock,:)=&
               BodynDenNuSpecies_I*& 
               exp(-(R_BLK(i,j,k,iBlock)-Rbody)&
               /HNuSpecies_I)
          
       else
          nDenNuSpecies_CBI(i,j,k,globalBLK,:)=0.0

       end if


       if(DoTestCell)then
          write(*,*)'nDenNuSpecies_CBI=',nDenNuSpecies_CBI(iTest,jTest,kTest,iBlock,:)
          write(*,*)'HNuSpecies_I=', HNuSpecies_I
       end if 


    end do;end do;end do

    !    call neutral_density_averages  !calculate averaged neutral density

    ! calculate optical depth and producation rate
  !  do k=1,nK; do j=1,nJ; do i=1,nI
  !     cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,iBlock)))*&
  !          x_BLK(i,j,k,iBlock)/max(R_BLK(i,j,k,iBlock),1.0e-3)&
  !          +5.0e-4


 !      Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,iBlock,1:MaxNuSpecies)*&
 !           CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
 !           6.0e-3)/cosSZA        
 !      if( Optdep<11.5 .and. x_BLK(i,j,k,iBlock) > 0.0) then 
 !         Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-5)
 !      else
 !         Productrate_CB(i,j,k,iBlock) = 1.0e-5
 !      end if

!    end do; end do; end do

 do k=1,nK; do j=1,nJ; do i=1,nI
             !cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
             !     x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-4)&
             !     +1.0e-4
             !Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,globalBLK,1:MaxNuSpecies)*&
             !     CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
             !     6.0e-3)/cosSZA
             !if( Optdep<11.5 .and. x_BLK(i,j,k,globalBLK) > 0.0) then
             !   Productrate_CB(i,j,k,globalBLK) = max(exp(-Optdep), 1.0e-30)
             !else
             !   Productrate_CB(i,j,k,globalBLK) = 1.0e-30
             !end if

                       Optdep=HNuSpecies_I(CO2_)*nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)*&
                  CrossSection_I(CO2_)

             cosSZA = X_BLK(i,j,k,iBlock) /R_BLK(i,j,k,iBlock)
             if(UseChapman)then            
                if(Optdep > 13.8) then 
                   Productrate_CB(i,j,k,iBlock) = 1.0e-6*&
                        max(cosSZA, 1.0e-6)                
                else
                   Xp=R_BLK(i,j,k,iBlock)/HNuSpecies_I(1)
                   chap_y= sqrt(0.5*Xp)*abs(cosSZA)
                   if(cosSZA .gt.0.0)then   !SZA<90 deg (equation 13)
                      if(chap_y<8.0)then
                         chap = sqrt(3.1415926/2.0*Xp)*&
                              (1.0606963+0.5564383*chap_y)/&
                              (1.0619896+1.7245609*chap_y+chap_y*chap_y)
                      elseif(chap_y<100.0)then
                         chap = sqrt(3.1415926/2.0*Xp)*&
                              0.56498823/(0.6651874+chap_y)
                      else
                         chap=0.0
                      end if
                   elseif(cosSZA > -0.5) then
                      ! 120 >SZA > 90 deg (equation 15) Smith and Smith, 1972
                      sinSZA = sqrt(1.0 - cosSZA**2)
                      if(chap_y<8.0)then
                         chap =sqrt(2.0*3.1415926*Xp)* &
                              (sqrt(sinSZA)*exp(Xp*(1.0-sinSZA)) &
                              -0.5*(1.0606963+0.5564383*chap_y)/ &
                              (1.0619896+1.7245609*chap_y+chap_y*chap_y))
                      elseif(chap_y<100.0)then
                         chap =sqrt(2.0*3.1415926*Xp)* &
                              (sqrt(sinSZA)*exp(Xp*(1.0-sinSZA)) &
                              -0.5*0.56498823/(0.6651874+chap_y))
                      else
                         chap =0.0
                      end if
                   else
                      chap =1.0e10
                   end if
                   
                   Optdep=Optdep*chap      
                   Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-6)
                end if
             else  !old way of optical depth calculation, 
                !production rate drops to 1.0e-5 for SZA > 90
                
                cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,iBlock)))*&
                     x_BLK(i,j,k,iBlock)/max(R_BLK(i,j,k,iBlock),1.0e-5)&
                     +1.0e-5
                Optdep =max( sum(nDenNuSpecies_CBI(i,j,k,iBlock,1:MaxNuSpecies)*&
                     CrossSection_I(1:MaxNuSpecies)*HNuSpecies_I(1:MaxNuSpecies)),&
                     6.0e-3)/cosSZA
                if( Optdep<11.5 .and. x_BLK(i,j,k,iBlock) > 0.0) then 
                   Productrate_CB(i,j,k,iBlock) = max(exp(-Optdep), 1.0e-5)
                else
                   Productrate_CB(i,j,k,iBlock) = 1.0e-30
                end if
             end if
             
          end do; end do; end do


    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell= DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

          nu_BLK(i,j,k,iBlock)=(nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)+&
               nDenNuSpecies_CBI(i,j,k,iBlock,O_))*nu0

           nu_BLK(i,j,k,iBlock)=&
               sum(nDenNuSpecies_CBI(i,j,k,iBlock,:))*nu0

             nDenNuSpecies_CBI(i,j,k,iBlock,O_)= &
               nDenNuSpecies_CBI(i,j,k,iBlock,O_)+ &
               nDenNuSpecies_CBI(i,j,k,iBlock,Oh_)+&
               nDenNuSpecies_CBI(i,j,k,iBlock,Ohx_)

        !  nDenNuSpecies_CBI(i,j,k,iBlock,H_)= 1.0e-5

          Ionizationrate_CBI(i,j,k,iBlock,O_)= &
               Rate_I(O_hv__Op_em_)&
               *nDenNuSpecies_CBI(i,j,k,iBlock,O_)&
               *Productrate_CB(i,j,k,iBlock)

          Ionizationrate_CBI(i,j,k,iBlock,CO2_)= &
               Rate_I(CO2_hv__CO2p_em_)&
               *nDenNuSpecies_CBI(i,j,k,iBlock,CO2_)&
                *Productrate_CB(i,j,k,iBlock)


          if(DoTestCell)write(*,*)'In set_neutral_density, Ionizationrate_CBI=',&
               Ionizationrate_CBI(iTest,jTest,kTest,BLKTest,:)/No2Io_V(UnitT_)*No2Io_V(UnitN_)

    end do; end do; end do 

  end subroutine set_neutral_density
 !===================================================================

end module ModUser
