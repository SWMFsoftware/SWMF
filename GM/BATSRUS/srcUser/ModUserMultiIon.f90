!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! do electro-magnetic force term in multi-ion point-implicitly

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_init_point_implicit,        &
       IMPLEMENTED4 => user_set_ics

  use ModMultiFluid

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'POINT IMPLICIT MULTI-ION SOURCE, Toth and Ma'


  real :: CollisionRate = 0.001
  logical :: IsAnalytic = .true.

contains

  !==========================================================================
  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#FRICTION")
          call read_var('CollisionRate',CollisionRate)
       case('#DSDU')
          call read_var('IsAnalytic',IsAnalytic)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserPointImplicit: unknown command='//&
               NameCommand)
       end select
    end do

  end subroutine user_read_inputs
  !==========================================================================
  subroutine user_calc_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, Test_String
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModAdvance, ONLY: B0XCell_BLK, B0YCell_BLK, B0ZCell_BLK
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB
    use ModPhysics, ONLY: ElectronCharge, inv_gm1
    use ModMain,    ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: cross_product

    ! Variables for multi-ion MHD
    real    :: InvCharge, NumDens, InvNumDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, u_D, uPlus_D, uPlusHallU_D, Force_D
    real, dimension(3) :: Current_D
    real, dimension(nIonFluid) :: NumDens_I, InvRho_I, Ux_I, Uy_I, Uz_I, &
         Temp_I

    integer :: iBlock, i, j, k
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, Heating
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    ! Add source term n_s*(- u_+ - w_H + u_s )xB for multi-ions
    ! where u_+ is the number density weighted average ion velocity,
    ! and w_H = -J/(e n_e) is the Hall velocity. Here
    ! e is the electron charge and n_e is the electron number density.

    InvCharge = 1.0/ElectronCharge

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,globalBLK)

       ! Total magnetic field
       FullB_D = State_V(Bx_:Bz_) + (/ &
            B0xCell_BLK(i,j,k,globalBLK),&
            B0yCell_BLK(i,j,k,globalBLK),&
            B0zCell_BLK(i,j,k,globalBLK) /)

       ! calculate number densities
       NumDens_I  = State_V(iRhoIon_I) / MassFluid_I(1:nIonFluid)
       NumDens    = sum(NumDens_I)
       InvNumDens = 1.0/NumDens

       Temp_I     = State_V(iPIon_I)/NumDens_I
       AverageTemp= sum(State_V(iPIon_I))*InvNumDens

       InvRho_I = 1.0/State_V(iRhoIon_I)
       Ux_I  = InvRho_I*State_V(iUxIon_I)
       Uy_I  = InvRho_I*State_V(iUyIon_I)
       Uz_I  = InvRho_I*State_V(iUzIon_I)

       ! calculate the average positive charge velocity
       uPlus_D(x_) = InvNumDens* sum(NumDens_I*Ux_I)
       uPlus_D(y_) = InvNumDens* sum(NumDens_I*Uy_I)
       uPlus_D(z_) = InvNumDens* sum(NumDens_I*Uz_I)

       ! Add the Hall velocity -J/(e n)
       if(index(Test_String,'newj') > 0)then
          Current_D = vInv_CB(i,j,k,globalBLK)*&
               ( bCrossArea_DX(:,i+1,j,k) - bCrossArea_DX(:,i,j,k) &
               + bCrossArea_DY(:,i,j+1,k) - bCrossArea_DY(:,i,j,k) &
               + bCrossArea_DZ(:,i,j,k+1) - bCrossArea_DZ(:,i,j,k))
       else
          call get_current(i,j,k,GlobalBlk,Current_D)
       end if
       uPlusHallU_D = uPlus_D - InvNumDens*InvCharge*Current_D

       ! Calculate the source term for all the ion fluids
       do iFluid = 1, nIonFluid
          call select_fluid
          uIon_D = (/ Ux_I(iFLuid),  Uy_I(iFluid), Uz_I(iFluid) /)
          u_D    = uIon_D - uPlusHallU_D

          Force_D = &
               ElectronCharge*NumDens_I(iFluid)*cross_product(u_D, FullB_D) 

          ! Add simplified collisional term
          Force_D = Force_D + CollisionRate*(uPlus_D - uIon_D) &
               * MassFluid_I(iFluid) * NumDens_I(iFluid) &
               * (NumDens - NumDens_I(iFluid))

          Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) = &
               Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) + Force_D

          Heating = 2*CollisionRate*(AverageTemp - Temp_I(iFluid)) &
               * NumDens_I(iFluid)* (NumDens - NumDens_I(iFluid))

          Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) + Heating

          Source_VC(Energy_-1+iFluid,i,j,k) = &
               Source_VC(Energy_-1+iFluid,i,j,k) + sum(Force_D*uIon_D) &
               + inv_gm1*Heating

          if(.not.IsAnalytic) CYCLE

          Coef   = ElectronCharge/MassFluid_I(iFluid) &
               *(1.0 - InvNumDens*NumDens_I(iFluid))
          CoefBx = Coef*FullB_D(x_)
          CoefBy = Coef*FullB_D(y_)
          CoefBz = Coef*FullB_D(z_)
          DsDu_VVC(iRhoUx, iRhoUy, i,j,k) =   CoefBz
          DsDu_VVC(iRhoUx, iRhoUz, i,j,k) = - CoefBy
          DsDu_VVC(iRhoUy, iRhoUz, i,j,k) =   CoefBx
          DsDu_VVC(iRhoUy, iRhoUx, i,j,k) = - CoefBz
          DsDu_VVC(iRhoUz, iRhoUx, i,j,k) =   CoefBy
          DsDu_VVC(iRhoUz, iRhoUy, i,j,k) = - CoefBx

          Coef = 0.5*CollisionRate*(NumDens - NumDens_I(iFluid))
          DsDu_VVC(iP,    iP,     i,j,k) = -2*Coef
          DsDu_VVC(iRhoUx,iRhoUx, i,j,k) = -Coef
          DsDu_VVC(iRhoUy,iRhoUy, i,j,k) = -Coef
          DsDu_VVC(iRhoUz,iRhoUz, i,j,k) = -Coef

       end do
       if( .not. IsAnalytic) CYCLE
       Coef   = -ElectronCharge/MassFluid_I(1) &
            *InvNumDens*NumDens_I(1)
       CoefBx = Coef*FullB_D(x_)
       CoefBy = Coef*FullB_D(y_)
       CoefBz = Coef*FullB_D(z_)
       DsDu_VVC(RhoUx_, OpRhoUy_, i,j,k) =   CoefBz
       DsDu_VVC(RhoUx_, OpRhoUz_, i,j,k) = - CoefBy
       DsDu_VVC(RhoUy_, OpRhoUz_, i,j,k) =   CoefBx
       DsDu_VVC(RhoUy_, OpRhoUx_, i,j,k) = - CoefBz
       DsDu_VVC(RhoUz_, OpRhoUx_, i,j,k) =   CoefBy
       DsDu_VVC(RhoUz_, OpRhoUy_, i,j,k) = - CoefBx

       Coef   = -ElectronCharge/MassFluid_I(2) &
            *InvNumDens*NumDens_I(2)
       CoefBx = Coef*FullB_D(x_)
       CoefBy = Coef*FullB_D(y_)
       CoefBz = Coef*FullB_D(z_)
       DsDu_VVC(OpRhoUx_, RhoUy_, i,j,k) =   CoefBz
       DsDu_VVC(OpRhoUx_, RhoUz_, i,j,k) = - CoefBy
       DsDu_VVC(OpRhoUy_, RhoUz_, i,j,k) =   CoefBx
       DsDu_VVC(OpRhoUy_, RhoUx_, i,j,k) = - CoefBz
       DsDu_VVC(OpRhoUz_, RhoUx_, i,j,k) =   CoefBy
       DsDu_VVC(OpRhoUz_, RhoUy_, i,j,k) = - CoefBx

       Coef = -0.5*CollisionRate*NumDens_I(1)
       DsDu_VVC(P_,     OpP_,     i, j, k) = -2*Coef
       DsDu_VVC(RhoUx_, OpRhoUx_, i, j, k) = -Coef
       DsDu_VVC(RhoUy_, OpRhoUy_, i, j, k) = -Coef
       DsDu_VVC(RhoUz_, OpRhoUz_, i, j, k) = -Coef

       Coef = -0.5*CollisionRate*NumDens_I(2)
       DsDu_VVC(OpP_,     P_,     i, j, k) = -2*Coef
       DsDu_VVC(OpRhoUx_, RhoUx_, i, j, k) = -Coef
       DsDu_VVC(OpRhoUy_, RhoUy_, i, j, k) = -Coef
       DsDu_VVC(OpRhoUz_, RhoUz_, i, j, k) = -Coef

    end do; end do; end do

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !------------------------------------------------------------------------

    ! All ion momenta are implicit
    allocate(iVarPointImpl_I(4*nIonFluid))

    do iFluid = 1, nIonFluid
       iVarPointImpl_I(4*iFluid-3) = iRhoUx_I(iFluid)
       iVarPointImpl_I(4*iFluid-2) = iRhoUy_I(iFluid)
       iVarPointImpl_I(4*iFluid-1) = iRhoUz_I(iFluid)
       iVarPointImpl_I(4*iFluid)   = iP_I(iFluid)
    end do

    IsPointImplMatrixSet = IsAnalytic

  end subroutine user_init_point_implicit

  !=====================================================================
  subroutine user_set_ics
    use ModMain,     ONLY: globalBLK
    use ModGeometry, ONLY: r_BLK
    use ModVarIndexes
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: rBody, BodyRho_I, BodyP_I

    implicit none

    integer :: iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    where(r_BLK(:,:,:,iBlock) < 3*rBody)
       State_VGB(Rho_  ,:,:,:,iBlock) = BodyRho_I(1)
       State_VGB(OpRho_,:,:,:,iBlock) = BodyRho_I(2)
       State_VGB(RhoUx_,:,:,:,iBlock) = 0.0
       State_VGB(RhoUy_,:,:,:,iBlock) = 0.0
       State_VGB(RhoUz_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUx_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUy_,:,:,:,iBlock) = 0.0
       State_VGB(OpRhoUz_,:,:,:,iBlock) = 0.0
       State_VGB(P_,      :,:,:,iBlock) = BodyP_I(1)
       State_VGB(OpP_,    :,:,:,iBlock) = BodyP_I(2)
    end where

  end subroutine user_set_ics

end module ModUser
