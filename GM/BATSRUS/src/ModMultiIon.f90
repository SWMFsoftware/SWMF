module ModMultiIon

  ! Calculate source terms for multi-ion MHD. 
  ! For sake of numerical stability this should be done with 
  ! a point-implicit scheme. 
  ! Allow for extra point-implicit sources in ModUser

  use ModMultiFluid
  use ModMain, ONLY: UseUserSource
  use ModUser, ONLY: user_calc_sources, user_init_point_implicit

  implicit none

  private

  public:: multi_ion_sources
  public:: multi_ion_init_point_impl

  real    :: CollisionCoefDim = 1.0, CollisionCoef = 0.0

contains
  !===========================================================================
  subroutine multi_ion_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModProcMH,  ONLY: iProc
    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplMatrixSet
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, Test_String, BlkTest, ProcTest
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModAdvance, ONLY: B0XCell_BLK, B0YCell_BLK, B0ZCell_BLK
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB
    use ModPhysics, ONLY: ElectronCharge, gm1, inv_gm1, &
         Si2No_V, No2Si_V, UnitTemperature_, UnitT_
    use ModMain,    ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: cross_product

    ! Variables for multi-ion MHD
    real    :: InvCharge, NumDens, InvNumDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, uIon2_D, u_D, uPlus_D, uPlusHallU_D
    real, dimension(3) :: Current_D, Force_D
    real, dimension(nIonFluid) :: NumDens_I, InvRho_I, Ux_I, Uy_I, Uz_I, Temp_I

    integer :: iBlock, i, j, k, jFluid, iFirstIons
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, TemperatureCoef, Heating
    real :: CollisionRate_II(nIonFluid, nIonFluid), CollisionRate

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    logical :: DoTest, DoTestMe
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    ! Add user defined point implicit source terms here
    ! Explicit user sources are added in calc_sources
    if(UsePointImplicit .and. UseUserSource) call user_calc_sources

    ! Add source term n_s*(- u_+ - w_H + u_s )xB for multi-ions
    ! where u_+ is the number density weighted average ion velocity,
    ! and w_H = -J/(e n_e) is the Hall velocity. Here
    ! e is the electron charge and n_e is the electron number density.

    InvCharge = 1.0/ElectronCharge

    ! Rate = n*CoefDim / T^1.5 with T [K], n [/cc] and Rate [1/s]
    CollisionCoef = CollisionCoefDim &
         /No2Si_V(UnitTemperature_)**1.5/Si2No_V(UnitT_)

    do jFluid = 1, nIonFluid
       do iFluid = 1, nIonFluid
          CollisionRate_II(iFluid, jFluid) = CollisionCoef* &
               MassFluid_I(iFluid)*MassFluid_I(jFluid) &
               /(MassFluid_I(iFluid)+MassFluid_I(jFluid))
       end do
    end do

    ! Do not add
    iFirstIons = 1
    if(TypeFluid_I(1) == 'ion')iFirstIons = 2

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,globalBLK)

       if(TypeFluid_I(1) == 'ion')then
          ! Get first fluid quantities
          State_V(Rho_) = State_V(Rho_) &
               - sum(State_V(iRhoIon_I(2:nIonFluid)))
          State_V(RhoUx_) = State_V(RhoUx_) &
               - sum(State_V(iRhoUxIon_I(2:nIonFluid)))
          State_V(RhoUy_) = State_V(RhoUy_) &
               - sum(State_V(iRhoUyIon_I(2:nIonFluid)))
          State_V(RhoUz_) = State_V(RhoUz_) &
               - sum(State_V(iRhoUzIon_I(2:nIonFluid)))
          State_V(P_) = State_V(P_) &
               - sum(State_V(iPIon_I(2:nIonFluid)))
       end if

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

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       CollisionRate = CollisionCoef

       ! Calculate the source term for all the ion fluids
       do iFluid = iFirstIons, nIonFluid
          call select_fluid
          uIon_D = (/ Ux_I(iFLuid),  Uy_I(iFluid), Uz_I(iFluid) /)
          u_D    = uIon_D - uPlusHallU_D

          Force_D = &
               ElectronCharge*NumDens_I(iFluid)*cross_product(u_D, FullB_D) 

          Heating = 0.0

          if(.false.)then
             do jFluid = 1, nIonFluid
                if(jFluid == iFluid) CYCLE
                
                ! Add collisional term
                uIon2_D = (/ Ux_I(jFLuid),  Uy_I(jFluid), Uz_I(jFluid) /)

                if(  InvNumDens*NumDens_I(iFluid) < 0.01 .or. &
                     InvNumDens*NumDens_I(jFluid) < 0.01) then
                   CollisionRate = 100.0 * NumDens_I(iFluid) * NumDens_I(jFluid)
                else
                   CollisionRate = CollisionRate_II(iFluid, jFluid) * &
                        NumDens_I(iFluid) * NumDens_I(jFluid) &
                        * TemperatureCoef
                end if

                Force_D = Force_D + CollisionRate*(uIon2_D - uIon_D)

                Heating = Heating + CollisionRate* &
                     ( 2*(Temp_I(jFluid) - Temp_I(iFluid)) &
                     + gm1*sum((uIon2_D - uIon_D)**2) )
             end do
          end if

          Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) = &
               Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) + Force_D

          Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) + Heating

          Source_VC(Energy_-1+iFluid,i,j,k) = &
               Source_VC(Energy_-1+iFluid,i,j,k) + sum(Force_D*uIon_D) &
               + inv_gm1*Heating

       end do
    end do; end do; end do

    if(DoTestMe)then
       write(*,*)NameSub,' CollisionCoef=',CollisionCoef
       write(*,*)NameSub,' CollisionRate=',CollisionRate
       write(*,*)NameSub,' AverageTemp  =',AverageTemp
       write(*,*)NameSub,' AverageTempDim=', &
            AverageTemp*No2Si_V(UnitTemperature_)
    end if

  end subroutine multi_ion_sources
  !===========================================================================
  subroutine multi_ion_init_point_impl

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet

    logical :: IsPointImpl_V(nVar)
    integer :: iFirst, iVar, iPointImplVar, nPointImplVar
    !------------------------------------------------------------------------

    IsPointImpl_V = .false.

    if(UseUserSource)then
       call user_init_point_implicit
       if(allocated(iVarPointImpl_I)) &
            IsPointImpl_V(iVarPointImpl_I) = .true.
    end if

    IsPointImplMatrixSet = .false.

    ! All ion momenta are implicit
    if(TypeFluid_I(1) == 'ions')then
       iFirst = 1
    else
       iFirst = 2
    end if

    IsPointImpl_V(iRhoUxIon_I(iFirst:)) = .true.
    IsPointImpl_V(iRhoUyIon_I(iFirst:)) = .true.
    IsPointImpl_V(iRhoUzIon_I(iFirst:)) = .true.
    IsPointImpl_V(iRhoUxIon_I(iFirst:)) = .true.
    IsPointImpl_V(iPIon_I(iFirst:))     = .true.

    nPointImplVar = count(IsPointImpl_V)

    allocate(iVarPointImpl_I(nPointImplVar))

    iPointImplVar = 0
    do iVar = 1, nVar
       if(.not. IsPointImpl_V(iVar)) CYCLE
       iPointImplVar = iPointImplVar + 1
       iVarPointImpl_I(iPointImplVar) = iVar
    end do

  end subroutine multi_ion_init_point_impl

end module ModMultiIon
