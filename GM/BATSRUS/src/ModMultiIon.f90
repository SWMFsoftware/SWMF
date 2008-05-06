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
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, &
         iTest, jTest, kTest, Test_String, BlkTest, ProcTest
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModAdvance, ONLY: B0XCell_BLK, B0YCell_BLK, B0ZCell_BLK
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB, x_BLK, y_BLK, z_BLK
    use ModPhysics, ONLY: ElectronCharge, gm1, inv_gm1, &
         Si2No_V, No2Si_V, UnitTemperature_, UnitT_
    use ModMain,    ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: cross_product

    ! Variables for multi-ion MHD
    real    :: InvCharge, NumDens, InvNumDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, uIon2_D, u_D, uPlus_D, uPlusHallU_D
    real, dimension(3) :: Current_D, Force_D
    real, dimension(nIonFluid) :: NumDens_I, InvRho_I, Ux_I, Uy_I, Uz_I, Temp_I

    integer :: iBlock, i, j, k, iIon, jIon, iRhoUx, iRhoUz, iP, iEnergy
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, TemperatureCoef, Heating
    real :: CollisionRate_II(nIonFluid, nIonFluid), CollisionRate

    character (len=*), parameter :: NameSub = 'multi_ion_sources'
    logical :: DoTest, DoTestMe, DoTestCell
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if
    DoTestCell = .false.

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

    do jIon = 1, nIonFluid
       do iIon = 1, nIonFluid
          CollisionRate_II(iIon, jIon) = CollisionCoef* &
               MassIon_I(iIon)*MassIon_I(jIon) &
               /(MassIon_I(iIon)+MassIon_I(jIon))
       end do
    end do

    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,iBlock)

       ! Total magnetic field
       FullB_D = State_V(Bx_:Bz_) + (/ &
            B0xCell_BLK(i,j,k,iBlock),&
            B0yCell_BLK(i,j,k,iBlock),&
            B0zCell_BLK(i,j,k,iBlock) /)

       ! calculate number densities
       NumDens_I  = State_V(iRhoIon_I) / MassIon_I
       NumDens    = sum(NumDens_I)
       InvNumDens = 1.0/NumDens

       Temp_I     = State_V(iPIon_I)/NumDens_I
       AverageTemp= sum(State_V(iPIon_I))*InvNumDens

       if(AverageTemp <= 0.0)then
          write(*,*)'ERROR: AverageTemp =',AverageTemp
          write(*,*)'i,j,k,iBlock,iProc =',i,j,k,iBlock,iProc
          write(*,*)'x,y,z              =', &
               x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
          write(*,*)'iRhoIon_I          =',iRhoIon_I
          write(*,*)'RhoIon_I           =',State_V(iRhoIon_I)
          write(*,*)'MassIon_I          =',MassIon_I
          write(*,*)'NumDens_I          =',NumDens_I
          write(*,*)'iPIon_I            =',iPIon_I
          write(*,*)'PIon_I             =',State_V(iPIon_I)
          write(*,*)'Temp_I             =',Temp_I
          call stop_mpi(NameSub//': non-positive average temperature')
       end if

       InvRho_I = 1.0/State_V(iRhoIon_I)
       Ux_I  = InvRho_I*State_V(iUxIon_I)
       Uy_I  = InvRho_I*State_V(iUyIon_I)
       Uz_I  = InvRho_I*State_V(iUzIon_I)

       ! calculate the average positive charge velocity
       uPlus_D(x_) = InvNumDens* sum(NumDens_I*Ux_I)
       uPlus_D(y_) = InvNumDens* sum(NumDens_I*Uy_I)
       uPlus_D(z_) = InvNumDens* sum(NumDens_I*Uz_I)

       ! Add the Hall velocity -J/(e n)
       ! old version: call get_current(i,j,k,iBlock,Current_D)
       Current_D = vInv_CB(i,j,k,iBlock)*&
            ( bCrossArea_DX(:,i+1,j,k) - bCrossArea_DX(:,i,j,k) &
            + bCrossArea_DY(:,i,j+1,k) - bCrossArea_DY(:,i,j,k) &
            + bCrossArea_DZ(:,i,j,k+1) - bCrossArea_DZ(:,i,j,k))

       uPlusHallU_D = uPlus_D - InvNumDens*InvCharge*Current_D

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       CollisionRate = CollisionCoef

       ! Calculate the source term for all the ion fluids
       do iIon = 1, nIonFluid
          ! call select_fluid
          uIon_D = (/ Ux_I(iIon),  Uy_I(iIon), Uz_I(iIon) /)
          u_D    = uIon_D - uPlusHallU_D

          Force_D = &
               ElectronCharge*NumDens_I(iIon)*cross_product(u_D, FullB_D) 

          if(DoTestCell)then
             write(*,*) NameSub,' iIon =', iIon
             write(*,*)'bxO_DX(:,i+1)=', bCrossArea_DX(:,i+1,j,k)
             write(*,*)'bxO_DX(:,i  )=', bCrossArea_DX(:,i,j,k)
             write(*,*)'bxO_DY(:,j+1)=', bCrossArea_DY(:,i,j+1,k)
             write(*,*)'bxO_DY(:,j  )=', bCrossArea_DY(:,i,j,k)
             write(*,*)'bxO_DZ(:,k+1)=', bCrossArea_DZ(:,i,j,k+1)
             write(*,*)'bxO_DZ(:,k  )=', bCrossArea_DZ(:,i,j,k)
             write(*,*) NameSub,' uPlus_D  =', uPlus_D
             write(*,*) NameSub,' uIon_D   =', uIon_D
             write(*,*) NameSub,' Current_D=', Current_D
             write(*,*) NameSub,' uPlusHall=', uPlusHallU_D
             write(*,*) NameSub,' Force_D  =', Force_D
          end if

          Heating = 0.0

          if(.false.)then
             do jIon = 1, nIonFluid
                if(jIon == iIon) CYCLE
                
                ! Add collisional term
                uIon2_D = (/ Ux_I(jIon),  Uy_I(jIon), Uz_I(jIon) /)

                if(  InvNumDens*NumDens_I(iIon) < 0.01 .or. &
                     InvNumDens*NumDens_I(jIon) < 0.01) then
                   CollisionRate = 100 * NumDens_I(iIon) * NumDens_I(jIon)
                else
                   CollisionRate = CollisionRate_II(iIon, jIon) * &
                        NumDens_I(iIon) * NumDens_I(jIon) &
                        * TemperatureCoef
                end if

                Force_D = Force_D + CollisionRate*(uIon2_D - uIon_D)

                Heating = Heating + CollisionRate* &
                     ( 2*(Temp_I(jIon) - Temp_I(iIon)) &
                     + gm1*sum((uIon2_D - uIon_D)**2) )
             end do

             iP = iPIon_I(iIon)
             Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) + Heating

          end if

          iRhoUx = iRhoUxIon_I(iIon); iRhoUz = iRhoUzIon_I(iIon)
          Source_VC(iRhoUx:iRhoUz,i,j,k) = Source_VC(iRhoUx:iRhoUz,i,j,k) &
               + Force_D

          iEnergy = Energy_-2+iIon+IonFirst_
          Source_VC(iEnergy,i,j,k) = Source_VC(iEnergy,i,j,k) &
               + sum(Force_D*uIon_D) &
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
    integer :: iVar, iPointImplVar, nPointImplVar
    !------------------------------------------------------------------------

    IsPointImpl_V = .false.

    if(UseUserSource)then
       call user_init_point_implicit
       if(allocated(iVarPointImpl_I)) &
            IsPointImpl_V(iVarPointImpl_I) = .true.
    end if

    IsPointImplMatrixSet = .false.

    ! All ion momenta and pressures are implicit
    IsPointImpl_V(iRhoUxIon_I) = .true.
    IsPointImpl_V(iRhoUyIon_I) = .true.
    IsPointImpl_V(iRhoUzIon_I) = .true.
    IsPointImpl_V(iPIon_I)     = .true.

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
