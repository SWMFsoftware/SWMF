module ModUser
  ! This is the user module for Mars 

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics,                    &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_init_point_implicit

  use ModMultiFluid

  use ModNumConst, ONLY: cSqrtHalf, cDegToRad

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: NameUserModule = &
       'Mars 4 Fluids MHD code, Dalal Najib'

  real    :: CollisionCoefDim = 1.0, CollisionCoef

contains

  subroutine user_calc_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModProcMH,  ONLY: iProc
    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, &
         Test_String, iTest, jTest, kTest, BlkTest, ProcTest, n_step
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModAdvance, ONLY: B0XCell_BLK, B0YCell_BLK, B0ZCell_BLK
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB
    use ModPhysics, ONLY: ElectronCharge, gm1, inv_gm1, &
         Si2No_V, No2Si_V, UnitTemperature_, UnitT_, BodyRho_I, BodyP_I
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

    real :: RhoMin = 0.0, pMin = 0.0

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    logical :: DoTest, DoTestMe, DoTestCell
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

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
    if(TypeFluid_I(1) == 'ion')then
       iFirstIons = 2
       RhoMin = BodyRho_I(1) - sum(BodyRho_I(2:nFluid))
       pMin   = BodyP_I(1) - sum(BodyP_I(2:nFluid))
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest
      
       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,globalBLK)

       if(DoTestCell)write(*,*) NameSub,' orig State_V=',State_V

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

          if(State_V(Rho_) < RhoMin .or. State_V(P_)< pMin)then
             State_V(Rho_)=RhoMin
             State_V(P_)  =pMin
          end if
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

       if(DoTestCell .or. AverageTemp<=0.0)then
          write(*,*) NameSub,' i,j,k,iBlock,iProc=',i,j,k,iBlock,iProc
          write(*,*) NameSub,' State_V  =',State_V
          write(*,*) NameSub,' iRhoIon_I=',iRhoIon_I
          write(*,*) NameSub,' Rho(ions)=',State_V(iRhoIon_I)
          write(*,*) NameSub,' MassFluid=',MassFluid_I
          write(*,*) NameSub,' NumDens_I=',NumDens_I
          write(*,*) NameSub,' Ux_I     =',Ux_I
          write(*,*) NameSub,' Uy_I     =',Uy_I
          write(*,*) NameSub,' Uz_I     =',Uz_I
          write(*,*) NameSub,' Temp_I   =',Temp_I
          write(*,*) NameSub,' uPlus_D    =',uPlus_D
          write(*,*) NameSub,' Current_D  =',Current_D
          write(*,*) NameSub,' AverageTemp=',AverageTemp
          if(AverageTemp<=0.0) &
               call stop_mpi(NameSub//': average temperature is non-positive')
       end if

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       CollisionRate = CollisionCoef

       ! Calculate the source term for all the ion fluids
       do iFluid = iFirstIons, nIonFluid
          call select_fluid
          uIon_D = (/ Ux_I(iFLuid),  Uy_I(iFluid), Uz_I(iFluid) /)
          u_D    = uIon_D - uPlusHallU_D

          Force_D = &
               ElectronCharge*NumDens_I(iFluid)*cross_product(u_D, FullB_D) 

          if(DoTestCell)then
             write(*,*) NameSub,' iFluid=',iFluid
             write(*,*) NameSub,' Hall Force=',Force_D
          end if

          Heating = 0.0

          do jFluid = 1, nIonFluid
             if(jFluid == iFluid) CYCLE

             ! Add collisional term
             uIon2_D = (/ Ux_I(jFLuid),  Uy_I(jFluid), Uz_I(jFluid) /)

             ! In low density regions make the fluids stick together
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


             if(DoTestCell)then
                write(*,*) NameSub,'    jFluid=',jFluid
                write(*,*) NameSub,'    CollisionRate=',CollisionRate
                write(*,*) NameSub,'    Friction     =', &
                     CollisionRate*(uIon2_D - uIon_D)
                write(*,*) NameSub,'    Heat exchange=', &
                     ( 2*(Temp_I(jFluid) - Temp_I(iFluid)) &
                     + gm1*sum((uIon2_D - uIon_D)**2) )
                write(*,*) NameSub,'    New Force_D  =',Force_D
                write(*,*) NameSub,'    New Heating  =',Heating
             end if

          end do

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

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !------------------------------------------------------------------------

    ! All ion momenta are implicit
    if(TypeFluid_I(1) == 'ions')then
       allocate(iVarPointImpl_I(4*nIonFluid))

       do iFluid = 1, nIonFluid
          iVarPointImpl_I(4*iFluid-3) = iRhoUx_I(iFluid)
          iVarPointImpl_I(4*iFluid-2) = iRhoUy_I(iFluid)
          iVarPointImpl_I(4*iFluid-1) = iRhoUz_I(iFluid)
          iVarPointImpl_I(4*iFluid)   = iP_I(iFluid)
       end do
    else
       allocate(iVarPointImpl_I(4*(nIonFluid-1)))
       do iFluid = 1, nIonFluid-1
          iVarPointImpl_I(4*iFluid-3) = iRhoUx_I(iFluid+1)
          iVarPointImpl_I(4*iFluid-2) = iRhoUy_I(iFluid+1)
          iVarPointImpl_I(4*iFluid-1) = iRhoUz_I(iFluid+1)
          iVarPointImpl_I(4*iFluid)   = iP_I(iFluid+1)
       end do
    end if

    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

  !===========================================================================
  subroutine user_set_ICs
    use ModProcMH, ONLY : iProc
    use ModMain
    use ModAdvance
    use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK,true_cell
    use ModIO, ONLY : restart
    use ModPhysics
    use ModNumConst




    real :: CosSZA
    !real::O2pRho_dim,OpRho_dim,CO2pRho_dim,HpRho_dim
    integer :: i,j,k
    integer :: iBoundary

    !--------------------------------------------------------------------------
    !do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
    !   if (R_BLK(i,j,k,globalBLK)< Rbody) then
    !      !cosSZA=(cHalf+sign(cHalf,x_BLK(i,j,k,globalBLK)))*&
    !      ! x_BLK(i,j,k,globalBLK)/max(R_BLK(i,j,k,globalBLK),1.0e-3)+&
    !      ! 1.0e-3
    !      State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
    !
    !   else
    !      State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)
    !     ! State_VGB(Ux_:bz_,i,j,k,globalBLK)   =0.0          
    !   end if
    !end do;end do; end do;


    do k=1,nK; do j=1,nJ; do i=1,nI

       if (true_cell(i,j,k,globalBLK).and. &
            R_BLK(i,j,k,globalBLK)<1.5*Rbody) then

          State_VGB(iRho_I,i,j,k,globalBLK) = BodyRho_I
          State_VGB(iP_I  ,i,j,k,globalBLK) = BodyP_I
          State_VGB(iRhoUx_I,i,j,k,globalBLK) = 0.0
          State_VGB(iRhoUy_I,i,j,k,globalBLK) = 0.0
          State_VGB(iRhoUz_I,i,j,k,globalBLK) = 0.0

!!$          O2pRho_dim= State_VGB(O2pRho_ ,i,j,k,globalBLK)*No2Io_V(UnitRho_)
!!$          OpRho_dim= State_VGB(OpRho_ ,i,j,k,globalBLK)*No2Io_V(UnitRho_)
!!$          CO2pRho_dim= State_VGB(CO2pRho_,i,j,k,globalBLK)&
!!$                *No2Io_V(UnitRho_)         
!!$          HpRho_dim=(State_VGB(Rho_ ,i,j,k,globalBLK)&
!!$               - State_VGB(O2pRho_ ,i,j,k,globalBLK)&
!!$               - State_VGB(OpRho_ ,i,j,k,globalBLK)&
!!$               - State_VGB(CO2pRho_ ,i,j,k,globalBLK) )*No2Io_V(UnitRho_)
!!$          
!!$          write(*,*)'No2Io_V(UnitRho_)=',No2Io_V(UnitRho_)
!!$          write(*,*)'State_VGB(O2pRho_ ,i,j,k,globalBLK)=',&
!!$               State_VGB(O2pRho_ ,i,j,k,globalBLK)
!!$          write(*,*)'State_VGB(OpRho_ ,i,j,k,globalBLK)=',&
!!$               State_VGB(OpRho_ ,i,j,k,globalBLK)
!!$          write(*,*)'State_VGB(CO2pRho_ ,i,j,k,globalBLK)=',&
!!$               State_VGB(CO2pRho_ ,i,j,k,globalBLK)
!!$          write(*,*)'State_VGB(Rho_ ,i,j,k,globalBLK)=',&
!!$               State_VGB(Rho_ ,i,j,k,globalBLK)
!!$          write(*,*)'O2pRho_dim', O2pRho_dim
!!$          write(*,*)'OpRho_dim', OpRho_dim
!!$          write(*,*)'CO2pRho_dim', CO2pRho_dim
!!$          write(*,*)' HpRho_dim ',  HpRho_dim


       end if !(true_cell?)

    end do; end do; end do

  end subroutine user_set_ICs

end module ModUser
