module ModThreadedLC
  use ModFieldLineThread, ONLY: &
       BoundaryThreads, BoundaryThreads_B
  use ModCoronalHeating, ONLY:PoyntingFluxPerBSi, PoyntingFluxPerB, &
                              LPerpTimesSqrtBSi, LPerpTimesSqrtB
  use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
  use ModPhysics,    ONLY: AverageIonCharge
  implicit none
  real :: TeFraction, TiFraction 
  real,allocatable :: Te_G(:,:,:)

contains
  !=========================================================================
  subroutine init_threaded_lc
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModMultifluid,   ONLY: MassIon_I
    use ModFieldLineThread, ONLY: check_tr_table
    allocate(Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)); Te_G = 0.0

    ! TeFraction is used for ideal EOS:
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TiFraction = MassIon_I(1)
       TeFraction = MassIon_I(1)/AverageIonCharge
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge)
       TeFraction = TiFraction
    end if
    call check_tr_table
  end subroutine init_threaded_lc
  !================================
  !\
  ! Main routine:
  ! solves MHD equations along thread, to link the state above the
  ! inner boundary of the global solar corona to the photosphere
  !/
  subroutine solve_boundary_thread(j, k, iBlock, TeSiIn, USiIn, AMinorIn,&
       DTeOverDsSiOut, PAvrSiOut, AMajorOut)
    !\
    ! USE:
    !/
    use ModFieldLineThread, ONLY: HeatCondParSi
    use ModPhysics,      ONLY: inv_gm1, No2Si_V, UnitX_,Si2No_V, &
                               UnitEnergyDens_, UnitTemperature_, UnitB_
    use ModLookupTable,  ONLY: i_lookup_table, interpolate_lookup_table
    use ModImplicit,   ONLY: StateSemi_VGB, iTeImpl 
    !\
    !The initial version: pressure is constant along the thread,
    !reflection and dissipation of the major wave is ignored
    !/
    !INPUT:
    !\
    !Cell and block indexes for the boundary point
    !/
    integer,intent(in):: j, k, iBlock 
    !\
    ! Parameters of the state in the true cell near the boundary:
    ! TeSiIn: Temperature in K 
    ! USiIn:  Velocity progection on the magnetic field direction.
    ! It is positive if the wind blows outward the Sun.
    ! AMinorIn: for the wave propagating toward the Sun 
    !            EnergyDensity = (\Pi/B)\sqrt{\rho} AMinor**2 
    !/
    real,   intent(in):: TeSiIn, USiIn, AMinorIn
    !\
    !OUTPUT:
    !DTeOverDsSiOut: Temperature derivative along the thread, at the end point
    !                Used to find the electron temperature in the ghostcell
    !PAvrSiOut: The geometric mean of electron and ion pressure (\sqrt{Pe*Pi})
    !AMajorOut: For the wave propagating outward the Sun
    !            EnergyDensity = (\Pi/B)\sqrt{\rho} AMajor**2  
    !/
    real,  intent(out):: DTeOverDsSiOut, PAvrSiOut, AMajorOut

    !\
    ! Two components arrays to use lookup table
    !/ 
    real    :: Value_V(2), AWValue_V(2), Length, RhoNoDim, Heating
    integer :: iTable, iTableAW
    !------------------!
    iTable = i_lookup_table('TR')
    if(iTable<=0)call CON_stop('TR table is not set')
    call interpolate_lookup_table(iTable, TeSiIn, 1.0e8, Value_V, &
         DoExtrapolate=.false.)
    !\
    ! First value is now the product of the thread length in meters times
    ! a geometric mean pressure, so that
    !/
    PAvrSiOut = Value_V(1)/( BoundaryThreads_B(iBlock)% Length_III(0,j,k) * &
         No2Si_V(UnitX_) )

    RhoNoDim = (PAvrSiOut*Si2No_V(UnitEnergyDens_)/sqrt(AverageIonCharge))*&
          TeFraction/(TeSiIn*Si2No_V(UnitTemperature_))

    !Dimmensionless length (related to the wave dissipation length)
    Length = BoundaryThreads_B(iBlock)% Length2SqrtB_III(0,j,k)*&
         sqrt(sqrt(RhoNoDim)*PoyntingFluxPerB/LperpTimesSqrtBSi**2)

    !Calculate Alfven waves for the given thread length and BC for ingoing wave 
    iTable = i_lookup_table('AW_TR')
    if(iTable<=0)call CON_stop('AW_TR table is not set')
    call interpolate_lookup_table(iTable, Length, AMinorIn, AWValue_V, &
         DoExtrapolate=.false.)

    Heating = AWValue_V(1)*PoyntingFluxPerBSi*&
         BoundaryThreads_B(iBlock)% Length_III(0,j,k)*No2Si_V(UnitB_)
    !\
    ! Heat flux equals PAvr * UHeat (spent for radiation) +
    ! Pi * U * (5/2) + Pe * U * (5/2) (carried away by outflow), 
    ! the pressure ratio being  Pe = Pi * AverageIonCharge
    ! 
    ! Temperature gradient equals the heat flux divided by kappa0*T**2.5
    !/
    DTeOverDsSiOut = ( PAvrSiOut * Value_V(2) - & !Radiation losses
         Heating                                & !AW Heating
         +USiIn * (PAvrSiOut/sqrt(AverageIonCharge)) *(inv_gm1 +1) * & !5/2*U*Pi
         (1 + AverageIonCharge) ) /&
         (HeatCondParSi * TeSiIn**2.50)

    AMajorOut = AWValue_V(2) 
  end subroutine solve_boundary_thread
  !=========================================================================
  subroutine set_field_line_thread_bc(nGhost, iBlock, nVarState, State_VG, &
               iImplBlock)
    use ModAdvance,      ONLY: State_VGB
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_, UnitU_, UnitX_, OmegaBody
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_, Bx_, Bz_, RhoUx_, RhoUz_
    use ModSize,         ONLY: nDim
    use ModConst,        ONLY: cTolerance
    use ModImplicit,     ONLY: iTeImpl
    use ModWaves
    use ModGeometry,     ONLY: Xyz_DGB
    use ModMain,         ONLY: jTest, kTest, BlkTest, UseRotatingFrame, x_, z_
    integer, intent(in):: nGhost
    integer, intent(in):: iBlock
    integer, intent(in):: nVarState
    real, intent(inout):: State_VG(nVarState,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Optional arguments when called by semi-implicit scheme
    integer, optional, intent(in):: iImplBlock
    
    logical:: IsNewBlock
    integer :: i, j, k, iP, Major_, Minor_
    real :: FaceGrad_D(3), TeSi, BDir_D(3), FaceCoord_D(3), U_D(3), B_D(3)
    real :: PAvrSI, U, AMinor, AMajor, DTeOverDsSi, GradTeDotB0Dir
    !-------------
    IsNewBlock = .true.

    !\
    ! Start from floating boundary values
    !/
    do k = MinK, MaxK; do j = MinJ, maxJ; do i = 1 - nGhost, 0
       State_VG(:, i,j,k) = State_VG(:,1, j, k)
    end do; end do; end do
    !\
    ! Fill in the temperature array
    !/
    if(present(iImplBlock))then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Te_G(i, j, k) = State_VG(iTeImpl,i,j,k)
       end do; end do; end do
    else
       if(UseMultiIon)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = State_VG(Pe_,i,j,k) &
                  /sum(ChargeIon_I*State_VG(iRhoIon_I,i,j,k)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          iP = p_
          if(UseElectronPressure) iP = Pe_
          
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction*State_VG(iP,i,j,k) &
                  /State_VG(Rho_,i,j,k)
          end do; end do; end do
       else
          call CON_stop('Generic EOS is not applicable with threads')
       end if
    end if
       
    do k = 1, nK; do j = 1, nJ
       !\
       ! Gradient across the boundary face
       !/
       call get_face_gradient(1, 1, j, k, iBlock, &
            IsNewBlock, Te_G, FaceGrad_D, &
            UseFirstOrderBcIn=.true.)
       !\
       !Store B1_D field in the physical cell
       !/
       B_D(x_:z_) = State_VGB(Bx_:Bz_, 1, j, k, iBlock)
       
       BDir_D = B_D + &
            BoundaryThreads_B(iBlock) % B0Face_DII(:, j, k)
       BDir_D = BDir_D/max(sqrt(sum(BDir_D**2)), cTolerance)
       if(BoundaryThreads_B(iBlock) % SignBr_II(j, k) < 0.0)then
          BDir_D = -BDir_D
          Major_ = WaveLast_
          Minor_ = WaveFirst_
       else
          Major_ = WaveFirst_
          Minor_ = WaveLast_
       end if
       !\
       ! Calculate input parameters for solving the thread
       !/
       TeSi = Te_G(1, j, k) * No2Si_V(UnitTemperature_)
       AMinor = sqrt(State_VGB(Minor_, 1, j, k, iBlock)/&
            ( sqrt(State_VGB(Rho_, 1, j, k, iBlock))* &
            PoyntingFluxPerB)  )
       U_D = State_VGB(RhoUx_:RhoUz_, 1, j, k, iBlock)/&
            State_VGB(Rho_, 1, j, k, iBlock)
       U = sum(U_D*BDir_D)
       call solve_boundary_thread(j=j, k=k, iBlock=iBlock, &
            TeSiIn=TeSi, USiIn=U*No2Si_V(UnitU_), AMinorIn=AMinor, &
            DTeOverDsSiOut=DTeOverDsSi, PAvrSiOut=PAvrSi, AMajorOut=AMajor) 
       !\
       ! Calculate temperature in the ghost cell by adding the difference 
       ! between the required value DTeOverDs and the temperature gradient
       ! calculated with the floating BC 
       !/ 
       Te_G(0, j, k) = min(max(Te_G(0, j, k) +(&
            DTeOverDsSi * Si2No_V(UnitTemperature_)/Si2No_V(UnitX_) - &
            sum(FaceGrad_D*BDir_D))/&
            sum(BoundaryThreads_B(iBlock)% DGradTeOverGhostTe_DII(:, j, k) &
            * BDir_D), 0.60*Te_G(0, j, k)),1.30*Te_G(0, j, k))
       if(present(iImplBlock))then
          State_VG(iTeImpl, 0, j, k) = Te_G(0, j, k)
          CYCLE
       end if
      
       if(UseElectronPressure)then
          State_VG(Pe_, 1-nGhost:0, j, k) = PAvrSi*Si2No_V(UnitEnergyDens_)*&
               sqrt(AverageIonCharge)
          State_VG(p_, 1-nGhost:0, j, k) = State_VG(Pe_, 1-nGhost:0, j, k)/&
               AverageIonCharge
       else
          State_VG(p_, 1-nGhost:0, j, k) = PAvrSi*Si2No_V(UnitEnergyDens_)/&
               sqrt(AverageIonCharge)*(1 + AverageIonCharge) 
       end if
       State_VG(Rho_, 0, j, k) = State_VG(iP, 0, j, k)* &
            TeFraction/Te_G(0, j, k)
       !Prolong the density gradient further way
       State_VG(Rho_, 1-nGhost:-1, j, k) = 2*State_VG(Rho_, 0, j, k) -&
            State_VG(iP, 0, j, k)* TeFraction/Te_G(1, j, k)

       !\
       !Calculate radial component of the B1 field and
       !maintain the radial component of the field to be 0
       !while the tangential ones being floating
       !/
       FaceCoord_D(1:nDim) = Xyz_DGB(1:nDim,1,j,k,iBlock)
       FaceCoord_D = FaceCoord_D/sqrt(sum(FaceCoord_D**2))
       B_D = B_D - FaceCoord_D*sum(FaceCoord_D*B_D)
 
       do i = 1-nGhost, 0
          State_VG(Bx_:Bz_, i, j, k) = B_D
          State_VG(RhoUx_:RhoUz_, i, j, k) = State_VG(Rho_,  i, j, k) * &
               U*BDir_D
          State_VG(Major_, i, j, k) = AMajor**2 * PoyntingFluxPerB *&
               sqrt( State_VG(Rho_, i, j, k) )
       end do
    end do; end do
  end subroutine set_field_line_thread_bc
end module ModThreadedLC
