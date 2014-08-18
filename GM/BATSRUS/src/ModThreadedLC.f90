module ModThreadedLC
  use ModFieldLineThread, ONLY: &
       BoundaryThreads, BoundaryThreads_B, HeatCondParSi, check_tr_table
  use ModCoronalHeating, ONLY:PoyntingFluxPerBSi, PoyntingFluxPerB, &
       UseWaveReflection, UseWaveReflection, TypeCoronalHeating
  use ModAdvance,    ONLY: UseElectronPressure, UseIdealEos
  implicit none
  real :: TeFraction, TiFraction 
  real,allocatable :: Te_G(:,:,:)

contains
  !=========================================================================
  subroutine init_threaded_lc
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModPhysics,      ONLY: ElectronTemperatureRatio, AverageIonCharge
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
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       TeFraction = TiFraction*ElectronTemperatureRatio
    end if
    call check_tr_table
  end subroutine init_threaded_lc
  !\
  ! Main routine:
  ! solves MHD equations along thread, to link the state above the
  ! inner boundary of the global solar corona to the photosphere
  !/
  subroutine solve_boundary_thread(j, k, iBlock, TeSiIn, USiIn, AMinorIn,&
       DTeOverDsSiOut, PAvrSiOut, AMajorOut)
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
 
  end subroutine solve_boundary_thread
  !=========================================================================
  subroutine set_field_line_thread_bc(nGhost, iBlock, nVarState, State_VG, &
               iImplBlock)
    use ModAdvance,      ONLY: State_VGB
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModFaceGradient, ONLY: get_face_gradient
    use ModPhysics,      ONLY: inv_gm1, Si2No_V, UnitTemperature_, &
         UnitEnergyDens_
    use ModMultifluid,   ONLY: UseMultiIon, MassIon_I, ChargeIon_I, iRhoIon_I
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Pe_, Bx_, Bz_, RhoUx_, RhoUz_
    use ModSize,         ONLY: nDim
    use ModConst,        ONLY: cTolerance
    use ModWaves
    use ModLookupTable, ONLY: i_lookup_table, interpolate_lookup_table
    integer, intent(in):: nGhost
    integer, intent(in):: iBlock
    integer, intent(in):: nVarState
    real, intent(inout):: State_VG(nVarState,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Optional arguments when called by semi-implicit scheme
    integer, optional, intent(in):: iImplBlock
    
    logical:: IsNewBlock
    integer :: i, j, k, iP, iTable, Major_, Minor_
    real :: FaceGrad_D(3), TeSi, BDir_D(3), Value_V(2)
    real :: PeSI, U_D(3), GradTeDotB0Dir 
    !---------------------
    IsNewBlock = .true.
    if(.not.present(iImplBlock))then
       if(UseMultiIon)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock) &
                  /sum(ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
          end do; end do; end do
       elseif(UseIdealEos)then
          iP = p_
          if(UseElectronPressure) iP = Pe_
       
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Te_G(i,j,k) = TeFraction*State_VGB(iP,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)
          end do; end do; end do
       else
          call CON_stop('Generic EOS is not applicable with threads')
       end if
       !\
       ! Start from floating boundary values
       !/
       do k = MinK, MaxK; do j = MinJ, MaxJ
          Te_G( MinI:0, j, k) = Te_G( MinI:0, j, k)
       end do; end do

       iTable = i_lookup_table('TR')
       if(iTable<=0)call CON_stop('TR table is not set')

       do k = 1, nK; do j = 1, nJ
          call get_face_gradient(1, 1, j, k, iBlock, &
               IsNewBlock, Te_G, FaceGrad_D, &
               UseFirstOrderBcIn=.true.)
          BDir_D = State_VGB(Bx_:Bz_, 1, j, k, iBlock) + &
               BoundaryThreads_B(iBlock) % B0Face_DII(:, j, k)
          BDir_D = BDir_D/sqrt(sum(BDir_D**2))* &
               BoundaryThreads_B(iBlock) % SignBr_II(j, k)
          
       end do; end do
    end if
  end subroutine set_field_line_thread_bc
end module ModThreadedLC
