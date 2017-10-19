!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModFaceFlux

  use ModProcMH,     ONLY: iProc
  use ModSize,       ONLY:x_, y_, z_, nI, nJ, nK, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModMain,       ONLY: UseB, UseB0, cLimit, &
       iTest, jTest, kTest, ProcTest, BlkTest, DimTest
  use ModMain,       ONLY: UseRadDiffusion, UseHeatConduction, &
       UseIonHeatConduction, DoThinCurrentSheet
  use ModMain,       ONLY: UseBorisSimple, UseBoris => boris_correction
  use ModGeometry,   ONLY: true_cell
  use BATL_lib,      ONLY: IsCartesianGrid, IsCartesian, IsRzGeometry, &
       Xyz_DGB, CellSize_DB, CellFace_DB, CellFace_DFB, FaceNormal_DDFB, &
       UseHighFDGeometry, correct_face_value

  use ModB0, ONLY: B0_DX, B0_DY, B0_DZ, B0_DGB ! input: face/cell centered B0

  use ModAdvance, ONLY:&
       LeftState_VX,  LeftState_VY,  LeftState_VZ,  & ! input: left  face state
       RightState_VX, RightState_VY, RightState_VZ, & ! input: right face state
       Flux_VX, Flux_VY, Flux_VZ,        & ! output: flux*Area
       VdtFace_x, VdtFace_y, VdtFace_z,  & ! output: cMax*Area for CFL
       EDotFA_X, EDotFA_Y, EDotFA_Z,     & ! output: E.Area
       uDotArea_XI, uDotArea_YI, uDotArea_ZI,& ! output: U.Area for P source
       bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ,& ! output: B x Area for J
       UseIdealEos, UseElectronPressure, &
       eFluid_, &                        ! index for electron fluid (nFluid+1)
       UseEfield, &                      ! electric field
       UseFDFaceFlux, FluxCenter_VGD, &
       UseLowOrderOnly, UseLowOrder, &
       UseLowOrder_X, UseLowOrder_Y, UseLowOrder_Z

  use ModPhysics, ONLY: ElectronPressureRatio, PePerPtotal

  use ModMultiIon, ONLY: &
       Pe_X, Pe_Y, Pe_Z, & ! output: Pe for grad Pe in multi-ion MHD
       Pwave_X, Pwave_Y, Pwave_Z ! Pwave for grad Pwave in multi-ion MHD

  use ModHallResist, ONLY: UseHallResist, HallCmaxFactor, IonMassPerCharge_G, &
       IsNewBlockCurrent, HallFactor_DF, set_hall_factor_face, &
       set_ion_mass_per_charge, UseBiermannBattery

  use ModRadDiffusion, ONLY: IsNewBlockRadDiffusion, get_radiation_energy_flux
  use ModHeatConduction, ONLY: IsNewBlockHeatCond, IsNewBlockIonHeatCond, &
       get_heat_flux, get_ion_heat_flux

  use ModResistivity, ONLY: UseResistiveFlux, Eta_GB
  use ModIonElectron, ONLY: iVarUseCmax_I
  use ModVarIndexes
  use ModMultiFluid
  use ModNumConst
  use ModCoronalHeating, ONLY: IsNewBlockAlfven
  use ModViscosity, ONLY: UseViscosity, IsNewBlockViscosity, Visco_DDI,&
       get_viscosity_tensor, set_visco_factor_face, ViscoFactor_DF

  implicit none

  ! Number of fluxes including pressure and energy fluxes
  integer, parameter :: nFlux=nVar+nFluid

  ! Range of fluid and variable indexes for the current solver
  integer:: iFluidMin = 1, iFluidMax = nFluid
  integer:: iVarMin   = 1, iVarMax   = nVar
  integer:: iEnergyMin= nVar+1, iEnergyMax = nVar + nFluid

  ! Neutral fluids may use different flux function
  logical:: UseDifferentNeutralFlux = .false.
  character(len=10):: TypeFluxNeutral = 'default'

  ! Logicals so we don't need string comparisons
  logical :: DoSimple, DoLf, DoHll, DoLfdw, DoHlldw, DoHlld, &
       DoAw, DoRoeOld, DoRoe
  logical :: DoLfNeutral, DoHllNeutral, DoHlldwNeutral, DoLfdwNeutral, &
       DoAwNeutral, DoGodunovNeutral, DoHllcNeutral

  ! 1D Burgers' equation, works for Hd equations. 
  logical:: DoBurgers = .false.

  logical :: UseLindeFix
  logical :: DoTestCell
  logical :: IsBoundary

  ! Index of the block for this face
  integer :: iBlockFace

  ! Direction of the face
  integer :: iDimFace

  ! index of the face 
  integer :: iFace, jFace, kFace

  ! index of cell in the negative and positive directions from face
  integer :: iLeft,  jLeft, kLeft
  integer :: iRight, jRight, kRight

  real :: StateLeft_V(nVar) = 0.0, StateRight_V(nVar) = 0.0
  real :: FluxLeft_V(nFlux) = 0.0, FluxRight_V(nFlux) = 0.0
  real :: StateLeftCons_V(nFlux) = 0.0, StateRightCons_V(nFlux) = 0.0
  real :: DissipationFlux_V(nFlux) = 0.0
  real :: B0x = 0.0, B0y = 0.0, B0z = 0.0
  real :: DiffBb = 0.0 !     (1/4)(BnL-BnR)^2
  real :: DeltaBnL = 0.0, DeltaBnR = 0.0
  real :: Area = 0.0, Area2 = 0.0, AreaX, AreaY, AreaZ

  ! Allow diffusion across the pole
  logical:: UsePoleDiffusion = .false.

  ! Maximum speed for the Courant condition
  real :: CmaxDt = 0.0

  ! Normal velocities for all fluids plus electrons
  real :: Unormal_I(nFluid+1) = 0.0
  real :: UnLeft_I(nFluid+1)  = 0.0
  real :: UnRight_I(nFluid+1) = 0.0

  real :: bCrossArea_D(3) = (/ 0.0, 0.0, 0.0 /) !B x Area for current -> BxJ
  real :: Enormal = 0.0                         !normal electric field -> div E
  real :: Pe      = 0.0                         !electron pressure -> grad Pe
  real :: Pwave   = 0.0

  ! Variables for normal resistivity
  real :: EtaJx, EtaJy, EtaJz, Eta = 0.0

  ! Variables needed for Hall resistivity
  real :: InvDxyz, HallCoeff, HallJx, HallJy, HallJz

  ! Variables needed for Biermann battery term
  logical :: UseHallGradPe = .false., IsNewBlockGradPe = .true.
  real :: BiermannCoeff, GradXPeNe, GradYPeNe, GradZPeNe
  real, allocatable, save :: Pe_G(:,:,:)

  ! Variables needed by viscosity
  real :: ViscoCoeff = 0.0

  ! Variables for diffusion solvers (radiation diffusion, heat conduction)
  real :: DiffCoef = 0.0, EradFlux, RadDiffCoef, HeatFlux, IonHeatFlux, &
       HeatCondCoefNormal

  ! These are variables for pure MHD solvers (Roe and HLLD)
  ! Number of MHD fluxes including the pressure and energy fluxes
  integer, parameter :: nFluxMhd = 9

  ! Named conservative MHD variable indexes + pressure
  integer, parameter :: RhoMhd_=1, RhoUn_=2, RhoUt1_=3, RhoUt2_=4, &
       B1n_=5, B1t1_=6, B1t2_=7, eMhd_=8, pMhd_=9

  ! Variables for rotated coordinate system (n is normal to face)
  real :: Normal_D(3), NormalX, NormalY, NormalZ
  real :: Tangent1_D(3), Tangent2_D(3)   
  real :: B0n, B0t1, B0t2
  real :: UnL, Ut1L, Ut2L, B1nL, B1t1L, B1t2L
  real :: UnR, Ut1R, Ut2R, B1nR, B1t1R, B1t2R

  ! Limit propagation speeds to reduce numerical diffusion
  logical :: UseClimit = .false.
  real    :: ClimitDim = -1.0, rClimit = -1.0

  ! One of the two possible ways to treat the MHD-like systems
  ! (partially symmetrizable, following the Godunov definition).
  ! If UseRS7=.true. then use the 7-wave Riemann Solver (RS) with 
  ! continuous  normal component of the magnetic field across the face.
  ! The number of jumps in the physical variables across the face is equal
  ! to the number of waves, resulting in a well-posed Riemann problem. 
  ! This approach is an alternative to the 8-wave scheme.
  logical:: UseRS7 = .false.


  ! Local logical variables for various terms that may be switched off
  ! if they are treated with semi-implicit scheme.
  logical :: DoRadDiffusion = .false., DoHeatConduction = .false., &
       DoIonHeatConduction = .false., DoHallInduction = .false.

  logical :: DoClightWarning     = .true.
  real    :: FactorClightWarning = 2.0

  character(len=*), private, parameter :: NameMod="ModFaceFlux"

contains

  !=========================================================================
  subroutine face_flux_set_parameters(NameCommand)
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand
    !------------------------------------------------------------------------
    select case(NameCommand)
    case('#CLIMIT')
       call read_var('UseClimit', UseClimit)
       if(UseClimit)then
          call read_var('ClimitDim', ClimitDim)
          call read_var('rClimit',   rClimit)
       else
          ! Make sure Climit is negative (if it was set in a previous session)
          Climit = -1.0
       end if
    case('#CLIGHTWARNING')
       call read_var('DoClightWarning', DoClightWarning)
       if (DoClightWarning) &
            call read_var('FactorClightWarning', FactorClightWarning)
    end select

  end subroutine face_flux_set_parameters

  !=========================================================================
  subroutine rotate_state_vectors

    use ModCoordTransform, ONLY: cross_product

    ! Rotate the vector variables B0*, StateLeft_V(B*_), StateLeft_V(U*_)
    ! StateRight_V(B*_), StateRight_V(U*_) into normal and
    ! tangential components with respect to the face.
    ! Store the rotated vector components in scalar variables UnL, Ut1L, ....
    ! Also store the transformation for rotating back.
    ! Current implementation is for a single ion fluid.

    if(IsCartesianGrid)then
       select case (iDimFace)
       case (x_) ! x face
          ! B0 on the face
          B0n  = B0x
          B0t1 = B0y
          B0t2 = B0z
          ! Left face
          UnL   =  StateLeft_V(Ux_)
          Ut1L  =  StateLeft_V(Uy_)
          Ut2L  =  StateLeft_V(Uz_)
          B1nL  =  StateLeft_V(Bx_)
          B1t1L =  StateLeft_V(By_)
          B1t2L =  StateLeft_V(Bz_)
          ! Right face
          UnR   =  StateRight_V(Ux_)
          Ut1R  =  StateRight_V(Uy_)
          Ut2R  =  StateRight_V(Uz_)
          B1nR  =  StateRight_V(Bx_)
          B1t1R =  StateRight_V(By_)
          B1t2R =  StateRight_V(Bz_)
       case (y_) ! y face
          ! B0 on the face
          B0n  = B0y
          B0t1 = B0z
          B0t2 = B0x
          ! Left face
          UnL   =  StateLeft_V(Uy_)
          Ut1L  =  StateLeft_V(Uz_)
          Ut2L  =  StateLeft_V(Ux_)
          B1nL  =  StateLeft_V(By_)
          B1t1L =  StateLeft_V(Bz_)
          B1t2L =  StateLeft_V(Bx_)
          ! Right face
          UnR   =  StateRight_V(Uy_)
          Ut1R  =  StateRight_V(Uz_)
          Ut2R  =  StateRight_V(Ux_)
          B1nR  =  StateRight_V(By_)
          B1t1R =  StateRight_V(Bz_)
          B1t2R =  StateRight_V(Bx_)
       case (z_) ! z face
          ! B0 on the face
          B0n  = B0z
          B0t1 = B0x
          B0t2 = B0y
          ! Left face
          UnL   =  StateLeft_V(Uz_)
          Ut1L  =  StateLeft_V(Ux_)
          Ut2L  =  StateLeft_V(Uy_)
          B1nL  =  StateLeft_V(Bz_)
          B1t1L =  StateLeft_V(Bx_)
          B1t2L =  StateLeft_V(By_)
          ! Right face
          UnR   =  StateRight_V(Uz_)
          Ut1R  =  StateRight_V(Ux_)
          Ut2R  =  StateRight_V(Uy_)
          B1nR  =  StateRight_V(Bz_)
          B1t1R =  StateRight_V(Bx_)
          B1t2R =  StateRight_V(By_)
       end select
    else
       if(Normal_D(z_) < 0.5)then
          ! Tangent1 = Normal x (0,0,1)
          Tangent1_D(x_) =  Normal_D(y_)
          Tangent1_D(y_) = -Normal_D(x_)
          Tangent1_D(z_) = 0.0
       else
          ! Tangent1 = Normal x (1,0,0)
          Tangent1_D(x_) = 0.0
          Tangent1_D(y_) =  Normal_D(z_)
          Tangent1_D(z_) = -Normal_D(y_)
       end if
       ! Normalize Tangent1 vector
       Tangent1_D = Tangent1_D/sqrt(sum(Tangent1_D**2))
       ! Tangent2 = Normal x Tangent1
       Tangent2_D = cross_product(Normal_D, Tangent1_D)

       ! B0 on the face
       B0n   = sum(Normal_D  *(/B0x, B0y, B0z/))
       B0t1  = sum(Tangent1_D*(/B0x, B0y, B0z/))
       B0t2  = sum(Tangent2_D*(/B0x, B0y, B0z/))
       ! Left face
       UnL   = sum(Normal_D  *StateLeft_V(Ux_:Uz_))
       Ut1L  = sum(Tangent1_D*StateLeft_V(Ux_:Uz_))
       Ut2L  = sum(Tangent2_D*StateLeft_V(Ux_:Uz_))
       B1nL  = sum(Normal_D  *StateLeft_V(Bx_:Bz_))
       B1t1L = sum(Tangent1_D*StateLeft_V(Bx_:Bz_))
       B1t2L = sum(Tangent2_D*StateLeft_V(Bx_:Bz_))
       ! Right face
       UnR   = sum(Normal_D  *StateRight_V(Ux_:Uz_))
       Ut1R  = sum(Tangent1_D*StateRight_V(Ux_:Uz_))
       Ut2R  = sum(Tangent2_D*StateRight_V(Ux_:Uz_))
       B1nR  = sum(Normal_D  *StateRight_V(Bx_:Bz_))
       B1t1R = sum(Tangent1_D*StateRight_V(Bx_:Bz_))
       B1t2R = sum(Tangent2_D*StateRight_V(Bx_:Bz_))
    end if

  end subroutine rotate_state_vectors
  !==========================================================================
  subroutine rotate_flux_vector(FluxRot_V, Flux_V)
    real, intent(in)   :: FluxRot_V(:)
    real, intent(inout):: Flux_V(:)

    ! Rotate n,t1,t2 components back to x,y,z components
    if(IsCartesianGrid)then                             
       select case (iDimFace)
       case (x_)
          Flux_V(RhoUx_ ) = FluxRot_V(RhoUn_)
          Flux_V(RhoUy_ ) = FluxRot_V(RhoUt1_)
          Flux_V(RhoUz_ ) = FluxRot_V(RhoUt2_)
          Flux_V(Bx_    ) = FluxRot_V(B1n_)
          Flux_V(By_    ) = FluxRot_V(B1t1_)
          Flux_V(Bz_    ) = FluxRot_V(B1t2_)
       case (y_)
          Flux_V(RhoUx_ ) = FluxRot_V(RhoUt2_)
          Flux_V(RhoUy_ ) = FluxRot_V(RhoUn_)
          Flux_V(RhoUz_ ) = FluxRot_V(RhoUt1_)
          Flux_V(Bx_    ) = FluxRot_V(B1t2_)
          Flux_V(By_    ) = FluxRot_V(B1n_)
          Flux_V(Bz_    ) = FluxRot_V(B1t1_)
       case (z_)
          Flux_V(RhoUx_ ) = FluxRot_V(RhoUt1_)
          Flux_V(RhoUy_ ) = FluxRot_V(RhoUt2_)
          Flux_V(RhoUz_ ) = FluxRot_V(RhoUn_)
          Flux_V(Bx_    ) = FluxRot_V(B1t1_)
          Flux_V(By_    ) = FluxRot_V(B1t2_)
          Flux_V(Bz_    ) = FluxRot_V(B1n_)
       end select
    else
       Flux_V(RhoUx_) = Normal_D(x_)  *FluxRot_V(RhoUn_)  &
            +           Tangent1_D(x_)*FluxRot_V(RhoUt1_) &
            +           Tangent2_D(x_)*FluxRot_V(RhoUt2_)
       Flux_V(RhoUy_) = Normal_D(y_)  *FluxRot_V(RhoUn_)  &
            +           Tangent1_D(y_)*FluxRot_V(RhoUt1_) &
            +           Tangent2_D(y_)*FluxRot_V(RhoUt2_)
       Flux_V(RhoUz_) = Normal_D(z_)  *FluxRot_V(RhoUn_)  &
            +           Tangent1_D(z_)*FluxRot_V(RhoUt1_) &
            +           Tangent2_D(z_)*FluxRot_V(RhoUt2_)

       Flux_V(Bx_   ) = Normal_D(x_)  *FluxRot_V(B1n_)  &
            +           Tangent1_D(x_)*FluxRot_V(B1t1_) &
            +           Tangent2_D(x_)*FluxRot_V(B1t2_)
       Flux_V(By_   ) = Normal_D(y_)  *FluxRot_V(B1n_)  &
            +           Tangent1_D(y_)*FluxRot_V(B1t1_) &
            +           Tangent2_D(y_)*FluxRot_V(B1t2_)
       Flux_V(Bz_   ) = Normal_D(z_)  *FluxRot_V(B1n_)  &
            +           Tangent1_D(z_)*FluxRot_V(B1t1_) &
            +           Tangent2_D(z_)*FluxRot_V(B1t2_)
    end if

  end subroutine rotate_flux_vector
  !===========================================================================
  subroutine calc_face_flux(DoResChangeOnly, iBlock)

    use ModAdvance,  ONLY: TypeFlux => FluxType
    use ModParallel, ONLY: &
         neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth
    use ModMain,     ONLY: nIFace, nJFace, nKFace, &
         iMinFace, iMaxFace, jMinFace, jMaxFace, kMinFace, kMaxFace, &
         UseHyperbolicDivb
    use ModImplicit, ONLY: TypeSemiImplicit, UseSemiHallResist
    use ModWaves,    ONLY: UseWavePressure
    use ModViscosity,ONLY: UseArtificialVisco, AlphaVisco

    logical, intent(in) :: DoResChangeOnly
    integer, intent(in) :: iBlock
    logical :: DoTest, DoTestMe


    character(len=*), parameter:: NameSub = 'calc_face_flux'
    !--------------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if(DoTestMe)call print_values

    DoSimple = TypeFlux == 'Simple'
    DoLf     = TypeFlux == 'Rusanov'
    DoHll    = TypeFlux == 'Linde'
    DoLfdw   = TypeFlux == 'LFDW'
    DoHlldw  = TypeFlux == 'HLLDW'
    DoHlld   = TypeFlux == 'HLLD'
    DoAw     = TypeFlux == 'Sokolov'
    DoRoeOld = TypeFlux == 'RoeOld'
    DoRoe    = TypeFlux == 'Roe'

    DoLfNeutral      = TypeFluxNeutral == 'Rusanov'
    DoHllNeutral     = TypeFluxNeutral == 'Linde'
    DoLfdwNeutral    = TypeFluxNeutral == 'LFDW'
    DoHlldwNeutral   = TypeFluxNeutral == 'HLLDW'
    DoAwNeutral      = TypeFluxNeutral == 'Sokolov'
    DoGodunovNeutral = TypeFluxNeutral == 'Godunov'
    DoHllcNeutral    = TypeFluxNeutral == 'HLLC'

    UseRS7 = DoRoe  ! This is always true for the current implementation

    UseLindeFix = UseB .and. &
         (UseHyperbolicDivb .or. DoHll .or. DoLfdw .or. DoHlldw &
         .or. DoHlld .or. DoAw)

    ! Make sure that Hall MHD recalculates the magnetic field 
    ! in the current block that will be used for the Hall term
    IsNewBlockCurrent      = .true.
    IsNewBlockGradPe       = .true.
    IsNewBlockRadDiffusion = .true.
    IsNewBlockHeatCond     = .true.
    IsNewBlockIonHeatCond  = .true.
    IsNewBlockViscosity    = .true.
    IsNewBlockAlfven       = .true.

    DoRadDiffusion      = UseRadDiffusion .and.&
         .not. (index(TypeSemiImplicit,'radiation')>0 .or.&
         index(TypeSemiImplicit,'radcond')>0 .or.&
         index(TypeSemiImplicit,'cond')>0)

    DoHeatConduction    = UseHeatConduction .and.&
         .not.  index(TypeSemiImplicit,'parcond')>0

    DoHallInduction = UseHallResist .and. .not. UseSemiHallResist

    if(UseHallResist)then
       call set_hall_factor_face(iBlock)
    elseif(UseBiermannBattery)then
       call set_ion_mass_per_charge(iBlock)
    end if

    if(UseViscosity) call set_visco_factor_face(iBlock)

    if(UseFDFaceFlux) call calc_simple_cell_flux(iBlock)
    if (DoResChangeOnly) then
       if(neiLeast(iBlock) == 1) &
            call get_flux_x(1,1,1,nJ,1,nK)
       if(neiLwest(iBlock) == 1) &
            call get_flux_x(nIFace,nIFace,1,nJ,1,nK)
       if(nJ > 1 .and. neiLsouth(iBlock) == 1) &
            call get_flux_y(1,nI,1,1,1,nK)
       if(nJ > 1 .and. neiLnorth(iBlock) == 1) &
            call get_flux_y(1,nI,nJFace,nJFace,1,nK)
       if(nK > 1 .and. neiLbot(iBlock)   == 1) &
            call get_flux_z(1,nI,1,nJ,1,1)
       if(nK > 1 .and. neiLtop(iBlock)   == 1) &
            call get_flux_z(1,nI,1,nJ,nKFace,nKFace)
    else
       call get_flux_x(1, nIFace, jMinFace, jMaxFace, kMinFace, kMaxFace)
       if(nJ > 1) &
            call get_flux_y(iMinFace, iMaxFace, 1, nJFace, kMinFace ,kMaxFace)
       if(nK > 1) &
            call get_flux_z(iMinFace, iMaxFace, jMinFace, jMaxFace, 1, nKFace)
    end if

  contains
    !=========================================================================
    subroutine print_values
      integer :: iVar
      !---------------------------------------------------------------------
      if(DoResChangeOnly)then
         write(*,*)'calc_facefluxes for DoResChangeOnly'
         RETURN
      end if

      if(DimTest==x_ .or. DimTest==0)then
         write(*,*)&
              'Calc_facefluxes, left and right states at i-1/2 and i+1/2:'

         do iVar=1,nVar
            write(*,'(2a,4(1pe13.5))')NameVar_V(iVar),'=',&
                 LeftState_VX(iVar,iTest,jTest,kTest),&
                 RightState_VX(iVar,iTest,  jTest,kTest),&
                 LeftState_VX(iVar,iTest+1,jTest,kTest),&
                 RightState_VX(iVar,iTest+1,jTest,kTest)
         end do
         if(UseB0)then
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
                 B0_DX(x_,iTest,jTest,kTest),' ',&
                 B0_DX(x_,iTest+1,jTest,kTest)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
                 B0_DX(y_,iTest,jTest,kTest),' ',&
                 B0_DX(y_,iTest+1,jTest,kTest)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
                 B0_DX(z_,iTest,jTest,kTest),' ',&
                 B0_DX(z_,iTest+1,jTest,kTest)
         end if
      end if

      if(DimTest==y_ .or. DimTest==0)then
         write(*,*)&
              'Calc_facefluxes, left and right states at j-1/2 and j+1/2:'

         do iVar=1,nVar
            write(*,'(2a,4(1pe13.5))')NameVar_V(iVar),'=',&
                 LeftState_VY(iVar,iTest,jTest,kTest),&
                 RightState_VY(iVar,iTest,  jTest,kTest),&
                 LeftState_VY(iVar,iTest,jTest+1,kTest),&
                 RightState_VY(iVar,iTest,jTest+1,kTest)
         end do
         if(UseB0)then
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
                 B0_DY(x_,iTest,jTest,kTest),' ',&
                 B0_DY(x_,iTest,jTest+1,kTest)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
                 B0_DY(y_,iTest,jTest,kTest),' ',&
                 B0_DY(y_,iTest,jTest+1,kTest)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
                 B0_DY(z_,iTest,jTest,kTest),' ',&
                 B0_DY(z_,iTest,jTest+1,kTest)
         end if
      end if

      if(DimTest==z_ .or. DimTest==0)then
         write(*,*)&
              'Calc_facefluxes, left and right states at k-1/2 and k+1/2:'
         do iVar=1,nVar
            write(*,'(2a,4(1pe13.5))')NameVar_V(iVar),'=',&
                 LeftState_VZ(iVar,iTest,jTest,kTest),&
                 RightState_VZ(iVar,iTest,  jTest,kTest),&
                 LeftState_VZ(iVar,iTest,jTest,kTest+1),&
                 RightState_VZ(iVar,iTest,jTest,kTest+1)
         end do
         if(UseB0)then
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0x:',&
                 B0_DZ(x_,iTest,jTest,kTest),' ',&
                 B0_DZ(x_,iTest,jTest,kTest+1)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0y:',&
                 B0_DZ(y_,iTest,jTest,kTest),' ',&
                 B0_DZ(y_,iTest,jTest,kTest+1)
            write(*,'(a,1pe13.5,a13,1pe13.5)')'B0z:',&
                 B0_DZ(z_,iTest,jTest,kTest),' ',&
                 B0_DZ(z_,iTest,jTest,kTest+1)
         end if
      end if

    end subroutine print_values

    !==========================================================================

    subroutine get_flux_x(iMin,iMax,jMin,jMax,kMin,kMax)

      use ModAdvance, ONLY: State_VGB
      integer, intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer:: iFlux
      !-----------------------------------------------------------------------

      call set_block_values(iBlock, x_)

      do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax

         DoTestCell = DoTestMe &
              .and. (iFace == iTest .or. iFace == iTest+1) &
              .and. jFace == jTest .and. kFace == kTest

         call set_cell_values_x

         if(  .not. true_cell(iLeft,jLeft,kLeft,iBlock) .and. &
              .not. true_cell(iRight,jRight,kRight,iBlock)) then
            uDotArea_XI(iFace,jFace,kFace,:) = 0.0
            VdtFace_X(iFace,jFace,kFace) = 0.0
            CYCLE
         endif

         if(UseB0)then
            B0x = B0_DX(x_,iFace,jFace,kFace)
            B0y = B0_DX(y_,iFace,jFace,kFace)
            B0z = B0_DX(z_,iFace,jFace,kFace)
         end if

         if(UseRS7.and..not.IsBoundary)then
            DeltaBnR=sum((RightState_VX(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace,jFace,kFace,iBlockFace))*&
                 Normal_D)
            RightState_VX(Bx_:Bz_, iFace, jFace, kFace)=&
                 RightState_VX(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnR* Normal_D   
            DeltaBnL=sum((LeftState_VX(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace-1,jFace,kFace,iBlockFace))*&
                 Normal_D)
            LeftState_VX(Bx_:Bz_, iFace, jFace, kFace)=&
                 LeftState_VX(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnL* Normal_D         
         else
            DeltaBnL=cZero;DeltaBnR=cZero
         end if
         StateLeft_V  = LeftState_VX( :, iFace, jFace, kFace)
         StateRight_V = RightState_VX(:, iFace, jFace, kFace)
         
         call get_numerical_flux(Flux_VX(:,iFace, jFace, kFace))
         
         if(UseArtificialVisco) call add_artificial_viscosity( &
              Flux_VX(:,iFace, jFace, kFace))

         VdtFace_x(iFace, jFace, kFace)       = CmaxDt*Area

         ! Correct Unormal_I to make div(u) achieve 6th order.
         if(UseFDFaceFlux) call correct_u_normal(x_) 
         uDotArea_XI(iFace, jFace, kFace,:)   = Unormal_I*Area

         if(UseB .and. UseBoris) &
              EDotFA_X(iFace,jFace,kFace) = Enormal*Area

         if(UseB .and. UseMultiIon)then
            Pe_X(iFace, jFace, kFace) = Pe
            if(UseWavePressure) Pwave_X(iFace,jFace,kFace) = Pwave
         end if

         if(UseB .and. (UseMultiIon .or. .not.IsMhd)) &
              bCrossArea_DX(:,iFace,jFace,kFace) = bCrossArea_D
      end do; end do; end do

      if(UseFDFaceFlux .and. .not.UseLowOrderOnly) then
         ! For FD method, modify flux so that df/dx=(f(j+1/2)-f(j-1/2))/dx 
         ! is 6th order. 
         do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax
            if(UseLowOrder)then
               if(UseLowOrder_X(iFace,jFace,kFace)) CYCLE
            end if
            do iFlux = 1, nFlux
               Flux_VX(iFlux,iFace,jFace,kFace) = &
                    correct_face_value(Flux_VX(iFlux,iFace,jFace,kFace),&
                    FluxCenter_VGD(iFlux,iFace-2:iFace+1,jFace,kFace,1))
            enddo
         end do; end do; enddo
      endif
    end subroutine get_flux_x

    !==========================================================================

    subroutine get_flux_y(iMin,iMax,jMin,jMax,kMin,kMax)

      use ModAdvance, ONLY: State_VGB
      integer, intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer:: iFlux
      !------------------------------------------------------------------------

      call set_block_values(iBlock, y_)

      do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax

         DoTestCell = DoTestMe .and. iFace == iTest .and. &
              (jFace == jTest .or. jFace == jTest+1) .and. kFace == kTest

         call set_cell_values_y

         if(  .not. true_cell(iLeft,jLeft,kLeft,iBlock) .and. &
              .not. true_cell(iRight,jRight,kRight,iBlock)) then
            uDotArea_YI(iFace,jFace,kFace,:) = 0.0
            VdtFace_Y(iFace,jFace,kFace) = 0.0
            CYCLE
         endif
         
         if(UseB0)then
            B0x = B0_DY(x_,iFace, jFace, kFace)
            B0y = B0_DY(y_,iFace, jFace, kFace)
            B0z = B0_DY(z_,iFace, jFace, kFace)
         end if

         if(UseRS7.and..not.IsBoundary)then
            DeltaBnR=sum((RightState_VY(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace,jFace,kFace,iBlockFace))*&
                 Normal_D)
            RightState_VY(Bx_:Bz_, iFace, jFace, kFace)=&
                 RightState_VY(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnR* Normal_D
            DeltaBnL=sum((LeftState_VY(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace,jFace-1,kFace,iBlockFace))*&
                 Normal_D)
            LeftState_VY(Bx_:Bz_, iFace, jFace, kFace)=&
                 LeftState_VY(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnL* Normal_D         
         else
            DeltaBnL=cZero;DeltaBnR=cZero
         end if

         StateLeft_V  = LeftState_VY( :, iFace, jFace, kFace)
         StateRight_V = RightState_VY(:, iFace, jFace, kFace)

         call get_numerical_flux(Flux_VY(:, iFace, jFace, kFace))

         if(UseArtificialVisco) call add_artificial_viscosity( &
              Flux_VY(:,iFace, jFace, kFace))
         
         VdtFace_y(iFace, jFace, kFace)       = CmaxDt*Area

         if(UseFDFaceFlux) call correct_u_normal(y_)
         uDotArea_YI(iFace, jFace, kFace, :)  = Unormal_I*Area

         if(UseB .and. UseBoris) &
              EDotFA_Y(iFace,jFace,kFace) = Enormal*Area

         if(UseB .and. UseMultiIon)then
            Pe_Y(iFace,jFace,kFace) = Pe
            if(UseWavePressure) Pwave_Y(iFace,jFace,kFace) = Pwave
         end if

         if(UseB .and. (UseMultiIon .or. .not.IsMhd)) &
              bCrossArea_DY(:,iFace,jFace,kFace) = bCrossArea_D

      end do; end do; end do

      !For FD method, modify flux so that df/dx=(f(j+1/2)-f(j-1/2))/dx (x=xj)
      !is 6th order. 
      if(UseFDFaceFlux .and. .not.UseLowOrderOnly) then
         do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax
            if(UseLowOrder)then
               if(UseLowOrder_Y(iFace,jFace,kFace)) CYCLE
            end if
            do iFlux = 1, nFlux
               Flux_VY(iFlux,iFace,jFace,kFace) = &
                    correct_face_value(&
                    Flux_VY(iFlux,iFace,jFace,kFace),&
                    FluxCenter_VGD(iFlux,iFace,jFace-2:jFace+1,kFace,2))
            enddo
         end do; end do; enddo
      end if
    end subroutine get_flux_y

    !==========================================================================

    subroutine get_flux_z(iMin, iMax, jMin, jMax, kMin, kMax)

      use ModAdvance, ONLY: State_VGB
      integer, intent(in):: iMin, iMax, jMin, jMax, kMin, kMax
      integer:: iFlux
      !------------------------------------------------------------------------

      call set_block_values(iBlock, z_)

      do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax

         DoTestCell = DoTestMe .and. iFace == iTest .and. &
              jFace == jTest .and. (kFace == kTest .or. kFace == kTest+1)

         call set_cell_values_z

         if(  .not. true_cell(iLeft,jLeft,kLeft,iBlock) .and. &
              .not. true_cell(iRight,jRight,kRight,iBlock)) then
            uDotArea_ZI(iFace,jFace,kFace,:) = 0.0
            VdtFace_Z(iFace,jFace,kFace) = 0.0
            CYCLE
         endif
         
         if(UseB0)then
            B0x = B0_DZ(x_,iFace, jFace, kFace)
            B0y = B0_DZ(y_,iFace, jFace, kFace)
            B0z = B0_DZ(z_,iFace, jFace, kFace)
         end if
         if(UseRS7.and..not.IsBoundary)then
            DeltaBnR=sum((RightState_VZ(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace,jFace,kFace,iBlockFace))*&
                 Normal_D)
            RightState_VZ(Bx_:Bz_, iFace, jFace, kFace)=&
                 RightState_VZ(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnR* Normal_D   
            DeltaBnL=sum((LeftState_VZ(Bx_:Bz_, iFace, jFace, kFace)-&
                 State_VGB(Bx_:Bz_,iFace,jFace,kFace-1,iBlockFace))*&
                 Normal_D)
            LeftState_VZ(Bx_:Bz_, iFace, jFace, kFace)=&
                 LeftState_VZ(Bx_:Bz_, iFace, jFace, kFace)-&
                 DeltaBnL* Normal_D         
         else
            DeltaBnL=cZero;DeltaBnR=cZero
         end if

         StateLeft_V  = LeftState_VZ( :, iFace, jFace, kFace)
         StateRight_V = RightState_VZ(:, iFace, jFace, kFace)

         call get_numerical_flux(Flux_VZ(:, iFace, jFace, kFace))

         if(UseArtificialVisco) call add_artificial_viscosity( &
              Flux_VZ(:,iFace, jFace, kFace))

         VdtFace_z(iFace, jFace, kFace)       = CmaxDt*Area

         if(UseFDFaceFlux) call correct_u_normal(z_)  
         uDotArea_ZI(iFace, jFace, kFace, :)  = Unormal_I*Area

         if(UseB .and. UseBoris) &
              EDotFA_Z(iFace,jFace,kFace) = Enormal*Area

         if(UseB .and. UseMultiIon)then
            Pe_Z(iFace,jFace,kFace) = Pe
            if(UseWavePressure) Pwave_Z(iFace,jFace,kFace) = Pwave
         end if

         if(UseB .and. (UseMultiIon .or. .not.IsMhd)) &
              bCrossArea_DZ(:,iFace,jFace,kFace)= bCrossArea_D

      end do; end do; end do

      if(UseFDFaceFlux .and. .not.UseLowOrderOnly) then
         do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax
            if(UseLowOrder)then
               if(UseLowOrder_Z(iFace,jFace,kFace)) CYCLE
            end if
            do iFlux = 1, nFlux
               Flux_VZ(iFlux,iFace,jFace,kFace) = &
                    correct_face_value(Flux_VZ(iFlux,iFace,jFace,kFace),&
                    FluxCenter_VGD(iFlux,iFace,jFace,kFace-2:kFace+1,3))
            enddo
         end do; end do; enddo
      end if
    end subroutine get_flux_z
    !=====================================================================

    subroutine add_artificial_viscosity(Flux_V)
      use ModAdvance,ONLY: State_VGB, Energy_GBI
      real,     intent(inout) :: Flux_V(nFlux)
      
      if(all(true_cell(iLeft:iFace,jLeft:jFace,kLeft:kFace,iBlockFace))) then

         Flux_V(1:nVar) = Flux_V(1:nVar) - AlphaVisco*CmaxDt* &
              (State_VGB(1:nVar,iFace,jFace,kFace,iBlockFace) - &
              State_VGB(1:nVar,iLeft,jLeft,kLeft,iBlockFace))* &
              Area 

         Flux_V(nVar+1:nFlux) = Flux_V(nVar+1:nFlux) - AlphaVisco*CmaxDt* &
              (Energy_GBI(iFace,jFace,kFace,iBlockFace,:) - &
              Energy_GBI(iLeft,jLeft,kLeft,iBlockFace,:))* &
              Area 
      endif
    end subroutine add_artificial_viscosity

  end subroutine calc_face_flux

  !===========================================================================
  subroutine set_block_values(iBlock, iDim)
    integer, intent(in) :: iBlock, iDim

    iBlockFace = iBlock
    iDimFace   = iDim

    ! For generalized coordinates the values below are calculated cell by cell
    if(.not.IsCartesian .and. .not. UseHighFDGeometry) RETURN    

    ! Calculate face normal and area vectors for Cartesian grid
    Normal_D = 0.0; Normal_D(iDim) = 1.0
    NormalX = Normal_D(x_); NormalY = Normal_D(y_); NormalZ = Normal_D(z_)

    Area = CellFace_DB(iDim,iBlockFace)
    InvDxyz = 1./CellSize_DB(iDim,iBlockFace)
    select case(iDim)
    case(x_)
       AreaX   = Area; AreaY = 0.0; AreaZ = 0.0
    case(y_)
       AreaY   = Area; AreaX = 0.0; AreaZ = 0.0
    case(z_)
       AreaZ   = Area; AreaX = 0.0; AreaY = 0.0
    end select
    Area2 = Area**2

  end subroutine set_block_values
  !===========================================================================
  subroutine set_cell_values

    select case(iDimFace)
    case(x_)
       call set_cell_values_x
    case(y_)
       call set_cell_values_y
    case(z_)
       call set_cell_values_z
    end select

  end subroutine set_cell_values
  !===========================================================================
  subroutine set_cell_values_x

    iLeft = iFace - 1; jLeft = jFace; kLeft = kFace

    if(.not.IsCartesian)then

       if(IsRzGeometry)then
          AreaX = CellFace_DFB(1,iFace,jFace,kFace,iBlockFace)
          AreaY = 0.0
          AreaZ = 0.0
       else
          AreaX = FaceNormal_DDFB(x_, 1, iFace, jFace, kFace, iBlockFace)
          if(nJ > 1)then
             AreaY = FaceNormal_DDFB(y_, 1, iFace, jFace, kFace, iBlockFace)
          else
             AreaY = 0.0
          end if
          if(nK > 1)then
             AreaZ = FaceNormal_DDFB(z_, 1, iFace, jFace, kFace, iBlockFace)
          else
             AreaZ = 0.0
          end if
       end if
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 1e-30)then

          if(DoTestCell)write(*,*)'set_cell_values_x: Area2=', Area2

          ! The face is at the pole
          Normal_D = Xyz_DGB(:,iFace,jFace,kFace,iBlockFace) &
               -     Xyz_DGB(:,iLeft,jLeft,kLeft,iBlockFace)
          Normal_D=Normal_D/sqrt(sum(Normal_D**2))
          Area =0.0
          Area2=0.0
       else
          Area = sqrt(Area2)
          Normal_D=(/AreaX, AreaY, AreaZ/)/Area
       end if
    end if

    call set_cell_values_common

  end subroutine set_cell_values_x

  !===========================================================================
  subroutine set_cell_values_y

    iLeft = iFace; jLeft = jFace - 1; kLeft = kFace

    if(.not.IsCartesian)then

       if(IsRzGeometry)then
          AreaX = 0.0
          AreaY = CellFace_DFB(2,iFace,jFace,kFace,iBlockFace)
          AreaZ = 0.0
       else
          AreaX = FaceNormal_DDFB(x_, 2, iFace, jFace, kFace, iBlockFace)
          AreaY = FaceNormal_DDFB(y_, 2, iFace, jFace, kFace, iBlockFace)
          if(nK > 1)then
             AreaZ = FaceNormal_DDFB(z_, 2, iFace, jFace, kFace, iBlockFace)
          else
             AreaZ = 0.0
          end if
       end if
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 1e-30)then

          if(DoTestCell)write(*,*)'set_cell_values_y: Area2=', Area2

          !The face is at the pole
          Normal_D = Xyz_DGB(:,iFace,jFace,kFace,iBlockFace) &
               -     Xyz_DGB(:,iLeft,jLeft,kLeft,iBlockFace)
          Normal_D=Normal_D/sqrt(sum(Normal_D**2))
          Area = 0.0
          Area2= 0.0
       else
          Area = sqrt(Area2)
          Normal_D = (/AreaX, AreaY, AreaZ/)/Area

       end if

    end if

    call set_cell_values_common

  end subroutine set_cell_values_y
  !===========================================================================

  subroutine set_cell_values_z

    iLeft = iFace; jLeft = jFace; kLeft = kFace - 1

    if(.not.IsCartesian)then

       AreaX = FaceNormal_DDFB(x_, 3, iFace, jFace, kFace, iBlockFace)
       AreaY = FaceNormal_DDFB(y_, 3, iFace, jFace, kFace, iBlockFace)
       AreaZ = FaceNormal_DDFB(z_, 3, iFace, jFace, kFace, iBlockFace)
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 1e-30)then

          if(DoTestCell)write(*,*)'set_cell_values_z: Area2 = ', Area2

          !The face is at the pole
          Normal_D = Xyz_DGB(:,iFace,jFace,kFace,iBlockFace) &
               -     Xyz_DGB(:,iLeft,jLeft,kLeft,iBlockFace)
          Normal_D=Normal_D/sqrt(sum(Normal_D**2))

          Area = 0.0
          Area2= 0.0
       else
          Area = sqrt(Area2)
          Normal_D = (/AreaX, AreaY, AreaZ/)/Area
       end if
    end if

    call set_cell_values_common

  end subroutine set_cell_values_z

  !==========================================================================

  subroutine set_cell_values_common

    use ModPhysics, ONLY: Io2No_V, UnitU_
    use ModGeometry, ONLY: r_BLK
    real :: r
    !--------------------------------------------------------------------

    if(DoTestCell .and. .not.IsCartesian) &
         write(*,*)'Area2,AreaX,AreaY,AreaZ,Normal_D=', &
         Area2, AreaX, AreaY, AreaZ, Normal_D

    iRight= iFace; jRight = jFace; kRight = kFace

    IsBoundary = true_cell(iLeft, jLeft, kLeft, iBlockFace) &
         .neqv.  true_cell(iRight,jRight,kRight,iBlockFace)

    HallCoeff     = -1.0
    BiermannCoeff = -1.0
    if(UseHallResist .or. UseBiermannBattery)then
       if(UseHallResist)then
          HallCoeff = HallFactor_DF(iDimFace,iFace,jFace,kFace)
          BiermannCoeff = HallCoeff
       else
          BiermannCoeff = &
               0.5*( IonMassPerCharge_G(iLeft,jLeft,kLeft) &
               +     IonMassPerCharge_G(iRight,jRight,kRight) )
       end if

       ! No need to calculate anything related to Hall MHD if
       ! we use semi-implicit scheme to update the magnetic field
       ! and there is no electron pressure equation that needs Hall velocity.
       if(.not.UseElectronPressure .and. .not.DoHallInduction) HallCoeff = -1.0
    end if

    ViscoCoeff = 0.0
    if(UseViscosity) ViscoCoeff = ViscoFactor_DF(iDimFace,iFace,jFace,kFace)

    ! Calculate -grad(pe)/(n_e * e) term for Hall MHD if needed
    UseHallGradPe = BiermannCoeff > 0.0 .and. &
         (UseElectronPressure .or. ElectronPressureRatio > 0.0 .or. &
         .not.UseIdealEos)

    Eta       = 0.0
    if(UseResistiveFlux) Eta = 0.5* &
         ( Eta_GB(iLeft, jLeft  ,kLeft,iBlockFace) &
         + Eta_GB(iRight,jRight,kRight,iBlockFace))

    if(.not.IsCartesian)then
       NormalX = Normal_D(x_); NormalY = Normal_D(y_); NormalZ = Normal_D(z_)
       AreaX = Area*NormalX; AreaY = Area*NormalY; AreaZ = Area*NormalZ

       ! InvDxyz is needed for the time step limit of the explicit evaluation
       ! of the diffusion operator
       InvDxyz = 1.0/sqrt( sum( &
            ( Xyz_DGB(:,iRight,jRight,kRight, iBlockFace)          &
            - Xyz_DGB(:,iLeft, jLeft  ,kLeft, iBlockFace))**2) )
    end if

    if(UseClimit)then
       r = 0.5*(r_BLK(iLeft,  jLeft,  kLeft,  iBlockFace) &
            +   r_BLK(iRight, jRight, kRight, iBlockFace))
       if(r < rClimit)then
          Climit = Io2No_V(UnitU_)*ClimitDim
       else
          Climit = -1.0
       end if

       if(DoTestCell)write(*,*)'set_cell_values_common: Climit=', Climit

    end if

  end subroutine set_cell_values_common

  !==========================================================================

  subroutine get_numerical_flux(Flux_V)

    use ModWaves, ONLY: UseWavePressure
    use ModAdvance, ONLY: DoReplaceDensity, State_VGB, UseMultiSpecies
    use ModCharacteristicMhd, ONLY: get_dissipation_flux_mhd
    use ModCoordTransform, ONLY: cross_product
    use ModMain, ONLY: UseHyperbolicDivb, SpeedHyp, UseDtFixed
    use ModFaceGradient, ONLY: get_face_gradient, get_face_curl
    use ModPhysics,  ONLY: UnitTemperature_, UnitN_, Si2No_V, cLight
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModUserInterface ! user_material_properties

    real,    intent(out):: Flux_V(nFlux)

    real :: State_V(nVar)

    real :: Cmax
    real :: DiffBn_D(3), DiffE
    real :: EnLeft, EnRight, PeLeft, PeRight, PwaveLeft, PwaveRight, Jx, Jy, Jz
    real :: uLeft_D(3), uRight_D(3)
    real :: dB0_D(3), Current_D(3)

    real :: GradPe_D(3)
    real :: InvElectronDens
    integer :: i, j, k
    real :: NatomicSi, TeSi

    real, save :: b_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    character(len=*), parameter:: NameSub = 'get_numerical_flux'
    !-----------------------------------------------------------------------

    ! Initialize diffusion coefficient for time step restriction
    DiffCoef = 0.0

    if(UseMultiSpecies .and. DoReplaceDensity)then
       StateLeft_V (Rho_)=sum( StateLeft_V(SpeciesFirst_:SpeciesLast_) )
       StateRight_V(Rho_)=sum( StateRight_V(SpeciesFirst_:SpeciesLast_) )
    end if

    ! Calculate current for the face if needed for (Hall) resistivity
    if(HallCoeff > 0.0 .or. Eta > 0.0) then
       if(IsNewBlockCurrent) b_DG = State_VGB(Bx_:Bz_,:,:,:,iBlockFace)
       call get_face_curl(iDimFace, iFace, jFace, kFace, iBlockFace, &
            IsNewBlockCurrent, b_DG, Current_D)
       Jx = Current_D(1); Jy = Current_D(2); Jz = Current_D(3)
    end if

    ! Calculateing stress tensor for viscosity Visco_DDI
    if(ViscoCoeff > 0.0 ) then
       call get_viscosity_tensor(iDimFace, iFace, jFace, kFace,&
            iBlockFace,iFluidMin,iFluidMax,ViscoCoeff)
    end if

    if(Eta > 0.0)then
       EtaJx = Eta*Jx
       EtaJy = Eta*Jy
       EtaJz = Eta*Jz
    end if

    if(HallCoeff > 0.0)then
       HallJx = HallCoeff*Jx
       HallJy = HallCoeff*Jy
       HallJz = HallCoeff*Jz
    end if

    if(UseHallGradPe)then

       if(IsNewBlockGradPe)then

          if(.not.allocated(Pe_G)) &
               allocate(Pe_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

          ! Obtain electron pressure
          if(.not.UseIdealEos .and. .not.UseElectronPressure)then
             do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
                call user_material_properties(State_VGB(:,i,j,k,iBlockFace), &
                     i, j, k, iBlockFace, TeOut=TeSi, NatomicOut=NatomicSi)
                ! Single temperature mode: electron temperature is the same
                ! as the ion temperature.
                ! Subtract ion pressure from total pressure.
                Pe_G(i,j,k) = State_VGB(p_,i,j,k,iBlockFace) &
                     - NatomicSi*Si2No_V(UnitN_)*TeSi*Si2No_V(UnitTemperature_)
             end do; end do; end do 
          elseif(UseElectronPressure)then
             Pe_G = State_VGB(Pe_,:,:,:,iBlockFace)
          elseif(IsMhd)then
             Pe_G = State_VGB(p_,:,:,:,iBlockFace)*PePerPtotal
          else
             Pe_G = sum(State_VGB(iPIon_I,:,:,:,iBlockFace),DIM=1) &
                  *ElectronPressureRatio
          end if
       end if

       ! Calculate face centered grad(Pe)
       call get_face_gradient(iDimFace, iFace, jFace, kFace, iBlockFace, &
            IsNewBlockGradPe, Pe_G, GradPe_D)

       ! Calculate 1/(n_e * e)
       if(UseMultiIon)then
          InvElectronDens = BiermannCoeff/(0.5* &
               sum((StateLeft_V(iRhoIon_I) + StateRight_V(iRhoIon_I))&
               *ChargeIon_I/MassIon_I))
       else
          InvElectronDens = BiermannCoeff &
               /(0.5*(StateLeft_V(Rho_) + StateRight_V(Rho_)))
       end if

       ! Calculate grad(Pe)/(n_e * e)
       GradXPeNe = GradPe_D(1)*InvElectronDens
       GradYPeNe = GradPe_D(2)*InvElectronDens
       GradZPeNe = GradPe_D(3)*InvElectronDens

    end if

    if(DoRadDiffusion)then
       call get_radiation_energy_flux(iDimFace, iFace, jFace, kFace, &
            iBlockFace, StateLeft_V, StateRight_V, Normal_D, &
            RadDiffCoef, EradFlux)
       DiffCoef = DiffCoef + RadDiffCoef
    end if

    if(DoHeatConduction)then
       call get_heat_flux(iDimFace, iFace, jFace, kFace, iBlockFace, &
            StateLeft_V, StateRight_V, Normal_D, &
            HeatCondCoefNormal, HeatFlux)
       DiffCoef = DiffCoef + HeatCondCoefNormal
    end if

    if(DoIonHeatConduction)then
       call get_ion_heat_flux(iDimFace, iFace, jFace, kFace, iBlockFace, &
            StateLeft_V, StateRight_V, Normal_D, &
            HeatCondCoefNormal, IonHeatFlux)
       DiffCoef = DiffCoef + HeatCondCoefNormal
    end if

    if(UseB)then
       if(DoRoe)then
          if(IsBoundary)then
             uLeft_D  = StateLeft_V(Ux_:Uz_)
             uRight_D = StateRight_V(Ux_:Uz_)
          else
             !Since the divB source term is calculated using the
             !cell centered velocity, the numerical diffusion
             !for the normal magnetic field should be evaluated
             !in terms of the cell centered velocity too
             uLeft_D = &
                  State_VGB(RhoUx_:RhoUz_,iLeft,jLeft,kLeft,iBlockFace) &
                  /        State_VGB(Rho_,iLeft,jLeft,kLeft,iBlockFace)
             uRight_D = &
                  State_VGB(RhoUx_:RhoUz_,iRight,jRight,kRight,iBlockFace) &
                  /        State_VGB(Rho_,iRight,jRight,kRight,iBlockFace)
          end if

          if(UseB0)then
             dB0_D = B0_DGB(:, iLeft,  jLeft,  kLeft,  iBlockFace)  &
                  -  B0_DGB(:, iRight, jRight, kRight, iBlockFace)
          else
             dB0_D = 0.0
          end if

          call get_dissipation_flux_mhd(Normal_D,         &
               StateLeft_V, StateRight_V,                 &
               (/B0x,B0y,B0z/), dB0_D,                    &
               uLeft_D, uRight_D, DeltaBnL, DeltaBnR,     &
               IsBoundary, .false.,                       &
               DissipationFlux_V, cMax, Unormal_I(1))

          Unormal_I=Unormal_I(1)
       end if
       if(UseRS7 .or. UseLindeFix)then
          ! Sokolov's algorithm
          ! Calculate the jump in the normal magnetic field vector
          DiffBn_D = Normal_D* &
               0.5*sum((StateRight_V(Bx_:Bz_) - StateLeft_V(Bx_:Bz_))*Normal_D)

          ! Remove the jump in the normal magnetic field
          StateLeft_V(Bx_:Bz_)  =  StateLeft_V(Bx_:Bz_)  + DiffBn_D
          StateRight_V(Bx_:Bz_) =  StateRight_V(Bx_:Bz_) - DiffBn_D

          ! The energy jump is also modified by 
          ! 1/2(Br^2 - Bl^2) = 1/2(Br-Bl)*(Br+Bl)
          ! We store half of this in DiffE
          DiffE = &
               0.5*sum((StateRight_V(Bx_:Bz_) + StateLeft_V(Bx_:Bz_))*DiffBn_D)

          DiffBb = sum(DiffBn_D**2)
       end if
    end if

    ! Calculate average state (used by most solvers and also by bCrossArea_D)
    if(DoSimple)then
       State_V = StateLeft_V
       call get_physical_flux(State_V, B0x, B0y, B0z,&
            StateLeftCons_V, Flux_V, Unormal_I, Enormal, Pe, Pwave)
    else
       State_V = 0.5*(StateLeft_V + StateRight_V)
    end if

    if(DoLf .or. DoHll .or. DoLfdw .or. DoHlldw .or. DoAw .or. &
         DoRoe .or. DoRoeOld .or. &
         DoLfNeutral .or. DoHllNeutral .or. DoLfdwNeutral .or. &
         DoHlldwNeutral .or. DoAwNeutral .or. DoHllcNeutral)then
       ! These solvers use left and right fluxes
       call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
            StateLeftCons_V, FluxLeft_V, UnLeft_I, EnLeft, PeLeft, PwaveLeft)

       call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
            StateRightCons_V, FluxRight_V, UnRight_I, EnRight, PeRight, &
            PwaveRight)

       if(UseRS7)then
          call modify_flux(FluxLeft_V, UnLeft_I(1))
          call modify_flux(FluxRight_V, UnRight_I(1))
       end if
    end if

    if(UseB .and. (UseMultiIon .or. .not. IsMhd))then
       ! Calculate bCrossArea_D to be used for J in the J x B source term
       ! for the individual ion fluids in calc_sources.f90.
       ! The upwinded discretization of the current is J = sum(A x B) / V

       bCrossArea_D = cross_product(AreaX, AreaY, AreaZ, State_V(Bx_:Bz_))

       if(DoTestCell)then
          write(*,'(a,3es13.5)')'bCrossArea_D        =',bCrossArea_D
          write(*,'(a,3es13.5)')'AreaX, AreaY, AreaZ =',AreaX, AreaY, AreaZ
          write(*,'(a,3es13.5)')'State_V(Bx_:Bz_)    =',State_V(Bx_:Bz_)
       end if
    end if

    ! Initialize CmaxDt so we can take maximum of ion and neutral cmaxdt
    CmaxDt = 0.0

    ! Calculate ion fluxes first (all ion fluids together)
    if(UseB)then
       iFluidMin = 1; iFluidMax = IonLast_
       iVarMin = iRho_I(iFluidMin); iVarMax = iP_I(iFluidMax)
       iEnergyMin = nVar + iFluidMin; iEnergyMax = nVar + iFluidMax

       if(DoSimple)then
          call simple_flux
       elseif(DoLf)then
          call lax_friedrichs_flux
       elseif(DoHll)then
          call harten_lax_vanleer_flux
       elseif(DoLfdw .or. DoHlldw)then
          call dominant_wave_flux(DoLfdw)
       elseif(DoHlld)then
          call hlld_flux
       elseif(DoAw)then
          call artificial_wind
       elseif(DoRoeOld)then
          call roe_solver(Flux_V)
       elseif(DoRoe)then
          call roe_solver_new
       else
          call stop_mpi(NameSub//': Unknown flux type for ions')
       end if
    end if

    if(UseLindeFix)then

       if(UseHyperbolicDivb) then
          ! Overwrite the flux of the Hyp field with the Lax-Friedrichs flux
          Cmax = max(Cmax, SpeedHyp)
          Flux_V(Hyp_) = 0.5*(FluxLeft_V(Hyp_) + FluxRight_V(Hyp_) &
               - Cmax*(StateRight_V(Hyp_) - StateLeft_V(Hyp_)))
       end if

       if(.not.UseRS7)then
          ! Linde's idea: use Lax-Friedrichs flux for Bn
          Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) - Cmax*DiffBn_D

          ! Fix the energy diffusion
          Flux_V(Energy_) = Flux_V(Energy_) - Cmax*DiffE
       end if
    end if

    if (UseEfield) then
       Cmax = max(Cmax, clight)
       Flux_V(iVarUseCmax_I) = &
            0.5*(FluxLeft_V(iVarUseCmax_I) + FluxRight_V(iVarUseCmax_I) &
            - Cmax*(StateRightCons_V(iVarUseCmax_I)                     &
            -       StateLeftCons_V(iVarUseCmax_I)))
    end if

    ! Calculate neutral fluxes one-by-one
    do iFluidMin = NeutralFirst_, nFluid
       iFluidMax = iFluidMin
       iFluid = iFluidMin
       call select_fluid
       iVarMin = iRho; iVarMax = iP
       iEnergyMin = iEnergy; iEnergyMax = iEnergy
       if(DoLfNeutral)then
          call lax_friedrichs_flux
       elseif(DoHllNeutral)then
          call harten_lax_vanleer_flux
       elseif(DoLfdwNeutral .or. DoHlldwNeutral)then
          call dominant_wave_flux(DoLfdwNeutral)
       elseif(DoAwNeutral)then
          call artificial_wind
       elseif(DoGodunovNeutral)then
          call godunov_flux
       elseif(DoHllcNeutral)then
          call hllc_flux
       else
          call stop_mpi(NameSub//': Unknown flux type for neutrals')
       end if
    end do

    ! Restore
    iFluidMin = 1; iFluidMax = nFluid

    ! Multiply Flux by Area. This is needed in div Flux in update_states_MHD
    Flux_V = Flux_V*Area

    ! Increase maximum speed with the sum of diffusion speeds
    ! Resistivity, viscosity, heat conduction, radiation diffusion
    CmaxDt = CmaxDt + 2*(DiffCoef + ViscoCoeff + Eta)*InvDxyz

    ! Further limit timestep due to the hyperbolic cleaning equation
    if(UseHyperbolicDivb) CmaxDt = max(SpeedHyp, CmaxDt)

    if(DoTestCell)call write_test_info

  contains
    !==========================================================================
    subroutine modify_flux(Flux_V,Un)

      real, intent(in)   :: Un
      real, intent(inout):: Flux_V(nFlux)
      !----------------------------------------------------------------------
      Flux_V(RhoUx_:RhoUz_) = Flux_V(RhoUx_:RhoUz_) + 0.5*DiffBb*Normal_D
      Flux_V(Energy_)       = Flux_V(Energy_)       + Un*DiffBb

    end subroutine modify_flux
    !==========================================================================
    subroutine roe_solver_new

      Flux_V = 0.5*(FluxLeft_V + FluxRight_V) + DissipationFlux_V
      CmaxDt = Cmax

    end subroutine roe_solver_new
    !==========================================================================
    subroutine simple_flux

      real    :: Cmax_I(nFluid)
      !----------------------------------------------------------------------
      ! This is needed for the time step constraint only (CmaxDt)
      if(UseDtFixed)then 
         CmaxDt = 1.0
      else
         call get_speed_max(State_V, B0x, B0y, B0z, Cmax_I = Cmax_I)
      end if

    end subroutine simple_flux

    !==========================================================================
    subroutine lax_friedrichs_flux

      real    :: Cmax_I(nFluid)
      !----------------------------------------------------------------------
      call get_speed_max(State_V, B0x, B0y, B0z, Cmax_I = Cmax_I)

      Cmax = maxval(Cmax_I(iFluidMin:iFluidMax))
      Flux_V(iVarMin:iVarMax) = &
           0.5*(FluxLeft_V(iVarMin:iVarMax) + FluxRight_V(iVarMin:iVarMax) &
           - Cmax*(StateRightCons_V(iVarMin:iVarMax) &
           -       StateLeftCons_V(iVarMin:iVarMax)))
      ! energy flux
      Flux_V(iEnergyMin:iEnergyMax) = &
           0.5*(FluxLeft_V(iEnergyMin:iEnergyMax) + &
           FluxRight_V(iEnergyMin:iEnergyMax) &
           - Cmax*(StateRightCons_V(iEnergyMin:iEnergyMax) &
           -       StateLeftCons_V(iEnergyMin:iEnergyMax)))

      Unormal_I(iFluidMin:iFluidMax) = 0.5* &
           (UnLeft_I(iFluidMin:iFluidMax) + UnRight_I(iFluidMin:iFluidMax))

      ! These quantities should be calculated with the ion fluxes
      if(iFluidMin == 1)then
         Enormal   = 0.5*(EnLeft + EnRight)
         if(UseElectronPressure) &
              Unormal_I(eFluid_) = 0.5*(UnLeft_I(eFluid_) + UnRight_I(eFluid_))
         if(UseMultiIon)then
            Pe = 0.5*(PeLeft + PeRight)
            if(UseWavePressure) Pwave = 0.5*(PwaveLeft + PwaveRight)
         end if
      end if

      if (DoTestCell) then 
         write(*,*) 'iFluidMin, iFluidMax, Cmax =', iFluidMin, iFluidMax, Cmax
      end if

    end subroutine lax_friedrichs_flux
    !==========================================================================
    subroutine harten_lax_vanleer_flux

      real, dimension(nFluid) :: CleftStateLeft_I,   CleftStateHat_I, &
           Cmax_I, CrightStateRight_I, CrightStateHat_I
      real :: Cleft, Cright, WeightLeft, WeightRight, Diffusion
      !-----------------------------------------------------------------------

      call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
           Cleft_I =CleftStateLeft_I)

      call get_speed_max(StateRight_V, B0x, B0y, B0z, &
           Cright_I=CrightStateRight_I)

      call get_speed_max(State_V, B0x, B0y, B0z, &
           Cmax_I = Cmax_I, &
           Cleft_I = CleftStateHat_I, Cright_I = CrightStateHat_I)

      Cmax   =maxval(Cmax_I(iFluidMin:iFluidMax))
      Cleft  =min(0.0, &
           minval(CleftStateLeft_I(iFluidMin:iFluidMax)), &
           minval(CleftStateHat_I(iFluidMin:iFluidMax)))
      Cright =max(0.0, &
           maxval(CrightStateRight_I(iFluidMin:iFluidMax)), &
           maxval(CrightStateHat_I(iFluidMin:iFluidMax)))

      WeightLeft  = Cright/(Cright - Cleft)
      WeightRight = 1.0 - WeightLeft
      Diffusion   = Cright*WeightRight

      Flux_V(iVarMin:iVarMax) = &
           ( WeightRight*FluxRight_V(iVarMin:iVarMax)     &
           + WeightLeft*FluxLeft_V(iVarMin:iVarMax)       &
           - Diffusion*(StateRightCons_V(iVarMin:iVarMax) &
           -            StateLeftCons_V(iVarMin:iVarMax)) )

      ! Energy flux
      Flux_V(iEnergyMin:iEnergyMax) = &
           ( WeightRight*FluxRight_V(iEnergyMin:iEnergyMax)     &
           + WeightLeft*FluxLeft_V(iEnergyMin:iEnergyMax)       &
           - Diffusion*(StateRightCons_V(iEnergyMin:iEnergyMax) &
           -            StateLeftCons_V(iEnergyMin:iEnergyMax)) )

      ! Weighted average of the normal speed
      Unormal_I(iFluidMin:iFluidMax) = &
           WeightRight*UnRight_I(iFluidMin:iFluidMax) &
           + WeightLeft*UnLeft_I(iFluidMin:iFluidMax)

      ! These quantities should be calculated with the ion fluxes
      if(iFluidMin == 1)then
         Enormal   = WeightRight*EnRight + WeightLeft*EnLeft
         if(UseElectronPressure) Unormal_I(eFluid_) = &
              WeightRight*UnRight_I(eFluid_) + WeightLeft*UnLeft_I(eFluid_)
         if(UseMultiIon)then
            Pe = WeightRight*PeRight   + WeightLeft*PeLeft
            if(UseWavePressure) &
                 Pwave = WeightRight*PwaveRight + WeightLeft*PwaveLeft
         end if
      end if

    end subroutine harten_lax_vanleer_flux
    !==========================================================================
    subroutine dominant_wave_flux(DoLf)
      use ModPhysics,  ONLY: Gamma_I

      logical, intent(in):: DoLf

      real, dimension(nFluid):: CleftStateLeft_I, CleftStateHat_I, &
           Cmax_I, CrightStateRight_I, CrightStateHat_I

      real:: DeltaCons_V(nVar), DeltaFlux_V(nVar)
      real:: Cleft, Cright
      real:: WeightLeft=0.0, WeightRight=0.0, Diffusion=0.0, DiffusionDw
      real:: Nu1, Nu2, Nu, CsoundL, CsoundR, Cdw
      !-----------------------------------------------------------------------
      ! Get the max, left and right speeds for HLL (and DW?)
      call get_speed_max(State_V, B0x, B0y, B0z, &
           Cmax_I = Cmax_I, &
           Cleft_I = CleftStateHat_I, Cright_I = CrightStateHat_I)
      Cmax = maxval(Cmax_I(iFluidMin:iFluidMax))

      ! Andrea Mignone's hybridization parameters
      ! Pressure jump detector
      Nu1 = 1.0 - min(StateLeft_V(p_), StateRight_V(p_)) &
           /max(StateLeft_V(p_), StateRight_V(p_))

      ! Maybe the max speed from get_speed_max is better
      CsoundL = sqrt(Gamma_I(1)*StateLeft_V(p_)/StateLeft_V(Rho_))
      CsoundR = sqrt(Gamma_I(1)*StateRight_V(p_)/StateRight_V(Rho_))

      ! Rarefaction (and shock) detector ?!
      Nu2 = min(0.5*abs(UnRight_I(1) - UnLeft_I(1)) &
           /(CsoundL + CsoundR), 1.0)

      ! HLLE flux weight is Nu, dominant wave flux has weight (1-Nu)
      Nu = max(Nu1, Nu2)

      ! Round to 0 or 1 if close
      if(Nu < 1e-6)   Nu = 0.0
      if(Nu > 1-1e-6) Nu = 1.0

      if(Nu < 1.0 .or. (Nu > 0.0 .and. .not.DoLf))then
         call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
              Cleft_I =CleftStateLeft_I)

         call get_speed_max(StateRight_V, B0x, B0y, B0z, &
              Cright_I=CrightStateRight_I)

         Cleft  =min(0.0, &
              minval(CleftStateLeft_I(iFluidMin:iFluidMax)), &
              minval(CleftStateHat_I(iFluidMin:iFluidMax)))
         Cright =max(0.0, &
              maxval(CrightStateRight_I(iFluidMin:iFluidMax)), &
              maxval(CrightStateHat_I(iFluidMin:iFluidMax)))
      end if

      if(Nu > 0.0)then
         ! LF or HLLE scheme is needed
         if(DoLf)then
            WeightLeft  = 0.5
            WeightRight = 0.5
            Diffusion   = 0.5*Cmax
         else
            ! HLLE weights and diffusion
            WeightLeft  = Cright/(Cright - Cleft)
            WeightRight = 1.0 - WeightLeft
            Diffusion   = Cright*WeightRight
         end if
      end if

      if(Nu < 1.0)then
         ! DW scheme is needed
         ! Jump in conservative variables
         DeltaCons_V(iVarMin:iVarMax) = &
              StateRightCons_V(iVarMin:iVarMax) - &
              StateLeftCons_V(iVarMin:iVarMax)

         ! Overwrite pressure with energy
         DeltaCons_V(iP_I(iFluidMin:iFluidMax)) = &
              StateRightCons_V(iEnergyMin:iEnergyMax) - &
              StateLeftCons_V(iEnergyMin:iEnergyMax)

         ! Jump in flux
         DeltaFlux_V(iVarMin:iVarMax) = &
              FluxRight_V(iVarMin:iVarMax) - &
              FluxLeft_V( iVarMin:iVarMax)

         DeltaFlux_V(iP_I(iFluidMin:iFluidMax)) = &
              FluxRight_V(iEnergyMin:iEnergyMax) - &
              FluxLeft_V(iEnergyMin:iEnergyMax)

         ! Dominant wave diffusion coefficient 0.5*dF.dU/||dU||
         Cdw = dot_product(DeltaFlux_V(iVarMin:iVarMax), &
              DeltaCons_V(iVarMin:iVarMax)) &
              / max(sum(DeltaCons_V(iVarMin:iVarMax)**2), 1e-30)
         ! Dominant wave speed is limited between Cleft and Cright
         DiffusionDw = 0.5*abs(max(Cleft, min(Cright, Cdw)))

         ! Combine HLLE and DW weights and diffusion
         WeightLeft  = Nu*WeightLeft  + (1-Nu)*0.5
         WeightRight = Nu*WeightRight + (1-Nu)*0.5
         Diffusion   = Nu*Diffusion   + (1-Nu)*DiffusionDw
      end if

      ! LF-DW or HLL-DW flux
      Flux_V(iVarMin:iVarMax) = &
           ( WeightRight*FluxRight_V(iVarMin:iVarMax)     &
           + WeightLeft*FluxLeft_V(iVarMin:iVarMax)       &
           - Diffusion*(StateRightCons_V(iVarMin:iVarMax) &
           -            StateLeftCons_V(iVarMin:iVarMax)) )

      ! Energy flux
      Flux_V(iEnergyMin:iEnergyMax) = &
           ( WeightRight*FluxRight_V(iEnergyMin:iEnergyMax)     &
           + WeightLeft*FluxLeft_V(iEnergyMin:iEnergyMax)       &
           - Diffusion*(StateRightCons_V(iEnergyMin:iEnergyMax) &
           -            StateLeftCons_V(iEnergyMin:iEnergyMax)) )

      ! Weighted average of the normal speed
      Unormal_I(iFluidMin:iFluidMax) = &
           WeightRight*UnRight_I(iFluidMin:iFluidMax) &
           + WeightLeft*UnLeft_I(iFluidMin:iFluidMax)

      ! These quantities should be calculated with the ion fluxes
      if(iFluidMin == 1)then
         Enormal = WeightRight*EnRight + WeightLeft*EnLeft
         if(UseElectronPressure) Unormal_I(eFluid_) = &
              WeightRight*UnRight_I(eFluid_) + WeightLeft*UnLeft_I(eFluid_)
         if(UseMultiIon)then
            Pe = WeightRight*PeRight + WeightLeft*PeLeft
            if(UseWavePressure) &
                 Pwave = WeightRight*PwaveRight + WeightLeft*PwaveLeft
         end if
      end if

    end subroutine dominant_wave_flux
    !==========================================================================
    subroutine artificial_wind

      real, dimension(nFluid) :: Cleft_I, Cright_I, Cmax_I
      real :: Cleft, Cright, WeightLeft, WeightRight, Diffusion
      !-----------------------------------------------------------------------

      ! The propagation speeds are modified by the DoAw = .true. !
      call get_speed_max(State_V, B0x, B0y, B0z,  &
           Cleft_I = Cleft_I, Cright_I = Cright_I, Cmax_I = Cmax_I, &
           UseAwSpeedIn = .true.)

      Cmax   = maxval(Cmax_I(iFluidMin:iFluidMax))
      Cleft  = min(0.0, minval(Cleft_I(iFluidMin:iFluidMax)))
      Cright = max(0.0, maxval(Cright_I(iFluidMin:iFluidMax)))

      WeightLeft  = Cright/(Cright - Cleft)
      WeightRight = 1.0 - WeightLeft
      Diffusion   = Cright*WeightRight

      Flux_V(iVarMin:iVarMax) = &
           ( WeightRight*FluxRight_V(iVarMin:iVarMax)     &
           + WeightLeft*FluxLeft_V(iVarMin:iVarMax)       &
           - Diffusion*(StateRightCons_V(iVarMin:iVarMax) &
           -            StateLeftCons_V(iVarMin:iVarMax)) )

      ! Energy flux
      Flux_V(iEnergyMin:iEnergyMax) = &
           ( WeightRight*FluxRight_V(iEnergyMin:iEnergyMax)     &
           + WeightLeft*FluxLeft_V(iEnergyMin:iEnergyMax)       &
           - Diffusion*(StateRightCons_V(iEnergyMin:iEnergyMax) &
           -            StateLeftCons_V(iEnergyMin:iEnergyMax)) )

      ! Weighted average of the normal speed and electric field
      Unormal_I(iFluidMin:iFluidMax) = &
           WeightRight*UnRight_I(iFluidMin:iFluidMax) &
           + WeightLeft*UnLeft_I(iFluidMin:iFluidMax)

      ! These quantities should be calculated with the ion fluxes
      if(iFluidMin == 1)then
         Enormal   = WeightRight*EnRight + WeightLeft*EnLeft
         if(UseElectronPressure) Unormal_I(eFluid_) = &
              WeightRight*UnRight_I(eFluid_) + WeightLeft*UnLeft_I(eFluid_)
         if(UseMultiIon)then
            Pe = WeightRight*PeRight + WeightLeft*PeLeft
            if(UseWavePressure) &
                 Pwave = WeightRight*PwaveRight + WeightLeft*PwaveLeft
         end if
      end if

    end subroutine artificial_wind
    !==========================================================================
    subroutine hlld_flux

      ! The HLLD scheme works for single ion fluid only

      use ModPhysics, ONLY: InvGammaMinus1, GammaMinus1
      use ModNumConst, ONLY: cTiny

      implicit none

      ! Rotated flux for vector variables
      real :: FluxRot_V(nFluxMhd)

      ! Needed as an argument for get_physical_flux
      real :: StateCons_V(nFlux)

      ! Left and right state (scalars and extra variables only)
      real :: DsL, DsRhoL, RhoL, pL, eL, PbL, PtotL, uDotB1L, Bt1L, Bt2L
      real :: DsR, DsRhoR, RhoR, pR, eR, PbR, PtotR, uDotB1R, Bt1R, Bt2R

      ! First left and right intermediate states
      real :: sL1, InvSLUn, SqRhoL1, Ut1L1, Ut2L1, B1t1L1, B1t2L1, uDotB1L1
      real :: sR1, InvSRUn, SqRhoR1, Ut1R1, Ut2R1, B1t1R1, B1t2R1, uDotB1R1
      real :: SqRhoLR1

      ! Total pressure in all intermediate states
      real :: Ptot12

      ! Intermediate HLLD state
      real :: Rho, Un, Ut1, Ut2, B1n, B1t1, B1t2, p, e
      real :: RhoUn, uDotB1, Bn, Bt1, Bt2

      ! Resistivity
      real :: FluxBx, FluxBy, FluxBz, B1x, B1y, B1z

      real :: Tmp, B1n2, Bn2, SignBn

      real :: sL, CleftStateLeft_I(nFluid), CleftStateRight_I(nFluid)
      real :: sR, CrightStateLeft_I(nFluid), CrightStateRight_I(nFluid)

      !-----------------------------------------------------------------------

      ! This is the choice made in the hlld_tmp code. May not be the best.
      call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
           Cleft_I = CleftStateLeft_I, Cright_I = CrightStateLeft_I)

      call get_speed_max(StateRight_V, B0x, B0y, B0z, &
           Cleft_I = CleftStateRight_I, Cright_I = CrightStateRight_I)

      sL = min(CleftStateLeft_I(1),  CleftStateRight_I(1))
      sR = max(CrightStateLeft_I(1), CrightStateRight_I(1))

      Cmax   = max(sR, -sL)
      CmaxDt = Cmax

      if(DoTestCell)then
         write(*,*)'hlld: StateLeft =',StateLeft_V
         write(*,*)'hlld: StateRight=',StateRight_V
         write(*,*)'hlld: sL, sR    =',sL,sR
      endif

      if(sL >= 0.) then
         call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
              StateCons_V, Flux_V, Unormal_I, Enormal, Pe, Pwave)
         if(UseRs7)call modify_flux(Flux_V, Unormal_I(1))
         RETURN
      end if

      if(sR <= 0.) then 
         call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
              StateCons_V, Flux_V, Unormal_I, Enormal, Pe, Pwave)
         if(UseRs7)call modify_flux(Flux_V, Unormal_I(1))
         RETURN
      end if

      ! Scalar variables
      RhoL = StateLeft_V(Rho_)
      pL   = StateLeft_V(p_)
      RhoR = StateRight_V(Rho_)
      pR   = StateRight_V(p_)

      ! Rotate vector variables into a coordinate system orthogonal to the face
      call rotate_state_vectors

      ! Use average normal field
      B1n    = 0.5*(B1nL + B1nR)
      B1n2   = B1n**2

      ! Full magnetic field
      Bn     = B1n   + B0n
      Bt1L   = B1t1L + B0t1
      Bt2L   = B1t2L + B0t2
      Bt1R   = B1t1R + B0t1
      Bt2R   = B1t2R + B0t2
      Bn2    = Bn**2

      ! Sign of the total normal field
      SignBn = sign(1., Bn)

      ! Magnetic pressures
      PbL = 0.5*(B1n2 + B1t1L**2 + B1t2L**2)
      PbR = 0.5*(B1n2 + B1t1R**2 + B1t2R**2)

      ! Total pressure including B1.B0 term
      PtotL = pL + PbL + B1n*B0n + B1t1L*B0t1 + B1t2L*B0t2
      PtotR = pR + PbR + B1n*B0n + B1t1R*B0t1 + B1t2R*B0t2

      ! Propagation speed relative to bulk speed
      DsL = sL - UnL
      DsR = sR - UnR

      ! HLLD speed is used in all intermediate states for Un
      DsRhoL = DsL*RhoL
      DsRhoR = DsR*RhoR
      Tmp = 1./(DsRhoR - DsRhoL)
      Un  = (DsRhoR*UnR - DsRhoL*UnL - PtotR + PtotL)*Tmp

      ! Total pressure in all intermediate states
      Ptot12 =(DsRhoR*PtotL - DsRhoL*PtotR + DsRhoR*DsRhoL*(UnR - UnL))*Tmp

      if(DoTestCell)write(*,*)'hlld: Un = ',Un

      if(Un >= 0.)then
         ! Density and scalars for left intermediate states
         InvSLUn = 1.0/(sL-Un)
         Rho     = DsRhoL*InvSLUn

         ! Tangential velocity and magnetic field 
         ! for the outer intermediate left state
         Tmp = DsRhoL*(sL-Un) - Bn2
         if(Tmp < cTiny) then
            Ut1L1  = Ut1L
            Ut2L1  = Ut2L
            B1t1L1 = B1t1L
            B1t2L1 = B1t2L
         else
            Tmp    = 1.0/Tmp
            Ut1L1  = Ut1L - Bn*Bt1L*(Un - UnL)*Tmp
            Ut2L1  = Ut2L - Bn*Bt2L*(Un - UnL)*Tmp
            B1t1L1 = Bt1L*(DsRhoL*DsL - Bn2)*Tmp - B0t1
            B1t2L1 = Bt2L*(DsRhoL*DsL - Bn2)*Tmp - B0t2
         end if

         ! Calculate energy of outer left intermediate state
         eL = InvGammaMinus1*pL + PbL + 0.5*RhoL*(UnL**2 + Ut1L**2 + Ut2L**2)
         uDotB1L  = UnL*B1n + Ut1L *B1t1L  + Ut2L *B1t2L
         uDotB1L1 = Un *B1n + Ut1L1*B1t1L1 + Ut2L1*B1t2L1
         e = (DsL*eL - PtotL*UnL + Ptot12*Un + Bn*(uDotB1L - uDotB1L1))*InvSLUn

         ! Left going Alfven speed
         SqRhoL1 = sqrt(Rho)
         sL1     = Un - abs(Bn)/SqRhoL1

         if(DoTestCell)write(*,*)'hlld: sL1=',sL1

         ! Check sign of left going Alfven wave
         if(sL1 >= 0.) then
            ! Use first left intermediate state
            Ut1    = Ut1L1
            Ut2    = Ut2L1
            B1t1   = B1t1L1
            B1t2   = B1t2L1
            uDotB1 = uDotB1L1
         else
            ! Calculate some values for first right intermediate state
            SqRhoR1 = sqrt(DsRhoR/(sR-Un))
            Tmp = DsRhoR*(sR-Un) - Bn2
            if(Tmp < cTiny) then
               Ut1R1  = Ut1R
               Ut2R1  = Ut2R
               B1t1R1 = B1t1R
               B1t2R1 = B1t2R
            else
               Tmp    = 1.0/Tmp
               Ut1R1  = Ut1R - Bn*Bt1R*(Un-UnR)*Tmp
               Ut2R1  = Ut2R - Bn*Bt2R*(Un-UnR)*Tmp
               B1t1R1 = Bt1R*(DsRhoR*DsR - Bn2)*Tmp - B0t1
               B1t2R1 = Bt2R*(DsRhoR*DsR - Bn2)*Tmp - B0t2
            end if

            ! second left intermediate state
            Tmp      = 1./(SqRhoL1 + SqRhoR1)
            SqRhoLR1 = SqRhoL1*SqRhoR1*SignBn
            Ut1 =(SqRhoL1*Ut1L1  + SqRhoR1*Ut1R1  + (B1t1R1-B1t1L1)*SignBn)*Tmp
            Ut2 =(SqRhoL1*Ut2L1  + SqRhoR1*Ut2R1  + (B1t2R1-B1t2L1)*SignBn)*Tmp
            B1t1=(SqRhoL1*B1t1R1 + SqRhoR1*B1t1L1 + (Ut1R1-Ut1L1)*SqRhoLR1)*Tmp
            B1t2=(SqRhoL1*B1t2R1 + SqRhoR1*B1t2L1 + (Ut2R1-Ut2L1)*SqRhoLR1)*Tmp
            uDotB1 = Un*B1n + Ut1*B1t1 + Ut2*B1t2

            ! Modify energy density with difference between 1st and 2nd states
            e      = e - SqRhoL1*(uDotB1L1 - uDotB1)*SignBn
         end if
      else  ! Un < 0
         ! Density and scalars for right intermediate states
         InvSRUn = 1.0/(sR-Un)
         Rho     = DsRhoR*InvSRUn

         ! Tangential velocity and magnetic field 
         ! for the outer intermediate right state
         Tmp = DsRhoR*(sR-Un) - Bn2
         if(Tmp < cTiny) then
            Ut1R1  = Ut1R
            Ut2R1  = Ut2R
            B1t1R1 = B1t1R
            B1t2R1 = B1t2R
         else
            Tmp    = 1.0/Tmp
            Ut1R1  = Ut1R - Bn*Bt1R*(Un - UnR)*Tmp
            Ut2R1  = Ut2R - Bn*Bt2R*(Un - UnR)*Tmp
            B1t1R1 = Bt1R*(DsRhoR*DsR - Bn2)*Tmp - B0t1
            B1t2R1 = Bt2R*(DsRhoR*DsR - Bn2)*Tmp - B0t2
         end if

         ! Calculate energy of outer right intermediate state
         eR = InvGammaMinus1*pR + PbR + 0.5*RhoR*(UnR**2 + Ut1R**2 + Ut2R**2)
         uDotB1R  = UnR*B1n + Ut1R*B1t1R + Ut2R*B1t2R
         uDotB1R1 = Un*B1n + Ut1R1*B1t1R1 + Ut2R1*B1t2R1
         e = (DsR*eR - PtotR*UnR + Ptot12*Un + Bn*(uDotB1R - uDotB1R1))*InvSRUn

         ! Right going Alfven speed
         SqRhoR1 = sqrt(Rho)
         sR1     = Un + abs(Bn)/SqRhoR1

         if(DoTestCell)write(*,*)'hlld: sR1=',sR1

         if(sR1 <= 0.) then
            ! Use first right intermediate state
            Ut1    = Ut1R1
            Ut2    = Ut2R1
            B1t1   = B1t1R1
            B1t2   = B1t2R1
            uDotB1 = uDotB1R1
         else
            ! Calculate some values for the first left intermediate state
            SqRhoL1 = sqrt(DsRhoL/(sL-Un))
            Tmp = DsRhoL*(sL-Un) - Bn2
            if(Tmp < cTiny) then
               Ut1L1  = Ut1L
               Ut2L1  = Ut2L
               B1t1L1 = B1t1L
               B1t2L1 = B1t2L
            else
               Tmp  = 1.0/Tmp
               Ut1L1  = Ut1L - Bn*Bt1L*(Un - UnL)*Tmp
               Ut2L1  = Ut2L - Bn*Bt2L*(Un - UnL)*Tmp
               B1t1L1 = Bt1L*(DsRhoL*DsL - Bn2)*Tmp - B0t1
               B1t2L1 = Bt2L*(DsRhoL*DsL - Bn2)*Tmp - B0t2
            end if

            ! Second right intermediate state
            Tmp      = 1./(SqRhoL1 + SqRhoR1)
            SqRhoLR1 = SqRhoL1*SqRhoR1*SignBn
            Ut1 =(SqRhoL1*Ut1L1 + SqRhoR1*Ut1R1 + (B1t1R1-B1t1L1)*SignBn)*Tmp
            Ut2 =(SqRhoL1*Ut2L1 + SqRhoR1*Ut2R1 + (B1t2R1-B1t2L1)*SignBn)*Tmp
            B1t1=(SqRhoL1*B1t1R1 + SqRhoR1*B1t1L1 + (Ut1R1-Ut1L1)*SqRhoLR1)*Tmp
            B1t2=(SqRhoL1*B1t2R1 + SqRhoR1*B1t2L1 + (Ut2R1-Ut2L1)*SqRhoLR1)*Tmp
            uDotB1 = Un*B1n + Ut1*B1t1 + Ut2*B1t2

            ! Modify energy density with difference between 1st and 2nd states
            e   = e + SqRhoR1*(uDotB1R1 - uDotB1)*SignBn
         end if
      end if

      ! Calculate flux from HLLD state but use pTot12 instead of p+B^2/2
      ! Note that p derived from e is always positive, 
      ! but p derived from pTot12 may not be
      RhoUn  = Rho*Un
      Bt1 = B1t1 + B0t1
      Bt2 = B1t2 + B0t2
      p = GammaMinus1* &
           (e - 0.5*(B1n2 + B1t1**2 + B1t2**2 + Rho*(Un**2+Ut1**2+Ut2**2)))

      if(DoTestCell)write(*,*)'hlld: State,pTot12=',&
           Rho,Un,Ut1,Ut2,B1n,B1t1,B1t2,p,pTot12

      Flux_V(Rho_)       = RhoUn
      FluxRot_V(RhoUn_)  = RhoUn*Un  - B1n*Bn  - B0n*B1n + pTot12
      FluxRot_V(RhoUt1_) = RhoUn*Ut1 - B1n*Bt1 - B0n*B1t1
      FluxRot_V(RhoUt2_) = RhoUn*Ut2 - B1n*Bt2 - B0n*B1t2
      FluxRot_V(B1n_)    = 0.0
      FluxRot_V(B1t1_)   = Un*Bt1 - Ut1*Bn
      FluxRot_V(B1t2_)   = Un*Bt2 - Ut2*Bn
      Flux_V(p_)         = 0.5*(sR*UnL*pL - sL*UnR*pR  + sR*sL*(pR-pL))/(sR-sL)
      Flux_V(Energy_)    = Un*(e + pTot12) - Bn*uDotB1

      ! Rotate fluxes of vector variables back
      call rotate_flux_vector(FluxRot_V, Flux_V)

      ! Set normal velocity for all fluids (HLLD is for 1 fluid only)
      Unormal_I = Un

      if(UseRs7)call modify_flux(Flux_V, Unormal_I(1))

      if(Eta > 0.0)then
         ! Add flux corresponding to -curl Eta.J to induction equation
         FluxBx = NormalY*EtaJz - NormalZ*EtaJy
         FluxBy = NormalZ*EtaJx - NormalX*EtaJz
         FluxBz = NormalX*EtaJy - NormalY*EtaJx

         Flux_V(Bx_) = Flux_V(Bx_) + FluxBx
         Flux_V(By_) = Flux_V(By_) + FluxBy
         Flux_V(Bz_) = Flux_V(Bz_) + FluxBz

         ! Rotate back B1 of the HLLD state into the grid coordinates
         B1x = Normal_D(x_)*B1n + Tangent1_D(x_)*B1t1 + Tangent2_D(x_)*B1t2
         B1y = Normal_D(y_)*B1n + Tangent1_D(y_)*B1t1 + Tangent2_D(y_)*B1t2
         B1z = Normal_D(z_)*B1n + Tangent1_D(z_)*B1t1 + Tangent2_D(z_)*B1t2

         ! add B1.dB1/dt = div(B1 x EtaJ) term to the energy equation
         Flux_V(Energy_) = Flux_V(Energy_) &
              + B1x*FluxBx + B1y*FluxBy + B1z*FluxBz
      end if

    end subroutine hlld_flux

    !==========================================================================
    subroutine godunov_flux

      ! The Godunov flux works for hydro fluid (no magnetic field)
      ! Called for each fluid separately. Uses iFluid, iRho, ...

      use ModAdvance,  ONLY: UseElectronPressure
      use ModExactRS,  ONLY: wR, wL,  RhoL, RhoR, pL, pR, UnL, UnR, &
           UnStar, pStar, exact_rs_set_gamma, exact_rs_sample, exact_rs_pu_star
      use ModPhysics,  ONLY: InvGammaMinus1_I, Gamma_I, InvGammaMinus1
      use ModWaves,    ONLY: UseWavePressure, GammaWave

      real::Rho, Un, p, pTotal, e, StateStar_V(nVar)
      real::RhoSide,UnSide

      !The wave pressure: in the left state, right state and at the wall
      real::pWaveL=0.0,pWaveR=0.0, pWaveStar=0.0, pWaveSide=0.0
      real :: PeL, PeR, PeStar, PeSide
      real :: Adiabatic, Isothermal, GammaRatio, Factor
      integer::iVar
      !-----------------------------------------------------------------------
      ! Scalar variables
      RhoL = StateLeft_V(iRho)
      pL   = StateLeft_V(iP)
      RhoR = StateRight_V(iRho)
      pR   = StateRight_V(iP)

      UnL  = sum( StateLeft_V(iUx:iUz) *Normal_D )
      UnR  = sum( StateRight_V(iUx:iUz)*Normal_D )

      call exact_rs_set_gamma(Gamma_I(iFluid))

      if(UseWavePressure .and. iFluid==1)then
         ! Add the radiation/wave pressure to the total pressure
         ! This is for radiation pressure in CRASH applications.
         ! Increase maximum speed due to isotropic wave pressure
         pWaveL = (GammaWave - 1)*sum(StateLeft_V(WaveFirst_:WaveLast_))
         pWaveR = (GammaWave - 1)*sum(StateRight_V(WaveFirst_:WaveLast_))
         pL = pL + pWaveL
         pR = pR + pWaveR
      else
         !It is easier to add zero than using if(UseWavePressure)
         pWaveL = 0.0; pWaveR = 0.0
      end if
      if(UseElectronPressure .and. iFluid==1)then
         ! Add the electron pressure to the total pressure
         ! This is for ions coupled to electrons by collisions
         ! but there is no magnetic field (e.g. CRASH applications)
         ! GammaElectron=Gamma (ion) is assumed by the Godunov solver.
         PeL = StateLeft_V(Pe_)
         PeR = StateRight_V(Pe_)
         pL = pL + PeL
         pR = pR + PeR
      else
         ! StateLeft_V(iP) and StateRight_V(iP) already include the electron
         ! pressure
         PeL = 0.0; PeR = 0.0
      end if

      ! Take the parameters at the Contact Discontinuity (CD)
      call exact_rs_pu_star

      ! At strong shocks use the artificial wind scheme
      if((pStar > 2*pL .and. wL < 0.0).or.(pStar > 2*pR .and. wR > 0.0))then
         ! Temporary solution, should be the monotone numerical flux with 
         ! modified StateLeft_V and/or StateRight_V
         call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
              StateLeftCons_V, FluxLeft_V, UnLeft_I, EnLeft, PeLeft, PwaveLeft)

         call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
              StateRightCons_V, FluxRight_V, UnRight_I, EnRight, PeRight, &
              PwaveRight)

         call artificial_wind
         RETURN
      end if

      if(UnStar > 0.0)then
         ! The CD is to the right from the face
         ! The Left gas passes through the face
         RhoSide     = RhoL
         UnSide      = UnL
         StateStar_V = StateLeft_V
         pWaveSide   = pWaveL
         PeSide      = PeL
      else
         ! The CD is to the left from the face
         ! The Right gas passes through the face
         RhoSide     = RhoR
         UnSide      = UnR
         StateStar_V = StateRight_V
         pWaveSide   = pWaveR
         PeSide      = PeR
      end if

      ! Take the parameters at the face
      call exact_rs_sample(0.0, Rho, Un, p)

      ! In order the Riemann problem solution to be governed by the
      ! total pressure, the wave pressure should behave 
      ! adiabatically, with the same polytropic index as that for the gas

      Isothermal           = Rho/RhoSide
      Adiabatic            = Isothermal**Gamma_I(iFluid)
      pWaveStar            = pWaveSide*Adiabatic
      PeStar               = PeSide*Adiabatic

      ! Since the total pressure is not less than the adiabatic one
      ! the difference below is positive
      pTotal               = p
      p                    = pTotal - pWaveStar - PeStar

      StateStar_V(iRho)    = Rho
      StateStar_V(iUx:iUz) = StateStar_V(iUx:iUz) + (Un-UnSide)*Normal_D
      StateStar_V(P_)      = p
      do iVar=ScalarFirst_, ScalarLast_
         StateStar_V(iVar) = StateStar_V(iVar)*(Rho/RhoSide)
      end do

      ! Calculate flux 
      ! (1) calculate momenta
      StateStar_V(iRhoUx:iRhoUz) = StateStar_V(iUx:iUz) * Rho

      ! (2) take advective part of the flux
      Flux_V(1:nVar) = StateStar_V * Un 

      ! (3) add the pressure gradient
      ! also add the force due to the radiation and electron pressure gradient
      Flux_V(iRhoUx:iRhoUz) = Flux_V(iRhoUx:iRhoUz) + pTotal*Normal_D

      ! (4) energy flux: (e + p)*u
      ! also add the work done by the radiation and electron pressure gradient
      e = InvGammaMinus1_I(iFluid)*p &
           + 0.5*sum(StateStar_V(iRhoUx:iRhoUz)**2)/Rho
      Flux_V(iEnergyMin) = (e + pTotal)*Un

      Cmax                 = max(wR, -wL)
      CmaxDt               = Cmax
      Unormal_I(iFluid)    = Un

      if(iFluid == 1)then
         if(UseWavePressure)then
            GammaRatio = InvGammaMinus1*(GammaWave - 1)
            Factor = (1.0 - GammaRatio) + GammaRatio*Adiabatic/Isothermal
            do iVar = WaveFirst_, WaveLast_
               Flux_V(iVar) = Factor*StateStar_V(iVar)*Un
            end do
         end if
         if(UseElectronPressure) &
              Flux_V(Pe_) = (Adiabatic/Isothermal)*StateStar_V(Pe_)*Un

         if(DoRadDiffusion) Flux_V(Erad_) = Flux_V(Erad_) + EradFlux
      end if

    end subroutine godunov_flux

    !==========================================================================

    subroutine hllc_flux
      ! The HLLC scheme works for single ion fluid only
      ! HYDRO ONLY (NO MHD)
      implicit none

      real :: StateStarCons_V(nFlux)
      real :: UnStar

      ! Left and right state (scalars and extra variables only)
      real :: RhoL, totalPresL, sL
      real :: RhoR, totalPresR, sR

      real :: CleftStateLeft_I(nFluid), CleftStateRight_I(nFluid)
      real :: CrightStateLeft_I(nFluid), CrightStateRight_I(nFluid)
      !-----------------------------------------------------------------------

      call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
           Cleft_I = CleftStateLeft_I, Cright_I = CrightStateLeft_I)

      call get_speed_max(StateRight_V, B0x, B0y, B0z, &
           Cleft_I = CleftStateRight_I, Cright_I = CrightStateRight_I)

      sL = min(CleftStateLeft_I(1),  CleftStateRight_I(1))
      sR = max(CrightStateLeft_I(1), CrightStateRight_I(1))

      Cmax   = max(sR, -sL)
      CmaxDt = Cmax

      !=======================================================================
      ! Calculate intermediate states
      !=======================================================================

      RhoL = StateLeft_V(Rho_)
      RhoR = StateRight_V(Rho_)
      totalPresL = StateLeft_V(p_)  + PeLeft  + PwaveLeft
      totalPresR = StateRight_V(p_) + PeRight + PwaveRight

      ! Rotate vector variables into a coordinate system orthogonal to the face
      call rotate_state_vectors

      ! Normal velocity component
      ! UnStarL = UnStarR = UnStar
      UnStar = (RhoR*UnR*(sR-UnR) - RhoL*UnL*(sL-UnL) +totalPresL -totalPresR)&
           /(RhoR*(sR-UnR) - RhoL*(sL-UnL))

      if(sL >= 0.)then
         Flux_V    = FluxLeft_V
         Unormal_I = UnLeft_I
      elseif(sL < 0. .and. UnStar >= 0.)then
         StateStarCons_V = StateLeftCons_V
         StateStarCons_V(RhoUx_+iDimFace-x_) = StateLeftCons_V(Rho_)*UnStar
         StateStarCons_V(Energy_) = StateStarCons_V(Energy_) &
              + (UnStar-UnL)*(RhoL*UnStar + totalPresL/(sL-UnL))
         StateStarCons_V = StateStarCons_V*(sL-UnL)/(sL-UnStar)
         Flux_V    = FluxLeft_V  + sL*(StateStarCons_V - StateLeftCons_V)
         Unormal_I = UnStar
      elseif(UnStar < 0. .and. sR >= 0.)then
         StateStarCons_V = StateRightCons_V
         StateStarCons_V(RhoUx_+iDimFace-x_) = StateRightCons_V(Rho_)*UnStar
         StateStarCons_V(Energy_) = StateStarCons_V(Energy_) &
              + (UnStar-UnR)*(RhoR*UnStar + totalPresR/(sR-UnR))
         StateStarCons_V = StateStarCons_V*(sR-UnR)/(sR-UnStar)
         Flux_V    = FluxRight_V + sR*(StateStarCons_V - StateRightCons_V)
         Unormal_I = UnStar
      else
         Flux_V    = FluxRight_V
         Unormal_I = UnRight_I
      endif

    end subroutine hllc_flux

    !==========================================================================

    subroutine write_test_info
      integer :: iVar
      !--------------------------------------------------------------------
      write(*,'(1x,4(a,i4))')'Hat state for dir=',iDimFace,&
           ' at I=',iFace,' J=',jFace,' K=',kFace
      write(*,*)'rho=',0.5*(StateLeft_V(Rho_)+StateRight_V(Rho_))
      write(*,*)'Un =',0.5*(StateLeft_V(U_+iDimFace)+StateRight_V(U_+iDimFace))
      write(*,*)'P  =',0.5*(StateLeft_V(P_)+StateRight_V(P_))
      if(UseB)then
         write(*,*)'B  =', &
              0.5*(StateLeft_V(Bx_:Bz_) + StateRight_V(Bx_:Bz_)) + (/B0x,B0y,B0z/)
         write(*,*)'BB =', &
              sum( (0.5*(StateLeft_V(Bx_:Bz_) + StateRight_V(Bx_:Bz_)) &
              + (/B0x,B0y,B0z/))**2)
      end if
      write(*,'(1x,4(a,i4))') 'Fluxes for dir    =',iDimFace,&
           ' at I=',iFace,' J=',jFace,' K=',kFace

      write(*,'(1x,4(a,i4),a,es13.5)') 'Flux*Area for dir =',iDimFace,&
           ' at I=',iFace,' J=',jFace,' K=',kFace,' Area=',Area

      write(*,*)'Eigenvalue_maxabs=', Cmax
      write(*,*)'CmaxDt           =', CmaxDt
      do iVar = 1, nVar + nFluid
         write(*,'(a,a8,5es13.5)') 'Var,F,F_L,F_R,dU,c*dU/2=',&
              NameVar_V(iVar),&
              Flux_V(iVar), &
              FluxLeft_V(iVar)*Area, &
              FluxRight_V(iVar)*Area,&
              StateRightCons_V(iVar)-StateLeftCons_V(iVar),&
              0.5*Cmax*(StateRightCons_V(iVar)-StateLeftCons_V(iVar))*Area
      end do

    end subroutine write_test_info

  end subroutine get_numerical_flux

  !===========================================================================

  subroutine get_physical_flux(State_V, B0x, B0y, B0z, &
       StateCons_V, Flux_V, Un_I, En, Pe, Pwave)

    use ModMultiFluid
    use ModMain,     ONLY: UseHyperbolicDivb, SpeedHyp, UseResistivePlanet
    use ModPhysics,  ONLY: GammaMinus1, GammaElectronMinus1, GammaElectron
    use ModAdvance,  ONLY: UseElectronPressure, UseElectronEntropy
    use ModWaves
    use BATL_size,   ONLY: nDim
    use ModGeometry, ONLY: r_BLK

    real,    intent(in) :: State_V(nVar)       ! input primitive state
    real,    intent(in) :: B0x, B0y, B0z       ! B0
    real,    intent(out):: StateCons_V(nFlux)  !conservative states with energy
    real,    intent(out):: Flux_V(nFlux)       ! fluxes for all states
    real,    intent(out):: Un_I(nFluid+1)      ! normal velocities
    real,    intent(out):: En                  ! normal electric field
    real,    intent(out):: Pe                  ! electron pressure for multiion
    real,    intent(out):: Pwave               ! wave pressure for multiion

    real:: Hyp, Bx, By, Bz, FullBx, FullBy, FullBz, Bn, B0n, FullBn, Un, HallUn
    real:: FluxBx, FluxBy, FluxBz
    real:: FluxViscoX, FluxViscoY, FluxViscoZ

    integer:: iVar

    character(len=*), parameter:: NameSub = 'get_physical_flux'
    !--------------------------------------------------------------------------
    ! Calculate conservative state
    StateCons_V(1:nVar)  = State_V

    ! Make sure normal electric field is initialized
    En = 0.0

    ! Initialize wave pressure
    Pwave = 0.0

    if(UseMultiIon)then
       ! Pe has to be returned for multiion only
       if(UseElectronPressure)then
          Pe = State_V(Pe_)
       elseif(IsMhd)then
          Pe = State_V(p_)*PePerPtotal
       else
          Pe = sum(State_V(iPIon_I))*ElectronPressureRatio
       end if

       if(UseWavePressure)then
          if(UseWavePressureLtd)then
             Pwave = (GammaWave - 1)*State_V(Ew_)
          else
             Pwave = (GammaWave - 1)*sum(State_V(WaveFirst_:WaveLast_))
          end if
       end if
    else
       Pe = 0.0
    end if

    ! Set magnetic variables
    if(UseB)then
       Bx = State_V(Bx_)
       By = State_V(By_)
       Bz = State_V(Bz_)
       FullBx  = Bx + B0x
       FullBy  = By + B0y
       FullBz  = Bz + B0z
       Bn      = Bx*NormalX  + By*NormalY  + Bz*NormalZ
       B0n     = B0x*NormalX + B0y*NormalY + B0z*NormalZ
       FullBn  = B0n + Bn
    end if

    ! Make sure this is initialized
    HallUn = 0.0

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       if(iFluid == 1 .and. IsMhd)then
          ! Calculate MHD flux for first fluid
          if(UseBoris)then
             call get_boris_flux
          else
             call get_mhd_flux
          end if
       elseif(DoBurgers) then
          call get_burgers_flux
       else
          ! If there is no MHD fluid, calculate fluxes for magnetic field
          ! (and electric field) together with hydro fluxes for the first fluid
          if(iFluid == 1 .and. UseB)then
             if(UseEfield)then
                call get_electro_magnetic_flux
             else
                call get_magnetic_flux
             end if
          end if

          ! Calculate HD flux for individual ion and neutral fluids
          call get_hd_flux
       end if

       if(UseResistivePlanet .and. iFluid == 1)then
          ! Do not evolve magnetic field inside the body
          if(r_BLK(iLeft,jLeft,kLeft,iBlockFace)  < 1.0 .and. &
               r_BLK(iRight,jRight,kRight,iBlockFace) < 1.0) &
               Flux_V(Bx_:Bz_) = 0.0
       end if
       
       ! Store normal velocity (needed for source terms with div U)
       Un_I(iFluid) = Un

    end do

    ! The extra fluxes should be added at the same time as fluid 1 fluxes
    if(iFluidMin /= 1) RETURN

    ! Scalars advect with the first fluid's velocity
    do iVar = ScalarFirst_, ScalarLast_
       Flux_V(iVar) = Un_I(1)*State_V(iVar)
    end do

    ! Set flux for electron pressure
    if(UseElectronPressure)then
       if(UseElectronEntropy) &
            StateCons_V(Pe_) = State_V(Pe_)**(1/GammaElectron)
       Flux_V(Pe_) = HallUn*StateCons_V(Pe_)
    end if

    if(Ehot_ > 1) Flux_V(Ehot_) = HallUn*State_V(Ehot_)

    if(UseAlfvenWaves)then
       AlfvenSpeed = FullBn/sqrt(State_V(iRhoIon_I(1)))
       if(SignB_>1 .and. DoThinCurrentSheet)then
          if(State_V(SignB_) < 0.0) AlfvenSpeed = -AlfvenSpeed
       end if

       do iVar = AlfvenWavePlusFirst_, AlfvenWavePlusLast_
          Flux_V(iVar) = (Un_I(IonFirst_) + AlfvenSpeed)*State_V(iVar)
       end do

       do iVar = AlfvenWaveMinusFirst_, AlfvenWaveMinusLast_
          Flux_V(iVar) = (Un_I(IonFirst_) - AlfvenSpeed)*State_V(iVar)
       end do
    end if

    if(ViscoCoeff > 0.0 ) then
       do iFluid = 1, nFluid
          call select_fluid
          FluxViscoX     = sum(Normal_D(1:nDim)*Visco_DDI(:,x_,iFluid))
          Flux_V(iRhoUx) = Flux_V(iRhoUx) - State_V(iRho)*FluxViscoX
          Flux_V(Energy_)= Flux_V(Energy_) - &
               State_V(iRho)*State_V(iUx)*FluxViscoX
          if(nDim == 1) CYCLE
          FluxViscoY     = sum(Normal_D(1:nDim)*Visco_DDI(:,y_,iFluid))
          Flux_V(iRhoUy) = Flux_V(iRhoUy) - &
               State_V(iRho)*FluxViscoY
          Flux_V(Energy_)= Flux_V(Energy_) - &
               State_V(iRho)*State_V(iUy)*FluxViscoY
          if(nDim == 2) CYCLE
          FluxViscoZ     = sum(Normal_D(1:nDim)*Visco_DDI(:,z_,iFluid))
          Flux_V(iRhoUz) = Flux_V(iRhoUz) - State_V(iRho)*FluxViscoZ
          Flux_V(Energy_)= Flux_V(Energy_) - &
               State_V(iRho)*State_V(iUz)*FluxViscoZ
       end do
    end if

    if(UseB) then
       ! These terms are common for the induction equation
       ! If the first fluid is the total fluid, 
       ! the total energy density is also updated
       if(Eta > 0.0)then
          ! Add curl Eta.J to induction equation
          FluxBx = NormalY*EtaJz - NormalZ*EtaJy
          FluxBy = NormalZ*EtaJx - NormalX*EtaJz
          FluxBz = NormalX*EtaJy - NormalY*EtaJx

          Flux_V(Bx_) = Flux_V(Bx_) + FluxBx
          Flux_V(By_) = Flux_V(By_) + FluxBy
          Flux_V(Bz_) = Flux_V(Bz_) + FluxBz

          ! add B.dB/dt term to energy equation
          if(IsMhd) Flux_V(Energy_) = Flux_V(Energy_) &
               + Bx*FluxBx + By*FluxBy + Bz*FluxBz
       end if

       if(UseHallGradPe)then
          ! Add curl (-grad Pe/n e) to induction equation
          FluxBx = - (NormalY*GradZPeNe - NormalZ*GradYPeNe)
          FluxBy = - (NormalZ*GradXPeNe - NormalX*GradZPeNe)
          FluxBz = - (NormalX*GradYPeNe - NormalY*GradXPeNe)

          Flux_V(Bx_) = Flux_V(Bx_) + FluxBx
          Flux_V(By_) = Flux_V(By_) + FluxBy
          Flux_V(Bz_) = Flux_V(Bz_) + FluxBz

          ! add B.dB/dt term to energy equation
          if(IsMhd) Flux_V(Energy_) = Flux_V(Energy_) &
               + Bx*FluxBx + By*FluxBy + Bz*FluxBz
       end if

       if(UseHyperbolicDivb)then
          Hyp  = State_V(Hyp_)

          Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) + SpeedHyp*Normal_D*Hyp
          Flux_V(Hyp_)    = SpeedHyp*Bn

          if(IsMhd) Flux_V(Energy_) = Flux_V(Energy_) + SpeedHyp*Bn*Hyp
       elseif(Hyp_ > 1)then
          Flux_V(Hyp_) = 0.0
       end if
    end if

    if(DoRadDiffusion) Flux_V(Erad_) = Flux_V(Erad_) + EradFlux
    if(DoHeatConduction)then
       if(UseElectronPressure)then
          if(UseElectronEntropy) call stop_mpi(NameSub// &
               ' heat conduction for electron entropy is not implemented')
          Flux_V(Pe_) = Flux_V(Pe_) + GammaElectronMinus1*HeatFlux
       else
          Flux_V(p_) = Flux_V(p_) + GammaMinus1*HeatFlux
          Flux_V(Energy_) = Flux_V(Energy_) + HeatFlux
       end if
    end if
    if(DoIonHeatConduction)then
       Flux_V(p_) = Flux_V(p_) + GammaMinus1*IonHeatFlux
       Flux_V(Energy_) = Flux_V(Energy_) + IonHeatFlux
    end if

    ! Set the normal electron velocity used for Hall MHD and/or 
    ! the electron pressure source term
    Un_I(eFluid_) = HallUn

  contains

    subroutine get_boris_flux

      use ModPhysics, ONLY: InvGammaMinus1, Inv_C2light, InvClight
      use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e
      real :: B2, FullB2, pTotal, pTotal2, uDotB, DpPerB
      real :: Ex, Ey, Ez, E2Half
      !-----------------------------------------------------------------------

      ! Extract primitive variables
      Rho     = State_V(Rho_)
      Ux      = State_V(Ux_)
      Uy      = State_V(Uy_)
      Uz      = State_V(Uz_)
      p       = State_V(p_)

      B2      = Bx**2 + By**2 + Bz**2

      ! Electric field divided by speed of light: 
      ! E= - U x B / c = (B x U)/c
      Ex      = (FullBy*Uz - FullBz*Uy) * InvClight
      Ey      = (FullBz*Ux - FullBx*Uz) * InvClight
      Ez      = (FullBx*Uy - FullBy*Ux) * InvClight

      ! Electric field squared/c^2
      E2Half  = 0.5*(Ex**2 + Ey**2 + Ez**2)

      ! Calculate energy and total pressure
      e = InvGammaMinus1*p + 0.5*(Rho*(Ux**2 + Uy**2 + Uz**2) + B2)

      pTotal  = p + 0.5*B2 + B0x*Bx + B0y*By + B0z*Bz

      if(UseElectronPressure) pTotal = pTotal + State_V(Pe_)

      if(UseWavePressure)then
         if(UseWavePressureLtd)then
            pTotal = pTotal + (GammaWave-1)*State_V(Ew_)
         else
            pTotal = pTotal + (GammaWave-1)*sum(State_V(WaveFirst_:WaveLast_))
         end if
      end if

      ! pTotal = pperp + bb/2 = 3/2*p - 1/2*ppar + bb/2 
      !        = p + bb/2 + (p - ppar)/2
      if(UseAnisoPressure) pTotal = pTotal + 0.5*(p - State_V(Ppar_))

      pTotal2 = pTotal + E2Half

      ! The full momentum contains the ExB/c^2 term:
      ! rhoU_Boris = rhoU - ((U x B) x B)/c^2 = rhoU + (U B^2 - B U.B)/c^2
      uDotB   = Ux*FullBx + Uy*FullBy + Uz*FullBz
      FullB2  = FullBx**2 + FullBy**2 + FullBz**2
      StateCons_V(RhoUx_)  = Rho*Ux + (Ux*FullB2 - FullBx*uDotB)*inv_c2LIGHT
      StateCons_V(RhoUy_)  = Rho*Uy + (Uy*FullB2 - FullBy*uDotB)*inv_c2LIGHT
      StateCons_V(RhoUz_)  = Rho*Uz + (Uz*FullB2 - FullBz*uDotB)*inv_c2LIGHT

      ! The full energy contains the electric field energy
      StateCons_V(Energy_) = e + E2Half

      ! Normal direction
      Un     = Ux*NormalX + Uy*NormalY + Uz*NormalZ
      En     = Ex*NormalX + Ey*NormalY + Ez*NormalZ

      ! f_i[rho] = rho*u_i
      Flux_V(Rho_)   = Rho*Un

      ! f_i[rhou_k] = u_i*u_k*rho - b_k*b_i - B0_k*b_i - B0_i*b_k - E_i*E_k
      !          +n_i*[p + B0_j*b_j + 0.5*(b_j*b_j + E_j*E_j)]
      Flux_V(RhoUx_) = Un*Rho*Ux - Bn*FullBx - B0n*Bx - En*Ex + pTotal2*Normalx
      Flux_V(RhoUy_) = Un*Rho*Uy - Bn*FullBy - B0n*By - En*Ey + pTotal2*Normaly
      Flux_V(RhoUz_) = Un*Rho*Uz - Bn*FullBz - B0n*Bz - En*Ez + pTotal2*Normalz

      ! f_i[b_k]=u_i*(b_k+B0_k) - u_k*(b_i+B0_i)
      Flux_V(Bx_) = Un*FullBx - Ux*FullBn
      Flux_V(By_) = Un*FullBy - Uy*FullBn
      Flux_V(Bz_) = Un*FullBz - Uz*FullBn

      ! f_i[p]=u_i*p
      Flux_V(p_)  = Un*p

      ! f_i[e]=(u_i*(ptotal+e+(b_k*B0_k))-(b_i+B0_i)*(b_k*u_k))
      Flux_V(Energy_) = &
           Un*(pTotal + e) - FullBn*(Ux*Bx + Uy*By + Uz*Bz)

      if(UseAnisoPressure)then
         ! f_i[rhou_k] = f_i[rho_k] + (ppar - pperp)bb for anisopressure
         ! ppar - pperp = ppar - (3*p - ppar)/2 = 3/2*(ppar - p)
         DpPerB = 1.5*(State_V(Ppar_) - p)*FullBn/max(1e-30, FullB2)
         Flux_V(RhoUx_) = Flux_V(RhoUx_) + FullBx*DpPerB
         Flux_V(RhoUy_) = Flux_V(RhoUy_) + FullBy*DpPerB
         Flux_V(RhoUz_) = Flux_V(RhoUz_) + FullBz*DpPerB
         ! f_i[Ppar] = u_i*Ppar
         Flux_V(Ppar_)  = Un*State_V(Ppar_)
         Flux_V(Energy_) = Flux_V(Energy_) &
              + DpPerB*(Ux*FullBx + Uy*FullBy + Uz*FullBz)
      end if

      HallUn = Un

    end subroutine get_boris_flux

    !==========================================================================

    subroutine get_mhd_flux

      use ModPhysics, ONLY: InvGammaMinus1, inv_c2LIGHT
      use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e
      real :: HallUx, HallUy, HallUz, InvRho
      real :: B2, B0B1, FullB2, pTotal, DpPerB
      real :: Gamma2

      real :: InvElectronDens, UxPlus, UyPlus, UzPlus, UnPlus
      real, dimension(nIonFluid) :: ChargeDens_I, Ux_I, Uy_I, Uz_I, RhoUn_I
      !-----------------------------------------------------------------------

      ! Extract primitive variables
      Rho     = State_V(Rho_)
      Ux      = State_V(Ux_)
      Uy      = State_V(Uy_)
      Uz      = State_V(Uz_)
      p       = State_V(p_)

      B2      = Bx**2 + By**2 + Bz**2

      ! Calculate energy
      e = InvGammaMinus1*p + 0.5*(Rho*(Ux**2 + Uy**2 + Uz**2) + B2)

      ! Calculate some intermediate values for flux calculations
      B0B1    = B0x*Bx + B0y*By + B0z*Bz
      pTotal  = p + 0.5*B2 + B0B1

      if(UseElectronPressure) pTotal = pTotal + State_V(Pe_)

      if(UseWavePressure)then
         if(UseWavePressureLtd)then
            pTotal = pTotal + (GammaWave-1)*State_V(Ew_)
         else
            pTotal = pTotal + (GammaWave-1)*sum(State_V(WaveFirst_:WaveLast_))
         end if
      end if

      ! pTotal = pperp + bb/2 = 3/2*p - 1/2*ppar + bb/2 
      !        = p + bb/2 + (p - ppar)/2
      if(UseAnisoPressure) pTotal = pTotal + 0.5*(p - State_V(Ppar_))

      ! Calculate conservative state
      StateCons_V(RhoUx_)  = Rho*Ux
      StateCons_V(RhoUy_)  = Rho*Uy
      StateCons_V(RhoUz_)  = Rho*Uz
      StateCons_V(Energy_) = e

      ! Normal direction
      Un     = Ux*NormalX  + Uy*NormalY  + Uz*NormalZ

      if(UseMultiIon)then
         ! Calculate charge density averaged velocity U*Plus

         ! calculate charge densities
         ChargeDens_I    = ChargePerMass_I*State_V(iRhoIon_I)
         InvElectronDens = 1.0/sum(ChargeDens_I)

         Ux_I  = State_V(iUxIon_I)
         Uy_I  = State_V(iUyIon_I)
         Uz_I  = State_V(iUzIon_I)

         ! calculate the average positive charge velocity
         UxPlus = InvElectronDens*sum(ChargeDens_I*Ux_I)
         UyPlus = InvElectronDens*sum(ChargeDens_I*Uy_I)
         UzPlus = InvElectronDens*sum(ChargeDens_I*Uz_I)

         UnPlus = UxPlus*NormalX + UyPlus*NormalY + UzPlus*NormalZ
      else
         ! For single ion fluid the mass and charge average is the same
         UxPlus = Ux
         UyPlus = Uy
         UzPlus = Uz
         UnPlus = Un
      end if

      ! f_n[rho] = Rho*U_i
      Flux_V(Rho_) = Rho*Un

      if(UseMultiIon)then
         ! Add up the (rho u u) diads of the ion fluids:
         ! f_n[rhou_k] = sum_s(rho_s*u_n,s*u_k,s) 
         !               - b_n*(b_k + B0_k) - B0_n*b_k + Ptotal*n_k

         RhoUn_I = State_V(iRhoIon_I) &
              *(Ux_I*NormalX + Uy_I*NormalY + Uz_I*NormalZ)

         Flux_V(RhoUx_)=sum(RhoUn_I*Ux_I) - Bn*FullBx - B0n*Bx + pTotal*NormalX
         Flux_V(RhoUy_)=sum(RhoUn_I*Uy_I) - Bn*FullBy - B0n*By + pTotal*NormalY
         Flux_V(RhoUz_)=sum(RhoUn_I*Uz_I) - Bn*FullBz - B0n*Bz + pTotal*NormalZ
      else
         ! f_n[rhou_k] = u_n*u_k*rho - b_n*(b_k + B0_k) - B0_n*b_k + Ptotal*n_k
         Flux_V(RhoUx_) = Un*Rho*Ux - Bn*FullBx - B0n*Bx + pTotal*NormalX
         Flux_V(RhoUy_) = Un*Rho*Uy - Bn*FullBy - B0n*By + pTotal*NormalY
         Flux_V(RhoUz_) = Un*Rho*Uz - Bn*FullBz - B0n*Bz + pTotal*NormalZ
      end if

      ! f_n[b_k] = u_n*(b_k + B0_k) - u_k*(b_n + B0_n)
      if(HallCoeff > 0.0)then
         ! Note that the HallCoeff contains the ion mass/charge.
         InvRho = 1.0/Rho
         HallUx = UxPlus - HallJx*InvRho
         HallUy = UyPlus - HallJy*InvRho
         HallUz = UzPlus - HallJz*InvRho

         HallUn = NormalX*HallUx + NormalY*HallUy + NormalZ*HallUz
      else
         HallUn = UnPlus
      end if

      if(HallCoeff > 0.0 .and. DoHallInduction)then
         Flux_V(Bx_) = HallUn*FullBx - HallUx*FullBn
         Flux_V(By_) = HallUn*FullBy - HallUy*FullBn
         Flux_V(Bz_) = HallUn*FullBz - HallUz*FullBn
      else
         Flux_V(Bx_) = UnPlus*FullBx - UxPlus*FullBn
         Flux_V(By_) = UnPlus*FullBy - UyPlus*FullBn
         Flux_V(Bz_) = UnPlus*FullBz - UzPlus*FullBn
      end if

      ! f_n[p] = u_n*p
      Flux_V(p_) = Un*p

      if(UseAnisoPressure)then
         ! f_i[rhou_k] = f_i[rho_k] + (ppar - pperp)bb for anisopressure
         ! ppar - pperp = ppar - (3*p - ppar)/2 = 3/2*(ppar - p)
         FullB2 = FullBx**2 + FullBy**2 + FullBz**2
         DpPerB = 1.5*(State_V(Ppar_) - p)*FullBn/max(1e-30, FullB2)

         Flux_V(RhoUx_) = Flux_V(RhoUx_) + FullBx*DpPerB
         Flux_V(RhoUy_) = Flux_V(RhoUy_) + FullBy*DpPerB
         Flux_V(RhoUz_) = Flux_V(RhoUz_) + FullBz*DpPerB

         ! f_i[Ppar] = u_i*Ppar
         Flux_V(Ppar_)  = Un*State_V(Ppar_)
      end if

      ! de/dt = -div(E x B) = -B.curl E + E.curl B = B.dB/dt + E.j
      ! For semi-implicit MHD the energy gets updated at the same time 
      ! as the magnetic field, so de/dt = B.dB/dt = -B.curlE is taken
      ! into account. The other term is E.j = 0, since the 
      ! Hall related electric field is j x B/ne so orthogonal to j.
      if(HallCoeff > 0.0 .and. DoHallInduction) then
         Flux_V(Energy_) = &
              Un*(pTotal + e) &
              - FullBn*(HallUx*Bx + HallUy*By + HallUz*Bz)  &
              + (HallUn - Un)*(B2 + B0B1)
      elseif(UseMultiIon)then
         Flux_V(Energy_) = &
              Un*(pTotal + e) &
              - FullBn*(UxPlus*Bx + UyPlus*By + UzPlus*Bz)  &
              + (UnPlus - Un)*(B2 + B0B1)
      elseif(UseAnisoPressure)then
         ! can only work for single fluid without Hall 
         Flux_V(Energy_) = &
              Un*(pTotal + e) &
              + DpPerB*(Ux*FullBx + Uy*FullBy + Uz*FullBz) &
              - FullBn*(Ux*Bx + Uy*By + Uz*Bz)
      else
         ! f_i[e]=(u_i*(ptotal + e + (b_k*B0_k)) - (b_i+B0_i)*(b_k*u_k))
         Flux_V(Energy_) = &
              Un*(pTotal + e) - FullBn*(Ux*Bx + Uy*By + Uz*Bz)     
      end if

      ! Correct energy flux, so that the electron contribution to the energy
      ! flux is U_e*p_e. We add (U_e-U_ion)*p_e.
      if(UseElectronPressure .and. nIonFluid == 1 .and. iFluid == 1)then
         if(HallCoeff > 0) &
              Flux_V(Energy_) = Flux_V(Energy_) + (HallUn - Un)*State_V(Pe_)
      end if

      if(UseBorisSimple)then
         ! Correct the momentum using the (1+VA2/c^2)
         Gamma2 = 1.0 + (FullBx**2 + FullBy**2 + FullBz**2)/Rho*inv_c2LIGHT
         StateCons_V(RhoUx_:RhoUz_) = StateCons_V(RhoUx_:RhoUz_)*Gamma2
      end if

      if(DoTestCell)then
         write(*,*)'ChargeDens_I,InvRho         =', &
              ChargeDens_I, InvRho
         write(*,*)'UxyzPlus  =',UxPlus,UyPlus,UzPlus
         write(*,*)'HallUxyz  =',HallUx,HallUy,HallUz
         write(*,*)'FullBxyz  =',FullBx,FullBy,FullBz
         write(*,*)'B0x,y,z   =',B0x,B0y,B0z
         write(*,*)'Flux(Bxyz)=',Flux_V(Bx_:Bz_)
      end if

    end subroutine get_mhd_flux

    !==========================================================================

    subroutine get_magnetic_flux

      ! Calculate magnetic flux for multi-ion equations 
      ! without a global ion fluid

      real :: ChargeDens_I(nIonFluid), InvElectronDens
      real :: UxPlus, UyPlus, UzPlus, UnPlus
      real :: HallUx, HallUy, HallUz, InvRho
      !-----------------------------------------------------------------------

      ! calculate number densities
      ChargeDens_I    = ChargePerMass_I*State_V(iRhoIon_I)
      InvElectronDens = 1.0/sum(ChargeDens_I)

      ! calculate positive charge velocity
      UxPlus = InvElectronDens*sum(ChargeDens_I*State_V(iUxIon_I))
      UyPlus = InvElectronDens*sum(ChargeDens_I*State_V(iUyIon_I))
      UzPlus = InvElectronDens*sum(ChargeDens_I*State_V(iUzIon_I))
      UnPlus = UxPlus*NormalX + UyPlus*NormalY + UzPlus*NormalZ

      if(HallCoeff > 0.0)then
         ! The ion mass per charge that is contained in HallCoef (and HallJ*) 
         ! is normalized to be divided with the total mass density.
         InvRho = 1/sum(State_V(iRhoIon_I))
         HallUx = UxPlus - HallJx*InvRho
         HallUy = UyPlus - HallJy*InvRho
         HallUz = UzPlus - HallJz*InvRho
         HallUn = NormalX*HallUx + NormalY*HallUy + NormalZ*HallUz
      else
         HallUn = UnPlus
      end if

      if(HallCoeff > 0.0 .and. DoHallInduction)then
         Flux_V(Bx_) = HallUn*FullBx - HallUx*FullBn
         Flux_V(By_) = HallUn*FullBy - HallUy*FullBn
         Flux_V(Bz_) = HallUn*FullBz - HallUz*FullBn
      else
         Flux_V(Bx_) = UnPlus*FullBx - UxPlus*FullBn
         Flux_V(By_) = UnPlus*FullBy - UyPlus*FullBn
         Flux_V(Bz_) = UnPlus*FullBz - UzPlus*FullBn
      end if

      if(DoTestCell)then
         write(*,*)'ChargeDens_I,InvElectronDens,InvRho=', &
              ChargeDens_I, InvElectronDens,InvRho
         write(*,*)'UxyzPlus  =',UxPlus,UyPlus,UzPlus
         write(*,*)'HallUxyz  =',HallUx,HallUy,HallUz
         write(*,*)'FullBxyz  =',FullBx,FullBy,FullBz
         write(*,*)'B0x,y,z   =',B0x,B0y,B0z
         write(*,*)'Flux(Bxyz)=',Flux_V(Bx_:Bz_)
      end if

    end subroutine get_magnetic_flux

    !==========================================================================

    subroutine get_electro_magnetic_flux

      use ModPhysics, ONLY: Clight, C2light

      real:: Ex, Ey, Ez
      !----------------------------------------------------------------------
      Ex = State_V(Ex_); Ey = State_V(Ey_); Ez = State_V(Ez_)

      ! dB/dt + div F = 0
      ! div F = curl E = sum(Normal x E)/Volume
      Flux_V(Bx_) = NormalY*Ez - NormalZ*Ey
      Flux_V(By_) = NormalZ*Ex - NormalX*Ez
      Flux_V(Bz_) = NormalX*Ey - NormalY*Ex

      ! dE/dt + c^2(J - curl B) = 0   (curl B0 is assumed to be zero for now)
      Flux_V(Ex_) = -C2light*(NormalY*Bz - NormalZ*By)
      Flux_V(Ey_) = -C2light*(NormalZ*Bx - NormalX*Bz)
      Flux_V(Ez_) = -C2light*(NormalX*By - NormalY*Bx)

      ! dE/dt + c*grad PhiE ...
      Flux_V(Ex_:Ez_) = Flux_V(Ex_:Ez_) + Clight*Normal_D*State_V(HypE_)

      ! Flux part of dPhiE/dt + c*(div E - chargedensity/eps0) 
      Flux_V(HypE_) = Clight*(Ex*NormalX  + Ey*NormalY  + Ez*NormalZ)

      if(DoTestCell)then
         write(*,'(a,3es13.5)')'ChargeDens_I    =', &
              ChargePerMass_I*State_V(iRhoIon_I)
         write(*,'(a,3es13.5)')'Normal_D        =', Normal_D
         write(*,'(a,3es13.5)')'Bx,By,Bz        =', Bx,By,Bz
         write(*,'(a,3es13.5)')'Ex,Ey,Ez        =', Ex,Ey,Ez
         write(*,'(a,3es13.5)')'Flux_V(Bx_:Bz_) =', Flux_V(Bx_:Bz_)
         write(*,'(a,3es13.5)')'Flux_V(Ex_:Ez_) =', Flux_V(Ex_:Ez_)
         write(*,'(a,3es13.5)')'State_V(HypE_)  =', State_V(HypE_)
         write(*,'(a,3es13.5)')'Flux_V(HypE_)   =', Flux_V(HypE_)
      end if

    end subroutine get_electro_magnetic_flux
    !==========================================================================
    subroutine get_hd_flux

      use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure
      use ModPhysics, ONLY: InvGammaMinus1_I
      use ModWaves

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e, RhoUn, pTotal
      real :: DpPerB, FullB2
      !-----------------------------------------------------------------------
      ! Extract primitive variables
      Rho     = State_V(iRho)
      Ux      = State_V(iUx)
      Uy      = State_V(iUy)
      Uz      = State_V(iUz)
      p       = State_V(iP)

      ! Calculate energy
      e = InvGammaMinus1_I(iFluid)*p + 0.5*Rho*(Ux**2 + Uy**2 + Uz**2)

      pTotal = p

      if(nIonFluid == 1 .and. iFluid == 1)then
         if(UseElectronPressure) pTotal = pTotal + State_V(Pe_)

         if(UseWavePressure) &
              pTotal = pTotal +(GammaWave-1)*sum(State_V(WaveFirst_:WaveLast_))
      end if

      ! pTotal = pperp = 3/2*p - 1/2*ppar = p + (p - ppar)/2
      if(UseAnisoPressure .and. IsIon_I(iFluid)) &
           pTotal = pTotal + 0.5*(p - State_V(iPpar))

      ! Calculate conservative state
      StateCons_V(iRhoUx)  = Rho*Ux
      StateCons_V(iRhoUy)  = Rho*Uy
      StateCons_V(iRhoUz)  = Rho*Uz
      StateCons_V(iEnergy) = e

      ! Normal velocity
      Un     = Ux*NormalX  + Uy*NormalY  + Uz*NormalZ
      RhoUn  = Rho*Un

      ! f_i[rho] = rho*u_i
      Flux_V(iRho)   = RhoUn

      ! f_i[rhou_k] = u_i*rho*u_k + n_i*[ptotal]
      Flux_V(iRhoUx) = RhoUn*Ux + pTotal*NormalX
      Flux_V(iRhoUy) = RhoUn*Uy + pTotal*NormalY
      Flux_V(iRhoUz) = RhoUn*Uz + pTotal*NormalZ

      ! f_i[p] = u_i*p
      Flux_V(iP)  = Un*p

      Flux_V(iEnergy) = Un*(pTotal + e)

      if(UseAnisoPressure .and. IsIon_I(iFluid))then
         ! f_i[rhou_k] = f_i[rho_k] + (ppar - pperp)bb for anisopressure
         ! ppar - pperp = ppar - (3*p - ppar)/2 = 3/2*(ppar - p)
         FullB2 = FullBx**2 + FullBy**2 + FullBz**2
         DpPerB = 1.5*(State_V(iPpar) - p)*FullBn/max(1e-30, FullB2)

         Flux_V(iRhoUx) = Flux_V(iRhoUx) + FullBx*DpPerB
         Flux_V(iRhoUy) = Flux_V(iRhoUy) + FullBy*DpPerB
         Flux_V(iRhoUz) = Flux_V(iRhoUz) + FullBz*DpPerB

         ! f_i[Ppar] = u_i*Ppar
         Flux_V(iPpar)  = Un*State_V(iPpar)

         Flux_V(iEnergy) = Flux_V(iEnergy) &
              + DpPerB*(Ux*FullBx + Uy*FullBy + Uz*FullBz)
      end if

      ! Needed for adiabatic source term for electron pressure
      if(iFluid == 1 .and. .not.UseB) HallUn = Un

    end subroutine get_hd_flux

    !======================================================================
    subroutine get_burgers_flux
      Flux_V = 0.0
      Flux_V(iRho) = 0.5*State_V(iRho)**2
    end subroutine get_burgers_flux
    !=========================================================================

  end subroutine get_physical_flux

  !===========================================================================
  subroutine calc_cell_flux(iBlock)

    ! Calculate cell centered fluxes including ghost cells

    use ModAdvance, ONLY: nVar, State_VGB, FluxLeft_VGD, FluxRight_VGD
    use ModMultiFluid, ONLY: nFluid, iRho, iUx, iUz, iUx_I, iUz_I
    use ModMain, ONLY: UseDtFixed, Cfl, Dt
    use BATL_lib, ONLY: nDim

    integer, intent(in):: iBlock

    integer:: i, j, k, iDim, iFluid

    real:: Primitive_V(nVar), RhoInv, Flux_V(nFlux)
    real:: CmaxArea, CmaxAll, Cmax_I(nFluid), Conservative_V(nFlux)

    ! These are calculated but not used                                        
    real:: Un_I(nFluid+1), En, Pe, Pwave
    !------------------------------------------------------------------------

    if(.not.allocated(FluxLeft_VGD)) allocate( &
         FluxLeft_VGD(nFlux,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nDim), &
         FluxRight_VGD(nFlux,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nDim))

    UseHallGradPe = .false. !!! HallJx = 0; HallJy = 0; HallJz = 0
    DoTestCell = .false.
    do iDim = 1, nDim
       call set_block_values(iBlock, iDim)

       if(Dt > 0)then
          CmaxAll = CellSize_DB(iDim,iBlock)/Dt
          if(.not.UseDtFixed) CmaxAll = Cfl*CmaxAll
       else
          CmaxAll = 0.0
       end if

       do k = MinK, MaxK; kFace = k
          do j = MinJ, MaxJ; jFace = j
             do i = MinI, MaxI; iFace = i
                !!! call set_cell_values

                ! Get primitive variables used by get_physical_flux
                Primitive_V = State_VGB(:,i,j,k,iBlock)
                do iFluid = 1, nFluid
                   iRho = iRho_I(iFluid)
                   iUx = iUx_I(iFluid)
                   iUz = iUz_I(iFluid)
                   RhoInv = 1/Primitive_V(iRho)
                   Primitive_V(iUx:iUz) = RhoInv*Primitive_V(iUx:iUz)
                end do

                if(UseB0)then
                   B0x = B0_DGB(x_,i,j,k,iBlock)
                   B0y = B0_DGB(x_,i,j,k,iBlock)
                   B0z = B0_DGB(x_,i,j,k,iBlock)
                end if

                ! Get the flux
                call get_physical_flux(Primitive_V, B0x, B0y, B0z, &
                     Conservative_V, Flux_V, Un_I, En, Pe, Pwave)

                FluxLeft_VGD(:,i,j,k,iDim) = &
                     0.5*Area*(Flux_V + CmaxAll*Conservative_V)
                FluxRight_VGD(:,i,j,k,iDim)= &
                     0.5*Area*(Flux_V - CmaxAll*Conservative_V)

                ! Get the maximum speed
                call get_speed_max(Primitive_V, B0x, B0y, B0z, Cmax_I)
                CmaxArea = maxval(Cmax_I)*Area
                select case(iDim)
                case(1)
                   if(i>0 .and. i<=nI+1 .and. j>0 .and. j<=nJ .and. &
                        k>0 .and. k<=nK) &
                        VdtFace_X(i,j,k) = CmaxArea
                case(2)
                   if(i>0 .and. i<=nI .and. j>0 .and. j<=nJ+1 .and. &
                        k>0 .and. k<=nK) &
                        VdtFace_Y(i,j,k) = CmaxArea
                case(3)
                   if(i>0 .and. i<=nI .and. j>0 .and. j<=nJ .and. &
                        k>0 .and. k<=nK+1) &
                        VdtFace_Z(i,j,k) = CmaxArea
                end select

             end do
          end do
       end do
    end do

  end subroutine calc_cell_flux
  !===========================================================================
  subroutine calc_simple_cell_flux(iBlock)

    ! Calculate cell centered fluxes including ghost cells

    use ModImplicit, ONLY: UseSemiHallResist
    use ModAdvance, ONLY: nVar, State_VGB,FluxCenter_VGD
    use ModMultiFluid, ONLY: nFluid, iRho, iUx, iUz, iUx_I, iUz_I
    use BATL_lib, ONLY: nDim, Xi_, Eta_, Zeta_, CellCoef_DDGB

    integer, intent(in):: iBlock

    integer:: i, j, k, iDim, iFluid

    real:: Primitive_V(nVar), RhoInv, Flux_V(nFlux)
    real:: Conservative_V(nFlux)

    ! These are calculated but not used          
    real:: Un_I(nFluid+1), En, Pe, Pwave

    real, allocatable, save:: Flux_VD(:,:)
    integer:: iFlux
    character(len=*), parameter:: NameSub = 'calc_simple_cell_flux'
    !------------------------------------------------------------------------

    if(.not.allocated(FluxCenter_VGD)) allocate( &
         FluxCenter_VGD(nFlux,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nDim))
    if(.not.allocated(Flux_VD)) allocate(Flux_VD(nFlux,nDim))

    HallCoeff     = -1.0
    if(UseHallResist .and. .not.UseSemiHallResist) call stop_mpi(NameSub// &
         ": Hall Resistivity has not been implemented!")
    
    BiermannCoeff = -1.0
    if(UseBiermannBattery) call stop_mpi(NameSub// &
         ": BiermannBattery has not been implemented!")
    
    ViscoCoeff = 0.0
    if(UseViscosity) call stop_mpi(NameSub// &
         ": Viscosity has not been implemented!")

    UseHallGradPe = .false. !!! HallJx = 0; HallJy = 0; HallJz = 0         
    DoTestCell = .false.
    do iDim = 1, nDim
       call set_block_values(iBlock, iDim)

       do k = MinK, MaxK; kFace = k
          do j = MinJ, MaxJ; jFace = j
             do i = MinI, MaxI; iFace = i
                
                Eta       = 0.0
                if(UseResistiveFlux) Eta = Eta_GB(i,j,k,iBlock)

                if(UseClimit)  call stop_mpi(&
                     "Climit has not been implemented for cell fulx!")

                ! Get primitive variables used by get_physical_flux
                Primitive_V = State_VGB(:,i,j,k,iBlock)
                do iFluid = 1, nFluid
                   iRho = iRho_I(iFluid)
                   iUx = iUx_I(iFluid)
                   iUz = iUz_I(iFluid)
                   RhoInv = 1/Primitive_V(iRho)
                   Primitive_V(iUx:iUz) = RhoInv*Primitive_V(iUx:iUz)
                end do

                if(UseB0)then
                   B0x = B0_DGB(x_,i,j,k,iBlock)
                   B0y = B0_DGB(y_,i,j,k,iBlock)
                   B0z = B0_DGB(z_,i,j,k,iBlock)
                end if

                ! Get the flux
                call get_physical_flux(Primitive_V, B0x, B0y, B0z, &
                     Conservative_V, Flux_V, Un_I, En, Pe, Pwave)
                if(.not. UseHighFDGeometry) then
                   FluxCenter_VGD(:,i,j,k,iDim) = Flux_V*Area
                else
                   ! 'Area' is included in the transform coef: 
                   FluxCenter_VGD(:,i,j,k,iDim) = Flux_V                  
                endif
             end do
          end do
       end do
    end do

    if(UseHighFDGeometry) then 
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          do iFlux = 1, nFlux
             Flux_VD(iFlux,:) = FluxCenter_VGD(iFlux,i,j,k,:)
          enddo

          do iFlux = 1, nFlux
             FluxCenter_VGD(iFlux,i,j,k,x_) = &
                  sum(Flux_VD(iFlux,:)*CellCoef_DDGB(Xi_,:,i,j,k,iBlock))
             FluxCenter_VGD(iFlux,i,j,k,y_) = &
                  sum(Flux_VD(iFlux,:)*CellCoef_DDGB(Eta_,:,i,j,k,iBlock))
             if(nK > 1) FluxCenter_VGD(iFlux,i,j,k,z_) = &
                  sum(Flux_VD(iFlux,:)*CellCoef_DDGB(Zeta_,:,i,j,k,iBlock))
          enddo
       enddo; enddo; enddo
    endif
  end subroutine calc_simple_cell_flux
  !===========================================================================

  subroutine get_speed_max(State_V, B0x, B0y, B0z, cMax_I, cLeft_I, cRight_I,&
       UseAwSpeedIn)

    use ModMultiFluid, ONLY: select_fluid, iFluid, iRho, iUx, iUz, iP
    use ModWaves, ONLY: UseWavePressure, UseWavePressureLtd, &
         GammaWave, UseAlfvenWaves
    use ModMain, ONLY: time_accurate, Climit
    use ModPhysics, ONLY: Clight
    use ModAdvance, ONLY: State_VGB

    real,    intent(in) :: State_V(nVar)
    real,    intent(in) :: B0x, B0y, B0z
    real, optional, intent(out) :: Cmax_I(nFluid)   ! max speed relative to lab
    real, optional, intent(out) :: Cleft_I(nFluid)  ! maximum left speed
    real, optional, intent(out) :: Cright_I(nFluid) ! maximum right speed

    logical, optional, intent(in):: UseAwSpeedIn    ! Use AW speed definitions

    logical:: UseAwSpeed

    real :: CmaxDt_I(nFluid)
    real :: UnLeft, UnRight
    integer :: iError = -1
    !--------------------------------------------------------------------------
    UseAwSpeed = .false.
    if(present(UseAwSpeedIn)) UseAwSpeed = UseAwSpeedIn

    do iFluid = iFluidMin, iFluidMax
       call select_fluid

       if(iFluid == 1 .and. UseB)then
          if(UseAwSpeed)then
             ! For AW flux UnLeft_I,UnRight_I 
             ! are already set by get_physical_flux
             UnLeft = minval(UnLeft_I(1:nIonFluid))
             UnRight= maxval(UnRight_I(1:nIonFluid))
          end if

          if(UseBoris .or. (UseEfield .and. nTrueIon ==1))then
             ! In the five- and six-moment case, we should use Boris speed for
             ! numerical diffusion for single ion case. For multi-ion five- and
             ! six-moment we should call the UseBorisSimple correction in 
             ! get_mhd_speed. 
             ! In the five- and six-moment case, we should NOT call Boris 
             ! flux, as the Maxwell's equation already contains the 
             ! displacement current.
             call get_boris_speed             
          else
             call get_mhd_speed
          endif

       elseif(iFluid > 1 .and. iFluid <= IonLast_)then
          if(present(Cleft_I))  Cleft_I(iFluid)  = Cleft_I(1)
          if(present(Cright_I)) Cright_I(iFluid) = Cright_I(1)
          if(present(Cmax_I))   Cmax_I(iFluid)   = Cmax_I(1)
          CmaxDt_I(iFluid) = CmaxDt_I(1)
          CYCLE

       elseif(DoBurgers)then
          call get_burgers_speed
       else
          if(UseAwSpeed)then
             UnLeft = UnLeft_I(iFluid)
             UnRight= UnRight_I(iFluid)
          end if
          call get_hd_speed
       end if

    end do

    if(UseEfield .and. iFluidMin <= nIonFluid)then
       ! The light speed in the five-moment equations should exceed
       ! all the fluid wave speeds. Only Lax-Friedrichs scheme can be
       ! used because the left/right/max wave speeds are the same = Clight.
       if(present(Cmax_I)) then
          CmaxDt_I(iFluidMin:iFluidMax) = &
               max(cLight, Cmax_I(iFluidMin:iFluidMax))
          CmaxDt                        = maxval(CmaxDt_I(iFluidMin:iFluidMax))

          if(maxval(Cmax_I(iFluidMin:iFluidMax)) > Clight*factorClightWarning &
               .and. DoClightWarning) then
             write(*,'(a,i3)')      'dir     =', iDimFace
             write(*,'(a,10es15.6)')'Xyz_DGB =', &
                  Xyz_DGB(:,iFace,jFace,kFace,iBlockFace)
             write(*,'(a,10es15.6)')'Rho     =', &
                  State_VGB(iRho_I,iFace,jFace,kFace,iBLockFace)
             write(*,'(a,10es15.6)')'P       =', &
                  State_VGB(iP_I,iFace,jFace,kFace,iBLockFace)
             write(*,*) 'neighbor point: '
             select case(iDimFace)
             case(x_)
                write(*,'(a,10es15.6)')'Rho     =', &
                     State_VGB(iRho_I,iFace+1,jFace,kFace,iBLockFace)
                write(*,'(a,10es15.6)')'P       =', &
                     State_VGB(iP_I,iFace+1,jFace,kFace,iBLockFace)
             case(y_)
                write(*,'(a,10es15.6)')'Rho     =', &
                     State_VGB(iRho_I,iFace,jFace+1,kFace,iBLockFace)
                write(*,'(a,10es15.6)')'P       =', &
                     State_VGB(iP_I,iFace,jFace+1,kFace,iBLockFace)
             case(z_)
                write(*,'(a,10es15.6)')'Rho     =', &
                     State_VGB(iRho_I,iFace,jFace,kFace+1,iBLockFace)
                write(*,'(a,10es15.6)')'P       =', &
                     State_VGB(iP_I,iFace,jFace,kFace+1,iBLockFace)
             end select
             write(*,*)'cLight         =',cLight
             write(*,*)'maxval(Cmax_I) =',maxval(Cmax_I(iFluidMin:iFluidMax))
             write(*,*)'Cmax_I         =',Cmax_I(iFluidMin:iFluidMax)
          end if

          if (maxval(Cmax_I(iFluidMin:iFluidMax)) > Clight) then
             call error_report( &
                  'get_speed_max: Clihgt is smaller than maxval(Cmax_I)', &
                  maxval(Cmax_I(iFluidMin:iFluidMax)), iError, .true.)
             if(time_accurate) call stop_mpi &
                  ('get_speed_max: Clihgt is smaller than maxval(Cmax_I)')
          end if
       end if

       RETURN
    end if


    ! Take time step limit for the fluids that were calculated so far
    if (present(Cmax_I)) &
         CmaxDt = max(CmaxDt, maxval(CmaxDt_I(iFluidMin:iFluidMax)))

    ! Limit propagation speeds if required
    if (Climit > 0.0) then
       if(present(Cmax_I))  Cmax_I(iFluidMin:iFluidMax) &
            = min( Climit, Cmax_I(iFluidMin:iFluidMax))

       if(present(Cleft_I)) Cleft_I(iFluidMin:iFluidMax) &
            = max(-Climit, Cleft_I(iFluidMin:iFluidMax))

       if(present(Cright_I)) Cright_I(iFluidMin:iFluidMax) &
            = min( Climit, Cright_I(iFluidMin:iFluidMax))

       if(present(Cmax_I))then
          ! If Climit has reduced the diffusion, then the block has to be
          ! advanced with the implicit scheme, so set CmaxDt to a huge number
          if(Climit < CmaxDt) CmaxDt = 1e30
       end if

    end if

  contains

    !========================================================================
    subroutine get_boris_speed

      use ModPhysics, ONLY: Gamma, inv_c2LIGHT, GammaElectron
      use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure

      real :: InvRho, Sound2, FullBx, FullBy, FullBz, FullBn, FullB2
      real :: p, Ppar, Pperp, BnInvB2, GammaPe
      real :: Alfven2, Alfven2Normal, Un, Fast2, Discr, Fast, Slow
      real :: GammaA2, GammaU2
      real :: UnBoris, Sound2Boris, Alfven2Boris, Alfven2NormalBoris
      !-----------------------------------------------------------------------
      ! No explicit formula for multi-ion fluids
      if (nTrueIon > 1) call stop_mpi &
           ('get_boris_speed should not be called with multi-ion fluids')

      InvRho = 1.0/State_V(Rho_)
      ! iPIon_I = p_ for single ion MHD case. iPIon_I is need to add the 
      ! electron pressure(s) for single ion five- and six-moment case.
      p = sum(State_V(iPIon_I))
      FullBx = State_V(Bx_) + B0x
      FullBy = State_V(By_) + B0y
      FullBz = State_V(Bz_) + B0z
      FullB2 = FullBx**2 + FullBy**2 + FullBz**2
      FullBn = NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz

      ! Calculate sound speed squared
      if(UseAnisoPressure .and. FullB2 > 0)then
         ! iPparIon_I = Ppar_ for single ion MHD case. iPparIon_I is need to 
         ! add the electron pressure(s) for single ion six-moment case.
         Ppar  = sum(State_V(iPparIon_I))
         Pperp = (3*p - Ppar)/2.
         BnInvB2 = FullBn**2/FullB2
         Sound2 = InvRho*(2*Pperp + (2*Ppar - Pperp)*BnInvB2)
      else
         Sound2 = InvRho*Gamma*p
      end if

      ! Add contribution of electron pressure
      if(UseElectronPressure)then
         GammaPe = GammaElectron*State_V(Pe_)
         Sound2  = Sound2 + InvRho*GammaPe
      else
         ! For five- and six-moment, Pe should be 0 because electron pressure
         ! has already been added.
         GammaPe = 0.0
      end if

      ! Wave pressure = (GammaWave - 1)*WaveEnergy
      if(UseWavePressure) Sound2 = Sound2 + InvRho*GammaWave &
           * (GammaWave - 1)*sum(State_V(WaveFirst_:WaveLast_))

      Alfven2= FullB2*InvRho
      Alfven2Normal = InvRho*FullBn**2

      Un = State_V(Ux_)*NormalX + State_V(Uy_)*NormalY + State_V(Uz_)*NormalZ

      ! "Alfven Lorentz" factor
      GammaA2 = 1.0/(1.0 + Alfven2*inv_c2LIGHT) 

      ! 1-gA^2*Un^2/c^2
      GammaU2 = max(0.0, 1.0 - GammaA2*Un**2*inv_c2LIGHT) 

      ! Modified speeds
      Sound2Boris        = Sound2        *GammaA2*(1+Alfven2Normal*inv_c2LIGHT)
      Alfven2Boris       = Alfven2       *GammaA2*GammaU2
      Alfven2NormalBoris = Alfven2Normal *GammaA2*GammaU2

      ! Approximate slow and fast wave speeds
      Fast2  = Sound2Boris + Alfven2Boris

      if(UseAnisoPressure .and. FullB2 > 0)then 
         Discr = sqrt(max(0.0, Fast2**2  &
              + 4*((Pperp*InvRho)**2*BnInvB2*(1 - BnInvB2) &  
              - 3*Ppar*Pperp*InvRho**2*BnInvB2*(2 - BnInvB2) &
              + 3*Ppar*Ppar*(InvRho*BnInvB2)**2 &
              - (3*Ppar + GammaPe)*InvRho*Alfven2NormalBoris &
              + GammaPe*InvRho**2*(4*Ppar*BnInvB2 &
              - 3*Ppar - Pperp*BnInvB2)*BnInvB2)))

      else
         Discr = sqrt(max(0.0, Fast2**2 - 4.0*Sound2*Alfven2NormalBoris))
      end if

      ! Get fast and slow speeds multiplied with the face area
      Fast = sqrt( 0.5*(          Fast2 + Discr) )
      Slow = sqrt( 0.5*( max(0.0, Fast2 - Discr) ) )

      ! In extreme cases "slow" wave can be faster than "fast" wave
      ! so take the maximum of the two

      if(UseAwSpeed)then
         Un           = min(UnRight, UnLeft)
         Cleft_I(1)   = min(Un*GammaA2 - Fast, Un - Slow)
         Un           = max(UnLeft, UnRight)
         Cright_I(1)  = max(Un*GammaA2 + Fast, Un + Slow)
         Cmax_I(1)    = max(Cright_I(1), -Cleft_I(1))
         CmaxDt_I(1)  = Cmax_I(1)
      else
         UnBoris            = Un*GammaA2
         if(present(Cmax_I))then
            Cmax_I(1)   = max(abs(UnBoris) + Fast, abs(Un) + Slow)
            CmaxDt_I(1) = Cmax_I(1)
         end if
         if(present(Cleft_I))  Cleft_I(1)  = min(UnBoris - Fast, Un - Slow)
         if(present(Cright_I)) Cright_I(1) = max(UnBoris + Fast, Un + Slow)
      end if

    end subroutine get_boris_speed

    !========================================================================

    subroutine get_mhd_speed

      use ModB0,       ONLY: UseCurlB0
      use ModPhysics,  ONLY: Inv_C2Light, ElectronPressureRatio, &
           GammaElectron, GammaMinus1, Gamma_I
      use ModNumConst, ONLY: cPi
      use ModAdvance,  ONLY: State_VGB, eFluid_, UseElectronPressure, &
           UseAnisoPressure

      real:: UnMin, UnMax
      real:: Rho, InvRho, GammaPe, Pw, Sound2, Ppar, Ppar1, Pperp
      real:: FullBx, FullBy, FullBz, FullBn, FullB2, BnInvB2
      real:: Alfven2, Alfven2Normal, Un, Fast2, Discr, Fast, FastDt, cWhistler
      real:: dB1dB1

      real :: FullBt, Rho1, cDrift, cHall, HallUnLeft, HallUnRight, &
           B1B0L, B1B0R

      real :: MultiIonFactor, ChargeDens_I(nIonFluid)

      integer:: jFluid, iVar

      character(len=*), parameter:: NameSub=NameMod//'::get_mhd_speed'
      !------------------------------------------------------------------------

      if(DoTestCell) then
         write(*,*) NameSub,' nTrueIon     =', nTrueIon
         write(*,*) NameSub,' State_V, B0  ='
         do ivar=1,nVar
            write(*,'(a8,a1,es13.5)') NameVar_V(iVar),'=',State_V(iVar)
         end do
         write(*,'(a9,es13.5)') 'B0x  =', B0x
         write(*,'(a9,es13.5)') 'B0y  =', B0y
         write(*,'(a9,es13.5)') 'B0z  =', B0z
      end if

      Rho = State_V(iRhoIon_I(1))
      Sound2 = State_V(iPIon_I(1))*Gamma_I(1)/Rho
      Un = sum( State_V(iUxIon_I(1):iUzIon_I(1))*Normal_D )
      UnMin = Un
      UnMax = Un

      if(DoTestCell) write(*,*) & 
           NameSub,' Initial Sound2              =', Sound2

      do jFluid = 2, nTrueIon
         Rho1= State_V(iRhoIon_I(jFluid))
         Rho = Rho + Rho1
         ! The (approximate) fast speed fromula for multi-ion MHD
         ! contains the maximum of ion sound speeds squared
         Sound2 = max(Sound2, State_V(iPIon_I(jFluid))*Gamma_I(jFluid)/Rho1)
         Un = sum( State_V(iUxIon_I(jFluid):iUzIon_I(jFluid))*Normal_D )
         ! A reliable upper and lower estimate for wave speeds 
         ! uses the max and min of all ion bulk velocities.
         UnMin = min(Un, UnMin)
         UnMax = max(Un, UnMax)
      end do

      if(DoTestCell) then
         write(*,*) NameSub,' UnMin, UnMax                =', UnMin, UnMax
         write(*,*) NameSub,' After multi-fluid, Sound2   =', Sound2
      end if

      ! InvRho = 1/Sum(RhoIon_I)
      InvRho = 1.0/Rho

      if(UseMultiIon)then
         ! The Alfven velocity and the electron pressure are multiplied
         ! with a Factor >= 1 in multi-ion MHD.
         ! UseMultiIon = .false. for the five moment case.
         ChargeDens_I = ChargePerMass_I*State_V(iRhoIon_I)
         MultiIonFactor = &
              Rho*sum(ChargeDens_I**2/State_V(iRhoIon_I))/sum(ChargeDens_I)**2

         ! Add contribution of electron pressure=fraction of ion pressure
         if(.not.UseElectronPressure) Sound2 = Sound2 + MultiIonFactor &
              *GammaElectron*sum(State_V(iPIon_I))*ElectronPressureRatio*InvRho

         if(DoTestCell) then
            write(*,*) NameSub,' UseMultiIon, ChargeDens_I   =', ChargeDens_I
            write(*,*) NameSub,' UseMultiIon, MultiIonFactor =', MultiIonFactor
         end if
      end if

      if(UseElectronPressure)then
         ! UseElectronPressure = .false. for the five moment case
         GammaPe = GammaElectron*State_V(Pe_)
         if(UseMultiIon) GammaPe = GammaPe*MultiIonFactor
         Sound2 = Sound2 + GammaPe*InvRho
         if(DoTestCell) then
            write(*,*) NameSub,' GammaPe                     =', GammaPe
            write(*,*) NameSub,' after Pe correction, Sound2 =', Sound2
         end if
      else
         GammaPe = 0.0 ! possibly needed for aniso pressure
      endif

      if(UseEfield) then
         ChargeDens_I = ChargePerMass_I*State_V(iRhoIon_I)
         ! MultiIonFactor is used to correct Alfven2/Alfven2Normal, 
         ! it must be calculated first, even for single ion fluid.
         if(nTrueIon > 1)then
            MultiIonFactor = Rho*sum(              &
                 ChargeDens_I(1:nTrueIon)**2       &
                 /State_V(iRhoIon_I(1:nTrueIon)) ) &
                 /sum(ChargeDens_I(1:nTrueIon))**2
         else
            MultiIonFactor = 1.0
         end if

         ! Added electron pressure to Sound2 for the five moment equation.
         ! GammaPe = 0.0 for the six moment equation.
         if (.not. UseAnisoPressure) then
            GammaPe = sum(Gamma_I(Electron_:nIonFluid) &
                 *State_V(iPIon_I(Electron_:nIonFluid)) )
            GammaPe = GammaPe*MultiIonFactor
            Sound2  = Sound2 + GammaPe*InvRho
         end if

         if(DoTestCell) then
            write(*,*) NameSub,' UseEfield, MultiIonFactor   =', MultiIonFactor
            write(*,*) NameSub,' UseEfield, GammaPe          =', GammaPe
            write(*,*) NameSub,' after Pe correction, Sound2 =', Sound2
         end if
      end if

      if(UseRS7) Sound2 = Sound2 + GammaMinus1*DiffBb*InvRho

      if(UseWavePressure)then
         if(UseWavePressureLtd)then
            Sound2 = Sound2 + &
                 GammaWave * (GammaWave - 1)*&
                 max(StateLeft_V(Ew_)/StateLeft_V(Rho_),&
                 StateRight_V(Ew_)/StateRight_V(Rho_))
         else
            Pw = (GammaWave - 1)*sum(State_V(WaveFirst_:WaveLast_))
            Sound2 = Sound2 + GammaWave*Pw*InvRho
         end if
      end if

      FullBx = State_V(Bx_) + B0x
      FullBy = State_V(By_) + B0y
      FullBz = State_V(Bz_) + B0z
      if(UseAwSpeed)then
         ! According to I. Sokolov adding (Bright-Bleft)^2/4 to
         ! the average field squared (Bright+Bleft)^2/4 results in
         ! an upper estimate of the left and right Alfven speeds 
         ! max(Bleft^2/RhoLeft, Bright^2/RhoRight)/
         !
         ! For B0=Bleft=0 and Bright=1 RhoLeft=RhoRight=1 
         ! this is clearly not true.
         !
         dB1dB1 = 0.25*sum((StateRight_V(Bx_:Bz_)-StateLeft_V(Bx_:Bz_))**2)
         Alfven2= (FullBx**2 + FullBy**2 + FullBz**2 + dB1dB1)*InvRho
      else
         Alfven2= (FullBx**2 + FullBy**2 + FullBz**2)*InvRho
      end if
      if(UseCurlB0)then
         B1B0L = StateLeft_V(Bx_)*B0x &
              +  StateLeft_V(By_)*B0y &
              +  StateLeft_V(Bz_)*B0z
         B1B0R = StateRight_V(Bx_)*B0x &
              +  StateRight_V(By_)*B0y &
              +  StateRight_V(Bz_)*B0z
         Alfven2 = Alfven2 +(abs(B1B0L) - B1B0L + abs(B1B0R) - B1B0R)*InvRho
      end if

      FullBn = NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz
      Alfven2Normal = InvRho*FullBn**2

      if(UseMultiIon .or. UseEfield)then
         Alfven2 = Alfven2*MultiIonFactor
         Alfven2Normal = Alfven2Normal*MultiIonFactor
      end if

      if(DoTestCell) &
           write(*,*) NameSub,' Alfven2, Alfven2Normal      =', &
           Alfven2, Alfven2Normal

      ! Calculate fast speed for anisotropic ion pressure.
      ! Formulas refer to V. B. Baranov, 1970 and MAPLE calculation
      if(UseAnisoPressure) FullB2 = FullBx**2 + FullBy**2 + FullBz**2
      if(UseAnisoPressure .and. FullB2 > 0)then
         Ppar  = State_V(Ppar_)
         Pperp = (3*State_V(p_) - Ppar)/2.
         if(.not. IsMhd)then
!!! Most likely the parallel and perpendicular sound speeds should be added up here !!!
            do jFluid = IonFirst_+1, IonLast_
               Ppar1 = State_V(iPparIon_I(jFluid))
               Ppar  = Ppar + Ppar1
               Pperp = Pperp + 0.5*(3*State_V(iP_I(jFluid)) - Ppar1)
            end do
         end if
         BnInvB2 = FullBn**2/FullB2
         Sound2 = InvRho*(2*Pperp + (2*Ppar - Pperp)*BnInvB2 + GammaPe)
         Fast2 = Sound2 + Alfven2 
         Discr = sqrt(max(0.0, Fast2**2  &
              + 4*((Pperp*InvRho)**2*BnInvB2*(1 - BnInvB2) &  
              - 3*Ppar*Pperp*InvRho**2*BnInvB2*(2 - BnInvB2) &
              + 3*Ppar*Ppar*(InvRho*BnInvB2)**2 &
              - (3*Ppar + GammaPe)*InvRho*Alfven2Normal &
              + GammaPe*InvRho**2*(4*Ppar*BnInvB2 &
              - 3*Ppar - Pperp*BnInvB2)*BnInvB2)))
         if(DoTestCell) write(*,*) NameSub,' AnisoP, Sound2, Fast2, Discr=', &
              Sound2, Fast2, Discr
      else
         Fast2  = Sound2 + Alfven2
         Discr  = sqrt(max(0.0, Fast2**2 - 4*Sound2*Alfven2Normal))
         if(DoTestCell) write(*,*) NameSub,' Sound2, Fast2, Discr        =', &
              Sound2, Fast2, Discr
      endif

      if(Fast2 + Discr < 0.0)then
         write(*,*)NameSub, &
              ' negative fast speed squared, Fast2, Discr=', Fast2, Discr
         write(*,*) NameSub, &
              ' iFluid, rho, p(face)   =', iFluid, Rho, State_V(p_)
         if(UseAnisoPressure) write(*,*) NameSub, &
              ' Ppar, Perp             =', Ppar, Pperp
         if(UseElectronPressure) write(*,*) NameSub, &
              ' State_V(Pe_)           =', State_V(Pe_)
         if(UseWavePressure) write(*,*) NameSub, &
              ' GammaWave, State(Waves)=', &
              GammaWave, State_V(WaveFirst_:WaveLast_)

         write(*,*) NameSub, &
              ' Sound2, Alfven2       =', Sound2, Alfven2
         write(*,*) NameSub, &
              ' FullBn, Alfven2Normal =', FullBn, Alfven2Normal
         write(*,*) NameSub, &
              ' FullBx, FullBy, FullBz=', FullBx, FullBy, FullBz
         write(*,*) NameSub, &
              ' State_VGB(left)       =', &
              State_VGB(:,iLeft,jLeft,kLeft,iBlockFace)
         write(*,*) NameSub, &
              ' State_VGB(right)      =', &
              State_VGB(:,iRight,jRight,kRight,iBlockFace)
         write(*,*)NameSub, &
              ' Xyz_DGB(right)        =', &
              Xyz_DGB(:,iFace,jFace,kFace,iBlockFace)

         write(*,*) NameSub, &
              ' iDim,i,j,k,BlockFace,iProc=', &
              iDimFace, iFace, jFace, kFace, iBlockFace, iProc

         call stop_mpi(NameSub//' negative fast speed squared')
      end if

      ! Fast speed multipled by the face area
      if(UseBorisSimple .or. (UseEfield))then
         Fast = sqrt( 0.5*(Fast2 + Discr) &
              /       (1.0 + Alfven2*Inv_C2light) )
      else
         Fast = sqrt( 0.5*(Fast2 + Discr) )
      end if
      FastDt = Fast

      if(DoTestCell) &
           write(*,*) NameSub,' multipled face area, Fast   =', Fast

      ! Add whistler wave speed for the shortest wavelength 2 dx
      if(HallCoeff > 0.0 .and. DoHallInduction) then
         ! Tangential component of B
         FullBt = sqrt(max(0.0, &
              (FullBx**2 + FullBy**2 + FullBz**2) - FullBn**2))

         ! Calculate Ln = d ln(Rho)/dx = (dRho/dx) / Rho
         Rho1 = sum(State_VGB(iRhoIon_I,iLeft,jLeft,kLeft,iBlockFace))

         ! Calculate drift speed and whistler speed
         cDrift    = abs(FullBt)*2.0*abs(Rho1 - Rho)/(Rho1 + Rho)
         cWhistler = cPi*abs(FullBn)

         ! Take the faster speed
         cHall     = HallCoeff*InvDxyz*InvRho*max(cWhistler, cDrift)

         !cHall    = HallCoeff*InvDxyz*InvRho*cWhistler
         FastDt = Fast + cHall
         Fast   = Fast + HallCmaxFactor*cHall
         if(DoTestCell) then
            write(*,*) NameSub,' HallCoeff,HallCmaxFactor    =', &
                 HallCoeff,HallCmaxFactor
            write(*,*) NameSub,' after whistler wave, Fast   =', Fast
         end if
      end if

      HallUnLeft  = UnLeft_I(eFluid_)
      HallUnRight = UnRight_I(eFluid_)

      if(UseAlfvenWaves.and.UseAwSpeed) then
         !\
         ! In this case the propagation speed for 
         ! Alfven waves equal to the Alvfen speed
         ! may happen to be larger that the fast wave
         ! speed in the "hat" state
         !/
         FullBx = StateLeft_V(Bx_) + B0x
         FullBy = StateLeft_V(By_) + B0y
         FullBz = StateLeft_V(Bz_) + B0z
         FullBn = NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz
         Fast = max(Fast, sqrt( FullBn*FullBn / StateLeft_V(iRhoIon_I(1)) ))

         FullBx = StateRight_V(Bx_) + B0x
         FullBy = StateRight_V(By_) + B0y
         FullBz = StateRight_V(Bz_) + B0z
         FullBn = NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz
         Fast = max(Fast, sqrt( FullBn*FullBn / StateRight_V(iRhoIon_I(1)) ))

         if(DoTestCell) &
              write(*,*) NameSub,' after Alfven waves, Fast    =', Fast
      end if

      if(UseAwSpeed)then
         if(HallCoeff > 0.0)then
            Cleft_I(1)   = min(UnLeft, UnRight, HallUnLeft, HallUnRight)
            Cright_I(1)  = max(UnLeft, UnRight, HallUnLeft, HallUnRight)
            CmaxDt_I(1)  = max(Cright_I(1) + FastDt, - Cleft_I(1) - FastDt)
            Cleft_I(1)   = Cleft_I(1)  - Fast
            Cright_I(1)  = Cright_I(1) + Fast
            Cmax_I(1)    = max(Cright_I(1), -Cleft_I(1))
         else
            Cleft_I(1)   = min(UnLeft, UnRight) - Fast
            Cright_I(1)  = max(UnLeft, UnRight) + Fast
            Cmax_I(1)    = max(Cright_I(1), -Cleft_I(1))
            CmaxDt_I(1) = Cmax_I(1)
         end if
      else
         if(present(Cmax_I))then
            if(HallCoeff > 0.0)then
               Cmax_I(1)   = max(abs(UnMin), abs(UnMax), &
                    abs(HallUnLeft), abs(HallUnRight))
               CmaxDt_I(1) = Cmax_I(1) + FastDt
               Cmax_I(1)   = Cmax_I(1) + Fast
            else
               Cmax_I(1)   = max(abs(UnMin), abs(UnMax)) + Fast
               CmaxDt_I(1) = Cmax_I(1)
            end if
         end if
         if(present(Cleft_I))  Cleft_I(1)  = UnMin - Fast
         if(present(Cright_I)) Cright_I(1) = UnMax + Fast
      end if

      if(DoTestCell)then
         if(.not.IsCartesian)then
            write(*,*)NameSub,' AreaX,Y,Z      =', AreaX, AreaY, AreaZ
            write(*,*)NameSub,' Area,Area2     =', Area, Area2
         end if
         if(UseAwSpeed)then
            write(*,*)NameSub,' UnLeft         =', UnLeft
            write(*,*)NameSub,' UnRight        =', UnRight
         end if
         if(HallCoeff > 0.0) then
            write(*,*)NameSub,' HallCoeff      =', HallCoeff
            write(*,*)NameSub,' HallUnLeft     =', HallUnLeft
            write(*,*)NameSub,' HallUnRight    =', HallUnRight
         end if
         write(*,*)NameSub,' Un             =', Un
         write(*,*)NameSub,' Csound2        =', Sound2
         write(*,*)NameSub,' Cfast2         =', Fast2
         write(*,*)NameSub,' Discr2         =', Discr**2
         write(*,*)NameSub,' Calfven        =', sqrt(Alfven2)
         write(*,*)NameSub,' Calfven_normal =', sqrt(Alfven2Normal)
         write(*,*)NameSub,' Cfast          =', Fast
         if(present(Cmax_I)) write(*,*)NameSub,' Cmax_I(1)=',Cmax_I(1)
      end if

    end subroutine get_mhd_speed
    !========================================================================
    subroutine get_hd_speed

      use ModAdvance, ONLY: UseElectronPressure, State_VGB
      use ModPhysics, ONLY: Gamma_I, GammaElectron

      real :: InvRho, Sound2, Sound, Un, GammaP

      character(len=*), parameter:: NameSub=NameMod//'::get_hd_speed'
      !------------------------------------------------------------------------

      if(DoTestCell) then 
         write(*,'(1x,a,a,i3,i3)')    NameSub,' iRho, iP =',iRho, iP
         write(*,'(1x,a,a,30es13.5)') NameSub,' State_V  =',State_V(iRho:iP)
      end if

      ! Calculate sound speed and normal speed
      InvRho = 1.0/State_V(iRho)
      GammaP = Gamma_I(iFluid)*State_V(iP)

      ! If no ion fluids present, then electron pressure is added to
      ! the neutral pressure (CRASH applications)
      if(UseElectronPressure .and. iFluid==1 .and. .not.UseB) &
           GammaP = GammaP + GammaElectron*State_V(Pe_)

      ! Similarly for wave pressure
      ! Wave pressure = (GammaWave - 1)*WaveEnergy
      if(UseWavePressure .and. iFluid==1 .and. .not.UseB) GammaP = GammaP &
           + GammaWave*(GammaWave-1)*sum(State_V(WaveFirst_:WaveLast_))

      ! Sound speed
      Sound2 = GammaP*InvRho

      if(Sound2 <= 0.0)then
         write(*,*)NameSub,' negative Sound2=',Sound2
         write(*,*)NameSub,' iFluid, rho, gamma, p(face) =', &
              iFluid, State_V(iRho), Gamma_I(iFluid), State_V(iP)

         if(UseWavePressure)write(*,*)NameSub,' GammaWave, State(Waves):',&
              GammaWave, State_V(WaveFirst_:WaveLast_)

         if(UseElectronPressure) &
              write(*,*)NameSub,' GammaElectron, State_V(Pe_)=',&
              GammaElectron,State_V(Pe_)

         write(*,*)NameSub,' State(left cell) =', &
              State_VGB(:,iLeft,jLeft,kLeft,iBlockFace)
         write(*,*)NameSub,' State(right cell)=', &
              State_VGB(:,iRight,jRight,kRight,iBlockFace)
         write(*,*)NameSub,' idim,i,j,k,BlockFace,iProc=', &
              iDimFace, iFace, jFace, kFace, iBlockFace, iProc
         write(*,*)NameSub,' xyz(right)=', &
              Xyz_DGB(:,iFace,jFace,kFace,iBlockFace)
         call stop_mpi(NameSub//' negative soundspeed squared')
      end if

      Sound = sqrt(Sound2)
      Un    = sum(State_V(iUx:iUz)*Normal_D)

      if(DoAw)then
         Cleft_I(iFluid)  = min(UnLeft, UnRight) - Sound
         Cright_I(iFluid) = max(UnLeft, UnRight) + Sound
         Cmax_I(iFluid)   = max(Cright_I(iFluid), -Cleft_I(iFluid))
         CmaxDt_I(iFluid) = Cmax_I(iFluid)
      else
         if(present(Cmax_I))then
            Cmax_I(iFluid)   = abs(Un) + Sound
            CmaxDt_I(iFluid) = Cmax_I(iFluid)
         end if
         if(present(Cleft_I))  Cleft_I(iFluid)  = Un - Sound
         if(present(Cright_I)) Cright_I(iFluid) = Un + Sound
      end if

      if(DoTestCell)then
         write(*,*)NameSub,' Un     =',Un
         write(*,*)NameSub,' Csound =',Sound
         if(present(Cmax_I))write(*,*)NameSub,' Cmax   =',Cmax_I(iFluid)
      end if

    end subroutine get_hd_speed
    !========================================================================
    subroutine get_burgers_speed

      real:: Rho
      !---------------------------------------------------------------------
      Rho = State_V(Rho_)
      if(present(Cmax_I))then
         Cmax_I(1)   = Rho
         CmaxDt_I(1) = Rho
      end if
      if(present(Cleft_I))  Cleft_I(1)  = Rho
      if(present(Cright_I)) Cright_I(1) = Rho

    end subroutine get_burgers_speed

  end subroutine get_speed_max

  !============================================================================
  subroutine correct_u_normal(iDim)

    ! Make Unormal 6th order accuracte                 
    use ModMultiFluid, ONLY: iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
    use ModAdvance,    ONLY: State_VGB
    use BATL_lib,  ONLY: correct_face_value, CellCoef_DDGB, &
         Xi_, Eta_, Zeta_, nDim

    integer, intent(in):: iDim

    integer:: iFluid, iRho, iRhoUx, iRhoUy, iRhoUz
    real :: Ucell_I(4), Unormal, Ucell_D(3)
    real :: Normal0_D(3), Area0_D(3)
    integer:: iCell, iCount
    !--------------------------------------------------------------------------

    do iFluid = iFluidMin, iFluidMax

       Unormal = Unormal_I(iFluid)
       iRho = iRho_I(iFluid)
       iRhoUx = iRhoUx_I(iFluid)
       iRhoUz = iRhoUz_I(iFluid)

       if(.not. IsCartesian) then
          if(iDim == x_) then
             iCount = 1
             do iCell = iFace - 2, iFace + 1
                Area0_D(1:nDim) = CellCoef_DDGB( &
                     Xi_,:,iCell,jFace,kFace,iBlockFace)
                Normal0_D(1:nDim) = Area0_D(1:nDim)/sqrt(sum(Area0_D(1:nDim)**2))
                Ucell_D = &
                     State_VGB(iRhoUx:iRhoUz,iCell,jFace,kFace,iBlockFace)/&
                     State_VGB(iRho,iCell,jFace,kFace,iBlockFace)
                Ucell_I(iCount) = sum(Ucell_D(1:nDim)*Normal0_D(1:nDim))
                iCount = iCount + 1
             enddo

          elseif(iDim == y_) then
             iCount = 1
             do iCell = jFace - 2, jFace + 1
                Area0_D(1:nDim) = CellCoef_DDGB( &
                     Eta_,:,iFace,iCell,kFace,iBlockFace)
                Normal0_D(1:nDim) = Area0_D(1:nDim)/sqrt(sum(Area0_D(1:nDim)**2))

                Ucell_D = &
                     State_VGB(iRhoUx:iRhoUz,iFace,iCell,kFace,iBlockFace)/&
                     State_VGB(iRho,iFace,iCell,kFace,iBlockFace)
                Ucell_I(iCount) = sum(Ucell_D(1:nDim)*Normal0_D(1:nDim))
                iCount = iCount + 1
             enddo
          elseif(iDim == z_) then
             iCount = 1
             do iCell = kFace - 2, kFace + 1
                Area0_D = CellCoef_DDGB( &
                     Zeta_,:,iFace,jFace,iCell,iBlockFace)
                Normal0_D = Area0_D/sqrt(sum(Area0_D(:)**2))

                Ucell_D = &
                     State_VGB(iRhoUx:iRhoUz,iFace,jFace,iCell,iBlockFace)/&
                     State_VGB(iRho,iFace,jFace,iCell,iBlockFace)
                Ucell_I(iCount) = sum(Ucell_D*Normal0_D)
                iCount = iCount + 1
             enddo
          endif
       else
          if(iDim == x_) then
             iRhoUx = iRhoUx_I(iFluid)
             Ucell_I = &
                  State_VGB(iRhoUx,iFace-2:iFace+1,jFace,kFace,iBlockFace)/&
                  State_VGB(iRho,iFace-2:iFace+1,jFace,kFace,iBlockFace)
          elseif(iDim == y_) then
             iRhoUy = iRhoUy_I(iFluid)
             Ucell_I = &
                  State_VGB(iRhoUy,iFace,jFace-2:jFace+1,kFace,iBlockFace)/&
                  State_VGB(iRho,iFace,jFace-2:jFace+1,kFace,iBlockFace)
          else
             iRhoUz = iRhoUz_I(iFluid)
             Ucell_I = State_VGB(iRhoUz,iFace,jFace,kFace-2:kFace+1,iBlockFace)/&
                  State_VGB(iRho,iFace,jFace,kFace-2:kFace+1,iBlockFace)
          endif
       endif
       Unormal_I(iFluid) = correct_face_value(Unormal, Ucell_I)
    enddo

  end subroutine correct_u_normal

end module ModFaceFlux

!==============================================================================
subroutine roe_solver(Flux_V)

  use ModFaceFlux, ONLY: &
       nFlux, IsBoundary, Climit, &
       StateLeft_V,  StateRight_V, FluxLeft_V, FluxRight_V, &
       StateLeftCons_V, StateRightCons_V, CmaxDt, Unormal_I, &
       nFluxMhd, RhoMhd_, RhoUn_, RhoUt1_, RhoUt2_, &
       B1n_, B1t1_, B1t2_, eMhd_, pMhd_, B0n, B0t1, B0t2, &
       UnL, Ut1L, Ut2L, B1nL, B1t1L, B1t2L, &
       UnR, Ut1R, Ut2R, B1nR, B1t1R, B1t2R, &
       rotate_state_vectors, rotate_flux_vector

  use ModVarIndexes, ONLY: Rho_, p_, Energy_, ScalarFirst_, ScalarLast_

  use ModPhysics,  ONLY: Gamma,GammaMinus1,InvGammaMinus1
  use ModNumConst

  implicit none

  real,    intent(out):: Flux_V(nFlux)

  ! Number of MHD waves including the divB wave
  integer, parameter :: nWaveMhd=8

  ! Named MHD wave indexes
  integer, parameter :: EntropyW_=1, AlfvenRW_=2, AlfvenLW_=3, &
       SlowRW_=4, FastRW_=5, SlowLW_=6, FastLW_=7, DivBW_=8 

  ! Loop variables
  integer :: iFlux, iVar, iWave

  ! Left and right face
  real :: RhoL, BnL, Bt1L, Bt2L, BbL, pL, eL, aL, CsL, CfL
  real :: RhoR, BnR, Bt1R, Bt2R, BbR, pR, eR, aR, CsR, CfR

  ! Roe average (hat)
  real :: RhoH,UnH,Ut1H,Ut2H
  real :: BnH,Bt1H,Bt2H,BbH
  real :: B1nH,B1t1H,B1t2H,Bb1H
  real :: pH,eH,UuH
  real :: aH,CsH,CfH

  real :: BetaY, BetaZ, AlphaS, AlphaF

  real :: RhoInvL,RhoInvR,RhoInvH
  real :: RhoSqrtH,    RhoSqrtL,    RhoSqrtR, &
       RhoInvSqrtH, RhoInvSqrtL, RhoInvSqrtR

  ! Jump in the conservative state
  real, dimension(nWaveMhd) :: dCons_V

  ! Eigenvalues and jumps in characteristic variable
  real, dimension(nWaveMhd) :: Eigenvalue_V, DeltaWave_V 

  ! Eigenvectors
  real, dimension(nWaveMhd, nWaveMhd):: EigenvectorL_VV  ! Left  eigenvectors
  real, dimension(nWaveMhd, nFluxMhd):: EigenvectorR_VV  ! Right eigenvectors

  ! Fluxes
  real, dimension(nFluxMhd)       :: Diffusion_V      ! Diffusive fluxes

  ! Misc. scalar variables
  real :: SignBnH, Tmp1, Tmp2, Tmp3, Gamma1A2Inv
  !---------------------------------------------------------------------------
  ! Scalar variables
  RhoL  =  StateLeft_V(Rho_)
  pL    =  StateLeft_V(p_ )
  RhoR  =  StateRight_V(Rho_)
  pR    =  StateRight_V(p_  )

  ! Rotate vector variables into a coordinate system orthogonal to the face
  call rotate_state_vectors

  ! Jump in scalar conservative variables
  dCons_V(RhoMhd_) = RhoR      - RhoL
  dCons_V(RhoUn_ ) = RhoR*UnR  - RhoL*UnL
  dCons_V(RhoUt1_) = RhoR*Ut1R - RhoL*Ut1L
  dCons_V(RhoUt2_) = RhoR*Ut2R - RhoL*Ut2L
  dCons_V(B1n_   ) = B1nR      - B1nL
  dCons_V(B1t1_  ) = B1t1R     - B1t1L
  dCons_V(B1t2_  ) = B1t2R     - B1t2L
  dCons_V(eMhd_)   = StateRightCons_V(Energy_) - StateLeftCons_V(Energy_)

  ! Derived variables
  RhoInvL = 1./RhoL
  BnL  = B0n+B1nL
  Bt1L = B0t1+B1t1L
  Bt2L = B0t2+B1t2L
  BbL  = BnL**2 + Bt1L**2 + Bt2L**2
  aL   = Gamma*pL*RhoInvL

  RhoInvR = 1./RhoR
  BnR  = B0n+B1nR
  Bt1R = B0t1+B1t1R
  Bt2R = B0t2+B1t2R
  BbR  = BnR**2 + Bt1R**2 + Bt2R**2
  aR   = Gamma*pR*RhoInvR

  !\
  ! Hat face
  !/
  RhoH = 0.5*(RhoL + RhoR)
  RhoInvH = 1./RhoH
  UnH  = 0.5*(  UnL +   UnR)
  Ut1H = 0.5*( Ut1L +  Ut1R)
  Ut2H = 0.5*( Ut2L +  Ut2R)
  BnH  = 0.5*(  BnL +   BnR)
  Bt1H = 0.5*( Bt1L +  Bt1R)
  Bt2H = 0.5*( Bt2L +  Bt2R)
  B1nH = 0.5*( B1nL +  B1nR)
  B1t1H= 0.5*(B1t1L + B1t1R)
  B1t2H= 0.5*(B1t2L + B1t2R)
  pH   = 0.5*(   pL +    pR)
  BbH  = BnH**2  + Bt1H**2  + Bt2H**2

  Bb1H = B1nH**2 + B1t1H**2 + B1t2H**2
  aH   = Gamma*pH*RhoInvH

  !if(aL<0.0)then
  !   write(*,*)'NEGATIVE aL Me, iDir, i, j, k, iBlockFace',&
  !        aL,iProc,iDimFace,i,j,k,&
  !        Xyz_DGB(:,i,j,k,iBlockFace)
  !   call stop_mpi
  !end if

  aL=sqrt(aL)
  aR=sqrt(aR)
  aH=sqrt(aH)

  eL = aL*aL + BbL*RhoInvL
  CfL = max(0., (eL**2 - 4.*aL**2 * BnL**2 * RhoInvL))
  eR = aR**2 + BbR*RhoInvR
  CfR = max(0., (eR**2 - 4.*aR**2 * BnR**2 * RhoInvR))
  eH = aH**2 + BbH*RhoInvH
  CfH = max(0., (eH**2 - 4.*aH**2 * BnH**2 * RhoInvH))

  CfL=sqrt(CfL)
  CfR=sqrt(CfR)
  CfH=sqrt(CfH)

  CsL  = max(0.,0.5*(eL-CfL))
  CfL  = 0.5*(eL+CfL)

  CsR  = max(0.,0.5*(eR-CfR))
  CfR  = 0.5*(eR+CfR)

  CsH  = max(0.,0.5*(eH-CfH))
  CfH  = 0.5*(eH+CfH)

  UuH  = UnH**2 + Ut1H**2 + Ut2H**2
  eH   = pH*InvGammaMinus1 + 0.5*RhoH*UuH + 0.5*Bb1H

  CsL=sqrt(CsL)
  CsR=sqrt(CsR)
  CsH=sqrt(CsH)
  CfL=sqrt(CfL)
  CfR=sqrt(CfR)
  CfH=sqrt(CfH)

  CsL  = min(CsL,aL)
  CfL  = max(CfL,aL)
  CsR  = min(CsR,aR)
  CfR  = max(CfR,aR)
  CsH  = min(CsH,aH)
  CfH  = max(CfH,aH)

  !\
  ! Non-dimensional scaling factors
  !/
  Tmp1 = Bt1H**2 + Bt2H**2

  if(Tmp1 > 1e-8)then
     Tmp1=sqrt(1./Tmp1)
     BetaY = Bt1H*Tmp1
     BetaZ = Bt2H*Tmp1
  else
     BetaY = cSqrtHalf
     BetaZ = cSqrtHalf
  end if

  Tmp1 = CfH**2 - CsH**2
  if (Tmp1 > 1.0e-08) then
     AlphaF = max(0.0, (aH**2  - CsH**2)/Tmp1)
     AlphaS = max(0.0, (CfH**2 - aH**2 )/Tmp1)

     AlphaF = sqrt(AlphaF)
     AlphaS = sqrt(AlphaS)
  else if (BnH**2 * RhoInvH <= aH**2 ) then
     AlphaF = 1.0
     AlphaS = 0.0
  else
     AlphaF = 0.0
     AlphaS = 1.0
  endif

  !\
  ! Set some values that are reused over and over
  !/

  RhoSqrtH   =sqrt(RhoH)
  RhoSqrtL   =sqrt(RhoL)
  RhoSqrtR   =sqrt(RhoR)
  RhoInvSqrtH=1./RhoSqrtH
  RhoInvSqrtL=1./RhoSqrtL
  RhoInvSqrtR=1./RhoSqrtR


  SignBnH     = sign(1.,BnH)
  Gamma1A2Inv = GammaMinus1 / aH**2

  !\
  ! Eigenvalues
  !/
  Eigenvalue_V(EntropyW_) = UnH
  Eigenvalue_V(AlfvenRW_) = UnH + BnH*RhoInvSqrtH
  Eigenvalue_V(AlfvenLW_) = UnH - BnH*RhoInvSqrtH
  Eigenvalue_V(SlowRW_)   = UnH + CsH
  Eigenvalue_V(FastRW_)   = UnH + CfH
  Eigenvalue_V(SlowLW_)   = UnH - CsH
  Eigenvalue_V(FastLW_)   = UnH - CfH
  Eigenvalue_V(DivBW_)    = UnH

  !\
  ! Entropy fix for Eigenvalues
  !/
  Tmp1 = UnR - UnL
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(1)) < Tmp1*0.5) then
     Eigenvalue_V(1) = sign(1.,Eigenvalue_V(1))*   &
          ((Eigenvalue_V(1)*Eigenvalue_V(1)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR + BnR*RhoInvSqrtR) - (UnL + BnL*RhoInvSqrtL)
  Tmp1 = max(cTiny,4.*Tmp1)
  if (abs(Eigenvalue_V(2)) < Tmp1*0.5) then
     Eigenvalue_V(2) = sign(1.,Eigenvalue_V(2))*   &
          ((Eigenvalue_V(2)*Eigenvalue_V(2)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR - BnR*RhoInvSqrtR) - (UnL - BnL*RhoInvSqrtL)
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(3)) < Tmp1*0.5) then
     Eigenvalue_V(3) = sign(1.,Eigenvalue_V(3))*   &
          ((Eigenvalue_V(3)*Eigenvalue_V(3)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR + CsR) - (UnL + CsL)
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(4)) < Tmp1*0.5) then
     Eigenvalue_V(4) = sign(1.,Eigenvalue_V(4))*   &
          ((Eigenvalue_V(4)*Eigenvalue_V(4)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR+CfR) - (UnL+CfL)
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(5)) < Tmp1*0.5) then
     Eigenvalue_V(5) = sign(1.,Eigenvalue_V(5))*   &
          ((Eigenvalue_V(5)*Eigenvalue_V(5)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR-CsR) - (UnL-CsL)
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(6)) < Tmp1*0.5) then
     Eigenvalue_V(6) = sign(1.,Eigenvalue_V(6))*   &
          ((Eigenvalue_V(6)*Eigenvalue_V(6)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = (UnR-CfR) - (UnL-CfL)
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(7)) < Tmp1*0.5) then
     Eigenvalue_V(7) = sign(1.,Eigenvalue_V(7))*   &
          ((Eigenvalue_V(7)*Eigenvalue_V(7)/Tmp1) + Tmp1*0.25)
  end if

  Tmp1 = UnR - UnL
  Tmp1 = max(cTiny, 4.*Tmp1)
  if (abs(Eigenvalue_V(8)) < Tmp1*0.5) then
     Eigenvalue_V(8) = sign(1.,Eigenvalue_V(8))* &
          ((Eigenvalue_V(8)*Eigenvalue_V(8)/Tmp1) + Tmp1*0.25)
  end if

  !\
  ! Timur's divergence wave fix
  !/
  !original  version
  !      Eigenvalue_V(8)=abs(Eigenvalue_V(8))+aH
  !
  !Enhanced diffusion, the maximum eigenvalue
  Eigenvalue_V(DivbW_) = max( Eigenvalue_V(FastRW_), -Eigenvalue_V(FastLW_))

  ! The original version was proposed by Timur Linde for heliosphere
  ! simulations. Enhanced version was found to be more robust on 8/20/01
  ! The original version was commented out in versions 6x resulting in 
  ! worse stability for the Roe solver.

  ! At inner BC replace all eigenvalues with the enhanced eigenvalue
  ! of divB, which is the maximum eigenvalue
  if(IsBoundary) Eigenvalue_V(1:nWaveMhd) = Eigenvalue_V(DivBW_)

  !\
  ! Eigenvectors
  !/
  Tmp1=1./(2.*RhoH*aH**2)
  Tmp2=RhoInvH*cSqrtHalf
  Tmp3=RhoInvSqrtH*cSqrtHalf

  ! Left eigenvector for Entropy wave
  EigenvectorL_VV(1,1) = 1.-0.5*Gamma1A2Inv*UuH
  EigenvectorL_VV(2,1) = Gamma1A2Inv*UnH
  EigenvectorL_VV(3,1) = Gamma1A2Inv*Ut1H
  EigenvectorL_VV(4,1) = Gamma1A2Inv*Ut2H
  EigenvectorL_VV(5,1) = Gamma1A2Inv*B1nH
  EigenvectorL_VV(6,1) = Gamma1A2Inv*B1t1H
  EigenvectorL_VV(7,1) = Gamma1A2Inv*B1t2H
  EigenvectorL_VV(8,1) = -Gamma1A2Inv

  ! Left eigenvector for Alfven wave +
  EigenvectorL_VV(1,2) = (Ut1H*BetaZ-Ut2H*BetaY)*Tmp2
  EigenvectorL_VV(2,2) = 0.
  EigenvectorL_VV(3,2) = -(BetaZ*Tmp2)
  EigenvectorL_VV(4,2) = (BetaY*Tmp2)
  EigenvectorL_VV(5,2) = 0.
  EigenvectorL_VV(6,2) = (BetaZ*Tmp3)
  EigenvectorL_VV(7,2) = -(BetaY*Tmp3)
  EigenvectorL_VV(8,2) = 0.

  ! Left eigenvector for Alfven wave -
  EigenvectorL_VV(1,3) = (Ut1H*BetaZ-Ut2H*BetaY)*Tmp2
  EigenvectorL_VV(2,3) = 0.
  EigenvectorL_VV(3,3) = -(BetaZ*Tmp2)
  EigenvectorL_VV(4,3) = (BetaY*Tmp2)
  EigenvectorL_VV(5,3) = 0.
  EigenvectorL_VV(6,3) = -(BetaZ*Tmp3)
  EigenvectorL_VV(7,3) = (BetaY*Tmp3)
  EigenvectorL_VV(8,3) = 0.

  ! Left eigenvector for Slow magnetosonic wave +
  EigenvectorL_VV(1,4) = Tmp1* &
       (AlphaS*(GammaMinus1*UuH/2. - UnH*CsH) - &
       AlphaF*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,4) = Tmp1*(AlphaS*(-UnH*GammaMinus1+CsH))
  EigenvectorL_VV(3,4) = Tmp1*(-GammaMinus1*AlphaS*Ut1H+AlphaF*CfH*BetaY*SignBnH)
  EigenvectorL_VV(4,4) = Tmp1*(-GammaMinus1*AlphaS*Ut2H+AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorL_VV(5,4) = Tmp1*(-GammaMinus1*B1nH*AlphaS)
  EigenvectorL_VV(6,4) = Tmp1*(-AlphaF*BetaY*aH*RhoSqrtH-GammaMinus1*B1t1H*AlphaS)
  EigenvectorL_VV(7,4) = Tmp1*(-AlphaF*BetaZ*aH*RhoSqrtH-GammaMinus1*B1t2H*AlphaS)
  EigenvectorL_VV(8,4) = Tmp1*(GammaMinus1*AlphaS)

  ! Left eigenvector for Fast magnetosonic wave +
  EigenvectorL_VV(1,5) = Tmp1* &
       (AlphaF*(GammaMinus1*UuH/2. - UnH*CfH)+AlphaS*CsH*SignBnH* &
       (Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,5) = Tmp1*(AlphaF*(-UnH*GammaMinus1+CfH))
  EigenvectorL_VV(3,5) = Tmp1*(-GammaMinus1*AlphaF*Ut1H-AlphaS*CsH*BetaY*SignBnH)
  EigenvectorL_VV(4,5) = Tmp1*(-GammaMinus1*AlphaF*Ut2H-AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorL_VV(5,5) = Tmp1*(-GammaMinus1*B1nH*AlphaF)
  EigenvectorL_VV(6,5) = Tmp1*(AlphaS*BetaY*aH*RhoSqrtH-GammaMinus1*B1t1H*AlphaF)
  EigenvectorL_VV(7,5) = Tmp1*(AlphaS*BetaZ*aH*RhoSqrtH-GammaMinus1*B1t2H*AlphaF)
  EigenvectorL_VV(8,5) = Tmp1*(GammaMinus1*AlphaF)

  ! Left eigenvector for Slow magnetosonic wave -
  EigenvectorL_VV(1,6) = Tmp1* &
       (AlphaS*(GammaMinus1*UuH/2. + UnH*CsH) + &
       AlphaF*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,6) = Tmp1*(AlphaS*(-UnH*GammaMinus1-CsH))
  EigenvectorL_VV(3,6) = Tmp1*(-GammaMinus1*AlphaS*Ut1H-AlphaF*CfH*BetaY*SignBnH)
  EigenvectorL_VV(4,6) = Tmp1*(-GammaMinus1*AlphaS*Ut2H-AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorL_VV(5,6) = Tmp1*(-GammaMinus1*B1nH*AlphaS)
  EigenvectorL_VV(6,6) = Tmp1*(-AlphaF*BetaY*aH*RhoSqrtH-GammaMinus1*B1t1H*AlphaS)
  EigenvectorL_VV(7,6) = Tmp1*(-AlphaF*BetaZ*aH*RhoSqrtH-GammaMinus1*B1t2H*AlphaS)
  EigenvectorL_VV(8,6) = Tmp1*(GammaMinus1*AlphaS)

  ! Left eigenvector for Fast magnetosonic wave -
  EigenvectorL_VV(1,7) = Tmp1* &
       (AlphaF*(GammaMinus1*UuH/2. + UnH*CfH) - &
       AlphaS*CsH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,7) = Tmp1*(AlphaF*(-UnH*GammaMinus1-CfH))
  EigenvectorL_VV(3,7) = Tmp1*(-GammaMinus1*AlphaF*Ut1H+AlphaS*CsH*BetaY*SignBnH)
  EigenvectorL_VV(4,7) = Tmp1*(-GammaMinus1*AlphaF*Ut2H+AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorL_VV(5,7) = Tmp1*(-GammaMinus1*B1nH*AlphaF)
  EigenvectorL_VV(6,7) = Tmp1*(AlphaS*BetaY*aH*RhoSqrtH-GammaMinus1*B1t1H*AlphaF)
  EigenvectorL_VV(7,7) = Tmp1*(AlphaS*BetaZ*aH*RhoSqrtH-GammaMinus1*B1t2H*AlphaF)
  EigenvectorL_VV(8,7) = Tmp1*(GammaMinus1*AlphaF)

  ! Left eigenvector for Divergence wave
  EigenvectorL_VV(1,8) = 0.
  EigenvectorL_VV(2,8) = 0.
  EigenvectorL_VV(3,8) = 0.
  EigenvectorL_VV(4,8) = 0.
  EigenvectorL_VV(5,8) = 1.
  EigenvectorL_VV(6,8) = 0.
  EigenvectorL_VV(7,8) = 0.
  EigenvectorL_VV(8,8) = 0.

  !coefficient for pressure component of the Right vector
  Tmp1=Gamma*max(pL,pR) 

  !Pressure component is not linearly independent and obeys the 
  ! equation as follows:
  !EigenvectorR_VV(1:8,9)=(0.5*UuH*EigenvectorR_VV(1:8,1)-&
  !                     UnH*EigenvectorR_VV(1:8,2)-&
  !                     Ut1H*EigenvectorR_VV(1:8,3)-&
  !                     Ut2H*EigenvectorR_VV(1:8,4)-&
  !                     B1nH*EigenvectorR_VV(1:8,5)-&
  !                     B1t1H*EigenvectorR_VV(1:8,6)-&
  !                     B1t2H*EigenvectorR_VV(1:8,7)+
  !                     EigenvectorR_VV(1:8,8))*inv_gm1         

  ! Right eigenvector for Entropy wave
  EigenvectorR_VV(1,1) = 1.
  EigenvectorR_VV(1,2) = UnH
  EigenvectorR_VV(1,3) = Ut1H
  EigenvectorR_VV(1,4) = Ut2H
  EigenvectorR_VV(1,5) = 0.
  EigenvectorR_VV(1,6) = 0.
  EigenvectorR_VV(1,7) = 0.
  EigenvectorR_VV(1,eMhd_) = 0.5*UuH
  EigenvectorR_VV(1,pMhd_)=cZero

  ! Right eigenvector for Alfven wave +
  EigenvectorR_VV(2,1) = 0.
  EigenvectorR_VV(2,2) = 0.
  EigenvectorR_VV(2,3) = -BetaZ*RhoH*cSqrtHalf
  EigenvectorR_VV(2,4) = BetaY*RhoH*cSqrtHalf
  EigenvectorR_VV(2,5) = 0.
  EigenvectorR_VV(2,6) = BetaZ*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(2,7) = -BetaY*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(2,eMhd_) = (BetaY*Ut2H - BetaZ*Ut1H)*RhoH*cSqrtHalf &
       + (B1t1H*BetaZ - B1t2H*BetaY)*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(2,pMhd_)=cZero

  ! Right eigenvector for Alfven wave -
  EigenvectorR_VV(3,1) = 0.
  EigenvectorR_VV(3,2) = 0.
  EigenvectorR_VV(3,3) = -BetaZ*RhoH*cSqrtHalf
  EigenvectorR_VV(3,4) = BetaY*RhoH*cSqrtHalf
  EigenvectorR_VV(3,5) = 0.
  EigenvectorR_VV(3,6) = -BetaZ*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(3,7) = BetaY*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(3,eMhd_) = (BetaY*Ut2H - BetaZ*Ut1H)*RhoH*cSqrtHalf &
       - (B1t1H*BetaZ - B1t2H*BetaY)*RhoSqrtH*cSqrtHalf
  EigenvectorR_VV(3,pMhd_)=cZero

  ! Right eigenvector for Slow magnetosonic wave +
  EigenvectorR_VV(4,1) = RhoH*AlphaS
  EigenvectorR_VV(4,2) = RhoH*AlphaS*(UnH+CsH)
  EigenvectorR_VV(4,3) = RhoH*(AlphaS*Ut1H + AlphaF*CfH*BetaY*SignBnH)
  EigenvectorR_VV(4,4) = RhoH*(AlphaS*Ut2H + AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorR_VV(4,5) = 0.
  EigenvectorR_VV(4,6) = -AlphaF*aH*BetaY*RhoSqrtH
  EigenvectorR_VV(4,7) = -AlphaF*aH*BetaZ*RhoSqrtH
  EigenvectorR_VV(4,eMhd_) = &
       AlphaS*(RhoH*UuH*0.5 + Gamma*pH*InvGammaMinus1+RhoH*UnH*CsH) &
       - AlphaF*(aH*RhoSqrtH*(BetaY*B1t1H + BetaZ*B1t2H) &
       - RhoH*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorR_VV(4,pMhd_)=Tmp1*AlphaS

  ! Right eigenvector for Fast magnetosonic wave +
  EigenvectorR_VV(5,1) = RhoH*AlphaF
  EigenvectorR_VV(5,2) = RhoH*AlphaF* (UnH+CfH)
  EigenvectorR_VV(5,3) = RhoH* (AlphaF*Ut1H - AlphaS*CsH*BetaY*SignBnH)
  EigenvectorR_VV(5,4) = RhoH* (AlphaF*Ut2H - AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorR_VV(5,5) = 0.
  EigenvectorR_VV(5,6) = AlphaS*aH*BetaY*RhoSqrtH
  EigenvectorR_VV(5,7) = AlphaS*aH*BetaZ*RhoSqrtH
  EigenvectorR_VV(5,eMhd_) = &
       AlphaF*(RhoH*UuH*0.5 + Gamma*pH*InvGammaMinus1+RhoH*UnH*CfH) &
       + AlphaS*(aH*RhoSqrtH*(BetaY*B1t1H + BetaZ*B1t2H) &
       - RhoH*CsH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorR_VV(5,pMhd_)=Tmp1*AlphaF

  ! Right eigenvector for Slow magnetosonic wave -
  EigenvectorR_VV(6,1) = RhoH*AlphaS
  EigenvectorR_VV(6,2) = RhoH*AlphaS*(UnH-CsH)
  EigenvectorR_VV(6,3) = RhoH* (AlphaS*Ut1H - AlphaF*CfH*BetaY*SignBnH)
  EigenvectorR_VV(6,4) = RhoH* (AlphaS*Ut2H - AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorR_VV(6,5) = 0.
  EigenvectorR_VV(6,6) = - AlphaF*aH*BetaY*RhoSqrtH
  EigenvectorR_VV(6,7) = - AlphaF*aH*BetaZ*RhoSqrtH
  EigenvectorR_VV(6,eMhd_) = &
       AlphaS*(RhoH*UuH*0.5 + Gamma*pH*InvGammaMinus1-RhoH*UnH*CsH) &
       - AlphaF*(aH*RhoSqrtH*(BetaY*B1t1H + BetaZ*B1t2H) &
       + RhoH*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorR_VV(6,pMhd_)=Tmp1*AlphaS

  ! Right eigenvector for Fast magnetosonic wave -
  EigenvectorR_VV(7,1) = RhoH*AlphaF
  EigenvectorR_VV(7,2) = RhoH*AlphaF* (UnH-CfH)
  EigenvectorR_VV(7,3) = RhoH*(AlphaF*Ut1H + AlphaS*CsH*BetaY*SignBnH)
  EigenvectorR_VV(7,4) = RhoH*(AlphaF*Ut2H + AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorR_VV(7,5) = 0.
  EigenvectorR_VV(7,6) = AlphaS*aH*BetaY*RhoSqrtH
  EigenvectorR_VV(7,7) = AlphaS*aH*BetaZ*RhoSqrtH
  EigenvectorR_VV(7,eMhd_) = &
       AlphaF*(RhoH*UuH*0.5 + Gamma*pH*InvGammaMinus1-RhoH*UnH*CfH) &
       + AlphaS*(aH*RhoSqrtH*(BetaY*B1t1H + BetaZ*B1t2H) &
       + RhoH*CsH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorR_VV(7,pMhd_)=Tmp1*AlphaF

  ! Right eigenvector for Divergence wave
  EigenvectorR_VV(8,1) = 0.
  EigenvectorR_VV(8,2) = 0.
  EigenvectorR_VV(8,3) = 0.
  EigenvectorR_VV(8,4) = 0.
  EigenvectorR_VV(8,5) = 1.
  EigenvectorR_VV(8,6) = 0.
  EigenvectorR_VV(8,7) = 0.
  EigenvectorR_VV(8,eMhd_) = B1nH
  EigenvectorR_VV(8,pMhd_) = cZero

  !\
  ! Alphas (elemental wave strengths)
  !/
  ! matmul is slower than the loop for the NAG F95 compiler
  ! DeltaWave_V = matmul(dCons_V, EigenvectorL_VV)
  do iWave = 1, nWaveMhd
     DeltaWave_V(iWave) = sum(dCons_V*EigenvectorL_VV(:,iWave))
  end do

  ! Take absolute value of eigenvalues
  Eigenvalue_V = abs(Eigenvalue_V)

  ! Limit them if required
  if(Climit > 0.0) Eigenvalue_V = min(Climit, Eigenvalue_V)

  !\
  ! Calculate the Roe Interface fluxes 
  ! F = A * 0.5 * [ F_L+F_R - sum_k(|lambda_k| * alpha_k * r_k) ]
  !/
  ! First get the diffusion: sum_k(|lambda_k| * alpha_k * r_k)
  !  Diffusion_V = matmul(abs(Eigenvalue_V)*DeltaWave_V, EigenvectorR_VV)
  do iFlux = 1, nFluxMhd
     Diffusion_V(iFlux) = &
          sum(Eigenvalue_V*DeltaWave_V*EigenvectorR_VV(:,iFlux))
  end do

  ! Scalar variables
  Flux_V(Rho_   ) = Diffusion_V(RhoMhd_)
  Flux_V(P_     ) = Diffusion_V(pMhd_)
  Flux_V(Energy_) = Diffusion_V(eMhd_)

  ! Rotate fluxes of vector variables back
  call rotate_flux_vector(Diffusion_V, Flux_V)

  ! The diffusive flux for the advected scalar variables is simply
  ! 0.5*|Velocity|*(U_R - U_L)
  do iVar = ScalarFirst_, ScalarLast_
     Flux_V(iVar) = abs(UnH)*(StateRightCons_V(iVar) - StateLeftCons_V(iVar))
  end do

  ! Roe flux = average of left and right flux plus the diffusive flux
  Flux_V  = 0.5*(FluxLeft_V + FluxRight_V - Flux_V)

  ! Normal velocity and maximum wave speed
  Unormal_I = UnH
  CmaxDt    = abs(UnH) + CfH

end subroutine roe_solver

!===========================================================================

subroutine calc_electric_field(iBlock)

  ! Calculate the total electric field which includes numerical resistivity
  ! This estimate averages the numerical fluxes to the cell centers 
  ! for sake of simplicity.

  use ModSize,       ONLY: nI, nJ, nK
  use ModVarIndexes, ONLY: Bx_,By_,Bz_
  use ModAdvance,    ONLY: Flux_VX, Flux_VY, Flux_VZ, ExNum_CB, EyNum_CB, EzNum_CB
  use BATL_lib,      ONLY: CellFace_DB

  implicit none
  integer, intent(in) :: iBlock
  !------------------------------------------------------------------------
  ! E_x=(fy+fy-fz-fz)/4
  ExNum_CB(:,:,:,iBlock) = - 0.25*(                                    &
       ( Flux_VY(Bz_,1:nI,1:nJ  ,1:nK  )                            &
       + Flux_VY(Bz_,1:nI,2:nJ+1,1:nK  )) / CellFace_DB(2,iBlock) - &
       ( Flux_VZ(By_,1:nI,1:nJ  ,1:nK  )                            &
       + Flux_VZ(By_,1:nI,1:nJ  ,2:nK+1)) / CellFace_DB(3,iBlock) )

  ! E_y=(fz+fz-fx-fx)/4
  EyNum_CB(:,:,:,iBlock) = - 0.25*(                                    &
       ( Flux_VZ(Bx_,1:nI  ,1:nJ,1:nK  )                            &
       + Flux_VZ(Bx_,1:nI  ,1:nJ,2:nK+1)) / CellFace_DB(3,iBlock) - &
       ( Flux_VX(Bz_,1:nI  ,1:nJ,1:nK  )                            &
       + Flux_VX(Bz_,2:nI+1,1:nJ,1:nK  )) / CellFace_DB(1,iBlock) )

  ! E_z=(fx+fx-fy-fy)/4
  EzNum_CB(:,:,:,iBlock) = - 0.25*(                                    &
       ( Flux_VX(By_,1:nI  ,1:nJ  ,1:nK)                            &
       + Flux_VX(By_,2:nI+1,1:nJ  ,1:nK)) / CellFace_DB(1,iBlock) - &
       ( Flux_VY(Bx_,1:nI  ,1:nJ  ,1:nK)                            &
       + Flux_VY(Bx_,1:nI  ,2:nJ+1,1:nK)) / CellFace_DB(2,iBlock))

end subroutine calc_electric_field
