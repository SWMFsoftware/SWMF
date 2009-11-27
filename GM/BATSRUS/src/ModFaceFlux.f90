module ModFaceFlux

  use ModProcMH,     ONLY: iProc
  use ModMain,       ONLY: x_, y_, z_, nI, nJ, nK, UseB, UseB0, cLimit, &
       iTest, jTest, kTest, ProcTest, BlkTest, DimTest
  use ModMain,       ONLY: UseRadDiffusion, UseParallelConduction
  use ModMain,       ONLY: UseBorisSimple                 !^CFG IF SIMPLEBORIS
  use ModMain,       ONLY: UseBoris => boris_correction   !^CFG IF BORISCORR
  use ModMultiFluid, ONLY: UseMultiIon, nIonFluid
  use ModGeometry,   ONLY: fAx_BLK, fAy_BLK, fAz_BLK, dx_BLK, dy_BLK, dz_BLK
  use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, true_cell

  use ModGeometry,   ONLY: UseCovariant, &                
       FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB, &     
       FaceArea2MinI_B, FaceArea2MinJ_B, FaceArea2MinK_B  

  use ModB0, ONLY:B0_DX, B0_DY, B0_DZ, B0_DGB ! input: face/cell centered B0

  use ModAdvance, ONLY:&
       LeftState_VX,  LeftState_VY,  LeftState_VZ,  & ! input: left  face state
       RightState_VX, RightState_VY, RightState_VZ, & ! input: right face state
       Flux_VX, Flux_VY, Flux_VZ,        & ! output: flux*Area
       VdtFace_x, VdtFace_y, VdtFace_z,  & ! output: cMax*Area for CFL
       EDotFA_X, EDotFA_Y, EDotFA_Z,     & ! output: E.Area !^CFG IF BORISCORR
       uDotArea_XI, uDotArea_YI, uDotArea_ZI,& ! output: U.Area for P source
       bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ,& ! output: B x Area for J
       UseRS7, UseTotalSpeed, UseElectronPressure,  &
       eFluid_                           ! index for electron fluid (nFluid+1)

  use ModPhysics, ONLY: AverageIonCharge, ElectronTemperatureRatio

  use ModMultiIon, ONLY: &
       Pe_X, Pe_Y, Pe_Z ! output: Pe for grad Pe in multi-ion MHD

  use ModHallResist, ONLY: UseHallResist, HallCmaxFactor, IonMassPerCharge_G, &
       IsNewBlockHall, hall_factor, get_face_current, set_ion_mass_per_charge

  !^CFG IF IMPLICIT BEGIN
  use ModRadDiffusion, ONLY: IsNewBlockRadDiffusion, get_radiation_energy_flux
  use ModHeatConduction, ONLY: IsNewBlockHeatConduction, get_heat_flux
  !^CFG END IMPLICIT

  use ModResistivity, ONLY: UseResistivity, Eta_GB  !^CFG IF DISSFLUX

  use ModVarIndexes
  use ModMultiFluid
  use ModNumConst

  implicit none

  ! Number of fluxes including pressure and energy fluxes
  integer, parameter :: nFlux=nVar+nFluid

  logical :: DoLf                !^CFG IF RUSANOVFLUX
  logical :: DoHll               !^CFG IF LINDEFLUX
  logical :: DoHlld              !^CFG IF HLLDFLUX
  logical :: DoAw                !^CFG IF AWFLUX
  logical :: DoRoeOld            !^CFG IF ROEFLUX
  logical :: DoRoe               !^CFG IF ROEFLUX
  logical :: DoGodunov           
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
  real :: Area = 0.0, Area2 = 0.0, AreaX, AreaY, AreaZ, Area2Min = 0.0

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

  ! Variables for normal resistivity
  real :: EtaJx, EtaJy, EtaJz, Eta = -1.0

  ! Variables needed for Hall resistivity
  real :: InvDxyz, HallCoeff, HallJx, HallJy, HallJz
  logical :: UseHallGradPe = .false., IsNewBlockGradPe = .true.
  real :: GradXPeNe, GradYPeNe, GradZPeNe

  ! Variables for diffusion solvers (radiation diffusion, heat conduction)
  real :: DiffCoef, EradFlux, RadDiffCoef, HeatFlux, HeatCondCoefNormal

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
    end select

  end subroutine face_flux_set_parameters

  !=========================================================================
  subroutine rotate_state_vectors

    ! Rotate the vector variables B0*, StateLeft_V(B*_), StateLeft_V(U*_)
    ! StateRight_V(B*_), StateRight_V(U*_) into normal and
    ! tangential components with respect to the face.
    ! Store the rotated vector components in scalar variables UnL, Ut1L, ....
    ! Also store the transformation for rotating back.
    ! Current implementation is for a single ion fluid.

    if(UseCovariant)then                           
       if(Normal_D(z_) < 0.5)then
          ! Tangent1 = Normal x (0,0,1)
          Tangent1_D(x_) =  Normal_D(y_)
          Tangent1_D(y_) = -Normal_D(x_)
          Tangent1_D(z_) = 0.0
          ! Tangent2 = Normal x Tangent1
          Tangent2_D(x_) =  Normal_D(z_)*Normal_D(x_)
          Tangent2_D(y_) =  Normal_D(z_)*Normal_D(y_)
          Tangent2_D(z_) = -Normal_D(x_)**2 - Normal_D(y_)**2
       else
          ! Tangent1 = Normal x (1,0,0)
          Tangent1_D(x_) = 0.0
          Tangent1_D(y_) =  Normal_D(z_)
          Tangent1_D(z_) = -Normal_D(y_)
          ! Tangent2 = Normal x Tangent1
          Tangent2_D(x_) = -Normal_D(y_)**2 - Normal_D(z_)**2
          Tangent2_D(y_) =  Normal_D(x_)*Normal_D(y_)
          Tangent2_D(z_) =  Normal_D(x_)*Normal_D(z_)
       end if
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
    else                                                
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
    end if                                          
  end subroutine rotate_state_vectors
  !==========================================================================
  subroutine rotate_flux_vector(FluxRot_V, Flux_V)
    real, intent(in)   :: FluxRot_V(:)
    real, intent(inout):: Flux_V(:)

    ! Rotate n,t1,t2 components back to x,y,z components
    if(UseCovariant)then                             
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
    else                                             
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
    end if                                           

  end subroutine rotate_flux_vector
  !===========================================================================
  subroutine calc_face_flux(DoResChangeOnly, iBlock)

    use ModAdvance,  ONLY: UseRS7,TypeFlux => FluxType
    use ModParallel, ONLY: &
         neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth
    use ModMain, ONLY: nIFace, nJFace, nKFace, &
         jMinFaceX, jMaxFaceX, kMinFaceX, kMaxFaceX, &
         iMinFaceY, iMaxFaceY, kMinFaceY,kMaxFaceY, &
         iMinFaceZ,iMaxFaceZ, jMinFaceZ, jMaxFaceZ, &
         UseHyperbolicDivb

    logical, intent(in) :: DoResChangeOnly
    integer, intent(in) :: iBlock
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('calc_face_flux', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if(DoTestMe)call print_values

    DoLf     = TypeFlux == 'Rusanov'     !^CFG IF RUSANOVFLUX
    DoHll    = TypeFlux == 'Linde'       !^CFG IF LINDEFLUX
    DoHlld   = TypeFlux == 'HLLD'        !^CFG IF HLLDFLUX
    DoAw     = TypeFlux == 'Sokolov'     !^CFG IF AWFLUX
    DoRoeOld = TypeFlux == 'RoeOld'      !^CFG IF ROEFLUX
    DoRoe    = TypeFlux == 'Roe'         !^CFG IF ROEFLUX
    DoGodunov= TypeFlux == 'Godunov'

    UseRS7 = DoRoe  ! This is always true for the current implementation

    UseLindeFix = (UseHyperbolicDivb &
         .or. DoHll      &               !^CFG IF LINDEFLUX
         .or. DoHllD     &               !^CFG IF HLLDFLUX
         .or. DoAw       &               !^CFG IF AWFLUX
          ).and.UseB
    ! Make sure that Hall MHD recalculates the magnetic field 
    ! in the current block that will be used for the Hall term
    IsNewBlockHall   = .true.
    IsNewBlockGradPe = .true.
    IsNewBlockRadDiffusion = .true.      !^CFG IF IMPLICIT
    IsNewBlockHeatConduction = .true.    !^CFG IF IMPLICIT

    if(UseResistivity) call set_resistivity(iBlock)      !^CFG IF DISSFLUX

    if(UseHallResist .and. UseMultiSpecies) &
         call set_ion_mass_per_charge(iBlock)

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
       call get_flux_x(1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
       if(nJ > 1) &
            call get_flux_y(iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
       if(nK > 1) &
            call get_flux_z(iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
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
      !-----------------------------------------------------------------------
      call set_block_values(iBlock, x_)

      do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax

         call set_cell_values_x

         DoTestCell = DoTestMe &
              .and. (iFace == iTest .or. iFace == iTest+1) &
              .and. jFace == jTest .and. kFace == kTest

         if(UseB0)then
            B0x = B0_DX(x_,iFace, jFace, kFace)
            B0y = B0_DX(y_,iFace, jFace, kFace)
            B0z = B0_DX(z_,iFace, jFace, kFace)
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

         VdtFace_x(iFace, jFace, kFace)       = CmaxDt*Area
         uDotArea_XI(iFace, jFace, kFace,:)   = Unormal_I*Area
         bCrossArea_DX(:, iFace, jFace, kFace)= bCrossArea_D
         EDotFA_X(iFace, jFace, kFace)        = Enormal*Area !^CFG IF BORISCORR

         if(UseMultiIon) Pe_X(iFace, jFace, kFace) = Pe

      end do; end do; end do
    end subroutine get_flux_x

    !==========================================================================

    subroutine get_flux_y(iMin,iMax,jMin,jMax,kMin,kMax)

      use ModAdvance, ONLY: State_VGB
      integer, intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      call set_block_values(iBlock, y_)

      do kFace = kMin, kMax; do jFace = jMin, jMax; do iFace = iMin, iMax

         call set_cell_values_y

         DoTestCell = DoTestMe .and. iFace == iTest .and. &
              (jFace == jTest .or. jFace == jTest+1) .and. kFace == kTest

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

         VdtFace_y(iFace, jFace, kFace)       = CmaxDt*Area
         uDotArea_YI(iFace, jFace, kFace, :)  = Unormal_I*Area
         bCrossArea_DY(:, iFace, jFace, kFace)= bCrossArea_D
         EDotFA_Y(iFace, jFace, kFace)        = Enormal*Area !^CFG IF BORISCORR

         if(UseMultiIon) Pe_Y(iFace, jFace, kFace) = Pe

      end do; end do; end do

    end subroutine get_flux_y

    !==========================================================================

    subroutine get_flux_z(iMin, iMax, jMin, jMax, kMin, kMax)

      use ModAdvance, ONLY: State_VGB
      integer, intent(in):: iMin, iMax, jMin, jMax, kMin, kMax
      !------------------------------------------------------------------------
      call set_block_values(iBlock, z_)

      do kFace = kMin, kMax; 
         do jFace = jMin, jMax; do iFace = iMin, iMax

         call set_cell_values_z

         DoTestCell = DoTestMe .and. iFace == iTest .and. &
              jFace == jTest .and. (kFace == kTest .or. kFace == kTest+1)
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

         VdtFace_z(iFace, jFace, kFace)       = CmaxDt*Area
         uDotArea_ZI(iFace, jFace, kFace, :)  = Unormal_I*Area
         bCrossArea_DZ(:, iFace, jFace, kFace)= bCrossArea_D
         EDotFA_Z(iFace, jFace, kFace)        = Enormal*Area !^CFG IF BORISCORR

         if(UseMultiIon) Pe_Z(iFace, jFace, kFace) = Pe

      end do; end do; end do
    end subroutine get_flux_z

  end subroutine calc_face_flux

  !===========================================================================
  subroutine set_block_values(iBlock, iDim)
    integer, intent(in) :: iBlock, iDim

    iBlockFace = iBlock
    iDimFace   = iDim

    ! For generalized coordinates the values below are calculated cell by cell
    if(UseCovariant) RETURN    

    ! Calculate face normal and area vectors for Cartesian grid
    Normal_D = 0.0; Normal_D(iDim) = 1.0
    NormalX = Normal_D(x_); NormalY = Normal_D(y_); NormalZ = Normal_D(z_)

    select case(iDim)
    case(x_)
       Area    = fAx_BLK(iBlockFace)
       AreaX   = Area; AreaY = 0.0; AreaZ = 0.0
       InvDxyz = 1./dx_BLK(iBlockFace)
    case(y_)
       Area    = fAy_BLK(iBlockFace)
       AreaY   = Area; AreaX = 0.0; AreaZ = 0.0
       InvDxyz = 1./dy_BLK(iBlockFace)
    case(z_)
       Area    = fAz_BLK(iBlockFace)
       AreaZ   = Area; AreaX = 0.0; AreaY = 0.0
       InvDxyz = 1./dz_BLK(iBlockFace)
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

    use ModMain, ONLY: Test_String

    iLeft = iFace - 1; jLeft = jFace; kLeft = kFace

    if(UseCovariant)then                  
       AreaX = FaceAreaI_DFB(x_, iFace, jFace, kFace, iBlockFace)
       AreaY = FaceAreaI_DFB(y_, iFace, jFace, kFace, iBlockFace)
       AreaZ = FaceAreaI_DFB(z_, iFace, jFace, kFace, iBlockFace)
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 0.5*FaceArea2MinI_B(iBlockFace))then
          ! The face is at the pole
          Normal_D(x_)= &
               x_BLK(iFace, jFace, kFace, iBlockFace)-&
               x_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(y_)= &
               y_BLK(iFace, jFace, kFace, iBlockFace)-&
               y_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(z_)= &
               z_BLK(iFace, jFace, kFace, iBlockFace)-&
               z_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D = Normal_D / sqrt(sum(Normal_D**2))
          Area =0.0
          Area2=0.0
          ! Store this in case pole diffusion is required
          Area2Min = FaceArea2MinI_B(iBlockFace)
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

    if(UseCovariant)then                   
       AreaX = FaceAreaJ_DFB(x_, iFace, jFace, kFace, iBlockFace)
       AreaY = FaceAreaJ_DFB(y_, iFace, jFace, kFace, iBlockFace)
       AreaZ = FaceAreaJ_DFB(z_, iFace, jFace, kFace, iBlockFace)
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 0.5*FaceArea2MinJ_B(iBlockFace))then
          !The face is at the pole
          Normal_D(x_)=x_BLK(iFace, jFace, kFace, iBlockFace)-&
                       x_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(y_)=y_BLK(iFace, jFace, kFace, iBlockFace)-&
                       y_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(z_)=z_BLK(iFace, jFace, kFace, iBlockFace)-&
                       z_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D=Normal_D/&
               sqrt(Normal_D(x_)**2+Normal_D(y_)**2+Normal_D(z_)**2)
          Area = 0.0
          Area2= 0.0
          ! Store this in case pole diffusion is required
          Area2Min = FaceArea2MinJ_B(iBlockFace)
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

    if(UseCovariant)then                  
       AreaX = FaceAreaK_DFB(x_, iFace, jFace, kFace, iBlockFace)
       AreaY = FaceAreaK_DFB(y_, iFace, jFace, kFace, iBlockFace)
       AreaZ = FaceAreaK_DFB(z_, iFace, jFace, kFace, iBlockFace)
       Area2 = AreaX**2 + AreaY**2 + AreaZ**2
       if(Area2 < 0.5*FaceArea2MinK_B(iBlockFace))then
          !The face is at the pole
          Normal_D(x_)=x_BLK(iFace, jFace, kFace, iBlockFace)-&
                       x_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(y_)=y_BLK(iFace, jFace, kFace, iBlockFace)-&
                       y_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D(z_)=z_BLK(iFace, jFace, kFace, iBlockFace)-&
                       z_BLK(iLeft, jLeft, kLeft, iBlockFace)
          Normal_D=Normal_D/&
               sqrt(Normal_D(x_)**2+Normal_D(y_)**2+Normal_D(z_)**2)
          Area = 0.0
          Area2= 0.0
          ! Store this in case pole diffusion is required
          Area2Min = FaceArea2MinK_B(iBlockFace)
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

    iRight= iFace; jRight = jFace; kRight = kFace

    IsBoundary = true_cell(iLeft, jLeft, kLeft, iBlockFace) &
         .neqv.  true_cell(iRight,jRight,kRight,iBlockFace)

    HallCoeff = -1.0
    if(UseHallResist) HallCoeff =                              &
         0.5*hall_factor(iDimFace, iFace, jFace, kFace, iBlockFace)* &
         ( IonMassPerCharge_G(iLeft, jLeft  ,kLeft)            &
         + IonMassPerCharge_G(iRight,jRight,kRight) )

    ! Calculate -grad(pe)/(n_e * e) term for Hall MHD if needed
    UseHallGradPe = HallCoeff > 0.0 .and. &
         (UseElectronPressure .or. ElectronTemperatureRatio > 0.0)

    Eta       = -1.0                                !^CFG IF DISSFLUX BEGIN
    if(UseResistivity) Eta = 0.5* &
         ( Eta_GB(iLeft, jLeft  ,kLeft,iBlockFace) &
         + Eta_GB(iRight,jRight,kRight,iBlockFace))  !^CFG END DISSFLUX

    if(UseCovariant)then
       NormalX = Normal_D(x_); NormalY = Normal_D(y_); NormalZ = Normal_D(z_)
       AreaX = Area*NormalX; AreaY = Area*NormalY; AreaZ = Area*NormalZ

       if(HallCoeff > 0.0 .or. Eta > 0.0 .or. UseParallelConduction) &
            InvDxyz = 1.0/sqrt(  &
            ( x_BLK(iRight,jRight,kRight, iBlockFace)          &
            - x_BLK(iLeft, jLeft  ,kLeft, iBlockFace))**2 +    &
            ( y_BLK(iRight,jRight,kRight, iBlockFace)          &
            - y_BLK(iLeft, jLeft  ,kLeft, iBlockFace))**2 +    &
            ( z_BLK(iRight,jRight,kRight, iBlockFace)          &
            - z_BLK(iLeft, jLeft  ,kLeft, iBlockFace))**2 )
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

    use ModAdvance, ONLY: DoReplaceDensity, State_VGB
    use ModCharacteristicMhd, ONLY: get_dissipation_flux_mhd
    use ModCoordTransform, ONLY: cross_product
    use ModMain, ONLY: UseHyperbolicDivb, SpeedHyp
    use ModFaceGradient, ONLY: get_face_gradient
    use ModImplicit, ONLY: UseSemiImplicit  !^CFG IF IMPLICIT

    real,    intent(out):: Flux_V(nFlux)

    real :: State_V(nVar)

    real :: Cmax
    real :: DiffBn_D(3), DiffE
    real :: EnLeft, EnRight, PeLeft, PeRight, Jx, Jy, Jz
    real :: uLeft_D(3), uRight_D(3) !,cDivBWave
    real :: dB0_D(3)

    real, save :: Pe_G(-1:nI+2,-1:nJ+2,-1:nK+2)
    real       :: GradPe_D(3)
    real       :: InvNumDens, Coef
    !-----------------------------------------------------------------------

    if(UseMultiSpecies .and. DoReplaceDensity)then
       StateLeft_V (Rho_)=sum( StateLeft_V(SpeciesFirst_:SpeciesLast_) )
       StateRight_V(Rho_)=sum( StateRight_V(SpeciesFirst_:SpeciesLast_) )
    end if

    ! Calculate current for the face if needed for (Hall) resistivity
    if(HallCoeff > 0.0 .or. Eta > 0.0) &
         call get_face_current(iDimFace,iFace,jFace,kFace,iBlockFace,Jx,Jy,Jz)

    if(Eta > 0.0)then                  !^CFG IF DISSFLUX BEGIN
       EtaJx = Eta*Jx
       EtaJy = Eta*Jy
       EtaJz = Eta*Jz
    end if                             !^CFG END DISSFLUX

    if(HallCoeff > 0.0)then
       HallJx = HallCoeff*Jx
       HallJy = HallCoeff*Jy
       HallJz = HallCoeff*Jz
    end if

    if(UseHallGradPe)then

       if(IsNewBlockGradPe)then
          ! Obtain electron pressure
          if(UseElectronPressure)then
             Pe_G = State_VGB(Pe_,:,:,:,iBlockFace) 
          elseif(IsMhd)then
             Coef = AverageIonCharge*ElectronTemperatureRatio
             Pe_G = State_VGB(p_,:,:,:,iBlockFace)*Coef/(1 + Coef)
          else
             Pe_G = sum(State_VGB(iPIon_I,:,:,:,iBlockFace),DIM=1) &
                  *ElectronTemperatureRatio
          end if
       end if

       ! Calculate face centered grad(Pe)
       call get_face_gradient(iDimFace, iFace, jFace, kFace, iBlockFace, &
            IsNewBlockGradPe, Pe_G, GradPe_D)

       ! Calculate 1/(n_e * e)
       if(UseMultiIon)then
          InvNumDens = HallCoeff/(0.5* &
               sum((StateLeft_V(iRhoIon_I)+StateRight_V(iRhoIon_I))/MassIon_I))
       else
          InvNumDens = HallCoeff/(0.5*(StateLeft_V(Rho_) + StateRight_V(Rho_)))
       end if

       ! Calculate grad(Pe)/(n_e * e)
       GradXPeNe = GradPe_D(x_)*InvNumDens
       GradYPeNe = GradPe_D(y_)*InvNumDens
       GradZPeNe = GradPe_D(z_)*InvNumDens

    end if

    !^CFG IF IMPLICIT BEGIN
    if(.not.UseSemiImplicit)then
       ! Initialize diffusion coefficient for time step restriction
       DiffCoef = 0.0

       if(UseRadDiffusion)then
          call get_radiation_energy_flux(iDimFace, iFace, jFace, kFace, &
               iBlockFace, StateLeft_V, StateRight_V, Normal_D, &
               RadDiffCoef, EradFlux)
          DiffCoef = DiffCoef + RadDiffCoef
       end if

       if(UseParallelConduction)then
          call get_heat_flux(iDimFace, iFace, jFace, kFace, iBlockFace, &
               StateLeft_V, StateRight_V, Normal_D, &
               HeatCondCoefNormal, HeatFlux)
          DiffCoef = DiffCoef + HeatCondCoefNormal
       end if
    end if
    !^CFG END IMPLICIT

    if(DoRoe)then
       if(UseB)then
          if(IsBoundary)then
             uLeft_D=StateLeft_V(Ux_:Uz_); uRight_D=StateRight_V(Ux_:Uz_)
          else
             !Since the divB source term is calculated using the
             !cell centered velocity, the numerical diffusion
             !for the normal magnetic field should be evaluated
             !in terms of the cell centered velocity too
             uLeft_D = &
                  State_VGB(RhoUx_:RhoUz_, iLeft, jLeft, kLeft, iBlockFace)/&
                  State_VGB(Rho_, iLeft, jLeft, kLeft, iBlockFace)
             uRight_D = &
                  State_VGB(RhoUx_:RhoUz_, iRight,jRight,kRight,iBlockFace)/&
                  State_VGB(Rho_, iRight,jRight,kRight,iBlockFace)
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
    end if
    if(UseRS7 .or. UseLindeFix)then
       ! Sokolov's algorithm
       ! Calculate the jump in the normal magnetic field vector
       DiffBn_D = Normal_D* &
            0.5*sum( (StateRight_V(Bx_:Bz_) - StateLeft_V(Bx_:Bz_))*Normal_D )

       ! Remove the jump in the normal magnetic field
       StateLeft_V(Bx_:Bz_)  =  StateLeft_V(Bx_:Bz_)  + DiffBn_D
       StateRight_V(Bx_:Bz_) =  StateRight_V(Bx_:Bz_) - DiffBn_D

       ! The energy jump is also modified by 
       ! 1/2(Br^2 - Bl^2) = 1/2(Br-Bl)*(Br+Bl)
       ! We store half of this in DiffE
       DiffE = &
            0.5*sum( (StateRight_V(Bx_:Bz_) + StateLeft_V(Bx_:Bz_))*DiffBn_D )

       DiffBb = sum(DiffBn_D**2)
    end if

    ! Calculate average state (used by most solvers and also by bCrossArea_D)
    State_V = 0.5*(StateLeft_V + StateRight_V)

    if(.not.DoGodunov &
         .and. .not.DoHlld &                                !^CFG IF HLLDFLUX
         )then
       ! All solvers, except HLLD, use left and right fluxes and avarage state
       call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
            StateLeftCons_V, FluxLeft_V, UnLeft_I, EnLeft, PeLeft)

       call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
            StateRightCons_V, FluxRight_V, UnRight_I, EnRight, PeRight)

       if(UseRS7)then
          call modify_flux(FluxLeft_V,UnLeft_I(1))
          call modify_flux(FluxRight_V,UnRight_I(1))
       end if
    end if

    if(UseB .and. (UseMultiIon .or. .not. IsMhd))then
       ! Calculate bCrossArea_D to be used for J in the J x B source term
       ! for the individual ion fluids in calc_sources.f90.
       ! The upwinded discretization of the current is J = sum(A x B) / V

       bCrossArea_D = cross_product(AreaX, AreaY, AreaZ, State_V(Bx_:Bz_))

       if(DoTestCell)then
          write(*,*)'bCrossArea_D=',bCrossArea_D
          write(*,*)'AreaX, AreaY, AreaZ=',AreaX, AreaY, AreaZ
          write(*,*)'State_V(Bx_:Bz_)=',State_V(Bx_:Bz_)
       end if
    end if

    if(DoLf)     call lax_friedrichs_flux       !^CFG IF RUSANOVFLUX
    if(DoHll)    call harten_lax_vanleer_flux   !^CFG IF LINDEFLUX
    if(DoHlld)   call hlld_flux                 !^CFG IF HLLDFLUX
    if(DoAw)     call artificial_wind           !^CFG IF AWFLUX
    if(DoRoeOld) call roe_solver(Flux_V)        !^CFG IF ROEFLUX
    if(DoRoe)    call roe_solver_new            !^CFG IF ROEFLUX
    if(DoGodunov)call godunov_flux

    if(UseRS7 .and. .not.DoRoe &
         .and. .not.DoHlld &                    !^CFG IF HLLDFLUX
         )then
       call stop_mpi('Second order RS7 is implemented for Roe solver only')
       !cDivBWave=max(abs(AreaX*UL_D(x_)+AreaY*UL_D(y_)+AreaZ*UL_D(z_)),&
       !     abs(AreaX*UR_D(x_)+AreaY*UR_D(y_)+AreaZ*UR_D(z_)))
       !Flux_V(Bx_) = Flux_V(Bx_) - cDivBWave*DiffBx
       !Flux_V(By_) = Flux_V(By_) - cDivBWave*DiffBy
       !Flux_V(Bz_) = Flux_V(Bz_) - cDivBWave*DiffBz
       !
       ! Fix the energy diffusion
       !Flux_V(Energy_) = Flux_V(Energy_) - cDivBWave*DiffE
    end if
    if(UseLindeFix)then

       if(UseHyperbolicDivb) then
          ! Overwrite the Flux of the Hyp field with the Lax-Friedrichs Flux
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

    ! Multiply Flux by Area. This is needed in div Flux in update_states_MHD
    Flux_V = Flux_V*Area

    ! Add diffusive flux across the pole if required
    if(Area == 0.0 .and. UsePoleDiffusion) Flux_V = Flux_V &
         - 0.5*sqrt(Area2Min)*Cmax*(StateRightCons_V - StateLeftCons_V)

    ! Increase maximum speed with resistive diffusion speed if necessary
    if(Eta > 0.0) CmaxDt = CmaxDt + 2*Eta*InvDxyz !^CFG IF DISSFLUX

    ! Increase maximum speed with diffusion speed if necessary
    !^CFG IF IMPLICIT BEGIN
    if(.not. UseSemiImplicit)then
       if(UseParallelConduction .or. UseRadDiffusion) &
            CmaxDt = CmaxDt + 2.0*DiffCoef*InvDxyz
    end if
    !^CFG END IMPLICIT

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
      cMaxDt = cMax

    end subroutine roe_solver_new
    !^CFG IF RUSANOVFLUX BEGIN
    !==========================================================================
    subroutine lax_friedrichs_flux

      real    :: Cmax_I(nFluid)
      integer :: iVar
      !----------------------------------------------------------------------
      call get_speed_max(State_V, B0x, B0y, B0z, Cmax_I = Cmax_I)

      if(UseTotalSpeed)then
         Cmax = maxval(Cmax_I)
         Flux_V = 0.5*(FluxLeft_V + FluxRight_V &
              - Cmax*(StateRightCons_V - StateLeftCons_V))
      else
         do iFluid = 1, nFluid
            Cmax = Cmax_I(iFluid); iRho=iRho_I(iFluid); iP=iP_I(iFluid) 
            do iVar = iRho_I(iFluid), iP_I(iFluid)
               Flux_V(iVar) = 0.5*(FluxLeft_V(iVar) + FluxRight_V(iVar) &
                    - Cmax*(StateRightCons_V(iVar) - StateLeftCons_V(iVar)))
            end do
            iVar = nVar + iFluid ! energy index
            Flux_V(iVar) = 0.5*(FluxLeft_V(iVar) + FluxRight_V(iVar) &
                 - Cmax*(StateRightCons_V(iVar) - StateLeftCons_V(iVar)))
         end do
         Cmax = Cmax_I(1)
      end if

      Unormal_I = 0.5*(UnLeft_I + UnRight_I)
      Enormal   = 0.5*(EnLeft + EnRight)                !^CFG IF BORISCORR
      if(UseMultiIon) &
           Pe   = 0.5*(PeLeft + PeRight)

    end subroutine lax_friedrichs_flux
    !^CFG END RUSANOVFLUX
    !^CFG IF LINDEFLUX BEGIN
    !==========================================================================
    subroutine harten_lax_vanleer_flux

      real, dimension(nFluid) :: CleftStateLeft_I,   CleftStateHat_I, &
           Cmax_I, CrightStateRight_I, CrightStateHat_I
      real :: Cleft, Cright, WeightLeft, WeightRight, Diffusion
      integer :: iVar
      !-----------------------------------------------------------------------

      call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
           Cleft_I =CleftStateLeft_I)

      call get_speed_max(StateRight_V, B0x, B0y, B0z, &
           Cright_I=CrightStateRight_I)

      call get_speed_max(State_V, B0x, B0y, B0z, &
           Cmax_I = Cmax_I, &
           Cleft_I = CleftStateHat_I, Cright_I = CrightStateHat_I)


      if(UseTotalSpeed)then
         Cmax   = maxval(Cmax_I)
         Cleft  =min(0.0, minval(CleftStateLeft_I), minval(CleftStateHat_I))
         Cright =max(0.0, maxval(CrightStateRight_I), maxval(CrightStateHat_I))

         WeightLeft  = Cright/(Cright - Cleft)
         WeightRight = 1.0 - WeightLeft
         Diffusion   = Cright*WeightRight

         Flux_V = &
              (WeightRight*FluxRight_V + WeightLeft*FluxLeft_V &
              - Diffusion*(StateRightCons_V - StateLeftCons_V))

         ! Weighted average of the normal speed
         Unormal_I = WeightRight*UnRight_I + WeightLeft*UnLeft_I

      else
         Cmax   = Cmax_I(1)
         do iFluid = 1, nFluid
            Cleft =min(0.0, CleftStateLeft_I(iFluid), CleftStateHat_I(iFluid))
            Cright=max(0.0,CrightStateRight_I(iFluid),CrightStateHat_I(iFluid))

            WeightLeft  = Cright/(Cright - Cleft)
            WeightRight = 1.0 - WeightLeft
            Diffusion   = Cright*WeightRight

            do iVar = iRho_I(iFluid), iP_I(iFluid)
               Flux_V(iVar) = &
                    (WeightRight*FluxRight_V(iVar)+WeightLeft*FluxLeft_V(iVar)&
                    - Diffusion*(StateRightCons_V(iVar)-StateLeftCons_V(iVar)))
            end do
            iVar = nVar + iFluid ! energy index
            Flux_V(iVar) = &
                 (WeightRight*FluxRight_V(iVar)+WeightLeft*FluxLeft_V(iVar) &
                 - Diffusion*(StateRightCons_V(iVar)-StateLeftCons_V(iVar)))

            ! Weighted average of the normal speed
            Unormal_I(iFluid) = &
                 WeightRight*UnRight_I(iFluid) + WeightLeft*UnLeft_I(iFluid)

            if(iFluid==1) Unormal_I(eFluid_) = &
                 WeightRight*UnRight_I(eFluid_) + WeightLeft*UnLeft_I(eFluid_)

         end do
      end if

      Enormal   = WeightRight*EnRight   + WeightLeft*EnLeft !^CFG IF BORISCORR
      if(UseMultiIon) &
           Pe   = WeightRight*PeRight   + WeightLeft*PeLeft

    end subroutine harten_lax_vanleer_flux
    !^CFG END LINDEFLUX
    !^CFG IF AWFLUX BEGIN
    !==========================================================================
    subroutine artificial_wind

      real, dimension(nFluid) :: Cleft_I, Cright_I, Cmax_I
      real :: Cleft, Cright, WeightLeft, WeightRight, Diffusion
      !-----------------------------------------------------------------------

      ! The propagation speeds are modified by the DoAw = .true. !
      call get_speed_max(State_V, B0x, B0y, B0z,  &
           Cleft_I = Cleft_I, Cright_I = Cright_I, Cmax_I = Cmax_I)

      Cmax   = maxval(Cmax_I)
      Cleft  = min(0.0, minval(Cleft_I))
      Cright = max(0.0, maxval(Cright_I))

      WeightLeft  = Cright/(Cright - Cleft)
      WeightRight = 1.0 - WeightLeft
      Diffusion   = Cright*WeightRight

      Flux_V = &
           (WeightRight*FluxRight_V + WeightLeft*FluxLeft_V &
           - Diffusion*(StateRightCons_V - StateLeftCons_V))

      ! Weighted average of the normal speed and electric field
      Unormal_I = WeightRight*UnRight_I + WeightLeft*UnLeft_I
      Enormal   = WeightRight*EnRight   + WeightLeft*EnLeft !^CFG IF BORISCORR
      if(UseMultiIon) &
           Pe   = WeightRight*PeRight   + WeightLeft*PeLeft

    end subroutine artificial_wind
    !^CFG END AWFLUX
    !^CFG IF HLLDFLUX BEGIN
    !==========================================================================
    subroutine hlld_flux

      use ModPhysics, ONLY: Inv_Gm1, gm1
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

      ! Resistivity                                  !^CFG IF DISSFLUX
      real :: FluxBx, FluxBy, FluxBz, B1x, B1y, B1z  !^CFG IF DISSFLUX

      real :: Tmp, B1n2, Bn2, SignBn

      real :: sL, CleftStateLeft_I(nFluid), CleftStateRight_I(nFluid)
      real :: sR, CrightStateLeft_I(nFluid), CrightStateRight_I(nFluid)
      
      integer, parameter :: ScalarMax_ = max(ScalarFirst_, ScalarLast_)
      real :: Scalar_V(ScalarFirst_:ScalarMax_)
      integer :: iVar
      !-----------------------------------------------------------------------

      ! This is the choice made in the hlld_tmp code. May not be the best.
      call get_speed_max(StateLeft_V,  B0x, B0y, B0z, &
           Cleft_I = CleftStateLeft_I, Cright_I = CrightStateLeft_I)

      call get_speed_max(StateRight_V, B0x, B0y, B0z, &
           Cleft_I = CleftStateRight_I, Cright_I = CrightStateRight_I)

      sL = min(minval(CleftStateLeft_I),  minval(CleftStateRight_I)) 
      sR = max(maxval(CrightStateLeft_I), maxval(CrightStateRight_I))

      Cmax   = max(sR, -sL)
      CmaxDt = Cmax

      if(DoTestCell)then
         write(*,*)'hlld: StateLeft =',StateLeft_V
         write(*,*)'hlld: StateRight=',StateRight_V
         write(*,*)'hlld: sL, sR    =',sL,sR
      endif

      if(sL >= 0.) then
         call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
              StateCons_V, Flux_V, Unormal_I, Enormal, Pe)
         if(UseRs7)call modify_flux(Flux_V, Unormal_I(1))
         RETURN
      end if

      if(sR <= 0.) then 
         call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
              StateCons_V, Flux_V, Unormal_I, Enormal, Pe)
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
         do iVar = ScalarFirst_, ScalarLast_
            Scalar_V(iVar) = StateLeft_V(iVar)*DsL*InvSLUn
         end do

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
         eL = Inv_Gm1*pL + PbL + 0.5*RhoL*(UnL**2 + Ut1L**2 + Ut2L**2)
         UdotB1L  = UnL*B1n + Ut1L *B1t1L  + Ut2L *B1t2L
         UdotB1L1 = Un *B1n + Ut1L1*B1t1L1 + Ut2L1*B1t2L1
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
         do iVar = ScalarFirst_, ScalarLast_
            Scalar_V(iVar) = StateRight_V(iVar)*DsR*InvSRUn
         end do

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
         eR = Inv_Gm1*pR + PbR + 0.5*RhoR*(UnR**2 + Ut1R**2 + Ut2R**2)
         UdotB1R  = UnR*B1n + Ut1R*B1t1R + Ut2R*B1t2R
         UdotB1R1 = Un*B1n + Ut1R1*B1t1R1 + Ut2R1*B1t2R1
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
      p = gm1*(e - 0.5*(B1n2 + B1t1**2 + B1t2**2 + Rho*(Un**2+Ut1**2+Ut2**2)))

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

      if(Eta > 0.0)then                          !^CFG IF DISSFLUX BEGIN
         ! Add flux corresponding to curl Eta.J to induction equation
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
      end if                                     !^CFG END DISSFLUX

    end subroutine hlld_flux
    !^CFG END HLLDFLUX

    !==========================================================================

    subroutine godunov_flux

      use ModAdvance,  ONLY: UseElectronPressure
      use ModExactRS,  ONLY: wR, wL, sample, pu_star, RhoL, RhoR, &
           pL, pR, UnL, UnR, UnStar, pStar
      use ModImplicit, ONLY: UseSemiImplicit  !^CFG IF IMPLICIT
      use ModPhysics,  ONLY: inv_gm1, g
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
      RhoL = StateLeft_V(Rho_)
      pL   = StateLeft_V(p_)
      RhoR = StateRight_V(Rho_)
      pR   = StateRight_V(p_)

      UnL  = sum( StateLeft_V(Ux_:Uz_) *Normal_D )
      UnR  = sum( StateRight_V(Ux_:Uz_)*Normal_D )

      if(UseWavePressure)then
         ! Increase maximum speed due to isotropic wave pressure
         ! To achieve this the pressure is increased by the value of the
         ! wave pressure
         pWaveL = (GammaWave - 1)*sum(StateLeft_V(WaveFirst_:WaveLast_))
         pWaveR = (GammaWave - 1)*sum(StateRight_V(WaveFirst_:WaveLast_))
         pL = pL + pWaveL
         pR = pR + pWaveR
      else
         !It is easier to add zero than using if(UseWavePressure)
         pWaveL = 0.0; pWaveR = 0.0
      end if
      if(UseElectronPressure)then
         ! StateLeft_V(p_) and StateRight_V(p_) are the ion pressure
         ! Add the electron pressure to the total pressure
         PeL = StateLeft_V(Pe_)
         PeR = StateRight_V(Pe_)
         pL = pL + PeL
         pR = pR + PeR
      else
         ! StateLeft_V(p_) and StateRight_V(p_) already include the electron
         ! pressure
         PeL = 0.0; PeR = 0.0
      end if

      !Take the parameters at the Contact Discontinuity (CD)
      call pu_star

      ! At strong shocks use the artificial wind scheme
      if((pStar > 2*pL .and. wL < 0.0).or.(pStar > 2*pR .and. wR > 0.0))then
         ! Temporary solution, should be the monotone numerical flux with 
         ! modified StateLeft_V and/or StateRight_V
         call get_physical_flux(StateLeft_V, B0x, B0y, B0z,&
              StateLeftCons_V, FluxLeft_V, UnLeft_I, EnLeft, PeLeft)
      
         call get_physical_flux(StateRight_V, B0x, B0y, B0z,&
              StateRightCons_V, FluxRight_V, UnRight_I, EnRight, PeRight)
      
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
      call sample(0.0, Rho, Un, p)

      ! In order the Riemann problem solution to be governed by the
      ! total pressure, the wave pressure should behave 
      ! adiabatically, with the same polytropic index as that for the gas

      Isothermal           = Rho/RhoSide
      Adiabatic            = Isothermal**g
      pWaveStar            = pWaveSide*Adiabatic
      PeStar               = PeSide*Adiabatic

      ! Since the total pressure is not less than the adiabatic one
      ! the difference below is positive
      pTotal               = p
      p                    = pTotal - pWaveStar - PeStar

      StateStar_V(Rho_)    = Rho
      StateStar_V(Ux_:Uz_) = StateStar_V(Ux_:Uz_) + (Un-UnSide)*Normal_D
      StateStar_V(P_)      = p
      do iVar=ScalarFirst_, ScalarLast_
         StateStar_V(iVar) = StateStar_V(iVar)*(Rho/RhoSide)
      end do

      ! Calculate flux 
      ! (1) calculate momenta
      StateStar_V(RhoUx_:RhoUz_) = StateStar_V(Ux_:Uz_) * Rho

      ! (2) take advective part of the flux
      Flux_V(1:nVar) = StateStar_V * Un 

      ! (3) add the pressure gradient
      ! also add the force due to the radiation and electron pressure gradient
      Flux_V(RhoUx_:RhoUz_) = Flux_V(RhoUx_:RhoUz_) + pTotal*Normal_D

      ! (4) energy flux: (e + p)*u
      ! also add the work done by the radiation and electron pressure gradient
      e = inv_gm1*p + 0.5*sum(StateStar_V(RhoUx_:RhoUz_)**2)/Rho
      Flux_V(Energy_) = (e + pTotal)*Un

      Cmax      = max(wR, -wL)
      CmaxDt    = Cmax
      Unormal_I = Un

      if(UseWavePressure)then
         GammaRatio = inv_gm1*(GammaWave - 1)
         Factor = (1.0 - GammaRatio) + GammaRatio*Adiabatic/Isothermal
         do iVar = WaveFirst_, WaveLast_
            Flux_V(iVar) = Factor*StateStar_V(iVar)*Un
         end do
      end if
      if(UseElectronPressure) &
           Flux_V(Pe_) = (Adiabatic/Isothermal)*StateStar_V(Pe_)*Un

      !^CFG IF IMPLICIT BEGIN
      if(.not.UseSemiImplicit)then
         if(UseRadDiffusion) Flux_V(Erad_) = Flux_V(Erad_) + EradFlux
      end if
      !^CFG END IMPLICIT

    end subroutine godunov_flux

    !==========================================================================

    subroutine write_test_info
      integer :: iVar
      !--------------------------------------------------------------------
      write(*,*)'Hat state for face=',iDimFace,&
           ' at I=',iFace,' J=',jFace,' K=',kFace
      write(*,*)'rho=',0.5*(StateLeft_V(Rho_)+StateRight_V(Rho_))
      write(*,*)'Un =',0.5*(StateLeft_V(U_+iDimFace)+StateRight_V(U_+iDimFace))
      write(*,*)'P  =',0.5*(StateLeft_V(P_)+StateRight_V(P_))
      write(*,*)'B  =', &
           0.5*(StateLeft_V(Bx_:Bz_) + StateRight_V(Bx_:Bz_)) + (/B0x,B0y,B0z/)
      write(*,*)'BB =', &
           sum( (0.5*(StateLeft_V(Bx_:Bz_) + StateRight_V(Bx_:Bz_)) &
           + (/B0x,B0y,B0z/))**2)

      write(*,*)'Fluxes for dir=',iDimFace,&
           ' at I=',iFace,' J=',jFace,' K=',kFace

      write(*,*)'Eigenvalue_maxabs=', Cmax
      write(*,*)'CmaxDt           =', CmaxDt
      do iVar = 1, nVar + nFluid
         write(*,'(a,a8,4(1pe13.5))') 'Var,F,F_L,F_R,dU=',&
              NameVar_V(iVar),&
              Flux_V(iVar), &
              FluxLeft_V(iVar), &
              FluxRight_V(iVar),&
              StateRightCons_V(iVar)-StateLeftCons_V(iVar)
      end do

    end subroutine write_test_info

  end subroutine get_numerical_flux

  !===========================================================================

  subroutine get_physical_flux(State_V, B0x, B0y, B0z, &
       StateCons_V, Flux_V, Un_I, En, Pe)

    use ModMultiFluid
    use ModMain,     ONLY: UseHyperbolicDivb, SpeedHyp2
    use ModImplicit, ONLY: UseSemiImplicit  !^CFG IF IMPLICIT
    use ModPhysics,  ONLY: gm1

    real,    intent(in) :: State_V(nVar)       ! input primitive state
    real,    intent(in) :: B0x, B0y, B0z       ! B0
    real,    intent(out):: StateCons_V(nFlux)  !conservative states with energy
    real,    intent(out):: Flux_V(nFlux)       ! fluxes for all states
    real,    intent(out):: Un_I(nFluid+1)      ! normal velocities
    real,    intent(out):: En                  ! normal electric field
    real,    intent(out):: Pe                  ! electron pressure for multiion

    real:: Hyp, Bx, By, Bz, FullBx, FullBy, FullBz, Bn, B0n, FullBn, Un, HallUn
    real:: FluxBx, FluxBy, FluxBz, Coef
    !--------------------------------------------------------------------------

    ! Calculate conservative state
    StateCons_V(1:nVar)  = State_V

    ! Make sure normal electric field is initialized
    En = 0.0

    if(UseMultiIon)then
       ! Pe has to be returned for multiion only
       if(UseElectronPressure)then
          Pe = State_V(Pe_)
       elseif(IsMhd)then
          Coef = AverageIonCharge*ElectronTemperatureRatio
          Pe = State_V(p_)*Coef/(1 + Coef)
       else
          Pe = sum(State_V(iPIon_I))*ElectronTemperatureRatio
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

    iFluid = 1
    if(IsMhd)then
       ! single ion fluid MHD (possibly with extra neutrals)
       if(UseBoris)then           !^CFG IF BORISCORR BEGIN
          call get_boris_flux
       else                       !^CFG END BORISCORR
          call get_mhd_flux
       end if                     !^CFG IF BORISCORR
    else
       ! If there is no MHD fluid, calculate fluxes for magnetic field
       ! together with hydro fluxes for the first fluid
       if(UseB)call get_magnetic_flux
       call select_fluid
       call get_hd_flux
    end if
    Un_I(1) = Un

    ! Calculate hydro fluxes for individual ion and neutral fluids
    do iFluid = 2, nFluid
       call select_fluid
       call get_hd_flux
       Un_I(iFLuid) = Un
    end do

    ! These terms are common for the induction equation
    ! If the first fluid is the total fluid, 
    ! the total energy density is also updated
    if(Eta > 0.0)then                          !^CFG IF DISSFLUX BEGIN
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
    end if                                     !^CFG END DISSFLUX

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

    !^CFG IF  IMPLICIT BEGIN
    if(.not.UseSemiImplicit)then
       if(UseRadDiffusion) Flux_V(Erad_) = Flux_V(Erad_) + EradFlux
       if(UseParallelConduction)then
          if(UseElectronPressure)then
             Flux_V(Pe_) = Flux_V(Pe_) + gm1*HeatFlux
          else
             Flux_V(p_) = Flux_V(p_) + gm1*HeatFlux
             Flux_V(Energy_) = Flux_V(Energy_) + HeatFlux
          end if
       end if
    end if
    !^CFG END IMPLICIT

    if(UseHyperbolicDivb)then
       Hyp  = State_V(Hyp_)

       Flux_V(Bx_:Bz_) = Flux_V(Bx_:Bz_) + Normal_D*Hyp
       Flux_V(Hyp_)    = SpeedHyp2*FullBn

       if(IsMhd) Flux_V(Energy_) = Flux_V(Energy_) + Bn*Hyp
    elseif(Hyp_ > 1)then
       Flux_V(Hyp_) = 0.0
    end if

    ! Set the normal electron velocity used for Hall MHD and/or 
    ! the electron pressure source term
    Un_I(eFluid_) = HallUn

  contains

    !^CFG IF BORISCORR BEGIN
    subroutine get_boris_flux

      use ModPhysics, ONLY: inv_gm1, Inv_C2light, InvClight

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e
      real :: B2, FullB2, pTotal, pTotal2, UDotB
      real :: Ex, Ey, Ez, E2Half
      integer :: iVar
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

      ! Calculate energy
      e = inv_gm1*p + 0.5*(Rho*(Ux**2 + Uy**2 + Uz**2) + B2)

      ! The full momentum contains the ExB/c^2 term:
      ! rhoU_Boris = rhoU - ((U x B) x B)/c^2 = rhoU + (U B^2 - B U.B)/c^2
      UDotB   = Ux*FullBx + Uy*FullBy + Uz*FullBz
      FullB2  = FullBx**2 + FullBy**2 + FullBz**2
      StateCons_V(RhoUx_)  = Rho*Ux + (Ux*FullB2 - FullBx*UdotB)*inv_c2LIGHT
      StateCons_V(RhoUy_)  = Rho*Uy + (Uy*FullB2 - FullBy*UdotB)*inv_c2LIGHT
      StateCons_V(RhoUz_)  = Rho*Uz + (Uz*FullB2 - FullBz*UdotB)*inv_c2LIGHT

      ! The full energy contains the electric field energy
      StateCons_V(Energy_) = e + E2Half

      ! Calculate some intermediate values for flux calculations
      pTotal  = p + 0.5*B2 + B0x*Bx + B0y*By + B0z*Bz
      pTotal2 = pTotal + E2Half

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

      ! f_i[scalar] = Un*scalar
      do iVar = ScalarFirst_, ScalarLast_
         Flux_V(iVar) = Un*State_V(iVar)
      end do

      HallUn = Un

    end subroutine get_boris_flux

    !==========================================================================
    !^CFG END BORISCORR

    subroutine get_mhd_flux

      use ModPhysics, ONLY: inv_gm1, g, inv_c2LIGHT
      use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure
      use ModWaves

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e
      real :: HallUx, HallUy, HallUz, InvRho
      real :: pAlfven
      real :: B2, B0B1, FullB2, pTotal, DpPerB
      real :: Gamma2                           !^CFG IF SIMPLEBORIS
      integer :: iVar

      real :: InvNumDens, StateTmp_V(nVar), UxPlus, UyPlus, UzPlus, UnPlus
      real, dimension(nIonFluid) :: NumDens_I, InvRho_I, Ux_I, Uy_I, Uz_I
      !-----------------------------------------------------------------------

      ! Extract primitive variables
      Rho     = State_V(Rho_)
      Ux      = State_V(Ux_)
      Uy      = State_V(Uy_)
      Uz      = State_V(Uz_)
      p       = State_V(p_)

      B2      = Bx**2 + By**2 + Bz**2

      ! Calculate energy
      e = inv_gm1*p + 0.5*(Rho*(Ux**2 + Uy**2 + Uz**2) + B2)

      ! Calculate conservative state
      StateCons_V(RhoUx_)  = Rho*Ux
      StateCons_V(RhoUy_)  = Rho*Uy
      StateCons_V(RhoUz_)  = Rho*Uz
      StateCons_V(Energy_) = e

      ! Calculate some intermediate values for flux calculations
      B0B1    = B0x*Bx + B0y*By + B0z*Bz
      pTotal  = p + 0.5*B2 + B0B1     ! For anisopressure, p = Pperp

      if(UseWavePressure) &
           pTotal = pTotal + (GammaWave-1.0)*sum(State_V(WaveFirst_:WaveLast_))
      if(UseElectronPressure) pTotal = pTotal + State_V(Pe_)

      ! Normal direction
      Un     = Ux*NormalX  + Uy*NormalY  + Uz*NormalZ

      if(UseMultiIon)then
         ! Calculate charge density averaged velocity U*Plus

         ! calculate number densities
         NumDens_I  = State_V(iRhoIon_I) / MassIon_I
         InvNumDens = 1.0/sum(NumDens_I)

         InvRho_I = 1.0/State_V(iRhoIon_I)
         Ux_I  = State_V(iUxIon_I)
         Uy_I  = State_V(iUyIon_I)
         Uz_I  = State_V(iUzIon_I)

         ! calculate the average positive charge velocity
         UxPlus = InvNumDens* sum(NumDens_I*Ux_I)
         UyPlus = InvNumDens* sum(NumDens_I*Uy_I)
         UzPlus = InvNumDens* sum(NumDens_I*Uz_I)

         UnPlus = UxPlus*NormalX + UyPlus*NormalY + UzPlus*NormalZ
      else
         ! For single ion fluid the mass and charge average is the same
         UxPlus = Ux
         UyPlus = Uy
         UzPlus = Uz
         UnPlus = Un
      end if

      ! f_i[rho] = Rho*U_i
      Flux_V(Rho_) = Rho*Un

      ! f_i[scalar] = Un*scalar
      do iVar = ScalarFirst_, ScalarLast_
         Flux_V(iVar) = Un*State_V(iVar)
      end do

      ! f_i[rhou_k] = u_i*u_k*rho - b_k*b_i -B0_k*b_i - B0_i*b_k + Ptotal*n_i
      Flux_V(RhoUx_) = Un*Rho*Ux - Bn*FullBx - B0n*Bx + pTotal*NormalX
      Flux_V(RhoUy_) = Un*Rho*Uy - Bn*FullBy - B0n*By + pTotal*NormalY
      Flux_V(RhoUz_) = Un*Rho*Uz - Bn*FullBz - B0n*Bz + pTotal*NormalZ

      ! f_i[b_k]=u_i*(b_k+B0_k) - u_k*(b_i+B0_i)
      if(HallCoeff > 0.0)then
         InvRho = 1.0/Rho
         HallUx = UxPlus - HallJx*InvRho
         HallUy = UyPlus - HallJy*InvRho
         HallUz = UzPlus - HallJz*InvRho

         HallUn = NormalX*HallUx + NormalY*HallUy + NormalZ*HallUz

         Flux_V(Bx_) = HallUn*FullBx - HallUx*FullBn
         Flux_V(By_) = HallUn*FullBy - HallUy*FullBn
         Flux_V(Bz_) = HallUn*FullBz - HallUz*FullBn
      else
         HallUn = UnPlus

         Flux_V(Bx_) = UnPlus*FullBx - UxPlus*FullBn
         Flux_V(By_) = UnPlus*FullBy - UyPlus*FullBn
         Flux_V(Bz_) = UnPlus*FullBz - UzPlus*FullBn
      end if
  
      ! f_i[Pe] = u_e,i*p_e
      if(UseElectronPressure)Flux_V(Pe_) = HallUn*State_V(Pe_)

      ! f_i[p] = u_i*p
      Flux_V(p_) = Un*p

      if(UseAnisoPressure)then
         ! f_i[rhou_k] = f_i[rho_k] + (ppar - pperp)bb for anisopressure
         FullB2 = FullBx**2 + FullBy**2 + FullBz**2
         DpPerB = (State_V(Ppar_) - State_V(Pperp_))*FullBn/max(1e-30,FullB2)
         Flux_V(RhoUx_) = Flux_V(RhoUx_) + FullBx*DpPerB
         Flux_V(RhoUy_) = Flux_V(RhoUy_) + FullBy*DpPerB
         Flux_V(RhoUz_) = Flux_V(RhoUz_) + FullBz*DpPerB
         ! f_i[Ppar] = u_i*Ppar, f_i[Pperp] = u_i*Pperp
         Flux_V(Ppar_)  = Un*State_V(Ppar_)
         Flux_V(Pperp_) = Un*State_V(Pperp_)
      end if

      ! f_i[e]=(u_i*(ptotal+e+(b_k*B0_k))-(b_i+B0_i)*(b_k*u_k))
      if(HallCoeff > 0.0) then
         Flux_V(Energy_) = &
              Un*(pTotal + e) &
              - FullBn*(HallUx*Bx + HallUy*By + HallUz*Bz)  &
              + (HallUn-Un)*(B2 + B0B1)
      else if(UseMultiIon)then
         Flux_V(Energy_) = &
              Un*(pTotal + e) &
              - FullBn*(UxPlus*Bx + UyPlus*By + UzPlus*Bz)  &
              + (UnPlus-Un)*(B2 + B0B1)
      else
         Flux_V(Energy_) = &
              Un*(pTotal + e) - FullBn*(Ux*Bx + Uy*By + Uz*Bz)     
      end if

      if(UseAlfvenSpeed)then
         AlfvenSpeed = FullBn/sqrt(Rho)

         do iVar = AlfvenSpeedPlusFirst_, AlfvenSpeedPlusLast_
            Flux_V(iVar) = Flux_V(iVar) + AlfvenSpeed*State_V(iVar) !!PLUS
         end do

         do iVar = AlfvenSpeedMinusFirst_, AlfvenSpeedMinusLast_
            Flux_V(iVar) = Flux_V(iVar) - AlfvenSpeed*State_V(iVar) !!MINUS
         end do
      end if

      !^CFG IF SIMPLEBORIS BEGIN
      if(UseBorisSimple)then
         ! Correct the momentum using the (1+VA2/c^2)
         Gamma2 = 1.0 + (FullBx**2 + FullBy**2 + FullBz**2)/Rho*inv_c2LIGHT
         StateCons_V(RhoUx_:RhoUz_) = StateCons_V(RhoUx_:RhoUz_)*Gamma2
      end if
      !^CFG END SIMPLEBORIS

    end subroutine get_mhd_flux

    !==========================================================================

    subroutine get_magnetic_flux

      ! Calculate magnetic flux for multi-ion equations 
      ! without a global ion fluid

      real :: NumDens_I(nIonFluid), InvNumDens
      real :: UxPlus, UyPlus, UzPlus
      real :: HallUx, HallUy, HallUz
      !-----------------------------------------------------------------------

      ! calculate number densities
      NumDens_I  = State_V(iRhoIon_I) / MassIon_I
      InvNumDens = 1.0/sum(NumDens_I)

      ! calculate positive charge velocity
      UxPlus = InvNumDens*sum(NumDens_I*State_V(iUxIon_I))
      UyPlus = InvNumDens*sum(NumDens_I*State_V(iUyIon_I))
      UzPlus = InvNumDens*sum(NumDens_I*State_V(iUzIon_I))

      if(HallCoeff > 0.0)then
         HallUx = UxPlus - HallJx*InvNumDens
         HallUy = UyPlus - HallJy*InvNumDens
         HallUz = UzPlus - HallJz*InvNumDens
      else
         HallUx = UxPlus
         HallUy = UyPlus
         HallUz = UzPlus
      end if

      HallUn = NormalX*HallUx + NormalY*HallUy + NormalZ*HallUz

      Flux_V(Bx_) = HallUn*FullBx - HallUx*FullBn
      Flux_V(By_) = HallUn*FullBy - HallUy*FullBn
      Flux_V(Bz_) = HallUn*FullBz - HallUz*FullBn

      if(DoTestCell)then
         write(*,*)'NumDens_I,InvNumDens=',NumDens_I,InvNumDens
         write(*,*)'UxyzPlus  =',UxPlus,UyPlus,UzPlus
         write(*,*)'HallUxyz  =',HallUx,HallUy,HallUz
         write(*,*)'FullBxyz  =',FullBx,FullBy,FullBz
         write(*,*)'B0x,y,z   =',B0x,B0y,B0z
         write(*,*)'Flux(Bxyz)=',Flux_V(Bx_:Bz_)
      end if

    end subroutine get_magnetic_flux

    !==========================================================================
    subroutine get_hd_flux

      use ModAdvance, ONLY: UseElectronPressure
      use ModPhysics, ONLY: inv_gm1
      use ModWaves

      ! Variables for conservative state and flux calculation
      real :: Rho, Ux, Uy, Uz, p, e, RhoUn, pTotal
      integer :: iVar
      !-----------------------------------------------------------------------
      ! Extract primitive variables
      Rho     = State_V(iRho)
      Ux      = State_V(iUx)
      Uy      = State_V(iUy)
      Uz      = State_V(iUz)
      p       = State_V(iP)

      ! Calculate energy
      e = inv_gm1*p + 0.5*Rho*(Ux**2 + Uy**2 + Uz**2)

      ! Calculate conservative state
      StateCons_V(iRhoUx)  = Rho*Ux
      StateCons_V(iRhoUy)  = Rho*Uy
      StateCons_V(iRhoUz)  = Rho*Uz
      StateCons_V(iEnergy) = e

      pTotal = p

      if(UseWavePressure) &
           pTotal = pTotal + (GammaWave-1.0)*sum(State_V(WaveFirst_:WaveLast_))
      if(UseElectronPressure) pTotal = pTotal + State_V(Pe_)

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

      ! f_i[scalar] = Un*scalar
      do iVar = ScalarFirst_, ScalarLast_
         Flux_V(iVar) = Un*State_V(iVar)
      end do

      ! Needed for adiabatic source term for electron pressure
      HallUn = Un

    end subroutine get_hd_flux

  end subroutine get_physical_flux

  !===========================================================================

  subroutine get_speed_max(State_V, B0x, B0y, B0z, cMax_I, cLeft_I, cRight_I)

    use ModMultiFluid, ONLY: select_fluid, iFluid, iRho, iUx, iUy, iUz, iP
    use ModMain, ONLY: Climit
    use ModWaves, ONLY: UseWavePressure, GammaWave, WaveEnergy

    real,    intent(in) :: State_V(nVar)
    real,    intent(in) :: B0x, B0y, B0z
    real, optional, intent(out) :: Cmax_I(nFluid)   ! max speed relative to lab
    real, optional, intent(out) :: Cleft_I(nFluid)  ! maximum left speed
    real, optional, intent(out) :: Cright_I(nFluid) ! maximum right speed

    real :: CmaxDt_I(nFluid)
    real :: UnLeft, UnRight                         !^CFG IF AWFLUX
    integer :: iVar
    !--------------------------------------------------------------------------

    if(DoAW)then                                 !^CFG IF AWFLUX BEGIN
       ! For AW flux UnLeft_I,UnRight_I are already set by get_physical_flux
       UnLeft = minval(UnLeft_I(1:nIonFluid))
       UnRight= maxval(UnRight_I(1:nIonFluid))
    end if                                       !^CFG END AWFLUX
    if(UseB)then
       if(UseBoris)then                             !^CFG IF BORISCORR BEGIN
          call get_boris_speed
       else                                         !^CFG END BORISCORR
          call get_mhd_speed
       endif                                        !^CFG IF BORISCORR    
    else
       iFluid = 1
       call get_hd_speed
    end if

    if(nFluid > 1)then

       do iFluid = 2, nFluid
          call select_fluid
          if(DoAw)then                           !^CFG IF AWFLUX BEGIN
             UnLeft = UnLeft_I(iFluid)
             UnRight= UnRight_I(iFluid)
          end if                                 !^CFG END AWFLUX
          call get_hd_speed
       end do

       ! For ion fluids the maximum speed is taken to be 
       ! the maximum of the individual HD speed and the total MHD speed
       do iFluid = 2, IonLast_
          if(present(Cmax_I))   CmaxDt_I(iFluid) = &
               max(CmaxDt_I(1), CmaxDt_I(iFluid))
          if(present(Cmax_I))   Cmax_I(iFluid)   = &
               max(Cmax_I(1),   Cmax_I(iFluid))
          if(present(Cleft_I))  Cleft_I(iFluid)  = &
               min(Cleft_I(1),  Cleft_I(iFluid))
          if(present(Cright_I)) Cright_I(iFluid) = &
               max(Cright_I(1), Cright_I(iFluid))
       end do

    end if

    ! Take maximum time step limit for all the fluids
    if (present(Cmax_I)) CmaxDt=maxval(CmaxDt_I)

    ! Limit propagation speeds if required
    if (Climit > 0.0) then
       if(present(Cmax_I))   Cmax_I   = min( Climit, Cmax_I)
       if(present(Cleft_I))  Cleft_I  = max(-Climit, Cleft_I)
       if(present(Cright_I)) Cright_I = min( Climit, Cright_I)

       if(present(Cmax_I))then
          ! If Climit has reduced the diffusion, then the block has to be
          ! advanced with the implicit scheme, so set CmaxDt to a huge number
          if(Climit < CmaxDt) CmaxDt = 1e30
       end if

    end if

  contains

    !^CFG IF BORISCORR BEGIN
    !========================================================================
    subroutine get_boris_speed

      use ModPhysics, ONLY: g, inv_c2LIGHT

      real :: InvRho, Sound2, FullBx, FullBy, FullBz
      real :: Alfven2, Alfven2Normal, Un, Fast2, Discr, Fast, Slow
      real :: GammaA2, GammaU2
      real :: UnBoris, Sound2Boris, Alfven2Boris, Alfven2NormalBoris
      !-----------------------------------------------------------------------
      InvRho = 1.0/State_V(Rho_)
      Sound2 = g*State_V(p_)*InvRho
      FullBx = State_V(Bx_) + B0x
      FullBy = State_V(By_) + B0y
      FullBz = State_V(Bz_) + B0z
      Alfven2= (FullBx**2 + FullBy**2 + FullBz**2)*InvRho

      Un = State_V(Ux_)*NormalX + State_V(Uy_)*NormalY + State_V(Uz_)*NormalZ
      Alfven2Normal = &
           InvRho*(NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz)**2

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
      Discr  = sqrt(max(0.0, Fast2**2 - 4.0*Sound2*Alfven2NormalBoris))

      ! Get fast and slow speeds multiplied with the face area
      Fast = sqrt( 0.5*(          Fast2 + Discr) )
      Slow = sqrt( 0.5*( max(0.0, Fast2 - Discr) ) )

      ! In extreme cases "slow" wave can be faster than "fast" wave
      ! so take the maximum of the two

      if(DoAw)then                                       !^CFG IF AWFLUX BEGIN
         Un           = min(UnRight, UnLeft)
         Cleft_I(1)   = min(Un*GammaA2 - Fast, Un - Slow)
         Un           = max(UnLeft, UnRight)
         Cright_I(1)  = max(Un*GammaA2 + Fast, Un + Slow)
         Cmax_I(1)    = max(Cright_I(1), -Cleft_I(1))
         CmaxDt_I(1)  = Cmax_I(1)
      else                                                !^CFG END AWFLUX
         UnBoris            = Un*GammaA2
         if(present(Cmax_I))then
            Cmax_I(1)   = max(abs(UnBoris) + Fast, abs(Un) + Slow)
            CmaxDt_I(1) = Cmax_I(1)
         end if
         if(present(Cleft_I))  Cleft_I(1)  = min(UnBoris - Fast, Un - Slow)
         if(present(Cright_I)) Cright_I(1) = max(UnBoris + Fast, Un + Slow)
      end if                                              !^CFG IF AWFLUX

    end subroutine get_boris_speed
    !^CFG END BORISCORR
    !========================================================================

    subroutine get_mhd_speed

      use ModMain,    ONLY: UseCurlB0
      use ModPhysics, ONLY: g, Inv_C2Light, ElectronTemperatureRatio
      use ModNumConst, ONLY: cPi
      use ModAdvance, ONLY: State_VGB, eFluid_, UseElectronPressure, &
           UseAnisoPressure

      real :: RhoU_D(3)
      real :: Rho, p, InvRho, Sound2, FullBx, FullBy, FullBz, FullBn, FullB2
      real :: Ppar, Pperp, BnInvB2, GammaPe
      real :: Alfven2, Alfven2Normal, Un, Fast2, Discr, Fast, FastDt, cWhistler
      real :: dB1dB1                                     !^CFG IF AWFLUX

      real :: FullBt, Rho1, cDrift, cHall, HallUnLeft, HallUnRight, &
           B1B0L, B1B0R

      character(len=*), parameter:: NameSub=NameMod//'::get_mhd_speed'
      !------------------------------------------------------------------------

      if(DoTestCell)write(*,*) NameSub,' State_V, B0=',State_V, B0x, B0y, B0z

      Rho    = State_V(Rho_)
      p      = State_V(p_)
      RhoU_D = Rho*State_V(Ux_:Uz_)
      if(.not. IsMhd)then
         do iFluid = 2, nIonFluid
            Rho1= State_V(iRhoIon_I(iFluid))
            Rho = Rho + Rho1
            p   = p   + State_V(iPIon_I(iFluid))
            RhoU_D = RhoU_D + Rho1*State_V(iUxIon_I(iFluid):iUzIon_I(iFluid))
         end do
         if(.not.UseElectronPressure) p = p * (1 + ElectronTemperatureRatio)
      end if

      if(UseElectronPressure) p = p + State_V(Pe_)

      InvRho = 1.0/Rho
      if(UseRS7)then
         Sound2 = (g*p+(g-1)*DiffBb)*InvRho
      else
         Sound2 = g*p*InvRho
      end if
      if(UseWavePressure)&
         Sound2 = Sound2 + GammaWave * (GammaWave - 1)*&
         sum(State_V(WaveFirst_:WaveLast_))*InvRho
     
      Un     = InvRho*sum( RhoU_D*Normal_D )

      FullBx = State_V(Bx_) + B0x
      FullBy = State_V(By_) + B0y
      FullBz = State_V(Bz_) + B0z
      if(DoAw)then                                       !^CFG IF AWFLUX BEGIN
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
      else                                               !^CFG END AWFLUX
         Alfven2= (FullBx**2 + FullBy**2 + FullBz**2)*InvRho
      end if                                             !^CFG IF AWFLUX
      if(UseCurlB0)then
         B1B0L = StateLeft_V(Bx_)*B0x &
              +  StateLeft_V(By_)*B0y &
              +  StateLeft_V(Bz_)*B0z
         B1B0R = StateRight_V(Bx_)*B0x &
              +  StateRight_V(By_)*B0y &
              +  StateRight_V(Bz_)*B0z
         Alfven2 = Alfven2 +(abs(B1B0L)-B1B0L+abs(B1B0R)-B1B0R)*InvRho
      end if

      FullBn = NormalX*FullBx + NormalY*FullBy + NormalZ*FullBz
      Alfven2Normal = InvRho*FullBn**2

      ! Calculate Fast speed for anisopressure.
      ! Formulae refer to V. B. Baranov, 1970
      if(UseAnisoPressure) FullB2 = FullBx**2+FullBy**2+FullBz**2
      if(UseAnisoPressure .and. FullB2 > 0)then
         GammaPe = 0.0 
         if(UseElectronPressure) GammaPe = g*State_V(Pe_) ! considering Pe
         Ppar = State_V(Ppar_)
         Pperp = State_V(Pperp_)
         BnInvB2 = FullBn**2/FullB2
         Sound2 = InvRho*(2*Pperp + (2*Ppar - Pperp)*BnInvB2 &
              + GammaPe)                        ! define Sound2 in this way
         Fast2 = Sound2 + Alfven2 
         Discr = sqrt(max(0.0, Fast2**2  &
              + 4*((Pperp*InvRho)**2*BnInvB2*(1 - BnInvB2) &  
              - (3*Ppar + GammaPe)*Pperp*InvRho**2*BnInvB2*(2-BnInvB2) &
              + (3*Ppar + GammaPe)*Ppar*(InvRho*BnInvB2)**2 &
              - (3*Ppar + GammaPe)*InvRho*Alfven2Normal)))
      else
         Fast2  = Sound2 + Alfven2
         Discr  = sqrt(max(0.0, Fast2**2 - 4*Sound2*Alfven2Normal))
      endif

      ! Fast speed multipled by the face area
      if(UseBorisSimple)then                         !^CFG IF SIMPLEBORIS BEGIN
         Fast = sqrt( 0.5*(Fast2 + Discr) &
              /       (1.0 + Alfven2*Inv_C2light) )
      else                                           !^CFG END SIMPLEBORIS
         Fast = sqrt( 0.5*(Fast2 + Discr) )
      end if                                         !^CFG IF SIMPLEBORIS


      ! Add whistler wave speed for the shortest wavelength 2 dx
      if(HallCoeff > 0.0) then
         ! Tangential component of B
         FullBt = sqrt(max(0.0, &
              (FullBx**2+FullBy**2+FullBz**2) - FullBn**2))
         ! Calculate Ln = d ln(Rho)/dx = (dRho/dx) / Rho
         Rho1 = State_VGB(Rho_,iLeft,jLeft,kLeft,iBlockFace)

         ! Calculate drift speed and whistler speed
         cDrift    = abs(FullBt)*2.0*abs(Rho1 - Rho)/(Rho1 + Rho)
         cWhistler = cPi*abs(FullBn)
         ! Take the faster speed
         cHall     = HallCoeff*InvDxyz*InvRho*max(cWhistler, cDrift)
         !cHall     = HallCoeff*InvDxyz*InvRho*cWhistler
         FastDt = Fast + cHall
         Fast   = Fast + HallCmaxFactor*cHall
      end if

      HallUnLeft  = UnLeft_I(eFluid_)
      HallUnRight = UnRight_I(eFluid_)

      if(DoAw)then                                   !^CFG IF AWFLUX BEGIN
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
      else                                           !^CFG END AWFLUX
         if(present(Cmax_I))then
            if(HallCoeff > 0.0)then
               Cmax_I(1)   = max(abs(Un), abs(HallUnLeft), abs(HallUnRight))
               CmaxDt_I(1) = Cmax_I(1) + FastDt
               Cmax_I(1)   = Cmax_I(1) + Fast
            else
               Cmax_I(1)   = abs(Un) + Fast
               CmaxDt_I(1) = Cmax_I(1)
            end if
         end if
         if(present(Cleft_I))  Cleft_I(1)  = Un - Fast
         if(present(Cright_I)) Cright_I(1) = Un + Fast
      end if                                         !^CFG IF AWFLUX

      if(DoTestCell)then
         if(UseCovariant)then
            write(*,*)NameSub,' AreaX,Y,Z =',AreaX, AreaY, AreaZ
            write(*,*)NameSub,' Area,Area2=',Area, Area2
            write(*,*)NameSub,' InvRho, RhoU_D=',InvRho, RhoU_D
         end if
         if(DoAw)then                               !^CFG IF AWFLUX BEGIN
            write(*,*)NameSub,' UnLeft=',  UnLeft
            write(*,*)NameSub,' UnRight=', UnRight
         end if                                     !^CFG END AWFLUX
         if(HallCoeff > 0.0) then
            write(*,*)NameSub,' HallCoeff=',   HallCoeff
            write(*,*)NameSub,' HallUnLeft=',  HallUnLeft
            write(*,*)NameSub,' HallUnRight=', HallUnRight
         end if
         write(*,*)NameSub,' Un=',Un
         write(*,*)NameSub,' Csound2=',Sound2
         write(*,*)NameSub,' Cfast2=', Fast2
         write(*,*)NameSub,' Discr2=', Discr**2
         write(*,*)NameSub,' Calfven=',sqrt(Alfven2)
         write(*,*)NameSub,' Calfven_normal=',sqrt(Alfven2Normal)
         write(*,*)NameSub,' Cfast=',Fast
         if(present(Cmax_I)) write(*,*)NameSub,' Cmax_I(1)=',Cmax_I(1)
      end if

    end subroutine get_mhd_speed
    !========================================================================
    subroutine get_hd_speed

      use ModAdvance, ONLY: UseElectronPressure, State_VGB
      use ModPhysics, ONLY: g

      real :: InvRho, Sound2, Sound, Un, p

      character(len=*), parameter:: NameSub=NameMod//'::get_hd_speed'
      !------------------------------------------------------------------------

      if(DoTestCell)write(*,*) NameSub,' State_V=',State_V(iRho:iP)

      ! Calculate sound speed and normal speed
      InvRho = 1.0/State_V(iRho)
      p = State_V(iP)
      if(UseElectronPressure) p = p + State_V(Pe_)
      Sound2 = g*p*InvRho

      if(UseWavePressure) &
           Sound2 = Sound2 + GammaWave*(GammaWave-1)*InvRho &
           *sum(State_V(WaveFirst_:WaveLast_))

      if(Sound2 <= 0.0)then
         write(*,*)NameSub,' negative Sound2=',Sound2
         write(*,*)NameSub,' iFluid, rho, p(face) =', &
              iFluid, State_V(iRho), State_V(iP)

         if(UseWavePressure)write(*,*)NameSub,' GammaWave, State(Waves):',&
              GammaWave, State_V(WaveFirst_:WaveLast_)

         if(UseElectronPressure)write(*,*)NameSub,' g,State_V(Pe_)=',&
              g,State_V(Pe_)

         write(*,*)NameSub,' rho, p(left) =', &
              State_VGB( (/Rho_,p_/),iLeft,jLeft,kLeft,iBlockFace)
         write(*,*)NameSub,' rho, p(right)=', &
              State_VGB( (/Rho_,p_/),iRight,jRight,kRight,iBlockFace)


         write(*,*)NameSub,' idim,i,j,k,BlockFace,iProc=', &
              iDimFace, iFace, jFace, kFace, iBlockFace, iProc
         write(*,*)NameSub,' xyz(right)=', &
              x_BLK(iFace,jFace,kFace,iBlockFace), &
              y_BLK(iFace,jFace,kFace,iBlockFace), &
              z_BLK(iFace,jFace,kFace,iBlockFace)
         call stop_mpi(NameSub//' negative soundspeed squared')
      end if

      Sound = sqrt(Sound2)
      Un    = sum(State_V(iUx:iUz)*Normal_D)


      if(DoAw)then                                   !^CFG IF AWFLUX BEGIN
         Cleft_I(iFluid)  = min(UnLeft, UnRight) - Sound
         Cright_I(iFluid) = max(UnLeft, UnRight) + Sound
         Cmax_I(iFluid)   = max(Cright_I(iFluid), -Cleft_I(iFluid))
         CmaxDt_I(iFluid) = Cmax_I(iFluid)
      else                                           !^CFG END AWFLUX
         if(present(Cmax_I))then
            Cmax_I(iFluid)   = abs(Un) + Sound
            CmaxDt_I(iFluid) = Cmax_I(iFluid)
         end if
         if(present(Cleft_I))  Cleft_I(iFluid)  = Un - Sound
         if(present(Cright_I)) Cright_I(iFluid) = Un + Sound
      end if                                         !^CFG IF AWFLUX

      if(DoTestCell)then
         write(*,*)NameSub,' Un     =',Un
         write(*,*)NameSub,' Csound =',Sound
         if(present(Cmax_I))write(*,*)NameSub,' Cmax=',Cmax_I(iFluid)
      end if

    end subroutine get_hd_speed

  end subroutine get_speed_max

end module ModFaceFlux

!^CFG IF ROEFLUX BEGIN
!==============================================================================
subroutine roe_solver(Flux_V)

  use ModFaceFlux, ONLY: &
       nFlux, IsBoundary, &
       StateLeft_V,  StateRight_V, FluxLeft_V, FluxRight_V, &
       StateLeftCons_V, StateRightCons_V, CmaxDt, Unormal_I, &
       nFluxMhd, RhoMhd_, RhoUn_, RhoUt1_, RhoUt2_, &
       B1n_, B1t1_, B1t2_, eMhd_, pMhd_, B0n, B0t1, B0t2, &
       UnL, Ut1L, Ut2L, B1nL, B1t1L, B1t2L, &
       UnR, Ut1R, Ut2R, B1nR, B1t1R, B1t2R, &
       rotate_state_vectors, rotate_flux_vector

  use ModVarIndexes, ONLY: Rho_, p_, Energy_, ScalarFirst_, ScalarLast_

  use ModPhysics,  ONLY: g,gm1,inv_gm1
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
  aL   = g*pL*RhoInvL

  RhoInvR = 1./RhoR
  BnR  = B0n+B1nR
  Bt1R = B0t1+B1t1R
  Bt2R = B0t2+B1t2R
  BbR  = BnR**2 + Bt1R**2 + Bt2R**2
  aR   = g*pR*RhoInvR

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
  aH   = g*pH*RhoInvH

  !if(aL<0.0)then
  !   write(*,*)'NEGATIVE aL Me, iDir, i, j, k, iBlockFace',&
  !        aL,iProc,iDimFace,i,j,k,&
  !        x_BLK(i,j,k,iBlockFace),&
  !        y_BLK(i,j,k,iBlockFace),&
  !        z_BLK(i,j,k,iBlockFace)
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
  eH   = pH*inv_gm1 + 0.5*RhoH*UuH + 0.5*Bb1H

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
  Tmp1 = max(1.00e-08, Bt1H**2 + Bt2H**2)
  Tmp1=sqrt(1./Tmp1)

  if (Tmp1 < 1.0e04) then
     BetaY = Bt1H*Tmp1
     BetaZ = Bt2H*Tmp1
  else
     BetaY = cSqrtHalf
     BetaZ = cSqrtHalf
  end if

  Tmp1 = CfH**2 - CsH**2
  if (Tmp1 > 1.0e-08) then
     AlphaF = max(0.00,(aH**2  - CsH**2)/Tmp1)
     AlphaS = max(0.00,(CfH**2 - aH**2 )/Tmp1)

     AlphaF = sqrt(AlphaF)
     AlphaS = sqrt(AlphaS)
  else if (BnH**2 * RhoInvH <= aH**2 ) then
     AlphaF = 1.00
     AlphaS = 0.00
  else
     AlphaF = 0.00
     AlphaS = 1.00
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
  Gamma1A2Inv = gm1 / aH**2

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

  Tmp1 = (UnR) - (UnL)
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
       (AlphaS*(gm1*UuH/2. - UnH*CsH) - &
       AlphaF*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,4) = Tmp1*(AlphaS*(-UnH*gm1+CsH))
  EigenvectorL_VV(3,4) = Tmp1*(-gm1*AlphaS*Ut1H + AlphaF*CfH*BetaY*SignBnH)
  EigenvectorL_VV(4,4) = Tmp1*(-gm1*AlphaS*Ut2H + AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorL_VV(5,4) = Tmp1*(-gm1*B1nH*AlphaS)
  EigenvectorL_VV(6,4) = Tmp1*(-AlphaF*BetaY*aH*RhoSqrtH-gm1*B1t1H*AlphaS)
  EigenvectorL_VV(7,4) = Tmp1*(-AlphaF*BetaZ*aH*RhoSqrtH - gm1*B1t2H*AlphaS)
  EigenvectorL_VV(8,4) = Tmp1*(gm1*AlphaS)

  ! Left eigenvector for Fast magnetosonic wave +
  EigenvectorL_VV(1,5) = Tmp1* &
       (AlphaF*(gm1*UuH/2. - UnH*CfH)+AlphaS*CsH*SignBnH* &
       (Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,5) = Tmp1*(AlphaF*(-UnH*gm1+CfH))
  EigenvectorL_VV(3,5) = Tmp1*(-gm1*AlphaF*Ut1H - AlphaS*CsH*BetaY*SignBnH)
  EigenvectorL_VV(4,5) = Tmp1*(-gm1*AlphaF*Ut2H - AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorL_VV(5,5) = Tmp1*(-gm1*B1nH*AlphaF)
  EigenvectorL_VV(6,5) = Tmp1*(AlphaS*BetaY*aH*RhoSqrtH - gm1*B1t1H*AlphaF)
  EigenvectorL_VV(7,5) = Tmp1*(AlphaS*BetaZ*aH*RhoSqrtH - gm1*B1t2H*AlphaF)
  EigenvectorL_VV(8,5) = Tmp1*(gm1*AlphaF)

  ! Left eigenvector for Slow magnetosonic wave -
  EigenvectorL_VV(1,6) = Tmp1* &
       (AlphaS*(gm1*UuH/2. + UnH*CsH) + &
       AlphaF*CfH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,6) = Tmp1*(AlphaS*(-UnH*gm1-CsH))
  EigenvectorL_VV(3,6) = Tmp1*(-gm1*AlphaS*Ut1H - AlphaF*CfH*BetaY*SignBnH)
  EigenvectorL_VV(4,6) = Tmp1*(-gm1*AlphaS*Ut2H - AlphaF*CfH*BetaZ*SignBnH)
  EigenvectorL_VV(5,6) = Tmp1*(-gm1*B1nH*AlphaS)
  EigenvectorL_VV(6,6) = Tmp1*(-AlphaF*BetaY*aH*RhoSqrtH-gm1*B1t1H*AlphaS)
  EigenvectorL_VV(7,6) = Tmp1*(-AlphaF*BetaZ*aH*RhoSqrtH-gm1*B1t2H*AlphaS)
  EigenvectorL_VV(8,6) = Tmp1*(gm1*AlphaS)

  ! Left eigenvector for Fast magnetosonic wave -
  EigenvectorL_VV(1,7) = Tmp1* &
       (AlphaF*(gm1*UuH/2. + UnH*CfH) - &
       AlphaS*CsH*SignBnH*(Ut1H*BetaY + Ut2H*BetaZ))
  EigenvectorL_VV(2,7) = Tmp1*(AlphaF*(-UnH*gm1-CfH))
  EigenvectorL_VV(3,7) = Tmp1*(-gm1*AlphaF*Ut1H + AlphaS*CsH*BetaY*SignBnH)
  EigenvectorL_VV(4,7) = Tmp1*(-gm1*AlphaF*Ut2H + AlphaS*CsH*BetaZ*SignBnH)
  EigenvectorL_VV(5,7) = Tmp1*(-gm1*B1nH*AlphaF)
  EigenvectorL_VV(6,7) = Tmp1*(AlphaS*BetaY*aH*RhoSqrtH-gm1*B1t1H*AlphaF)
  EigenvectorL_VV(7,7) = Tmp1*(AlphaS*BetaZ*aH*RhoSqrtH-gm1*B1t2H*AlphaF)
  EigenvectorL_VV(8,7) = Tmp1*(gm1*AlphaF)

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
  Tmp1=g*max(pL,pR) 

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
       AlphaS*(RhoH*UuH*0.5 + g*pH*inv_gm1+RhoH*UnH*CsH) &
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
       AlphaF*(RhoH*UuH*0.5 + g*pH*inv_gm1+RhoH*UnH*CfH) &
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
       AlphaS*(RhoH*UuH*0.5 + g*pH*inv_gm1-RhoH*UnH*CsH) &
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
       AlphaF*(RhoH*UuH*0.5 + g*pH*inv_gm1-RhoH*UnH*CfH) &
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
  !\
  ! Calculate the Roe Interface fluxes 
  ! F = A * 0.5 * [ F_L+F_R - sum_k(|lambda_k| * alpha_k * r_k) ]
  !/
  ! First get the diffusion: sum_k(|lambda_k| * alpha_k * r_k)
  !  Diffusion_V = matmul(abs(Eigenvalue_V)*DeltaWave_V, EigenvectorR_VV)
  do iFlux = 1, nFluxMhd
     Diffusion_V(iFlux) = &
          sum(abs(Eigenvalue_V)*DeltaWave_V*EigenvectorR_VV(:,iFlux))
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
!^CFG END ROEFLUX

!===========================================================================

subroutine calc_electric_field(iBlock)

  ! Calculate the total electric field which includes numerical resistivity
  ! This estimate averages the numerical fluxes to the cell centers 
  ! for sake of simplicity.

  use ModSize,       ONLY: nI, nJ, nK
  use ModVarIndexes, ONLY: Bx_,By_,Bz_
  use ModAdvance,    ONLY: Flux_VX, Flux_VY, Flux_VZ, Ex_CB, Ey_CB, Ez_CB
  use ModGeometry,   ONLY: fAx_BLK, fAy_BLK, fAz_BLK 
 
  implicit none
  integer, intent(in) :: iBlock
  !------------------------------------------------------------------------
  ! E_x=(fy+fy-fz-fz)/4
  Ex_CB(:,:,:,iBlock) = - 0.25*(                              &
       ( Flux_VY(Bz_,1:nI,1:nJ  ,1:nK  )                      &
       + Flux_VY(Bz_,1:nI,2:nJ+1,1:nK  )) / fAy_BLK(iBlock) - &
       ( Flux_VZ(By_,1:nI,1:nJ  ,1:nK  )                      &
       + Flux_VZ(By_,1:nI,1:nJ  ,2:nK+1)) / fAz_BLK(iBlock) )

  ! E_y=(fz+fz-fx-fx)/4
  Ey_CB(:,:,:,iBlock) = - 0.25*(                              &
       ( Flux_VZ(Bx_,1:nI  ,1:nJ,1:nK  )                      &
       + Flux_VZ(Bx_,1:nI  ,1:nJ,2:nK+1)) / fAz_BLK(iBlock) - &
       ( Flux_VX(Bz_,1:nI  ,1:nJ,1:nK  )                      &
       + Flux_VX(Bz_,2:nI+1,1:nJ,1:nK  )) / fAx_BLK(iBlock) )

  ! E_z=(fx+fx-fy-fy)/4
  Ez_CB(:,:,:,iBlock) = - 0.25*(                              &
       ( Flux_VX(By_,1:nI  ,1:nJ  ,1:nK)                      &
       + Flux_VX(By_,2:nI+1,1:nJ  ,1:nK)) / fAx_BLK(iBlock) - &
       ( Flux_VY(Bx_,1:nI  ,1:nJ  ,1:nK)                      &
       + Flux_VY(Bx_,1:nI  ,2:nJ+1,1:nK)) / fAy_BLK(iBlock))

end subroutine calc_electric_field

