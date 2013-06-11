!This code is a copyright protected software (c) 2002- University of Michigan
module ModHallResist

  use ModSize, ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxDim

  implicit none

  private !except

  ! Logical for adding the Biermann battery term
  logical, public :: UseBiermannBattery = .false.

  ! Logical for adding hall resistivity
  logical, public:: UseHallResist=.false.
  logical, public:: IsNewBlockCurrent=.true.
  ! Coefficient for taking whistler wave speed into account
  real, public:: HallCmaxFactor = 1.0

  ! Non-diagonal part (Hall) resistivity with an arbitrary factor
  real, public:: HallFactorMax = 1.0

  ! Ion mass per charge may depend on space and time for multispecies
  real, public, allocatable:: IonMassPerCharge_G(:,:,:)
  real:: IonMassPerChargeCoef

  ! Arrays for the implicit preconditioning
  real, public, allocatable :: HallJ_CD(:,:,:,:), &
       BxPerN_G(:,:,:),ByPerN_G(:,:,:),BzPerN_G(:,:,:)

  ! Name/shape of the region where Hall effect is used
  character(len=20), public :: NameHallRegion ='all'

  ! Parameters for exclusion of pole region in spherical grids
  real, public :: PoleAngleHall=-1.0, dPoleAngleHall=-1.0

  ! Parameters for exclusion of inner boundary (if present)
  real, public :: rInnerHall=-1.0, DrInnerHall=0.0

  ! Parameters for the center location of Hall region
  real, public :: X0Hall=0.0, Y0Hall=0.0, Z0Hall=0.0

  ! Parameters for spherical Hall region
  real, public :: rSphereHall=-1.0, DrSphereHall=-1.0

  ! Parameters for box Hall region
  real, public :: xSizeBoxHall=-1.0, ySizeBoxHall=-1.0, zSizeBoxHall=-1.0
  real, public :: DxSizeBoxHall=-1.0, DySizeBoxHall=-1.0, DzSizeBoxHall=-1.0

  ! Jacobian matrix for general grid: Dgencoord/Dcartesian
  real, public :: DgenDxyz_DDFD(MaxDim,MaxDim,nI+1,nJ+1,nK+1,MaxDim)
  real, public :: DgenDxyz_DDC(MaxDim,MaxDim,nI,nJ,nK)

  ! Public methods
  public :: init_hall_resist, get_face_current, hall_factor, test_face_current
  public :: set_ion_mass_per_charge
  public :: set_ion_mass_per_charge_point
  public :: set_block_jacobian_cell           

  ! Magnetic field with 3rd order accurate ghost cells
  real:: b_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

  ! Inverse of cell size
  real :: InvDx, InvDy, InvDz

  ! Local variables for Hall regions
  real :: TanSqr1 = -1.0, TanSqr2 = -1.0                  ! pole
  real :: rSqrInner1 = -1.0, rSqrInner2 = -1.0            ! inner body
  real :: rSqrSphere1 = -1.0, rSqrSphere2 = -1.0          ! spherical region
  real :: xSizeBox1=-1.0, ySizeBox1=-1.0, zSizeBox1=-1.0  ! box region
  real :: xSizeBox2=-1.0, ySizeBox2=-1.0, zSizeBox2=-1.0  ! box region

  save

contains
  !============================================================================
  subroutine init_hall_resist
    use ModConst,   ONLY: cMu, cElectronCharge
    use ModSize,    ONLY: nI, nJ, nK, MaxDim
    use ModPhysics, ONLY: IonMassPerCharge, Si2No_V, No2Si_V, UnitX_, &
         UnitRho_, UnitT_, UnitB_, UnitMass_, UnitCharge_, UnitU_, UnitJ_

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_hall_resist'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if (DoTestMe) then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>>>>>>>> HALL Resistivity Parameters <<<<<<<<<<'
       write(*,*)
       write(*,*) 'HallFactorMax    = ', HallFactorMax
       write(*,*) 'HallCmaxFactor   = ', HallCmaxFactor
       write(*,*) 'IonMassPerCharge = ', IonMassPerCharge
       ! Omega_Bi=B0/IonMassPerCharge'
       write(*,*)
       write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
       write(*,*) ''
    end if

    if(.not.allocated(HallJ_CD)) allocate(&
         HallJ_CD(nI,nJ,nK,MaxDim), &
         BxPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         ByPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         BzPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         IonMassPerCharge_G(0:nI+1,0:nJ+1,0:nK+1) )

    IonMassPerCharge_G = IonMassPerCharge

    ! This is used in combination with normalized density
    ! divided by SI charge density.
    IonMassPerChargeCoef = &
         Si2No_V(UnitX_)**3 / (cElectronCharge*Si2No_V(UnitCharge_))

    rSqrInner1 = -1.0
    rSqrInner2 = -1.0
    TanSqr1 = -1.0
    TanSqr2 = -1.0

    if(PoleAngleHall > 0.0)then
       TanSqr1 = tan(PoleAngleHall)**2
       TanSqr2 = tan(PoleAngleHall+dPoleAngleHall)**2
    end if

    if(rInnerHall > 0.0)then
       rSqrInner1 = rInnerHall**2
       rSqrInner2 = (rInnerHall + DrInnerHall)**2
    endif

    select case(NameHallRegion)
    case('sphere')
       rSqrSphere1 = rSphereHall**2
       rSqrSphere2 = (rSphereHall+DrSphereHall)**2
    case('box')
       xSizeBox1 = xSizeBoxHall
       ySizeBox1 = ySizeBoxHall
       zSizeBox1 = zSizeBoxHall
       xSizeBox2 = xSizeBoxHall + 2*DxSizeBoxHall
       ySizeBox2 = ySizeBoxHall + 2*DySizeBoxHall
       zSizeBox2 = zSizeBoxHall + 2*DzSizeBoxHall
    case('all')
       ! Do nothing
    case default
       call stop_mpi('ERROR in init_hall_resist: NameHallRegion=' &
            //NameHallRegion)
    end select

  end subroutine init_hall_resist

  !=========================================================================
  subroutine set_ion_mass_per_charge(iBlock)

    use ModAdvance, ONLY: State_VGB

    ! Set IonMassPerCharge_G based on average mass
    integer, intent(in) :: iBlock

    integer :: i, j, k
    !-------------------------------------------------------------------------

    ! Multiply IonMassPerCharge_G by average ion mass = rho_total / n_total

    do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
       call set_ion_mass_per_charge_point(State_VGB(:,i,j,k,iBlock), &
            IonMassPerCharge_G(i,j,k))
    end do; end do; end do

  end subroutine set_ion_mass_per_charge

  !===========================================================================

  subroutine set_ion_mass_per_charge_point(State_V, IonMassPerChargeOut)

    use ModAdvance,    ONLY: UseIdealEos
    use ModVarIndexes, ONLY: nVar, Rho_, &
         UseMultiSpecies, SpeciesFirst_, SpeciesLast_, MassSpecies_V
    use ModMultiFluid, ONLY: UseMultiIon, iRhoIon_I, MassIon_I,ChargeIon_I
    use ModPhysics,    ONLY: IonMassPerCharge
    use ModUser,       ONLY: user_material_properties

    real, intent(in) :: State_V(nVar)
    real, intent(out):: IonMassPerChargeOut

    real :: zAverage, NatomicSi
    !--------------------------------------------------------------------------

    if(.not.UseIdealEos)then
       call user_material_properties(State_V, &
            AverageIonChargeOut = zAverage, NatomicOut = NatomicSi)

       ! Avoid using small zAverage, since then we will generate magnetic
       ! field with the Biermann Battery term based numerical errors.
       zAverage = max(zAverage, 1.0)

       IonMassPerChargeOut = IonMassPerChargeCoef*State_V(Rho_) &
            /(zAverage*NatomicSi)

    elseif(UseMultiSpecies)then
       IonMassPerChargeOut = IonMassPerCharge*State_V(Rho_) &
            / sum(State_V(SpeciesFirst_:SpeciesLast_)/MassSpecies_V)

    elseif(UseMultiIon)then
       ! Get mass density per total number denisity
       IonMassPerChargeOut = IonMassPerCharge*sum(State_V(iRhoIon_I)) &
            / sum(State_V(iRhoIon_I)*ChargeIon_I / MassIon_I)

    else
       IonMassPerChargeOut = IonMassPerCharge

    end if

  end subroutine set_ion_mass_per_charge_point

  !============================================================================

  subroutine set_block_jacobian_face(iBlock)
    use ModCoordTransform, ONLY: inverse_matrix
    use BATL_lib, ONLY: Xyz_DGB

    integer, intent(in):: iBlock

    ! Dxyz/Dgen matrix for one cell
    real:: DxyzDgen_DD(MaxDim,MaxDim)

    ! Transverse gradients
    real:: TransGrad_DDG(MaxDim,MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Cell center coordinates for this block
    real:: Xyz_DG(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Indexes
    integer:: i, j, k

    !coeff of Ui+2 and Ui+1 to get normal derivative
    real, parameter:: fp2 = -1./24.0, fp1 = 9.0/8.0 
    !coeff of Ui+2 and Ui+1 for transverse derivatives
    real, parameter:: dp2 = -1./12.0, dp1 = 2.0/3.0 
    !coeff to average transverse derivatives
    real, parameter:: ap2 = -1./16.0, ap1 = 9.0/16. 

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='set_block_jacobian_face'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dGencoord/dCartesian matrix

    Xyz_DG(:,:,:,:) = Xyz_DGB(:,:,:,:,iBlock)

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=1,nI
       TransGrad_DDG(:,1,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i+1,j,k) - Xyz_DG(:,i-1,j,k)) &
            + dp2* (Xyz_DG(:,i+2,j,k) - Xyz_DG(:,i-2,j,k)))
    end do; end do; end do

    do k=MinK,MaxK; do j=1,nJ; do i=MinI,MaxI
       TransGrad_DDG(:,2,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i,j+1,k) - Xyz_DG(:,i,j-1,k)) &
            + dp2* (Xyz_DG(:,i,j+2,k) - Xyz_DG(:,i,j-2,k)))
    end do; end do; end do

    do k=1,nK; do j=MinJ,MaxJ; do i=MinI,MaxI
       TransGrad_DDG(:,3,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i,j,k+1) - Xyz_DG(:,i,j,k-1)) &
            + dp2* (Xyz_DG(:,i,j,k+2) - Xyz_DG(:,i,j,k-2)))
    end do; end do; end do

    !gen1 face
    do k=1,nK; do j=1,nJ; do i=1,nI+1
       !dxyzdgen along gen1 face
       DxyzDgen_DD(:,1) = InvDx* &
            (  fp1*(Xyz_DG(:,i  ,j,k) - Xyz_DG(:,i-1,j,k)) &
            +  fp2*(Xyz_DG(:,i+1,j,k) - Xyz_DG(:,i-2,j,k)))
       DxyzDgen_DD(:,2) = InvDy* &
            ( ap1*( TransGrad_DDG(:,2,i  ,j,k) + TransGrad_DDG(:,2,i-1,j,k)) &
            + ap2*( TransGrad_DDG(:,2,i+1,j,k) + TransGrad_DDG(:,2,i-2,j,k)))
       DxyzDgen_DD(:,3) = InvDz* &
            ( ap1*( TransGrad_DDG(:,3,i  ,j,k) + TransGrad_DDG(:,3,i-1,j,k)) &
            + ap2*( TransGrad_DDG(:,3,i+1,j,k) + TransGrad_DDG(:,3,i-2,j,k)))
       DgenDxyz_DDFD(:,:,i,j,k,1) = &
            inverse_matrix(DxyzDgen_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

    !gen2 face
    do k=1,nK; do j=1,nJ+1; do i=1,nI
       !dxyzdgen along gen2 face
       DxyzDgen_DD(:,1) = InvDx* &
            ( ap1*( TransGrad_DDG(:,1,i,j  ,k) + TransGrad_DDG(:,1,i,j-1,k)) &
            + ap2*( TransGrad_DDG(:,1,i,j+1,k) + TransGrad_DDG(:,1,i,j-2,k)))
       DxyzDgen_DD(:,2) = InvDy* &
            (  fp1*(Xyz_DG(:,i,j  ,k) - Xyz_DG(:,i,j-1,k)) &
            +  fp2*(Xyz_DG(:,i,j+1,k) - Xyz_DG(:,i,j-2,k)))
       DxyzDgen_DD(:,3) = InvDz* &
            ( ap1*( TransGrad_DDG(:,3,i,j  ,k) + TransGrad_DDG(:,3,i,j-1,k)) &
            + ap2*( TransGrad_DDG(:,3,i,j+1,k) + TransGrad_DDG(:,3,i,j-2,k)))

       DgenDxyz_DDFD(:,:,i,j,k,2) = &
            inverse_matrix(DxyzDgen_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

    !gen3 face
    do k=1,nK+1; do j=1,nJ; do i=1,nI

       !dxyzdgen along gen1 face
       DxyzDgen_DD(:,1) = InvDx* &
            ( ap1*( TransGrad_DDG(:,1,i,j,k  ) + TransGrad_DDG(:,1,i,j,k-1)) &
            + ap2*( TransGrad_DDG(:,1,i,j,k+1) + TransGrad_DDG(:,1,i,j,k-2)))
       DxyzDgen_DD(:,2) = InvDy* &
            ( ap1*( TransGrad_DDG(:,2,i,j,k  ) + TransGrad_DDG(:,2,i,j,k-1)) &
            + ap2*( TransGrad_DDG(:,2,i,j,k+1) + TransGrad_DDG(:,2,i,j,k-2)))
       DxyzDgen_DD(:,3) = InvDz* &
            (  fp1*(Xyz_DG(:,i,j,k  ) - Xyz_DG(:,i,j,k-1)) &
            +  fp2*(Xyz_DG(:,i,j,k+1) - Xyz_DG(:,i,j,k-2)))

       DgenDxyz_DDFD(:,:,i,j,k,3) = &
            inverse_matrix(DxyzDgen_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

  end subroutine set_block_jacobian_face

  !============================================================================

  subroutine set_block_jacobian_cell(iBlock)
    use ModMain, ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: inverse_matrix
    use BATL_lib, ONLY: Xyz_DGB

    integer, intent(in):: iBlock
    real:: InvDx1Half, InvDx2Half, InvDx3Half
    real:: DxyzDgen_DD(MaxDim, MaxDim)
    integer:: i,j,k
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='set_block_jacobian'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dCartesian/dGencoord matrix
    
    InvDx1Half = InvDx*0.5
    InvDx2Half = InvDy*0.5
    InvDx3Half = InvDz*0.5

    ! Get the dCartesian/dGencoord matrix with finite differences
    do k=1,nK; do j=1,nJ; do i=1,nI
       DxyzDgen_DD(:,1) = InvDx1Half &
            *(Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock))

       DxyzDgen_DD(:,2) = InvDx2Half &
            *(Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock))

       DxyzDgen_DD(:,3) = InvDx3Half &
            *(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))

       DgenDxyz_DDC(:,:,i,j,k) = &
            inverse_matrix(DxyzDgen_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

  end subroutine set_block_jacobian_cell
 
  !============================================================================

  subroutine get_face_current(iDir, i, j, k, iBlock, Jx, Jy, Jz)

    use ModAdvance,      ONLY: State_VGB
    use ModFaceGradient, ONLY: set_block_field3
    use ModProcMH,       ONLY: iProc
    use ModMain,         ONLY: nI, nJ, nK, x_, y_, z_, &
         iTest, jTest, kTest, BlkTest, ProcTest
    use BATL_lib,        ONLY: IsCartesianGrid, IsRzGeometry, &
         CellSize_DB, Xyz_DGB, DiLevelNei_IIIB
    use ModParallel,     ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use ModVarIndexes,   ONLY: Bx_, Bz_
    implicit none

    real, parameter :: cTwoThird = 2.0/3.0, cFourThird = 4.0/3.0

    character(len=*), parameter :: NameSub='get_face_current'

    integer, intent(in):: iDir, i, j, k, iBlock
    real, intent(out)  :: Jx, Jy, Jz

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz

    real :: b1_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    logical :: DoTest, DoTestMe
    integer :: i1,j1,k1
    !-------------------------------------------------------------------------
    if(iProc==ProcTest.and.iBlock==BlkTest.and. &
             (i==iTest.or.i==iTest+1.and.iDir==x_) .and. &
             (j==jTest.or.j==jTest+1.and.iDir==y_) .and. &
             (k==kTest.or.k==kTest+1.and.iDir==z_))then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTestMe = .false.; DoTest =.false.
    end if

    Jx = 0.0; Jy = 0.0; Jz = 0.0

    InvDx = 1.0/CellSize_DB(x_,iBlock)
    InvDy = 1.0/CellSize_DB(y_,iBlock)
    InvDz = 1.0/CellSize_DB(z_,iBlock)

    if( IsNewBlockCurrent ) then
       b_DG = State_VGB(bx_:bz_,:,:,:,iBlock)
       call set_block_field3(iBlock, 3, b1_DG, b_DG)
       if(.not.IsCartesianGrid) call set_block_jacobian_face(iBlock)

       IsNewBlockCurrent = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1
    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB(-1,-1,0, iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB(-1, 1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB(-1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB(-1, 0, 1, iBlock)==-1))&
            ) then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB( 1,-1, 0, iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB( 1, 1, 0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB( 1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB( 1, 0, 1, iBlock)==-1))&
            ) then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    ! y direction
    if(nJ == 1)then
       ! 1D
       jR = j; jL = j
       Ay = 0.0; By = 0.0; Cy = 0.0
    else
       jR = j+1; jL = j-1
       Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
       if(j==1)then
          if(NeiLsouth(iBlock)==-1 &
               .or. (iDir==x_ .and. &
               (i==1    .and. DiLevelNei_IIIB(-1,-1,0, iBlock)==-1) .or. &
               (i==nI+1 .and. DiLevelNei_IIIB( 1,-1,0, iBlock)==-1)) &
               .or. (iDir==z_ .and. &
               (k==1    .and. DiLevelNei_IIIB( 0,-1,-1, iBlock)==-1) .or. &
               (k==nK+1 .and. DiLevelNei_IIIB( 0,-1, 1, iBlock)==-1))&
               )then
             jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
          end if
       elseif(j==nJ)then
          if(NeiLnorth(iBlock)==-1 &
               .or. (iDir==x_ .and. &
               (i==1    .and. DiLevelNei_IIIB(-1, 1,0, iBlock)==-1) .or. &
               (i==nI+1 .and. DiLevelNei_IIIB( 1, 1,0, iBlock)==-1)) &
               .or. (iDir==z_ .and. &
               (k==1    .and. DiLevelNei_IIIB( 0, 1,-1, iBlock)==-1) .or. &
               (k==nK+1 .and. DiLevelNei_IIIB( 0, 1, 1, iBlock)==-1))&
               )then
             jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
          end if
       end if
    end if

    ! z direction
    if(nK == 1)then
       ! 1D or 2D
       kR = k; kL = k
       Az = 0.0; Bz = 0.0; Cz = 0.0
    else
       kR = k+1; kL = k-1
       Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz
       if(k==1)then
          if(NeiLbot(iBlock)==-1 &
               .or. (iDir==x_ .and. &
               (i==1    .and. DiLevelNei_IIIB(-1,0,-1, iBlock)==-1) .or. &
               (i==nI+1 .and. DiLevelNei_IIIB( 1,0,-1, iBlock)==-1)) &
               .or. (iDir==y_ .and. &
               (j==1    .and. DiLevelNei_IIIB( 0,-1,-1, iBlock)==-1) .or. &
               (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1,-1, iBlock)==-1))&
               )then
             kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
          end if
       elseif(k==nK)then
          if(NeiLtop(iBlock)==-1 &
               .or. (iDir==x_ .and. &
               (i==1    .and. DiLevelNei_IIIB(-1,0, 1, iBlock)==-1) .or. &
               (i==nI+1 .and. DiLevelNei_IIIB( 1,0, 1, iBlock)==-1)) &
               .or. (iDir==y_ .and. &
               (j==1    .and. DiLevelNei_IIIB( 0,-1,1, iBlock)==-1) .or. &
               (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1,1, iBlock)==-1))&
               )then
             kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
          end if
       end if
    end if

    if(IsCartesianGrid)then
       call calc_cartesian_j
    else
       call calc_gencoord_j
    end if                             

    if(DoTestMe)then
       write(*,*)NameSub,': iDir,i,j,k,iBlock=',iDir,i,j,k,iBlock
       write(*,*)NameSub,': iL,jL,kL,iR,jR,kR=',iL,jL,kL,iR,jR,kR
       write(*,*)NameSub,': Ax,Bx,Cx=',Ax,Bx,Cx
       write(*,*)NameSub,': Ay,By,Cy=',Ay,By,Cy
       write(*,*)NameSub,': Az,Bz,Cz=',Az,Bz,Cz
       do k1=k-1,k+1; do j1=j-1,j+1; do i1=i-1,i+1
          write(*,*)NameSub,': i,j,k,b_DG(x_:z_)=',&
               i1,j1,k1,b_DG(x_:z_,i1,j1,k1)
       end do; end do; end do
       write(*,*)NameSub,': Jx=',Jx
       write(*,*)NameSub,': Jy=',Jy
       write(*,*)NameSub,': Jz=',Jz
    end if

  contains
    !==========================================================================
    subroutine calc_cartesian_j

      select case(iDir)
      case(x_)
         Jy = -InvDx* (b_DG(z_,i  ,j,k) - b_DG(z_,i-1,j,k)) &
              + Az*(b_DG(x_,i-1,j,kL)+b_DG(x_,i,j ,kL))     &
              + Bz*(b_DG(x_,i-1,j,k )+b_DG(x_,i,j ,k ))     &
              + Cz*(b_DG(x_,i-1,j,kR)+b_DG(x_,i,j ,kR))

         Jz = +InvDx* (b_DG(y_,i,j,k) - b_DG(y_,i-1,j,k)) &
              - Ay*(b_DG(x_,i-1,jL,k)+b_DG(x_,i,jL,k)) &
              - By*(b_DG(x_,i-1,j ,k)+b_DG(x_,i,j ,k)) &
              - Cy*(b_DG(x_,i-1,jR,k)+b_DG(x_,i,jR,k)) 

         Jx = + Ay*(b_DG(z_,i-1,jL,k)+b_DG(z_,i,jL,k )) &
              + By*(b_DG(z_,i-1,j ,k)+b_DG(z_,i,j ,k )) &
              + Cy*(b_DG(z_,i-1,jR,k)+b_DG(z_,i,jR,k )) &
              - Az*(b_DG(y_,i-1,j,kL)+b_DG(y_,i,j ,kL)) &
              - Bz*(b_DG(y_,i-1,j,k )+b_DG(y_,i,j ,k )) &
              - Cz*(b_DG(y_,i-1,j,kR)+b_DG(y_,i,j ,kR))

         ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
         if(IsRzGeometry) Jx = Jx + 0.5*(b_DG(z_,i,j,k)+b_DG(z_,i-1,j,k)) &
              /Xyz_DGB(2,i,j,k,iBlock)

      case(y_)
         Jx = + InvDy*(b_DG(z_,i,j,k) - b_DG(z_,i,j-1,k)) &
              - Az*(b_DG(y_,i,j-1,kL) + b_DG(y_,i,j ,kL)) &
              - Bz*(b_DG(y_,i,j-1,k ) + b_DG(y_,i,j ,k )) &
              - Cz*(b_DG(y_,i,j-1,kR) + b_DG(y_,i,j ,kR))

         Jz = - InvDy*(b_DG(x_,i,j,k) - b_DG(x_,i,j-1,k)) &
              + Ax*(b_DG(y_,iL,j-1,k) + b_DG(y_,iL,j ,k)) &
              + Bx*(b_DG(y_,i ,j-1,k) + b_DG(y_,i ,j ,k)) &
              + Cx*(b_DG(y_,iR,j-1,k) + b_DG(y_,iR,j ,k))

         Jy = + Az*(b_DG(x_,i,j-1,kL) + b_DG(x_,i,j ,kL)) &
              + Bz*(b_DG(x_,i,j-1,k ) + b_DG(x_,i,j ,k )) &
              + Cz*(b_DG(x_,i,j-1,kR) + b_DG(x_,i,j ,kR)) &
              - Ax*(b_DG(z_,iL,j-1,k) + b_DG(z_,iL,j ,k)) &
              - Bx*(b_DG(z_,i ,j-1,k) + b_DG(z_,i ,j ,k)) &
              - Cx*(b_DG(z_,iR,j-1,k) + b_DG(z_,iR,j ,k))

         ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
         if(IsRzGeometry)then
            if(Xyz_DGB(2,i,j-1,k,iBlock)<0.0)then
               ! Just for bookkeeping. It's effect is zeroed by zero face area
               Jx = Jx + b_DG(z_,i,j,k)/Xyz_DGB(2,i,j,k,iBlock)
            else
               Jx = Jx + (b_DG(z_,i,j,k)+b_DG(z_,i,j-1,k)) &
                    /(Xyz_DGB(2,i,j,k,iBlock) + Xyz_DGB(2,i,j-1,k,iBlock))
            end if
         end if

      case(z_)
         Jx = -InvDz*(b_DG(y_,i,j,k) - b_DG(y_,i,j,k-1)) & 
              + Ay*(b_DG(z_,i,jL,k-1) + b_DG(z_,i,jL,k))  &
              + By*(b_DG(z_,i,j ,k-1) + b_DG(z_,i,j ,k))  &
              + Cy*(b_DG(z_,i,jR,k-1) + b_DG(z_,i,jR,k))

         Jy = +InvDz*(b_DG(x_,i,j,k) - b_DG(x_,i,j,k-1)) &
              - Ax*(b_DG(z_,iL,j,k-1) + b_DG(z_,iL,j,k))  &
              - Bx*(b_DG(z_,i ,j,k-1) + b_DG(z_,i ,j,k))  &
              - Cx*(b_DG(z_,iR,j,k-1) + b_DG(z_,iR,j,k))

         Jz = + Ax*(b_DG(y_,iL,j,k-1) + b_DG(y_,iL,j,k))  &
              + Bx*(b_DG(y_,i ,j,k-1) + b_DG(y_,i ,j,k))  &
              + Cx*(b_DG(y_,iR,j,k-1) + b_DG(y_,iR,j,k))  &
              - Ay*(b_DG(x_,i,jL,k-1) + b_DG(x_,i,jL,k))  &
              - By*(b_DG(x_,i,j ,k-1) + b_DG(x_,i,j ,k))  &
              - Cy*(b_DG(x_,i,jR,k-1) + b_DG(x_,i,jR,k)) 

         ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
         if(IsRzGeometry) Jx = Jx + b_DG(z_,i,j,k)/Xyz_DGB(2,i,j,k,iBlock)

      case default
         write(*,*)'Error in get_face_current: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_j

    
    !==========================================================================
    subroutine calc_gencoord_j

      real :: DbDgen_DD(MaxDim,MaxDim)
      !-----------------------------------------------------------------------

      ! Calculate the partial derivatives dB/dCoord
      select case(iDir)
      case(x_)
         DbDgen_DD(:,1) = InvDx*(b_DG(:,i,j,k) - b_DG(:,i-1,j,k))
         DbDgen_DD(:,2) = Ay*(b_DG(:,i-1,jL,k) + b_DG(:,i,jL,k))    &
              +           By*(b_DG(:,i-1,j ,k) + b_DG(:,i,j ,k))    &
              +           Cy*(b_DG(:,i-1,jR,k) + b_DG(:,i,jR,k))
         DbDgen_DD(:,3) = Az*(b_DG(:,i-1,j,kL) + b_DG(:,i,j ,kL))   &
              +           Bz*(b_DG(:,i-1,j,k ) + b_DG(:,i,j ,k ))   &
              +           Cz*(b_DG(:,i-1,j,kR) + b_DG(:,i,j ,kR))
         
      case(y_)
         DbDgen_DD(:,1) = Ax*(b_DG(:,iL,j-1,k) + b_DG(:,iL,j ,k))   &
              +           Bx*(b_DG(:,i ,j-1,k) + b_DG(:,i ,j ,k))   &
              +           Cx*(b_DG(:,iR,j-1,k) + b_DG(:,iR,j ,k))

         DbDgen_DD(:,2) = InvDy*(b_DG(:,i,j,k) - b_DG(:,i,j-1,k))
         DbDgen_DD(:,3) = Az*(b_DG(:,i,j-1,kL) + b_DG(:,i,j ,kL))   &
              +           Bz*(b_DG(:,i,j-1,k ) + b_DG(:,i,j ,k ))   &
              +           Cz*(b_DG(:,i,j-1,kR) + b_DG(:,i,j ,kR))

      case(z_)
         DbDgen_DD(:,1) = Ax*(b_DG(:,iL,j,k-1) + b_DG(:,iL,j,k))    &
              +           Bx*(b_DG(:,i ,j,k-1) + b_DG(:,i ,j,k))    &
              +           Cx*(b_DG(:,iR,j,k-1) + b_DG(:,iR,j,k))
         DbDgen_DD(:,2) = Ay*(b_DG(:,i,jL,k-1) + b_DG(:,i,jL,k))    &
              +           By*(b_DG(:,i,j ,k-1) + b_DG(:,i,j ,k))    &
              +           Cy*(b_DG(:,i,jR,k-1) + b_DG(:,i,jR,k))
         DbDgen_DD(:,3) = InvDz*(b_DG(:,i,j,k) - b_DG(:,i,j,k-1))
      end select

      ! Jx = Dbz/Dy - Dby/Dz
      Jx =   sum(DbDgen_DD(z_,:)*DgenDxyz_DDFD(:,y_,i,j,k,iDir)) &
           - sum(DbDgen_DD(y_,:)*DgenDxyz_DDFD(:,z_,i,j,k,iDir))

      ! Jy = Dbx/Dz - Dbz/Dx
      Jy =   sum(DbDgen_DD(x_,:)*DgenDxyz_DDFD(:,z_,i,j,k,iDir)) &
           - sum(DbDgen_DD(z_,:)*DgenDxyz_DDFD(:,x_,i,j,k,iDir))

      ! Jz = Dby/Dx - Dbx/Dy
      Jz =   sum(DbDgen_DD(y_,:)*DgenDxyz_DDFD(:,x_,i,j,k,iDir)) &
           - sum(DbDgen_DD(x_,:)*DgenDxyz_DDFD(:,y_,i,j,k,iDir))

    end subroutine calc_gencoord_j

  end subroutine get_face_current

  !==========================================================================

  subroutine test_face_current

    use ModMain,     ONLY: nI, nJ, nK, nBlock, Unused_B, x_, y_, z_
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, nVar
    use ModParallel, ONLY: NeiLev
    use ModFaceValue,ONLY: correct_monotone_restrict
    use BATL_lib,    ONLY: message_pass_cell, IsCartesianGrid, IsRLonLat, &
         Xyz_DGB, CellSize_DB

    integer, parameter :: nTest = 2
    integer :: i,j,k,iBlock,iTest
    real :: Jx, Jy, Jz
    !------------------------------------------------------------------------

    write(*,*)'test_face_current starting !!!'

    do iTest = 1, nTest

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          select case(iTest)
          case(1)                          
             State_VGB(Bx_,:,:,:,iBlock) = &
                  + 1*Xyz_DGB(x_,:,:,:,iBlock) &
                  + 2*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 3*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 4*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 5*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 6*Xyz_DGB(y_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock)
             State_VGB(By_,:,:,:,iBlock) = &
                  + 10*Xyz_DGB(x_,:,:,:,iBlock) &
                  + 20*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 30*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 40*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 50*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 60*Xyz_DGB(y_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock)
             State_VGB(Bz_,:,:,:,iBlock) = &
                  + 100*Xyz_DGB(x_,:,:,:,iBlock) &
                  + 200*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 300*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 400*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(y_,:,:,:,iBlock) &
                  + 500*Xyz_DGB(x_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock) &
                  + 600*Xyz_DGB(y_,:,:,:,iBlock)*Xyz_DGB(z_,:,:,:,iBlock)
             if(IsRlonLat)then
                State_VGB(Bx_,:,:,:,iBlock) = &
                     Xyz_DGB(y_,:,:,:,iBlock)
                State_VGB(By_,:,:,:,iBlock) = Xyz_DGB(z_,:,:,:,iBlock)
                State_VGB(Bz_,:,:,:,iBlock) = 0.0
             end if                                    
          case(2)
             State_VGB(Bx_,:,:,:,iBlock) = 1.0 + &
                  0.01*Xyz_DGB(x_,:,:,:,iBlock)**2 + &
                  0.02*Xyz_DGB(y_,:,:,:,iBlock)**2 + &
                  0.03*Xyz_DGB(z_,:,:,:,iBlock)**2
             State_VGB(By_,:,:,:,iBlock) = 10.0 + &
                  0.1*Xyz_DGB(x_,:,:,:,iBlock)**2 + &
                  0.2*Xyz_DGB(y_,:,:,:,iBlock)**2 + &
                  0.3*Xyz_DGB(z_,:,:,:,iBlock)**2
             State_VGB(Bz_,:,:,:,iBlock) = 100.0 + &
                  1.0*Xyz_DGB(x_,:,:,:,iBlock)**2 + &
                  2.0*Xyz_DGB(y_,:,:,:,iBlock)**2 + &
                  3.0*Xyz_DGB(z_,:,:,:,iBlock)**2
          end select

       end do

       call message_pass_cell(nVar,State_VGB)

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          IsNewBlockCurrent = .true.

          call correct_monotone_restrict(iBlock)

          do k=1, nK; do j=1,nJ; do i=1,nI+1
             call get_face_current(x_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('x')
          end do; end do; end do

          do k=1, nK; do j=1,nJ+1; do i=1,nI
             call get_face_current(y_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('y')
          end do; end do; end do

          do k=1, nK+1; do j=1,nJ; do i=1,nI
             call get_face_current(z_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('z')
          end do; end do; end do

       end do

       write(*,*)'test_face_current: test ',iTest,' passed !!!'

    end do

    call timing_report

    call stop_mpi('test_face_current succeeded !!! ')

  contains

    subroutine check_error(NameDir)
      character, intent(in):: NameDir
      real :: x, y, z, JxGood, JyGood, JzGood, Tolerance, r

      !--------------------------------------------------------------------
      ! Face center coordinates
      select case(NameDir)
      case('x')
         if(i==1   .and.neiLEV(1,iBlock)==-1) RETURN
         if(i==nI+1.and.neiLEV(2,iBlock)==-1) RETURN
         x = 0.5*(Xyz_DGB(x_,i-1,j,k,iBlock) + Xyz_DGB(x_,i,j,k,iBlock))
         y = 0.5*(Xyz_DGB(y_,i-1,j,k,iBlock) + Xyz_DGB(y_,i,j,k,iBlock))
         z = 0.5*(Xyz_DGB(z_,i-1,j,k,iBlock) + Xyz_DGB(z_,i,j,k,iBlock))
     case('y')
         if(j==1   .and.neiLEV(3,iBlock)==-1) RETURN
         if(j==nJ+1.and.neiLEV(4,iBlock)==-1) RETURN
         x = 0.5*(Xyz_DGB(x_,i,j-1,k,iBlock) + Xyz_DGB(x_,i,j,k,iBlock))
         y = 0.5*(Xyz_DGB(y_,i,j-1,k,iBlock) + Xyz_DGB(y_,i,j,k,iBlock))
         z = 0.5*(Xyz_DGB(z_,i,j-1,k,iBlock) + Xyz_DGB(z_,i,j,k,iBlock))
      case('z')
         if(k==1   .and.neiLEV(5,iBlock)==-1) RETURN
         if(k==nK+1.and.neiLEV(6,iBlock)==-1) RETURN
         x = 0.5*(Xyz_DGB(x_,i,j,k-1,iBlock) + Xyz_DGB(x_,i,j,k,iBlock))
         y = 0.5*(Xyz_DGB(y_,i,j,k-1,iBlock) + Xyz_DGB(y_,i,j,k,iBlock))
         z = 0.5*(Xyz_DGB(z_,i,j,k-1,iBlock) + Xyz_DGB(z_,i,j,k,iBlock))
         !avoid pole
         if(Xyz_DGB(x_,i,j,k-2,iBlock)*Xyz_DGB(x_,i,j,k,iBlock) < 0.0 ) return
         if(Xyz_DGB(x_,i,j,k,iBlock)*Xyz_DGB(x_,i,j,k+1,iBlock) < 0.0 ) return

      end select

      select case(iTest)
      case(1)
         r = sqrt(x**2+y**2+z**2)
         JxGood = 100.0*y/r - 10.0*z/r -30.0
         JyGood = 3.0+ z/r - 100.0*x/r
         JzGood = 10.0*x/r - y/r

         JxGood = 200 + 400*x + 600*z -  30 -  50*x -  60*y
         JyGood =   3 +   5*x +   6*y - 100 - 400*y - 500*z
         JzGood =  10 +  40*y +  50*z -   2 -   4*x -   6*z

         if(IsCartesianGrid)then                
            Tolerance = 1.e-6
         else
            if(IsRlonLat)then
               ! This is an easier test
               JxGood = -1.0
               JyGood =  0.0
               JzGood = -1.0
            end if
            Tolerance = 5.e-3
         end if                              

      case(2)
         JxGood = 4.0 *y - 0.6 *z
         JyGood = 0.06*z - 2.0 *x
         JzGood = 0.2 *x - 0.04*y

         if(IsCartesianGrid)then                
            Tolerance = 5e-2
         else
            Tolerance = 1.e-1
         end if                              
      end select

      if(       abs(Jx-JxGood) > Tolerance*max(abs(JxGood),1.0)  &
           .or. abs(Jy-JyGood) > Tolerance*max(abs(JyGood),1.0)  &
           .or. abs(Jz-JzGood) > Tolerance*max(abs(JzGood),1.0)) then

         write(*,*)'Face=',NameDir
         write(*,*)'iTest, i,j,k,iBlock=',iTest, i,j,k,iBlock
         write(*,*)'x,y,z (i,j,k)=', Xyz_DGB(:,i,j,k,iBlock)
         write(*,*)'x,y,z (face)=', x,y,z
         write(*,*)'x,y,z (i,j,k-1)=', Xyz_DGB(:,i,j,k-1,iBlock)
         write(*,*)'dx,dy,dz=', CellSize_DB(:,iBlock)
         write(*,*)'Bad  Jx,Jy,Jz =',Jx,Jy,Jz
         write(*,*)'Good Jx,Jy,Jz =',JxGood,JyGood,JzGood
         write(*,*)'NeiLev=',NeiLev(:,iBlock)
         call stop_mpi('Error')
      end if

    end subroutine check_error

  end subroutine test_face_current
  
  !=========================================================================

  real function hall_factor(iDir, iFace, jFace, kFace , iBlock)

    use BATL_lib, ONLY: Xyz_DGB, x_, y_, z_

    integer, intent(in)::iDir, iFace, jFace, kFace, iBlock 

    real :: x,y,z,rSqr,TanSqr,Distance1,Distance2
    real :: HallFactor

    !--------------------------------------------------------------
    select case(iDir)
    case(0)  !for cell center
       x = Xyz_DGB(x_,iFace,jFace,kFace,iBlock)
       y = Xyz_DGB(y_,iFace,jFace,kFace,iBlock)
       z = Xyz_DGB(z_,iFace,jFace,kFace,iBlock)       
    case(1)
       x = 0.5*sum(Xyz_DGB(x_,iFace-1:iFace,jFace,kFace,iBlock))
       y = 0.5*sum(Xyz_DGB(y_,iFace-1:iFace,jFace,kFace,iBlock))
       z = 0.5*sum(Xyz_DGB(z_,iFace-1:iFace,jFace,kFace,iBlock))
    case(2)
       x = 0.5*sum(Xyz_DGB(x_,iFace,jFace-1:jFace,kFace,iBlock))
       y = 0.5*sum(Xyz_DGB(y_,iFace,jFace-1:jFace,kFace,iBlock))
       z = 0.5*sum(Xyz_DGB(z_,iFace,jFace-1:jFace,kFace,iBlock))
    case(3)
       x = 0.5*sum(Xyz_DGB(x_,iFace,jFace,kFace-1:KFace,iBlock))
       y = 0.5*sum(Xyz_DGB(y_,iFace,jFace,kFace-1:KFace,iBlock))
       z = 0.5*sum(Xyz_DGB(z_,iFace,jFace,kFace-1:KFace,iBlock))
    end select

    HallFactor = HallFactorMax

    rSqr = (x**2 + y**2 + z**2)
    if(rSqr < rSqrInner1)then
       hall_factor=0.0
       RETURN
    else if(rSqr < rSqrInner2)then
       HallFactor = HallFactor*(rSqr-rSqrInner1)/(rSqrInner2 - rSqrInner1)
    endif

    if(TanSqr1 > 0.0 .and. abs(z)>0.0)then
       TanSqr = (x**2+y**2)/z**2
       if(TanSqr < TanSqr1)then
          hall_factor=0.0
          RETURN
       else if(TanSqr < TanSqr2)then
          HallFactor = HallFactor*(TanSqr-TanSqr1)/(TanSqr2-TanSqr1)
       end if
    end if

    select case(NameHallRegion)
    case('all')
       ! Do nothing
    case('user')
       ! hall_factor= &
       !   user_hall_factor(x, y, z, iDir, iFace, jFace, kFace, iBlock)
    case('sphere')
       rSqr = (x-x0Hall)**2 + (y-y0Hall)**2 + (z-z0Hall)**2
       if(rSqr > rSqrSphere2)then
          hall_factor=0.0
          RETURN
       else if(rSqr > rSqrSphere1)then
          HallFactor = HallFactor*(rSqrSphere2-rSqr)/(rSqrSphere2-rSqrSphere1)
       end if
    case('box')
       Distance2 = max( &
            abs(x-x0Hall)/xSizeBox2, &
            abs(y-y0Hall)/ySizeBox2, &
            abs(z-z0Hall)/zSizeBox2 )
       Distance1 = max( &
            abs(x-x0Hall)/xSizeBox1, &
            abs(y-y0Hall)/ySizeBox1, &
            abs(z-z0Hall)/zSizeBox1 )
       if(Distance2 > 0.5)then
          hall_factor=0.0
          RETURN
       else if(Distance1 > 0.5)then
          HallFactor = HallFactor*(0.5-Distance2)/(Distance1-Distance2)
       end if
    case default
       call stop_mpi("Unknown value for NameHallRegion="//NameHallRegion)
    end select

    hall_factor = HallFactor

  end function hall_factor

end module ModHallResist
