!^CFG COPYRIGHT UM
module ModHallResist

  use ModSize, ONLY: nI, nJ, nK, nDim

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
  real, public:: HallFactor = 1.0

  ! Ion mass per charge may depend on space and time for multispecies
  real, public, allocatable:: IonMassPerCharge_G(:,:,:)

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

  ! Jacobian matrix for covariant grid: Dcovariant/Dcartesian
  real, public :: DgenDxyz_DDFD(nDim, nDim,1:nI+1, 1:nJ+1, 1:nK+1, nDim)
  real, public :: DgenDxyz_DDC(nDim, nDim, nI, nJ, nK)

  ! Public methods
  public :: init_hall_resist, get_face_current, hall_factor, test_face_current
  public :: set_ion_mass_per_charge
  public :: set_block_jacobian_cell           

  ! Local variables
  real :: b_DG(3,-1:nI+2,-1:nJ+2,-1:nK+2)

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
    use ModSize,    ONLY: nI, nJ, nK, nDim
    use ModPhysics, ONLY: IonMassPerCharge

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_hall_resist'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    IonMassPerCharge = HallFactor*IonMassPerCharge

    if (DoTestMe) then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>>>>>>>> HALL Resistivity Parameters <<<<<<<<<<'
       write(*,*)
       write(*,*) 'HallFactor       = ', HallFactor
       write(*,*) 'HallCmaxFactor   = ', HallCmaxFactor
       write(*,*) 'IonMassPerCharge = ', IonMassPerCharge
       ! Omega_Bi=B0/IonMassPerCharge'
       write(*,*)
       write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
       write(*,*) ''
    end if

    if(.not.allocated(HallJ_CD)) allocate(&
         HallJ_CD(nI,nJ,nK,nDim), &
         BxPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         ByPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         BzPerN_G(0:nI+1,0:nJ+1,0:nK+1),&
         IonMassPerCharge_G(0:nI+1,0:nJ+1,0:nK+1) )

    IonMassPerCharge_G = IonMassPerCharge

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

    use ModAdvance, ONLY: State_VGB, Rho_, UseIdealEos
    Use ModConst, ONLY: cMu, cElectronCharge
    use ModVarIndexes, ONLY: &
         UseMultiSpecies, SpeciesFirst_, SpeciesLast_, MassSpecies_V, nVar
    use ModMultiFluid, ONLY: UseMultiIon, iRhoIon_I, MassIon_I
    use ModPhysics, ONLY: IonMassPerCharge, No2Si_V, Si2No_V, UnitX_, &
         UnitRho_, UnitB_, UnitT_
    use ModUser, ONLY: user_material_properties

    ! Set IonMassPerCharge_G based on average mass
    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: State_V(nVar)

    real :: Zav, NatomicSi, Coeff
    !-------------------------------------------------------------------------

    ! Multiply IonMassPerCharge_G by average ion mass = rho_total / n_total

    if(.not.UseIdealEos)then
       Coeff = No2Si_V(UnitRho_)/(cMu*cElectronCharge) &
            *Si2No_V(UnitX_)**2*Si2No_V(UnitRho_) &
            /(Si2No_V(UnitB_)*Si2No_V(UnitT_))

       do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               AverageIonChargeOut = Zav, NatomicOut = NatomicSi)

          IonMassPerCharge_G(i,j,k) = Coeff*State_VGB(Rho_,i,j,k,iBlock) &
               /(Zav*NatomicSi)
       end do; end do; end do
       
    elseif(UseMultiSpecies)then
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          IonMassPerCharge_G(i,j,k) = IonMassPerCharge * &
               State_VGB(Rho_,i,j,k,iBlock) / &
               sum(State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,iBlock) &
               /MassSpecies_V)
       end do; end do; end do
    elseif(UseMultiIon)then
       ! Get mass density per total number denisity
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          IonMassPerCharge_G(i,j,k) = IonMassPerCharge * &
               sum(State_VGB(iRhoIon_I, i, j, k, iBlock)) / &
               sum(State_VGB(iRhoIon_I, i, j, k, iBlock) / MassIon_I)
       end do; end do; end do
    end if

  end subroutine set_ion_mass_per_charge

  !===========================================================================

  subroutine set_block_jacobian_face(iBlock)
    use ModMain, ONLY: x_, y_, z_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModCoordTransform, ONLY: inverse_matrix

    integer, intent(in):: iBlock

    ! Dxyz/Dgen matrix for one cell
    real:: DxyzDgen_DD(nDim, nDim)

    ! Transverse gradients
    real:: TransGrad_DDG(nDim, nDim, -1:nI+2, -1:nJ+2, -1:nK+2)

    ! Cell center coordinates for this block
    real:: Xyz_DG(nDim, -1:nI+2, -1:nJ+2, -1:nK+2)

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

    ! Calculate the dCovariant/dCartesian matrix

    Xyz_DG(x_,:,:,:) = x_BLK(:,:,:,iBlock)
    Xyz_DG(y_,:,:,:) = y_BLK(:,:,:,iBlock)
    Xyz_DG(z_,:,:,:) = z_BLK(:,:,:,iBlock)

    do k=-1,nK+2; do j=-1,nJ+2; do i=1,nI
       TransGrad_DDG(:,1,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i+1,j,k) - Xyz_DG(:,i-1,j,k)) &
            + dp2* (Xyz_DG(:,i+2,j,k) - Xyz_DG(:,i-2,j,k)))
    end do; end do; end do
    do k=-1,nK+2; do j=1,nJ; do i=-1,nI+2
       TransGrad_DDG(:,2,i,j,k)=  &
            ( dp1* (Xyz_DG(:,i,j+1,k) - Xyz_DG(:,i,j-1,k)) &
            + dp2* (Xyz_DG(:,i,j+2,k) - Xyz_DG(:,i,j-2,k)))
    end do; end do; end do

    do k=1,nK; do j=-1,nJ+2; do i=-1,nI+2
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
  subroutine set_block_sph_jacobian_face(iBlock)
    use ModMain,     ONLY: x_, y_, z_, r_, Theta_, Phi_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, TypeGeometry

    integer, intent(in):: iBlock

    ! Dxyz/Dgen matrix for one cell
    real:: DgenDxyz_DD(nDim, nDim)

    ! Indexes
    integer:: i, j, k, Di, Dj, Dk, iL, jL, kL, iDir

    ! All coordinates refer to the face
    real :: x, y, z
    real :: r, InvR, InvR2, zInvR2
    real :: Rxy, InvRxy, Rxy2, InvRxy2

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='set_block_sph_jacobian_face'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    DIRECTION: do iDir = 1, nDim
       select case(iDir)
       case(x_)
          Di = 1; Dj = 0; Dk = 0
       case(y_)
          Di = 0; Dj = 1; Dk = 0
       case(z_)
          Di = 0; Dj = 0; Dk = 1
       end select

       ! This element is always zero
       DgenDxyz_DD(Phi_,z_) =  0.0

       THETA: do k = 1, nK+1
          kL = k - Dk

          RADIUS: do i = 1, nI+1
             iL = i - Di

             ! Vertical coordinate depends on R (i) and Theta (k) only
             z = 0.5*(z_BLK(iL,1,kL,iBlock) + z_BLK(i,1,k,iBlock))

             ! The radial distance of the face depends on R and Theta only
             ! so use an arbitrary index for Phi (j=1):
             x = 0.5*(x_BLK(iL,1-Dj,kL,iBlock) + x_BLK(i,1,k,iBlock))
             y = 0.5*(y_BLK(iL,1-Dj,kL,iBlock) + y_BLK(i,1,k,iBlock))

             r = sqrt(x**2 + y**2 + z**2)

             InvR  = 1.0/r
             InvR2 = InvR**2
             zInvR2 = z*InvR2

             Rxy2   = r**2 - z**2
             if(Rxy2 > 0.0)then
                Rxy    = sqrt(Rxy2)
                InvRxy = 1.0/Rxy
             else
                ! On the Z axis InvRxy is not needed
                Rxy    = 0.0
                InvRxy = 1.0
             end if
             InvRxy2 = InvRxy**2

             ! These elements do not depend on Phi:
             DgenDxyz_DD(r_    ,z_) = z*InvR
             DgenDxyz_DD(Theta_,z_) = InvR2*Rxy
             
             PHI: do j = 1, nJ+1
                jL = j - Dj

                x = 0.5*(x_BLK(iL,jL,kL,iBlock) + x_BLK(i,j,k,iBlock))
                y = 0.5*(y_BLK(iL,jL,kL,iBlock) + y_BLK(i,j,k,iBlock))

                ! dRadius/..
                DgenDxyz_DD(r_, x_)     = x*InvR
                DgenDxyz_DD(r_, y_)     = y*InvR

                ! If radial coordinate is logarithmic, multiply by extra 1/r
                if(TypeGeometry=='spherical_lnr') &
                     DgenDxyz_DD(r_, :) = DgenDxyz_DD(r_, :)*InvR
                
                ! dPhi/..
                DgenDxyz_DD(Phi_, x_)   = -y*InvRxy2
                DgenDxyz_DD(Phi_, y_)   =  x*InvRxy2

                ! dTheta/..
                DgenDxyz_DD(Theta_, x_) = -x*InvRxy*zInvR2
                DgenDxyz_DD(Theta_, y_) = -y*InvRxy*zInvR2

                ! Store the matrix
                DgenDxyz_DDFD(:,:,i,j,k,iDir) = DgenDxyz_DD      

             end do PHI
          end do RADIUS
       end do THETA
    end do DIRECTION
 
  end subroutine set_block_sph_jacobian_face

  !============================================================================

  subroutine set_block_jacobian_cell(iBlock)
    use ModMain, ONLY: x_, y_, z_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModCoordTransform, ONLY: inverse_matrix

    integer, intent(in):: iBlock
    real:: InvDx1Half, InvDx2Half, InvDx3Half
    real:: DxyzDgen_DD(nDim, nDim)
    integer:: i,j,k
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='set_block_jacobian'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dCovariant/dCartesian matrix
    
    InvDx1Half = InvDx*0.5
    InvDx2Half = InvDy*0.5
    InvDx3Half = InvDz*0.5

    ! Get the dCartesian/dCovariant matrix with finite differences
    do k=1,nK; do j=1,nJ; do i=1,nI
       DxyzDgen_DD(x_,1) = InvDx1Half &
            *(x_BLK(i+1,j,k,iBlock) - x_BLK(i-1,j,k,iBlock))
       DxyzDgen_DD(y_,1) = InvDx1Half &
            *(y_BLK(i+1,j,k,iBlock) - y_BLK(i-1,j,k,iBlock))
       DxyzDgen_DD(z_,1) = InvDx1Half &
            *(z_BLK(i+1,j,k,iBlock) - z_BLK(i-1,j,k,iBlock))
       DxyzDgen_DD(x_,2) = InvDx2Half &
            *(x_BLK(i,j+1,k,iBlock) - x_BLK(i,j-1,k,iBlock))
       DxyzDgen_DD(y_,2) = InvDx2Half &
            *(y_BLK(i,j+1,k,iBlock) - y_BLK(i,j-1,k,iBlock))
       DxyzDgen_DD(z_,2) = InvDx2Half &
            *(z_BLK(i,j+1,k,iBlock) - z_BLK(i,j-1,k,iBlock))
       DxyzDgen_DD(x_,3) = InvDx3Half &
            *(x_BLK(i,j,k+1,iBlock) - x_BLK(i,j,k-1,iBlock))
       DxyzDgen_DD(y_,3) = InvDx3Half &
            *(y_BLK(i,j,k+1,iBlock) - y_BLK(i,j,k-1,iBlock))
       DxyzDgen_DD(z_,3) = InvDx3Half &
            *(z_BLK(i,j,k+1,iBlock) - z_BLK(i,j,k-1,iBlock))

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
    use ModGeometry,     ONLY: Dx_BLK, Dy_BLK, Dz_BLK, y_BLK
    use ModCovariant,    ONLY: UseCovariant, IsRzGeometry
    use ModParallel,     ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, BlkNeighborLev
    use ModVarIndexes,   ONLY: Bx_, Bz_
    implicit none

    real, parameter :: cTwoThird = 2.0/3.0, cFourThird = 4.0/3.0

    character(len=*), parameter :: NameSub='get_face_current'

    integer, intent(in):: iDir, i, j, k, iBlock
    real, intent(out)  :: Jx, Jy, Jz

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz

    real :: b1_DG(3,-1:nI+2,-1:nJ+2,-1:nK+2)

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

    InvDx = 1.0/dx_Blk(iBlock)
    InvDy = 1.0/dy_Blk(iBlock)
    InvDz = 1.0/dz_Blk(iBlock)

    if( IsNewBlockCurrent ) then
       b_DG = State_VGB(bx_:bz_,:,:,:,iBlock)
       call set_block_field3(iBlock, 3, b1_DG, b_DG)
       if(UseCovariant .and. .not.IsRzGeometry)then
          !call timing_start('set_block_jac')
          !call set_block_jacobian_cell(iBlock) ! Fast but not accurate
          !if(TypeGeometry=='spherical'.or.TypeGeometry=='spherical_lnr') then
          !   call set_block_sph_jacobian_face(iBlock) !optimized for spherical
          !else
             call set_block_jacobian_face(iBlock)      !general
          !end if
          !call timing_stop('set_block_jac')
       end if                                          

       IsNewBlockCurrent = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev(-1,-1,0, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev(-1, 1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev(-1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev(-1, 0, 1, iBlock)==-1))&
            ) then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 1,-1, 0, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 1, 1, 0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 1, 0, 1, iBlock)==-1))&
            ) then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,-1,0, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,-1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0,-1,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0,-1, 1, iBlock)==-1))&
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif(j==nJ)then
       if(NeiLnorth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 1,0, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0, 1,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0, 1, 1, iBlock)==-1))&
            )then
          jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
       end if
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,0,-1, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,0,-1, iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,-1, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,-1, iBlock)==-1))&
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif(k==nK)then
       if(NeiLtop(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,0, 1, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,0, 1, iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,1, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,1, iBlock)==-1))&
            )then
          kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
       end if
    end if

    if(UseCovariant .and. .not.IsRzGeometry)then
       call calc_covariant_j
    else                              
       call calc_cartesian_j
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
              /y_BLK(i,j,k,iBlock)

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
            if(y_BLK(i,j-1,k,iBlock)<0.0)then
               ! Just for bookkeeping. It's effect is zeroed by zero face area
               Jx = Jx + b_DG(z_,i,j,k)/y_BLK(i,j,k,iBlock)
            else
               Jx = Jx + (b_DG(z_,i,j,k)+b_DG(z_,i,j-1,k)) &
                    /(y_BLK(i,j,k,iBlock)+y_BLK(i,j-1,k,iBlock))
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
         if(IsRzGeometry) Jx = Jx + b_DG(z_,i,j,k)/y_BLK(i,j,k,iBlock)

      case default
         write(*,*)'Error in get_face_current: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_j

    
    !==========================================================================
    subroutine calc_covariant_j

      real :: DbDgen_DD(nDim, nDim)
      !-----------------------------------------------------------------------

      ! Calculate the partial derivatives dB/dCovariant
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

    end subroutine calc_covariant_j

  end subroutine get_face_current

  !==========================================================================

  subroutine test_face_current

    use ModMain,     ONLY: nI, nJ, nK, nBlock, UnusedBlk, x_, y_, z_, &
         east_, west_, south_, north_, bot_, top_
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, nVar
    use ModCovariant,ONLY: UseCovariant                
    use ModGeometry, ONLY: TypeGeometry                
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK, dz_BLK
    use ModParallel, ONLY: NeiLev
    use ModFaceValue,ONLY: correct_monotone_restrict

    integer, parameter :: nTest = 2
    integer :: i,j,k,iBlock,iTest
    real :: Jx, Jy, Jz
    !------------------------------------------------------------------------

    write(*,*)'test_face_current starting !!!'

    do iTest = 1, nTest

       do iBlock = 1, nBlock
          if(UnusedBlk(iBlock)) CYCLE

          select case(iTest)
          case(1)                          
             State_VGB(Bx_,:,:,:,iBlock) = &
                  + 1*x_BLK(:,:,:,iBlock) &
                  + 2*y_BLK(:,:,:,iBlock) &
                  + 3*z_BLK(:,:,:,iBlock) &
                  + 4*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 5*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 6*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
             State_VGB(By_,:,:,:,iBlock) = &
                  + 10*x_BLK(:,:,:,iBlock) &
                  + 20*y_BLK(:,:,:,iBlock) &
                  + 30*z_BLK(:,:,:,iBlock) &
                  + 40*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 50*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 60*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
             State_VGB(Bz_,:,:,:,iBlock) = &
                  + 100*x_BLK(:,:,:,iBlock) &
                  + 200*y_BLK(:,:,:,iBlock) &
                  + 300*z_BLK(:,:,:,iBlock) &
                  + 400*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 500*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 600*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
             if(TypeGeometry == 'spherical' .or.&     
                  TypeGeometry == 'spherical_lnr')then
                State_VGB(Bx_,:,:,:,iBlock) = &
                     y_BLK(:,:,:,iBlock)
                State_VGB(By_,:,:,:,iBlock) = z_BLK(:,:,:,iBlock)
                State_VGB(Bz_,:,:,:,iBlock) = 0.0
             end if                                    
          case(2)
             State_VGB(Bx_,:,:,:,iBlock) = 1.0 + &
                  0.01*x_BLK(:,:,:,iBlock)**2 + &
                  0.02*y_BLK(:,:,:,iBlock)**2 + &
                  0.03*z_BLK(:,:,:,iBlock)**2
             State_VGB(By_,:,:,:,iBlock) = 10.0 + &
                  0.1*x_BLK(:,:,:,iBlock)**2 + &
                  0.2*y_BLK(:,:,:,iBlock)**2 + &
                  0.3*z_BLK(:,:,:,iBlock)**2
             State_VGB(Bz_,:,:,:,iBlock) = 100.0 + &
                  1.0*x_BLK(:,:,:,iBlock)**2 + &
                  2.0*y_BLK(:,:,:,iBlock)**2 + &
                  3.0*z_BLK(:,:,:,iBlock)**2
          end select

       end do

       call message_pass_cells8(.false., .false., .true.,nVar,State_VGB)

       do iBlock = 1, nBlock
          if(UnusedBlk(iBlock)) CYCLE

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
         if(i==1   .and.neiLEV(east_,iBlock)==-1) RETURN
         if(i==nI+1.and.neiLEV(west_,iBlock)==-1) RETURN
         x = 0.5*(x_BLK(i-1,j,k,iBlock) + x_BLK(i,j,k,iBlock))
         y = 0.5*(y_BLK(i-1,j,k,iBlock) + y_BLK(i,j,k,iBlock))
         z = 0.5*(z_BLK(i-1,j,k,iBlock) + z_BLK(i,j,k,iBlock))
     case('y')
         if(j==1   .and.neiLEV(south_,iBlock)==-1) RETURN
         if(j==nJ+1.and.neiLEV(north_,iBlock)==-1) RETURN
         x = 0.5*(x_BLK(i,j-1,k,iBlock) + x_BLK(i,j,k,iBlock))
         y = 0.5*(y_BLK(i,j-1,k,iBlock) + y_BLK(i,j,k,iBlock))
         z = 0.5*(z_BLK(i,j-1,k,iBlock) + z_BLK(i,j,k,iBlock))
      case('z')
         if(k==1   .and.neiLEV(bot_,iBlock)==-1) RETURN
         if(k==nK+1.and.neiLEV(top_,iBlock)==-1) RETURN
         x = 0.5*(x_BLK(i,j,k-1,iBlock) + x_BLK(i,j,k,iBlock))
         y = 0.5*(y_BLK(i,j,k-1,iBlock) + y_BLK(i,j,k,iBlock))
         z = 0.5*(z_BLK(i,j,k-1,iBlock) + z_BLK(i,j,k,iBlock))
         !avoid pole
         if(x_BLK(i,j,k-2,iBlock)*x_BLK(i,j,k,iBlock) < 0.0 ) return
         if(x_BLK(i,j,k,iBlock)*x_BLK(i,j,k+1,iBlock) < 0.0 ) return

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

         if(UseCovariant)then                
            if(  TypeGeometry == 'spherical' .or. &
                 TypeGeometry == 'spherical_lnr')then
               ! This is an easier test
               JxGood = -1.0
               JyGood =  0.0
               JzGood = -1.0
            end if
            Tolerance = 5.e-3
         else                               
            Tolerance = 1.e-6
         end if                              

      case(2)
         JxGood = 4.0 *y - 0.6 *z
         JyGood = 0.06*z - 2.0 *x
         JzGood = 0.2 *x - 0.04*y

         if(UseCovariant)then                
            Tolerance = 1.e-1
         else                               
            Tolerance = 5e-2
         end if                              
      end select

      if(       abs(Jx-JxGood) > Tolerance*max(abs(JxGood),1.0)  &
           .or. abs(Jy-JyGood) > Tolerance*max(abs(JyGood),1.0)  &
           .or. abs(Jz-JzGood) > Tolerance*max(abs(JzGood),1.0)) then

         write(*,*)'Face=',NameDir
         write(*,*)'iTest, i,j,k,iBlock=',iTest, i,j,k,iBlock
         write(*,*)'x,y,z (i,j,k)=', &
              x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
         write(*,*)'x,y,z (face)=', &
              x,y,z
         write(*,*)'x,y,z (i,j,k-1)=', &
              x_BLK(i,j,k-1,iBlock),y_BLK(i,j,k-1,iBlock),z_BLK(i,j,k-1,iBlock)
         write(*,*)'dx,dy,dz=',&
              dx_BLK(iBlock), dy_BLK(iBlock), dz_BLK(iBlock)
         write(*,*)'Bad  Jx,Jy,Jz =',Jx,Jy,Jz
         write(*,*)'Good Jx,Jy,Jz =',JxGood,JyGood,JzGood
         write(*,*)'NeiLev=',NeiLev(:,iBlock)
         call stop_mpi('Error')
      end if

    end subroutine check_error

  end subroutine test_face_current
  
  !=========================================================================

  real function hall_factor(iDir, iFace, jFace, kFace , iBlock)

    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK

    integer, intent(in)::iDir, iFace, jFace, kFace, iBlock 

    real :: x,y,z,rSqr,TanSqr,Distance1,Distance2
    real :: HallFactor

    !--------------------------------------------------------------
    select case(iDir)
    case(0)  !for cell center
       x = x_BLK(iFace,jFace,kFace,iBlock)
       y = y_BLK(iFace,jFace,kFace,iBlock)
       z = z_BLK(iFace,jFace,kFace,iBlock)       
    case(1)
       x = 0.5*sum(x_BLK(iFace-1:iFace,jFace,kFace,iBlock))
       y = 0.5*sum(y_BLK(iFace-1:iFace,jFace,kFace,iBlock))
       z = 0.5*sum(z_BLK(iFace-1:iFace,jFace,kFace,iBlock))
    case(2)
       x = 0.5*sum(x_BLK(iFace,jFace-1:jFace,kFace,iBlock))
       y = 0.5*sum(y_BLK(iFace,jFace-1:jFace,kFace,iBlock))
       z = 0.5*sum(z_BLK(iFace,jFace-1:jFace,kFace,iBlock))
    case(3)
       x = 0.5*sum(x_BLK(iFace,jFace,kFace-1:KFace,iBlock))
       y = 0.5*sum(y_BLK(iFace,jFace,kFace-1:KFace,iBlock))
       z = 0.5*sum(z_BLK(iFace,jFace,kFace-1:KFace,iBlock))
    end select

    HallFactor = 1.0

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
