!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module ModViscosity

  use ModProcMH, ONLY: iProc

  implicit none

  private ! except

  public viscosity_read_param
  public viscosity_init
  public viscosity_clean
  public viscosity_factor
  public get_viscosity_tensor

  ! Use viscosity or not
  logical, public :: UseViscosity=.false.

  ! Visosity tensor for each fluid
  real, public, allocatable:: Visco_DDI(:,:,:)

  ! Logical needed to fil up ghost cells only once per  block
  logical, public :: IsNewBlockViscosity = .true.

  ! Local variables

  ! Velosity vector for each block and fluid
  real, allocatable:: u_DGI(:,:,:,:,:)

  ! Gradient of velocity centered for faces
  real, allocatable:: GradU_DDI(:,:,:)

  real,save :: ViscoCoeff = 0.0, ViscoCoeffSi = 0.0

  ! Name/shape of the region where Visco effect is used
  character(len=20):: NameViscoRegion ='all'

  ! Parameters for exclusion of inner boundary (if present)
  real :: rInnerVisco=-1.0, DrInnerVisco=0.0

  ! Parameters for exclusion of pole region in spherical grids
  real :: PoleAngleVisco=-1.0, dPoleAngleVisco=-1.0

  ! Parameters for the center location of Visco region
  real :: X0Visco=0.0, Y0Visco=0.0, Z0Visco=0.0

  ! Parameters for spherical Visco region
  real :: rSphereVisco=-1.0, DrSphereVisco=-1.0

  ! Parameters for box Visco region
  real :: xSizeBoxVisco=-1.0, ySizeBoxVisco=-1.0, zSizeBoxVisco=-1.0
  real :: DxSizeBoxVisco=-1.0, DySizeBoxVisco=-1.0, DzSizeBoxVisco=-1.0

  ! Local variables for Visco regions
  real :: TanSqr1 = -1.0, TanSqr2 = -1.0                  ! pole
  real :: rSqrInner1 = -1.0, rSqrInner2 = -1.0            ! inner body
  real :: rSqrSphere1 = -1.0, rSqrSphere2 = -1.0          ! spherical region
  real :: xSizeBox1=-1.0, ySizeBox1=-1.0, zSizeBox1=-1.0  ! box region
  real :: xSizeBox2=-1.0, ySizeBox2=-1.0, zSizeBox2=-1.0  ! box region

  real:: rZeroVisco = -1.0
  real:: rFullVisco = 1e30

  ! Artificial viscosity.
  logical, public :: UseArtificialVisco = .false.
  real, public :: alphaVisco 

  save 

contains

  !=========================================================================

  subroutine viscosity_read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter :: NameSub = "viscosity_read_param" 
    integer :: i
    !------------------------------------------------------------------------


    select case(NameCommand)
    case('#VISCOSITY')
       call read_var('UseViscosity',      UseViscosity)
       if(UseViscosity) &
            call read_var('ViscosityCoeffSi',  ViscoCoeffSi)
    case('#ARTIFICIALVISCOSITY')
       call read_var('UseArtificialVisco',      UseArtificialVisco)
       if(UseArtificialVisco) then 
          call read_var('alphaVisco', alphaVisco)
       endif
    case('#VISCOSITYREGION')
       call read_var('NameViscosityRegion', NameViscoRegion)

       i = index(NameViscoRegion, '0')
       if(i < 1 .and. NameViscoRegion /= 'all' .and. &
            NameViscoRegion /= 'mask')then
          call read_var("x0Visco", x0Visco)
          call read_var("y0Visco", y0Visco)
          call read_var("z0Visco", z0Visco)
       else
          x0Visco = 0.0; y0Visco = 0.0; z0Visco = 0.0
          if(i>1)NameViscoRegion = &
               NameViscoRegion(1:i-1)//NameViscoRegion(i+1:len(NameViscoRegion))
       end if

       select case(NameViscoRegion)
       case("all")
       case("sphere")
          call read_var("rSphereVisco",rSphereVisco)
          call read_var("DrSphereVisco",DrSphereVisco)
       case("box")
          call read_var("xSizeBoxVisco ",xSizeBoxVisco)
          call read_var("DxSizeBoxVisco",DxSizeBoxVisco)
          call read_var("ySizeBoxVisco ",ySizeBoxVisco)
          call read_var("DySizeBoxVisco",DySizeBoxVisco)
          call read_var("zSizeBoxVisco ",zSizeBoxVisco)
          call read_var("DzSizeBoxVisco",DzSizeBoxVisco)
       case('mask')
          call read_var('rZeroVisco', rZeroVisco)
          call read_var('rFullVisco', rFullVisco)
       case default
          call stop_mpi(NameSub//': unknown NameViscoRegion='&
               //NameViscoRegion)
       end select
    end select

  end subroutine viscosity_read_param

  !=========================================================================

  subroutine  viscosity_init

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nDim
    use ModMultiFluid, ONLY: nFluid
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitT_
    !------------------------------------------------------------------------


    if(.not.UseViscosity) RETURN

    ViscoCoeff = ViscoCoeffSi*Si2No_V(UnitX_)**2/Si2No_V(UnitT_)

    if(allocated(u_DGI))     deallocate(u_DGI)
    if(allocated(GradU_DDI)) deallocate(GradU_DDI)
    if(allocated(Visco_DDI)) deallocate(Visco_DDI)

    allocate(u_DGI(nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nFluid))
    u_DGI = 0.0

    allocate(GradU_DDI(nDim,nDim,nFluid))
    GradU_DDI = 0.0

    allocate(Visco_DDI(nDim,nDim,nFluid))
    Visco_DDI = 0.0

    rSqrInner1 = -1.0
    rSqrInner2 = -1.0
    TanSqr1 = -1.0
    TanSqr2 = -1.0

    if(PoleAngleVisco > 0.0)then
       TanSqr1 = tan(PoleAngleVisco)**2
       TanSqr2 = tan(PoleAngleVisco+dPoleAngleVisco)**2
    end if

    if(rInnerVisco > 0.0)then
       rSqrInner1 = rInnerVisco**2
       rSqrInner2 = (rInnerVisco + DrInnerVisco)**2
    endif

    select case(NameViscoRegion)
    case('sphere')
       rSqrSphere1 = rSphereVisco**2
       rSqrSphere2 = (rSphereVisco+DrSphereVisco)**2
    case('box')
       xSizeBox1 = xSizeBoxVisco
       ySizeBox1 = ySizeBoxVisco
       zSizeBox1 = zSizeBoxVisco
       xSizeBox2 = xSizeBoxVisco + 2*DxSizeBoxVisco
       ySizeBox2 = ySizeBoxVisco + 2*DySizeBoxVisco
       zSizeBox2 = zSizeBoxVisco + 2*DzSizeBoxVisco
    case('all','mask')
       ! Do nothing
    case default
       call stop_mpi('ERROR in init_Visco_resist: NameViscoRegion=' &
            //NameViscoRegion)
    end select

  end subroutine viscosity_init

  !=========================================================================
  subroutine  viscosity_clean

    if(.not.UseViscosity) RETURN

    if(allocated(u_DGI))     deallocate(u_DGI)
    if(allocated(GradU_DDI)) deallocate(GradU_DDI)
    if(allocated(Visco_DDI)) deallocate(Visco_DDI)

    rSqrInner1 = -1.0
    rSqrInner2 = -1.0
    TanSqr1 = -1.0
    TanSqr2 = -1.0

    NameViscoRegion ='all'

    UseViscosity=.false.

  end subroutine viscosity_clean

  !=========================================================================

  real function viscosity_factor(iDir, iFace, jFace, kFace , iBlock)

    use BATL_lib,    ONLY: Xyz_DGB, x_, y_, z_
    use ModGeometry, ONLY: r_BLK

    integer, intent(in)::iDir, iFace, jFace, kFace, iBlock 

    real :: x,y,z,r,rSqr,Distance1,Distance2
    real :: ViscoFactor
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

    ViscoFactor = 1.0

    rSqr = (x**2 + y**2 + z**2)

    select case(NameViscoRegion)
    case('all')
       ! Do nothing
       ViscoFactor = 1.0
    case('sphere')
       rSqr = (x-x0Visco)**2 + (y-y0Visco)**2 + (z-z0Visco)**2
       if(rSqr > rSqrSphere2)then
          ViscoFactor=0.0
       else if(rSqr > rSqrSphere1)then
          ViscoFactor = ViscoFactor*(rSqrSphere2-rSqr) &
               / (rSqrSphere2 - rSqrSphere1)
       end if
    case('box')
       Distance2 = max( &
            abs(x-x0Visco)/xSizeBox2, &
            abs(y-y0Visco)/ySizeBox2, &
            abs(z-z0Visco)/zSizeBox2 )
       Distance1 = max( &
            abs(x-x0Visco)/xSizeBox1, &
            abs(y-y0Visco)/ySizeBox1, &
            abs(z-z0Visco)/zSizeBox1 )
       if(Distance2 > 0.5)then
          ViscoFactor=0.0
       else if(Distance1 > 0.5)then
          ViscoFactor = ViscoFactor*(0.5-Distance2)/(Distance1-Distance2)
       end if
    case('mask')
       ! Check if the block is fully outside the masked region
       if(rSqr <= rFullVisco**2) then

          r = r_BLK(iFace,jFace,kFace,iBlock)
          if(r < rZeroVisco)then
             ViscoFactor=0.0
          else if(r < rFullVisco)then
             ViscoFactor =  &
                  (r - rZeroVisco)/(rFullVisco - rZeroVisco)
          end if
       end if
       !    ViscoFactor = 0.0
    case default
       call stop_mpi("Unknown value for NameViscoRegion="//NameViscoRegion)
    end select

    viscosity_factor = ViscoCoeff*ViscoFactor

  end function viscosity_factor

  !=========================================================================

  subroutine get_viscosity_tensor(iDimFace, iFace, jFace, &
       kFace,iBlockFace,iFluidMin,iFluidMax,ViscoCoeff)   

    use ModAdvance, ONLY: State_VGB
    use BATL_size, ONLY: nDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use BATL_lib,  ONLY: x_, y_, z_
    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iRho, iRhoUx
    use ModFaceGradient, ONLY: get_face_gradient_field

    integer, intent(in) :: iDimFace, iFace, jFace,kFace,iBlockFace
    integer, intent(in) :: iFluidMin,iFluidMax
    real,    intent(in) :: ViscoCoeff

    real    :: Diag
    logical :: IsNewBlock = .true.
    real, parameter :: TraceCoeff = 2.0/3.0
    integer :: i,j,k
    !------------------------------------------------------------------------    

    ! Get velocity vector for the block, only done ones per block
    if(IsNewBlockViscosity) then
       do iFluid = iFluidMin, iFluidMax
          call select_fluid
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(State_VGB(iRho,i,j,k,iBlockFace) > 0.0) then 
                u_DGI(:,i,j,k,iFluid) = &
                     State_VGB(iRhoUx:iRhoUx+nDim-1,i,j,k,iBlockFace) / &
                     State_VGB(iRho,i,j,k,iBlockFace)
             else
                u_DGI(:,i,j,k,iFluid) = 0.0
             end if
          end do; end do; end do
       end do
    end if

    ! Get the velocity gradient on the faces. 
    ! Fill in ghost cells for the block for each fluid once.
    IsNewBlock = IsNewBlockViscosity
    do iFluid = iFluidMin, iFluidMax
       call get_face_gradient_field(iDimFace, iFace, jFace, kFace, &
            iBlockFace, nDim,  &
            IsNewBlock, u_DGI(:,:,:,:,iFluid), GradU_DDI(:,:,iFluid))
       ! so ghost cell for all fluids are updated
       IsNewBlock = IsNewBlockViscosity 
    end do
    IsNewBlockViscosity = .false.        

    ! Get the viscosity tensor
    do iFluid = 1, nFluid
       Diag              =        GradU_DDI(x_,x_,iFluid)
       if(nDim > 1) Diag = Diag + GradU_DDI(y_,y_,iFluid)
       if(nDim > 2) Diag = Diag + GradU_DDI(z_,z_,iFluid)
       Diag = -TraceCoeff*Diag

       ! Diagonal
       Visco_DDI(x_,x_,iFluid)             = 2*GradU_DDI(x_,x_,iFluid) + Diag
       if(nDim > 1)Visco_DDI(y_,y_,iFluid) = 2*GradU_DDI(y_,y_,iFluid) + Diag
       if(nDim > 2)Visco_DDI(z_,z_,iFluid) = 2*GradU_DDI(z_,z_,iFluid) + Diag

       ! Off-diagonal terms are symmetrized
       if(nDim > 1)then
          Visco_DDI(x_,y_,iFluid) = &
               GradU_DDI(x_,y_,iFluid) + GradU_DDI(y_,x_,iFluid)
          Visco_DDI(y_,x_,iFluid) = Visco_DDI(x_,y_,iFluid)
       end if
       if(nDim > 2)then
          Visco_DDI(x_,z_,iFluid) = &
               GradU_DDI(x_,z_,iFluid) + GradU_DDI(z_,x_,iFluid)
          Visco_DDI(z_,x_,iFluid) = Visco_DDI(x_,z_,iFluid)
          Visco_DDI(y_,z_,iFluid) = &
               GradU_DDI(y_,z_,iFluid) + GradU_DDI(z_,y_,iFluid)
          Visco_DDI(z_,y_,iFluid) = Visco_DDI(y_,z_,iFluid)
       end if

       Visco_DDI(:,:,iFluid) = ViscoCoeff*Visco_DDI(:,:,iFluid)          
    end do

  end subroutine get_viscosity_tensor

end module ModViscosity
