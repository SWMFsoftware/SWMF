!^CFG COPYRIGHT UM
!==============================================================================
module ModViscosity

  use ModProcMH, ONLY: iProc

  implicit none

  private ! except

  public Viscosity_set_parameters
  public Viscosity_init
  public Viscosity_clean
  public Viscosity_factor

  ! Velosity vector for each block and fluid
  real, public, allocatable, dimension(:,:,:,:,:) :: U_DGI
  ! Visosity tensor for each fluid
  real, public, allocatable, dimension(:,:,:)     :: Visco_DDI

  logical, public :: UseViscosity=.false.
  logical, public :: IsNewBlockViscosity = .true.

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

  save 

contains

  !=========================================================================

  subroutine Viscosity_set_parameters(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand

    character(len=*), parameter :: NameSub = "Viscosity_set_parameters" 
    integer :: i
    !------------------------------------------------------------------------


    select case(NameCommand)
    case('#VISCOSITY')
       call read_var('UseViscosity',  UseViscosity)
       call read_var('VisosityCoeffisient',  ViscoCoeffSi)
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

  end subroutine Viscosity_set_parameters

  !=========================================================================

  subroutine  Viscosity_init()

    use BATL_size, ONLY: MaxBlock, MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
         nDim, maxDim
    use ModMultiFluid, ONLY: nFluid
    use ModPhysics,  ONLY: Si2No_V, UnitX_, UnitT_, UnitJ_
    !------------------------------------------------------------------------


    if(.not.UseViscosity) RETURN

    ViscoCoeff = ViscoCoeffSi*Si2No_V(UnitX_)**2/Si2No_V(UnitT_)

    if(allocated(U_DGI)) deallocate (U_DGI)
    if(allocated(Visco_DDI)) deallocate (Visco_DDI)

    allocate(U_DGI(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nFluid))
    U_DGI = 0.0

    allocate(Visco_DDI(MaxDim,MaxDim,nFluid))
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

  end subroutine  Viscosity_init

  !=========================================================================
  subroutine  Viscosity_clean()

    if(.not.UseViscosity) RETURN

    if(allocated(U_DGI)) deallocate (U_DGI)
    if(allocated(Visco_DDI)) deallocate (Visco_DDI)

    rSqrInner1 = -1.0
    rSqrInner2 = -1.0
    TanSqr1 = -1.0
    TanSqr2 = -1.0

    NameViscoRegion ='all'

    UseViscosity=.false.

  end subroutine  Viscosity_clean

  !=========================================================================

  real function Viscosity_factor(iDir, iFace, jFace, kFace , iBlock)

    use BATL_lib, ONLY: Xyz_DGB, x_, y_, z_
    use ModPhysics,   ONLY: Si2No_V,No2Si_V,UnitX_,UnitT_   
    use ModGeometry, ONLY: rMin_BLK, r_BLK

    integer, intent(in)::iDir, iFace, jFace, kFace, iBlock 

    real :: x,y,z,r,rSqr,TanSqr,Distance1,Distance2
    real :: ViscoFactor
    logical,save :: doOnes=.true.
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
          ViscoFactor = ViscoFactor*(rSqrSphere2-rSqr)/(rSqrSphere2-rSqrSphere1)
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

    Viscosity_factor = ViscoCoeff*ViscoFactor

  end function Viscosity_factor

end module ModViscosity
