!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHallResist

  use ModSize, ONLY: nI, nJ, nK, MaxDim, j0_, nJp1_, k0_, nKp1_

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
  real, public, allocatable :: HallJ_CD(:,:,:,:)

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

  ! Public methods
  public :: init_hall_resist, hall_factor
  public :: set_ion_mass_per_charge
  public :: set_ion_mass_per_charge_point

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

    use ModConst,   ONLY: cElectronCharge
    use ModPhysics, ONLY: IonMassPerCharge, Si2No_V, UnitX_, UnitCharge_

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

    if(.not.allocated(HallJ_CD)) allocate(              &
         HallJ_CD(nI,nJ,nK,MaxDim),                     &
         IonMassPerCharge_G(0:nI+1,j0_:nJp1_,k0_:nKp1_) )

    HallJ_CD = 0.0

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

    use ModAdvance, ONLY: State_VGB, UseIdealEos
    use ModVarIndexes, ONLY: UseMultiSpecies
    use ModMultiFluid, ONLY: UseMultiIon

    ! Set IonMassPerCharge_G based on average mass
    integer, intent(in) :: iBlock

    integer :: i, j, k
    !-------------------------------------------------------------------------

    ! Check if IonMassPerCharge_G varies at all. Return if it is constant.
    if(.not.UseMultiIon .and. .not.UseMultiSpecies .and. UseIdealEos) RETURN

    ! Set IonMassPerCharge_G to the average ion mass = rho_total / n_total
    ! including 1 layer of ghost cells
    do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
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

    real, intent(in) :: State_V(nVar)
    real, intent(out):: IonMassPerChargeOut

    real :: zAverage, NatomicSi
    !--------------------------------------------------------------------------

    if(.not.UseIdealEos)then
       call user_material_z_n(State_V, zAverage, NatomicSi)

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
