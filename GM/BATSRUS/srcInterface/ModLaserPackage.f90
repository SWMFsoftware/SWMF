!=================================THE LASER PACKAGE STARTS HERE=================!
!Our intent is to make ray_path routine re-usable,
!which is a reason to put the package here.
!===============================================================================!
module ModLaserPulse
  implicit none
  !The laser pulse is assumed to have the linear raise front
  !and linear decay with the constant irradiance within the
  !pulse,  at  t_raise < t < tPulse - tDecay
  real:: tRaise = 1.0e-10 ![s]
  real:: tDecay = 1.0e-10 ![s]
  real:: tPulse = 1.1e-9  ![s]
  real:: IrradianceSi = 3.8e12 ![J/s]
contains
  real function irradiance_t(TimeSi)
    real,intent(in)::TimeSi
    !----------------------
    irradiance_t = IrradianceSi *&
         max(0.0, min(1.0,       &
         TimeSi/tRaise,          &
         (tPulse -TimeSi)/tDecay))
         
  end function irradiance_t
end module ModLaserPulse
!===========================
subroutine read_laser_pulse_param
  use ModReadparam
  use ModLaserPulse
  use ModMain, ONLY: UseLaserPackage
  implicit none
  !-------------
  !Usage
  !#LASERPULSE
  !T                UseLaserPackage
  !3.8e12           IrradianceSI
  !1.1e-9           tPulse
  !1.0e-10          tRaise
  !1.0e-10          tDecay
  !-------------
  call read_var('UseLaserPackage', UseLaserPackage)
  call read_var('IrradianceSI'   , IrradianceSI   )
  call read_var('tPulse'         , tPulse         )
  call read_var('tRaise'         , tRaise         )
  call read_var('tDecay'         , tDecay         )
  
end subroutine read_laser_pulse_param
!===========================
module ModLaserPackage
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  use ModRadioWaveRaytracing, ONLY: ray_path
  implicit none

  real, allocatable, save:: SourceE_CB(:,:,:,:)
  real, dimension(:,:), pointer:: Position_DI
  logical, save:: DoInit = .true.
  integer, save:: nRay = -1
  !------------
contains
  !================================
  subroutine init_laser_package
    use ModSize, ONLY: MaxBlock, nI, nJ, nK
    use ModDensityAndGradient, ONLY: Density_I, GradDensity_DI,DeltaSNew_I, &
         NameVector
    use ModAbsorption, ONLY: AbsorptionCoeff_I
    use ModMain,ONLY:NameThisComp
    use CON_global_vector, ONLY: allocate_vector, associate_with_global_vector
    !--------------------------------------
    if(nRay < 0)call stop_mpi(&
         'get_rays should be called before the first use of laser package')

    !Allocate arrays for rays
    allocate(SourceE_CB(nI,nJ,nK,MaxBlock))

    if(nRay >0)then
       allocate(Density_I(nRay), GradDensity_DI(3,nRay),DeltaSNew_I(nRay))
       allocate(AbsorptionCoeff_I(nRay))
    end if

    NameVector=NameThisComp//'_Rays_DI'
    call allocate_vector(NameVector, 3, nRay)
    call associate_with_global_vector(Position_DI, NameVector)
    
  end subroutine init_laser_package
  !================================
  subroutine get_impl_energy_source
    if(DoInit) call init_laser_package
    SourceE_CB(:,:,:,:) = 0.0
    
  end subroutine get_impl_energy_source
end module ModLaserPackage

subroutine add_laser_energy_deposition
  !This routine should add the laser energy deposition to the "Explicit residual"
  !within the semi-implicit scheme. Effectively, in this way the heat conduction
  !equation, which otherwise is solved within the semi-implicit scheme WITH ZERO
  !right-hand-side, is solved with the right-hand-side including the laser energy
  !deposition
  use ModSize, ONLY: nI, nJ, nK
  use ModImplicit, ONLY: ResExpl_VCB, iTeImpl, nImplBlk, impl2iblk
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  use ModPhysics, ONLY: Si2No_V, UnitEnergydens_, UnitX_, UnitT_
  use ModLaserPackage, ONLY: get_impl_energy_source, SourceE_CB
  use ModMain, ONLY: Time_Simulation
  use ModLaserPulse, ONLY: irradiance_t
  implicit none
  real:: Irradiance
  integer :: iImplBlock, iBlock, i, j, k, iVar

  character(len=*), parameter:: NameSub = 'add_laser_energy_deposition'
  !----------------
  Irradiance = irradiance_t(Time_Simulation) * &
       Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  call get_impl_energy_source

  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)
     do k = 1, nK; do j = 1, nJ; do i = 1, nI
        ResExpl_VCB(iTeImpl,i,j,k,iImplBlock) = &
             ResExpl_VCB(iTeImpl,i,j,k,iImplBlock) + &
             SourceE_CB(i,j,k,iBlock) * Irradiance
     end do; end do; end do
  end do
  
end subroutine add_laser_energy_deposition
!==========================================
