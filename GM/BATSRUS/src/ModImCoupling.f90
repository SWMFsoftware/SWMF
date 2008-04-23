!^CFG COPYRIGHT UM
!^CFG FILE RCM
!==========================================================================
module ModImPressure

  use ModUtilities, ONLY: check_allocate

  save

  ! The number of IM pressures obtained so far
  integer :: iNewPIm = 0

  ! The size of the IM grid
  integer :: iSize, jSize

  real, dimension(:), allocatable   :: RCM_lat, RCM_lon
  real, dimension(:,:), allocatable :: RCM_p, RCM_dens

contains

  subroutine im_pressure_init(iSizeIn,jSizeIn)
    integer :: iSizeIn, jSizeIn, iError

    iSize = iSizeIn
    jSize = jSizeIn
    allocate(&
         RCM_lat(iSize), &
         RCM_lon(jSize), &
         RCM_p(iSize,jSize), &
         RCM_dens(iSize,jSize), &
         stat=iError)

    call check_allocate(iError, 'GM_im_pressure_init: RCM arrays')

  end subroutine im_pressure_init

end module ModImPressure

!==========================================================================

subroutine get_im_pressure(iBlock, pIm_C, dIm_C)

  use ModImPressure
  use ModMain,     ONLY : nI, nJ, nK, DoFixPolarRegion, rFixPolarRegion
  use ModRaytrace, ONLY : ray
  use ModPhysics,  ONLY : Si2No_V, UnitP_, UnitRho_, PolarRho_I, PolarP_I
  use ModGeometry, ONLY : R_BLK, z_BLK
  use ModAdvance,  ONLY : State_VGB, RhoUz_
  implicit none

  integer, intent(in)  :: iBlock
  real,    intent(out) :: pIm_C(1:nI, 1:nJ, 1:nK)
  real,    intent(out) :: dIm_C(1:nI, 1:nJ, 1:nK)

  integer :: i,j,k, n, i1,i2

  real :: rLat,rLon
  !--------------------------------------------------------------------------
  !\
  ! Check to see if cell centers are on closed fieldline
  !/
  do k=1,nK; do j=1,nJ; do i=1,nI
     if(nint(ray(3,1,i,j,k,iBlock)) == 3) then
        rLat=ray(1,1,i,j,k,iBlock)
        rLon=ray(2,1,i,j,k,iBlock)

        !NOTE: RCM_lat in decending order
        i1=1
        do n=2,isize
           if(rLat < 0.5*(RCM_lat(n-1)+RCM_lat(n))) i1 = n
        end do

        !NOTE: RCM_lon in accending order
        i2=1
        do n=2,jsize
           if(rLon > 0.5*(RCM_lon(n-1)+RCM_lon(n))) i2 = n
        end do
        if(rLon > 0.5*(RCM_lon(jsize)+RCM_lon(1)+360.)) i2=1

        pIm_C(i,j,k) = RCM_p(i1,i2)*Si2No_V(UnitP_)
        dIm_C(i,j,k) = RCM_dens(i1,i2)*Si2No_V(UnitRho_)
     else if(DoFixPolarRegion .and. R_BLK(i,j,k,iBlock) < rFixPolarRegion &
        .and. z_BLK(i,j,k,iBlock)*State_VGB(RhoUz_,i,j,k,iBlock) > 0.0)then
        pIm_C(i,j,k) = PolarP_I(1)
        dIm_C(i,j,k) = PolarRho_I(1)
     else
        pIm_C(i,j,k) = -1.0
        dIm_C(i,j,k) = -1.0
     end if
  end do; end do; end do

end subroutine get_im_pressure

!==========================================================================

subroutine apply_im_pressure

  use ModMain,    ONLY: nI, nJ, nK, nBlock, iNewGrid, TauCoupleIm, &
       time_accurate, Dt, DoCoupleImPressure,DoCoupleImDensity, unusedBLK
  use ModAdvance, ONLY: State_VGB, Energy_GBI, &
       Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, p_
  use ModPhysics, ONLY: Si2No_V, UnitT_, inv_gm1
  use ModImPressure

  implicit none

  real :: Factor
  real :: pIm_C(1:nI, 1:nJ, 1:nK)
  real :: dIm_C(1:nI, 1:nJ, 1:nK)

  integer :: iLastPIm = -1, iLastGrid = -1

  integer :: iBlock
  character (len=*), parameter :: NameSub='apply_im_pressure'
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  if(iNewPIm < 1) RETURN ! No IM pressure has been obtained yet
  if(.not.DoCoupleImPressure .and. .not.DoCoupleImDensity) RETURN  ! Nothing to do

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! redo ray tracing if necessary 
  ! (load_balance takes care of new decomposition)
  if(iNewPIm > iLastPIm .or. iNewGrid > iLastGrid) then
     if(DoTestMe)write(*,*)'GM_apply_im_pressure: call ray_trace ',&
          'iNewPIm,iLastPIm,iNewGrid,iLastGrid=',&
          iNewPIm,iLastPIm,iNewGrid,iLastGrid
     call ray_trace
  end if

  ! Remember this call
  iLastPIm = iNewPIm; iLastGrid = iNewGrid

  ! Now use the pressure from the IM to nudge the pressure in the MHD code.  
  ! This will happen only on closed magnetic field lines.
  ! Determining which field lines are closed is done by using the ray
  ! tracing.

  if(time_accurate)then
     ! Ramp up is based on physical time: p' = p + dt/tau * (pIM - p)
     ! A typical value might be 5, to get close to the RCM pressure 
     ! in 10 seconds

     Factor = 1.0/(TauCoupleIM*Si2No_V(UnitT_))

  else
     ! Ramp up is based on number of iterations: p' = (ntau*p + pIm)/(1+ntau)
     ! A typical value might be 20, to get close to the RCM pressure 
     ! in 20 iterations
     
     Factor = 1.0/(1.0+TauCoupleIM)

  end if

  do iBlock = 1, nBlock
     if(unusedBLK(iBlock)) CYCLE

     call get_im_pressure(iBlock, pIm_C, dIm_C)

     if(all(pIm_C < 0.0)) CYCLE  ! Nothing to do

     !Put velocity into momentum temporarily when density is changed
     if(DoCoupleImDensity)then
        where(dIm_C > 0.0)
           State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBlock)/ &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
           State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBlock)/ &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
           State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBlock)/ &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
        end where
     end if

     if(time_accurate)then
        if(DoCoupleImPressure)then
           where(pIm_C > 0.0) &
                State_VGB(P_,1:nI,1:nJ,1:nK,iBlock) = &
                State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)   &
                + min(1.0, Factor * Dt) &
                * (pIm_C - State_VGB(P_,1:nI,1:nJ,1:nK,iBlock))
        end if
        if(DoCoupleImDensity)then
           where(dIm_C > 0.0) &
                State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock) = &
                State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)   &
                + min(1.0, Factor * Dt) &
                * (dIm_C - State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock))
        end if
     else
        if(DoCoupleImPressure)then
           where(pIm_C > 0.0) &
                State_VGB(P_,1:nI,1:nJ,1:nK,iBlock) = Factor* &
                (TauCoupleIM*State_VGB(P_,1:nI,1:nJ,1:nK,iBlock) + pIm_C)
        end if
        if(DoCoupleImDensity)then
           where(dIm_C > 0.0) &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock) = Factor* &
                (TauCoupleIM*State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock) + dIm_C)
        end if
     end if

     !Convert back to momentum
     if(DoCoupleImDensity)then
        where(dIm_C > 0.0)
           State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBlock)* &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
           State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBlock)* &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
           State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBlock)= &
                State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBlock)* &
                State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)
        end where
     end if

     ! Now get the energy that corresponds to these new values
     Energy_GBI(1:nI,1:nJ,1:nK,iBlock,1) = &
          inv_gm1*State_VGB(P_,1:nI,1:nJ,1:nK,iBlock) + 0.5*( &
          ( State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBlock)**2 &
          + State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBlock)**2 &
          + State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBlock)**2 &
          )/State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock)  &
          + State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)**2 &
          + State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)**2 &
          + State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)**2)
  end do

end subroutine apply_im_pressure
