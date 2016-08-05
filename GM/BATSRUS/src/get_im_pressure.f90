!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan

!==========================================================================
module ModImPressure

  use ModUtilities, ONLY: check_allocate
  use ModMain, ONLY: DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling

  save

  ! The number of IM pressures obtained so far
  integer :: iNewPIm = 0

  ! The size of the IM grid
  integer :: iSize, jSize

  real, dimension(:), allocatable   :: IM_lat, IM_lon
  real, dimension(:,:), allocatable :: IM_p, IM_dens, IM_ppar, IM_bmin, &
       IM_Hpp, IM_Opp, IM_Hpdens, IM_Opdens

contains

  subroutine im_pressure_init(iSizeIn,jSizeIn)
    integer :: iSizeIn, jSizeIn, iError

    iSize = iSizeIn
    jSize = jSizeIn
    allocate(&
         IM_lat(iSize), &
         IM_lon(jSize), &
         IM_p(iSize,jSize), &
         IM_dens(iSize,jSize), &
         stat=iError)

    call check_allocate(iError, 'GM_im_pressure_init: IM arrays')

    if(DoMultiFluidIMCoupling)then
        allocate(&
             IM_Hpp(iSize,jSize), &
             IM_Opp(iSize,jSize), &
             IM_Hpdens(iSize,jSize), &
             IM_Opdens(iSize,jSize), &
             stat = ierror)
        call check_allocate(ierror, 'GM_im_pressure_init: IM arrays')
     end if

     if(DoAnisoPressureIMCoupling)then
        allocate(&
             IM_ppar(iSize,jSize), &
             IM_bmin(iSize,jSize), &
             stat = ierror)
        call check_allocate(ierror, 'GM_im_pressure_init: IM arrays')
     end if

  end subroutine im_pressure_init

end module ModImPressure

!==========================================================================

subroutine get_im_pressure(iBlock, pIm_IC, RhoIm_IC, TauCoeffIm_C, PparIm_C)

  use ModImPressure
  use ModMain,     ONLY : nI, nJ, nK, DoFixPolarRegion, rFixPolarRegion, &
       dLatSmoothIm, UseB0
  use ModRaytrace, ONLY : ray
  use ModPhysics,  ONLY : Si2No_V, UnitB_, UnitP_, UnitRho_, PolarRho_I, PolarP_I
  use ModGeometry, ONLY : R_BLK, Xyz_DGB, z_
  use ModAdvance,  ONLY : State_VGB, RhoUz_, Bx_, Bz_
  use ModB0,       ONLY: B0_DGB
  use ModMultiFluid, ONLY : IonFirst_, IonLast_, iFluid
  implicit none

  integer, intent(in)  :: iBlock

  real,    intent(out) :: pIm_IC(3,1:nI, 1:nJ, 1:nK)
  real,    intent(out) :: RhoIm_IC(3,1:nI, 1:nJ, 1:nK)
  real,    intent(out) :: TauCoeffIm_C(1:nI, 1:nJ, 1:nK)
  real,    intent(out) :: PparIm_C(1:nI, 1:nJ, 1:nK)

  real    :: BminIm_C(1:nI, 1:nJ, 1:nK), b_D(3)

  integer :: i,j,k, n, iLat1,iLat2, iLon1,iLon2

  real :: Lat,Lon, LatWeight1,LatWeight2, LonWeight1,LonWeight2
  real :: LatMaxIm
  ! variables for anisotropic pressure coupling
  real :: Pperp, PperpInvPpar, Coeff

  integer :: iIonSecond, nIons
  !--------------------------------------------------------------------------
  iIonSecond = min(IonFirst_+1, IonLast_)
  if (DoMultiFluidIMCoupling) then
     nIons = iIonSecond
  else
     nIons = 1
  end if

  TauCoeffIm_C = 1.0

  ! Maximum latitude (ascending or descending) of the IM grid 
  LatMaxIm = max(IM_lat(1), IM_lat(iSize))

  ! Check to see if cell centers are on closed fieldline
  do k=1,nK; do j=1,nJ; do i=1,nI

     ! Default is negative, which means that do not nudge GM values
     pIm_IC(:,i,j,k) = -1.0
     RhoIm_IC(:,i,j,k) = -1.0
     PparIm_C(i,j,k) = -1.0
     BminIm_C(i,j,k) = -1.0

     ! For closed field lines nudge towards IM pressure/density
     if(nint(ray(3,1,i,j,k,iBlock)) == 3) then

        ! Map the point down to the IM grid 
        ! Note: ray values are in SM coordinates!
        Lat = ray(1,1,i,j,k,iBlock)

        ! Do not modify pressure along field lines outside the IM grid
        if(Lat > LatMaxIm) CYCLE

        Lon = ray(2,1,i,j,k,iBlock)

        if (IM_lat(1) > IM_lat(2)) then
           ! IM_lat is in descending order
           do iLat1 = 2, iSize
              if(Lat > IM_lat(iLat1)) EXIT
           end do
           iLat2 = iLat1-1
           LatWeight1 = (Lat - IM_lat(iLat2))/(IM_lat(iLat1) - IM_lat(iLat2))
           LatWeight2 = 1 - LatWeight1
        else
           ! IM lat is in ascending order
           do iLat1 = 2, iSize
              if(Lat < IM_lat(iLat1)) EXIT
           end do
           iLat2 = iLat1-1
           LatWeight1 = &
                (Lat - IM_lat(iLat2))/(IM_lat(iLat1) - IM_lat(iLat2))
           LatWeight2 = 1 - LatWeight1              
        endif

        ! Note: IM_lon is in ascending order
        if(Lon < IM_lon(1)) then
           ! periodic before 1
           iLon1 = 1
           iLon2 = jSize
           LonWeight1 =     (Lon           + 360 - IM_lon(iLon2)) &
                /           (IM_lon(iLon1) + 360 - IM_lon(iLon2))
        elseif(Lon > IM_lon(jSize)) then
           ! periodic after jSize
           iLon1 = 1
           iLon2 = jSize
           LonWeight1 = (Lon                 - IM_lon(iLon2)) &
                /       (IM_lon(iLon1) + 360 - IM_lon(iLon2))
        else
           do iLon1 = 2, jSize
              if(Lon < IM_lon(iLon1)) EXIT
           end do
           iLon2 = iLon1-1
           LonWeight1 = (Lon           - IM_lon(iLon2)) &
                /       (IM_lon(iLon1) - IM_lon(iLon2))
        end if
        LonWeight2 = 1 - LonWeight1

        if(all( IM_p( (/iLat1,iLat2/), (/iLon1, iLon2/) ) > 0.0 ))then
           pIm_IC(1,i,j,k) = Si2No_V(UnitP_)*( &
                LonWeight1 * ( LatWeight1*IM_p(iLat1,iLon1) &
                +              LatWeight2*IM_p(iLat2,iLon1) ) + &
                LonWeight2 * ( LatWeight1*IM_p(iLat1,iLon2) &
                +              LatWeight2*IM_p(iLat2,iLon2) ) )
           RhoIm_IC(1,i,j,k) = Si2No_V(UnitRho_)*( &
                LonWeight1 * ( LatWeight1*IM_dens(iLat1,iLon1) &
                +              LatWeight2*IM_dens(iLat2,iLon1) ) + &
                LonWeight2 * ( LatWeight1*IM_dens(iLat1,iLon2) &
                +              LatWeight2*IM_dens(iLat2,iLon2) ) )
           if(DoMultiFluidIMCoupling)then
              pIm_IC(2,i,j,k) = Si2No_V(UnitP_)*( &
                   LonWeight1 * ( LatWeight1*IM_Hpp(iLat1,iLon1) &
                   +              LatWeight2*IM_Hpp(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_Hpp(iLat1,iLon2) &
                   +              LatWeight2*IM_Hpp(iLat2,iLon2) ) )
              RhoIm_IC(2,i,j,k) = Si2No_V(UnitRho_)*( &
                   LonWeight1 * ( LatWeight1*IM_Hpdens(iLat1,iLon1) &
                   +              LatWeight2*IM_Hpdens(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_Hpdens(iLat1,iLon2) &
                   +              LatWeight2*IM_Hpdens(iLat2,iLon2) ) )
              pIm_IC(3,i,j,k) = Si2No_V(UnitP_)*( &
                   LonWeight1 * ( LatWeight1*IM_Opp(iLat1,iLon1) &
                   +              LatWeight2*IM_Opp(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_Opp(iLat1,iLon2) &
                   +              LatWeight2*IM_Opp(iLat2,iLon2) ) )
              RhoIm_IC(3,i,j,k) = Si2No_V(UnitRho_)*( &
                   LonWeight1 * ( LatWeight1*IM_Opdens(iLat1,iLon1) &
                   +              LatWeight2*IM_Opdens(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_Opdens(iLat1,iLon2) &
                   +              LatWeight2*IM_Opdens(iLat2,iLon2) ) )
           end if
           ! ppar at minimum B
           if(DoAnisoPressureIMCoupling)then
              PparIm_C(i,j,k) = Si2No_V(UnitP_)*( &
                   LonWeight1 * ( LatWeight1*IM_ppar(iLat1,iLon1) &
                   +              LatWeight2*IM_ppar(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_ppar(iLat1,iLon2) &
                   +              LatWeight2*IM_ppar(iLat2,iLon2) ) ) 
              BminIm_C(i,j,k) = Si2No_V(UnitB_)*( &
                   LonWeight1 * ( LatWeight1*IM_bmin(iLat1,iLon1) &
                   +              LatWeight2*IM_bmin(iLat2,iLon1) ) + &
                   LonWeight2 * ( LatWeight1*IM_bmin(iLat1,iLon2) &
                   +              LatWeight2*IM_bmin(iLat2,iLon2) ) ) 
           end if

           if(.not. DoAnisoPressureIMCoupling)then
              ! If coupled with RCM; Or if GM is not using anisop: set ppar = p.
              PparIm_C(i,j,k) = pIm_IC(1,i,j,k)
           else
              ! Anisotropic pressure coupling, not dealing with multifluid for now 
              ! Pperp at minimum B
              Pperp = (3.0*pIm_IC(1,i,j,k) - PparIm_C(i,j,k))/2.0
              PperpInvPpar = Pperp/PparIm_C(i,j,k)
              b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
              if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)
              Coeff = 1/(PperpInvPpar &
                   + min(1.0, BminIm_C(i,j,k)/sqrt(sum(b_D**2)))*(1 - PperpInvPpar))
              ! pressures and density at an arbitrary location of a field line.
              pIm_IC(1,i,j,k) = pIm_IC(1,i,j,k)*Coeff &
                   + 2.0*Pperp*(Coeff - 1)*Coeff/3.0
              PparIm_C(i,j,k) = PparIm_C(i,j,k)*Coeff
              RhoIm_IC(1,i,j,k) = RhoIm_IC(1,i,j,k)*Coeff
           end if

           if(dLatSmoothIm > 0.0)then
              ! Go from low to high lat and look for first unset field line
!!! WHAT ABOUT ASCENDING VS DESCENDING ORDER FOR LAT???
              do n = iSize,1,-1
                 if(IM_p(n,iLon1) < 0.0) EXIT
              enddo
              ! Make sure n does not go below 1
              n = max(1, n)
              ! Set TauCoeff as a function of lat distance from unset lines
              ! No adjustment at the unset line, full adjustment if latitude
              ! difference exceeds dLatSmoothIm
              TauCoeffIm_C(i,j,k) = &
                   min( abs(IM_lat(n) - IM_lat(iLat1))/dLatSmoothIm, 1.0 )
           end if
        end if
     end if

     ! If the pressure is not set by IM, and DoFixPolarRegion is true
     ! and the cell is within radius rFixPolarRegion and flow points outward
     ! then nudge the pressure (and density) towards the "polarregion" values
     if(pIm_IC(1,i,j,k) < 0.0 .and. DoFixPolarRegion .and. &
          R_BLK(i,j,k,iBlock) < rFixPolarRegion .and. &
          Xyz_DGB(z_,i,j,k,iBlock)*State_VGB(RhoUz_,i,j,k,iBlock) > 0)then
        do iFluid = 1, nIons
           pIm_IC(iFluid,i,j,k) = PolarP_I(iFluid)
           RhoIm_IC(iFluid,i,j,k) = PolarRho_I(iFluid)
        end do
     end if

  end do; end do; end do

end subroutine get_im_pressure

!==========================================================================

subroutine apply_im_pressure

  use ModMain, ONLY: nI, nJ, nK, nBlock, Unused_B, iNewGrid, TauCoupleIm, &
       time_accurate, Dt, DoCoupleImPressure, DoCoupleImDensity, RhoMinDimIm
  use ModAdvance, ONLY: State_VGB, UseAnisoPressure
  use ModVarIndexes, ONLY: Ppar_
  use ModPhysics, ONLY: Io2No_V, UnitT_, UnitRho_
  use ModImPressure
  use ModMultiFluid, ONLY : IonFirst_, IonLast_, iRho_I, iP_I, &
       iRhoUx_I, iRhoUy_I, iRhoUz_I, iFluid
  use ModEnergy, ONLY: calc_energy_cell

  implicit none

  real :: Factor

  real :: RhoIm_IC(3,nI,nJ,nK), RhoMinIm
  real :: pIm_IC(3,nI,nJ,nK)
  real :: TauCoeffIm_C(nI,nJ,nK)
  real :: PparIm_C(nI,nJ,nK)

  integer :: iLastPIm = -1, iLastGrid = -1
  integer :: iIonSecond, nIons
  integer :: iBlock

  character (len=*), parameter :: NameSub='apply_im_pressure'
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------
  if(iNewPIm < 1) RETURN ! No IM pressure has been obtained yet

  ! Are we coupled at all?
  if(.not.DoCoupleImPressure .and. .not.DoCoupleImDensity) RETURN  

  call set_oktest(NameSub, DoTest, DoTestMe)

  iIonSecond = min(IonFirst_+1, IonLast_)

  ! Set density floor in normalized units
  if(DoCoupleImDensity) RhoMinIm = Io2No_V(UnitRho_)*RhoMinDimIm

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
     ! A typical value might be 5, to get close to the IM pressure 
     ! in 10 seconds

     Factor = 1.0/(TauCoupleIM*Io2No_V(UnitT_))

  else
     ! Ramp up is based on number of iterations: p' = (ntau*p + pIm)/(1+ntau)
     ! A typical value might be 20, to get close to the IM pressure 
     ! in 20 iterations

     Factor = 1.0/(1 + TauCoupleIM)

  end if

  if (DoMultiFluidIMCoupling)then
     nIons = iIonSecond
  else
     nIons = 1
  end if

  do iBlock = 1, nBlock
     if(Unused_B(iBlock)) CYCLE
     
     call get_im_pressure(iBlock, pIm_IC, RhoIm_IC, TauCoeffIm_C, PparIm_C)
     if(all(pIm_IC < 0.0)) CYCLE  ! Nothing to do

     !Put velocity into momentum temporarily when density is changed
     if(DoCoupleImDensity)then
        do iFluid = 1, nIons
           where(RhoIm_IC(iFluid,:,:,:) > 0.0)
              State_VGB(iRhoUx_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUx_I(iFluid),1:nI,1:nJ,1:nK,iBlock)/ &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
              State_VGB(iRhoUy_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUy_I(iFluid),1:nI,1:nJ,1:nK,iBlock)/ &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
              State_VGB(iRhoUz_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUz_I(iFluid),1:nI,1:nJ,1:nK,iBlock)/ &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
           end where
        end do
     end if

     if(time_accurate)then
        if(DoCoupleImPressure)then
           do iFluid = 1, nIons
              where(pIm_IC(iFluid,:,:,:) > 0.0) &
                   State_VGB(iP_I(iFluid),1:nI,1:nJ,1:nK,iBlock) = &
                   State_VGB(iP_I(iFluid),1:nI,1:nJ,1:nK,iBlock)   &
                   + min(1.0, Factor * Dt) * TauCoeffIm_C &
                   * (pIm_IC(iFluid,:,:,:) - &
                   State_VGB(iP_I(iFluid),1:nI,1:nJ,1:nK,iBlock))
           end do
           if(UseAnisoPressure)then
              where(PparIm_C(:,:,:) > 0.0) &
                   State_VGB(Ppar_,1:nI,1:nJ,1:nK,iBlock) = &
                   State_VGB(Ppar_,1:nI,1:nJ,1:nK,iBlock)   &
                   + min(1.0, Factor * Dt) * TauCoeffIm_C &
                   * (PparIm_C(:,:,:) - &
                   State_VGB(Ppar_,1:nI,1:nJ,1:nK,iBlock))
           end if
        end if
        if(DoCoupleImDensity)then
           do iFluid = 1, nIons
              where(RhoIm_IC(iFluid,:,:,:) > 0.0) &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock) = &
                   max( RhoMinIm, &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)   &
                   + min(1.0, Factor * Dt) * TauCoeffIm_C &
                   * (RhoIm_IC(iFluid,:,:,:) - &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)))
           end do
        end if
     else
        if(DoCoupleImPressure)then
           do iFluid = 1, nIons
              where(pIm_IC(iFluid,:,:,:) > 0.0) &
                   State_VGB(iP_I(iFluid),1:nI,1:nJ,1:nK,iBlock) = Factor* &
                   (TauCoupleIM*State_VGB(iP_I(iFluid),1:nI,1:nJ,1:nK,iBlock)+&
                   pIm_IC(iFluid,:,:,:))
           end do
           if(UseAnisoPressure)then
              where(PparIm_C(:,:,:) > 0.0) &
                   State_VGB(Ppar_,1:nI,1:nJ,1:nK,iBlock) = Factor* &
                   (TauCoupleIM*State_VGB(Ppar_,1:nI,1:nJ,1:nK,iBlock) + &
                   PparIm_C(:,:,:))
           end if
        end if
        if(DoCoupleImDensity)then
           do iFluid = 1, nIons
              where(RhoIm_IC(iFluid,:,:,:) > 0.0) &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock) = &
                   max(RhoMinIm, Factor*( &
                   TauCoupleIM*State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)&
                   + RhoIm_IC(iFluid,:,:,:)))
           end do
        end if
     end if
     !Convert back to momentum
     if(DoCoupleImDensity)then
        do iFluid = 1, nIons
           where(RhoIm_IC(iFluid,:,:,:) > 0.0)
              State_VGB(iRhoUx_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUx_I(iFluid),1:nI,1:nJ,1:nK,iBlock)* &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
              State_VGB(iRhoUy_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUy_I(iFluid),1:nI,1:nJ,1:nK,iBlock)* &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
              State_VGB(iRhoUz_I(iFluid),1:nI,1:nJ,1:nK,iBlock)= &
                   State_VGB(iRhoUz_I(iFluid),1:nI,1:nJ,1:nK,iBlock)* &
                   State_VGB(iRho_I(iFluid),1:nI,1:nJ,1:nK,iBlock)
           end where
        end do
     end if

     ! Now get the energy that corresponds to these new values
     call calc_energy_cell(iBlock)

  end do

end subroutine apply_im_pressure
