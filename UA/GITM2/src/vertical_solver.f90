!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine advance_vertical_1d

  use ModVertical
  use ModGITM, ONLY : Dt, iCommGITM, iProc
  use ModInputs, only: UseBarriers, iDebugLevel
  implicit none
  !-----------------------------------------------------------

  integer :: iError

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc

  ! Fill in ghost cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)

  ! Copy input state into New state
  NewLogNS  = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp   = Temp
  NewVertVel = VertVel

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) write(*,*) "========> stage 1", iproc

  ! Do the half step: U^n+1/2 = U^n + (Dt/2) * R(U^n)
!  Dt = Dt/2
!
!  write(*,*) "vv, before av1s 1: ", VertVel(3, :),LogNS(3,1)
!
!  call advance_vertical_1stage(&
!       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
!       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!
!  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
!  if (iDebugLevel > 7) write(*,*) "========> vertical bcs 2", iproc
!
!  ! Fill in ghost cells for U^n+1/2 state
!  call set_vertical_bcs(NewLogRho, NewLogNS, NewVel_GD, NewTemp,NewLogINS, &
!       IVel, NewVertVel)
!
!  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
!  if (iDebugLevel > 7) write(*,*) "========> stage 2", iproc
!
!  write(*,*) "vv, before av1s 2: ", NewVertVel(3, :),NewLogNS(3,1)
!
!  ! Do full step U^n+1 = U^n + Dt * R(U^n+1/2)
!  Dt = 2*Dt
  call advance_vertical_1stage(&
       NewLogRho, NewLogNS, NewVel_GD, NewTemp, LogRho, LogNS, Vel_GD, Temp, &
       NewLogINS, LogINS, IVel, NewVertVel, VertVel)
  
!  write(*,*) "vv, after av1s 2: ", VertVel(3, :),LogNS(3,1)

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) write(*,*) "========> vertical bcs 3", iproc

  ! Fill in ghost cells for updated U^n+1 state
  call set_vertical_bcs(LogRho, LogNS, Vel_GD, Temp, LogINS, IVel, VertVel)

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d

!=============================================================================
subroutine advance_vertical_1stage( &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: RadialDistance, &
       Dt, iO_, dAlt, Gravity, Altitude, iEast_, iNorth_, iUp_, TempUnit
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       Heating, KappaNS, nDensityS, KappaTemp, Centrifugal, Coriolis
  use ModTime
  use ModInputs
  use ModConstants
  implicit none

  real, intent(in) :: LogRho(-1:nAlts+2)
  real, intent(in) :: LogNS(-1:nAlts+2,nSpecies)
  real, intent(in) :: LogINS(-1:nAlts+2,nIonsAdvect)
  real, intent(in) :: Vel_GD(-1:nAlts+2,3)
  real, intent(in) :: IVel(-1:nAlts+2,3)
  real, intent(in) :: Temp(-1:nAlts+2)
  real, intent(in) :: VertVel(-1:nAlts+2,nSpecies)

  real, intent(inout) :: NewLogRho(-1:nAlts+2)
  real, intent(inout) :: NewLogNS(-1:nAlts+2,nSpecies)
  real, intent(inout) :: NewLogINS(-1:nAlts+2,nIonsAdvect)
  real, intent(inout) :: NewVel_GD(-1:nAlts+2,3)
  real :: NewVel2_G(-1:nAlts+2)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)
  real :: NS(-1:nAlts+2,nSpecies)
  real :: Rho(-1:nAlts+2)

  real :: TempKoM(-1:nAlts+2), AveMass(-1:nAlts+2), LogNum(-1:nAlts+2)

  real, dimension(1:nAlts)    :: GradLogRho, DivVel, GradTemp, GradTempKoM, &
       DiffLogRho, DiffTemp, GradTmp, DiffTmp, DiffLogNum, GradLogNum
  real, dimension(1:nAlts,3) :: GradVel_CD, DiffVel_CD

  real, dimension(1:nAlts,nSpecies)    :: GradLogNS, DiffLogNS, &
       GradVertVel, DiffVertVel, DivVertVel
  real, dimension(1:nAlts,nIonsAdvect) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho, rat

  integer :: iAlt, iSpecies, jSpecies, iDim
  !--------------------------------------------------------------------------

  NS = exp(LogNS)
  Rho = exp(LogRho)
  LogNum = alog(sum(NS,dim=2))
  AveMass = Rho/sum(NS,dim=2)
  TempKoM = Temp*TempUnit * Boltzmanns_Constant / AveMass

  call calc_rusanov(LogRho ,GradLogRho,  DiffLogRho)
  call calc_rusanov(LogNum ,GradLogNum,  DiffLogNum)
  call calc_rusanov(Temp   ,GradTemp,    DiffTemp)
  call calc_rusanov(TempKoM,GradTempKoM, DiffTmp)
  do iDim = 1, 3
     call calc_rusanov(Vel_GD(:,iDim), GradVel_CD(:,iDim), DiffVel_CD(:,iDim))
  end do

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)/RadialDistance(1:nAlts)

  do iSpecies=1,nSpecies
     call calc_rusanov(LogNS(:,iSpecies),GradTmp, DiffTmp)
     GradLogNS(:,iSpecies) = GradTmp
     DiffLogNS(:,iSpecies) = DiffTmp
     call calc_rusanov(VertVel(:,iSpecies),GradTmp, DiffTmp)
     GradVertVel(:,iSpecies) = GradTmp
     DiffVertVel(:,iSpecies) = DiffTmp
     DivVertVel(:,iSpecies) = GradVertVel(:,iSpecies) + &
          2*VertVel(1:nAlts,iSpecies)/RadialDistance(1:nAlts)
  enddo

  do iSpecies=1,nIonsAdvect
     call calc_rusanov(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  do iAlt = 1,nAlts

     NewLogRho(iAlt) = NewLogRho(iAlt) - Dt * &
          (DivVel(iAlt) + Vel_GD(iAlt,iUp_) * GradLogRho(iAlt) ) &
          + Dt * DiffLogRho(iAlt)

     do iSpecies=1,nSpecies

!        NewLogNS(iAlt,iSpecies) = NewLogNS(iAlt,iSpecies) - Dt * &
!             (DivVel(iAlt) + Vel_GD(iAlt,iUp_) * GradLogNS(iAlt,iSpecies) ) &
!             + Dt * DiffLogNS(iAlt,iSpecies)
        NewLogNS(iAlt,iSpecies) = LogNS(iAlt,iSpecies) - Dt * &
             (DivVertVel(iAlt,iSpecies) + &
             VertVel(iAlt,iSpecies) * GradLogNS(iAlt,iSpecies) ) &
             + Dt * DiffLogNS(iAlt,iSpecies)
     enddo

     do iSpecies=1,nIonsAdvect

        NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - Dt * &
             (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
             + Dt * DiffLogINS(iAlt,iSpecies)

     enddo

!     ! dVr/dt = -[ (V grad V)_r + grad T + T grad ln Rho - g ]
!     ! and V grad V contains the centripetal acceleration 
!     ! (Vphi**2+Vtheta**2)/R
!
!     NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) - Dt * &
!          (Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iUp_) &
!          - (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
!          / RadialDistance(iAlt) &
!          - Gravity(iAlt)) &
!          + Dt * DiffVel_CD(iAlt,iUp_)

!     NewVel2_G(iAlt) = NewVel_GD(iAlt,iUp_) - Dt * &
!          (GradTempKoM(iAlt) + TempKoM(iAlt)*GradLogRho(iAlt))
!     NewVel_GD(iAlt, iUp_) = NewVel2_G(iAlt) 

     NewVel_GD(iAlt,iUp_) = 0.0

!     rat = Rho(iAlt) / Rho(1)
!     rat = 0.0

     do iSpecies=1,nSpecies

!        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) - Dt * &
!             (exp(LogNS(iAlt,iSpecies)) / (exp(LogRho(iAlt)) / Mass(1))) * &
!             (GradTemp(iAlt) + Temp(iAlt)*GradLogNS(iAlt,iSpecies))

        NewVertVel(iAlt, iSpecies) = VertVel(iAlt, iSpecies) - Dt * &
             (VertVel(iAlt,iSpecies)*GradVertVel(iAlt,iSpecies) &
             - (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
             / RadialDistance(iAlt) &
             - Gravity(iAlt) + &
             Mass(1)/Mass(iSpecies) * &
             (GradTemp(iAlt) + Temp(iAlt)*GradLogNS(iAlt,iSpecies))) &
             + Dt * DiffVertVel(iAlt,iSpecies)

        NewVertVel(iAlt, iSpecies) = max(-500.0, NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( 500.0, NewVertVel(iAlt, iSpecies))

        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies) * &
             (Mass(iSpecies) * NS(iAlt,iSpecies) / Rho(iAlt))

     enddo

     if (UseCoriolis) then
        ! Need to add vertical velocity here!!!!
        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + Dt * ( &
             Centrifugal * RadialDistance(iAlt) + &
             Coriolis * Vel_GD(iAlt,iEast_))
     endif

     ! dVphi/dt = - (V grad V)_phi
     NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iEast_) &
          + Dt * DiffVel_CD(iAlt,iEast_)

     ! dVtheta/dt = - (V grad V)_theta
     NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) - Dt * &
          Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iNorth_) &
          + Dt * DiffVel_CD(iAlt,iNorth_)

     ! dT/dt = -(V.grad T + (gamma - 1) T div V

     NewTemp(iAlt)   = NewTemp(iAlt) - Dt * &
          (Vel_GD(iAlt,iUp_)*GradTemp(iAlt) + &
          (Gamma - 1.0) * Temp(iAlt)*DivVel(iAlt))&
          + Dt * DiffTemp(iAlt)

  end do

!  ! Add molecular diffusion
!  if(UseDiffusion) call add_diffusion

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

!contains

!  subroutine add_diffusion
!
!    ! Colgrove, JGR, 1966, p2229
!
!!    ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
!!    real, parameter, dimension(5, 5) :: Diff0 = 0.0001 * reshape( (/ &
!!         ! 0      02     N2      N     NO
!!         !---------------------------------+
!!         0.00,  0.260, 0.260, 0.260, 0.260, &            ! O
!!         0.26,  0.000, 0.181, 0.181, 0.181, &            ! O2
!!         0.26,  0.181, 0.000, 0.181, 0.181, &            ! N2
!!         0.26,  0.181, 0.181, 0.000, 0.260, &            ! N
!!         0.26,  0.181, 0.181, 0.260, 0.000  /), (/5,5/) )  ! NO
!!
!!    ! These are the exponents - 1
!!    real, parameter, dimension(5, 5) :: DiffExp = reshape( (/ &
!!         ! 0      02     N2
!!         !---------------------------------+
!!         0.00,  0.75,  0.75, 0.75, 0.75, &             ! O
!!         0.75,  0.00,  0.75, 0.75, 0.75, &             ! O2
!!         0.75,  0.75,  0.00, 0.75, 0.75, &             ! N2
!!         0.75,  0.75,  0.75, 0.00, 0.75, &             ! N
!!         0.75,  0.75,  0.75, 0.75, 0.00  /), (/5,5/) )  ! NO
!
!    ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
!    real, parameter, dimension(4, 4) :: Diff0 = 0.0001 * reshape( (/ &
!         ! 0      02     N2      N     NO
!         !---------------------------------+
!         0.00,  0.260, 0.260, 0.260, &            ! O
!         0.26,  0.000, 0.181, 0.181, &            ! O2
!         0.26,  0.181, 0.000, 0.181, &            ! N2
!         0.26,  0.181, 0.181, 0.000 /), (/4,4/) )  ! N
!
!    ! These are the exponents - 1
!    real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
!         ! 0      02     N2
!         !---------------------------------+
!         0.00,  0.75,  0.75, 0.75, &             ! O
!         0.75,  0.00,  0.75, 0.75, &             ! O2
!         0.75,  0.75,  0.00, 0.75, &             ! N2
!         0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N
!
!    real :: NS(nAlts, nSpecies), N0
!    real :: Diffusion(nAlts)
!    real :: Rate(nAlts, nSpecies)
!    !-------------------------------------------------------------------------
!
!    NS = exp(LogNS(1:nAlts, : ))  
!    
!    N0 = sum(NS(1,:))
!
!    Rate = 0
!    do iSpecies = 1, nSpecies
!
!       ! Rate_i = Sum n_j/(N D_ij)
!
!       do jSpecies = 1, nSpecies
!
!          if(jSpecies == iSpecies) CYCLE
!
!          Rate(:,iSpecies) = Rate(:,iSpecies) + &
!               NS(:,jSpecies) /  &
!               (Diff0(iSpecies, jSpecies)                               &
!               * (Temp(1:nAlts)/Temp(1)) ** DiffExp(iSpecies, jSpecies) &
!               * N0 &
!               )
!
!       end do
!    end do
!
!    do iSpecies=1,nSpecies
!
!       Diffusion = &
!            ( GradTemp/Temp(1:nAlts) + GradLogNS(1:nAlts,iSpecies) &
!            - Gravity(1:nAlts)/Temp(1:nAlts)*(Mass(iSpecies)/Mass(1)) &
!            )/ Rate(:,iSpecies)
!
!       do iAlt = 2, nAlts-1
!
!!          write(*,*) "diff : ", iAlt, iSpecies, &
!!               0.5*(Diffusion(iAlt+1)-Diffusion(iAlt-1))/dAlt(iAlt), &
!!               Diffusion(iAlt) * GradLogNS(iAlt, iSpecies), NewLogNS(iAlt,iSpecies)
!
!          NewLogNS(iAlt, iSpecies) = NewLogNS(iAlt,iSpecies) - Dt * &
!              ( 0.5*(Diffusion(iAlt+1)-Diffusion(iAlt-1))/dAlt(iAlt) &
!              + Diffusion(iAlt) * GradLogNS(iAlt, iSpecies))
!       end do
!
!       iAlt = nAlts
!       NewLogNS(iAlt, iSpecies) = NewLogNS(iAlt,iSpecies) + Dt * &
!            ( 0.5*(Diffusion(iAlt)-Diffusion(iAlt-2))/dAlt(iAlt-1) &
!            + Diffusion(iAlt) * GradLogNS(iAlt, iSpecies))
!
!    end do
!
!    !        write(*,*)'NewLogNS=',iSpecies,NewLogNS(50,iSpecies)
!    ! Fix total density
!    !do iAlt = 1, nAlts
!    !   NewLogRho(iAlt) = alog(sum(Mass(1:nSpecies)*&
!    !        exp(NewLogNS(iAlt,1:nSpecies))))
!    !end do
!
!  end subroutine add_diffusion

end subroutine advance_vertical_1stage

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModGITM, only : dAlt, Altitude
  use ModVertical, only : cmax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

!  call start_timing("rusanov")

  call calc_GITM_facevalues(nAlts, Altitude(-1:nAlts+2), Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt(1:nAlts)

!  call end_timing("rusanov")

end subroutine calc_rusanov

!\
! ------------------------------------------------------------
! calc_GITM_facevalues
! ------------------------------------------------------------
!/

subroutine calc_GITM_facevalues(n, Location, Var, VarLeft, VarRight)

  use ModInputs, only: TypeLimiter, UseMinMod, UseMC

  implicit none
  
  integer, intent(in) :: n
  real, intent(in) :: Var(-1:n+2), Location(-1:n+2)
  real, intent(out):: VarLeft(1:n+1), VarRight(1:n+1)

  real :: dVarUp, dVarDown, dVarLimited(0:n+1)

  integer :: i

  logical :: IsFirstWarning = .true.

!  call start_timing("facevalues")

!!! write(*,*)'Var(19:23)=',Var(19:23)

  do i=0,n+1

     dVarUp            = (Var(i+1) - Var(i))  / (Location(i+1) - Location(i))
     dVarDown          = (Var(i)   - Var(i-1))/ (Location(i)   - Location(i-1))

     if (UseMinMod) dVarLimited(i) = Limiter_minmod(dVarUp, dVarDown)

     if (UseMC) dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

!!! if(i>=20.and. i<23)write(*,*)'dVarUp,dVarDown,dVarLimited=',&
!!!     dVarUp,dVarDown,dVarLimited(i)

  end do

  do i=1,n+1

     VarLeft(i)  = Var(i-1) + 0.5*dVarLimited(i-1) * &
          (Location(i) - Location(i-1))
     VarRight(i) = Var(i)   - 0.5*dVarLimited(i) * &
          (Location(i) - Location(i-1))

  end do

!  call end_timing("facevalues")

contains

  !---------------------------------------------------
  !---------------------------------------------------
  real function Limiter_minmod(dUp, dDown)

    real :: dUp, dDown

    Limiter_minmod = (sign(0.5,dUp) + sign(0.5,dDown))*min(abs(dUp),abs(dDown))

  end function Limiter_minmod

  !---------------------------------------------------
  !---------------------------------------------------
  real function Limiter_mc(dUp, dDown)

    real :: dUp, dDown

    if (dUp > 0.0) then
       if (dDown > 0.0) then
          Limiter_mc = min(2*dUp,2*dDown,(dUp+dDown)*0.5)
       else
          Limiter_mc = 0.0
       endif
    else
       if (dDown < 0.0) then
          Limiter_mc = max(2*dUp,2*dDown,(dUp+dDown)*0.5)
       else
          Limiter_mc = 0.0
       endif
    endif

  end function Limiter_mc

end subroutine calc_GITM_facevalues

