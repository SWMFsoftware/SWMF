subroutine calc_conduction(iBlock, Quantity, Diff, MulFac, dTdt_cond)

  use ModSizeGitm
  use ModGITM, only: dAlt_GB, Latitude, Longitude, dt, Altitude_GB
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  real, intent(in) :: Quantity(nLons, nLats, -1:nAlts+2)
  real, intent(in) :: Diff(nLons, nLats, 0:nAlts+1)
  real, intent(in) :: MulFac(nLons, nLats,-1:nAlts+2)
  real, intent(out) :: dTdt_cond(nLons, nLats, nAlts)
  real :: solution(nLons, nLats, 0:nAlts+1)

  real :: beta(0:nAlts+1),  tempold(-1:nAlts+2), temp(0:nAlts+1)
  real, dimension(nAlts+2) :: a,b,c,r,u, gam
  real :: tempsave(0:nAlts+1), one_2kb(0:nAlts+1), bet
  integer :: iLon, iLat, iAlt, nIters
  real :: tolerance = 0.01, residual, maxresidual
  real :: theta 
  real :: F, dau, dal

  real :: eta(0:nAlts+1), betaP(0:nAlts+1), betaN(0:nAlts+1)
  real :: zeta(0:nAlts+1), gammaP(0:nAlts+1), gammaN(0:nAlts+1)

  logical :: NewSolver = .true., UseDebug = .false.

  call start_timing("conduction")
  call report("calc_conduction",3)

  !C-----------------------------------------------------------------------
  !      Subroutine tridag(a,b,c,r,u,n)
  !      Implicit Double Precision (A-z)
  !      Integer j, n, NMAX
  !      real :: a(n), b(n), c(n), r(n), u(n)
  !      real :: gam(n)
  !      Parameter (NMAX=500)
  !      !
  !      If (b(1).eq.0.D0) Pause 'tridag: rewrite equations'
  !      bet = b(1)
  !      u(1) = r(1)/bet
  !      Do 11 j = 2,n
  !         gam(j) = c(j-1)/bet
  !         bet = b(j)-a(j)*gam(j)
  !         If (bet.eq.0.D0) Pause 'tridag failed'
  !         u(j) = (r(j)-a(j)*u(j-1))/bet
  !   11 Continue
  !      Do 12 j = n-1,1,-1
  !         u(j) = u(j)-gam(j+1)*u(j+1)
  !   12 Continue
  !      Return
  !      End
  !C-----------------------------------------------------------------------


  !        write(*,*) "MaxTemp In : ", maxval(quantity)
  !        write(*,*) ' Testing in Calc_Conduction'

  theta = 1.0

  do iLon = 1, nLons
     do iLat = 1, nLats

!!! First, Set the "Old Data"
        do iAlt = -1, nAlts +2
           tempold(iAlt) = Quantity(iLon, iLat, iAlt)
        enddo 

        do iAlt = 0, nAlts +1

           dau =Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt,iBlock)

           dal =Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt-1,iBlock)

           F = (dau + dal)/2.0

!!!! Beta is now the Kt*Dt/(cp*rho) * (2.0/(dX+ + dX-))

           eta(iAlt) =          theta*(Dt*Diff(iLon,iLat,iAlt)) / MulFac(iLon, iLat, iAlt) / F
           zeta(iAlt) = (1.0 - theta)*(Dt*Diff(iLon,iLat,iAlt)) / MulFac(iLon, iLat, iAlt) / F

!           eta(iAlt) =          0.5*(Dt*Diff(iLon,iLat,iAlt)) / (F* MulFac(iLon, iLat, iAlt) )
!           zeta(iAlt) =         0.5*(Dt*Diff(iLon,iLat,iAlt)) / (F*MulFac(iLon, iLat, iAlt) )

           betaP(iAlt) =  eta(iAlt)/dau
           betaN(iAlt) =  eta(iAlt)/dal

           gammaP(iAlt) =  zeta(iAlt)/dau
           gammaN(iAlt) =  zeta(iAlt)/dal

           one_2kb(iAlt) = (1.0 + BetaP(iAlt) + BetaN(iAlt))

           a(iAlt) = -betaN(iAlt)

           c(iAlt) = -betaP(iAlt)

           b(iAlt) = one_2kb(iAlt)

           r(iAlt) = gammaP(iAlt)*tempold(iAlt+1) + &
                    (1.0 - gammaP(iAlt) - gammaN(iAlt) )*tempold(iAlt) + &
                     gammaN(iAlt)*tempold(iAlt-1)

        enddo

        ! Boundary Conditions:

        r(0) = tempold(0)
        a(0) = 0.0
        b(0) = 1.0
        c(0) = 0.0

        a(nAlts+1) = 0.0
        b(nAlts+1) = 1.0
        c(nAlts+1) = 0.0
        r(nAlts+1) = tempold(nAlts+1)

        ! This solver comes out of Numerical Recipies, page 48.

        bet = b(0)
        u(0) = r(0)/bet

        do iAlt = 1,nAlts+1
           gam(iAlt) = c(iAlt-1)/bet
           bet = b(iAlt)-a(iAlt)*gam(iAlt)
          If (bet.eq.0.0) then
                call stop_gitm("Error in tridiaginal solver in calc_cond")
           endif
           u(iAlt) = (r(iAlt)-a(iAlt)*u(iAlt-1))/bet
        enddo
        Do iAlt = nAlts,0,-1
           u(iAlt) = u(iAlt)-gam(iAlt+1)*u(iAlt+1)
        enddo

        temp(0:nAlts+1) = u(0:nAlts+1)

        dTdt_cond(iLon,iLat,1:nAlts) = temp(1:nAlts) - &
             Quantity(iLon,iLat,1:nAlts)

     enddo
  enddo

  if (UseDebug) then

     solution(:,:,1:nAlts) = dTdt_cond + Quantity(:,:,1:nAlts)
     solution(:,:,0)       = Quantity(:,:,0)
     solution(:,:,nAlts+1) = solution(:,:,nAlts)

     do iLon = 1, nLons
        do iLat = 1, nLats

           do iAlt = 1, nAlts

              tempold(iAlt) = Quantity(iLon,iLat,iAlt) + &
                   (0.5*(solution(iLon,iLat,iAlt+1)-solution(iLon,iLat,iAlt-1)) + &
                   Diff(iLon,iLat,iAlt)* &
                   (solution(iLon,iLat,iAlt+1) &
                   + solution(iLon,iLat,iAlt-1) &
                   - 2*solution(iLon,iLat,iAlt))) *  &
                   Dt * MulFac(iLon, iLat, iAlt) / &
                   dAlt_GB(iLon, iLat, iAlt, iBlock)**2

              if (abs(solution(iLon,iLat,iAlt) - tempold(iAlt)) > 1.0e-2) then
                 write(*,*) "Error : ", iLon, iLat, iAlt, &
                      solution(iLon,iLat,iAlt), tempold(iAlt), &
                      abs(solution(iLon,iLat,iAlt) - tempold(iAlt))
                 call stop_gitm("Must stop in calc_conduction")
              endif

           enddo
        enddo
     enddo

  endif

  call end_timing("conduction")

end subroutine calc_conduction
