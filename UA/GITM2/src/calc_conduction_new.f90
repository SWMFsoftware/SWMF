subroutine calc_conduction(iBlock, Quantity, Diff, MulFac, dTdt_cond)

  use ModSizeGitm
  use ModGITM, only: dAlt_GB, Latitude, Longitude, dt, Altitude_gb
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  real, intent(in) :: Quantity(nLons, nLats, -1:nAlts+2)
  real, intent(in) :: Diff(nLons, nLats, 0:nAlts+1)
  real, intent(in) :: MulFac(nLons, nLats, nAlts)
  real, intent(out) :: dTdt_cond(nLons, nLats, nAlts)
  real :: solution(nLons, nLats, 0:nAlts+1)

  real :: beta(nAlts), dk(nAlts), tempold(nAlts), temp(0:nAlts+1)
  real, dimension(nAlts) :: a,b,c,r,u, gam
  real :: tempsave(0:nAlts+1), one_2kb(nAlts), bet
  integer :: iLon, iLat, iAlt, nIters
  real :: tolerance = 0.01, residual, maxresidual

  real :: f, dau, dal

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


  !  write(*,*) "MaxTemp In : ", maxval(quantity)

  do iLon = 1, nLons
     do iLat = 1, nLats

        do iAlt = 1, nAlts

           ! F is defined as the ratio between the upper and lower 
           ! grid spacing.
           
           dau =Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt,iBlock)
           dal =Altitude_GB(iLon,iLat,iAlt,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt-1,iBlock)
           f = dau/dal

           ! Then we use dal (delta-altitude lower) in all of the eqns.

           beta(iAlt) = (dt / MulFac(iLon, iLat, iAlt)) / dal**2

           dk(iAlt) =   Diff(iLon,iLat,iAlt+1) - &
                      f*Diff(iLon,iLat,iAlt-1) - &
                (1-f*f)*Diff(iLon,iLat,iAlt)

           ! dk(iAlt) = 0.0

           tempold(iAlt) = Quantity(iLon, iLat, iAlt)
           one_2kb(iAlt) = &
                1.0 + beta(iAlt)*(&
                (1+f)*Diff(iLon,iLat,iAlt) + &
                (1-f*f)/(f*(1+f))*dk(iAlt))

           a(iAlt) = -beta(iAlt) * &
                (f*Diff(iLon,iLat,iAlt) - dk(iAlt)/(f+1))

           c(iAlt) = -beta(iAlt) * &
                (Diff(iLon,iLat,iAlt) + dk(iAlt)/(f*(f+1)))

           b(iAlt) = one_2kb(iAlt)

           r(iAlt) = tempold(iAlt)

        enddo

        ! Boundary Conditions:
        r(1) = r(1) - a(1) * Quantity(iLon, iLat, 0)
        a(1) = 0.0

        b(nAlts) = b(nAlts) + c(nAlts)
        c(nAlts) = 0.0

        ! This solver comes out of Numerical Recipies, page 48.

        bet = b(1)
        u(1) = r(1)/bet
        do iAlt = 2,nAlts
           gam(iAlt) = c(iAlt-1)/bet
           bet = b(iAlt)-a(iAlt)*gam(iAlt)
          If (bet.eq.0.0) then
                call stop_gitm("Error in tridiaginal solver in calc_cond")
           endif
           u(iAlt) = (r(iAlt)-a(iAlt)*u(iAlt-1))/bet
        enddo
        Do iAlt = nAlts-1,1,-1
           u(iAlt) = u(iAlt)-gam(iAlt+1)*u(iAlt+1)
        enddo

        temp(1:nAlts) = u

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

              dk(iAlt) = ( &
                   Diff(iLon,iLat,iAlt+1) - &
                   Diff(iLon,iLat,iAlt-1))/2

              tempold(iAlt) = Quantity(iLon,iLat,iAlt) + &
                   (dk(iAlt) * 0.5 * &
                   (solution(iLon,iLat,iAlt+1)-solution(iLon,iLat,iAlt-1)) + &
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
