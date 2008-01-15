
subroutine solve

  use ModRIM
  use ModParamRIM
  use ModNumConst, only : cPi
  use ModLinearSolver, ONLY: gmres, prehepta, Uhepta, Lhepta

  implicit none

  real, dimension(0:nLons+1,nLats) :: sinTheta
  real, dimension(:), allocatable :: x, y, rhs, b

  integer :: iLat, iLon, nTotalSolve, iI
  logical :: IsLowLat, IsHighLat, DoTest

  real :: Residual
  integer :: nIteration, iError

  external :: matvec_RIM

  if (.not. DoSolve) return

  sinTheta = cPi - Latitude

  nLatsSolve = 0

  DoTouchNorthPole = .false.
  DoTouchSouthPole = .false.

!!!  if (LowLatBoundary /= 0) then
!!!
!!!     ! We are doing two distinctly different regions (like old solver)
!!!
!!!     do iLat = nLats/2, nLats
!!!        if ( Latitude(1,iLat) >= LowLatBoundary .and. &
!!!             Latitude(1,iLat) <= HighLatBoundary) nLatsSolve = nLatsSolve + 1
!!!     enddo
!!!
!!!     nLatsSolve = nLatsSolve * 2
!!!
!!!     SolveType = SolveWithOutEquator_
!!!
!!!     if (HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats)) then
!!!        DoTouchNorthPole = .true.
!!!        DoTouchSouthPole = .true.
!!!     endif
!!!
!!!  else
!!!
!!!     ! We are doing at least the low latitude region
!!!
!!!     if (.not.DoFold) then 

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              nLatsSolve=nLatsSolve+1
!           else
!              Potential(:,iLat) = 10.0
           endif
        enddo

        SolveType = SolveAcrossEquator_

        write(*,*) "hlb : ",HighLatBoundary, Latitude(1,nLats), &
             dLatitude(1,nLats), Latitude(1,nLats)+dLatitude(1,nLats)
        if (HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats) .and. &
             HighLatBoundary > Latitude(1,nLats)) then
           DoTouchNorthPole = .true.
           DoTouchSouthPole = .true.
        endif

        write(*,*) "Across the equator : ", nLatsSolve

!!!     else
!!!
!!!        ! Hmmm....  I am not sure how to define OCFLB yet...
!!!
!!!        SolveType = SolveWithFold_
!!!
!!!        if (HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats)) then
!!!           DoTouchNorthPole = .true.
!!!           DoTouchSouthPole = .true.
!!!        endif
!!!
!!!     endif
!!!
!!!  endif

  ! Don't need Ghostcells here (I think)
  nTotalSolve = nLatsSolve*nLons

  allocate( x(nTotalSolve), y(nTotalSolve), rhs(nTotalSolve), &
       b(nTotalSolve), d_I(nTotalSolve), e_I(nTotalSolve), &
       e1_I(nTotalSolve), f_I(nTotalSolve), f1_I(nTotalSolve) )

  iI = 0

!!!  select case(SolveType)
!!!
!!!     case(SolveAcrossEquator_)

        ! Northern Hemisphere
        IsLowLat = .true.
        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              IsHighLat = .false.
              if (iLat == nLats) then
                 IsHighLat = .true.
              else
                 if (Latitude(1,iLat+1) >= HighLatBoundary) then
                    IsHighLat = .true.
                 endif
              endif
              do iLon = 1, nLons
                 call fill
              enddo
              IsLowLat = .false.
           endif
        enddo

!!!     case(SolveWithOutEquator_)
!!!
!!!        ! We have to do two separate solves - but can do this at once.
!!!
!!!        ! Southern Hemisphere
!!!        IsLowLat = .true.
!!!        do iLat = 1, nLats/2
!!!           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
!!!                Latitude(1,iLat) <= -LowLatBoundary) then
!!!              IsHighLat = .false.
!!!              if (Latitude(1,iLat+1) >= -LowLatBoundary) then
!!!                 IsHighLat = .true.
!!!              endif
!!!              do iLon = 1, nLons
!!!                 call fill
!!!              enddo
!!!              IsLowLat = .false.
!!!           endif
!!!        enddo
!!!
!!!        ! Northern Hemisphere
!!!        IsLowLat = .true.
!!!        do iLat = nLats/2, nLats
!!!           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
!!!                Latitude(1,iLat) <= HighLatBoundary) then
!!!              IsHighLat = .false.
!!!              if (iLat == nLats) then
!!!                 IsHighLat = .true.
!!!              else
!!!                 if (Latitude(1,iLat+1) >= HighLatBoundary) then
!!!                    IsHighLat = .true.
!!!                 endif
!!!              endif
!!!              do iLon = 1, nLons
!!!                 call fill
!!!              enddo
!!!              IsLowLat = .false.
!!!           endif
!!!        enddo
!!!
!!!     case(SolveWithFold_)
!!!        write(*,*) "Doesn't Work!"
!!!
!!!  end select

        call ridley_solve

!!!  write(*,*) "iI : ", iI, nTotalSolve
!!!
!!!  Rhs = b
!!!  if (DoPrecond) then
!!!     ! A -> LU
!!!
!!!     write(*,*) "precon"
!!!     call prehepta(nTotalSolve,1,nLatsSolve,nTotalSolve,-0.5,d_I,e_I,f_I,e1_I,f1_I)
!!!
!!!     ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs
!!!
!!!     ! rhs'=U^{-1}.L^{-1}.rhs
!!!     call Lhepta(       nTotalSolve,1,nLatsSolve,nTotalSolve,b,d_I,e_I,e1_I)
!!!     call Uhepta(.true.,nTotalSolve,1,nLatsSolve,nTotalSolve,b,    f_I,f1_I)
!!!     
!!!  end if
!!!
!!!  if (iDebugLevel > 2) &
!!!       write(*,*)'after precond: sum(b,abs(b),x,d,e,f,e1,f1)=',&
!!!       sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
!!!       sum(e1_I),sum(f1_I)
!!!
!!!  ! Solve A'.x = rhs'
!!!  Residual    = Tolerance
!!!  nIteration  = MaxIteration
!!!  if (.not.UseInitialGuess) x = 0.0
!!!  if (iDebugLevel > 1) DoTest = .true.
!!!
!!!  call gmres(matvec_RIM,b,x,UseInitialGuess,nTotalSolve,&
!!!       MaxIteration,Residual,'rel',nIteration,iError,DoTest)
!!!  if (iError /= 0 .and. iError /=3) &
!!!       write(*,*)'iono_solve: iter, resid, iError=',&
!!!       nIteration, Residual, iError
!!!  if (iError /= 0 .and. iError /=3)then
!!!     write(*,*)'IE_ERROR in iono_solve: gmres failed !!!'
!!!     if(iError < 0) &
!!!          call CON_stop('IE_ERROR in iono_solve: residual did not decrease')
!!!  end if
!!!
!!!  iI = 0
!!!  select case(SolveType)
!!!
!!!     case(SolveAcrossEquator_)
!!!
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
!!!                Latitude(1,iLat) <= HighLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 iI = iI + 1
!!!                 Potential(iLon, iLat) = x(iI)
!!!                 if (iLon == 25) write(*,*) "x(iI) : ",iI, iLon, iLat, x(iI)
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!     case(SolveWithOutEquator_)
!!!
!!!        ! South
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
!!!                Latitude(1,iLat) <= -LowLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 iI = iI + 1
!!!                 Potential(iLon, iLat) = x(iI)
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!        ! North
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
!!!                Latitude(1,iLat) <= HighLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 iI = iI + 1
!!!                 Potential(iLon, iLat) = x(iI)
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!     case(SolveWithFold_)
!!!        write(*,*) "Doesn't Work!"
!!!  end select
!!!
!!!  ! If we include poles, then the pole solution is the average of all
!!!  ! the cells around the pole:
!!!
!!!  if (DoTouchNorthPole) &
!!!       Potential(1:nLons,nLats) = sum(Potential(1:nLons,nLats-1))/nLons
!!!  if (DoTouchSouthPole) &
!!!       Potential(1:nLons,    1) = sum(Potential(1:nLons,      2))/nLons
!!!
!!!  ! Need to do message passing here......
!!!
!!!  ! Periodic Boundary Conditions:
!!!
!!!  Potential(      0,:) = Potential(nLons,:)
!!!  Potential(nLons+1,:) = Potential(    1,:)

  if(allocated(b)) deallocate(x, y, b, rhs, d_I, e_I, f_I, e1_I, f1_I)

contains

  subroutine fill

    iI = iI + 1
    b(iI)    = Jr(iLon,iLat)*(Radius*sinTheta(iLon,iLat))**2
    x(iI)    = OldPotential(iLon,iLat)
    d_I(iI)  = SolverA(iLon,iLat)
    e_I(iI)  = SolverB(iLon,iLat)
    f_I(iI)  = SolverC(iLon,iLat)
    e1_I(iI) = SolverD(iLon,iLat)
    f1_I(iI) = SolverE(iLon,iLat)

    if (IsLowLat  ) e_I(iI)  = 0.0
    if (IsHighLat ) f_I(iI)  = 0.0
    if (iLon == 1    ) e1_I(iI) = 0.0
    if (iLon == nLons) f1_I(iI) = 0.0

  end subroutine fill

end subroutine solve

!-----------------------------------------------------------------------
! matvec routine
!-----------------------------------------------------------------------

subroutine matvec_RIM(x_I, y_I, n)

  use ModRIM
  use ModParamRIM
  use ModLinearsolver, ONLY: Uhepta, Lhepta

  implicit none

  integer, intent(in) :: n          ! number of unknowns
  real, intent(in) :: x_I(n)        ! vector of unknowns
  real, intent(out):: y_I(n)        ! y = A.x

  real :: x_G(0:nLons+1, nLats)

  integer :: iI, iLon, iLat

  iI = 0

  x_G = Potential

!!!  select case(SolveType)
!!!
!!!     case(SolveAcrossEquator_)

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 x_G(iLon, iLat) = x_I(iI)
              enddo
           endif
        enddo

!!!     case(SolveWithOutEquator_)
!!!
!!!        ! South
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
!!!                Latitude(1,iLat) <= -LowLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 iI = iI + 1
!!!                 x_G(iLon, iLat) = x_I(iI)
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!        ! North
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
!!!                Latitude(1,iLat) <= HighLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 iI = iI + 1
!!!                 x_G(iLon, iLat) = x_I(iI)
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!     case(SolveWithFold_)
!!!        write(*,*) "Doesn't Work!"
!!!  end select

  ! If we include poles, then the pole solution is the average of all
  ! the cells around the pole:

!  if (DoTouchNorthPole) x_G(1:nLons,nLats) = sum(x_G(1:nLons,nLats-1))/nLons
!  if (DoTouchSouthPole) x_G(1:nLons,    1) = sum(x_G(1:nLons,      2))/nLons

  ! Need to do message passing here......

  ! Periodic Boundary Conditions:

  x_G(      0,:) = x_G(nLons,:)
  x_G(nLons+1,:) = x_G(    1,:)

  iI = 0
  
!!!  select case(SolveType)
!!!
!!!     case(SolveAcrossEquator_)

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 call fill
              enddo
           endif
        enddo

        write(*,*) iI

!!!     case(SolveWithOutEquator_)
!!!
!!!        ! South
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
!!!                Latitude(1,iLat) <= -LowLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 call fill
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!        ! North
!!!        do iLat = 1, nLats
!!!           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
!!!                Latitude(1,iLat) <= HighLatBoundary) then
!!!              do iLon = 1, nLons
!!!                 call fill
!!!              enddo
!!!           endif
!!!        enddo
!!!
!!!     case(SolveWithFold_)
!!!        write(*,*) "Doesn't Work!"
!!!  end select

  ! Preconditioning: y'= U^{-1}.L^{-1}.y
  if(DoPrecond)then
     call Lhepta(       n,1,nLatsSolve,n,y_I,d_I,e_I,e1_I)
     call Uhepta(.true.,n,1,nLatsSolve,n,y_I,    f_I,f1_I)
  end if

contains

  subroutine fill

    iI = iI + 1
    
    y_I(iI) = &
         SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
         SolverB(iLon, iLat)*x_G(iLon,  iLat-1) + &
         SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
         SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
         SolverE(iLon, iLat)*x_G(iLon+1,iLat  )

    if (iLon == 25) write(*,*) "iLat : ",iLat, &
         x_G(iLon,  iLat-1), x_G(iLon,  iLat), x_G(iLon,  iLat+1), &
         SolverA(iLon, iLat), SolverB(iLon, iLat), SolverC(iLon, iLat), &
         SolverD(iLon, iLat), SolverE(iLon, iLat), y_I(iI)

  end subroutine fill

end subroutine matvec_RIM

subroutine ridley_solve

  use ModRIM
  use ModParamRIM

  implicit none

  integer :: iLon, iLat, nIters
  real :: OldSolution(0:nLons+1,nLats), Solution(0:nLons+1,nLats)
  real :: Residual

  ! Don't worry about high latitude (i.e. Pole) BCs yet

  Solution = Potential
  Solution(0,iLat) = Solution(nLons,iLat)
  Solution(nLons+1,iLat) = Solution(1,iLat)
  OldSolution = Solution
  Residual = 1.0e32
  nIters = 0

  do while (Residual > 500.0)

     do iLat = 2, nLats-1
        if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
             Latitude(1,iLat) <= HighLatBoundary) then
           do iLon = 1, nLons
              Solution(iLon,iLat) = &
                   (SolverB(iLon, iLat)*Solution(iLon,  iLat-1) + &
                    SolverC(iLon, iLat)*Solution(iLon,  iLat+1) + &
                    SolverD(iLon, iLat)*Solution(iLon-1,iLat  ) + &
                    SolverE(iLon, iLat)*Solution(iLon+1,iLat  )) / &
                    (-SolverA(iLon, iLat))

!              if (iLon == 25) write(*,*) "sol : ",iLat,&
!                   Solution(iLon,iLat)

           enddo
           Solution(0,iLat) = Solution(nLons,iLat)
           Solution(nLons+1,iLat) = Solution(1,iLat)
        endif
     enddo

     Residual = sum(abs(Solution-OldSolution))

     OldSolution = Solution
     nIters = nIters + 1

     write(*,*) "Residual : ", Residual, nIters

  enddo

  Potential = Solution

end subroutine ridley_solve
