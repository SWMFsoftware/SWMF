
subroutine solve

  use ModRIM
  use ModParamRIM
  use ModNumConst,     only: cPi
  use ModProcIE,       only: iComm
  use ModLinearSolver, only: gmres, prehepta, Uhepta, Lhepta

  implicit none

  real, dimension(0:nLons+1,nLats) :: sinTheta
  real, dimension(:), allocatable :: x, y, rhs, b
  real, dimension(0:nLons+1) :: nPotential, sPotential

  integer :: iLat, iLon, nTotalSolve, iI
  logical :: IsLowLat, IsHighLat, DoTest, DoIdealTest=.false.

  real :: Residual, r
  integer :: nIteration, iError

  external :: matvec_RIM

  if (.not. DoSolve) return


!  DoIdealTest = .true.

  if (DoIdealTest) then
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "               IDEAL TEST"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     jr = -1.0e-6*exp(-abs((abs(Latitude)-75.0*cDegToRad))/(5.0*cDegToRad)) &
          * sin(Longitude)
  endif

  sinTheta = sin(cPi/2 - Latitude)

  nLatsSolve = 0

  DoTouchNorthPole = .false.
  DoTouchSouthPole = .false.

  if (LowLatBoundary /= 0) then

     ! We are doing two distinctly different regions (like old solver)

     do iLat = nLats/2, nLats
        if ( Latitude(1,iLat) >= LowLatBoundary .and. &
             Latitude(1,iLat) <= HighLatBoundary) nLatsSolve = nLatsSolve + 1
     enddo

     nLatsSolve = nLatsSolve * 2

     SolveType = SolveWithOutEquator_

     if ((HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats)) .and. &
          HighLatBoundary > Latitude(1,nLats)) then
        DoTouchNorthPole = .true.
        DoTouchSouthPole = .true.
     endif

  else

     ! We are doing at least the low latitude region

     if (.not.DoFold) then 

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              nLatsSolve=nLatsSolve+1
!           else
!              Potential(:,iLat) = 10.0
           endif
        enddo

        SolveType = SolveAcrossEquator_

        if (HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats) .and. &
             HighLatBoundary > Latitude(1,nLats)) then
           DoTouchNorthPole = .true.
           DoTouchSouthPole = .true.
        endif

     else

        ! Hmmm....  I am not sure how to define OCFLB yet...

        SolveType = SolveWithFold_

        if (HighLatBoundary < Latitude(1,nLats)+dLatitude(1,nLats) .and. &
             HighLatBoundary > Latitude(1,nLats)) then
           ! In this case we touch the pole, and we have to do the OCFLB
           ! Complicated....
           DoTouchNorthPole = .true.
           DoTouchSouthPole = .true.

           do iLat = 1, nLats
              if (Latitude(1,iLat)<=-minval(OCFLB) .or. Latitude(1,iLat)>=0) &
                   nLatsSolve = nLatsSolve + 1
           enddo

        else
           ! In this case, we have selected a high latitude boundary, and
           ! can just assume that we are solving below this latitude
           ! using a folded over conductance pattern.  This should be 
           ! relatively easy to do. Maybe.
           ! Only need to count northern hemishere, since we are doing a
           ! fold over.
           do iLat = nLats/2, nLats
              if ( Latitude(1,iLat) <= HighLatBoundary) &
                   nLatsSolve = nLatsSolve + 1
              if ( Latitude(1,iLat) >= HighLatBoundary .and. &
                   Latitude(1,iLat) < HighLatBoundary+OCFLBBuffer) then
                 r = (1-(Latitude(1,iLat) - HighLatBoundary)/OCFLBBuffer)/2
                 nPotential = Potential(:,iLat)
                 sPotential = Potential(:,nLats-iLat+1)
                 Potential(:,iLat) = nPotential*(1-r) + sPotential*r
                 Potential(:,nLats-iLat+1) = nPotential*r + sPotential*(1-r)
              endif
           enddo

        endif

     endif

  endif

  ! Don't need Ghostcells here (I think)
  nTotalSolve = nLatsSolve*nLons

  allocate( x(nTotalSolve), y(nTotalSolve), rhs(nTotalSolve), &
       b(nTotalSolve), d_I(nTotalSolve), e_I(nTotalSolve), &
       e1_I(nTotalSolve), f_I(nTotalSolve), f1_I(nTotalSolve) )

  iI = 0

  select case(SolveType)

     case(SolveAcrossEquator_)

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

     case(SolveWithOutEquator_)

        ! We have to do two separate solves - but can do this at once.

        ! Southern Hemisphere
        IsLowLat = .true.
        do iLat = 1, nLats/2
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= -LowLatBoundary) then
              IsHighLat = .false.
              if (Latitude(1,iLat+1) >= -LowLatBoundary) then
                 IsHighLat = .true.
              endif
              do iLon = 1, nLons
                 call fill
              enddo
              IsLowLat = .false.
           endif
        enddo

        ! Northern Hemisphere
        IsLowLat = .true.
        do iLat = nLats/2, nLats
           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
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

     case(SolveWithFold_)

        ! Do the same as for the Northern Hemisphere above, but
        ! LowLatBoundary is equator, so we can remove the conditional
        if (.not.DoTouchNorthPole) then
           IsLowLat = .true.
           do iLat = nLats/2, nLats
              if (Latitude(1,iLat) <= HighLatBoundary) then
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
        else
           IsLowLat = .true.
           do iLat = 1, nLats
              IsHighLat = .false.
              if (iLat == nLats) IsHighLat = .true.
              if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                   Latitude(1,iLat)>=0) then
                 do iLon = 1, nLons
                    call fill
                 enddo
              endif
              IsLowLat = .false.
           enddo
        endif

  end select

!  call ridley_solve

  Rhs = b
  if (DoPrecond) then
     ! A -> LU

     call prehepta(nTotalSolve,1,nLatsSolve,nTotalSolve,-0.5, &
          d_I,e_I,f_I,e1_I,f1_I)

     ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs

     ! rhs'=U^{-1}.L^{-1}.rhs
     call Lhepta(       nTotalSolve,1,nLatsSolve,nTotalSolve,b,d_I,e_I,e1_I)
     call Uhepta(.true.,nTotalSolve,1,nLatsSolve,nTotalSolve,b,    f_I,f1_I)
     
  end if

  if (iDebugLevel > 2) &
       write(*,*)'after precond: sum(b,abs(b),x,d,e,f,e1,f1)=',&
       sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
       sum(e1_I),sum(f1_I)

  ! Solve A'.x = rhs'
  Residual    = Tolerance
  nIteration  = MaxIteration
  if (.not.UseInitialGuess) x = 0.0
  if (iDebugLevel > 1) DoTest = .true.
  DoTest = .false.

  call gmres(matvec_RIM,b,x,UseInitialGuess,nTotalSolve,&
       MaxIteration,Residual,'abs',nIteration,iError,DoTest,iComm)
  if (iError /= 0 .and. iError /=3)then
     write(*,*)'IE_ERROR in iono_solve: gmres failed !!!'
     write(*,*)'iono_solve: iter, resid, iError=',&
          nIteration, Residual, iError
     if(iError < 0) &
          call CON_stop('IE_ERROR in iono_solve: residual did not decrease')
  end if

  iI = 0
  select case(SolveType)

     case(SolveAcrossEquator_)

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 Potential(iLon, iLat) = x(iI)
              enddo
           endif
        enddo

     case(SolveWithOutEquator_)

        ! South
        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= -LowLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 Potential(iLon, iLat) = x(iI)
              enddo
           endif
        enddo

        ! North
        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 Potential(iLon, iLat) = x(iI)
              enddo
           endif
        enddo

     case(SolveWithFold_)

        if (.not. DoTouchNorthPole) then
           ! North with mirrored south
           do iLat = nLats/2, nLats
              if (Latitude(1,iLat) <= HighLatBoundary) then
                 do iLon = 1, nLons
                    iI = iI + 1
                    Potential(iLon, iLat) = x(iI)
                    Potential(iLon, nLats-iLat+1) = x(iI)
                 enddo
              endif
           enddo
        else
           do iLat = 1, nLats
              if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                   Latitude(1,iLat)>=0) then
                 do iLon = 1, nLons
                    iI = iI + 1
                    Potential(iLon, iLat) = x(iI)
                    if ( Latitude(1,iLat) > 0 .and. &
                         Latitude(1,iLat) <= minval(OCFLB)) &
                         Potential(iLon, nLats-iLat+1) = x(iI)
                 enddo
              endif
           enddo
        endif

  end select

  ! If we include poles, then the pole solution is the average of all
  ! the cells around the pole:

!  if (DoTouchNorthPole) &
!       Potential(1:nLons,nLats) = sum(Potential(1:nLons,nLats-1))/nLons
!  if (DoTouchSouthPole) &
!       Potential(1:nLons,    1) = sum(Potential(1:nLons,      2))/nLons

  ! Need to do message passing here......

  ! Periodic Boundary Conditions:

  Potential(      0,:) = Potential(nLons,:)
  Potential(nLons+1,:) = Potential(    1,:)

  if(allocated(b)) deallocate(x, y, b, rhs, d_I, e_I, f_I, e1_I, f1_I)

  OldPotential = Potential

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

    if (IsLowLat  ) then 
       e_I(iI)  = 0.0
       if (iLat > 1 .and. .not.DoFold) &
            b(iI)=b(iI)-SolverB(iLon,iLat)*Potential(iLon,iLat-1)
    endif

    if (IsHighLat ) then
       f_I(iI)  = 0.0
       if (iLat < nLats) b(iI)=b(iI)-SolverC(iLon,iLat)*Potential(iLon,iLat+1)
    endif

!    if (IsLowLat  ) e_I(iI)  = 0.0
!    if (IsHighLat ) f_I(iI)  = 0.0
!    if (iLon == 1    ) e1_I(iI) = 0.0
!    if (iLon == nLons) f1_I(iI) = 0.0

  end subroutine fill

end subroutine solve

!-----------------------------------------------------------------------
! matvec routine
!-----------------------------------------------------------------------

subroutine matvec_RIM(x_I, y_I, n)

  use ModRIM
  use ModParamRIM
  use ModLinearsolver, ONLY: Uhepta, Lhepta
  use ModMpi
  use ModProcIE

  implicit none

  integer, intent(in) :: n          ! number of unknowns
  real, intent(in) :: x_I(n)        ! vector of unknowns
  real, intent(out):: y_I(n)        ! y = A.x

  real :: x_G(0:nLons+1, nLats), SouthPolePotential, NorthPolePotential
  real :: BufferIn(nLats), BufferOut(nLats), r

  integer :: iI, iLon, iLat, iLh

  logical :: IsHighLat, IsLowLat
  integer :: iProcFrom, iProcTo, iError

  ! MPI status variable
  integer :: iStatus_I(MPI_STATUS_SIZE)

  iI = 0

!  if (DoTouchNorthPole) Potential(:,nLats/2:nLats) = 0.0
!  if (DoTouchSouthPole) Potential(:,1:nLats/2) = 0.0

  x_G = Potential

  select case(SolveType)

     case(SolveAcrossEquator_)

        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 x_G(iLon, iLat) = x_I(iI)
              enddo
           endif
        enddo

     case(SolveWithOutEquator_)

        ! South
        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= -LowLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 x_G(iLon, iLat) = x_I(iI)
              enddo
           endif
        enddo

        ! North
        do iLat = 1, nLats
           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
                Latitude(1,iLat) <= HighLatBoundary) then
              do iLon = 1, nLons
                 iI = iI + 1
                 x_G(iLon, iLat) = x_I(iI)
              enddo
           endif
        enddo

     case(SolveWithFold_)

        ! North
        if (.not. DoTouchNorthPole) then
           do iLat = nLats/2, nLats
              if ( Latitude(1,iLat) <= HighLatBoundary) then
                 do iLon = 1, nLons
                    iI = iI + 1
                    x_G(iLon, iLat) = x_I(iI)
                 enddo
              endif
           enddo
        else
           do iLat = 1, nLats
              if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                   Latitude(1,iLat)>=0) then
                 do iLon = 1, nLons
                    iI = iI + 1
                    x_G(iLon, iLat) = x_I(iI)
                 enddo
              endif
           enddo
           ! Fill in missing southern hemisphere stuff
           do iLat = 1, nLats/2
              if ( Latitude(1,iLat) > -minval(OCFLB)) then
                 do iLon = 1, nLons
                    x_G(iLon, iLat) = x_G(iLon, nLats-iLat+1)
                 enddo
              endif
           enddo
        endif
  end select

  ! This is not really correct.

  if (DoTouchNorthPole) NorthPolePotential = sum(x_G(1:nLons,nLats-1))/nLons
  if (DoTouchSouthPole) SouthPolePotential = sum(x_G(1:nLons,2))/nLons

  ! Periodic Boundary Conditions:

  if (nProc > 1) then

     x_G(      0,:) = 0.0
     x_G(nLons+1,:) = 0.0

     ! Counterclockwise
     ! try isend and irecv
     do iProcFrom = 0, nProc-1
        iProcTo = mod(iProcFrom+1,nProc)
        if (iProc == iProcFrom) then
           BufferOut = x_G(nLons,:)
           call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo  ,1,iComm,iError)
        endif
        if (iProc == IProcTo) then
           call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                iStatus_I,iError)
           x_G(0,:) = BufferIn
        endif
     enddo

     ! Clockwise
     do iProcFrom = 0, nProc-1
        iProcTo = iProcFrom-1
        if (iProcTo == -1) iProcTo = nProc-1
        if (iProc == iProcFrom) then
           BufferOut = x_G(1,:)
           call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo  ,1,iComm,iError)
        endif
        if (iProc == IProcTo) then
           call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                iStatus_I,iError)
           x_G(nLons+1,:) = BufferIn
        endif
     enddo

  else
     x_G(      0,:) = x_G(nLons,:)
     x_G(nLons+1,:) = x_G(    1,:)
  endif

  iI = 0
  
  select case(SolveType)

     case(SolveAcrossEquator_)

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
                 if (DoTouchNorthPole .and. DoTouchSouthPole .and. &
                      iLon == 1 .and. iLat == nLats/2) then
                    iI = iI + 1
                    y_I(iI) = 0.0
                 else
                    call fill
                 endif
              enddo
              IsLowLat = .false.
           endif
        enddo

     case(SolveWithOutEquator_)

        ! Southern Hemisphere
        IsLowLat = .true.
        do iLat = 1, nLats/2
           if ( Latitude(1,iLat) >= -HighLatBoundary .and. &
                Latitude(1,iLat) <= -LowLatBoundary) then
              IsHighLat = .false.
              if (Latitude(1,iLat+1) >= -LowLatBoundary) then
                 IsHighLat = .true.
              endif
              do iLon = 1, nLons
                 call fill
              enddo
              IsLowLat = .false.
           endif
        enddo

        ! Northern Hemisphere
        IsLowLat = .true.
        do iLat = nLats/2, nLats
           if ( Latitude(1,iLat) >= LowLatBoundary .and. &
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

     case(SolveWithFold_)

        if (.not. DoTouchNorthPole) then
           ! Northern Hemisphere
           IsLowLat = .true.
           do iLat = nLats/2, nLats
              if (Latitude(1,iLat) <= HighLatBoundary) then
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
        else

           IsLowLat = .true.
           IsHighLat = .false.
           do iLat = 1, nLats
              if (iLat == nLats) IsHighLat = .true.
              if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                   Latitude(1,iLat)>=0) then
                 do iLon = 1, nLons
                    if (iLon == 1 .and. iLat == nLats/2) then
                       iI = iI + 1
                       y_I(iI) = 0.0
                    else
                       call fill
                    endif
                 enddo
              endif
              IsLowLat = .false.
           enddo

           ! Now, mirror the north to the south

           ! Need to know updated North and South potentials, so lets
           ! move the solution back into a 2D array.
           iI = 0
           do iLat = 1, nLats
              if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                   Latitude(1,iLat)>=0) then
                 do iLon = 1, nLons
                    iI = iI + 1
                    x_G(iLon,iLat) = y_I(iI)
                 enddo
              endif
           enddo

           ! For the Southern Hemisphere, make sure to fill in values
           do iLat = 1, nLats/2
              if (Latitude(1,iLat) >= -minval(OCFLB)) then
                 do iLon = 1, nLons
                    x_G(iLon,iLat) = x_G(iLon, nLats-iLat+1)
                 enddo
              endif
           enddo

!           Potential = x_G

           iI = 0
           do iLat = 1, nLats
              do iLon = 1, nLons
                 if ( Latitude(1,iLat)<=-minval(OCFLB) .or. &
                      Latitude(1,iLat)>=0) then
                    iI = iI + 1
                    ! We only care about the region within the band 
                    ! between the OCFLB and the OCFLB Buffer.
                    if (abs(Latitude(iLon,iLat)) < &
                         minval(OCFLB)+OCFLBBuffer .and. &
                         abs(Latitude(iLon,iLat)) >= minval(OCFLB)) then
                       r =  ( 1.0 - &
                            (abs(Latitude(iLon,iLat))-minval(OCFLB)) &
                            /OCFLBBuffer) / 2.0
                       ! iLh = latitude of other hemisphere
                       iLh = nLats - iLat + 1
                       y_I(iI) = (1-r) * x_G(iLon,iLat) + r * x_G(iLon,iLh)
                    endif
                 endif
              enddo
           enddo
        endif

  end select

!  do iLat = 1,nLats
!     write(*,*) 'p: ',iLat, x_G(19,iLat),x_G(20,iLat),x_G(21,iLat)
!  enddo

  ! Preconditioning: y'= U^{-1}.L^{-1}.y
  if(DoPrecond)then
     call Lhepta(       n,1,nLatsSolve,n,y_I,d_I,e_I,e1_I)
     call Uhepta(.true.,n,1,nLatsSolve,n,y_I,    f_I,f1_I)
  end if

contains

  subroutine fill

    iI = iI + 1
    
    if (iLat /= 1 .and. iLat /= nLats) then
       if (iLat == nLats/2 .and. DoFold) then
          ! The boundary condition for the fold is that the potential
          ! just below the equator is the same as the potential just
          ! above the equator.
          y_I(iI) = &
               SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
               SolverB(iLon, iLat)*x_G(iLon,  iLat+1) + &
               SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
               SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
               SolverE(iLon, iLat)*x_G(iLon+1,iLat  )
       else
          y_I(iI) = &
               SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
               SolverB(iLon, iLat)*x_G(iLon,  iLat-1) + &
               SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
               SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
               SolverE(iLon, iLat)*x_G(iLon+1,iLat  )
       endif
    else
       if (iLat == 1) then
          y_I(iI) = &
               SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
               SolverB(iLon, iLat)*SouthPolePotential + &
               SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
               SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
               SolverE(iLon, iLat)*x_G(iLon+1,iLat  )
       endif
       if (iLat == nLats) then
          y_I(iI) = &
               SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
               SolverB(iLon, iLat)*x_G(iLon,  iLat-1) + &
               SolverC(iLon, iLat)*NorthPolePotential + &
               SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
               SolverE(iLon, iLat)*x_G(iLon+1,iLat  )
       endif
    endif

    if (IsHighLat .and. iLat < nLats) then
       y_I(iI) = y_I(iI) - SolverC(iLon, iLat)*x_G(iLon,  iLat+1)
    endif
    if (IsLowLat .and. iLat > 1 .and. .not. DoFold) &
         y_I(iI) = y_I(iI) - SolverB(iLon, iLat)*x_G(iLon,  iLat-1)

  end subroutine fill

end subroutine matvec_RIM

