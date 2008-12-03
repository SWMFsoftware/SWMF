subroutine solve

  use ModRIM
  use ModParamRIM
  use ModNumConst,     only: cPi
  use ModProcIE,       only: iComm, iProc, nProc
  use ModLinearSolver, only: gmres, prehepta, Uhepta, Lhepta
  use ModMpi

  implicit none

  real, dimension(0:nLons+1,nLats) :: sinTheta
  real, dimension(:), allocatable :: x, y, rhs, b
  real, dimension(0:nLons+1) :: nPotential, sPotential

  integer :: iLat, iLon, nTotalSolve, iI
  logical :: IsLowLat, IsHighLat, DoTest, DoIdealTest=.false.

  real :: Residual, r
  integer :: nIteration, iError

  real :: BufferIn(nLats), BufferOut(nLats)
  integer :: iProcFrom, iProcTo

  ! MPI status variable
  integer :: iStatus_I(MPI_STATUS_SIZE)

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

  if (DoFold) then

     if (iDebugLevel > 1) write(*,*) "RIM==> Using Ridley Solver"
     call ridley_solve

  else

     if (iDebugLevel > 1) write(*,*) "RIM==> Using Linear GMRES Solver"

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

       Rhs = b
       if (UsePreconditioner) then
          ! A -> LU

          call prehepta(nTotalSolve,1,nLatsSolve,nTotalSolve,-0.5, &
               d_I,e_I,f_I,e1_I,f1_I)

          ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs
          
          ! rhs'=U^{-1}.L^{-1}.rhs
          call Lhepta(       nTotalSolve,1,nLatsSolve,&
               nTotalSolve,b,d_I,e_I,e1_I)
          call Uhepta(.true.,nTotalSolve,1,nLatsSolve,&
               nTotalSolve,b,    f_I,f1_I)
          
       end if

       if (iDebugLevel > 2) &
            write(*,*)'RIM===> after precond: sum(b,abs(b),x,d,e,f,e1,f1)=',&
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
       if (iError /= 0 .and. iError /=3 .and. iDebugLevel > -1)then
          write(*,*)'IE_ERROR in iono_solve: gmres failed !!!'
          write(*,'(a,i5,1p,e12.3,i4)')&
               'iono_solve: iter, resid (kV), iError=',&
               nIteration, Residual/1000.0, iError
          if(iError < 0) &
               call CON_stop('IE_ERROR in iono_solve: residual not decreasing')
       else
          if (iDebugLevel > 2) &
               write(*,"(a,i5,1p,e12.3)") " RIM===>nIter : ", &
               nIteration, Residual
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

  endif

  ! If we include poles, then the pole solution is the average of all
  ! the cells around the pole:

  if (DoTouchNorthPole) &
       Potential(1:nLons,nLats) = sum(Potential(1:nLons,nLats-1))/nLons
  if (DoTouchSouthPole) &
       Potential(1:nLons,    1) = sum(Potential(1:nLons,      2))/nLons

  ! Periodic Boundary Conditions:

  if (nProc > 1) then

     Potential(      0,:) = 0.0
     Potential(nLons+1,:) = 0.0

     ! Counterclockwise
     ! try isend and irecv
     do iProcFrom = 0, nProc-1
        iProcTo = mod(iProcFrom+1,nProc)
        if (iProc == iProcFrom) then
           BufferOut = Potential(nLons,:)
           call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo  ,1,iComm,iError)
        endif
        if (iProc == IProcTo) then
           call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                iStatus_I,iError)
           Potential(0,:) = BufferIn
        endif
     enddo
          
     ! Clockwise
     do iProcFrom = 0, nProc-1
        iProcTo = iProcFrom-1
        if (iProcTo == -1) iProcTo = nProc-1
        if (iProc == iProcFrom) then
           BufferOut = Potential(1,:)
           call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo  ,1,iComm,iError)
        endif
        if (iProc == IProcTo) then
           call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                iStatus_I,iError)
           Potential(nLons+1,:) = BufferIn
        endif
     enddo

  else
     Potential(      0,:) = Potential(nLons,:)
     Potential(nLons+1,:) = Potential(    1,:)
  endif

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

    if (iLat == nLats/2 .and. iProc == 0) then
       b(iI)=b(iI)-SolverB(iLon,iLat)*Potential(iLon,iLat-1)
    endif

  end subroutine fill

  !----------------------

  subroutine ridley_solve

    integer :: nIters
    real :: Old(0:nLons+1,nLats), LocalVar, NorthPotential, SouthPotential
    logical :: IsDone

    nIters = 0

    if (.not.UseInitialGuess) Potential = 0.0

    IsDone = .false.

    do while (.not.IsDone)

       Old = Potential

       do iLat = 2, nLats-1
          do iLon = 1, nLons

             ! Start at the Southern Pole and move northwards.
             ! Skip over the solve for the north potential on the closed
             ! field-lines and average between north and south with a 
             ! weighting for the slush region.

             if ( Latitude(iLon,iLat) < -LowLatBoundary .or. &
                  Latitude(iLon,iLat) > OCFLB(1,iLon)) then
                Potential(iLon,iLat) =  &
                     (Jr(iLon,iLat)*(Radius*sinTheta(iLon,iLat))**2 - &
                     (SolverB(iLon,iLat)*Old(iLon,iLat+1) + &
                      SolverC(iLon,iLat)*Old(iLon,iLat-1) + &
                      SolverD(iLon,iLat)*Old(iLon-1,iLat) + &
                      SolverE(iLon,iLat)*Old(iLon+1,iLat)) ) / &
                      SolverA(iLon,iLat)
             else
                if (OCFLB(1,iLon)-abs(Latitude(iLon,iLat)) < OCFLBBuffer) then
                   r = 1.0 - 0.5*&
                        (OCFLB(1,iLon)-abs(Latitude(iLon,iLat)))/OCFLBBuffer
                else
                   r = 0.5
                endif
                Potential(iLon,iLat) = &
                     (1.0 - r) * Potential(iLon,nLats-iLat+1) + &
                     (      r) * Old(iLon,iLat)
             endif 
          enddo
       enddo

       Potential(1,nLats/2) = 0.0
       Potential(1,nLats/2+1) = 0.0

       if (DoTouchNorthPole) then
          Potential(1:nLons,nLats) = sum(Potential(1:nLons,nLats-1))/nLons
          LocalVar = Potential(1,nLats)
          NorthPotential = 0.0
          call MPI_REDUCE(LocalVar, NorthPotential, 1, MPI_REAL, &
               MPI_SUM, 0, iComm, iError)
          NorthPotential = NorthPotential/nProc
          call MPI_Bcast(NorthPotential,1,MPI_Real,0,iComm,iError)
          Potential(1:nLons,nLats) = NorthPotential
       endif

       if (DoTouchSouthPole) then
          Potential(1:nLons,1) = sum(Potential(1:nLons,2))/nLons
          LocalVar = Potential(1,1)
          SouthPotential = 0.0
          call MPI_REDUCE(LocalVar, SouthPotential, 1, MPI_REAL, &
               MPI_SUM, 0, iComm, iError)
          SouthPotential = SouthPotential/nProc
          call MPI_Bcast(SouthPotential,1,MPI_Real,0,iComm,iError)
          Potential(1:nLons,1) = SouthPotential
       endif

       ! Periodic Boundary Conditions:

       if (nProc > 1) then

          Potential(      0,:) = 0.0
          Potential(nLons+1,:) = 0.0

          ! Counterclockwise
          ! try isend and irecv
          do iProcFrom = 0, nProc-1
             iProcTo = mod(iProcFrom+1,nProc)
             if (iProc == iProcFrom) then
                BufferOut = Potential(nLons,:)
                call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo,1,iComm,iError)
             endif
             if (iProc == IProcTo) then
                call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                     iStatus_I,iError)
                Potential(0,:) = BufferIn
             endif
          enddo
          
          ! Clockwise
          do iProcFrom = 0, nProc-1
             iProcTo = iProcFrom-1
             if (iProcTo == -1) iProcTo = nProc-1
             if (iProc == iProcFrom) then
                BufferOut = Potential(1,:)
                call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo,1,iComm,iError)
             endif
             if (iProc == IProcTo) then
                call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                     iStatus_I,iError)
                Potential(nLons+1,:) = BufferIn
             endif
          enddo

       else
          Potential(      0,:) = Potential(nLons,:)
          Potential(nLons+1,:) = Potential(    1,:)
       endif

       Old = Potential

       do iLat = nLats-1, 2, -1

          do iLon = 1, nLons

             if ( Latitude(iLon,iLat) > LowLatBoundary .or. &
                  Latitude(iLon,iLat) < -OCFLB(2,iLon)) then
                Potential(iLon,iLat) =  &
                     (Jr(iLon,iLat)*(Radius*sinTheta(iLon,iLat))**2 - &
                     (SolverB(iLon,iLat)*Old(iLon,iLat+1) + &
                      SolverC(iLon,iLat)*Old(iLon,iLat-1) + &
                      SolverD(iLon,iLat)*Old(iLon-1,iLat) + &
                      SolverE(iLon,iLat)*Old(iLon+1,iLat)) ) / &
                      SolverA(iLon,iLat)
             else
                if (OCFLB(2,iLon)-abs(Latitude(iLon,iLat)) < OCFLBBuffer) then
                   r = 1.0 - 0.5* &
                        (OCFLB(2,iLon)-abs(Latitude(iLon,iLat)))/OCFLBBuffer
                else
                   r = 0.5
                endif
                Potential(iLon,iLat) = &
                     (1.0 - r) * Potential(iLon,nLats-iLat+1) + &
                     (      r) * Old(iLon,iLat)
             endif 
          enddo

       enddo

!       do iLat = 2, nLats-1
!          do iLon = 1, nLons
!             if ( abs(Latitude(iLon,iLat)) < minval(OCFLB)) &
!                  Potential(iLon,iLat) = &
!                  Potential(iLon,nLats-iLat+1) 
!          enddo
!       enddo

       Potential(1,nLats/2) = 0.0
       Potential(1,nLats/2+1) = 0.0

       if (DoTouchNorthPole) then
          Potential(1:nLons,nLats) = sum(Potential(1:nLons,nLats-1))/nLons
          LocalVar = Potential(1,nLats)
          NorthPotential = 0.0
          call MPI_REDUCE(LocalVar, NorthPotential, 1, MPI_REAL, &
               MPI_SUM, 0, iComm, iError)
          NorthPotential = NorthPotential/nProc
          call MPI_Bcast(NorthPotential,1,MPI_Real,0,iComm,iError)
          Potential(1:nLons,nLats) = NorthPotential
       endif

       if (DoTouchSouthPole) then
          Potential(1:nLons,1) = sum(Potential(1:nLons,2))/nLons
          LocalVar = Potential(1,1)
          SouthPotential = 0.0
          call MPI_REDUCE(LocalVar, SouthPotential, 1, MPI_REAL, &
               MPI_SUM, 0, iComm, iError)
          SouthPotential = SouthPotential/nProc
          call MPI_Bcast(SouthPotential,1,MPI_Real,0,iComm,iError)
          Potential(1:nLons,1) = SouthPotential
       endif

       ! Periodic Boundary Conditions:

       if (nProc > 1) then

          Potential(      0,:) = 0.0
          Potential(nLons+1,:) = 0.0

          ! Counterclockwise
          ! try isend and irecv
          do iProcFrom = 0, nProc-1
             iProcTo = mod(iProcFrom+1,nProc)
             if (iProc == iProcFrom) then
                BufferOut = Potential(nLons,:)
                call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo,1,iComm,iError)
             endif
             if (iProc == IProcTo) then
                call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                     iStatus_I,iError)
                Potential(0,:) = BufferIn
             endif
          enddo
          
          ! Clockwise
          do iProcFrom = 0, nProc-1
             iProcTo = iProcFrom-1
             if (iProcTo == -1) iProcTo = nProc-1
             if (iProc == iProcFrom) then
                BufferOut = Potential(1,:)
                call MPI_send(BufferOut,nLats,MPI_REAL,iProcTo,1,iComm,iError)
             endif
             if (iProc == IProcTo) then
                call MPI_recv(BufferIn ,nLats,MPI_REAL,iProcFrom,1,iComm, &
                     iStatus_I,iError)
                Potential(nLons+1,:) = BufferIn
             endif
          enddo

       else
          Potential(      0,:) = Potential(nLons,:)
          Potential(nLons+1,:) = Potential(    1,:)
       endif

       Residual = sum((Old-Potential)**2)

       nIters = nIters + 1

       LocalVar = Residual
       call MPI_REDUCE(localVar, Residual, 1, MPI_REAL, MPI_SUM, &
            0, iComm, iError)
       call MPI_Bcast(Residual,1,MPI_Real,0,iComm,iError)
       Residual = sqrt(Residual)

       if (Residual < Tolerance) IsDone = .true.
       if (nIters >= MaxIteration) IsDone = .true.

    enddo

    if (iDebugLevel > 1) &
         write(*,*) "RIM==> Final Residual : ", nIters, Residual

  end subroutine ridley_solve

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
              if (iProc == 0 .and. iLat == nLats/2 + 1) x_G(iLon, iLat) = 0.0
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

  ! This is not really correct, since it only works on a single processor...

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
                 call fill
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
                 if (abs(Latitude(iLon,iLat)-LowLatBoundary) <= OCFLBBuffer) then
                    y_i(iI) = y_i(iI) * &
                         (abs(Latitude(iLon,iLat)-LowLatBoundary)/OCFLBBuffer)**2
                 endif
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
                    call fill
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

  ! Preconditioning: y'= U^{-1}.L^{-1}.y
  if(UsePreconditioner)then
     call Lhepta(       n,1,nLatsSolve,n,y_I,d_I,e_I,e1_I)
     call Uhepta(.true.,n,1,nLatsSolve,n,y_I,    f_I,f1_I)
  end if

contains

  subroutine fill

    real :: factor1, factor2
    iI = iI + 1
    
    if (iLat /= 1 .and. iLat /= nLats) then
       if (iLat == nLats/2+1 .and. DoFold) then
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
!          y_I(iI) = &
!               SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
!               SolverB(iLon, iLat)*x_G(iLon,  iLat-1) + &
!               SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
!               SolverD(iLon, iLat)*x_G(iLon-1,iLat  ) + &
!               SolverE(iLon, iLat)*x_G(iLon+1,iLat  )

          if (abs(Latitude(iLon,iLat))-LowLatBoundary <= OCFLBBuffer) then

             factor1 = sin(((abs(Latitude(iLon,iLat+1))-LowLatBoundary)/OCFLBBuffer)*3.141592/2)
             factor2 = sin(((abs(Latitude(iLon,iLat-1))-LowLatBoundary)/OCFLBBuffer)*3.141592/2)

             y_I(iI) = &
                  SolverA(iLon, iLat)*x_G(iLon,  iLat  ) + &
                  factor2*SolverB(iLon, iLat)*x_G(iLon,  iLat-1) + &
                  factor1*SolverC(iLon, iLat)*x_G(iLon,  iLat+1) + &
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

    ! Enforce 0 potential at the equator at midnight
!    if (iLat == nLats/2+1 .and. iLon == 1 .and. iProc == 0) then
    if (iLat == nLats/2+1 .and. iProc == 0) then
!       y_I(iI) = 0.0
         y_I(iI) = y_I(iI) - SolverB(iLon, iLat)*x_G(iLon,  iLat-1)
!       y_I(iI) = y_I(iI) - SolverA(iLon, iLat)*x_G(iLon,  iLat)
    endif

  end subroutine fill

end subroutine matvec_RIM

