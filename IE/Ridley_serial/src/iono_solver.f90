!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine ionosphere_solver(iBlock, Jr, &
     SigmaThTh, SigmaThPs, SigmaPsPs, &
     dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
     dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
     Theta, Psi, dTheta, dPsi, Phi_C)

  !
  ! This subroutine solves for the ionospheric potential PHI
  ! using the field aligned currents Jr and the conductivity tensor Sigma
  ! and its various derivatives as input.
  !
  ! The idea is that the Jr current is balanced by the divergence of the
  ! horizontal currents. The horizontal currents are Sigma . E, 
  ! the height integrated conductivity Sigma is 2x2 antisymmetric matrix
  ! acting on the Theta and Phi components of the electric field.
  ! The electric field is assumed to be a potential field: E = Grad Phi.
  !
  ! We are solving 
  !
  !    Div( Sigma . Grad Phi ) = - Jr
  ! 
  ! in spherical coordinates. 
  !
  ! This leads to a penta-diagonal linear equation:
  !  
  !  C_A*Phi(i,j) + C_B*Phi(i-1,j) + C_C*Phi(i+1,j) +
  !                 C_D*Phi(i,j-1) + C_E*Phi(i,j+1) = RHS
  !
  ! To avoid division by zero at the poles, the equation is 
  ! multiplied by (sin(Theta)*Radius)**2, thus the 
  ! RHS = Jr * (sin(Theta)*Radius)**2
  !
  ! The linear elliptic PDE for the
  ! electric field potential PHI is defined on the domain 
  ! 0 < Theta < ThetaMax  for the northern hemisphere, and
  ! ThetaMin < Theta < PI for the southern hemisphere), and
  ! 0 < Psi < 2 PI.
  !
  ! The following boundary conditions are applied:
  !
  !      PHI(ThetaMax,Psi) = PHI(ThetaMin,Psi) = 0,
  !
  !      PHI(Theta,0) = PHI(Theta,2*PI).
  ! 
  ! There is no boundary at the poles, but to avoid numerical 
  ! difficulties, the cell value at the pole is replaced with the 
  ! average of the first neighbors:
  !
  !      PHI(0,Psi)  = average( PHI(dTheta, Psi) )
  !
  !      PHI(PI,Psi) = average( PHI(PI-dTheta, Psi) )
  !
  ! where dTheta is the grid resolution at the poles.
  !/

  use ModIonosphere
  use IE_ModMain, ONLY: DoCoupleUaCurrent, LatBoundaryGm, LatBoundary, &
       NameSolver, UsePreconditioner, UseInitialGuess, Tolerance, MaxIteration
  use IE_ModIo, ONLY: write_prefix, iUnitOut
  use ModLinearSolver, ONLY: gmres, bicgstab, prehepta, Uhepta, Lhepta

  implicit none

  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1

  ! INPUT ARGUMENTS:
  integer, intent(in) :: iBlock
  real, intent(in), dimension(nTheta,nPsi) ::  &
       Jr, &
       SigmaThTh, SigmaThPs, SigmaPsPs, &
       dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
       dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
       Theta, Psi
  real, intent(in) :: dTheta(nTheta), dPsi(nPsi)

  ! OUTPUT ARGUMENTS:
  real, intent(out) :: Phi_C(nTheta,nPsi) ! Potential

  ! Local variables
  real :: dTheta2(nTheta), dPsi2(nPsi)
  real :: SinTheta_I(nTheta), CosTheta_I(nTheta)
  real :: sn, cs, sn2
  real :: TermTheta2, TermTheta1, TermPsi2, TermPsi1

  real :: PhiOld_CB(nTheta, nPsi, 2) = 0.0  ! Previous solution

  integer :: i, j, iMin, iMax, iI

  ! 1D vectors for the linear solver
  real, dimension(:), allocatable :: x, y, rhs, b

  real :: PhiMin, PhiMax

  logical :: DoTest, DoTestMe

  real    :: Residual
  integer :: nIteration, iError

  external :: matvec_ionosphere

  SAVE
  character(len=*), parameter:: NameSub = 'ionosphere_solver'
  !-------------------------------------------------------------------------
  call CON_set_do_test('ionosphere',DoTest,DoTestMe)
  call timing_start('iono_solve')
  if(DoTest)write(*,*)'iono_solve starting'

  if(DoTest)write(*,*)'North, LatBoundaryGm=',north,LatBoundaryGm, LatBoundary

  if (LatBoundary > LatBoundaryGm - 5.1 * cDegToRad) &
       LatBoundary = LatBoundaryGm - 5.1 * cDegToRad
!  if (.not. DoCoupleUaCurrent) LatBoundary = LatBoundaryGm - 5.1 * cDegToRad

  ! Count the points above the latitude boundary
  nThetaUsed = count(abs(cHalfPi-Theta(1:nTheta,1)) > LatBoundary)

  ! The pole is a boundary point
  if(north)then
     iMin = 2; iMax = nThetaUsed+1  
  else
     iMin = nTheta - nThetaUsed; iMax = nTheta-1
  end if

  ! Calculate (Delta Psi)^2 (note that dPsi = 2 Delta Psi)
  dPsi2 = (dPsi/2.0)**2

  ! Calculate (Delta Theta)^2  (dTheta = 2 Delta Theta for i=1 and nTheta)
  dTheta2        = (dTheta/2.0)**2
  dTheta2(1)     = 4*dTheta2(1)
  dTheta2(nTheta)= 4*dTheta2(nTheta)

  SinTheta_I = sin(Theta(:,1))
  CosTheta_I = cos(Theta(:,1))

  if(DoTest) write(*,*) 'iono_solver: iMin, iMax=',iMin, iMax
  if(DoTest .and. north) &
       write(*,*)'north:nThetaUsed, LatBoundary, Lat(iMax)=',&
       nThetaUsed,LatBoundary*cRadToDeg,&
       abs(cHalfPi-Theta(iMax,1))*cRadToDeg
  if(DoTest .and. .not. north) &
       write(*,*)'south:nThetaUsed, LatBoundary, Lat(iMin)=',&
       nThetaUsed,LatBoundary*cRadToDeg,&
       abs(cHalfPi-Theta(iMin,1))*cRadToDeg

  nX = nPsiUsed*nThetaUsed

  allocate( x(nX), y(nX), rhs(nX), b(nX), &
       d_I(nX), e_I(nX), e1_I(nX), f_I(nX), f1_I(nX) )

  do j = 1, nPsiUsed
     do i= iMin, iMax

        sn  = SinTheta_I(i)
        cs  = CosTheta_I(i)
        sn2 = sn**2

        ! Central difference coefficients for second and first derivatives
        TermTheta2 = SigmaThTh(i,j)*sn2/dTheta2(i)
        TermTheta1 = ( SigmaThTh(i,j)*sn*cs   &
             + dSigmaThTh_dTheta(i,j)*sn2 &
             - dSigmaThPs_dPsi(i,j)*sn) / dTheta(i)
        TermPsi2 = SigmaPsPs(i,j) / dPsi2(j)
        TermPsi1 = (dSigmaThPs_dTheta(i,j)*sn + dSigmaPsPs_dPsi(i,j)) / dPsi(j)

        ! Form the complete matrix
        C_A(i,j) = -2.0 * (TermTheta2 + TermPsi2)
        C_B(i,j) =         TermTheta2 - TermTheta1
        C_C(i,j) =         TermTheta2 + TermTheta1
        C_D(i,j) =         TermPsi2   - TermPsi1
        C_E(i,j) =         TermPsi2   + TermPsi1

     enddo
  enddo

  ! Fill in the diagonal vectors
  iI = 0
  do j = 1, nPsiUsed; do i=iMin,iMax
     iI = iI + 1

     ! The right-hand-side is Jr * Radius^2 sin^2(Theta).
     b(iI)    = Jr(i,j)*(Radius*SinTheta_I(i))**2
     x(iI)    = PhiOld_CB(i,j,iBlock)
     d_I(iI)  = C_A(i,j)
     e_I(iI)  = C_B(i,j)
     f_I(iI)  = C_C(i,j)
     e1_I(iI) = C_D(i,j)
     f1_I(iI) = C_E(i,j)

     if(i == iMin)   e_I(iI)  = 0.0
     if(i == iMax)   f_I(iI)  = 0.0
     if(j == 1     ) e1_I(iI) = 0.0
     if(j == nPsiUsed) f1_I(iI) = 0.0

  end do; end do

  if(DoTest)write(*,*)'north,sum(b,abs(b),x,d,e,f,e1,f1)=',&
       north,sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
       sum(e1_I),sum(f1_I)

  DoPrecond = UsePreconditioner
  Rhs = b
  if(DoPrecond)then
     ! A -> LU
     call prehepta(nX,1,nThetaUsed,nX,-0.5,d_I,e_I,f_I,e1_I,f1_I)

     ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs

     ! rhs'=U^{-1}.L^{-1}.rhs
     call Lhepta(       nX,1,nThetaUsed,nX,b,d_I,e_I,e1_I)
     call Uhepta(.true.,nX,1,nThetaUsed,nX,b,    f_I,f1_I)

  end if

  if(DoTest)write(*,*)'after precond: north,sum(b,abs(b),x,d,e,f,e1,f1)=',&
       north,sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
       sum(e1_I),sum(f1_I)

  ! Solve A'.x = rhs'
  Residual    = Tolerance
  if(.not.UseInitialGuess) x = 0.0

  select case(NameSolver)
  case('gmres')
     nIteration = MaxIteration
     call gmres(matvec_ionosphere, b, x, UseInitialGuess, nX, MaxIteration, &
          Residual, 'abs', nIteration, iError, DoTestMe)
  case('bicgstab')
     nIteration = 3*MaxIteration
     call bicgstab(matvec_ionosphere, b, x, UseInitialGuess, nX, &
          Residual, 'abs', nIteration, iError, DoTestMe)
  case default
     call CON_stop(NameSub//': unknown NameSolver='//NameSolver)
  end select

  if(DoTest .or. (iError /= 0 .and. iError /=3) ) &
       write(*,*)'iono_solve: north, iter, resid, iError=',&
       north, nIteration, Residual, iError
  if(iError /= 0 .and. iError /=3)then
     write(*,*)'IE_ERROR in iono_solve: solver failed !!!'
     if(iError < 0) &
          call CON_stop('IE_ERROR in iono_solve: residual did not decrease')
  end if

  !Check solution:
  if(DoTest)then
     DoPrecond = .false.
     call check_solution(x, y, rhs, nX)
  end if

  ! Put solution vector into 2D array
  Phi_C(:,:) = 0.0
  iI = 0
  do j = 1, nPsiUsed; do i = iMin, iMax
     iI = iI + 1
     Phi_C(i,j) = x(iI)
  enddo; enddo
  PhiMax = maxval(x(1:iI))
  PhiMin = minval(x(1:iI))

  if (north) then
     ! Apply average condition at north pole
     Phi_C(1,:) = sum(Phi_C(2,1:nPsiUsed))/nPsiUsed
     cpcp_north = (PhiMax - PhiMin)/1000.0
     call write_prefix; write(iUnitOut,'(a,G14.6,a)') &
          "iono_solver: Northern Cross Polar Cap Potential=",&
          cpcp_north," kV"
  else
     ! Apply average condition at south pole
     Phi_C(nTheta,:) = sum(Phi_C(nTheta-1,1:nPsiUsed))/nPsiUsed
     cpcp_south = (PhiMax - PhiMin)/1000.0
     call write_prefix; write(iUnitOut,'(a,G14.6,a)') &
          "iono_solver: Southern Cross Polar Cap Potential=",&
          cpcp_south," kV"
  endif
  ! Apply periodic boundary condition in Psi direction
  Phi_C(:,nPsi) = Phi_C(:,1)

  if(allocated(b)) deallocate(x, y, b, rhs, d_I, e_I, f_I, e1_I, f1_I)

  ! Save the solution for next time
  PhiOld_CB(:,:,iBlock) = Phi_C

  call timing_stop('iono_solve')

contains
  !============================================================================
  subroutine check_solution(x_I, y_I, b_I, n)

    use IE_ModMain,   ONLY: Time_Simulation, nSolve
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file

    ! Calculate y = A.x where A is the (pentadiagonal) matrix

    integer, intent(in) :: n          ! number of unknowns
    real, intent(in) :: x_I(n)        ! vector of unknowns
    real, intent(out):: y_I(n)        ! y = A.x
    real, intent(in) :: b_I(n)        ! rhs

    integer :: iTheta, iPsi, i, iError(1)
    !-------------------------------------------------------------------------
    call matvec_ionosphere(x_I, y_I, n)

    if(north)then
       call open_file(FILE='iono_north.out')
    else
       call open_file(FILE='iono_south.out')
    end if
    write(UnitTmp_,'(a79)')'Iono Solver_var22'
    write(UnitTmp_,'(i7,es13.5,3i3)') nSolve,Time_Simulation,2,1,13
    write(UnitTmp_,'(2i4)') nThetaUsed, nPsiUsed
    write(UnitTmp_,'(es13.5)') 0.0
    write(UnitTmp_,'(a79)') &
         'Theta Psi A B C D E d e f e1 f1 Solution Rhs Error Param'

    i = 0
    iError = maxloc(abs(y_I-b_I))
    do iPsi = 1, nPsiUsed; do iTheta = iMin, iMax; i = i+1
       if(i==iError(1))write(*,*) &
            '!!! Check: max(abs(y-b)), iPsi, iTheta, nPtsTheta, north=',&
            iPsi, iTheta, nThetaUsed, abs(y_I(iError)-b_I(iError)), north
       write(UnitTmp_,'(100es18.10)') &
            cRadToDeg*Theta(iTheta,1), cRadToDeg*Psi(1,iPsi), &
            C_A(iTheta, iPsi),  C_B(iTheta, iPsi),  C_C(iTheta, iPsi), &
            C_D(iTheta, iPsi),  C_E(iTheta, iPsi),  &
            d_I(i), e_I(i), f_I(i), e1_I(i), f1_I(i), &
            y_I(i), b_I(i), y_I(i)-b_I(i)
    end do; end do
    call close_file

  end subroutine check_solution

end subroutine ionosphere_solver

!============================================================================
subroutine matvec_ionosphere(x_I, y_I, n)

  use ModIonosphere, ONLY: IONO_nPsi, IONO_nTheta, nThetaUsed, &
       north, DoPrecond, d_I, e_I, f_I, e1_I, f1_I, C_A, C_B, C_C, C_D, C_E
  use ModLinearsolver, ONLY: Uhepta, Lhepta

  implicit none

  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi, nPsiUsed=nPsi-1

  ! Calculate y = A.x where A is the (pentadiagonal) matrix

  integer, intent(in) :: n          ! number of unknowns
  real, intent(in) :: x_I(n)        ! vector of unknowns
  real, intent(out):: y_I(n)        ! y = A.x

  integer :: iTheta, iTheta2, iPsi, i, iMin, iMax
  real :: x_G(nTheta, 0:nPsi) ! 2D array with ghost cells
  !-------------------------------------------------------------------------

  if(north)then
     iMin = 2; iMax = nThetaUsed + 1
  else
     iMin = IONO_nTheta - nThetaUsed; iMax = IONO_nTheta-1
  end if

  ! Put 1D vector into 2D solution
  i = 0;
  do iPsi = 1, IONO_nPsi-1; do iTheta = iMin, iMax; i = i+1
     x_G(iTheta, iPsi) = x_I(i)
  enddo; enddo

  if(north)then
     ! Apply pole boundary condition at the north pole
     x_G(iMin-1,:) = sum(x_G(iMin,1:IONO_nPsi-1))/(IONO_nPsi-1)
     
     ! Apply 0 value below the lowest latitude
     x_G(iMax+1,:) = 0.0
  else
     ! Apply pole boundary condition at the south pole
     x_G(iMax+1,:) = sum(x_G(iMax,1:IONO_nPsi-1))/(IONO_nPsi-1)

     ! Apply 0 value above the highest latitude
     x_G(iMin-1,:) = 0.0
  end if

  ! Apply periodic boundary conditions in Psi direction
  x_G(:,nPsi) = x_G(:,1)
  x_G(:,0)    = x_G(:,nPsiUsed)

  i = 0
  do iPsi = 1, IONO_nPsi-1; do iTheta = iMin, iMax; i = i+1
     y_I(i) = &
          C_A(iTheta, iPsi)*x_G(iTheta,   iPsi)   + &
          C_B(iTheta, iPsi)*x_G(iTheta-1, iPsi)   + &
          C_C(iTheta, iPsi)*x_G(iTheta+1, iPsi)   + &
          C_D(iTheta, iPsi)*x_G(iTheta,   iPsi-1) + &
          C_E(iTheta, iPsi)*x_G(iTheta,   iPsi+1)
  end do; end do

  ! Preconditioning: y'= U^{-1}.L^{-1}.y
  if(DoPrecond)then
     call Lhepta(       n,1,nThetaUsed,n,y_I,d_I,e_I,e1_I)
     call Uhepta(.true.,n,1,nThetaUsed,n,y_I,    f_I,f1_I)
  end if

end subroutine matvec_ionosphere

