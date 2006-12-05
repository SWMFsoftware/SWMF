!^CFG COPYRIGHT UM

!*************************************************************************
!
! SOLVER Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! ionospheric_solver
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_solver(PHI, &
     SigmaThTh, SigmaThPs, SigmaPsPs, &
     dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
     dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
     Theta, Psi, Radius, nTheta, nPsi, &
     dTheta, dPsi,                     &
     ncycle, iboundary)

  !\
  ! This subroutine solves for the ionospheric potential PHI
  ! using the field aligned currents and the conductivities as input.
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
  use IE_ModMain, ONLY: DoCoupleUaCurrent, LatBoundary
  use IE_ModIo, ONLY: write_prefix, iUnitOut
  use ModLinearSolver, ONLY: gmres, prehepta, Uhepta, Lhepta

  implicit none

  integer :: nTheta, nPsi, ncycle, iboundary
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
       PHI, &
       SigmaThTh, SigmaThPs, SigmaPsPs, &
       dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
       dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
       Theta, Psi, &
       kappa_Theta2, kappa_Theta1, kappa_Psi2, kappa_Psi1, &
       sn, cs, sn2, cs2

  real, dimension(1:IONO_nTheta) :: dTheta, dTheta2
  real, dimension(1:IONO_nPsi)   :: dPsi, dPsi2

  integer :: nThetaC, nPsiC, nThetaF, nPsiF
  integer :: i, j, i1, i2, j2, ind_i, ind_j, k, jj, n, npts
  integer :: saved_npts_Theta=-1

  real, dimension(:), allocatable :: x, y, rhs, b

  real :: ave, min_pot, max_pot, cpcp

  logical :: oktest, oktest_me, UseOldSolver

  real    :: Tolerance
  integer :: iter, info

  external :: matvec_ionosphere

  SAVE
  !-------------------------------------------------------------------------

  call CON_set_do_test('ionosphere',oktest,oktest_me)
  call timing_start('iono_solve')
  if(oktest)write(*,*)'iono_solve starting'

  north = .false.

  if (Theta(1,1) < cPi/4.0) north = .true.

  if(oktest)write(*,*)'North=',north

  if (.not. DoCoupleUaCurrent) then
     ! Set the latitude boundary 5 degrees below lowest non-zero FAC

     ! Set a value just in case there are no currents at all
     LatBoundary = 45.0 * cDegToRad

     if (north) then
        do i=1,nTheta
           if (PHI(i,nPsi/4) /= 0.0) LatBoundary = abs(cHalfPi-Theta(i,1))
        enddo
     else
        do i=nTheta,1,-1
           if (PHI(i,nPsi/4) /= 0.0) LatBoundary = abs(cHalfPi-Theta(i,1))
        enddo
     endif
     LatBoundary = max(LatBoundary - 5.0*cDegToRad,0.0)
  endif

!!! write(*,*)'PHI(:,nPsi/4), LatBoundary=',PHI(:,nPsi/4), LatBoundary !!!

  if(oktest)write(*,*)'sum(abs(PHI),LatBoundary=',sum(abs(PHI)),LatBoundary

  do j = 1, nPsi
     dPsi2(j) = (dPsi(j)/2.0)**2
  enddo

  do i = 2, nTheta-1
     dTheta2(i) = (dTheta(i)/2.0)**2
  enddo
  dTheta2(1)     = dTheta(1)**2
  dTheta2(nTheta)= dTheta(nTheta)**2

  ! Multiply the right-hand-side source terms
  ! (in effect the field-aligned current) by
  ! Radius^2 sin^2(Theta).

  do j = 1, nPsi
     do i = 1, nTheta
        if (north) then
           if (abs(cHalfPi-Theta(i,j)) > LatBoundary+5.0*cDegToRad) then
              PHI(i,j) = PHI(i,j)*(Radius*sin(Theta(i,j)))**2
           else
              PHI(i,j) = 0.0
           endif
        else
           if (abs(Theta(i,j)-cHalfPi) > LatBoundary+5.0*cDegToRad) then
              PHI(i,j) = PHI(i,j)*(Radius*sin(Theta(i,j)))**2
           else
              PHI(i,j) = 0.0
           endif
        endif
     end do
  end do

  ! The Allocate Matrix routine creates a sparse matrix structure of the
  ! size specified.  The over all matrix size is "nPsi*nTheta^2", but we are
  ! only going to have (at most) "5*(nPsi*nTheta)" values.
  !

  !     write(6,*) '=> ALLOCATING MATRIX S in ionosphere'

  ! Start with 0:

  npts = 0
  npts_Theta = 0

  do j=1,nPsi

     ! Add the Pole:

     npts = npts + 5
     if (j.eq.1) npts_Theta = npts_Theta + 1

     ! Add the Equator:

     npts = npts + 5
     if (j.eq.1) npts_Theta = npts_Theta + 1

     ! Add the Central Points

     do i=2,nTheta-1
        if (abs(cHalfPi-Theta(i,j)) > LatBoundary) then
           npts = npts + 5
           if (j.eq.1) npts_Theta = npts_Theta + 1
        endif
     enddo

  enddo

  npts = 5*npts_Theta*nPsi

  if(oktest)write(*,*)'npts_Theta=',npts_Theta

  if(npts_Theta /= saved_npts_Theta)then
     if(allocated(b)) deallocate(b)
     if(allocated(x)) deallocate(x)
     if(allocated(y)) deallocate(y)
  end if
  saved_npts_Theta = npts_Theta

  nX = nPsi*npts_Theta

  if (.NOT.allocated(x)) allocate( x(nX), y(nX), rhs(nX), b(nX), &
       d_I(nX), e_I(nX), e1_I(nX), f_I(nX), f1_I(nX) )

  !
  ! We basically are solving the equation:
  !  
  !  C_A*phi(i,j) + C_B*phi(i-1,j) + C_C*phi(i+1,j) +
  !                 C_D*phi(i,j-1) + C_E*phi(i,j+1) = S
  !
  ! When we come into the subroutine, Clinton has stored
  ! S in phi, so that is why we are setting S = phi(i,j)!
  !
  ! We are using a sparse matrix solver. So we are going to shove
  !  A, B, C, D, and E into the sparse matrix at the right locations,
  !  and shove S into the RHS of the matrix equation.
  !

  do j = 1, nPsi
     do i= 1, nTheta

        sn(i,j)  = sin(Theta(i,j))
        cs(i,j)  = cos(Theta(i,j)) 
        sn2(i,j) = sn(i,j)*sn(i,j) 
        cs2(i,j) = cs(i,j)*cs(i,j) 

        kappa_Theta2(i,j) = SigmaThTh(i,j)*sn2(i,j) 
        kappa_Theta1(i,j) = SigmaThTh(i,j)*sn(i,j)*cs(i,j)   &
             + dSigmaThTh_dTheta(i,j)*sn2(i,j)               &
             - dSigmaThPs_dPsi(i,j)*sn(i,j) 
        kappa_Psi2(i,j)   = SigmaPsPs(i,j)
        kappa_Psi1(i,j)   = dSigmaThPs_dTheta(i,j)*sn(i,j)   &
             + dSigmaPsPs_dPsi(i,j)

        C_A(i,j) = -2.0 * (&
             kappa_Theta2(i,j)/dTheta2(i) + kappa_Psi2(i,j)/dPsi2(j))

        C_B(i,j) = &
             kappa_Theta2(i,j)/dTheta2(i) - kappa_Theta1(i,j)/dTheta(i)

        C_C(i,j) = &
             kappa_Theta2(i,j)/dTheta2(i) + kappa_Theta1(i,j)/dTheta(i)

        C_D(i,j) = &
             kappa_Psi2(i,j)/dPsi2(j) - kappa_Psi1(i,j)/dPsi(j)

        C_E(i,j) = &
             kappa_Psi2(i,j)/dPsi2(j) + kappa_Psi1(i,j)/dPsi(j)

     enddo
  enddo

  !  write(6,*) "=> Filling MSR Array"

  if (north) then   

     do j = 1, nPsi

        ! we are at the pole

        i = 1
        ind_i = (j-1)*npts_Theta + i
        b(ind_i) = phi(i,j)
        x(ind_i) = 0.0
        d_I(ind_i)  = C_A(i,j)
        e_I(ind_i)  = C_B(i,j)
        f_I(ind_i)  = C_C(i,j)
        e1_I(ind_i) = C_D(i,j)
        f1_I(ind_i) = C_E(i,j)


        i = 2

        do while (abs(cHalfPi-Theta(i,j)) > LatBoundary)

           ind_i = (j-1)*npts_Theta + i
           b(ind_i) = phi(i,j)
           x(ind_i) = 0.0
           d_I(ind_i)  = C_A(i,j)
           e_I(ind_i)  = C_B(i,j)
           f_I(ind_i)  = C_C(i,j)
           e1_I(ind_i) = C_D(i,j)
           f1_I(ind_i) = C_E(i,j)


           i = i + 1

        enddo

        ind_i = (j-1)*npts_Theta + i
        b(ind_i) = 0.0
        x(ind_i) = 0.0
        d_I(ind_i)  = C_A(i,j)
        e_I(ind_i)  = C_B(i,j)
        f_I(ind_i)  = C_C(i,j)
        e1_I(ind_i) = C_D(i,j)
        f1_I(ind_i) = C_E(i,j)

     enddo

  else

     n = nTheta - npts_Theta + 1

     do j = 1, nPsi

        i = n

        ind_i = (j-1)*npts_Theta + (i-n+1)

        b(ind_i) = phi(i,j)
        x(ind_i) = 0.0
        d_I(ind_i)  = C_A(i,j)
        e_I(ind_i)  = C_B(i,j)
        f_I(ind_i)  = C_C(i,j)
        e1_I(ind_i) = C_D(i,j)
        f1_I(ind_i) = C_E(i,j)

        do i = n+1, nTheta-1

           ind_i = (j-1)*npts_Theta + (i-n+1)

           b(ind_i) = phi(i,j)
           x(ind_i) = 0.0
           d_I(ind_i)  = C_A(i,j)
           e_I(ind_i)  = C_B(i,j)
           f_I(ind_i)  = C_C(i,j)
           e1_I(ind_i) = C_D(i,j)
           f1_I(ind_i) = C_E(i,j)

        enddo

        i = nTheta
        ind_i = (j-1)*npts_Theta + i - n + 1
        b(ind_i) = phi(i,j)
        x(ind_i) = 0.0
        d_I(ind_i)  = C_A(i,j)
        e_I(ind_i)  = C_B(i,j)
        f_I(ind_i)  = C_C(i,j)
        e1_I(ind_i) = C_D(i,j)
        f1_I(ind_i) = C_E(i,j)

     enddo

  endif

  iter = 100
  DoPrecond = .true.
  rhs = b
  if(DoPrecond)then
     ! A -> LU
     call prehepta(nX,1,npts_theta,nX,-0.5,d_I,e_I,f_I,e1_I,f1_I)
     ! rhs'=L^{-1}.rhs
     call Lhepta(nX,1,npts_theta,nX,b,d_I,e_I,e1_I)
  end if
  ! Solve A'.x' = rhs'
  Tolerance = 1.0e-4
  call gmres(matvec_ionosphere,b,x,.false.,nPsi*npts_Theta,&
       100,Tolerance,'rel',Iter,info,.true.)
  write(*,*)'!!! gmres: iter, info=',iter, Tolerance, info
  if(DoPrecond)then
     ! x = U^{-1}.x'
     call Uhepta(.true.,nX,1,npts_theta,nX,x,f_I,f1_I)
  end if

  !Check solution:
  DoPrecond = .false.
  call check_solution(x, y, rhs, nPsi*npts_Theta)

  phi(:,:) = 0.0

  min_pot = 9999.0
  max_pot = -9999.0

  if (north) then

     do j = 1, nPsi
        do i = 1, npts_Theta

           ind_i = (j-1)*npts_Theta + i
           phi(i,j) = x(ind_i)
           if (x(ind_i) > max_pot) max_pot = x(ind_i)
           if (x(ind_i) < min_pot) min_pot = x(ind_i)

        enddo
     enddo

  else

     do j = 1, nPsi
        do i = 1, npts_Theta

           ind_i = (j-1)*npts_Theta + i
           phi(i + (nTheta - npts_Theta),j) = x(ind_i)
           if (x(ind_i) > max_pot) max_pot = x(ind_i)
           if (x(ind_i) < min_pot) min_pot = x(ind_i)

        enddo
     enddo

  endif

  cpcp = max_pot - min_pot

  i = nTheta-1
  if (north) i = 2
  ave = 0.0
  do j = 1, nPsi
     ave = ave + phi(i,j)
  enddo
  ave = ave/nPsi

  call write_prefix
  if (north) then
     cpcp_north = cpcp/1000.0
     write(iUnitOut,*) &
          "iono_solver: Northern Cross Polar Cap Potential=",&
          cpcp_north," kV"
  else
     cpcp_south = cpcp/1000.0
     write(iUnitOut,*) &
          "iono_solver: Southern Cross Polar Cap Potential=",&
          cpcp_south," kV"
  endif

  i = nTheta
  if (north) i = 1
  do j = 1, nPsi
     phi(i,j) = ave
  enddo

  if(allocated(b)) deallocate(x, y, b, rhs, d_I, e_I, f_I, e1_I, f1_I)

  call timing_stop('iono_solve')

contains
  !============================================================================
  subroutine check_solution(x_I, y_I, b_I, n)

    use ModIoUnit, ONLY: UnitTmp_

    ! Calculate y = A.x where A is the (pentadiagonal) matrix

    integer, intent(in) :: n          ! number of unknowns
    real, intent(in) :: x_I(n)        ! vector of unknowns
    real, intent(out):: y_I(n)        ! y = A.x
    real, intent(in) :: b_I(n)        ! rhs

    integer :: iTheta, iPsi, i, iError(1)
    !-------------------------------------------------------------------------
    call matvec_ionosphere(x_I, y_I, n)

    if(north)then
       open(UnitTmp_, FILE='iono_north.out')
    else
       open(UnitTmp_, FILE='iono_south.out')
    end if
    write(UnitTmp_,*)'Iono Solver_var22'
    write(UnitTmp_,*)'0.0 0 2 1 3'
    write(UnitTmp_,*) npts_Theta, nPsi
    write(UnitTmp_,*)'0.0'
    write(UnitTmp_,*)'Theta Psi Solution Rhs Error Param'

    i = 0
    iError = maxloc(abs(y_I-b_I))
    do iPsi = 1, nPsi; do iTheta = 1, npts_Theta; i = i+1
       if(i==iError(1))write(*,*) &
            '!!! Check: max(abs(y-b)), iPsi, iTheta, nPtsTheta, north=',&
            iPsi, iTheta, npts_Theta, abs(y_I(iError)-b_I(iError)), north
       write(UnitTmp_,*)iTheta, iPsi, y_I(i), b_I(i), y_I(i)-b_I(i)
    end do; end do
    close(UnitTmp_)

  end subroutine check_solution

end subroutine ionosphere_solver

!============================================================================
subroutine matvec_ionosphere(x_I, y_I, n)

  use ModIonosphere, ONLY: IONO_nPsi, IONO_nTheta, npts_Theta, &
       north, DoPrecond, d_I, e_I, f_I, e1_I, f1_I, C_A, C_B, C_C, C_D, C_E
  use ModLinearsolver, ONLY: Uhepta, Lhepta

  implicit none
  ! Calculate y = A.x where A is the (pentadiagonal) matrix

  integer, intent(in) :: n          ! number of unknowns
  real, intent(in) :: x_I(n)        ! vector of unknowns
  real, intent(out):: y_I(n)        ! y = A.x

  integer :: iTheta, iTheta2, iPsi, i
  real :: x_G(0:IONO_nTheta+1, 0:IONO_nPsi+1) ! 2D array with ghost cells
  !-------------------------------------------------------------------------

  ! Preconditioning: x' = U^{-1}.x stored in y_I
  y_I = x_I

  if(DoPrecond) call Uhepta(.true.,n,1,npts_Theta,n,y_I,f_I,f1_I)

  ! Put 1D vector into 2D solution
  i = 0;
  do iPsi = 1, IONO_nPsi; do iTheta = 1, npts_Theta; i = i+1
     x_G(iTheta, iPsi) = y_I(i)
  enddo; enddo

  ! Apply periodic boundary conditions in Psi direction
  x_G(:,0)           = x_G(:,IONO_nPsi-1)
  x_G(:,IONO_nPsi+1) = x_G(:,2)

  if(north)then
     ! Apply pole boundary condition at the north pole
     do iPsi = 1, IONO_nPsi
        x_G(0, iPsi) = x_G(2, mod(iPsi+IONO_nPsi/2, IONO_nPsi))
     end do
     ! Apply 0 value at the lowest latitude
     x_G(npts_Theta+1,:) = 0.0
  else
     ! Apply 0 value at the highest latitude
     x_G(0,:) = 0.0
     ! Apply pole boundary condition at the south pole
     do iPsi = 1, IONO_nPsi
        x_G(npts_Theta+1, iPsi) = &
             x_G(npts_Theta-1, mod(iPsi+IONO_nPsi/2, IONO_nPsi))
     end do
  end if

  i = 0
  if(north)then
     do iPsi = 1, IONO_nPsi; do iTheta = 1, npts_Theta; i = i+1
        y_I(i) = &
             C_A(iTheta, iPsi)*x_G(iTheta,   iPsi)   + &
             C_B(iTheta, iPsi)*x_G(iTheta-1, iPsi)   + &
             C_C(iTheta, iPsi)*x_G(iTheta+1, iPsi)   + &
             C_D(iTheta, iPsi)*x_G(iTheta,   iPsi-1) + &
             C_E(iTheta, iPsi)*x_G(iTheta,   iPsi+1)

     end do; end do
  else
     do iPsi = 1, IONO_nPsi; do iTheta = 1, npts_Theta; i = i+1
        iTheta2 = iTheta + IONO_nTheta - npts_Theta
        y_I(i) = &
             C_A(iTheta2, iPsi)*x_G(iTheta,   iPsi)   + &
             C_B(iTheta2, iPsi)*x_G(iTheta-1, iPsi)   + &
             C_C(iTheta2, iPsi)*x_G(iTheta+1, iPsi)   + &
             C_D(iTheta2, iPsi)*x_G(iTheta,   iPsi-1) + &
             C_E(iTheta2, iPsi)*x_G(iTheta,   iPsi+1)

     end do; end do
  end if

  ! Preconditioning: y'=L^{-1}y
  if(DoPrecond)call Lhepta(n,1,npts_theta,n,y_I,d_I,e_I,e1_I) 

end subroutine matvec_ionosphere

