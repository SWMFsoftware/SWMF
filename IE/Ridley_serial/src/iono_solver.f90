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
     dTheta, dPsi,                                                           &
     ncycle, iboundary)

  !\
  ! This subroutine applies a full linear multigrid solution
  ! algorithm to the problem of determining the solution of the
  ! boundary value problem (BVP) for the linear elliptic equation
  ! that governs planetary ionospheric currents on a spherical
  ! surface of radius Radius.  The linear elliptic PDE for the
  ! electric field potential PHI is defined on the domain 
  ! 0 < Theta < PI/2 (this is for the northern hemisphere, 
  ! PI/2 < Theta < PI for southern hemisphere), 0 < Psi < 2 PI, 
  ! and subject to the Dirichlet boundary data
  !
  !      PHI(0,Psi) = PHI(PI,Psi) = 0,
  !
  !      PHI(PI/2,Psi) = 0,
  !
  !      PHI(Theta,0) = PHI(Theta,2*PI).
  ! 
  ! A red-black Gauss-Seidel scheme is used as the smoothing
  ! operator, bilinear interpolation is used for the prolongation
  ! operator, and half-weighting is used for the restriction
  ! operator.
  !/

  use ModIonosphere
  use IE_ModIo, ONLY: write_prefix, iUnitOut

  USE MSR_Module
  USE Matrix_Arithmetic_Module
  USE ILU_Module
  USE GMRES_Module
  USE Iteration_Defaults_Module

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
                  C_A, C_B, C_C, C_D, C_E, sn, cs, sn2, cs2

  real, dimension(1:IONO_nTheta) :: dTheta, dTheta2
  real, dimension(1:IONO_nPsi)   :: dPsi, dPsi2

  integer :: nThetaC, nPsiC, nThetaF, nPsiF
  integer :: i, j, i2, j2, ind_i, ind_j, k, jj, n, npts, npts_Theta

  TYPE(MSR)                             :: S
  TYPE(LU_CSR_MSR)                      :: PC
  REAL(prec), DIMENSION(:), allocatable :: x, b
  REAL(prec)                            :: tmp
  REAL(prec), DIMENSION(1:5)            :: row
  INTEGER, DIMENSION(1:5)               :: cols
  REAL(prec), DIMENSION(1:1)            :: row1
  INTEGER, DIMENSION(1:1)               :: cols1

  REAL :: ave, lat_boundary, min_pot, max_pot, cpcp
  LOGICAL :: north
  LOGICAL, PARAMETER :: IsManualLatBoundary = .false.
  
  logical :: oktest, oktest_me
  
  call CON_set_do_test('ionosphere',oktest,oktest_me)
  call timing_start('iono_solve')
  if(oktest)write(*,*)'iono_solve starting'

  lat_boundary = 30.0 * IONO_PI/180.0

  north = .false.

  if (Theta(1,1) < IONO_PI/4.0) north = .true.

  if(oktest)write(*,*)'North=',north

  if (IsManualLatBoundary) then
     lat_boundary = 60.0 * IONO_PI/180.0
  else
     if (north) then
        do i=1,nTheta
           if (PHI(i,nPsi/4).ne.0.0) lat_boundary = abs(IONO_PI/2.0-Theta(i,1))
        enddo
     else
        do i=nTheta,1,-1
           if (PHI(i,nPsi/4).ne.0.0) lat_boundary = abs(IONO_PI/2.0-Theta(i,1))
        enddo
     endif
     lat_boundary = lat_boundary - 5.0*IONO_PI/180.0
  endif

  !!! write(*,*)'PHI(:,nPsi/4), lat_boundary=',PHI(:,nPsi/4), lat_boundary !!!

  if(oktest)write(*,*)'sum(abs(PHI),lat_boundary=',sum(abs(PHI)),lat_boundary

  do j = 1, nPsi
     dPsi2(j) = (dPsi(j)/2.0)*(dPsi(j)/2.0)
  enddo

  do i = 2, nTheta-1
     dTheta2(i) = (dTheta(i)/2.0)*(dTheta(i)/2.0)
  enddo
  dTheta2(1)     = dTheta(1)*dTheta(1)
  dTheta2(nTheta)= dTheta(nTheta)*dTheta(nTheta)

!  dTheta_l=(Theta(nTheta,1)-Theta(1,1))/real(nTheta-1)
!  dPsi_l=(Psi(1,nPsi)-Psi(1,1))/real(nPsi-1)
!  dTheta_l2=dTheta_l*dTheta_l
!  dPsi_l2=dPsi_l*dPsi_l
!  dTheta_lX2=2.0*dTheta_l
!  dPsi_lX2=2.0*dPsi_l
!  dd=dPsi_l
!  dd2=dPsi_l2

  ! Multiply the right-hand-side source terms
  ! (in effect the field-aligned current) by
  ! Radius^2 sin^2(Theta).
  
  do j = 1, nPsi
     do i = 1, nTheta
        if (north) then
           if (abs(IONO_PI/2.0-Theta(i,j)).gt.lat_boundary+5.0*IONO_PI/180.0) then
              PHI(i,j) = PHI(i,j)*Radius*Radius* &
                   sin(Theta(i,j))*sin(Theta(i,j))
           else
              PHI(i,j) = 0.0
           endif
        else
           if (abs(Theta(i,j)-IONO_PI/2.0).gt.lat_boundary+5.0*IONO_PI/180.0) then
              PHI(i,j) = PHI(i,j)*Radius*Radius* &
                   sin(Theta(i,j))*sin(Theta(i,j))
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
        if (abs(IONO_PI/2.0-Theta(i,j)).gt.lat_boundary) then
           npts = npts + 5
           if (j.eq.1) npts_Theta = npts_Theta + 1
        endif
     enddo

  enddo

  if(oktest)write(*,*)'npts_Theta=',npts_Theta

  allocate(b(npts_Theta*nPsi))
  allocate(x(npts_Theta*nPsi))

  npts = 5*npts_Theta*nPsi

  CALL Nullify_Matrix ( S )
  CALL Allocate_Matrix ( S, nPsi*npts_Theta, npts)
  if (.NOT.allocated_matrix(S)) &
     call CON_stop("Error allocating MSR structure S in ionosphere solver")


!     write(6,*) '=> ALLOCATING MATRIX PC in ionosphere'

  CALL Nullify_Matrix ( PC )
  CALL Allocate_Matrix (PC, nPsi*npts_Theta,npts,npts)
  if (.NOT.allocated_matrix(PC)) &
     call CON_stop("Error allocating PC in ionosphere solver")
!
! We basically are solving the equation:
!  
!  C_A*phi(i,j) + C_B*phi(i-1,j) + C_C*phi(i+1,j) +
!                 C_D*phi(i,j-1) + C_E*phi(i,j-1) = S
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
        kappa_Theta1(i,j) = SigmaThTh(i,j)*sn(i,j)*cs(i,j)             &
             + dSigmaThTh_dTheta(i,j)*sn2(i,j)               &
             - dSigmaThPs_dPsi(i,j)*sn(i,j) 
        kappa_Psi2(i,j)   = SigmaPsPs(i,j)
        kappa_Psi1(i,j)   = dSigmaThPs_dTheta(i,j)*sn(i,j)        &
             + dSigmaPsPs_dPsi(i,j)
           
        C_A(i,j) = -2.0 * (kappa_Theta2(i,j)/dTheta2(i) +              &
             kappa_Psi2(i,j)/dPsi2(j))
        
        C_B(i,j) = kappa_Theta2(i,j)/dTheta2(i) -                      &
             kappa_Theta1(i,j)/dTheta(i)

        C_C(i,j) = kappa_Theta2(i,j)/dTheta2(i) +                      &
             kappa_Theta1(i,j)/dTheta(i)

        C_D(i,j) = kappa_Psi2(i,j)/dPsi2(j) -                          &
             kappa_Psi1(i,j)/dPsi(j)

        C_E(i,j) = kappa_Psi2(i,j)/dPsi2(j) +                          &
             kappa_Psi1(i,j)/dPsi(j)

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

        j2 = j-1
        i2 = i+1
        if (j2 < 1) j2 = nPsi-1
        ind_j = (j2-1)*npts_Theta + i2
        cols(1) = ind_j
        row(1) = C_D(i,j)
        
        j2 = j+(nPsi-1)/2
        if (j2.gt.nPsi) j2 = j2 - nPsi
        i2 = i+1
        ind_j = (j2-1)*npts_Theta + i2
        cols(2) = ind_j
        row(2) = C_B(i,j)

        i2 = i
        j2 = j
        ind_j = (j2-1)*npts_Theta + i2
        cols(3) = ind_j
        row(3) = C_A(i,j)
        
        j2 = j
        i2 = i+1
        ind_j = (j2-1)*npts_Theta + i2
        cols(4) = ind_j
        row(4) = C_C(i,j)

        j2 = j+1
        i2 = i+1
        if (j2 > nPsi) j2 = 2
        ind_j = (j2-1)*npts_Theta + i2
        cols(5) = ind_j
        row(5) = C_E(i,j)

        CALL SetRow(S,ind_i,cols,row)

        i = 2

        do while (abs(IONO_PI/2.0-Theta(i,j)).gt.lat_boundary)

           ind_i = (j-1)*npts_Theta + i

           b(ind_i) = phi(i,j)
           x(ind_i) = 0.0

           j2 = j-1
           i2 = i
           if (j2 < 1) j2 = nPsi-1
           ind_j = (j2-1)*npts_Theta + i2
           cols(1) = ind_j
           row(1) = C_D(i,j)

           j2 = j
           i2 = i-1
           ind_j = (j2-1)*npts_Theta + i2
           cols(2) = ind_j
           row(2) = C_B(i,j)

           j2 = j
           i2 = i
           ind_j = (j2-1)*npts_Theta + i2
           cols(3) = ind_j
           row(3) = C_A(i,j)
        
           j2 = j
           i2 = i+1
           ind_j = (j2-1)*npts_Theta + i2
           cols(4) = ind_j
           row(4) = C_C(i,j)

           j2 = j+1
           i2 = i
           if (j2 > nPsi) j2 = 2
           ind_j = (j2-1)*npts_Theta + i2
           cols(5) = ind_j
           row(5) = C_E(i,j)
        
           call setrow(S,ind_i,cols,row)

           i = i + 1

        enddo

        ind_i = (j-1)*npts_Theta + i
        b(ind_i) = 0.0
        x(ind_i) = 0.0

        cols1(1) = ind_i
        row1(1) = 1.0
        call setrow(S,ind_i,cols1,row1)

     enddo

  else

     n = nTheta - npts_Theta + 1

     do j = 1, nPsi

        i = n

        ind_i = (j-1)*npts_Theta + (i-n+1)

        b(ind_i) = phi(i,j)
        x(ind_i) = 0.0
        cols1(1) = ind_i
        row1(1) = 1.0
        call setrow(S,ind_i,cols1,row1)

        do i = n+1, nTheta-1

           ind_i = (j-1)*npts_Theta + (i-n+1)

           b(ind_i) = phi(i,j)
           x(ind_i) = 0.0

           j2 = j-1
           i2 = i - n + 1
           if (j2 < 1) j2 = nPsi-1
           ind_j = (j2-1)*npts_Theta + i2
           cols(1) = ind_j
           row(1) = C_D(i,j)

           j2 = j
           i2 = i-1 - n + 1
           ind_j = (j2-1)*npts_Theta + i2
           cols(2) = ind_j
           row(2) = C_B(i,j)

           j2 = j
           i2 = i - n + 1
           ind_j = (j2-1)*npts_Theta + i2
           cols(3) = ind_j
           row(3) = C_A(i,j)
        
           j2 = j
           i2 = i+1 - n + 1
           ind_j = (j2-1)*npts_Theta + i2
           cols(4) = ind_j
           row(4) = C_C(i,j)

           j2 = j+1
           i2 = i - n + 1
           if (j2 > nPsi) j2 = 2
           ind_j = (j2-1)*npts_Theta + i2
           cols(5) = ind_j
           row(5) = C_E(i,j)

           call setrow(S,ind_i,cols,row)

        enddo

        i = nTheta
        ind_i = (j-1)*npts_Theta + i - n + 1
        b(ind_i) = phi(i,j)

        j2 = j-1
        i2 = i-1 - n + 1
        if (j2 < 1) j2 = nPsi-1
        ind_j = (j2-1)*npts_Theta + i2
        cols(1) = ind_j
        row(1) = C_D(i,j)

        j2 = j
        i2 = i-1 - n + 1
        ind_j = (j2-1)*npts_Theta + i2
        cols(2) = ind_j
        row(2) = C_B(i,j)

        i2 = i - n + 1
        j2 = j
        ind_j = (j2-1)*npts_Theta + i2
        cols(3) = ind_j
        row(3) = C_A(i,j)
        
        j2 = j+(nPsi-1)/2
        if (j2.gt.nPsi) j2 = j2 - nPsi
        i2 = i-1 - n + 1
        ind_j = (j2-1)*npts_Theta + i2
        cols(4) = ind_j
        row(4) = C_C(i,j)

        j2 = j+1
        i2 = i-1 - n + 1
        if (j2 > nPsi) j2 = 2
        ind_j = (j2-1)*npts_Theta + i2
        cols(5) = ind_j
        row(5) = C_E(i,j)

        call setrow(S,ind_i,cols,row)

     enddo

  endif

  if(oktest)write(*,*)'MSR array filled'

  if (.NOT.is_ok(S)) &
     call CON_stop("Error in MSR structure S in ionosphere solver")

!  write(6,*) '=> Going into Symbolic Cancellation.'

  CALL Symbolic_Cancellation_ilu (1, S, PC)

  if(oktest)write(*,*)'Symbolic_Cancellation_ilu done'

  if (.NOT.is_ok(PC)) &
     call CON_stop("Error in structure PC in ionosphere solver")

!  write(6,*) '=> Going into ILU.'

  CALL ILU(S,PC)

  if(oktest)write(*,*)'ILU done'

  if (.NOT.is_ok(S)) &
     call CON_stop("Error in MSR structure S in ionosphere solver")
  if (.NOT.is_ok(PC)) &
     call CON_stop("Error in MSR structure PC in ionosphere solver")

!  write(6,*) '=> Going into ionospheric GMRES Solver.'

  CALL GMRES_MSR(10,x,S,b,PC)

  if(oktest)write(*,*)'GMRES_MSR done'

!  write(6,*) '=> Deallocating Matrix S and PC'

  CALL Deallocate_Matrix(S)
  CALL Deallocate_Matrix(PC)

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

  deallocate(b)
  deallocate(x)

  call timing_stop('iono_solve')

end subroutine ionosphere_solver
