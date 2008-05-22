
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

subroutine ionosphere_solver(PHI, Jr, &
     SigmaThTh, SigmaThPs, SigmaPsPs, &
     dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
     dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
     Theta, Psi, Radius, nTheta, nPsi, &
     dTheta, dPsi,                                                           &
     ncycle, iboundary)

  use ModIonosphere

  USE MSR_Module
  USE Matrix_Arithmetic_Module
  USE ILU_Module
  USE GMRES_Module
  USE Iteration_Defaults_Module
  use ModMpi
  implicit none


  integer :: nTheta, nPsi, ncycle, iboundary
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, Jr, Source, &
                  SigmaThTh, SigmaThPs, SigmaPsPs, &
                  dSigmaThTh_dTheta, dSigmaThPs_dTheta, dSigmaPsPs_dTheta, &
                  dSigmaThTh_dPsi, dSigmaThPs_dPsi, dSigmaPsPs_dPsi, &
                  Theta, Psi, &
                  kappa_Theta2, kappa_Theta1, kappa_Psi2, kappa_Psi1, &
                  C_A, C_B, C_C, C_D, C_E, sn, cs, sn2, cs2

  real, dimension(1:IONO_nTheta) :: dTheta, dTheta2
  real, dimension(1:IONO_nPsi)   :: dPsi, dPsi2

  integer :: nThetaC, nPsiC, nThetaF, nPsiF
  integer :: i, j, i2, j2, ii, ind_i, ind_j, k, jj, n, npts, npts_Theta
  integer :: iStart, iEnd
  integer, save :: saved_npts_Theta=-1

  TYPE(MSR)                             :: S
  TYPE(LU_CSR_MSR)                      :: PC
  REAL(prec), DIMENSION(:), allocatable :: x, b
  REAL(prec)                            :: tmp
  REAL(prec), DIMENSION(1:5)            :: row
  INTEGER, DIMENSION(1:5)               :: cols
  REAL(prec), DIMENSION(1:1)            :: row1
  INTEGER, DIMENSION(1:1)               :: cols1

  REAL :: ave, LowLatBoundary, min_pot, max_pot, cpcp, HighLatBoundary
  LOGICAL :: north
  
  real*8 :: StartTimer,EndTimer
  real*8 :: StartTimer2,EndTimer2

  LowLatBoundary = 35.0 * IONO_PI/180.0

  north = .false.
!  if (Theta(1,1) < 2.0*IONO_Theta_0) north = .true.
  if (Theta(1,1) < IONO_PI/4.0) north = .true.

  StartTimer = mpi_wtime()

  HighLatBoundary = 0.0

  if (north) then
     do i=1,nTheta
        if (Jr(i,nPsi/4) /= 0.0) LowLatBoundary = abs(IONO_PI/2.0-Theta(i,1))
        if (Jr(i,nPsi/4) /= 0.0 .and. HighLatBoundary == 0.0) &
             HighLatBoundary = abs(IONO_PI/2.0-Theta(i,1))
     enddo
  else
     do i=nTheta,1,-1
        if (Jr(i,nPsi/4) /= 0.0) LowLatBoundary = abs(IONO_PI/2.0-Theta(i,1))
        if (Jr(i,nPsi/4) /= 0.0 .and. HighLatBoundary == 0.0) &
             HighLatBoundary = abs(IONO_PI/2.0-Theta(i,1))
     enddo
  endif

!  LowLatBoundary = LowLatBoundary - 5.0*IONO_PI/180.0
!  HighLatBoundary = HighLatBoundary + 5.0*IONO_PI/180.0
!!! Using set latitudes for the low and high boundary location sometimes:
  LowLatBoundary = 30.0*IONO_PI/180.0
  HighLatBoundary = 67.0*IONO_PI/180.0

  if (north) then
     do i=1,nTheta
        if (abs(IONO_PI/2.0-Theta(i,1)) > HighLatBoundary) iEnd = i
        if (abs(IONO_PI/2.0-Theta(i,1)) > LowLatBoundary)  iStart = i
     enddo
  else
     do i=nTheta,1,-1
        if (abs(IONO_PI/2.0-Theta(i,1)) > HighLatBoundary) iEnd = i
        if (abs(IONO_PI/2.0-Theta(i,1)) > LowLatBoundary)  iStart = i
     enddo
  endif

!!!  write(*,*) 'iStart, iEnd (1) : ', iStart, iEnd

  if (iStart > iEnd) then
     i = iStart
     iStart = iEnd
     iEnd = i
  endif

  write(*,*)'LowLatBoundary=', &
       LowLatBoundary/IONO_PI*180.0, HighLatBoundary/IONO_PI*180.0
!!!  write(*,*) 'iStart, iEnd : ', iStart, iEnd

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
        Source(i,j) = Jr(i,j)*Radius*Radius* &
                      sin(Theta(i,j))*sin(Theta(i,j))
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

!  do j=1,nPsi
!
!! Add the Pole:
!
!     npts = npts + 5
!     if (j.eq.1) npts_Theta = npts_Theta + 1
!
!! Add the Equator:
!
!     npts = npts + 5
!     if (j.eq.1) npts_Theta = npts_Theta + 1
!
!! Add the Central Points
!
!     do i=2,nTheta-1
!        if (abs(IONO_PI/2.0-Theta(i,j)).gt.LowLatBoundary) then
!           npts = npts + 5
!           if (j.eq.1) npts_Theta = npts_Theta + 1
!        endif
!     enddo
!
!  enddo


  !\
  ! Start Test Code
  !/
  npts_Theta = iEnd - iStart + 1 + 2
  !\
  ! End Test Code
  !/

  if(npts_Theta /= saved_npts_Theta)then
     if(allocated(b)) deallocate(b)
     if(allocated(x)) deallocate(x)
!     if(allocated_matrix(S))  call Deallocate_Matrix(S)
!     if(allocated_matrix(PC)) call Deallocate_Matrix(PC)
  end if
  saved_npts_Theta = npts_Theta

  if (.NOT.allocated(b)) allocate( b(nPsi*npts_Theta) )
  if (.NOT.allocated(x)) allocate( x(nPsi*npts_Theta) )

  npts = 5*npts_Theta*nPsi

  CALL Nullify_Matrix ( S )
  CALL Allocate_Matrix ( S, nPsi*npts_Theta, npts)
  if (.NOT.allocated_matrix(S)) then
     write(6,*) "Error allocating MSR structure S in ionosphere solver"
     stop
  endif

!     write(6,*) '=> ALLOCATING MATRIX PC in ionosphere'

  CALL Nullify_Matrix ( PC )
  CALL Allocate_Matrix (PC, nPsi*npts_Theta,npts,npts)
  if (.NOT.allocated_matrix(PC)) then
     write(6,*) "Error allocating PC in ionosphere solver"
     stop
  endif

!  write(*,*) 'Done with the allocation'

!
! We basically are solving the equation:
!  
!  C_A*phi(i,j) + C_B*phi(i-1,j) + C_C*phi(i+1,j) +
!                 C_D*phi(i,j-1) + C_E*phi(i,j-1) = S
!
! We are using a sparse matrix solver. So we are going to shove
!  A, B, C, D, and E into the sparse matrix at the right locations,
!  and shove S into the RHS of the matrix equation.
!

!  write(*,*) 'Entering nPsi and nTheta loops for coefficient setup'

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

!  write(*,*) 'Done with initial setup. Entering nPsi loop for more setup.'

  do j = 1, nPsi

     i = iStart - 1
     ind_i = (j-1)*npts_Theta + (i - iStart + 2)

     b(ind_i) = Phi(i,j)

!     if (north) then 
!        b(ind_i) = -5000.0 * sin(Psi(1,j))
!     else
!        b(ind_i) = 0.0
!     endif

     x(ind_i) = 0.0

     cols1(1) = ind_i
     row1(1) = 1.0

     call setrow(S,ind_i,cols1,row1)

     do i = iStart, iEnd

        ii = i - iStart + 2

        ind_i = (j-1)*npts_Theta + ii

        b(ind_i) = Source(i,j)
        x(ind_i) = 0.0

        j2 = j-1
        i2 = ii
        if (j2 < 1) j2 = nPsi-1
        ind_j = (j2-1)*npts_Theta + i2
        cols(1) = ind_j
        row(1) = C_D(i,j)

        j2 = j
        i2 = ii-1
        ind_j = (j2-1)*npts_Theta + i2
        cols(2) = ind_j
        row(2) = C_B(i,j)

        j2 = j
        i2 = ii
        ind_j = (j2-1)*npts_Theta + i2
        cols(3) = ind_j
        row(3) = C_A(i,j)
        
        j2 = j
        i2 = ii+1
        ind_j = (j2-1)*npts_Theta + i2
        cols(4) = ind_j
        row(4) = C_C(i,j)

        j2 = j+1
        i2 = ii
        if (j2 > nPsi) j2 = 2
        ind_j = (j2-1)*npts_Theta + i2
        cols(5) = ind_j
        row(5) = C_E(i,j)
        
        call setrow(S,ind_i,cols,row)

     enddo

     i = iEnd + 1 
     ii = i - iStart + 2
     ind_i = (j-1)*npts_Theta + ii

     b(ind_i) = Phi(i,j)
!     if (north) then 
!        b(ind_i) = 0.0
!     else
!        b(ind_i) = -5000.0 * sin(Psi(1,j))
!     endif

     x(ind_i) = 0.0

     cols1(1) = ind_i
     row1(1) = 1.0
     call setrow(S,ind_i,cols1,row1)

  enddo

!  write(*,*) 'Done with second setup.'

!  write(6,*) "=> Filling MSR Array"

!  if (north) then   
!
!     do j = 1, nPsi
!
!! we are at the pole
!
!        i = 1
!        ind_i = (j-1)*npts_Theta + i
!        b(ind_i) = phi(i,j)
!        x(ind_i) = 0.0
!
!        j2 = j-1
!        i2 = i+1
!        if (j2 < 1) j2 = nPsi-1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(1) = ind_j
!        row(1) = C_D(i,j)
!        
!        j2 = j+(nPsi-1)/2
!        if (j2.gt.nPsi) j2 = j2 - nPsi
!        i2 = i+1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(2) = ind_j
!        row(2) = C_B(i,j)
!
!        i2 = i
!        j2 = j
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(3) = ind_j
!        row(3) = C_A(i,j)
!        
!        j2 = j
!        i2 = i+1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(4) = ind_j
!        row(4) = C_C(i,j)
!
!        j2 = j+1
!        i2 = i+1
!        if (j2 > nPsi) j2 = 2
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(5) = ind_j
!        row(5) = C_E(i,j)
!
!        CALL SetRow(S,ind_i,cols,row)
!
!        i = 2
!
!        do while (abs(IONO_PI/2.0-Theta(i,j)).gt.LowLatBoundary)
!
!           ind_i = (j-1)*npts_Theta + i
!
!           b(ind_i) = phi(i,j)
!           x(ind_i) = 0.0
!
!           j2 = j-1
!           i2 = i
!           if (j2 < 1) j2 = nPsi-1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(1) = ind_j
!           row(1) = C_D(i,j)
!
!           j2 = j
!           i2 = i-1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(2) = ind_j
!           row(2) = C_B(i,j)
!
!           j2 = j
!           i2 = i
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(3) = ind_j
!           row(3) = C_A(i,j)
!        
!           j2 = j
!           i2 = i+1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(4) = ind_j
!           row(4) = C_C(i,j)
!
!           j2 = j+1
!           i2 = i
!           if (j2 > nPsi) j2 = 2
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(5) = ind_j
!           row(5) = C_E(i,j)
!        
!           call setrow(S,ind_i,cols,row)
!
!           i = i + 1
!
!        enddo
!
!        ind_i = (j-1)*npts_Theta + i
!        b(ind_i) = 0.0
!        x(ind_i) = 0.0
!
!        cols1(1) = ind_i
!        row1(1) = 1.0
!        call setrow(S,ind_i,cols1,row1)
!
!     enddo
!
!  else
!
!     n = nTheta - npts_Theta + 1
!
!     do j = 1, nPsi
!
!        i = n
!
!        ind_i = (j-1)*npts_Theta + (i-n+1)
!
!        b(ind_i) = phi(i,j)
!        x(ind_i) = 0.0
!        cols1(1) = ind_i
!        row1(1) = 1.0
!        call setrow(S,ind_i,cols1,row1)
!
!        do i = n+1, nTheta-1
!
!           ind_i = (j-1)*npts_Theta + (i-n+1)
!
!           b(ind_i) = phi(i,j)
!
!           x(ind_i) = 0.0
!
!           j2 = j-1
!           i2 = i - n + 1
!           if (j2 < 1) j2 = nPsi-1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(1) = ind_j
!           row(1) = C_D(i,j)
!
!           j2 = j
!           i2 = i-1 - n + 1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(2) = ind_j
!           row(2) = C_B(i,j)
!
!           j2 = j
!           i2 = i - n + 1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(3) = ind_j
!           row(3) = C_A(i,j)
!        
!           j2 = j
!           i2 = i+1 - n + 1
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(4) = ind_j
!           row(4) = C_C(i,j)
!
!           j2 = j+1
!           i2 = i - n + 1
!           if (j2 > nPsi) j2 = 2
!           ind_j = (j2-1)*npts_Theta + i2
!           cols(5) = ind_j
!           row(5) = C_E(i,j)
!
!           call setrow(S,ind_i,cols,row)
!
!        enddo
!
!        i = nTheta
!        ind_i = (j-1)*npts_Theta + i - n + 1
!        b(ind_i) = phi(i,j)
!
!        j2 = j-1
!        i2 = i-1 - n + 1
!        if (j2 < 1) j2 = nPsi-1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(1) = ind_j
!        row(1) = C_D(i,j)
!
!        j2 = j
!        i2 = i-1 - n + 1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(2) = ind_j
!        row(2) = C_B(i,j)
!
!        i2 = i - n + 1
!        j2 = j
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(3) = ind_j
!        row(3) = C_A(i,j)
!        
!        j2 = j+(nPsi-1)/2
!        if (j2.gt.nPsi) j2 = j2 - nPsi
!        i2 = i-1 - n + 1
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(4) = ind_j
!        row(4) = C_C(i,j)
!
!        j2 = j+1
!        i2 = i-1 - n + 1
!        if (j2 > nPsi) j2 = 2
!        ind_j = (j2-1)*npts_Theta + i2
!        cols(5) = ind_j
!        row(5) = C_E(i,j)
!
!        call setrow(S,ind_i,cols,row)
!
!     enddo
!
!  endif

  if (.NOT.is_ok(S)) then
     write(6,*) "Error in MSR structure S in ionosphere solver"
     stop
  endif

!  write(6,*) '=> Going into Symbolic Cancellation.'

  StartTimer2 = mpi_wtime()

  CALL Symbolic_Cancellation_ilu (1, S, PC)

  if (.NOT.is_ok(PC)) then
     write(6,*) "Error in structure PC in ionosphere solver"
     stop
  endif

!  write(6,*) '=> Going into ILU.'

  CALL ILU(S,PC)

  EndTimer2 = mpi_wtime()

!  write(*,*) "-------------------ILU Time : ",EndTimer2-StartTimer2


  if (.NOT.is_ok(S)) then
     write(6,*) "Error in MSR structure S in ionosphere solver"
     stop
  endif
  if (.NOT.is_ok(PC)) then
     write(6,*) "Error in MSR structure PC in ionosphere solver"
     stop
  endif

!  write(6,*) '=> Going into ionospheric GMRES Solver.'

  CALL GMRES_MSR(10,x,S,b,PC)

!  write(6,*) '=> Deallocating Matrix S and PC'

  if(allocated_matrix(S))  call Deallocate_Matrix(S)
  if(allocated_matrix(PC)) call Deallocate_Matrix(PC)

  min_pot = 9999.0
  max_pot = -9999.0

  if (north) then

     do j = 1, nPsi
        do i = iStart-1, iEnd+1
           ii = i - iStart + 2
           ind_i = (j-1)*npts_Theta + ii
           phi(i,j) = x(ind_i)
           if (x(ind_i) > max_pot) max_pot = x(ind_i)
           if (x(ind_i) < min_pot) min_pot = x(ind_i)
        enddo
     enddo

  else

     do j = 1, nPsi
        do i = iStart-1, iEnd+1
           ii = i - iStart + 2
           ind_i = (j-1)*npts_Theta + ii
           phi(i,j) = x(ind_i)
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

  write(6,*) "=> Cross Polar Cap Potential : ",cpcp/1000.0," kV"
  write(6,*) "=> Potential at pole : ",ave/1000.0," kV"

  i = nTheta
  if (north) i = 1
  do j = 1, nPsi
     phi(i,j) = ave
  enddo

  if(allocated(b)) deallocate(b)
  if(allocated(x)) deallocate(x)

  EndTimer = mpi_wtime()

!  write(*,*) "Timing in iono_solver : ",EndTimer-StartTimer

end subroutine ionosphere_solver
