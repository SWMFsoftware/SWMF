!^CFG FILE IMPLICIT
subroutine gmres(matvec,Rhs,Sol,IsInit,n,nKrylov,Tol,TypeStop,Iter,info,&
     oktest)

  !*************************************************************
  ! This code was initially written by Youcef Saad
  ! then revised by Henk A. van der Vorst and Mike Botchev (oct. 1996)
  ! Rewritten into F90 and parallelized for the BATSRUS code (may 2002) 
  ! by Gabor Toth
  !*************************************************************
  ! gmres algorithm . simple version .  (may 23, 1985)

  implicit none

  interface
     subroutine matvec(a,b,n)
       ! Calculate b = M.a where M is the matrix
       integer, intent(in) :: n
       real, intent(in) ::  a(n)
       real, intent(out) :: b(n)
     end subroutine matvec
  end interface
  !        subroutine for matrix vector multiplication 

  integer, intent(in) :: n
  !        on input:  number of unknowns.

  integer, intent(in) :: nKrylov
  !        on input:  size of krylov subspace

  real, intent(in) :: Rhs(n)
  !        on input:  right hand side

  real, intent(inout) :: Sol(n)
  !        on input : initial guess if IsInit is .true.
  !        on output: solution vector

  logical, intent(in):: IsInit
  !        on input : true  if Sol contains initial guess
  !                   false if initial guess is zero

  real, intent(inout) :: Tol
  !        on input : required (relative) 2-norm or maximum norm of residual
  !        on output: achieved (relative) 2-norm or maximum norm of residual
  ! eps    == tolerance for stopping criterion. process is stopped
  !           as soon as ( ||.|| is the euclidean norm):
  !           || current residual||/||initial residual|| <= eps
  !           on OUTPUT: actual achieved norm residual (if iabs.ne.0)
  !           or achieved relative residual norm reduction

  character (len=3), intent(in) :: TypeStop
  !      Determine stopping criterion (||.|| denotes the 2-norm):
  !      typestop='rel'    -- relative stopping crit.: ||res|| <= Tol*||res0||
  !      typestop='abs'    -- absolute stopping crit.: ||res|| <= Tol
  !      typestop='max'    -- maximum  stopping crit.: max(abs(res)) <= Tol

  integer, intent(inout) :: Iter
  !       on input:  maximum number of iterations to be performed.
  !       on output: actual  number of iterations done.

  integer, intent(out)   :: info
  !       Gives reason for returning:
  !     abs(info)=  0 - solution found satisfying given tolerance.
  !                 2 - no convergence within maximum number of iterations.
  !                 3 - initial guess satisfies the stopping criterion.
  !    sign(info)=  + - residual decreased
  !                 - - residual did not reduce

  logical, intent(in)    :: oktest
  !       write debug info if true

  ! Local variables
  integer :: i,i1,its,j,k,k1
  real :: coeff,Tol1,epsmac,gam,ro,ro0,t,tmp

  ! Automatic array for vectors in Krylov subspace
  real :: Krylov_II(n,nKrylov+2)

  ! Automatic arrays (Hessenberg matrix and some vectors)
  real :: hh(nKrylov+1,nKrylov),c(nKrylov),s(nKrylov),rs(nKrylov+1)
  !-----------------------------------------------------------------------
  if(oktest)write(*,*)'GMRES tol,iter:',Tol,Iter

  if(range(1.0)>100)then
     epsmac=0.0000000000000001
  else
     epsmac=0.00000001
  endif

  its = 0
  !-------------------------------------------------------------
  ! **  outer loop starts here..
  !-------------- compute initial residual vector --------------

  RESTARTLOOP: do
     !
     !           Krylov_II(1):=A*Sol
     !
     if(IsInit.or.its>0)then
        call matvec(Sol,Krylov_II,n)
        Krylov_II(:,1)=Rhs - Krylov_II(:,1)
     else
        ! Save a matvec when starting from zero initial condition
        Krylov_II(:,1)=Rhs
     endif
     !-------------------------------------------------------------
     ro = sqrt( dot_product_mpi(Krylov_II(:,1),Krylov_II(:,1)) )
     if (ro == 0.0) then
        if(its == 0)then
           info=3
        else
           info = 0
        endif
        Tol = ro
        Iter = its 
        RETURN
     end if

     ! set Tol1 for stopping criterion
     if (its == 0) then
        ro0 = ro
        if(oktest) print *,'initial rnrm:',ro0
        if (TypeStop=='abs') then
           Tol1=Tol
           if (ro <= Tol1) then ! quit if accurate enough
              info = 3
              Tol  = ro
              Iter = its
              if(oktest) print *,'GMRES: nothing to do. info = ',info
              return
           end if
        else
           Tol1=Tol*ro
        end if
     end if

     coeff = 1.0 / ro
     Krylov_II(:,1)=coeff*Krylov_II(:,1)

     ! initialize 1-st term  of rhs of hessenberg system
     rs(1) = ro
     i = 0
     KRYLOVLOOP: do
        i=i+1
        its = its + 1
        i1 = i + 1
        !
        !           Krylov_II(i1):=A*Krylov_II(i)
        !
        call matvec(Krylov_II(1,i),Krylov_II(1,i1),n) 
        !-----------------------------------------
        !  modified gram - schmidt...
        !-----------------------------------------
        do j=1, i
           t = dot_product_mpi(Krylov_II(:,j),Krylov_II(:,i1))
           hh(j,i) = t
           Krylov_II(:,i1) = Krylov_II(:,i1) - t*Krylov_II(:,j)
        end do
        t = sqrt( dot_product_mpi(Krylov_II(:,i1),Krylov_II(:,i1)) )

        hh(i1,i) = t
        if (t /= 0.0)then
           t = 1.0 / t
           Krylov_II(:,i1) = t*Krylov_II(:,i1)
        endif
        !--------done with modified gram schmidt and arnoldi step

        !-------- now  update factorization of hh

        !-------- perform previous transformations  on i-th column of h
        do k=2,i
           k1 = k-1
           t = hh(k1,i)
           hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
           hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
        end do
        gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
        if (gam == 0.0) gam = epsmac
        !-----------#  determine next plane rotation  #-------------------
        c(i) = hh(i,i)/gam
        s(i) = hh(i1,i)/gam
        rs(i1) = -s(i)*rs(i)
        rs(i) =  c(i)*rs(i)
        !---determine residual norm and test for convergence-
        hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
        ro = abs(rs(i1))
        if (oktest) then
           select case(TypeStop)
           case('rel')
              write(*,*) its,' matvecs, ',' ||rn||/||r0|| =',ro/ro0
           case('abs')
              write(*,*) its,' matvecs, ',' ||rn|| =',ro
           end select
        end if
        if (i >= nKrylov .or. (ro <= Tol1)) exit KRYLOVLOOP
     enddo KRYLOVLOOP

     !
     ! now compute solution. first solve upper triangular system.
     !
     ! rs := hh(1:i,1:i) ^-1 * rs

     do j=i,1,-1
        if (rs(j)/=0.0) then
           rs(j)=rs(j)/hh(j,j)
           tmp = rs(j)
           do k=j-1,1,-1
              rs(k) = rs(k) - tmp*hh(k,j)
           enddo
        endif
     enddo

     ! done with back substitution..
     ! now form linear combination to get solution
     do j=1, i
        t = rs(j)
        Sol = Sol + t*Krylov_II(:,j)
     end do

     ! exit from outer loop if converged or too many iterations
     if (ro <= Tol1 .or. its >= Iter) exit RESTARTLOOP
  end do RESTARTLOOP

  Iter=its
  Tol=Tol/Tol1*ro ! (relative) tolerance achieved
  if(ro < Tol1)then
     info = 0
  elseif(ro < ro0)then
     info =  2
  else
     info = -2
  endif
contains
  !============================================================================
  real function dot_product_mpi(a,b)

    use ModProcMH
    use ModMpi
    implicit none
    real, intent(in)    :: a(:), b(:)

    real :: local_dot_product, global_dot_product
    integer :: iError
    !--------------------------------------------------------------------------

    local_dot_product = dot_product(a,b)

    call MPI_allreduce(local_dot_product, global_dot_product, 1, MPI_REAL, &
         MPI_SUM, iComm, iError)

    dot_product_mpi = global_dot_product

  end function dot_product_mpi
  !============================================================================
end subroutine gmres

