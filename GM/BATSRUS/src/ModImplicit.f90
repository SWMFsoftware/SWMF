!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
module ModImplicit

  use ModSize

  implicit none
  SAVE

  !\
  ! Logical determining if the scheme is part or full implicit (spatially)
  ! and whether we make it fully conservative
  !/
  LOGICAL :: UsePartImplicit, UseFullImplicit, UseConservativeImplicit

  !\
  ! Logical arrays determining if a block is handled implicitly
  !/
  LOGICAL, DIMENSION(nBLK) :: ImplicitBLK, SkippedBLK

  !\
  ! Parameters for selecting implicit blocks
  !/
  CHARACTER (LEN=10) :: ImplCritType
  REAL               :: Rimplicit

  !\
  ! Parameters for the Implicit Time Stepping
  !/

  ! Conversion between BATSRUS and VAC index order for variables
  INTEGER, PARAMETER, DIMENSION(8):: bat2vac=(/1,2,3,4,8,5,6,7/)

  INTEGER, PARAMETER:: nw=8

  ! Number of cells and unknowns per block
  INTEGER, PARAMETER:: nwIJK = nw*nIJK

  ! Maximum number of unknowns per processor
  INTEGER, PARAMETER:: MaxImplVar=nwIJK*MaxImplBLK

  REAL,    PARAMETER:: smalldouble=1.E-30, bigdouble=1.E+30

  !-- Common variables:

  !Kronecker delta
  INTEGER :: kr(3,3)
  DATA kr /1,0,0,0,1,0,0,0,1/

  ! Actual number of implicitly treated blocks for a processor
  INTEGER :: nImplBLK

  ! Actual number of implicit variables (unknowns) per processor, total
  INTEGER :: nimpl, nimpl_total

  ! Test variables
  INTEGER :: implBLKtest, implVARtest

  ! Indirect index array from implicit block index to general block index
  INTEGER :: impl2iBLK(MaxImplBLK)

  ! Parameters for the implicit techniques

  ! Implicit scheme parameters
  REAL   :: explCFL, implCFL, ImplCoeff0, ImplCoeff
  LOGICAL:: UseBDF2
  LOGICAL:: ImplSource

  REAL :: impldwlimit !!!

  ! Newton iteration parameters
  LOGICAL :: UseNewton=.false.
  LOGICAL :: NewMatrix=.true.
  INTEGER :: NewtonIterMax=1

  ! Jacobian parameters
  CHARACTER (LEN=10) :: JacobianType='prec'
  REAL               :: JacobianEps =1.E-12

  ! Preconditioner parameters
  REAL               :: GustafssonPar=0.5
  CHARACTER (LEN=10) :: PrecondSide
  CHARACTER (LEN=10) :: PrecondType

  ! Krylov scheme parameters
  CHARACTER (LEN=10) :: KrylovType     ='bicgstab'
  CHARACTER (LEN=10) :: KrylovInitType ='nul'
  INTEGER            :: KrylovMatvecMax=100
  REAL               :: KrylovErrorMax =0.001
  INTEGER            :: nKrylovVector  =100

  ! Implicit scheme
  INTEGER            :: nORDER_impl
  CHARACTER (LEN=10) :: FluxTypeImpl

  LOGICAL:: sourceunsplit=.false.
  LOGICAL:: compactres=.true.

  ! Second norm of the variables and the residual
  REAL:: wnrm(nw), residual

  ! Previous time step, explicit time step, ratio of explicit and implicit dt
  REAL:: dtexpl, dt_prev, dtcoeff

  ! The k-th Newton iterate, one layer of ghost cells for the Jacobian
  REAL, DIMENSION(0:nI+1,0:nJ+1,0:nK+1,nw,MaxImplBLK) :: w_k

  ! Previous (n-1) state
  REAL, DIMENSION(nI,nJ,nK,nw,nBLK) :: w_prev

  ! low and high order residuals
  REAL, DIMENSION(nI,nJ,nK,nw,MaxImplBLK) :: RES_impl, RES_expl

  ! Time step when the previous state was stored
  INTEGER :: n_prev=-100

  ! Heptadiagonal Jacobian matrix
  INTEGER, PARAMETER:: nstencil=2*ndim+1
  REAL :: MAT(nw,nw,1:nI,1:nJ,1:nK,nstencil,MaxImplBLK)

  ! Right hand side for UseNewton, right hand side, dw=w_n+1 - w_n
  REAL, DIMENSION(MaxImplVar) :: rhs0, rhs, dw

  ! Is dw=0 initiallly
  LOGICAL :: non0dw

  ! Counters for reports
  INTEGER:: nexpl=0,nnewton=0,niterimpl=0,nmatvec=0
  
  INTEGER,PARAMETER::E_=8

end module ModImplicit
