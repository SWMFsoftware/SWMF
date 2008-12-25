!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
module ModImplicit

  use ModSize
  use ModVarIndexes, ONLY: nVar, P_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  SAVE

  !\
  ! Number of variables and the index of energy (replacing pressure)
  !/
  integer :: nw
  integer, parameter :: E_ = P_

  !\
  ! Logicals determining if and what implcit scheme is used
  !/
  logical :: UseImplicit = .false.

  ! Use fully implicit scheme (all blocks)
  logical :: UseFullImplicit=.false.  

  ! Use part implicit scheme (some blocks)
  logical :: UsePartImplicit=.false. 

  ! Use temporally 2nd order part implicit by doing explicit predictor step
  ! in all the blocks to get second order flux at expl/impl interfaces
  logical :: UsePartImplicit2=.false. 

  ! Do a conservative update at the end of the iterative solver 
  ! This does not work well without full Newton iteration
  logical :: UseConservativeImplicit=.false.  

  ! Use semi-implicit scheme
  logical :: UseSemiImplicit = .false.
  character(len=40) :: TypeSemiImplicit = 'radiation'

  !\
  ! Parameters for selecting implicit blocks
  !/
  character (len=10) :: ImplCritType
  real               :: Rimplicit

  !\
  ! Parameters for the Implicit Time Stepping
  !/
  ! Number of cells and unknowns per block
  integer:: nwIJK

  ! Maximum number of unknowns per processor
  integer:: MaxImplVar

  real,    parameter:: smalldouble=1.E-30, bigdouble=1.E+30

  !-- Common variables:

  !Kronecker delta
  integer, parameter :: kr(3,3) = reshape( (/1,0,0,0,1,0,0,0,1/), (/3,3/) )

  ! Actual number of implicitly treated blocks for a processor
  integer :: nImplBLK=0

  ! Actual number of implicit variables (unknowns) per processor, total
  integer :: nimpl=0, nimpl_total=0

  ! Test variables
  integer :: implBLKtest=1, implVARtest=1

  ! Indirect index array from implicit block index to general block index
  integer :: impl2iBLK(MaxImplBLK)

  ! Parameters for the implicit techniques

  ! Implicit scheme parameters
  real   :: explCFL=0.8, implCFL=100.0, ImplCoeff0=1.0, ImplCoeff=1.0
  logical:: UseBDF2    = .false.
  logical:: ImplSource = .false.

  ! Newton iteration parameters
  logical :: UseNewton     = .false.
  logical :: NewMatrix     = .true.
  integer :: NewtonIterMax = 1

  ! Jacobian parameters
  character (len=10) :: JacobianType = 'prec'
  real               :: JacobianEps  = 1.E-12

  ! Preconditioner parameters
  real               :: GustafssonPar = 0.5
  character (len=10) :: PrecondSide   = 'symmetric'
  character (len=10) :: PrecondType   = 'MBILU'

  ! Krylov scheme parameters
  character (len=10) :: KrylovType     ='gmres'
  character (len=10) :: KrylovInitType ='nul'
  integer            :: KrylovMatvecMax=100
  real               :: KrylovErrorMax =0.001
  integer            :: nKrylovVector  =100

  ! Implicit scheme
  integer            :: nOrder_Impl  = 1
  character (len=10) :: FluxTypeImpl = 'default' ! same as explicit

  logical:: sourceunsplit = .false.
  logical:: compactres    = .true.

  ! Previous time step, explicit time step, ratio of explicit and implicit dt
  real:: dtexpl=-1.0, dt_prev=-1.0, dtcoeff=-1.0

  ! Second norm of the variables and the residual
  real:: residual=-1.0
  real, allocatable:: wnrm(:)

  ! The k-th Newton iterate, one layer of ghost cells for the Jacobian
  real, allocatable :: w_k(:,:,:,:,:)

  ! Previous (n-1) state
  real, allocatable :: w_prev(:,:,:,:,:)

  ! low and high order residuals
  real, allocatable :: RES_impl(:,:,:,:,:)
  real, allocatable :: RES_expl(:,:,:,:,:)

  ! semi-implicit state
  real, allocatable :: StateSemi_VGB(:,:,:,:,:)

  ! Time step when the previous state was stored
  integer :: n_prev=-100

  ! Heptadiagonal Jacobian matrix
  integer, parameter:: nstencil=2*ndim+1
  real, allocatable :: MAT(:,:,:,:,:,:,:)

  ! Right hand side for UseNewton, right hand side, dw=w_n+1 - w_n
  real, allocatable :: rhs0(:)
  real, allocatable :: rhs(:)
  real, allocatable :: dw(:)

  ! Is dw=0 initiallly
  logical :: non0dw = .false.

  ! Counters for reports
  integer:: nExpl=0, nNewton=0, nIterImpl=0, nMatvec=0

  ! Update check
  real :: &
       RejectStepLevel   = 0.3, RejectStepFactor   = 0.50, &
       ReduceStepLevel   = 0.6, ReduceStepFactor   = 0.95, &
       IncreaseStepLevel = 0.8, IncreaseStepFactor = 1.05

contains
  !============================================================================
  subroutine read_implicit_param(NameCommand)

    use ModReadParam,     ONLY: read_var
    use ModMain,          ONLY: DoSplitDb0Dt
    use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B

    character(len=*), intent(in) :: NameCommand
    character(len=*), parameter:: NameSub = 'read_implicit_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#IMPLICIT")                                 !^CFG IF IMPLICIT BEGIN
       call read_var('UsePointImplicit', UsePointImplicit)
       call read_var('UsePartImplicit',  UsePartImplicit)
       call read_var('UseFullImplicit',  UseFullImplicit)

       UseImplicit = UseFullImplicit .or. UsePartImplicit
       UsePointImplicit_B = UsePointImplicit

       ! For implicit scheme it is better to use unsplit dB0/dt evaluation
       DoSplitDb0Dt = .not.UseImplicit

       if(UsePartImplicit  .and. UseFullImplicit) call stop_mpi(&
            'Only one of UsePartImplicit and UseFullImplicit can be true')

       if(UseImplicit)call read_var('ImplCFL',ImplCFL)

    case("#IMPLCRITERIA", "#IMPLICITCRITERIA", "#STEPPINGCRITERIA")
       call read_var('TypeImplCrit', ImplCritType)
       select case(ImplCritType)
       case('R','r')
          call read_var('rImplicit', rImplicit)
       case('test','dt')
       case default
          if(iProc==0)then
             write(*,'(a)')NameSub// &
                  ' WARNING: invalid ImplCritType='//trim(ImplCritType)// &
                  ' !!!'
             write(*,*)NameSub//' Setting ImplCritType=dt'
          end if
          ImplCritType = 'dt'
       end select

    case("#PARTIMPL", "#PARTIMPLICIT")
       call read_var('UsePartImplicit2',UsePartImplicit2)

    case("#SEMIIMPLICIT", "#SEMIIMPL")
       call read_var('UseSemiImplicit', UseSemiImplicit)
       if(UseSemiImplicit)call read_var('TypeSemiImplicit', TypeSemiImplicit)

    case("#IMPLSCHEME", "#IMPLICITSCHEME")
       call read_var('nOrderImpl',nOrder_Impl)
       call read_var('TypeFluxImpl',FluxTypeImpl,IsUpperCase=.true.)

    case("#IMPLSTEP", "#IMPLICITSTEP")
       call read_var('ImplCoeff ',ImplCoeff0)
       call read_var('UseBDF2   ',UseBdf2)
       call read_var('ImplSource',ImplSource)

    case("#IMPLCHECK", "#IMPLICITCHECK")
       call read_var('RejectStepLevel' ,   RejectStepLevel)
       call read_var('RejectStepFactor',   RejectStepFactor)
       call read_var('ReduceStepLevel' ,   ReduceStepLevel)
       call read_var('ReduceStepFactor',   ReduceStepFactor)
       call read_var('IncreaseStepLevel' , IncreaseStepLevel)
       call read_var('IncreaseStepFactor', IncreaseStepFactor)

    case("#NEWTON")
       call read_var('UseConservativeImplicit',UseConservativeImplicit)
       call read_var('UseNewton',UseNewton)
       if(UseNewton)then
          call read_var('UseNewMatrix ',NewMatrix)
          call read_var('MaxIterNewton',NewtonIterMax)
       endif

    case("#JACOBIAN")
       call read_var('TypeJacobian',JacobianType)
       call read_var('JacobianEps', JacobianEps)

    case("#PRECONDITIONER")
       call read_var('TypePrecondSide',PrecondSide)
       call read_var('TypePrecond'    ,PrecondType)
       call read_var('GustafssonPar'  ,GustafssonPar)
       if(GustafssonPar<=0.)GustafssonPar=0.
       if(GustafssonPar>1.)GustafssonPar=1.

    case("#KRYLOV")
       call read_var('TypeKrylov'     ,KrylovType)
       call read_var('TypeInitKrylov' ,KrylovInitType)
       call read_var('ErrorMaxKrylov' ,KrylovErrorMax)
       call read_var('MaxMatvecKrylov',KrylovMatvecMax)
       nKrylovVector = KrylovMatvecMax

    case("#KRYLOVSIZE")
       call read_var('nKrylovVector',nKrylovVector)

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_implicit_param
  !============================================================================
  subroutine init_mod_implicit

    use ModUtilities, ONLY: check_allocate
    integer :: iError
    !----------------------------------------------------------------------
    if(allocated(w_k)) return

    if(UseSemiImplicit)then
       select case(TypeSemiImplicit)
       case('radiation')
          nw = 2
       end select
    else
       nw = nVar
    end if
    nwIJK = nw*nIJK
    MaxImplVar=nwIJK*MaxImplBLK

    allocate(w_k(0:nI+1,0:nJ+1,0:nK+1,nw,MaxImplBLK))
    allocate(w_prev(nI,nJ,nK,nw,nBLK))
    allocate(RES_impl(nI,nJ,nK,nw,MaxImplBLK))
    allocate(RES_expl(nI,nJ,nK,nw,MaxImplBLK))
    allocate(StateSemi_VGB(nw,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
    allocate(MAT(nw,nw,1:nI,1:nJ,1:nK,nstencil,MaxImplBLK))
    allocate(rhs0(MaxImplVar))
    allocate(rhs(MaxImplVar))
    allocate(dw(MaxImplVar))
    allocate(wnrm(nw))

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_implicit allocated arrays'
    end if

    wnrm = 0.0

  end subroutine init_mod_implicit
  !============================================================================
  subroutine clean_mod_implicit

    if(.not.allocated(w_k)) return
    deallocate(w_k)
    deallocate(w_prev)
    deallocate(RES_impl)
    deallocate(RES_expl)
    deallocate(MAT)
    deallocate(rhs0)
    deallocate(rhs)
    deallocate(dw)
    deallocate(wnrm)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_implicit deallocated arrays'
    end if

  end subroutine clean_mod_implicit

end module ModImplicit
