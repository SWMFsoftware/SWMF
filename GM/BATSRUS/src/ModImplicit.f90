!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
module ModImplicit

  use ModSize
  use BATL_size,     ONLY: nDimBatl => nDim
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

  ! Use energy or pressure as the implicit variable
  logical :: UseImplicitEnergy = .true.

  ! Do a conservative update at the end of the iterative solver 
  ! This does not work well without full Newton iteration
  logical :: UseConservativeImplicit=.false.  

  ! Use semi-implicit scheme
  logical :: UseSemiImplicit = .false.
  character(len=40) :: TypeSemiImplicit = 'radiation'

  logical:: UseSplitSemiImplicit = .false.
  ! Index range for split semi implicit scheme
  integer:: iVarSemiMin, iVarSemiMax
  ! Number and index of implicit variables per splitting
  integer :: nVarSemi, iVarSemi

  ! Named indices for semi-implicit variables
  integer :: iTeImpl=0, iTrImplFirst=0, iTrImplLast=0, iEradImpl=0

  ! Increase accurary of radiation solver by switching to third order
  ! interpolation at the resolution change
  logical, public :: UseAccurateRadiation = .false.

  ! Shall we zero out contribution from ghost cells
  logical :: UseNoOverlap = .true.

  !\
  ! Parameters for selecting implicit blocks
  !/
  character (len=10) :: ImplCritType = 'dt'
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
  real               :: PrecondParam  = -0.5
  character (len=10) :: PrecondSide   = 'symmetric'
  character (len=10) :: PrecondType   = 'MBILU'
  integer            :: DnInitHypreAmg= 1

  ! Krylov scheme parameters
  character (len=10) :: KrylovType     ='GMRES'
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
  real, allocatable :: Impl_VGB(:,:,:,:,:)

  ! Previous (n-1) state
  real, allocatable :: ImplOld_VCB(:,:,:,:,:)

  ! low and high order residuals
  real, allocatable :: ResImpl_VCB(:,:,:,:,:)
  real, allocatable :: ResExpl_VCB(:,:,:,:,:)

  ! semi-implicit state
  real, allocatable :: StateSemi_VGB(:,:,:,:,:)

  real, allocatable :: DconsDsemi_VCB(:,:,:,:,:)

  real, allocatable, dimension(:,:,:,:,:) :: &
       FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB

  ! Time step when the previous state was stored
  integer :: n_prev=-100

  ! Heptadiagonal Jacobian matrix
  integer, parameter:: nStencil = 2*nDimBatl + 1
  integer, parameter:: Stencil1_ = 1, Stencil2_ = 2, Stencil3_ = 3, &
       Stencil4_ = min(4,nStencil), Stencil5_ = min(5,nStencil), &
       Stencil6_ = min(6,nStencil), Stencil7_ = min(7,nStencil)

  real, allocatable :: MAT(:,:,:,:,:,:,:)
  real, allocatable :: JacobiPrec_I(:)

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
    use ModLinearSolver,  ONLY: &
         Jacobi_, BlockJacobi_, GaussSeidel_, Dilu_, Bilu_

    character(len=*), intent(in) :: NameCommand

    integer :: i
    character(len=*), parameter:: NameSub = 'read_implicit_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#IMPLICIT')                                 !^CFG IF IMPLICIT BEGIN
       call read_var('UsePointImplicit', UsePointImplicit)
       call read_var('UsePartImplicit',  UsePartImplicit)
       call read_var('UseFullImplicit',  UseFullImplicit)

       UseImplicit = UseFullImplicit .or. UsePartImplicit
       UsePointImplicit_B = UsePointImplicit

       ! For implicit scheme it is better to use unsplit dB0/dt evaluation
!!!       DoSplitDb0Dt = .not.UseImplicit

       if(UsePartImplicit  .and. UseFullImplicit) call stop_mpi(&
            'Only one of UsePartImplicit and UseFullImplicit can be true')

       if(UseImplicit)call read_var('ImplCFL',ImplCFL)

    case('#IMPLCRITERIA', '#IMPLICITCRITERIA', '#STEPPINGCRITERIA')
       call read_var('TypeImplCrit', ImplCritType, IsLowerCase=.true.)
       select case(ImplCritType)
       case('r')
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

    case('#PARTIMPL', '#PARTIMPLICIT')
       call read_var('UsePartImplicit2', UsePartImplicit2)

    case('#IMPLENERGY', '#IMPLICITENERGY')
       call read_var('UseImplicitEnergy', UseImplicitEnergy)

    case('#SEMIIMPLICIT', '#SEMIIMPL')
       call read_var('UseSemiImplicit', UseSemiImplicit)
       if(UseSemiImplicit)then
          call read_var('TypeSemiImplicit', TypeSemiImplicit)
          i = index(TypeSemiImplicit,'split')
          UseSplitSemiImplicit = i > 0
          if(UseSplitSemiImplicit) &
               TypeSemiImplicit = TypeSemiImplicit(1:i-1)
       end if

    case('#IMPLSCHEME', '#IMPLICITSCHEME')
       call read_var('nOrderImpl', nOrder_Impl)
       call read_var('TypeFluxImpl', FluxTypeImpl, IsUpperCase=.true.)

    case('#IMPLSTEP', '#IMPLICITSTEP')
       call read_var('ImplCoeff ',ImplCoeff0)
       call read_var('UseBDF2   ',UseBdf2)
       call read_var('ImplSource',ImplSource)

    case('#IMPLCHECK', '#IMPLICITCHECK')
       call read_var('RejectStepLevel' ,   RejectStepLevel)
       call read_var('RejectStepFactor',   RejectStepFactor)
       call read_var('ReduceStepLevel' ,   ReduceStepLevel)
       call read_var('ReduceStepFactor',   ReduceStepFactor)
       call read_var('IncreaseStepLevel' , IncreaseStepLevel)
       call read_var('IncreaseStepFactor', IncreaseStepFactor)

    case('#NEWTON')
       call read_var('UseConservativeImplicit',UseConservativeImplicit)
       call read_var('UseNewton',UseNewton)
       if(UseNewton)then
          call read_var('UseNewMatrix ', NewMatrix)
          call read_var('MaxIterNewton', NewtonIterMax)
       endif

    case('#JACOBIAN')
       call read_var('TypeJacobian',JacobianType, IsLowerCase=.true.)
       call read_var('JacobianEps', JacobianEps)

    case('#PRECONDITIONER')
       UseNoOverlap = .true.
       call read_var('TypePrecondSide',PrecondSide, IsLowerCase=.true.)
       call read_var('TypePrecond'    ,PrecondType, IsUpperCase=.true.)
       select case(PrecondType)
       case('HYPRE')
          PrecondSide = 'left'
          UseNoOverlap = .false.
          call read_var('DnInitHypreAmg', DnInitHypreAmg)
       case('JACOBI')
          PrecondParam = Jacobi_
          PrecondSide  = 'left'
       case('BLOCKJACOBI')
          PrecondParam = BlockJacobi_
          PrecondSide  = 'left'
       case('GS')
          PrecondParam = GaussSeidel_
       case('DILU')
          PrecondParam = Dilu_
          PrecondSide  = 'left'
       case('BILU')
          PrecondParam = Bilu_
       case('MBILU')
          call read_var('GustafssonPar', PrecondParam)
          PrecondParam = -PrecondParam
       case default
          call stop_mpi(NameSub//' invalid TypePrecond='//PrecondType)
       end select

    case('#KRYLOV')
       call read_var('TypeKrylov'     ,KrylovType, IsUpperCase=.true.)
       call read_var('TypeInitKrylov' ,KrylovInitType, IsLowerCase=.true.)
       call read_var('ErrorMaxKrylov' ,KrylovErrorMax)
       call read_var('MaxMatvecKrylov',KrylovMatvecMax)
       nKrylovVector = KrylovMatvecMax

    case('#KRYLOVSIZE')
       call read_var('nKrylovVector',nKrylovVector)

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_implicit_param
  !============================================================================
  subroutine init_mod_implicit

    use ModUtilities,  ONLY: check_allocate
    use ModVarIndexes, ONLY: nWave
    use ModLinearSolver, ONLY: Bilu_

    integer :: iError

    character(len=*), parameter:: NameSub = 'init_mod_implicit'
    !----------------------------------------------------------------------
    if(allocated(Impl_VGB)) return

    if(UseSemiImplicit)then
       select case(TypeSemiImplicit)
       case('radiation')
          if(nWave == 1 .or. UseSplitSemiImplicit)then
             ! Waves with point implicit temperature
             nw = nWave
          else
             ! (electron) temperature and waves
             nw = 1 + nWave
          end if
       case('cond', 'parcond')
          nw = 1
       case('radcond')
          nw = nWave + 1
       case('resistivity')
          nw = MaxDim
       case default
          call stop_mpi(NameSub//': nw unknown for'//TypeSemiImplicit)
       end select

    else
       ! Solve for all variables implicitly. nVarSemi has to be set too.
       nw       = nVar
    end if

    if(UseSplitSemiImplicit)then
       ! Split semi-implicit scheme solves 1 implicit variable at a time
       nVarSemi = 1
    else
       ! Unsplit (semi-)implicit scheme solves all impl. variables together
       nVarSemi = nw
    end if
    ! Set range of (semi-)implicit variables for unsplit scheme.
    ! For split scheme this will be overwritten.
    iVarSemiMin = 1
    iVarSemiMax = nVarSemi

    if(UseSemiImplicit)then
       if(TypeSemiImplicit == 'parcond' .or. &
            TypeSemiImplicit == 'resistivity' .or. UseAccurateRadiation)then
          allocate( &
               FluxImpl_VXB(nVarSemi,nJ,nK,2,MaxBlock), &
               FluxImpl_VYB(nVarSemi,nI,nK,2,MaxBlock), &
               FluxImpl_VZB(nVarSemi,nI,nJ,2,MaxBlock) )
          FluxImpl_VXB = 0.0
          FluxImpl_VYB = 0.0
          FluxImpl_VZB = 0.0
       end if
    end if

    nwIJK      = nVarSemi*nIJK
    MaxImplVar = nwIJK*MaxImplBLK

    ! Fix preconditioner type from DILU to BILU for a single variable
    ! because BILU is optimized for scalar equation.
    if(nVarSemi == 1 .and. PrecondType == 'DILU')then
       PrecondType  = 'BILU'
       PrecondParam = Bilu_
    end if

    ! Arrays for all implicit variables
    allocate(Impl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBLK))
    allocate(ImplOld_VCB(nw,nI,nJ,nK,nBLK))
    allocate(wnrm(nw))
    if(UseSemiImplicit) allocate(DconsDsemi_VCB(nw,nI,nJ,nK,MaxImplBLK))

    ! Arrays for split implicit variables
    allocate(StateSemi_VGB(nVarSemi,-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock))
    allocate(ResExpl_VCB(nVarSemi,nI,nJ,nK,MaxImplBLK))
    allocate(ResImpl_VCB(nVarSemi,nI,nJ,nK,MaxImplBLK))

    allocate(MAT(nVarSemi,nVarSemi,nI,nJ,nK,nStencil,MaxImplBLK))
    allocate(rhs0(MaxImplVar))
    allocate(rhs(MaxImplVar))
    allocate(dw(MaxImplVar))
    if(PrecondType == 'JACOBI') allocate(JacobiPrec_I(MaxImplVar))

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_implicit allocated arrays'
    end if

    wnrm = 0.0

  end subroutine init_mod_implicit
  !============================================================================
  subroutine clean_mod_implicit

    if(.not.allocated(Impl_VGB)) return
    deallocate(Impl_VGB)
    deallocate(ImplOld_VCB)
    deallocate(ResImpl_VCB)
    deallocate(ResExpl_VCB)
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
