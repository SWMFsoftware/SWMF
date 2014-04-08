!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSemiImplVar

  implicit none

  ! Semi-implicit variables that cannot be in ModSemiImplicit
  ! because it would create circular dependencies

  logical, public:: UseSemiImplicit = .false.
  logical, public:: UseSplitSemiImplicit = .false.
  character(len=40), public:: TypeSemiImplicit = 'none'

  ! Number of all semi-implicit variables
  integer, public:: nVarSemiAll

  ! Index range of semi-implicit variables solved together
  integer, public:: iVarSemiMin, iVarSemiMax

  ! Number of semi-implicit variables solved together
  integer, public:: nVarSemi

  ! Named indices for semi-implicit variables
  integer, public :: iTeImpl=0, iTrImplFirst=0, iTrImplLast=0, iEradImpl=0

  ! Number of semi-implicit grid blocks
  integer, public:: nBlockSemi

  ! Conversion from compact block index iBlockSemi to normal index iBlock
  integer, public, allocatable:: iBlockFromSemi_B(:)

  ! Arrays for flux correction at resolution changes
  real, public, allocatable, dimension(:,:,:,:,:) :: &
       FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB

  ! Increase accurary of radiation solver by switching to third order
  ! interpolation at the resolution change
  logical, public:: UseAccurateRadiation = .false.

  ! Coefficient of the implicit part: 1 - Backward Euler, 0.5 - trapezoidal
  real, public:: SemiImplCoeff = 1.0

  ! This array is indexed with normal block index and has nG ghost cells
  real, allocatable:: SemiState_VGB(:,:,:,:,:)  ! Semi-implicit vars

end module ModSemiImplVar

!============================================================================

module ModImplicit

  use ModSize
  use BATL_size,     ONLY: nDim
  use ModVarIndexes, ONLY: nVar, P_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc, iComm
  use ModLinearSolver, ONLY: LinearSolverParamType

  ! for sake of fewer changes in the user modules
  use ModSemiImplVar, StateSemi_VGB => SemiState_VGB


  implicit none
  SAVE

  ! Heptadiagonal Jacobian matrix
  integer, parameter:: nStencil = 2*nDim + 1
  integer, parameter:: Stencil1_ = 1, Stencil2_ = 2, Stencil3_ = 3, &
       Stencil4_ = min(4,nStencil), Stencil5_ = min(5,nStencil), &
       Stencil6_ = min(6,nStencil), Stencil7_ = min(7,nStencil)

  !\
  ! Number of variables and the index of energy (replacing pressure)
  !/
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

  ! Shall we zero out contribution from ghost cells
  logical :: UseNoOverlap = .true.

  !\
  ! Parameters for selecting implicit blocks
  !/
  character (len=10) :: ImplCritType = 'dt'
  real               :: Rimplicit

  ! Second order in time implicit scheme (uses previous time step)
  logical:: UseBDF2    = .false.

  ! Time step when the previous state was stored
  integer :: n_prev=-100

  ! Previous (n-1) state
  real, allocatable :: ImplOld_VCB(:,:,:,:,:)

  ! Previous time step
  real::  dt_prev = -1.0

  ! CFL number for explicit scheme
  real   :: ExplCfl = 0.8

  ! Implicit scheme
  character (len=10) :: FluxTypeImpl = 'default' ! same as explicit

end module ModImplicit
