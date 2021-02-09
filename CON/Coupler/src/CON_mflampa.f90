!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_mflampa
  ! allocate the grid used in this model
  use ModUtilities,      ONLY: check_allocate
  implicit none
    ! Number of variables in the state vector and the identifications
  integer, public, parameter :: nMHData = 13,                 &
       LagrID_     = 0, & ! Lagrangian id           ^saved/   ^set to 0
       X_          = 1, & !                         |read in  |in copy_
       Y_          = 2, & ! Cartesian coordinates   |restart  |old_stat
       Z_          = 3, & !                         v/        |saved to
       Rho_        = 4, & ! Background plasma density         |mhd1
       T_          = 5, & ! Background temperature            |
       Ux_         = 6, & !                                   |may be
       Uy_         = 7, & ! Background plasma bulk velocity   |read from
       Uz_         = 8, & !                                   |mhd1
       Bx_         = 9, & !                                   |or
       By_         =10, & ! Background magnetic field         |received
       Bz_         =11, & !                                   |from
       Wave1_      =12, & !\                                  |coupler
       Wave2_      =13    ! Alfven wave turbulence            v
  !
  ! State vector is a pointer, which is joined to a target array
  ! For stand alone version the target array is allocated here
  !
  real, allocatable, target:: Target_VIB(:,:,:)
  ! Grid integer parameters:
  integer :: nParticleMF, nBlockMF
  ! Misc:
  integer:: iError
contains
  !============================================================================
  subroutine set_state_pointer(State_VIB, nBlock, nParticle)
    real, intent(inout), pointer :: State_VIB(:,:,:)
    integer, intent(in)          :: nBlock, nParticle
    integer :: iParticle
    character(len=*), parameter:: NameSub = 'set_state_pointer'
    !--------------------------------------------------------------------------
    !
    ! Store nBlock and nParticleMax
    nBlockMF    = nBlock
    nParticleMF = nParticle
    allocate(Target_VIB(LagrID_:nMHData, 1:nParticleMF, 1:nBlockMF), &
         stat=iError)
    call check_allocate(iError, NameSub//'Target_VIB')
    State_VIB => Target_VIB
    !
    State_VIB = 0.0
    !
    ! reset lagrangian ids
    !
    do iParticle = 1, nParticleMF
       State_VIB(LagrID_, iParticle, 1:nBlock) = real(iParticle)
    end do
  end subroutine set_state_pointer
  !============================================================================
end module CON_mflampa

