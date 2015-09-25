!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_particles
  use BATL_interpolate_amr_wrapper
  implicit none
  SAVE
  logical  ::UseParticles = .false.
  !\
  ! Number of different sorts of particles
  !/ 
  integer, parameter:: nPType = 1

  type particles
     !\
     ! The number of parameters which characterize the 'state' of
     ! the particles
     !/
     integer:: nVar 
     integer:: nParticle
     !\
     ! nVar*nParticleMax array. The econd index numerates 'particles'.
     !First nDim components are the cartezian coordinates the other 
     !can be are velocities, the momentum components, gyrokinetic
     !parameters or the physical particles, or, alternatively,
     !the magnetic field line ID and the number of the Lagrangian
     !fluid particle along the magnetic field line in application to
     !M-FLAMPA
     real,    pointer  :: State_VI(:,:)
     !\
     ! Array of the length of nParticleMax with the index enumerating
     ! particles. It stores the number of block which posesses an 
     ! information, sufficient to interpolate (with the possible use 
     ! of ghost cells) an information 
     integer,  pointer :: iBlock_I(:)
  end type particles
  type(particles),dimension(nPType) ::Of_I
  
contains
  !================
  subroutine set_pointer_to_particles(&
       iSortParticle, State_VI, iBlock_I, nVar, nParticle)
    integer,          intent(in)    :: iSortParticle
    real,    pointer, intent(inout) :: State_VI(:,:)
    integer, pointer, intent(inout) :: iBlock_I(:)
    integer,          intent(out)   :: nVar 
    integer,          intent(out) :: nParticle
    !------------
    nullify(State_VI); nullify(iBlock_I)
    State_VI  => Of_I(iSortParticle)%State_VI
    iBlock_I  => Of_I(iSortParticle)%iBlock_I
    nVar      =  Of_I(iSortParticle)%nVar
    nParticle =  Of_I(iSortParticle)%nParticle
  end subroutine set_pointer_to_particles
  !=================
  subroutine message_pass_particles
    integer          :: iSortParticle
    real,    pointer :: State_VI(:,:)
    integer, pointer :: iBlock_I(:)
    integer          :: nVar 
    integer          :: nParticle
    !--------------
    do iSortParticle = 1, nPType
       call set_pointer_to_particles(&
            iSortParticle, State_VI, &
            iBlock_I, nVar,nParticle)
       call pass_this_sort
       Of_I(iSortParticle)%nParticle = nParticle
    end do
  contains
    subroutine pass_this_sort
      use BATL_mpi
      integer:: iParticle
    end subroutine pass_this_sort
  end subroutine message_pass_particles
end module BATL_particles
