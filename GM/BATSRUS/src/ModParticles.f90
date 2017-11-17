!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModParticles

  use BATL_lib, ONLY: &
       test_start, test_stop
  use BATL_particles, ONLY: & 
       Particle_I, message_pass_particles, remove_undefined_particles, &
       mark_undefined, check_particle_location, put_particles, trace_particles
  use ModBatlInterface, ONLY: interpolate_grid_amr_gc
  use BATL_size,        ONLY: BatlNKind=>nKindParticle
  implicit none
  SAVE
  integer, private :: nKindParticle = 0
contains
  !============================================================================
  subroutine allocate_particles(iKindParticle, nVar, nIndex, nParticleMax)
    integer, intent(inout) :: iKindParticle
    integer, intent(in)    :: nVar, nIndex, nParticleMax

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'allocate_particles'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(iKindParticle > 0) then
       !\
       ! This sort of particles is already allocated
       !/
       if(nVar /= Particle_I(iKindParticle)%nVar .or. &
            nIndex/=Particle_I(iKindParticle)%nIndex)call stop_mpi(&
            NameSub//':'//&
            'Particles are already allocated for different purpose')
       if(nParticleMax <= Particle_I(iKindParticle)%nParticleMax)&
            RETURN
       !\
       ! Deallocate particles to be allocated with larger nParticles
       !/
       deallocate(Particle_I(iKindParticle)%State_VI)
       deallocate(Particle_I(iKindParticle)%iIndex_II)
    else
       !\
       ! Particles are allocated from scratch
       !/
       nKindParticle = nKindParticle + 1
       if(nKindParticle > BatlNKind)call stop_mpi(NameSub//':'//&
            'Too many sorts of particles are allocated')
       iKindParticle = nKindParticle
       Particle_I(iKindParticle)%nVar   = nVar
       Particle_I(iKindParticle)%nIndex = nIndex
    end if
    Particle_I(iKindParticle)%nParticleMax = nParticleMax

    nullify( Particle_I(iKindParticle)%State_VI)
    allocate(Particle_I(&
         iKindParticle)%State_VI(1:nVar,1:nParticleMax))

    nullify( Particle_I(iKindParticle)%iIndex_II)
    allocate(Particle_I(&
         iKindParticle)%iIndex_II(0:nIndex,1:nParticleMax))

    Particle_I(iKindParticle)%State_VI( :,:) = 0.0
    Particle_I(iKindParticle)%iIndex_II(:,:) = 0
    Particle_I(iKindParticle)%nParticle      = 0
    call test_stop(NameSub, DoTest)
  end subroutine allocate_particles
  !================================
  subroutine deallocate_particles(iKindParticle)
    integer, intent(inout) :: iKindParticle
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'deallocate_particles'
    !--------------
    call test_start(NameSub, DoTest)
    if(iKindParticle <= 0) call stop_mpi(NameSub//':'//&
         'Particles are not allocated, cannot be deallocated')
    if(iKindParticle /= nKindParticle)call stop_mpi(NameSub//':'//&
         'Only just allocated particles can be deallocated')
    Particle_I(iKindParticle)%nVar        = 0
    Particle_I(iKindParticle)%nIndex      = 0
    Particle_I(iKindParticle)%nParticle   = 0
    Particle_I(iKindParticle)%nParticleMax= 0
    deallocate(Particle_I(iKindParticle)%State_VI)
    deallocate(Particle_I(iKindParticle)%iIndex_II)
    call test_stop(NameSub, DoTest)
  end subroutine deallocate_particles
  !============================================================================
end module ModParticles
!==============================================================================
