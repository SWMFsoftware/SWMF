!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModParticles
  use BATL_lib, ONLY: &
       iProc, MaxDim, nDim, check_interpolate_amr_gc, Particle_I,      &
       Particle_I, message_pass_particles, remove_undefined_particles, &
       mark_undefined, check_particle_location, put_particles, trace_particles
  use ModBatlInterface, ONLY: interpolate_grid_amr_gc
  use BATL_size,        ONLY: BatlNKind=>nKindParticle
  implicit none
  SAVE
  integer, private :: nKindParticle = 0
contains
  subroutine allocate_particles(iKindParticle, nVar, nIndex, nParticleMax)
    integer, intent(inout) :: iKindParticle
    integer, intent(in)    :: nVar, nIndex, nParticleMax
    !---------------
    if(iKindParticle > 0) RETURN
    nKindParticle = nKindParticle + 1
    if(nKindParticle > BatlNKind)&
         call stop_mpi('Too many sorts of particles are allocated')
    iKindParticle = nKindParticle
    Particle_I(iKindParticle)%nVar   = nVar 
    Particle_I(iKindParticle)%nIndex = nIndex
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
  end subroutine allocate_particles
end module ModParticles
