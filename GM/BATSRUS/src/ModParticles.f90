!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModParticles
  use BATL_lib, ONLY: &
       iProc, MaxDim, nDim, check_interpolate_amr_gc, Particle_I,      &
       Particle_I, message_pass_particles, remove_undefined_particles, &
       mark_undefined, check_particle_location, put_particles, trace_particles
  use ModBatlInterface, ONLY: interpolate_grid_amr_gc
  use BATL_particles,    ONLY: BATL_allocate=>allocate_particles
  implicit none
  SAVE
  integer, private, parameter :: nKindParticle = 2
contains
  subroutine allocate_particles
    integer :: iKind
    !---------------
    do iKind = 1, nKindParticle
       call BATL_allocate(iKind)
    end do
  end subroutine allocate_particles
end module ModParticles
