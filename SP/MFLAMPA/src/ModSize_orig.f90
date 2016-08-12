  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModSize

  implicit none

  private ! except
  public:: &
       nDim, nVar, &
       ROrigin, RSc, &
       Particle_, OriginLat_, OriginLon_, &
       R_, Lat_, Lon_, Bx_, By_, Bz_, &
       nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle

  ! Dimensionality
  integer, parameter:: nDim = 3

  ! Number of variables in the state vector
  integer, parameter:: nVar = 6

  ! Starting position of field lines in Rs
  real, parameter:: ROrigin = 2.5

  ! Boundary of the solar corona in Rs
  real, parameter:: RSc =21.0

  ! Indices of internal cooridnates
  integer, parameter:: &
       Particle_  = 1, & ! Index along  field line
       OriginLon_ = 2, & ! Latitude of field line origin
       OriginLat_ = 3    ! Longitude  of field line origin

  ! Indices of variables in state vector
  integer, parameter:: &
       R_   = 1,    & ! Radial coordinate
       Lon_ = 2,    & ! Longitude
       Lat_ = 3,    & ! Latitude
       Bx_  = 4,    & !
       By_  = 5,    & ! Magnetic field
       Bz_  = 6       !
       

  ! Min and Max possible index of a particle along a field line,
  ! both are set by Config.pl
  integer, parameter:: iParticleMin = 0
  integer, parameter:: nParticle    = 1000
  integer, parameter:: iParticleMax = iParticleMin + nParticle - 1

  ! angular grid at origin surface
  integer, parameter:: nLon  = 4
  integer, parameter:: nLat  = 4
  integer, parameter:: nNode = nLat * nLon

end module ModSize
