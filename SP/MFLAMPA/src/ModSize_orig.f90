  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModSize

  implicit none

  private ! except
  public:: &
       nDim, nMomentum, nPitchAngle, IsPitchAngleAveraged, &
       Particle_,  OriginLon_, OriginLat_, &
       nLon, nLat, nNode, nParticleMax

  ! Dimensionality
  integer, parameter:: nDim = 3

  ! Indices of internal cooridnates
  integer, parameter:: &
       Particle_  = 1, & ! Index along  field line
       OriginLon_ = 2, & ! Longitude of field line origin
       OriginLat_ = 3    ! Latitude  of field line origin

      

  ! Max possible index of a particle on a line set by Config.pl
  integer, parameter:: nParticleMax = 1000

  ! angular grid at origin surface
  integer, parameter:: nLon  = 4
  integer, parameter:: nLat  = 4
  integer, parameter:: nNode = nLon*nLat

  ! number of points along the phase coords (see ModAdvance);
  integer, parameter:: nMomentum   = 100
  integer, parameter:: nPitchAngle = 1  

  ! whether to use pitch-angle averaged equations
  logical, parameter:: IsPitchAngleAveraged = nPitchAngle == 1

end module SP_ModSize
