  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModSize

  implicit none

  private ! except
  public:: &
       nDim, nMomentumBin, nPitchAngleBin, IsPitchAngleAveraged, &
       ROrigin, RSc, &
       Particle_, OriginLat_, OriginLon_, &
       nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle

  ! Dimensionality
  integer, parameter:: nDim = 3

  ! Starting position of field lines in Rs
  real, parameter:: ROrigin = 2.5

  ! Indices of internal cooridnates
  integer, parameter:: &
       Particle_  = 1, & ! Index along  field line
       OriginLon_ = 2, & ! Latitude of field line origin
       OriginLat_ = 3    ! Longitude  of field line origin

      

  ! Min and Max possible index of a particle along a field line,
  ! both are set by Config.pl
  integer, parameter:: iParticleMin = 0
  integer, parameter:: nParticle    = 1000
  integer, parameter:: iParticleMax = iParticleMin + nParticle - 1

  ! angular grid at origin surface
  integer, parameter:: nLon  = 4
  integer, parameter:: nLat  = 4
  integer, parameter:: nNode = nLat * nLon

  ! number of bins in the distribution (see ModGrid);
  integer, parameter:: nMomentumBin   = 100
  integer, parameter:: nPitchAngleBin = 1  

  ! whether to use pitch-angle averaged equations
  logical, parameter:: IsPitchAngleAveraged = nPitchAngleBin == 1

end module ModSize
