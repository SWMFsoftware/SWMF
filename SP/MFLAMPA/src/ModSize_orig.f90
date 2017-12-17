  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModSize

  implicit none

  private ! except
  public:: &
       nDim, nMomentumBin, nPitchAngleBin, IsPitchAngleAveraged, &
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

  ! number of bins in the distribution (see ModGrid);
  integer, parameter:: nMomentumBin   = 100
  integer, parameter:: nPitchAngleBin = 1  

  ! whether to use pitch-angle averaged equations
  logical, parameter:: IsPitchAngleAveraged = nPitchAngleBin == 1

end module SP_ModSize
