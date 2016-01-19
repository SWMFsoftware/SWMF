  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModSize

  implicit none

  private ! except
  public:: Id_, Lon_, Lat_, iIdMin, iIdMax

  ! Indices of cooridnates
  integer, parameter:: &
       Id_  = 1, & ! Index along field line
       Lon_ = 2, & ! Longitude
       Lat_ = 3    ! Latitude

  ! Min and Max possible index of a particle along a field line,
  ! both are set by Config.pl
  integer, parameter:: iIdMin = 0
  integer, parameter:: iIdMax = 1000

end module ModSize
