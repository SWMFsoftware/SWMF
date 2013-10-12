!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine set_vertical_grid
  use ModCommonPlanet,   ONLY: rLowerBoundary, rPlanet 
  use ModCommonVariables,ONLY: DrBnd, ALTMIN, ALTMAX, ALTD, &
       RAD, RBOUND ,DAREA,  NEXP, CellVolume_C, CellVolumeTop, &
       AR12, AR23, AR12top,AR23top
  use ModPWOM, ONLY: nAlt

  implicit none
  
  integer :: iAlt
  real    :: rLower
  !----------------------------------------------------------------------------
  
  ! DrBnd=altitude step in cm
  ! rLower=lower boundary of the simulation
  ! RAD=radial distance of cell centers?      
  ! RBOUND=radial distance of lower boundary of cell     
  ! ALTD = same as RAD but distance is from surface, not center of planet
  rLower=rLowerBoundary+0.5*DrBnd
  do iAlt=1,nAlt+1
     RBOUND(iAlt)=rLower+(iAlt-1)*DrBnd
  enddo
  do iAlt=1,nAlt
     RAD(iAlt)=0.5*(RBOUND(iAlt)+RBOUND(iAlt+1))
     ALTD(iAlt)=RAD(iAlt)-rPlanet
  enddo
  ALTMIN=ALTD(1)-DrBnd
  ALTMAX=ALTD(nAlt)+DrBnd
  
  !     READ THE EXPONENT OF THE  A(R)=R**NEXP  AREA FUNCTION            C
  
  NEXP=3
  
  ! AR stands for area function. 12 is the lower boundary of the cell
  ! and 23 is the upper boundary
  
  do iAlt=1,nAlt
     DAREA(iAlt)=NEXP/RAD(iAlt)
     AR12(iAlt)=(RBOUND(iAlt)/RAD(iAlt))**NEXP
     AR23(iAlt)=(RBOUND(iAlt+1)/RAD(iAlt))**NEXP
     
     !     For the cell volume in the Rusanov Solver we use a cell volume 
     !     calculated by assuming an area function in the form A=alpha r^3
     !     and then assuming each cell is a truncated cone.
     !     So CellVolume_C is the volume of cell j divided by the crossesction
     !     of cell j, A(j). 
     
     CellVolume_C(iAlt)=1.0/3.0 * DrBnd *&
          ( Ar12(iAlt) + Ar23(iAlt) + ( Ar12(iAlt)*Ar23(iAlt) )**0.5 ) 
     
  enddo
  ! area and volume for ghost cell
  AR12top(1)=((rLower+nAlt*DrBnd)/(rLower+nAlt*DrBnd+DrBnd*0.5))**NEXP
  AR23top(1)=((rLower+(nAlt+1.0)*DrBnd)/(rLower+(nAlt+1.0)*DrBnd+DrBnd*0.5))**NEXP
  AR12top(2)= AR23top(1)
  AR23top(2)=&
       ((rLower+(nAlt+2.0)*DrBnd)/(rLower+(nAlt+2.0)*DrBnd+DrBnd*0.5))**NEXP
  
  CellVolumeTop(1) = 1.0/3.0 * DrBnd *&
       ( Ar12top(1) + Ar23top(1) + ( Ar12top(1)*Ar23top(1) )**0.5 ) 
  CellVolumeTop(2) = 1.0/3.0 * DrBnd * &
       ( Ar12top(2) + Ar23top(2) + ( Ar12top(2)*Ar23top(2) )**0.5 ) 
  
end subroutine set_vertical_grid
