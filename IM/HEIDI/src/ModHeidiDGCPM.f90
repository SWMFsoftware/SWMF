!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModHeidiDGCPM
  !\
  ! DGCPM interfacing variable definition module for the HEIDI program.
  !/
  
  use ModHeidiSize, ONLY: nthetacells,nphicells
  
  ! Define DGCPM interfacing variables
  real :: vthetacells(nthetacells),vphicells(nphicells)
  real :: vlzcells(nthetacells),vmltcells(nphicells)
  real :: potdgcpm(nthetacells,nphicells)
  real :: dendgcpm(nthetacells,nphicells)
  real :: gridx(nthetacells,nphicells)
  real :: gridy(nthetacells,nphicells)
  real :: gridoc(nthetacells,nphicells)
  
end Module ModHeidiDGCPM
