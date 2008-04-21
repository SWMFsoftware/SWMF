Module ModHeidiDGCPM
  !\
  ! DGCPM interfacing variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/

	use ModHeidiSize

! Define DGCPM interfacing variables
! Formerly: Common block CDGCPM
        real vthetacells(nthetacells),vphicells(nphicells)
	real vlzcells(nthetacells),vmltcells(nphicells)
        real potdgcpm(nthetacells,nphicells)
        real dendgcpm(nthetacells,nphicells)
        real gridx(nthetacells,nphicells)
        real gridy(nthetacells,nphicells)
        real gridoc(nthetacells,nphicells)

end Module ModHeidiDGCPM
