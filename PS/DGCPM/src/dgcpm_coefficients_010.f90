! File name: dgcpm_coefficients_010.f90
!
! Contains: drift and diffusion coefficient definition routines for DGCPM
!	MAGCONV
!
! Last Modified: December 2006, Mike Liemohn
!
! **********************************************************************
!                             MAGCONV
!	Calculates magnetospheric convection strength and puts it
!	into the drift terms 
!	(added to already-calculated corotation and gradient-curvature)
!  Calculation options:
!	IA=1  Kp driven V-S with Maynard and Chen activity dependence
! **********************************************************************


SUBROUTINE MAGCONV(I3)

  use ModSizeDGCPM
  use ModMainDGCPM
  use ModIoDGCPM
  use ModHeidiDGCPM
  use ModConstants

  implicit none

  integer, intent(in) :: i3
  integer j,i

  !  Calculate base convection electric field

  !CC Shielded Volland-Stern

  !!	IF (ABASE(IA+1).EQ.1) THEN

  do J=1,nphicells   ! Fill in DGCPM potentials

!     write(*,*) j, vphicells(j)

     do I=1,nthetacells

!        if (j == nphicells) write(*,*) i, vlzcells(i)

        potdgcpm(i,j)=a*re*vlzcells(i)**(lamgam)*sin(dtor*vphicells(J))
     enddo
  enddo

  write(*,*) "cpcp : ", maxval(potdgcpm)-minval(potdgcpm)

  return

end subroutine magconv


!
! End of subroutine MAGCONV
!


