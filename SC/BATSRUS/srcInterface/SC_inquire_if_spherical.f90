subroutine SC_inquire_if_spherical(IsSphericalSc)
  use SC_ModGeometry,ONLY:UseCovariant,TypeGeometry
  implicit none
  logical,intent(out)::IsSphericalSc
  
  character(len=*), parameter :: NameSub='SC_inquire_if_spherical'
  if(UseCovariant)then
     IsSphericalSc=index(TypeGeometry,'spherical')>0
  else
     IsSphericalSc=.false.
  end if
end subroutine SC_inquire_if_spherical
!===================================================================!
