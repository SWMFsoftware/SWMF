subroutine IM_plot_fieldline(nAlt,iLat,iLon,Length_I,RadialDist_I,Bfield_I)

  use ModIoUnit, ONLY: UnitTmp_
  use ModNumConst, ONLY:cDegToRad
  implicit none
  
  integer, intent(in) :: nAlt,iLat,iLon
  real   , intent(in) :: Length_I(nAlt),RadialDist_I(nAlt),Bfield_I(nAlt)
  integer :: i
  real :: Gamma = 5.0/3.0 
  Character(len=100) :: NameFieldLine
!-----------------------------------------------------------------------------


  write(NameFieldLine,"(a,i2.2,i2.2,a)") &
       'IM/FieldLine',iLat,iLon,'.idl'
  open(UnitTmp_,FILE=NameFieldLine)

  write (UnitTmp_,"(a79)") 'RB_line'
  write (UnitTmp_,"(i8,1pe13.5,3i3)") 1,0.0,1,1,2
  write (UnitTmp_,"(3i4)") nAlt
  write (UnitTmp_,"(100(1pe13.5))") Gamma
  write (UnitTmp_,"(a79)")&
       'Length Radius Bfield g'
  do i=1,nAlt  
     WRITE (UnitTmp_,"(100es18.10)") &
          Length_I(i),RadialDist_I(i),Bfield_I(i)
  enddo
  close(UnitTmp_)

end subroutine IM_plot_fieldline
