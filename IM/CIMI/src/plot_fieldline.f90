subroutine IM_plot_fieldline(nAlt, iLat, iLon, tSimulation, &
     Length_I,RadialDist_I,Bfield_I)

  use ModIoUnit, ONLY: UnitTmp_

  implicit none
  
  integer, intent(in) :: nAlt,iLat,iLon
  real   , intent(in) :: tSimulation
  real   , intent(in) :: Length_I(nAlt),RadialDist_I(nAlt),Bfield_I(nAlt)

  integer :: i,iTimeOut
  character(len=100) :: NameFieldLine
  !---------------------------------------------------------------------------

  iTimeOut=int(tSimulation)
  write(NameFieldLine,"(a,i2.2,i2.2,a,i8.8,a)") &
       'IM/FieldLine',iLat,iLon,'_',iTimeOut,'.idl'
  open(UnitTmp_,FILE=NameFieldLine)
  write (UnitTmp_,"(a79)") 'RB_line'
  write (UnitTmp_,"(i8,es13.5,3i3)") iTimeOut, tSimulation, 1, 0, 2
  write (UnitTmp_,"(i4)") nAlt
  write (UnitTmp_,"(a79)") 'Length Radius Bfield'
  do i = 1, nAlt  
     write(UnitTmp_,"(3es18.10)") &
          Length_I(i), RadialDist_I(i), Bfield_I(i)
  enddo
  close(UnitTmp_)

end subroutine IM_plot_fieldline
