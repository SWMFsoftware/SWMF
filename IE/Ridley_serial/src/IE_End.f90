subroutine IE_End(iError)

  use ModIE_Interface

  integer, intent(out) :: iError

  iError = 0
  if (allocated(IEr3_HaveMLTs)) deallocate(IEr3_HaveMLTs, stat = iError)
  if (allocated(IEr3_HaveLats)) deallocate(IEr3_HaveLats, stat = iError)
  if (allocated(IEr3_HavePotential)) &
       deallocate(IEr3_HavePotential, stat = iError)
  if (allocated(IEr3_HaveEFlux)) deallocate(IEr3_HaveEFlux, stat = iError)
  if (allocated(IEr3_HaveAveE)) deallocate(IEr3_HaveAveE, stat = iError)

  stop

end subroutine IE_End
