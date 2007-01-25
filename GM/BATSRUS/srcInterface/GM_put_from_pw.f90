subroutine GM_put_from_pw(Buffer_VI, nFieldLine, nVar)

  implicit none
  character (len=*),parameter :: NameSub='GM_put_from_pw'
  
  integer, intent(in)           :: nVar, nFieldLine
  real, intent(out)             :: Buffer_VI(nVar, nFieldLine)

  !----------------------------------------------------------------------------

  write(*,*) NameSub,': !!!!!!!!!!! Couple GM and PW !!!!!!!!!!!!!!!!'



end subroutine GM_put_from_pw
