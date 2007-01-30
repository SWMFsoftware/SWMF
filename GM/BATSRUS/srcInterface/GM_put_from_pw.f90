subroutine GM_put_from_pw(Buffer_VI, nFieldLine, nVar, Name_V)

  use ModVarIndexes, ONLY: UseMultiSpecies
  use ModUtilities, ONLY: lower_case
  use ModPhysics, ONLY: UnitSi_Rho, UnitSi_RhoU

  implicit none
  character (len=*),parameter :: NameSub='GM_put_from_pw'
  
  integer, intent(in)           :: nVar, nFieldLine
  real, intent(in)              :: Buffer_VI(nVar, nFieldLine)
  character (len=*), intent(in) :: Name_V(nVar)

  integer, parameter :: Theta_=1, Phi_=2, Rho_=1, RhoUb_=2
  integer :: iRhoPwFirst=-1, iRhoPwLast=-1, iUPwFirst=-1, iUPwLast=-1
  logical :: DoInitialize = .true.

  real, allocatable :: CoordPw_VI(:,:), StatePw_VI(:,:)

  integer :: nVarPw=-1, nLinePw=-1
  integer :: iVar
  character (len=40):: NameVar
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoInitialize)then
     nLinePw = nFieldLine
     nVarPw  = nVar
     allocate(CoordPw_VI(2, nFieldLine))
     do iVar = 1, nVar
        NameVar = Name_V(iVar)
        call lower_case(NameVar)
        if(iRhoPwFirst == -1 .and. NameVar(1:1) == 'd') iRhoPwFirst = iVar
        if(NameVar(1:1) == 'd')                         iRhoPwLast  = iVar
        if(iUPwFirst == -1 .and. NameVar(1:1) == 'v')   iUPwFirst   = iVar
        if(NameVar(1:1) == 'v')                         iUPwLast    = iVar
     end do
  end if
  if(nLinePw /= nFieldLine .or. nVarPw /= nVar)then
     write(*,*)NameSub,' ERROR: nLinePw=',nLinePw,'/= nFieldLine=',nFieldLine,&
          ' or nVarPw=',nVarPw,' /= nVar=',nVar
     call CON_stop(NameSub,' nFieldLine or nVar has changed')
  end if

  CoordPw_VI(Theta_,:) = Buffer_VI(Theta_,:)
  CoordPw_VI(Phi_,:)   = Buffer_VI(Phi_,:)

  if(UseMultiSpecies)then
  else
     if(.not. allocated(StatePw_VI)) allocate(StatePw_VI(2, nFieldLine))

     StatePw_VI(Rho_,:)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:), dim=1) &
          / UnitSI_rho
     StatePw_VI(RhoUb_,:)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
          *Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
          / UnitSI_rhoU

     if(DoTestMe)then
        write(*,*)'!!! nVarPw, nLinePw=',nVarPw, nLinePw
        write(*,*)'!!! Buffer_VI=',Buffer_VI
        write(*,*)'!!! CoordPw_VI=',CoordPw_VI
        write(*,*)'!!! StatePw_VI=',StatePw_VI
     end if
  end if

end subroutine GM_put_from_pw
