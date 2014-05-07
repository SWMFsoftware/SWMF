!=============================================================================

real function logvar_ionosphere(NameLogvar)
  use ModProcMH,  ONLY: iProc
  use ModIO,      ONLY: write_myname
  use ModIonoPotential, ONLY: nThetaIono, IonoPotential_II
  implicit none
  character (len=*), intent(in) :: NameLogvar
  integer :: nWarn = 0 ! warn multiple times to catch multiple log variables
  !---------------------------------------------------------------------------
  if(.not.allocated(IonoPotential_II))then
     logvar_ionosphere = 0.0
     return
  endif

  select case(NameLogvar)
  case('cpcpn','cpcp_n','cpcp_north','cpcpnorth')
     logvar_ionosphere = &
          maxval(IonoPotential_II(1:(nThetaIono+1)/2,:)) - &
          minval(IonoPotential_II(1:(nThetaIono+1)/2,:))
  case('cpcps','cpcp_s','cpcp_south','cpcpsouth')
     logvar_ionosphere = &
          maxval(IonoPotential_II((nThetaIono+1)/2:nThetaIono,:)) - &
          minval(IonoPotential_II((nThetaIono+1)/2:nThetaIono,:))
  case default
     if(nWarn < 2 .and. iProc==0)then
        call write_myname;
        write(*,'(a)')'WARNING in logvar_ionosphere: unknown NameLogvar='//&
             trim(NameLogvar)//', returning -777.77'
        nWarn = nWarn + 1
     end if
     logvar_ionosphere = -777.77
  end select

end function logvar_ionosphere

