!^CFG COPYRIGHT UM
program saha
  use ModStatSum
  implicit NONE
  real :: &
       Nao = 1.00e18,  &  ! 25 , &  ! 25, &     ! cm-3
       vTe,Na 

  real    ::  dTe,dLogN
  integer :: iT,nT=1000,iN
  integer,parameter::nN=5
  real    ::z_I(0:nN),z2_I(0:nN)
  character(LEN=*),parameter,dimension(0:nN)::Separator_I='|'
  character(LEN=*),parameter,dimension(0:nN)::Separator1_I='/'
  logical::IsDegenerated
  dTe = 2.0; dLogN=log(10.0)

  call set_element( 54 )
  call mod_init

  do  iT = 1,nT
     vTe = dTe * iT
     do iN=0,nN
        Na  = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,Na*1000000.0,IsDegenerated)
        Z_I(iN) = z_averaged() 
        Z2_I(iN)= z2_averaged()/Z_I(iN)
        if(IsDegenerated)then
           Z_I(iN)=-1.0
           Z2_I(iN)=-1.0
        end if
     end do
     write(*,'(a,f5.0,a,6(f4.1,a,f4.1,a))')'|',vTe,'|',&
          (Z_I(iN), Separator1_I(iN), Z2_I(iN), Separator_I(iN), iN=0,nN )
  end do
end program saha

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test


