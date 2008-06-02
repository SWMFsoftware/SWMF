program caxa     ! test of get_ionization()
  use ModSaha    
  use ModIonizPotential

  implicit NONE

  integer          :: i 

  call get_ioniz_potential( 54 , U_I)  ! new, 1-30 materials in use

  vNatomII    = 1.d025!19 ! 25!<=std TEST         ! [cm^{-3}]  <== set for default test
                               ! deb:> vNatomII    = vNatomII *0.0001
  OCTPOBA:   do i= 5,300,25    ! 5,300,25  <== set for default test
     vTe = 1.0*i               ! [eV]

     call ConcNafter

     write (*,*) "No=", vNatomII," Te=",vTe," Z=",Z," Ne=", vNe," SUi=",CUiTotal
     write (*,*) " ooooooooooooo\ "      
     write (*,*) " ooooooo  ooooo\  Te[K]=", vTe*cEvToK 
     write (*,*) " ooooooooooooooo\  "

     write (*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
  end do OCTPOBA
  !.............................................
end program caxa





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


