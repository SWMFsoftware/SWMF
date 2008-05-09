program caxa
  use ModSaha    ! CAXA
  implicit NONE

     integer          :: i
      
     real :: Uxe(11)=  &  !<> Xenon
                  (/ 12.1299, 21.21 ,  32.10,  52.42, 65.31, 89.80, 103.41,   & !1-7
                     176.88,  204.09, 225.86, 255.79                        /)   !8-11

     real  :: Ube(4)= &    !<> Beryllium 
                  (/ 15.5  ,  25.5 , 35.5 , 45.5 /)    
  

        Ube(1) = 14.4  !<->NO {-dusty -w}                  

      do i=1,11
         U_I(i)=Uxe(i)
      end do


      NatomII    = 1.d025  !1/cm^3


   OCTPOGA:   do i=5,300,25

      vTe = 1.0*i     ! eV

      call ConcNafter

     write (*,*) "No=", NatomII," Te=",vTe," Z=",Z," Ne=", vNe," SUi=",SUi
     write (*,*) "oooooooooooooo\ "      
     write (*,*) "oooooooo  ooooo\  Te[K]=", vTe*cEvToK 
     write (*,*) "oooooooooooooooo\  "
     
     write (*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

   end do OCTPOGA
!............................................
     ! contains      ! 
     ! <insert  ALL internal procedures HERE>
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

