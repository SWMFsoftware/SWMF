program caxa     ! test of get_ionization()
  use ModSaha    
  use ModIonizPotential
  implicit NONE

  integer          :: i, NZ

  !  real :: Uxe(11)=  &  !<> Xenon
  !       (/ 12.1299, 21.21 ,  32.10,  52.42, 65.31, 89.80, 103.41,   &  !1-7
  !       176.88,  204.09, 225.86, 255.79                        /)   !8-11
  !**************************************  debug perposes, only
  !  real  :: Ube(4)= &    !<> Beryllium 
  !       (/ 15.5  ,  25.5 , 35.5 , 45.5 /)    
  !  Ube(1) = 14.4  !<~> {-dusty -w} is  NOT required          
  !      do i=1,11        ! 
  !         U_I(i)=Uxe(i) !  version for XENON, keep it, temporary
  !      end do           !


  call get_ioniz_potential( 54 , U_I)  ! new, 1-30 materials in use



  vNatomII    = 1.d025         ! [cm^{-3}]  <== set for default test
                               ! deb:> vNatomII    = vNatomII *0.0001
  OCTPOBA:   do i= 5,300,25    ! 5,300,25  <== set for default test

     vTe = 1.0*i               ! [eV]

     call ConcNafter

     write (*,*) "No=", vNatomII," Te=",vTe," Z=",Z," Ne=", vNe," SUi=",CUiTotal
     write (*,*) " ooooooooooooo\ "      
     write (*,*) " ooooooo  ooooo\  Te[K]=", vTe*cEvToK 
     write (*,*) " ooooooooooooooo\  "

     write (*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
  end do OCTPOBA
  !............................................
  ! contains      ! 
  ! <insert  ALL internal procedures HERE>
  ! vDens  <== density, kg/m^{3}
  ! vTmas  <== temperature, K
  ! vnUi   <== max ionization      for given chemical element (nZ)
  ! vmUi   <== ioniz.potent.array  ------/------
  ! nMat   <== simple or compl. mat [0,1]
  !            if complex : nUi & mUi complex
  !                        common vNe, ROTation on mCi<~>mZi

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


