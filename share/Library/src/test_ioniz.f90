program caxaTe
  use CAXA
  implicit NONE

     integer          :: i
     double precision :: Uxe(11)=  &  !<> Xenon
                  (/ 12.1299, 21.21 ,  32.10,  52.42, 65.31, 89.80, 103.41,   & !1-7
                     176.88,  204.09, 225.86, 255.79                        /)   !8-11

     double precision :: Ube(4)= &    !<> Beryllium 
                  (/ 15.5  ,  25.5 , 35.5 , 45.5 /)    
  

        Ube(1) = 14.4  !dlia PERBAK                    

      do i=1,11
         Ui(i)=Uxe(i)
      end do


      NatomII    = 1.d027  !1/cm^3


  OCTPOGA:   do i=5,300,25

      Telectrons = 1.0*i     ! eV

      call ConcNafter

     
	   write (*,*) "Na=", NatomII, " Te=", Telectrons, " Z=", Z, " Ne=", Neplasma
	   write (*,*) "oooooooooooooo   ooooooooooooooooooo\ "      
	   write (*,*) "oooooooooooo       oooooooo Te[K]= ", Telectrons/Tev  
      write (*,*) "oooooooooooooo   ooooooooooooooooooo/ "
     
      write (*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

   end do OCTPOGA
!............................................
     ! contains      ! 
     ! <insert  ALL internal procedures HERE>
!.............................................
 end program caxaTe 

