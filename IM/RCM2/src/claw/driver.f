c     =================================================================
      program driver
c
c  Generic driver routine for claw2
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c
      implicit double precision (a-h,o-z)

      real, dimension(-1:182, -1:32) :: qIN,aIN,bIN

      T1=0.d0
      T2=1.d0
      
      call claw2ez(.true., T1, T2,
     &     iterations, 2, 180, 30,
     &     qIN, aIN, bIN)

      stop 
      end
c     =================================================================
