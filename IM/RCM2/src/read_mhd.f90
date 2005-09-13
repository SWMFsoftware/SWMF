   PROGRAM Read_file
   USE Rcm_variables
   USE Rcm_io
   IMPLICIT NONE
   REAL :: x,y,lon,lat,lat_bnd
   INTEGER :: i,j
!
!
   OPEN (UNIT=1, FILE='MHDoutput.dat', STATUS='OLD', ACTION='READ')
   READ (1,'(14/)')
   DO i = 1, isize
   DO j = 1, jsize + 1
      READ (1,*) x,y,lon,lat,lat_bnd,xmin(i,j),ymin(i,j),vm(i,j),density(i,j),&
                 temperature(i,j),v(i,j)
      IF (density(i,j) < 0.0 .OR. vm(i,j) < 0.0) THEN
         density(i,j) = -1.0
      ELSE
         density (i,j) = density (i,j) / (vm(i,j)**(-1.5)) / xmass(2) /1.0E+6 ! in cm-3
      END IF
     if (j==24) print*,i,xmin(i,j),density(i,j),vm(i,j)
      temperature (i,j) = temperature (i,j) * 1.38E-23 / 1.6E-19
   END DO
   END DO
   CLOSE (UNIT=1)
   CALL Wrap_around_ghostcells (v, isize, jsize, n_gc)
   CALL Wrap_around_ghostcells (vm, isize, jsize, n_gc)
   CALL Wrap_around_ghostcells (xmin, isize, jsize, n_gc)
   CALL Wrap_around_ghostcells (ymin, isize, jsize, n_gc)
   CALL Wrap_around_ghostcells (density, isize, jsize, n_gc)
   CALL Wrap_around_ghostcells (temperature, isize, jsize, n_gc)
!
!
   CALL Write_array ('rcmxmin_inp', 1, label, ARRAY_2D=xmin, SETUP=.TRUE., ASCI=asci_flag)
   CALL Write_array ('rcmymin_inp', 1, label, ARRAY_2D=ymin, SETUP=.TRUE., ASCI=asci_flag)
   CALL Write_array ('rcmvm_inp', 1, label, ARRAY_2D=vm, SETUP=.TRUE., ASCI=asci_flag)
!
   CALL Write_array ('rcmv_inp', 1, label, ARRAY_2D=v, SETUP=.TRUE., ASCI=asci_flag)
   CALL Write_array ('rcmdensity_inp', 1, label, ARRAY_2D=density, SETUP=.TRUE., ASCI=asci_flag)
   CALL Write_array ('rcmtemperature_inp', 1, label, ARRAY_2D=temperature, SETUP=.TRUE.,ASCI=asci_flag)
!
   STOP 
   END PROGRAM Read_file
