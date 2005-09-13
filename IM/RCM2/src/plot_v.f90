      PROGRAM Plot_v
      USE Rcm_variables
      USE Rcm_io
!     use f90_unix
      IMPLICIT NONE
!     include 'mpif.h'
      integer :: pgopen, add_cond, minind, maxind, iv, color_table, dev_type_index
      INTEGER (iprec), PARAMETER :: wsize = jsize, &
                                    j_shift = (jsize)/2, &
                                    npanel_size=6, ncut_size=3
      INTEGER (iprec) :: I, j, n_rec, n_rec_1, n_rec_2, itime_init,&
                         n_rec_depth, j_midnight, im, Iswitch, ncut, npanel, &
                         Do_page, i_cut (ncut_size), j_cut( ncut_size),      &
                         n_levels
      REAL (rprec), TARGET :: v_rlb (isize, wsize),                          &
                      aloct_rlb (isize, wsize), bndloc_rlb (wsize),          &
                      e_east_rlb (isize, wsize), e_south_rlb (isize,wsize),  &
                      e_south (isize,jsize), e_east (isize,jsize),           &
                      x_mark_text, y_mark_text, max_y, max_x, min_x,         &
                      x1_vsize, x2_vsize, y1_vsize, y2_vsize,                &
                      colat_low, colat_high, x_b, y_b, R_max,                &
                      colat_cut (ncut_size), mlt_cut (ncut_size),            &
                      aloct_cut (ncut_size), lat_cut (ncut_size),            &
                      x_pt, y_pt, x_prime, y_prime, scale, &
                      bright, contra, a1, a2, av, tr(6), xpoly(4), ypoly(4), &
                      eeta_to_plot (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                      x_plot(2), y_plot(2), x_wedge(2), y_wedge(2), r, p
       REAL (rprec) :: Gntrp_2d, Gntrp_2d_ang, Bndy
      REAL (rprec), PARAMETER :: x_width   = 3.0,   &  ! 3 reg, 3.6 for poster
                                 x_spacing = 0.2,   &  ! 1.5 reg., 0.2 for poster
                                 x_start   = 1.5,   &  ! 1.5 reg., 1.2 for poster
                                 y_height  = 1.5,   &  ! 1.5 reg., 1.7 for poster
                                 y_spacing = 0.2        ! all in inches
      REAL (rprec), DIMENSION (npanel_size) :: &
                      x1plot = (/ x_start, x_start, x_start,            &
                                  x_start+x_width+x_spacing,    &
                                  x_start+x_width+x_spacing,    &
                                  x_start+x_width+x_spacing /), &
                      x2plot = (/ x_start + x_width, &
                                  x_start + x_width, &
                                  x_start + x_width, &
                                  x_start + 2*x_width+x_spacing, &
                                  x_start + 2*x_width+x_spacing, &
                                  x_start + 2*x_width+x_spacing /),&
                      y1plot  = (/ 1.0, 1.0 + y_height+y_spacing, &
                                   1.0 + 2*(y_height+y_spacing), &
                                   1.0, 1.0 + y_height+y_spacing, &
                                   1.0 + 2*(y_height+y_spacing) /),&
                      y2plot  = (/ 1.0+y_height, 1.0 + 2*y_height+y_spacing, &
                                   1.0 + 2*(y_height+y_spacing)+y_height, &
                                   1.0+y_height, 1.0 + 2*y_height+y_spacing, &
                                   1.0 + 2*(y_height+y_spacing)+y_height /)
      REAL (rprec), ALLOCATABLE :: data (:,:), xdata (:), workspace (:), v_levels (:)
      REAL (rprec), POINTER :: slice_x (:), slice_y (:)
!
      CHARACTER (LEN=1) :: choice
      CHARACTER (LEN=10):: time_string
      CHARACTER (LEN=5) :: v_label
      CHARACTER (LEN=62) :: info_string      
      CHARACTER (LEN=80) :: lats_char='', dev_name, dev_type, plot_labeL
      CHARACTER (LEN=8) :: char_cut (ncut_size)
      CHARACTER (LEN=5), DIMENSION (7), PARAMETER :: page_name = &
         (/ 'v_mlt', 'v_eqt', 'v_ion', 'v_ica', 'v_icp', 'e_mlt', 'e_lat' /)
      CHARACTER (LEN=6) :: ut_string
      CHARACTER (LEN=11), PARAMETER :: date_string_all(4) = (/ &
!        'Oct-18-1998', 'Oct-19-1998', 'Oct-20-1998', 'Oct-21-1998' /)
         '', '', '', '' /)
      CHARACTER (LEN=11) :: date_string
!

      INTERFACE
          SUBROUTINE Plot_v_io (visble, x, y, z)
          USE Rcm_variables, ONLY : iprec, rprec
          IMPLICIT NONE
          REAL (rprec), INTENT (IN) :: x, y, z
          INTEGER (iprec), INTENT (IN) :: visble
          END SUBROUTINE Plot_v_io
      END INTERFACE

      INTERFACE
          SUBROUTINE Plot_v_equat (visble, x, y, z)
          USE Rcm_variables, ONLY : iprec, rprec
          IMPLICIT NONE
          REAL (rprec), INTENT (IN) :: x, y, z
          INTEGER (iprec), INTENT (IN) :: visble
          END SUBROUTINE Plot_v_equat
      END INTERFACE

      INTERFACE
          SUBROUTINE Plot_v_io_cart (visble, x, y, z)
          USE Rcm_variables, ONLY : iprec, rprec
          IMPLICIT NONE
          REAL (rprec), INTENT (IN) :: x, y, z
          INTEGER (iprec), INTENT (IN) :: visble
          END SUBROUTINE Plot_v_io_cart
      END INTERFACE
!
!
!
CALL Getarg (1, info_string)
READ (info_string, '(I5)') Do_page
CALL Getarg (2, info_string)
READ (info_string,'(I5)') n_rec
CALL Getarg (3, dev_type)
CALL Getarg (4, info_string)
READ (info_string,'(I2)') add_cond
CALL Getarg (5, info_string)
READ (info_string,'(I3)') color_table
CALL Getarg (6, plot_label)
plot_label=TRIM(plot_label)
!
!
      IF (dev_type(1:LEN_TRIM(dev_type)) == 'GIF' .OR. &
          dev_type(1:LEN_TRIM(dev_type)) == 'gif') THEN
          dev_type_index = 1
      ELSE
          dev_type_index = 0
      END IF
!
!
      Iswitch = 1 ! 1 is forward derivs, 2 is central derivs
      lat_cut = (/ 15.0, 40.0, 55.0 /)
      mlt_cut = (/ 18.0, 24.0, 06.0 /)
!
!
      CALL Read_grid ()
      CALL Read_plasma ()
!
!
      j_midnight = -999
      Search_loop_1: DO j = 1, jsize
         IF (ABS(aloct(1,j)-pi) < 10.0*EPSILON(1.0_rprec)) THEN
            j_midnight = j
            EXIT Search_loop_1
         END IF
      END DO Search_loop_1
      IF (j_midnight == -999) STOP 'NO MIDNIGHT'
!
!
      CALL Read_array ('rcmbndloc', n_rec, label, ARRAY_1D = bndloc)
imin_j = CEILING(bndloc)
      CALL Read_array ('rcmv',      n_rec, label, ARRAY_2D = v)
      CALL Read_array ('rcmxmin',   n_rec, label, ARRAY_2D = xmin)
      CALL Read_array ('rcmymin',   n_rec, label, ARRAY_2D = ymin)
      CALL Read_array ('rcmpedpsi',   n_rec, label, ARRAY_2D = pedpsi)
!do j = 1, jsize
!xmin(1:imin_j(j)-1,j) = xmin(imin_j(j),j)
!ymin(1:imin_j(j)-1,j) = ymin(imin_j(j),j)
!print*,j,imin_j(j),xmin(imin_j(j),j),ymin(imin_j(j),j)
!end do
!pause
!call outp (birk, 1, isize,1,1,jsize,1,0.0,label%intg,'birk',101,80)
      IF (Do_page >4) CALL Compute_efield (e_south, e_east, Iswitch)
!
      date_string = date_string_all (1+FLOOR(label%intg(6)/3600.0/24.00))

      WRITE (time_string,'(I2.2,A1,I2.2,A1,I2.2)') &
            label%intg(3),':',label%intg(4),':',label%intg(5)
      ut_string = time_string(1:2)//time_string(4:5)//time_string(7:8)
      WRITE (time_string(1:2),'(I2.2)') MODULO(label%intg(3),24)
      WRITE (*,'(A,A10,A)',ADVANCE='NO') 'TIME IS: ',time_string
!
!
!
     IF (Do_page > 4) THEN
!
!    Relabel RCM arrays in local time so that j=1 is LT=0,
!    j=wsize is LT=24:
!
      stop 'disabled'
!
      END IF
!     Construct name of plotting device:
!      -- if Xwindows window, pass only device type to PGPLOT
!      -- if a GIF file, make a file with .gif extension
!      -- if any of 4 types of postscript, make a .ps file
!
!
      IF (dev_type(1:LEN_TRIM(dev_type)) == 'xwin' .OR. &
          dev_type(1:LEN_TRIM(dev_type)) == 'XWIN' ) THEN
         dev_name = '/'//dev_type(1:LEN_TRIM(dev_type))
      ELSE IF (dev_type(1:LEN_TRIM(dev_type)) == 'GIF' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'gif' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'PS'  .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'ps' ) THEN
         dev_name = page_name(do_page)//'_'//ut_string//'.'// &
                    dev_type(1:LEN_TRIM(dev_type))//'/'//&
                    dev_type(1:LEN_TRIM(dev_type))
      ELSE IF (dev_type(1:LEN_TRIM(dev_type)) == 'VPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'vps' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'CPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'cps' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'VCPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'vcps' ) THEN
         dev_name = page_name(do_page)//'_'//ut_string//'.'// &
                    'ps'//'/'//&
                    dev_type(1:LEN_TRIM(dev_type))
      ELSE
         STOP 'DEVICE-FILENAME NOT IMPLEMENTED'
      END IF
      IF (Pgopen(dev_name) <= 0) STOP 'could not open graphics'
      CALL Pgslw (2)  ! make lines thicker
      CALL Pgscf (2)  ! roman font will be used
!
      CALL Pgqvsz (1, x1_vsize, x2_vsize, y1_vsize, y2_vsize)
      IF (.NOT.(ABS(y2_vsize-7.8)<0.2.AND.ABS(x2_vsize-10.5)<0.2)) THEN
!
!        this is not letter postscript page, so re-scale:
!
         WRITE (6,'(A,F5.1,A,F5.1)') 'RESCALING FROM 11x8.5 PS TO', x2_vsize, ' BY', y2_vsize
         x1plot = x1plot * (x2_vsize/11.0)
         x2plot = x2plot * (x2_vsize/11.0)
         y1plot = y1plot * (y2_vsize/8.5)
         y2plot = y2plot * (y2_vsize/8.5)
      END IF
!
   IF (Do_page == 1) THEN !potential vs. MLT for each I-value (overlap curves):
!
      CALL Pgpage ()
      CALL Pgvstd ()
      CALL Pgswin (0.0, 24.0, -50.0, +50.0)
      CALL Pgbox ('ABCNTS', 3.0, 3, 'ABCNTS', 0.0_rprec, 0)
      CALL Pglab ('LT, hrs', 'V', 'RCM')
      DO i = INT ( MINVAL (bndloc(:))), isize, 1
         slice_x => aloct_rlb (i, :)
         slice_y => v_rlb     (i, :)
         slice_y = slice_y /1000.0_rprec
         CALL Pgline ( SIZE(slice_x), slice_x, slice_y )
      END DO
!
!
   ELSE IF (Do_page ==2) THEN !contours of potential in the equatorial plane:
!
      IF (dev_type == 'GIF' .OR. dev_type == 'gif'.OR.&
          dev_type == 'XWIN' .OR. dev_type == 'xwin' ) THEN
         CALL Pgpap (6.0, 1.0) !set view surface to square 6.0 inches
      END IF
      CALL Setup_contour_levels ()
      CALL Pgpage ()
!     IF (dev_type == 'GIF' .OR. dev_type == 'gif' .OR.&
!         dev_type == 'XWIN' .OR. dev_type == 'xwin') THEN
         CALL Pgvsiz (0.8, 5.3, 0.8, 5.3)
!     ELSE
!        CALL Pgvsiz (1.5, 7.0, 1.5, 7.0)
!     END IF
!     CALL Pgwnad (+15.0, -60.0, +15.0, -15.0) !window in eq. plane GSM
!     CALL Pgwnad (+15.0, -25.0, +20.0, -20.0) !window in eq. plane GSM
      CALL Pgwnad (+8.0, -12.0, +10.0, -10.0) !window in eq. plane GSM
      CALL Pgbox ('abctsn', 0.0, 0, 'abctsnv', 0.0, 0)
      CALL Pgmtxt ('B', 3.0, 0.5, 0.5, 'X\dGSM\u')
      CALL Pgmtxt ('L', 3.0, 0.5, 0.5, 'Y\dGSM\u')
      WRITE (info_string,'(I5.5)') label%intg(1)
      CALL Pgmtxt ('T', 1.0, 0.0, 0.0, 'RCM: UT='//time_string//' '//date_string)
      CALL Pgmtxt ('T', 2.2, 0.0, 0.0, plot_label)
      CALL Pgsave ()
      CALL Pgstbg (15)
      CALL Pgmtxt ('T', 1.0, 1.0, 1.0,'   \gF-\gF\dcorot\u')
      CALL Pgunsa ()
!     CALL Pgmtxt ('T', 2.5, 0.0, 0.0, label%char)
!
!
!     Plot color-coded conductances if requested:
!
      IF (add_cond /= 0) THEN
!
         eeta_to_plot = 0.5*pedpsi
         CALL Pgsave ()
         bright = 0.5
         contra = 1.0
         CALL Palett (color_table, contra, bright)
         CALL Pgqcir (minind, maxind)
         a2 = 10.0
         a1 = 0.01*a2
         IF (a1 >= a2) STOP 'STOP: A1 >= A2'
         IF (minind >= maxind) STOP 'NOT ENOUGH COLORS'
         DO j = 1, jsize
         DO i = 1, isize
            IF (i < Bndy(bndloc,jsize, REAL(j,rprec))) CYCLE
            av = eeta_to_plot (i,j)
            IF (a2 > a1) THEN
               av = MIN (a2, MAX(a1,av))
            ELSE
               av = MIN (a1, MAX(a2,av))
            END IF
            iv = NINT ( (minind*(a2-av) + maxind*(av-a1)) / (a2-a1))
!        
!           Define a polygon:
!
            r = Gntrp_2d (xmin, isize, jsize, REAL(i-0.5,rprec), REAL (j-0.5, rprec))
            p = Gntrp_2d (ymin, isize, jsize, REAL(i-0.5,rprec), REAL (j-0.5, rprec))
            xpoly (1) = r 
            ypoly (1) = p
!
            r = Gntrp_2d (xmin, isize, jsize, REAL(i-0.5,rprec), REAL (j+0.5, rprec))
            p = Gntrp_2d (ymin, isize, jsize, REAL(i-0.5,rprec), REAL (j+0.5, rprec))
            xpoly (2) = r 
            ypoly (2) = p
!
            r = Gntrp_2d (xmin, isize, jsize, REAL(i+0.5,rprec), REAL (j+0.5, rprec))
            p = Gntrp_2d (ymin, isize, jsize, REAL(i+0.5,rprec), REAL (j+0.5, rprec))
            xpoly (3) = r 
            ypoly (3) = p
!
            r = Gntrp_2d (xmin, isize, jsize, REAL(i+0.5,rprec), REAL (j-0.5, rprec))
            p = Gntrp_2d (ymin, isize, jsize, REAL(i+0.5,rprec), REAL (j-0.5, rprec))
            xpoly (4) = r 
            ypoly (4) = p
!
            CALL Pgsci (iv)
            CALL Pgpoly (4, xpoly, ypoly)
         END DO
         END DO
         CALL Pgunsa ()
!        CALL Pgvsiz (x_wedge(1), x_wedge(2), y_wedge(1), y_wedge(2))
         tr = (/ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
!        a1 = eeta_min / eeta_max * 100.0
!        a2 = 100.0
!        DO i = 1, SIZE(wedge)
!           wedge(i) =  (SIZE(wedge)*a1-a2)/(SIZE(wedge)-1) + &
!                       (a2-a1)/(SIZE(wedge)-1)*REAL(i,rprec)
!        END DO
!        CALL Pgimag (wedge, SIZE(wedge), 1, 1, SIZE(wedge), 1, 1, a1, a2, tr)
         CALL Pgsave ()
         CALL Pgsch (1.5)
         CALL Pgwedg ('RI',  0.5, 3.0, a1, a2, '\gS, S')
         CALL Pgunsa ()
      END IF
!
!
      CALL Pgsch (0.5_rprec)
      DO i = 1, n_levels
         WRITE (v_label,'(F5.1)') v_levels(i)/1000.0
         CALL Pgconx (v(1:isize,1:jsize+1), isize, jsize+1, 1, isize, 1, jsize+1, v_levels(i), -1, Plot_v_equat)
      END DO
      CALL pgsch (1.0_rprec)
!
!     Write contour step size and max/min values on the plot:
!
      CALL Pgsave ()
      CALL Pgscf (1)
      CALL Pgsch (0.8)
      CALL Pgstbg (0)
      WRITE (info_string,'(SP,(A,F6.1,A3))') '\gDV=',(v_levels(2)-v_levels(1))/1000.0, ' kV'
      CALL Pgmtxt ('T', -1.5, 0.7, 0.0, info_string)
      WRITE (info_string,'(SP,(A,F6.1,A3))') 'V\dmax\u=', MAXVAL(v)/1000.0, ' kV'
      CALL Pgmtxt ('T', -2.5, 0.7, 0.0, info_string)
      WRITE (info_string,'(SP,(A,F6.1,A3))') 'V\dmin\u=', MINVAL(v)/1000.0, ' kV'
      CALL Pgmtxt ('T', -3.5, 0.7, 0.0, info_string)
      CALL Pgunsa ()
!
!     Draw earth with white to erase possible lines inside inner boundary:
      CALL Pgsfs(1)
      CALL Pgsci(0)  ! set default color to "white"
      CALL Pgmove (xmin(isize,1), ymin(isize,1))
      DO j = 1, jsize+1
         CALL Pgdraw (xmin(isize,j),ymin(isize,j))
      END DO
!     CALL Pgcirc (0.0, 0.0, 0.99)
      CALL Pgsci(1)  ! set default color back to "black"
      CALL Pgsfs (2)
      CALL Pgmove (xmin(isize,1), ymin(isize,1))
      DO j = 1, jsize
         CALL Pgdraw (xmin(isize,j)-0.01*cos(aloct(isize,j)), &
                      ymin(isize,j)-0.01*sin(aloct(isize,j)))
      END DO
!     CALL Pgcirc (0.0, 0.0, 0.99)
!
!
!     Draw RCM high- and low-lat. boundaries as thick lines:
      CALL Pgsave ()
      CALL Pgslw (6)
      j = 1
!     x_b = xmin (imin_j(j),j)
!     y_b = ymin (imin_j(j),j)
      x_b = Gntrp_2d (xmin, isize, jsize, bndloc(j), REAL(j,rprec))
      y_b = Gntrp_2d (ymin, isize, jsize, bndloc(j), REAL(j,rprec))
      CALL Pgmove (x_b, y_b)
      DO j = 2, jsize+1
      x_b = Gntrp_2d (xmin, isize, jsize, bndloc(j), REAL(j,rprec))
      y_b = Gntrp_2d (ymin, isize, jsize, bndloc(j), REAL(j,rprec))
!        x_b = xmin (imin_j(j),j)
!        y_b = ymin (imin_j(j),j)
         CALL Pgdraw (x_b, y_b)
      END DO
      CALL Pgunsa ()
!
!
   ELSE IF (Do_page == 3) THEN !contours of V in the ionosphere:
!
!     IF (dev_type_index == 1 ) THEN
         CALL Pgpap (6.0, 1.0) !set view surface to square 6.0 inches
!     END IF
      CALL Setup_contour_levels ()
      CALL Pgpage ()
!     IF (dev_type_index == 1) THEN
         x_plot = (/ 0.8, 5.3 /)
         y_plot = (/ 0.8, 5.3 /)
         x_wedge = (/ 5.4, 5.7 /)
         y_wedge = (/ 0.8, 5.3 /)
!     ELSE
!        x_plot = (/ 1.5, 6.0 /)
!        y_plot = (/ 1.5, 6.0 /)
!        x_wedge = (/ 6.1, 6.4 /)
!        y_wedge = (/ 1.5, 6.0 /)
!     END IF
      CALL Pgvsiz (x_plot(1), x_plot(2), y_plot(1), y_plot(2))
      R_max = 45.0
      CALL Pgwnad (-R_max, +R_max, -R_max, +R_max)
!     CALL Pgswin (-R_max, +R_max, -R_max, +R_max)
!
!
!     Plot color-coded conductances if requested:
!
      IF (add_cond /= 0) THEN
!
         eeta_to_plot = 0.5*pedpsi
         CALL Pgsave ()
         bright = 0.5
         contra = 1.0
         CALL Palett (color_table, contra, bright)
         CALL Pgqcir (minind, maxind)
         a2 = 15.0
         a1 = 0.01*a2
         IF (a1 >= a2) STOP 'STOP: A1 >= A2'
         IF (minind >= maxind) STOP 'NOT ENOUGH COLORS'
         DO j = 1, jsize
         DO i = 1, isize
            IF (i < Bndy(bndloc,jsize, REAL(j,rprec))) CYCLE
            av = eeta_to_plot (i,j)
            IF (a2 > a1) THEN
               av = MIN (a2, MAX(a1,av))
            ELSE
               av = MIN (a1, MAX(a2,av))
            END IF
            iv = NINT ( (minind*(a2-av) + maxind*(av-a1)) / (a2-a1))
!        
!           Define a polygon:
!
            r = RTD * Gntrp_2d (colat, isize, jsize, REAL(i-0.5,rprec), REAL (j-0.5, rprec))
            p =       Gntrp_2d (aloct, isize, jsize, REAL(i-0.5,rprec), REAL (j-0.5, rprec))
            xpoly (1) = r * COS(p+pi)
            ypoly (1) = r * SIN(p+pi)
!
            r = RTD * Gntrp_2d (colat, isize, jsize, REAL(i-0.5,rprec), REAL (j+0.5, rprec))
            p =       Gntrp_2d (aloct, isize, jsize, REAL(i-0.5,rprec), REAL (j+0.5, rprec))
            xpoly (2) = r * COS(p+pi)
            ypoly (2) = r * SIN(p+pi)
!
            r = RTD * Gntrp_2d (colat, isize, jsize, REAL(i+0.5,rprec), REAL (j+0.5, rprec))
            p =       Gntrp_2d (aloct, isize, jsize, REAL(i+0.5,rprec), REAL (j+0.5, rprec))
            xpoly (3) = r * COS(p+pi)
            ypoly (3) = r * SIN(p+pi)
!
            r = RTD * Gntrp_2d (colat, isize, jsize, REAL(i+0.5,rprec), REAL (j-0.5, rprec))
            p =       Gntrp_2d (aloct, isize, jsize, REAL(i+0.5,rprec), REAL (j-0.5, rprec))
            xpoly (4) = r * COS(p+pi)
            ypoly (4) = r * SIN(p+pi)
!
            CALL Pgsci (iv)
            CALL Pgpoly (4, xpoly, ypoly)
         END DO
         END DO
         CALL Pgunsa ()
!        CALL Pgvsiz (x_wedge(1), x_wedge(2), y_wedge(1), y_wedge(2))
         tr = (/ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
!        a1 = eeta_min / eeta_max * 100.0
!        a2 = 100.0
!        DO i = 1, SIZE(wedge)
!           wedge(i) =  (SIZE(wedge)*a1-a2)/(SIZE(wedge)-1) + &
!                       (a2-a1)/(SIZE(wedge)-1)*REAL(i,rprec)
!        END DO
!        CALL Pgimag (wedge, SIZE(wedge), 1, 1, SIZE(wedge), 1, 1, a1, a2, tr)
         CALL Pgsave ()
         CALL Pgsch (1.5)
         CALL Pgwedg ('RI',  0.5, 3.0, a1, a2, '\gS, S')
         CALL Pgunsa ()
      END IF
!
!
!     Write contour step size and max/min values on the plot:
      WRITE (info_string,'(I5.5)') label%intg(1)
      CALL Pgmtxt ('T', 1.0, 0.0, 0.0, 'RCM: UT='//time_string//' '//date_string)
      CALL Pgmtxt ('T', 2.1, 0.0, 0.0, plot_label)
      CALL Pgsave ()
      CALL Pgstbg (15)
      CALL Pgmtxt ('T', 1.0, 1.0, 1.0,'   \gF\dIONOSPHERE\u')
      CALL Pgunsa ()
!     CALL Pgmtxt ('T', 2.5, 0.0, 0.0, label%char)
!
!
      CALL Pgsave ()
      CALL Pgscf (1)
      CALL Pgsch (0.8)
      CALL Pgsci (0)
      CALL Pgstbg (1)
      WRITE (info_string,'(SP,(A,F6.1,A3))') '\gDV=',(v_levels(2)-v_levels(1))/1000.0, ' kV'
      CALL Pgmtxt ('T', -1.5, 0.7, 0.0, info_string)
      WRITE (info_string,'(SP,(A,F6.1,A3))') 'V\dmax\u=', MAXVAL(v)/1000.0, ' kV'
      CALL Pgmtxt ('T', -2.8, 0.7, 0.0, info_string)
      WRITE (info_string,'(SP,(A,F6.1,A3))') 'V\dmin\u=', MINVAL(v)/1000.0, ' kV'
      CALL Pgmtxt ('T', -4.1, 0.7, 0.0, info_string)
      CALL Pgsci (1)
      CALL Pgstbg (0)
      CALL Pgunsa ()
!
!     IF (dev_type_index /= 1) CALL Pgsci (0)
      DO i = 1, n_levels
         CALL Pgconx (v(1:isize,1:jsize), isize, jsize, 1, isize, 1, jsize, v_levels(i), +1, Plot_v_io)
      END DO
!     IF (dev_type_index /= 1) CALL Pgsci (1)
!
!
!     CALL Pgsave ()
         CALL Pgstbg (15)
         CALL Pgptxt (-R_max*0.95, 0.0,        +90.0, 0.5, 'MLT=12')
         CALL Pgptxt (+R_max*0.95, 0.0,        -90.0, 0.5, 'MLT=24')
         CALL Pgptxt (0.0,    -R_max*0.95,   0.0, 0.5, 'MLT=18')
!     CALL Pgunsa ()
!
!
!     Plot three circles:
      CALL Pgsave ()
         CALL Pgsfs (2) !fill-area style is outline
         CALL Pgsls (4) !dotted lines
         CALL Pgslw (1)
         CALL Pgcirc (0.0, 0.0, 10.0)
         CALL Pgcirc (0.0, 0.0, 20.0)
         CALL Pgcirc (0.0, 0.0, 30.0)
         CALL Pgptxt (-10.0/SQRT(2.0), 10.0/SQRT(2.0), 45.0, 0.5, '80')
         CALL Pgptxt (-20.0/SQRT(2.0), 20.0/SQRT(2.0), 45.0, 0.5, '70')
         CALL Pgptxt (-30.0/SQRT(2.0), 30.0/SQRT(2.0), 45.0, 0.5, '60')
      CALL Pgunsa ()
!
!     Plot the boundary:
!     x_b = Gntrp (colat, bndloc (jwrap), REAL (jwrap), 0) * RTD*COS (aloct(1,jwrap)+pi)
!     y_b = Gntrp (colat, bndloc(jwrap), REAL (jwrap), 0) * RTD*SIN (aloct(1,jwrap)+pi)
!     CALL Pgmove (x_b, y_b)
!     CALL Pgsave ()
!     CALL Pgslw (7)
!     DO j = jwrap + 1, jsize
!        x_b = Gntrp (colat, bndloc(j), REAL (j), 0) * RTD*COS (aloct(1,j)+pi)
!        y_b = Gntrp (colat, bndloc(j), REAL (j), 0) * RTD*SIN (aloct(1,j)+pi)
!        CALL Pgdraw (x_b, y_b)
!     END DO
!     CALL Pgunsa ()
!
!
   ELSE IF (Do_page == 4) THEN
!
!    Plot electrostatic potential contours in the ionosphere
!    in the rectangular coordinates.
!
     CALL Setup_contour_levels ()
!
     CALL Pgpage ()
     CALL Pgsave ()
     CALL Pgsch (0.7)
     CALL Pgvsiz (2.0, 6.0, 2.0, 5.0)
     CALL Pgswin (0.0, 24.0, 79.0, 10.0)
     CALL Pgbox ('ABCNTS', 3.0, 3, 'ABCNTSV', 10.0, 5)
     CALL Pgmtxt ('B', 2.5, 0.5, 0.5, 'MLT, hrs')
     CALL Pgmtxt ('L', 4.0, 0.5, 0.5, '\gh, deg')
     DO i = 1, n_levels
        CALL Pgconx (v, isize, jsize, 1, isize, 1 , j_midnight-1, &
                     v_levels (i), -1, plot_v_io_cart)
     END DO
     DO i = 1, n_levels
        CALL Pgconx (v, isize, jsize,  1, isize, j_midnight , jsize-1, &
                     v_levels (i), -1, plot_v_io_cart)
     END DO
!
     CALL Pgsave ()
     CALL Pgslw (6)
!
     j = 1
     y_b = RTD*Gntrp_2d (colat, isize,jsize,bndloc(j), REAL (j))
     x_b =     Gntrp_2d (aloct, isize,jsize,bndloc(j), REAL (j))
       x_b = MODULO (x_b*RTH+12.0, 24.0)
     CALL Pgmove (x_b,y_b)
     DO j = 1, j_midnight
       y_b = RTD*Gntrp_2d (colat, isize, jsize, bndloc(j), REAL(j)) 
       x_b =     Gntrp_2d (aloct, isize, jsize, bndloc(j), REAL(j))
         x_b = MODULO (x_b*RTH+12.0, 24.0)
         IF (j == j_midnight) x_b = 24.0
       CALL Pgdraw (x_b, y_b)
     END DO
!
     j = j_midnight
     y_b = RTD*Gntrp_2d (colat, isize, jsize, bndloc(j), REAL (j))
     x_b = Gntrp_2d_ang (aloct, isize, jsize, bndloc(j), REAL (j))
       x_b = MODULO (x_b*RTH+12.0, 24.0)
     CALL Pgmove (x_b,y_b)
     DO j = j_midnight+1, jsize
       y_b = RTD*Gntrp_2d (colat, isize, jsize, bndloc(j), REAL(j)) 
       x_b = Gntrp_2d_ang (aloct, isize, jsize, bndloc(j), REAL(j))
         x_b = MODULO (x_b*RTH+12.0, 24.0)
         IF (j == j_midnight) x_b = 24.0
       CALL Pgdraw (x_b, y_b)
     END DO
     CALL Pgunsa ()
!
!    Put labels on plot:
!
     WRITE (info_string,'(I5.5)') label%intg(1)
     CALL pgmtxt ('T', 1.0, 0.0, 0.0, ' T='//time_string&
                  &//' REC='//TRIM(info_string) )
     CALL pgmtxt ('T', 2.5, 0.0, 0.0, label%char)
     WRITE (info_string,'(SP,3(A,F6.1))') 'STEP=',(v_levels(2)-v_levels(1))/1000.,&
            '  V_max=',MAXVAL(v)/1000., '  V_min=', MINVAL(v)/1000.
     CALL pgmtxt ('R', 2.0, 0.5, 0.5, info_string)
!
     CALL Pgunsa ()
!
!
   ELSE IF (Do_page == 5) THEN
!
!    Same as above, but in a smaller subauroral region:
!    Plot electrostatic potential contours in the ionosphere
!    in the rectangular coordinates.
!
     CALL Setup_contour_levels ()
     WRITE (6,*) 'FOR THIS PLOT, NEED 2 COLAT LIMITS'
     WRITE (6,'(A)',ADVANCE='NO') 'COLAT_LOW (DEG): '
     READ (5,*) colat_low
     WRITE (6,'(A)',ADVANCE='NO') 'COLAT_HIGH (DEG):'
     READ (5,*) colat_high
!
     CALL Pgpage ()
     CALL Pgsch (0.7)
     CALL Pgvsiz (2.0, 6.0, 2.0, 5.0)
     CALL Pgswin (0.0, 24.0, colat_high, colat_low)
     CALL Pgbox ('ABCNTS', 3.0, 3, 'ABCNTSV', 10.0, 5)
     CALL Pgmtxt ('B', 2.5, 0.5, 0.5, 'MLT, hrs')
     CALL Pgmtxt ('L', 4.0, 0.5, 0.5, '\gh, deg')
     DO i = 1, n_levels
        CALL Pgconx (v, isize, jsize, 1, isize, 1 , j_midnight-1, &
                     v_levels (i), 1, plot_v_io_cart)
     END DO
     DO i = 1, n_levels
        CALL Pgconx (v, isize, jsize, 1, isize, j_midnight , jsize-1, &
                     v_levels (i), -1, plot_v_io_cart)
     END DO
!
!
!    Put labels on plot:
!
     WRITE (info_string,'(I5.5)') label%intg(1)
     CALL pgmtxt ('T', 1.0, 0.0, 0.0, 'RCM: T='//time_string&
                  &//' REC='//TRIM(info_string) )
     CALL pgmtxt ('T', 2.5, 0.0, 0.0, label%char)
     WRITE (info_string,'(SP,3(A,F6.1))') 'STEP=',(v_levels(2)-v_levels(1))/1000.,&
            '  V_max=',MAXVAL(v)/1000., '  V_min=', MINVAL(v)/1000.
     CALL pgmtxt ('R', +2.0, 0.5, 0.5, info_string)
!
!
!
   ELSE IF (Do_page == 6) THEN
!
!  Page_4: ionospheric electric field:
!
      CALL Pgsave ()
      CALL Pgsch (0.8)
      CALL Pgsclp (0)
      efield_loop_1: DO
!
         WRITE (6,'(A,3(TR2,F6.1),A)') 'lats:', &
               lat_cut(1:3), "TYPE 3 NEW LATS, C-NO CHANGE, Q-QUIT: "
         READ (5,'(A)') lats_char
!
         IF (lats_char(1:1) == 'q'.OR. lats_char(1:1) == 'Q') THEN
            EXIT efield_loop_1
         ELSE IF (lats_char(1:1) /= 'C'.AND. lats_char(1:1) /= 'c') THEN
            READ (lats_char,*) lat_cut (1:ncut_size)
         END IF
         DO ncut = 1, ncut_size
            colat_cut (ncut) = (90.-lat_cut (ncut)) * DTR
            WRITE (char_cut(ncut), '(A5,I2,A1)')&
                           ' \gL=',NINT(lat_cut(ncut)),ACHAR(176)
            i_cut (ncut) = Get_i_for_colat (colat_cut(ncut), colat)
         END DO
!
         CALL Pgpage ()
!
!        Do left side stack of panels:
!
         DO ncut = 1, 3
            npanel = ncut
            max_y = get_max_y (lat_cut (ncut), 'EW', &
                            e_east_rlb (i_cut(ncut),:), jsize)
            x_mark_text = 18.
            y_mark_text = max_y - 0.2*2*max_y/(y2plot(npanel)-y1plot(npanel))
            CALL pgvsiz (x1plot(npanel), x2plot(npanel), y1plot(npanel), y2plot(npanel))
            CALL pgswin (0.0, 24.0, -max_y, + max_y)
    !
            CALL Pgmove (0.0,  0.0)    ! draw horizontal line of zero level
            CALL Pgdraw (24.0, 0.0)    ! but without axis tics
            IF (npanel == 1) THEN
               CALL Pgbox ('B NTS', 3.0, 3, 'ABCNTSV', 0.0, 0)
               CALL Pglab ('MLT', ' ', ' ')
            ELSE IF (npanel == 2) THEN
               CALL Pgbox ('   TS',3.0, 3, 'ABCNTSV', 0.0, 0)
            ELSE
               CALL Pgbox (' C TS', 3.0, 3, 'ABCNTSV', 0.0, 0)
            END IF

            CALL pgtext (x_mark_text, y_mark_text, char_cut (ncut))
            IF (npanel == 3) THEN
               CALL pgmtxt ('T', 1.0, 0.5, 0.5, 'E\deast\u, [mV/m]')
               WRITE (info_string,'(I5.5)') label%intg(1)
               CALL pgmtxt ('T', 2.5, 0.0, 0.0, 'RCM: T='//time_string&
                            &//'REC='//TRIM(info_string) )
               CALL pgmtxt ('T', 4.0, 0.0, 0.0, label%char)
            END IF
!
            CALL Pgmove (aloct_rlb (1,1), e_east_rlb (i_cut(ncut),1))
            DO j = 2 , wsize
               CALL Pgdraw (aloct_rlb (1,j), e_east_rlb (i_cut(ncut),j))
               IF (REAL(i_cut(ncut)) < bndloc (j)) &
                  CALL Pgpt1 (aloct_rlb (1,j), &
                              e_east_rlb (i_cut(ncut),j), 11)
            END DO
    !
         END DO
!
!        Do right side stack of panels:
!
         DO ncut = 1, ncut_size
            npanel = ncut+3
            max_y = Get_max_y (lat_cut (ncut), 'NS', &
                               e_south_rlb (i_cut(ncut),:), jsize)
            x_mark_text = 18.
            y_mark_text = max_y - 0.2*2*max_y/(y2plot(npanel)-y1plot(npanel))
            CALL Pgvsiz (x1plot(npanel), x2plot(npanel), y1plot(npanel), y2plot(npanel))
            CALL Pgswin (0., 24., -max_y, + max_y)
    !
            CALL Pgmove (0.0,  0.0)    ! draw horizontal line of zero level
            CALL Pgdraw (24.0, 0.0)    ! but without axis tics
            IF (npanel == 4) THEN
               CALL Pgbox ('B NTS',3.,3,'ABCMTSV',0.0,0)
               CALL Pglab ('MLT',' ',' ')
            ELSE IF (npanel == 5) THEN
               CALL Pgbox ('   TS',3.,3,'ABCMTSV',0.0,0)
            ELSE
               CALL Pgbox (' C TS',3.,3,'ABCMTSV',0.0,0)
               CALL Pgmtxt ('T', 1.2, 0.5, 0.5, 'E\dsouth\u, [mV/m]')
            END IF
            CALL Pgtext (x_mark_text, y_mark_text, char_cut (ncut))
!
            CALL Pgmove (aloct_rlb (1,1), e_south_rlb (i_cut(ncut),1))
            DO j = 2 , wsize
               CALL Pgdraw (aloct_rlb (1,j), e_south_rlb (i_cut(ncut),j))
               IF (REAL(i_cut(ncut)) < bndloc (j)) &
                  CALL Pgpt1 (aloct_rlb (1,j), &
                              e_south_rlb (i_cut(ncut),j), 11)
            END DO
    !
         END DO
!
      END DO efield_loop_1
      CALL Pgunsa ()
!
!
   ELSE IF (Do_page == 7) THEN
!     Page_one more: latituninal cuts:
!
!
      CALL Pgsave ()
      CALL Pgsch (0.8)
      CALL Pgsclp (0)
!
   efield_loop_2: DO
!
      WRITE (6,'(A,3(2X,F6.1),A)',ADVANCE='NO') '3 MLTs to plot are:', &
             mlt_cut (1:ncut_size),"  TYPE 3 NEW Times, 'C-NO CHANGE, Q-QUIT: " 
      READ (5,'(A)') lats_char
!
      IF (lats_char(1:1) == 'q'.OR. &
          lats_char(1:1) == 'Q') EXIT efield_loop_2
      IF (lats_char(1:1) /= 'C'.AND. lats_char(1:1) /= 'c') THEN
         READ (lats_char,*) mlt_cut (1:ncut_size)
      END IF
      DO ncut = 1, ncut_size
         WRITE (char_cut (ncut), '(A4,F4.1)')'MLT=',mlt_cut(ncut)
         aloct_cut(ncut) = MODULO (mlt_cut(ncut)/12.*pi + pi, 2.*pi)
         j_cut (ncut) = Get_j_for_loct (aloct_cut(ncut), aloct, isize, jsize )
      END DO
!
      CALL pgpage ()
      min_x = 15.0
      max_x = 35.0
      im = INT (MINVAL(bndloc))
!
!     Do left side stack of panels (east-west component):
!
      DO ncut = 1, ncut_size
         npanel = ncut
         max_y = MAXVAL ( ABS (e_east (im:isize, j_cut (ncut)) ) )
         x_mark_text = 28.0
         y_mark_text = max_y - 0.2*2*max_y/(y2plot(npanel)-y1plot(npanel))
         CALL pgvsiz (x1plot(npanel), x2plot(npanel), y1plot(npanel), y2plot(npanel))
         CALL pgswin (min_x, max_x, -max_y, + max_y)
!
         IF (npanel == 1) THEN
            CALL pgbox ('AB NTS', 0.0, 0, 'ABCNTSV', 0.0, 0)
            CALL pglab ('Colat',' ',' ')
         ELSE IF (npanel == 2) THEN
            CALL pgbox ('A   TS', 0.0, 0, 'ABCNTSV', 0.0, 0)
         ELSE
            CALL pgbox ('AC  TS', 0.0, 0, 'ABCNTSV', 0.0, 0)
             CALL pgmtxt ('T', 1.0, 0.5, 0.5, 'E\deast\u, [mV/m]')
             WRITE (info_string,'(I5.5)') label%intg(1)
             CALL pgmtxt ('T', 2.5, 0.0, 0.0, 'RCM: T='//time_string&
                          &//'REC='//TRIM(info_string) )
             CALL pgmtxt ('T', 4.0, 0.0, 0.0, label%char)
         END IF
         CALL pgtext (x_mark_text, y_mark_text, char_cut(ncut))
!
         CALL pgmove (colat (im, j_cut(ncut))/pi*180., &
                      e_east (im,j_cut(ncut)))
         DO i = im+1, isize
            CALL pgdraw (colat (i, j_cut(ncut))/pi*180., &
                         e_east (i, j_cut(ncut)))
            IF (REAL(i) < bndloc (j_cut (ncut))) &
                CALL pgpt1 (colat(i,j_cut (ncut))/pi*180.,&
                            e_east(i,j_cut (ncut)), 11)
         END DO
      END DO
!
!     Do right side stack of panels (north-south component):
!
      DO ncut = 1, ncut_size
         npanel = ncut+3
         max_y = MAXVAL ( ABS (e_south (im:isize, j_cut (ncut)) ) )
         x_mark_text = 28.0
         y_mark_text = max_y - 0.2*2*max_y/(y2plot(npanel)-y1plot(npanel))
         CALL pgvsiz (x1plot(npanel), x2plot(npanel), y1plot(npanel), y2plot(npanel))
         CALL pgswin (min_x, max_x, -max_y, + max_y)
!
         IF (npanel == 4) THEN
            CALL pgbox ('AB NTS', 0.0, 0, 'ABCNTSV', 0.0, 0)
            CALL pglab ('Colat',' ',' ')
         ELSE IF (npanel == 5) THEN
            CALL pgbox ('A   TS', 0.0, 0, 'ABCNTSV', 0.0, 0)
         ELSE 
            CALL pgbox ('AC  TS', 0.0, 0, 'ABCNTSV', 0.0, 0)
            CALL pgmtxt ('T', 1.2, 0.5, 0.5, 'E\dsouth\u, [mV/m]')
         END IF
         CALL pgtext (x_mark_text, y_mark_text, char_cut(ncut))
!
         CALL pgmove (colat (im, j_cut(ncut))/pi*180., &
                      e_south (im,j_cut(ncut)))
         DO i = im+1, isize
            CALL pgdraw (colat (i, j_cut(ncut))/pi*180., &
                         e_south (i, j_cut (ncut)))
            IF (REAL(i) < bndloc (j_cut (ncut))) &
                CALL pgpt1 (colat(i,j_cut (ncut))/pi*180.,&
                            e_south (i,j_cut (ncut)), 11)
         END DO
      END DO
!
   END DO efield_loop_2
   CALL Pgunsa ()
!
   END IF
!
!
!
   GO TO 102
!
!
102   CALL Pgend ()
      STOP
      CONTAINS
!
!
         SUBROUTINE Setup_contour_levels ()
         IMPLICIT NONE
!
         INTEGER (iprec) :: i
         REAL (rprec) :: v_low, v_high, v_step
!
!        Set up contour levels for plotting. Plot contours from
!        -50 to +50 kV at intervals of DELTA_C_SPACING, given
!        from command line. Should be no more than 100 contours,
!        so minimum spacing at present is 1 kV.
!
         v_low =  -100.0 * 1.0E+3
         v_high = +100.0 * 1.0E+3
         v_step = +2.0   * 1.0E+3
         IF (v_low >= v_high .AND. (v_high-v_low)/v_step <= 1) THEN
            WRITE (*,'(T2,A)') ' SETUP OF CONTOUR LEVELS IS WRONG'
            STOP
         END IF
!
         IF (ALLOCATED (v_levels)) DEALLOCATE (v_levels)
         n_levels = CEILING ((v_high-v_low)/v_step)
         ALLOCATE (v_levels (n_levels))
!
         v_levels (1) = v_low
         DO i = 2, n_levels
            v_levels (i) = v_levels (i-1) + v_step
         END DO
         RETURN
         END SUBROUTINE Setup_contour_levels
!
!
      FUNCTION Get_max_y (latitude, component, array_1d, nmax)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN)   :: nmax
      REAL (rprec), INTENT (IN)      :: latitude, array_1d(nmax)
      CHARACTER (LEN=*), INTENT (IN) :: component
      REAL                           :: Get_max_y
!
!     Get y-limits (plus/minus get_max_y) for plotting
!     electric field in the ionosphere.
!
      IF (component == 'EW') THEN
!
         IF (latitude <= 30.) THEN
            get_max_y = 1.0
         ELSE IF (latitude <= 35.) THEN
            get_max_y = 1.5
         ELSE IF (latitude <= 45.) THEN
            get_max_y = 2.0
         ELSE IF (latitude <= 60.0 .AND. latitude > 50.0) THEN
            get_max_y = 5.0
         ELSE 
            get_max_y = 15.0
         END IF
!
      ELSE IF (component == 'NS' ) THEN
!
         IF (latitude <= 30.) THEN
            get_max_y = 1.0
         ELSE IF (latitude <= 35.) THEN
            get_max_y = 1.5
         ELSE IF (latitude <= 45.) THEN
            get_max_y = 2.0
         ELSE IF (latitude <= 60.0 .AND. latitude>50.0) THEN
            get_max_y = 10.0
         ELSE 
            get_max_y = 50.0
         END IF
!
      ELSE
         WRITE (*,*) array_1d(:)
         STOP 'component in GET_MAX_Y is not one of options'
      END IF
      RETURN
      END FUNCTION Get_max_y
!
!
      FUNCTION Get_i_for_colat (colat_1, colat)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN)    :: colat_1, colat (:,:)
      INTEGER (iprec)              :: Get_i_for_colat
!
!     Find the I-value such as to minimize
!     ABS(colat(i,:)-colat) (that is, such that colat(i,:) is
!     the closest to a given value COLATITUDE)
!
      INTEGER (iprec) :: i, j, i_max
!
      i_max = SIZE (colat, DIM = 1)
      j = 1
!
      DO i = 1, i_max - 1
        IF (colat (i,j) <= colat_1 .AND. &
            colat (i+1,j) > colat_1) THEN
           IF ( ABS (colat(i,j)-colat_1) < ABS (colat(i+1,j)-colat_1)) THEN
                Get_i_for_colat = i
           ELSE
                Get_i_for_colat = i + 1
           END IF
           RETURN
        END IF
      END DO 
!
      IF (colat_1 >= colat (i_max, j)) THEN
         WRITE (*,*) 'colatitude > colat_max !!! PRESS ENTER TO CONTINUE...'
         READ (*,*)
         Get_i_for_colat = i_max
         RETURN
      ELSE
         STOP 'could not find I in GET_I_FOR_COLAT'
      END IF
!
      RETURN
      END FUNCTION Get_i_for_colat
!
!
      FUNCTION Get_j_for_loct (loctime, aloct, i_max, j_max)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: i_max, j_max
      REAL (rprec), INTENT (IN) :: loctime, aloct (i_max,j_max)
      INTEGER (iprec) :: Get_j_for_loct
!
!     Find the J-value such as to minimize
!     ABS(aloct(:,j)-loctime) (that is, such that aloct(:,J) is
!     the closest to a given value LOCTIME). In this function,
!     ALOCT is meant to be relabeled for MLT. (j=1 is LT=0,
!     j=j_max=wsize  is LT=24).
!
      INTEGER (iprec) :: i, j
!
      i = 1
!
      DO j = 1, j_max - 1
        IF (aloct(i,j) <= loctime .AND. &
            aloct (i,j+1) > loctime) THEN
           IF ( ABS (aloct(i,j)-loctime) < &
                ABS (aloct(i,j+1)-loctime)) THEN
                get_j_for_loct = j
           ELSE
                get_j_for_loct = j + 1
           END IF
           RETURN
        END IF
      END DO 
!
      IF (loctime >= aloct (i, j_max)) THEN
         WRITE (6,*) 'loctime > loctime_max !!!'
         get_j_for_loct = j_max
         RETURN
      ELSE
         STOP 'could not find J in GET_J_FOR_LOCT'
      END IF
!
      RETURN
      END FUNCTION Get_j_for_loct
!
!
      SUBROUTINE Compute_efield (e_south, e_east, Iswitch)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: Iswitch
      REAL (rprec), INTENT(IN OUT):: e_south (isize,jsize), &
                                     e_east (isize,jsize)
!
!------------------------------------------------------------
!     Subroutine to differentiate RCM electric potential on 
!     the ionospheric grid to compute the two components of
!     the electric field. "e_east" is the easward electric
!     field (horizontal, positive east); "e_south" is the
!     meridional component of electric field perpendicular to
!     B-field (positive south). Since the RCM potential is in
!     [Volts], and distances are in [km], electric fields are 
!     computed in [mV/m]. Reference altitude for E-field is
!     R_i (ionosphere).
!
!     In computing E-field, we assume that V wraps around as
!     in RCM, that is, first phi-point (noon) is j=jwrap, 
!     last point (just before noon after going around once) is
!     j = jmax - 1, point j=jmax wraps onto j=jwrap, etc.
!
!     Iswitch = 1 means do central differences
!             = 2 means two-point (forward) differences
!
!------------------------------------------------------------
!
      INTEGER (iprec) :: i,j
      REAL (rprec) :: d_phi, d_theta, d_V
!
!
!
!  1. Compute fi-component (eastward) of electric field:
!     e_east = - (1./Ri/sin(theta)) * dV / dfi
!
   IF (Iswitch == 1) THEN
!
      DO j = 1, jsize
      DO i = 1,isize
!
         IF (aloct(i,j+1) < aloct(i,j-1)) THEN
            d_phi = aloct(i,j+1) + 2.*pi - aloct (i,j-1)
         ELSE
            d_phi = aloct (i,j+1) - aloct (i,j-1)
         END IF
         IF (d_phi > pi/8. .OR. d_phi <= 0.) &
            STOP 'problem in COMPUTE_EFIELD'
!
         d_v = v(i,j+1) - v(i,j-1)
         e_east(i,j) = - d_v / d_phi / Ri / SIN(colat(i,j))
!
      END DO
      END DO
!
   ELSE IF (Iswitch == 2) THEN
!
      DO j = 1, jsize
      DO i = 1,isize
!
         IF (aloct(i,j+1) < aloct(i,j)) THEN
            d_phi = aloct(i,j+1) + 2.*pi - aloct (i,j)
         ELSE
            d_phi = aloct (i,j+1) - aloct (i,j)
         END IF
         IF (d_phi > pi/8. .OR. d_phi <= 0.) &
            STOP 'problem in COMPUTE_EFIELD'
!
         d_v = v(i,j+1) - v(i,j)
         e_east(i,j) = - d_v / d_phi / Ri / SIN(colat(i,j))
!
      END DO
      END DO
!
   ELSE
!
      STOP 'Iswitch not implemented in compute_efield'
!
   END IF
!
      CALL Wrap_around_ghostcells (e_east, isize, jsize, n_gc)
!
!
!  2.1 Compute theta-component of electric field:
!      e_south = - (1./Ri) * dV / dtheta
!
!
   IF (Iswitch == 1) THEN
!
      DO j = 1, jsize
         DO i = 1, isize 
!
            IF (i == 1) THEN
               d_v     = v    (i+1,j) - v    (i,j)
               d_theta = colat(i+1,j) - colat(i,j)
            ELSE IF ( i == isize) THEN
               d_v     = v    (i,j) - v    (i-1,j)
               d_theta = colat(i,j) - colat(i-1,j)
            ELSE
               d_v     = v    (i+1,j) - v    (i-1,j)
               d_theta = colat(i+1,j) - colat(i-1,j)
            END IF
!
            e_south (i,j) = - d_v / d_theta / Ri / sini(i,j)
!
         END DO
      END DO
!
   ELSE IF (Iswitch == 2) THEN
!
      DO j = 1, jsize
         DO i = 1, isize 
!
            IF (i == 1) THEN
               d_v     = v    (i+1,j) - v    (i,j)
               d_theta = colat(i+1,j) - colat(i,j)
            ELSE IF ( i == isize) THEN
               d_v     = v    (i,j) - v    (i-1,j)
               d_theta = colat(i,j) - colat(i-1,j)
            ELSE
               d_v     = v    (i+1,j) - v    (i,j)
               d_theta = colat(i+1,j) - colat(i,j)
            END IF
!
            e_south (i,j) = - d_v / d_theta / Ri / sini(i,j)
!
         END DO
      END DO
!
   ELSE 
! 
      STOP 'Iswitch not implemented in compute_efield'
!
   END IF
!
   CALL Wrap_around_ghostcells (e_south, isize, jsize, n_gc)
!
      RETURN
      END SUBROUTINE Compute_efield
!

!
!
      END PROGRAM Plot_v
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
              SUBROUTINE Plot_v_io (visble, x, y, z)
              USE Rcm_variables, ONLY : iprec, rprec, pi, colat, aloct,&
                    isize, jsize, RTD, bndloc
              IMPLICIT NONE
              REAL (rprec), INTENT (IN) :: x, y, z
              INTEGER (iprec), INTENT (IN) :: visble

              REAL (rprec) :: xworld, yworld, r_world, p_world, Gntrp_2d, Gntrp_2d_ang, Bndy
              LOGICAL :: erase
!
!             X is colatitude (radius)
!             Y is local time angle (phi)
!
              erase = .FALSE.
              IF (x < Bndy (bndloc, jsize, y) - 0.001) erase = .TRUE.
              r_world = RTD * Gntrp_2d (colat, isize, jsize, x, y)
              p_world = Gntrp_2d (aloct, isize, jsize, x, y)
!
              xworld = r_world *COS(p_world+pi)
              yworld = r_world *SIN(p_world+pi)
!
!             IF (z < 0.0_rprec) THEN
!                CALL Pgsls (4)
!             END IF
              IF (visble == 0 .OR. erase) THEN
                 CALL Pgmove (xworld, yworld)
              ELSE
                 CALL Pgdraw (xworld, yworld)
              END IF
              CALL Pgsls (1)
              RETURN
              END SUBROUTINE Plot_v_io
!
!
!-------------------------------------------
!
!
              SUBROUTINE Plot_v_equat (visble, x, y, z)
              USE Rcm_variables, ONLY : iprec, rprec, isize,jsize,xmin, ymin, bndloc
              IMPLICIT NONE
              REAL (rprec), INTENT (IN) :: x,y,z
              INTEGER (iprec), INTENT (IN) :: visble

              REAL (rprec) :: xworld, yworld, Gntrp_2d, Gntrp_2d_ang, Bndy
              LOGICAL ::erase
!
!             X is GSM x in Re
!             Y is GSM y in Re
!
!             Determine if we are outside the boundary:
!
              erase = .FALSE.
              IF (x < Bndy (bndloc, jsize, y) + 0.001) erase = .TRUE.
!
              xworld = Gntrp_2d (xmin, isize, jsize, x, y)
              yworld = Gntrp_2d (ymin, isize, jsize, x, y)
!
              IF (z <= 0.0_rprec) THEN
!                CALL Pgsls (4)
                 CALL Pgsci (4)
              ELSE
                 CALL Pgsci (2)
              END IF
              IF (visble == 0 .OR. erase) THEN
                 CALL Pgmove (xworld, yworld)
              ELSE
                 CALL Pgdraw (xworld, yworld)
              END IF
!             CALL Pgsls (1)
              CALL Pgsci (1)
              RETURN
              END SUBROUTINE Plot_v_equat
!
!
!
              SUBROUTINE Plot_v_io_cart (visble, x, y, z)
              USE Rcm_variables, ONLY : iprec, rprec, colat, aloct, &
                  isize, jsize, RTD, RTH, bndloc
!
!             For making contour plots of RCM V array in the
!             ionosphere in cartesian coordinates; x is MLT in hrs,
!              y is colat in deg.
!
              IMPLICIT NONE
              REAL (rprec), INTENT (IN) :: x, y, z
              INTEGER (iprec), INTENT (IN) :: visble
              LOGICAL :: erase
!
              REAL (rprec) :: xworld, yworld, Gntrp_2d, Gntrp_2d_ang, Bndy
!
!
!             X is GSM x in Re
!             Y is GSM y in Re
!
!             Determine if we are outside the boundary:
!
              erase = .FALSE.
              IF (x < Bndy (bndloc, jsize, y) - 0.001) erase = .TRUE.
!
              xworld = Gntrp_2d_ang (aloct, isize, jsize, x, y)
              yworld = Gntrp_2d (colat, isize, jsize, x, y)
              xworld = MODULO (xworld*RTH + 12.0, 24.0)

              yworld = yworld*RTD
!
              IF (z < 0.0_rprec) THEN
                 CALL Pgsls (4)
              END IF
              IF (visble == 0 .OR. erase) THEN
                 CALL Pgmove (xworld, yworld)
              ELSE
                 CALL Pgdraw (xworld, yworld)
              END IF
              CALL Pgsls (1)
              RETURN
              END SUBROUTINE Plot_v_io_cart

      SUBROUTINE Palett (type, contra, bright)
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
      INTEGER, INTENT (IN) :: type
      REAL, INTENT (IN) :: contra, bright
!
      REAL :: gl(2), gr(2), gg(2), gb(2)
      REAL :: rl(9), rr(9), rg(9), rb(9)
      REAL :: hl(5), hr(5), hg(5), hb(5)
      REAL :: wl(10), wr(10), wg(10), wb(10)
      REAL :: al(20), ar(20), ag(20), ab(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
               0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
               0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (type == 1) THEN !        -- gray scale
         CALL Pgctab (gl, gr, gg, gb, 2, contra, bright)
      ELSE IF (type == 2) THEN !        -- rainbow
         CALL Pgctab (rl, rr, rg, rb, 9, contra, bright)
      ELSE IF (type == 3) THEN !        -- heat
         CALL Pgctab (hl, hr, hg, hb, 5, contra, bright)
      ELSE IF (type == 4) THEN !        -- weird IRAF
         CALL Pgctab (wl, wr, wg, wb, 10, contra, bright)
      ELSE IF (type == 5) THEN !        -- AIPS
         CALL Pgctab (al, ar, ag, ab, 20, contra, bright)
      END IF
      END SUBROUTINE Palett
