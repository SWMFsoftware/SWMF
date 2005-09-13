      PROGRAM Plot_eeta
      USE Rcm_variables
      USE Rcm_io
      IMPLICIT NONE
      CHARACTER (LEN=1) :: prefix,answer
      CHARACTER (LEN=10) :: time_string
      CHARACTER (LEN=62) :: info_string
      CHARACTER (LEN=05) :: itrack_string
      CHARACTER (LEN=80) :: dev_name, dev_type, plot_label, file_printout_name
      CHARACTER (LEN=6) :: ut_string
      CHARACTER (LEN=11), PARAMETER :: date_string_all(4) = (/ &
!        'Oct-18-1998', 'Oct-19-1998', 'Oct-20-1998', 'Oct-21-1998' /)
         '', '', '', '' /)
      CHARACTER (LEN=11) :: date_string
!     
      REAL :: x_bnd(jsize), y_bnd(jsize), &
              eeta_to_plot(1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
              xgsm_max, xgsm_min, ygsm_max, ygsm_min, xgsm_step, ygsm_step, xstep_tmp, ystep_tmp,&
              bright, contra, a1, a2, av, xpoly(4), ypoly(4), x_wedge(2), y_wedge(2), &
              wedge(1000), tr(6), x_plot(2), y_plot(2), x_win(2), y_win(2), &
              eeta_max_array (kcsize), x_b, y_b, x_b_n, y_b_n, &
              x1_new_win, x2_new_win, y1_new_win, y2_new_win, gntrp_2d, Bndy
!
      INTEGER :: array(ksize),k_count,kc_beg,kc_end , c1, c2, kc, &
                 type_of_plot, type_of_frame
!
      INTEGER :: n, k, i,  pgopen, getcwd, npts, j, n_rec, pgid, minind, maxind, iv
      INTEGER :: mpts (1000)
      CHARACTER (LEN=2) :: label_pts_number(10000),l_pts_number
      CHARACTER (LEN=5) :: label_pts_itrack(10000),l_pts_itrack
      REAL :: xpts(1000), ypts(1000)
!
!                                   different edges.
      CHARACTER (LEN=40), DIMENSION(80) :: label_pts
      CHARACTER (LEN=80) :: c_format
!
      LOGICAL :: Do_plot_gridlines = .FALSE., p_in, p_in_n
      INTEGER :: n_contours, color_use
      REAL :: eeta_min, eeta_max
      REAL :: cont_level, cont_step
      EXTERNAL  plot_equat, Plot_ij
!
      INTEGER (iprec) :: nx_rg, ny_rg
      REAL (rprec), ALLOCATABLE :: x_rg (:), y_rg (:), z_rg (:,:)
!
!
      CALL Getarg (1, info_string)
      READ (info_string,'(I6)') n_rec
      CALL Getarg (2, info_string)
      READ (info_string, '(I10)') type_of_plot
      CALL Getarg (3, info_string)
      READ (info_string, '(L)') do_plot_gridlines
      CALL Getarg (4, info_string)
      READ (info_string, '(A80)') dev_type
      CALL Getarg (5, info_string) 
      READ (info_string, '(I2)') color_use
      IF (color_use /= 0 .AND. color_use /= 1 .AND. color_use /= 2) STOP 'ILLEGAL COLOR_USE'
      CALL Getarg (6, info_string)
      READ (info_string, '(A80)') plot_label
      plot_label = TRIM(plot_label)
      CALL Getarg (7, info_string)
      READ (info_string, '(I)') type_of_frame
      IF (type_of_frame < 1 .OR. type_of_frame > 3) STOP 'TYPE_OF_FRAME'
!
      CALL Read_grid ()
! 
!
         CALL Read_array ('rcmbndloc', n_rec, label, ARRAY_1D = bndloc)
         imin_j = CEILING(bndloc)
!
         CALL Read_array ('rcmxmin',   n_rec, label, ARRAY_2D = xmin)
         CALL Read_array ('rcmymin',   n_rec, label, ARRAY_2D = ymin)
      IF (type_of_plot == 1) THEN
         CALL Read_array ('rcmbirk',   n_rec, label, ARRAY_2D = birk)
      ELSE
         stop
         CALL Read_array ('rcmbirkavg',   n_rec, label, ARRAY_2D = birk)
      END IF
!do j = 1 , jsize
!  xmin(1:imin_j(j)-1,j) = xmin(imin_j(j),j)
!  ymin(1:imin_j(j)-1,j) = ymin(imin_j(j),j)
!end do
!
         date_string = date_string_all (1+FLOOR(label%intg(6)/3600.0/24.00))
!
         time_string = Get_time_char_string (label)
         WRITE (time_string,'(I2.2,A1,I2.2,A1,I2.2)') &
               label%intg(3),':',label%intg(4),':',label%intg(5)
         WRITE (6,'(A,A10,A,I4)') 'TIME IS: ',time_string
         ut_string = time_string(1:2)//time_string(4:5)//time_string(7:8)
         WRITE (time_string(1:2),'(I2.2)') MODULO(label%intg(3),24)
!
      IF (dev_type(1:LEN_TRIM(dev_type)) == 'xwin' .OR. &
          dev_type(1:LEN_TRIM(dev_type)) == 'XWIN' ) THEN
         dev_name = '/'//'xwin'
      ELSE IF (dev_type(1:LEN_TRIM(dev_type)) == 'GIF' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'gif' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'PS'  .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'ps' ) THEN
         dev_name = 'jbirk_'//ut_string//'.'//&
                    dev_type(1:LEN_TRIM(dev_type))//'/'//&
                    dev_type(1:LEN_TRIM(dev_type))
      ELSE IF (dev_type(1:LEN_TRIM(dev_type)) == 'VPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'vps' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'CPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'cps' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'VCPS' .OR. &
               dev_type(1:LEN_TRIM(dev_type)) == 'vcps' ) THEN
         dev_name = 'jbirk_'//ut_string//'.'// &
                    'ps'//'/'//&
                    dev_type(1:LEN_TRIM(dev_type))
      ELSE
         STOP 'DEVICE-FILENAME NOT IMPLEMENTED'
      END IF
!
!
      pgid = Pgopen (dev_name)
      IF (pgid < 1) STOP 'Could not open graphics device'
      CALL Pgslw (2)
      CALL Pgscf (2)  ! roman font
!
!
      IF (dev_type == 'GIF' .OR. dev_type == 'gif' .OR. &
          dev_type == 'XWIN' .OR. dev_type == 'xwin' ) THEN
         CALL Pgpap (6.5, 1.0) !set view surface to square 6.0 inches
      END IF
      CALL Pgpage () 
      IF (dev_type == 'GIF' .OR. dev_type == 'gif'.OR.dev_type=='xwin'.OR.dev_type=='XWIN') THEN
         x_plot = (/ 0.8, 5.3 /)
         y_plot = (/ 0.8, 5.3 /)
         x_wedge = (/ 5.4, 5.7 /)
         y_wedge = (/ 0.8, 5.3 /)
      ELSE
         x_plot = (/ 1.5, 6.0 /)
         y_plot = (/ 1.5, 6.0 /)
         x_wedge = (/ 6.1, 6.4 /)
         y_wedge = (/ 1.5, 6.0 /)
      END IF
      CALL Pgvsiz (x_plot(1), x_plot(2), y_plot(1), y_plot(2))
!
!
    IF (type_of_frame == 1) THEN
!     Here coordinate system is GSM, so X runs from left to
!     right, and Y runs from top to bottom:
      x_win = (/ +15.0, -25.0 /)
      y_win = (/ +20.0, -20.0 /)
      CALL Pgwnad (x_win(1), x_win(2), y_win(1), y_win(2)) !window in eq. plane GSM
    ELSE IF (type_of_frame == 2) THEN
      x_win = (/ 20., 80. /)
      y_win = (/ REAL(1), REAL(52) /)
      CALL Pgswin (x_win(1), x_win(2), y_win(1), y_win(2)) !window in eq. plane GSM
    ELSE IF (type_of_frame == 3) THEN
      x_win = (/ 40., 80. /)
      y_win = (/ 53., REAL(jsize) /)
      CALL Pgswin (x_win(1), x_win(2), y_win(1), y_win(2)) !window in eq. plane GSM
    ELSE
      STOP 'TYPE_OF_FRAME'
    END IF
      CALL Pgqvp (1, x1_new_win, x2_new_win, y1_new_win, y2_new_win)
!
      CALL Pgsave ()
      CALL Pgstbg (15)
   IF (type_of_plot == 1) THEN
      WRITE (info_string,'(A)') 'J\dbirk\u'
   ELSE IF (type_of_plot == 2) THEN
      WRITE (info_string,'(A)') 'AVG J\dbirk\u'
   END IF
      CALL Pgmtxt ('T', 1.0, 1.0, 1.0, info_string)
      CALL Pgunsa ()
      CALL Pgmtxt ('T', 2.2, 0.0, 0.0, plot_label)
!
!
!        Compute and plot the location of the outer boundary
!        as thin line:
!
         CALL pgslw (5)
         j=1
      IF (type_of_frame == 1 ) THEN
         x_bnd (j) = gntrp_2d (xmin, isize, jsize, bndloc(j), REAL(j))
         y_bnd (j) = gntrp_2d (ymin, isize, jsize, bndloc(j), REAL(j))
!x_bnd(j) = xmin (imin_j(j),j)
!y_bnd(j) = ymin (imin_j(j),j)
      ELSE IF (type_of_frame == 2 .OR. type_of_frame == 3) THEN
         x_bnd(j) = REAL(imin_j(j))
         y_bnd(j) = REAL(j)
      ELSE
         STOP 'TYPE_OF_FRAME'
      END IF
         CALL Pgmove (x_bnd(j), y_bnd(j))
         DO j = 2, jsize
          IF (type_of_frame == 1) THEN
            x_bnd (j) = gntrp_2d (xmin, isize, jsize, bndloc(j), REAL(j))
            y_bnd (j) = gntrp_2d (ymin, isize, jsize, bndloc(j), REAL(j))
!x_bnd(j) = xmin (imin_j(j),j)
!y_bnd(j) = ymin (imin_j(j),j)
          ELSE IF (type_of_frame == 2 .OR. type_of_frame ==3) THEN
            x_bnd (j) = REAL(imin_j(j))
            y_bnd (j) = REAL (j)
          ELSE
            STOP 'TYPE_OF_FRAME'
          END IF
            CALL Pgdraw (x_bnd(j), y_bnd(j))
         END DO
         CALL pgslw (2)
!
! 
!>>>>>>>>>>>>>>>>>>
!
!    Make contour plot of EETA:
!
!    1. Figure out what quantity to plot:
!
        eeta_to_plot = birk
        eeta_max = +2.0
        eeta_min = -2.0
!
!
!     First, color plotting:
!     
      IF (color_use == 1) THEN
         CALL Pgsave ()
!        CALL Pgqcir (c1, c2)
!        WRITE (*,*) 'NUMBER OF COLOR INDICES:', MAX(0,c2-c1+1)
         bright = 0.5
         contra = 1.0
         CALL Palett (2, contra, bright)
         CALL Pgqcir (minind, maxind)
         a2 = eeta_max
         a1 = eeta_min
         IF (a1 >= a2) STOP 'STOP: A1 >= A2'
         IF (minind >= maxind) STOP 'NOT ENOUGH COLORS'
         DO j = 1, jsize
         DO i = 1, isize
            IF (i < imin_j(j)) CYCLE
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
          IF (type_of_frame == 1) THEN
            xpoly (1) = Gntrp_2d (xmin, isize, jsize, REAL(i-0.5,rprec), REAL(j-0.5, rprec))
            ypoly (1) = Gntrp_2d (ymin, isize, jsize, REAL(i-0.5,rprec), REAL(j-0.5, rprec))
            xpoly (2) = Gntrp_2d (xmin, isize, jsize, REAL(i-0.5,rprec), REAL(j+0.5, rprec))
            ypoly (2) = Gntrp_2d (ymin, isize, jsize, REAL(i-0.5,rprec), REAL(j+0.5, rprec))
            xpoly (3) = Gntrp_2d (xmin, isize, jsize, REAL(i+0.5,rprec), REAL(j+0.5, rprec))
            ypoly (3) = Gntrp_2d (ymin, isize, jsize, REAL(i+0.5,rprec), REAL(j+0.5, rprec))
            xpoly (4) = Gntrp_2d (xmin, isize, jsize, REAL(i+0.5,rprec), REAL(j-0.5, rprec))
            ypoly (4) = Gntrp_2d (ymin, isize, jsize, REAL(i+0.5,rprec), REAL(j-0.5, rprec))
!
          ELSE IF (type_of_frame == 2 .OR. type_of_frame == 3) THEN
            xpoly (1) = REAL(i-0.5,rprec)
            ypoly (1) = REAL(j-0.5, rprec)
            xpoly (2) = REAL(i-0.5,rprec)
            ypoly (2) = REAL(j+0.5, rprec)
            xpoly (3) = REAL(i+0.5,rprec)
            ypoly (3) = REAL(j+0.5, rprec)
            xpoly (4) = REAL(i+0.5,rprec)
            ypoly (4) = REAL(j-0.5, rprec)
          ELSE
            STOP 'TYPE_OF_FRAME'
          END IF
            CALL Pgsci (iv)
            CALL Pgpoly (4, xpoly, ypoly)
         END DO
         END DO
         CALL Pgunsa ()
!        CALL Pgvsiz (x_wedge(1), x_wedge(2), y_wedge(1), y_wedge(2))
         CALL Pgvsiz (x_wedge(1), x_wedge(2), y1_new_win, y2_new_win)
         CALL Pgswin (0.9, 1.1, 1.0, REAL(SIZE(wedge)))
         tr = (/ 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 /)
!        a1 = eeta_min / eeta_max * 100.0
!        a2 = 100.0
         a1 = eeta_min 
         a2 = eeta_max
!        DO i = 1, SIZE(wedge)
!           wedge(i) =  (SIZE(wedge)*a1-a2)/(SIZE(wedge)-1) + &
!                       (a2-a1)/(SIZE(wedge)-1)*REAL(i,rprec)
!        END DO
!        CALL Pgimag (wedge, SIZE(wedge), 1, 1, SIZE(wedge), 1, 1, a1, a2, tr)
         CALL Pgsave ()
         CALL Pgsch (1.5)
         CALL Pgwedg ('RI', -1.0, 3.0, a1, a2, 'J\dbirk\u, \gmA/m\u2\d')
         CALL Pgunsa ()
         CALL Pgvsiz (x_plot(1), x_plot(2), y_plot(1), y_plot(2))
        IF (type_of_frame == 1) THEN
         CALL Pgwnad (x_win(1), x_win(2), y_win(1), y_win(2)) !window in eq. plane GSM
        ELSE IF (type_of_frame == 2 .OR. type_of_frame == 3) THEN
         CALL Pgswin (x_win(1), x_win(2), y_win(1), y_win(2)) !window in eq. plane GSM
        ELSE
         STOP 'TYPE_OF_FRAME'
        END IF
      END IF
!
     IF (color_use == 0) THEN
      cont_step = 0.5
      cont_level = eeta_max + cont_step
      DO 
         cont_level = cont_level - cont_step
         IF (cont_level < eeta_min) EXIT
        IF (type_of_frame == 1) THEN
         CALL Pgconx (eeta_to_plot , isize, jsize, 1, isize, 1, jsize, &
                      cont_level, +1, Plot_equat)
        ELSE IF (type_of_frame == 2 .OR. type_of_frame == 3) THEN
         CALL Pgconx (eeta_to_plot , isize, jsize, 1, isize, 1, jsize, &
                      cont_level, +1, Plot_ij)
        ELSE
           STOP 'TYPE_OF_FRAME'
        END IF
      END DO
     END IF
!
!
!     Write contour step and info on this energy channel on the plot:
!
      CALL Pgsave ()
      CALL Pgscf (1)
      CALL Pgsch (0.75)
      CALL Pgstbg (0)
      WRITE (info_string,'(A,SP,F6.2)') 'MAX=', MAXVAL(birk)
      CALL Pgmtxt ('T', -3.0, 0.25, 0.0, info_string)
!
      WRITE (info_string,'(A,SP,F6.2)') 'MIN=', MINVAL(birk)
      CALL Pgmtxt ('T', -1.5, 0.25, 0.0, info_string)
!
      WRITE (info_string,'(A,F5.1)') '\gDJ\dbirk\u=',cont_step
      CALL Pgmtxt ('T', -1.5, 0.65, 0.0, info_string)
      CALL Pgunsa ()
!      
      CALL Pgbox ('abctsn', 0.0, 0, 'abctsnv', 0.0, 0)
    IF (type_of_frame == 1) THEN
      CALL Pgmtxt ('B', 3.0, 0.5, 0.5, 'X\dGSM\u')
      CALL Pgmtxt ('L', 3.0, 0.5, 0.5, 'Y\dGSM\u')
    ELSE IF (type_of_frame == 2 .OR. type_of_frame == 3) THEN
      CALL Pgmtxt ('B', 3.0, 0.5, 0.5, 'I-grid')
      CALL Pgmtxt ('L', 3.0, 0.5, 0.5, 'J-grid')
    ELSE
      STOP 'TYPE_OF_FRAME'
    END IF
         CALL Read_array ('rcmbndloc', n_rec, label, ARRAY_1D = bndloc) !get LABEL again
      WRITE (info_string,'(I5.5)') label%intg(1)
      CALL Pgmtxt ('T', 1.0, 0.0, 0.0, 'RCM: UT='//time_string//' '//date_string)
!
!
!     Plot the grid lines mapped out to the equator.plane if requested:
!
      IF (Do_plot_gridlines) THEN
        CALL Pgsave ()
        CALL Pgsls (1)
        DO i = isize, 1, -1
           j = 1
           x_b_n = xmin(i,j)
           y_b_n = ymin(i,j)
           p_in_n = .TRUE.
           IF (i<imin_j(j)) THEN
              p_in_n = .FALSE.
           END IF
           CALL Pgmove (x_b, y_b)
           DO j = 1+1, jsize
              x_b = x_b_n
              y_b = y_b_n
              p_in = p_in_n
              x_b_n = xmin(i,j)
              y_b_n = ymin(i,j)
              p_in_n = .TRUE.
              IF (i<imin_j(j)) THEN
                 p_in_n = .FALSE.
              END IF
              IF (p_in .AND. p_in_n) THEN
                 CALL Pgmove (x_b,y_b)
                 CALL Pgdraw (x_b_n,y_b_n)
              END IF
           END DO
        END DO
        CALL Pgsls (2)
        DO j = 1, jsize - 1
           CALL Pgline (isize-imin_j(j)+1, xmin(imin_j(j):,j), ymin(imin_j(j):,j))
        END DO
        CALL Pgunsa ()
      END IF
!
!
      CALL Pgend ()
      STOP
      CONTAINS
!
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
!
      END PROGRAM plot_eeta
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Plot_equat (visble, x, y, z)
      USE Rcm_variables
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: x,y,z
      INTEGER (iprec), INTENT (IN) :: visble
!
      INTEGER (iprec) :: i, j
      REAL (rprec) :: x_world,y_world, gntrp_2d, bndy
      LOGICAL :: erase
!
!
      erase = .FALSE.
      IF (x < Bndy(bndloc,y)) erase = .TRUE.
      x_world = Gntrp_2d (xmin, isize, jsize, x, y)
      y_world = Gntrp_2d (ymin, isize, jsize, x, y)
!
      IF (visble == 0 .or. erase) THEN
         CALL Pgmove (x_world, y_world)
      ELSE
         CALL Pgdraw (x_world, y_world)
      END IF
      RETURN
      END SUBROUTINE Plot_equat
!
!
      SUBROUTINE Plot_ij (visble, x, y, z)
      USE Rcm_variables
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: x,y,z
      INTEGER (iprec), INTENT (IN) :: visble
!
      INTEGER (iprec) :: i, j
      REAL (rprec) :: x_world,y_world, Bndy
      LOGICAL :: erase
!
!
      erase = .FALSE.
      IF (x < Bndy(bndloc,y)) erase = .TRUE.
      x_world = x
      y_world = y
!
      IF (visble == 0 .or. erase) THEN
         CALL Pgmove (x_world, y_world)
      ELSE
         CALL Pgdraw (x_world, y_world)
      END IF
      RETURN
      END SUBROUTINE Plot_ij
