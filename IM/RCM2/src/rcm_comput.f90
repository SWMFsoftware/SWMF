      SUBROUTINE Comput (jtime, dt )
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
      REAL (rprec),    INTENT (IN) :: dt
!
      INTEGER (iprec) :: j, i, k
      REAL (rprec)  ::  a(3), b(3), dx(3), dy(3), deqdt
!
!
      CALL Interp_by_marktimes (SIZE(ivtime), ivtime, vinput,       jtime, vdrop)
      CALL Interp_by_marktimes (SIZE(ivtime), ivtime, vinput_phase, jtime, vdrop_phase)
      CALL Interp_by_marktimes (SIZE(ikptime), ikptime, kpinput,     jtime, Kp)
      CALL Get_bfield  (jtime)
      CALL Get_boundary (boundary, bndloc)
      IF (i_eta_bc == 1) THEN
         CALL    Get_eta_on_bndy (jtime)
      ELSE IF (i_eta_bc == 2) THEN
         CALL Rcm_plasma_bc (i_eta_bc, 1)
      ELSE IF (i_eta_bc == 3) THEN
         CALL Rcm_plasma_bc (i_eta_bc, 1)
      ELSE 
         call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
             'ILLEGAL VALUE OF I_eta_bc')
      END IF
      CALL Get_jbirk (jtime)
      CALL Get_vparallel ()
!
!
      SELECT CASE (icond)   ! compute auroral conductances enhancements
      CASE (1)
         CALL Get_active_cond ( )
      CASE (2)
         CALL Get_hardy_cond ()
      CASE (11) ! for coupling, set S_h=0, S_p=4*2, no auroral enh.
         hall = 0.0
         pedpsi = 8.0
         pedlam = 8.0 / sini**2
         CALL Wrap_around_ghostcells (pedlam, isize, jsize, n_gc)
         CALL Wrap_around_ghostcells (pedpsi, isize, jsize, n_gc)
         CALL Wrap_around_ghostcells (hall, isize, jsize, n_gc)
      CASE DEFAULT
         call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
             'COMPUT: ILLEGAL icond')
      END SELECT

      RETURN 
      END SUBROUTINE Comput                         
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Interp_by_marktimes (nmax, mark_times, mark_values, jtime, value)
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: nmax, jtime, mark_times(nmax)
      REAL (rprec), INTENT (IN)  :: mark_values (nmax) 
      REAL (rprec), INTENT (OUT) :: value
!                                                                       
!-------------------------------------------------------------
!     Subroutine to interpolate an input value using "marktimes".
!     Used to get polar cap potential drop, its "phase", Kp, etc.
!     If jtime <= mark_times(1) then value = mark_values(1)
!     If jtime >  mark_times(nmax) then value = mark_values(nmax)
!     all other cases--interpolated.
!-------------------------------------------------------------
!
      INTEGER (iprec) :: n
      REAL (rprec)    :: f
!
      DO n = 1, nmax 
         IF (jtime <= mark_times (n) ) THEN 
            IF (n == 1) THEN 
               value = mark_values (1)
               RETURN 
            ELSE 
               f = REAL(jtime-mark_times(n-1),rprec) / &
                   REAL(mark_times(n)-mark_times(n-1), rprec)
               value = (1.0_rprec - f) * mark_values(n-1) + f * mark_values(n)
               RETURN 
            END IF 
         END IF 
      END DO 
      value = mark_values (nmax)
!
      RETURN 
      END SUBROUTINE Interp_by_marktimes
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_bfield (jtime )
      USE Rcm_variables
      USE Rcm_io
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
!
!_____________________________________________________________________________
!
!     Subroutine to find magnetic field arrays r,p,be and vm 
!     at time itime by interpolating in time between 
!     precomputed magnetic field
!     models.  nbf is the number of precomputed bfield models
!     and ibtime is a vector giving the event times associated
!     with each of the models.
!
!     rws   3/19/97; stanislav 5/01/98
!    
!     Stanislav: nold_bf is initialized to 0; this is actually
!     enough to "save" its value, but just in case add SAVE 
!     attribute. nold_bf is 
!     incremented by 1 only when the appropriate sets of 
!     B-models are read in from the files for interpolation;
!     "appropriate" here is:
!      | first B-model if itime <= ibtime(1)
!      | last  B-model if itime >  ibtime(nbf)
!      | two B-models for ibtime(n-1) and ibtime(n) where n
!      |  is such that ibtime(n-1) < itime <= ibtime(n)
!_____________________________________________________________________________
!
!
      INTEGER (iprec), SAVE :: nold_bf = 0
      INTEGER (iprec) :: n, nn, Lrec, nbf, i,j
      REAL    (rprec) :: f, fstoff1, fstoff2, fdst1, fdst2, fmeb1, fmeb2, &
                         fclps1, fclps2
!
!
 IF (itype_bf == 2) THEN ! use friction code results
    Lrec = 1
    CALL Read_array ('input/rcmxmin_inp', LREC , label, ARRAY_2D=xmin, ASCI=asci_flag)
    CALL Read_array ('input/rcmymin_inp', LREC , label, ARRAY_2D=ymin, ASCI=asci_flag)
    CALL Read_array ('input/rcmvm_inp',   LREC , label, ARRAY_2D=vm  , ASCI=asci_flag)
    CALL Read_array ('input/rcmbmin_inp', LREC , label, ARRAY_2D=bmin, ASCI=asci_flag)
    CALL Read_array ('input/rcmrmin_inp', LREC , label, ARRAY_2D=rmin, ASCI=asci_flag)
    CALL Read_array ('input/rcmpmin_inp', LREC , label, ARRAY_2D=pmin, ASCI=asci_flag)
    RETURN
ELSE IF (itype_bf == 1) THEN ! interpolate HV
!
  nbf = SIZE (ibtime(:))
!
  IF (jtime <= ibtime(1)) THEN
!
     IF (nold_bf < 1) THEN
        Lrec = 1
        CALL Read_array ('input/rcmxmin_inp', LREC , label, &
                          ARRAY_2D = xmin, ASCI=asci_flag)
        CALL Read_array ('input/rcmymin_inp', LREC , label, &
                          ARRAY_2D = ymin, ASCI=asci_flag)
        CALL Read_array ('input/rcmvm_inp',   LREC , label, &
                          ARRAY_2D = vm  , ASCI=asci_flag)
        CALL Read_array ('input/rcmbmin_inp', LREC , label, &
                          ARRAY_2D = bmin, ASCI=asci_flag)
!
        fmeb   = label%real (12)
        fstoff = label%real (13)
        fdst   = label%real (14)
        fclps  = label%real (15)
        nold_bf   = 1
     END IF
!
  ELSE IF (jtime > ibtime(nbf)) THEN
!
     IF (nold_bf < nbf+1) THEN
!
       Lrec = nbf
       CALL Read_array ('input/rcmxmin_inp', LREC , label, ARRAY_2D = xmin, ASCI=asci_flag)
       CALL Read_array ('input/rcmymin_inp', LREC , label, ARRAY_2D = ymin, ASCI=asci_flag)
       CALL Read_array ('input/rcmvm_inp',   LREC , label, ARRAY_2D = vm  , ASCI=asci_flag)
       CALL Read_array ('input/rcmbmin_inp', LREC , label, ARRAY_2D = bmin, ASCI=asci_flag)
!
       fmeb   = label%real (12)
       fstoff = label%real (13)
       fdst   = label%real (14)
       fclps  = label%real (15)
       nold_bf   = nbf + 1
     END IF
!
  ELSE 

     nn = -999
     find_loop: DO n = 2, nbf
        IF (jtime <= ibtime(n)) THEN
           nn = n
           EXIT find_loop
        END IF
     END DO find_loop
     IF (nn == -999) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
         'ibtime screwed up, stop in bfield')
!
     IF (nn /= nold_bf) THEN
!
        LREC = nn - 1
        CALL Read_array ('input/rcmxmin_inp', LREC , label, &
                          ARRAY_2D = xmin_1, ASCI=asci_flag)
        CALL Read_array ('input/rcmymin_inp', LREC , label, &
                          ARRAY_2D = ymin_1, ASCI=asci_flag)
        CALL Read_array ('input/rcmvm_inp',   LREC , label, &
                          ARRAY_2D = vm_1, ASCI=asci_flag )
        CALL Read_array ('input/rcmbmin_inp', LREC , label, &
                          ARRAY_2D = bmin_1, ASCI=asci_flag)
!
        fmeb1   = label%real (12)
        fstoff1 = label%real (13)
        fdst1   = label%real (14)
        fclps1  = label%real (15)
!
        LREC = nn
        CALL Read_array ('input/rcmxmin_inp', LREC , label, &
                          ARRAY_2D = xmin_2, ASCI=asci_flag)
        CALL Read_array ('input/rcmymin_inp', LREC , label, &
                          ARRAY_2D = ymin_2, ASCI=asci_flag)
        CALL Read_array ('input/rcmvm_inp',   LREC , label, &
                          ARRAY_2D = vm_2, ASCI=asci_flag )
        CALL Read_array ('input/rcmbmin_inp', LREC , label, &
                          ARRAY_2D = bmin_2, ASCI=asci_flag)
!
        fmeb2   = label%real (12)
        fstoff2 = label%real (13)
        fdst2   = label%real (14)
        fclps2  = label%real (15)
!
        nold_bf = nn
     END IF
!
     f = REAL(jtime-ibtime(nn-1), rprec) / &
         REAL(ibtime(nn)-ibtime(nn-1), rprec)
     xmin   = (1.0_rprec-f)*xmin_1 + f*xmin_2
     ymin   = (1.0_rprec-f)*ymin_1 + f*ymin_2
     bmin   = (1.0_rprec-f)*bmin_1 + f*bmin_2
     vm     = (1.0_rprec-f)*vm_1 + f*vm_2
     fstoff = (1.0_rprec-f)*fstoff1+f*fstoff2
     fmeb   = (1.0_rprec-f)*fmeb1+f*fmeb2
     fdst   = (1.0_rprec-f)*fdst1+f*fdst2
     fclps  = (1.0_rprec-f)*fclps1+f*fclps2
!
  END IF
!
  rmin = SQRT (xmin**2+ymin**2)
  pmin = ATAN2 (ymin, xmin)
!
  RETURN
ELSE IF (itype_bf == 3) THEN
! Get Bfield from somewhere

ELSE
print*,itype_bf
   call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
       'ILLEGAL BFIELD TYPE IN GET_BFIELD')
END IF
      END SUBROUTINE Get_bfield
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_eta_on_bndy (jtime )
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
!_____________________________________________________________________________
!
!_____________________________________________________________________________
!
!
      INTEGER (iprec), SAVE :: nold_t = 0
      INTEGER (iprec) :: n, nn, Lrec, n_t, kc, j, i, i_least
      REAL    (rprec) :: f
      REAL    (rprec), SAVE :: etac_1 (kcsize), etac_2 (kcsize)
!
!
      n_t = SIZE (itime_etac(:))
!
      IF (jtime <= itime_etac(1)) THEN
!
         etac = etac_inp (:,1)
!
      ELSE IF (jtime > itime_etac(n_t)) THEN
!
         etac = etac_inp (:,n_t)
!
      ELSE 
         nn = -999
         find_loop: DO n = 2, n_t
            IF (jtime <= itime_etac(n)) THEN
               nn = n
               EXIT find_loop
            END IF
         END DO find_loop
         IF (nn == -999) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
             'ibtime screwed up, stop in get_eta_on_bndy')
!
         f = REAL(jtime-itime_etac(nn-1), rprec) / &
             REAL(itime_etac(nn)-itime_etac(nn-1), rprec)
         etac   = (1.0_rprec-f)*etac_inp(:,nn-1) + f*etac_inp(:,nn)
!
      END IF
!
!
      i_least = MIN(MINVAL(imin_j)-2,1)
!
      DO kc = 1, kcsize
         DO j = 1, jsize
            DO i = i_least, imin_j(j)
               eeta(i,j,kc) = etac(kc)
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE Get_eta_on_bndy
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Get_boundary (boundary, bndloc)
    USE Rcm_variables, junk_boundary=>boundary, junk_bndloc=>bndloc
    IMPLICIT NONE
    TYPE (ellipse_def), DIMENSION (2), INTENT (OUT) :: boundary
    REAL (rprec), INTENT (OUT) :: bndloc (1-n_gc:jsize+n_gc)
!_____________________________________________________________________________
!   This function must be called AFTER calling BFIELD.
!
!   IBND_TYPE is an RCM control parameter, read from input file. It specifies
!             how to set the RCM high-latitude boundary.
!
!   Note that subroutine returns the parameters of ellipse,
!   but for different IMODE they have different meaning, namely:
!   IBND = 1:  AA, BB, XC, YC are in Re units, with XC, YC
!              being GSM coordinates, so XC>0 on dayside,
!              YC positive on duskside
!   IBND = 2:  AA, BB, XC, YC are in degrees measured from
!              north pole, displacements being >0 toward noon,
!              and dusk (same convention as in EFIELD).
!_____________________________________________________________________________
!
    INTEGER (iprec) :: j, jind (3), n, ierr,i,lp,lpMax
    REAL (rprec) :: a, b, c, f_c, colat_bnd, theta (3), bndloc_tmp, &
                    Gntrp_2d, Thet, R_12, R_24
    REAL (rprec), PARAMETER :: phi (3) = (/ 0.0, pi_by_two, pi /), &
!                                                  |     |        |
!                                                noon, dusk,    midnight
!
                       fdist (3) = (/0.95_rprec, 1.45_rprec, 1.9_rprec/)


    bndloc = -1.0
    imin_j = -1

    SELECT CASE (ibnd_type)

       CASE (1)

         !  This gives ellipse extending to 0.95*fstoff at noon,
         !  2*fstoff at midnight, and 1.5*fstoff at dawn and dusk:

         boundary(1)%aa =  1.475_rprec * fstoff
         boundary(1)%bb =  1.500_rprec * fstoff
         boundary(1)%xx = -0.525_rprec * fstoff
         boundary(1)%yy =  0.0_rprec


         DO j = 1-n_gc, jsize+n_gc ! place ellipse on the grid:

            a = 1.0_rprec
            b = REAL (isize, rprec)
            DO lp=1,20
               c   = 0.5 * (a+b)
!              IF ( ABS ( Fequat_of_x (c, REAL(j,rprec)) ) < 100*EPSILON(1.0_rprec)) EXIT
!              IF (       Fequat_of_x (c, REAL(j,rprec)) < 0.0_rprec) THEN
               IF (       Fequat_of_x (c, REAL(j)) < 0.0_rprec) THEN
                   b = c
               ELSE
                   a = c
               END IF
            END DO

            bndloc(j) = c

         END DO

         imin_j = CEILING (bndloc) ! first grid point inside modeling region.


       CASE (4)
      

         ! weird ellipse for SWMF

         ! The maximum distance for the tip of the ellipse towards the midnight meridian (-X)
         R_24 = 10.0
  
         ! The maximum distance for the tip of the ellipse towards the noon meridian (+X)
         R_12 = maxval(xmin(1:iSize,1:jSize),MASK=vm(1:iSize,1:jSize)>0.0)

         IF (R_12 < 1.0) THEN
            write(*,*)'R_24=',R_24,' R_12=',R_12
            write(*,*)'vm(:,1)=',vm(:,1)
            write(*,*)'xmin(:,1)=',xmin(:,1)
            call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:WRONG R_12')
         END IF

         boundary(1)%aa =  (R_12+R_24)/2.0
         boundary(1)%bb =  1.0 * min(R_12,R_24)
         boundary(1)%xx =  (R_12 - R_24) / 2.0
         boundary(1)%yy =  0.0_rprec

         ! Check if ellipse will fit in grid, shrink if neccessary
         lpMax=100
         ChkLoop: do lp=1,lpMax
            do j=1-n_gc, jsize+n_gc
               do i=isize,1,-1
                  if(vm(i,j)>0.0 .AND. vm(i-1,j)<0.0) then
                     IF ( (xmin(i,j)-boundary(1)%xx)**2 / boundary(1)%aa**2 + &
                          (ymin(i,j)-boundary(1)%yy)**2 / boundary(1)%bb**2 &
                          <= 1.0_rprec) then
                        if(lp==lpMax)then
                           write(*,*)'Error in grid, shrinking ellipse:','   i,j=',i,j
                           write(*,*)' Point=',xmin(i,j),ymin(i,j),vm(i,j)
                           write(*,*)' Ellipse:', &
                                boundary(1)%xx,boundary(1)%yy,boundary(1)%aa,boundary(1)%bb
                           call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90')
                        end if
                        if(xmin(i,j)>0.)then
                           r_12=0.95*r_12
                        else
                           r_24=0.95*r_24
                        end if
                        boundary(1)%aa =  (R_12+R_24)/2.0
                        boundary(1)%bb =  1.0 * min(R_12,R_24)
                        boundary(1)%xx =  (R_12 - R_24) / 2.0
                        boundary(1)%yy =  0.0_rprec
                        CYCLE ChkLoop
                     end IF
                     EXIT
                  end if
               end do
            end do
            call IM_write_prefix
            write(iUnitOut,*)' Ellipse size adjusted ',lp,' times. size=',R_12,R_24
            EXIT ChkLoop
         end do ChkLoop


         DO j = 1-n_gc, jsize+n_gc ! place ellipse on the grid:

            a = 1.0_rprec
            b = REAL (isize, rprec)
            DO lp=1,20
               c   = 0.5 * (a+b)
!              IF ( ABS ( Fequat_of_x (c, REAL(j,rprec)) ) < 100*EPSILON(1.0_rprec)) EXIT
!              IF (       Fequat_of_x (c, REAL(j,rprec)) < 0.0_rprec) THEN
               IF (       Fequat_of_x (c, REAL(j)) < 0.0_rprec) THEN
                   b = c
               ELSE
                   a = c
               END IF
            END DO

            bndloc(j) = c

         END DO

         imin_j = CEILING (bndloc) ! first grid point inside modeling region.

       CASE (3)

         !  boundary at L=6.6 in the magn. equat. plane (10/24/200)

         boundary(1)%aa =  6.6_rprec
         boundary(1)%bb =  6.6_rprec
         boundary(1)%xx =  0.0_rprec
         boundary(1)%yy =  0.0_rprec


         DO j = 1-n_gc, jsize+n_gc ! place ellipse on the grid:

            a = 1.0_rprec
            b = REAL (isize, rprec)
            DO lp=1,20
               c   = 0.5 * (a+b)
!              IF ( ABS ( Fequat_of_x (c, REAL(j,rprec)) ) < 100*EPSILON(1.0_rprec)) EXIT
!              IF (       Fequat_of_x (c, REAL(j,rprec)) < 0.0_rprec) THEN
               IF (       Fequat_of_x (c, REAL(j)) < 0.0_rprec) THEN
                   b = c
               ELSE
                   a = c
               END IF
            END DO
            
            bndloc(j) = c

         END DO

         imin_j = CEILING (bndloc) ! first grid point inside modeling region.

       CASE (2)

         !  Set up boundary in the ionosphere to an ellipse such that
         !  it maps out to the equatorial plane to given locations at
         !  noon, midnight and dusk.

         ! First, find colatitudes of 3 specified points in eq. plane:

         DO n = 1, 3 ! Since the ellipse will be in the ionosphere, map the three points:

            jind (n) = -1
            DO j = 1, jsize
              IF ( ABS(aloct(1,j)-phi(n)) < EPSILON(1.0_rprec)) THEN
                 jind (n) = j
              END IF
            END DO

            IF (jind(n) == -1) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
              'UNABLE TO LOCATE ONE OF PNTS IN GETBND')


            a = REAL (isize,rprec)
            b = REAL (1.0_rprec, rprec)
            DO
              c = 0.5 * (a+b)
              f_c = Gntrp_2d (rmin, isize, jsize, c, REAL(jind(n),rprec)) -fstoff*fdist(n)
              IF (ABS (f_c) < 10.0*EPSILON(1.0_rprec)) EXIT
              IF (f_c < 0.0_rprec) THEN
                 a = c
              ELSE
                 b = c
              END IF
            END DO
            theta (n) = Gntrp_2d (colat, isize, jsize, c, REAL(j,rprec))
         END DO

         !  Compute parameters of ellipse:

         boundary(2)%bb = theta(2)
         boundary(2)%xx = 0.5 * (theta(1) - theta(3))
         boundary(2)%aa = theta(1) - boundary(2)%xx
         boundary(2)%yy = 0.0_rprec


         !  From colatitudes of boundary points, estimate their I-values:

         DO j = 1-n_gc, jsize+n_gc
           colat_bnd = Thet (boundary(2)%aa, boundary(2)%bb, &
                             boundary(2)%xx, boundary(2)%yy, aloct(1,j) )
           a = 1.0_rprec
           b = REAL (isize,rprec)
           DO
             c = 0.5 * (a+b)
             f_c = Gntrp_2d (colat, isize, jsize, c, REAL(j,rprec)) - colat_bnd
             IF (ABS (f_c) < 10.8*EPSILON(1.0_rprec)) EXIT
             IF (f_c < 0.0_rprec) THEN
               a = c
             ELSE
               b = c
             END IF
           END DO
           bndloc(j) = c
         END DO

         imin_j = CEILING (bndloc) ! first grid point inside modeling region.

         boundary(2)%aa = boundary(2)%aa * RTD
         boundary(2)%bb = boundary(2)%bb * RTD
         boundary(2)%xx = boundary(2)%xx * RTD
         boundary(2)%yy = boundary(2)%yy * RTD

       CASE (5)

         ! set boundary close to open/closed field line boundary:
  
         DO j = 1 - n_gc, jsize + n_gc
            do_i_loop: DO i = isize - 1, 1, -1
               if (vm(i,j) > 0.0 .AND. vm(i-1,j) < 0.0) THEN
                  imin_j(j) = i
                  bndloc(j) = i
                  EXIT do_i_loop
               end if 
            END DO do_i_loop

            if (imin_j(j) < 1) then
              call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:&
                            &boundary failed to set')
            end if
         END DO

         imin_j = CEILING (bndloc) ! first grid point inside modeling region.

       CASE DEFAULT

          call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
           'ILLEGAL VALUE OF IBND_TYPE') 

    END SELECT

    RETURN

    CONTAINS
!
      FUNCTION Fequat_of_x (bi, bj)
      IMPLICIT NONE
      REAL (KIND=rprec), INTENT (IN) :: bi, bj
      REAL (KIND=rprec) :: Fequat_of_x
      REAL (KIND=rprec) :: xx, yy, Gntrp_2d
      xx  = Gntrp_2d (xmin, isize, jsize, bi, bj)
      yy  = Gntrp_2d (ymin, isize, jsize, bi, bj)
      Fequat_of_x = (xx-boundary(1)%xx)**2 / boundary(1)%aa**2 + &
                    (yy-boundary(1)%yy)**2 / boundary(1)%bb**2 - 1.0_rprec
      RETURN
      END FUNCTION Fequat_of_x
!
    END SUBROUTINE Get_boundary
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_jbirk (jtime)
      USE Rcm_variables, junk_ss => ss
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
!__________________________________________________________________________
!                                                                       
!  Program written by: r.w. spiro        
!  last update:
!     04-05-88          
!     01-29-96 frt                 - added ain,min_j arr
!  Algorithm by: r.a. wolf                                              
!                                                                       
!  This subroutine computes birk(i,j) given inner edge
!  locations 
!  modified 04-05-88 to include effects of gradients in eta.            
!  see raw document re including eeta in computation of jbirk           
!    dated feb 6, 1988.                                                 
!  birk is current density (2 hemispheres) in units of
!  microamp/m**2    
!
!  birk(i,j) here is the field-aligned current density per
!  unit of ionospheric area, so that it already includes
!  the factor sin(I); this is J_parallel*sin(I) in the RHS
!  of the Vasyliunas equation.
!
!  Issues with non-integer boundary (Stanislav's notes):
!  for BIRK from inner edge segments, this is not an issue
!  (except that if a segment is entirely outside the bndry,
!  then we don't compute its contribution); of course, we 
!  have to care about this somewhere else where motion of 
!  test particles is computed. For BIRK from gradients of 
!  EETA,  
!                                                                       
!______________________________________________________________________________
!
      REAL (rprec), PARAMETER :: cf1 = 1.0_rprec / pi_two, &
                cf2 =  - (3.0_rprec/4.0_rprec)*( (2.0_rprec / pi) - 0.5_rprec)
!
      INTEGER (iprec) :: i, j, k, kc, klbeg, klend, kl, klnext, &
                 ibmin, ibmax, jbmin, jbmax, jb1, jb2,  &
                 ig, jj, jindex, ib1, ib2, &
                 Bjmod_int
      REAL (rprec), DIMENSION (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc) :: &
                      detadi, detadj, &
                      dvmdi, dvmdj
      REAL (rprec) :: dbirk, &
                      vmkl, vmnext, sum, b1, b2, x, y, el, umax, umin, ss, &
                      z, dg1, dg2, dg3, qmin, qmax, qn, qx,  &
                      denom, a1, a2, bjm, range, bim, gkl (5000), &
                      Gntrp_2d, Bndy
!                                                                       
!
      birk  = 0.0_rprec
!
      DO_20: DO k = 1, ksize
!
        IF (ABS(eta(k)) < TINY(1.0_rprec)) CYCLE DO_20  ! in case edge is "empty"
!
        klbeg = mpoint (k) 
        klend = klbeg + npoint (k) - 1 
!                                                                       
        vmnext = Gntrp_2d (vm,  isize, jsize, bi (klbeg), bj (klbeg))
!                                                                       
!
        DO_30: DO kl = klbeg, klend   ! goes over pts on edge K
!
          vmkl = vmnext 
          klnext = KL + 1 
          IF (klnext > klend) klnext = klbeg 
          vmnext = Gntrp_2d (vm, isize, jsize, bi (klnext), bj (klnext))
!                                                             
!
!         Determine nearby grid pts:                                   
!
          ib1   = INT (bi (kl))
          ib2   = INT (bi (klnext))
          ibmin = MAX (MIN (ib1, ib2), 1 )
          ibmax = MIN (MAX (ib1 + 1, ib2 + 1), isize )
!
!         Skip this segment if it is outside the modeling region:
          IF (ibmax       <= 1                     .OR. &
              (bi(kl)     <= Bndy (bndloc, jsize, bj(kl)) .OR. &
               bi(klnext) <= Bndy (bndloc, jsize, bj (klnext) ) )) CYCLE DO_30
!                                                                       
!
          CALL Adjust_bj_2 (bj (kl), bj (klnext), b1, b2 )
          jb1 = b1 
          jb2 = b2 
          jbmin = MIN (jb1, jb2) 
          jbmax = MAX (jb1 + 1, jb2 + 1) 
!
          sum = 0.0_rprec
          ig = 0 
!                                                                       
!
!
          DO i = ibmin, ibmax ! loop over all nearby grid points
          DO j = jbmin, jbmax 
!
            ig = ig + 1 
            IF (ig > SIZE(gkl)) THEN
               WRITE (*,'(T2,A)') 'JBIRK: ig.gt.igdim'
               call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90')
            END IF
!
!           Consider straight line  bi=a1*bj+a2
!                                                                       
            denom = b1 - b2 
            IF (ABS (denom) < 1E-5) denom = b1-b2-1E-5
!
!           slope=a1; intercept=a2:
            a1 = (bi (kl) - bi (klnext) ) / denom  
            a2 = (b1 * bi (klnext) - bi (kl) * b2) / denom   
!                                                                       
!
!           bjm= j value of point that is closest to (i,j)
!           and lies on the line that passes through kl and klnext:
!                                                                       
            bjm = (j - a1 * (a2 - i) ) / (1.0_rprec + a1**2)
!                                                                       
!           gkl(ig)="nearness parameter" used to compute             
!              weighting factors for nearby grid points              
!                                                                       
!           Let                                                      
!           l = infinitely long straight line thru points kl, klnext.
!           p = point on line l that is closest to grid point i,j.   
!           x = distance from p to midpoint of segment (kl,klnext).  
!           y = distance from p to grid point (i,j).                 
!           el = length of segment (kl, klnext).                     
!                                                                       
!           j parallel from the segment is distributed among         
!           grid points that are within distance 'range' of segment.
!                                                                       
!           bim = bi-value on line l that corresponds to bj=bjm.     
!                                                                       
            range = 1.0_rprec
            bim = a1 * bjm + a2 
            x  = (SQRT((bim - 0.5 * (bi(kl)+bi(klnext)))**2 + &
                       (bjm - 0.5 * (b1+b2))**2) ) / range
            y  = (SQRT((REAL(i,rprec)-bim)**2+(REAL(j,rprec)-bjm)**2))/range
            el = (SQRT((bi(klnext)-bi(kl))**2+(b1-b2)**2))/range
!                                                                       
            IF (y >= 1.0_rprec) THEN
               umax = 0.0_rprec
               umin = 0.0_rprec
            ELSE 
               ss   = SQRT (1.0_rprec - y**2)
               umax = MIN ( (0.5 * el - x), ss)
               umin = MAX ( ( -0.5 * el - x), - ss)
               z    = MAX (y, 0.00001_rprec)
            END IF 
!                                                                       
            IF (umax <= umin) THEN 
               gkl (ig) = 0.0_rprec
            ELSE 
               dg1 = cf1 * ( &
                  -umax*LOG(umax**2 + z**2) + umin * LOG (umin**2 + z**2) &
                  + 2.0_rprec*(umax - umin) - &
                    2.0_rprec*y*(ATAN (umax/z) - ATAN (umin/z)))
!                                                                       
               dg2 = cf2 * ( &
                 (1.0_rprec - y**2) * (umax - umin) - (umax**3 - umin**3) / 3.0)
!                                                                       
                     qmax = MAX (0.0_rprec, 1.0_rprec - y**2 - umax**2)
                     qmin = MAX (0.0_rprec, 1.0_rprec - y**2 - umin**2)
                     IF (ss <= 0.0_rprec)  ss = 1.0E-6
                     qx = umax / ss 
                     IF (ABS (qx) > 1.0_rprec) qx = SIGN (1.0_rprec, qx)
                     qn = umin / ss 
                     IF (ABS(qn) >  1.0_rprec) qn = SIGN (1.0_rprec, qn)
               dg3 = cf1 * (umax * SQRT (qmax) - umin * SQRT (qmin) &
                      + ss**2 * (ASIN (qx) - ASIN (qn) ) )
               gkl (ig) = dg1 + dg2 + dg3 
            END IF 
!                                                                       
            IF (gkl (ig) < 0.) then 
!!7/13/99      WRITE (*,'(A,E10.2,4I5)') &
!!7/13/99     'gkl negative. gkl,ig,i,j,k =',gkl(ig),ig,i,j,k
            END IF 
!
            sum = sum + gkl (ig) 
!                                                                       
          END DO
          END DO
!                                                                       
          IF (ABS(sum) < TINY(1.0_rprec)) CYCLE DO_30
          ig = 0 
!
          DO i = ibmin, ibmax 
          DO jindex = jbmin, jbmax 
             ig = ig + 1 
             j = jindex 
!
             jj = Bjmod_int (j, jsize)
!
             IF (ABS(gkl (ig)) > 0.0_rprec) THEN
!
!              Compute contribution to jbirk(i,j)due to segment kl:klnext:
!
               dbirk = charge_e * signbe * (vmnext - vmkl) *  &
                       eta(k) * ABS(alam(k))*sgn(k)* &
                       (gkl (ig) / sum) / &
                       (alpha (i, jj) * beta (i, jj) * dlam * dpsi * ri**2)
!                                                                       
!              Birkeland current caused by the centrifugal force here:
!              dbirk = dbirk -  &
!                      0.5e6*signbe*romeca**2.0* &
!                      (rnext-rkl)*(rnext+rkl)* 
!                      xmass(ikflav(k))*eta(k)*sgn(k)*(gkl(ig)/sum)/ 
!                      (alpha(i,jj)*beta(i,jj)*dlam*dpsi)
!
               birk (i, jj) = birk (i, jj) + dbirk 
!
             END IF 
          END DO
          END DO
!                                                                       
        END DO DO_30
      END DO DO_20
!                                                                       
!
!
!
!     Compute J_parallel due to continuous channel:
!                                                                       
      CALL Deriv_i (vm, isize, jsize, imin_j, dvmdi)
      CALL Deriv_j (vm, isize, jsize, imin_j, 1.0E+25, dvmdj)
      WHERE (ABS(dvmdj) > 1.0E+24)  ! to prevent artificial inflows on bndy
          dvmdi = 0.0_rprec
          dvmdj = 0.0_rprec
      END WHERE 
!
      DO kc = 1, kcsize
!
         CALL Deriv_i (eeta (:,:,kc), isize, jsize, imin_j, detadi)
         CALL Deriv_j (eeta (:,:,kc), isize, jsize, imin_j, 1.0E+32, detadj)
         WHERE (ABS(detadj) > 1.0E+31)
           detadi = 0.0
           detadj = 0.0
         END WHERE
         DO  j = 1, jsize
         DO  i = 1, isize
            IF (i < imin_j(j) + 1) CYCLE
               dbirk  = charge_e * signbe * ABS(alamc(kc)) * &
                        (detadj(i,j) * dvmdi(i,j) - detadi(i,j) * dvmdj(i,j)) /  &
                        (alpha(i,j)*beta(i,j)*dlam*dpsi*Ri**2)
               birk (i, j) = birk (i, j) + dbirk 
            END DO 
         END DO 
      END DO 
!
!                                                                       
     CALL Wrap_around_ghostcells (birk, isize, jsize, n_gc)
!
      RETURN 
      END SUBROUTINE Get_jbirk
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_vparallel ()
      USE Rcm_variables
      IMPLICIT NONE
!______________________________________________________________________________
!  last update: 
!     05-05-87       by:rws                              
!     02-10-96          frt - added arrays ain,min_j     
!                                                                       
!  Birk is sum of current densities into both hemispheres.              
!  (micro amp/m**2).  Before activating parallel potential drop        
!  we need to check if birk is being used correctly in
!  this routine.   
!
!  Stanislav: VPAR is computed inside IE loop (for both
!             negative and positive particles), and will
!             be the one for the largest IE value. Which
!             is nonsense.
!  Stanislav: this subroutine needs grid-based formulation
!             of plasma (EETA). Before it was done by
!             computing EETA for electrons from the inner
!             edges of electrons, then it was changed to
!             use directly grid-based population. In the
!             latter case, array PVEC returned by this
!             routine is the electron pressure (without
!             the factor of 2/3) and is the same as what
!             routine PV returns as array PVGAM. If the
!             electrons are on the grid only, as in my case,
!             then we call PV in rcm main program to compute
!             the ion pressure, and we use PVEC from this
!             routine for the electron pressure. (04/20/99)
!  Stanislav, may 18,99: make all loops over electrons only,
!             by using iedim_local and setting it to 1.
!______________________________________________________________________________
!
      INTEGER (iprec) :: i, j, ie, iedim_local, kc
      REAL (rprec)    :: en, ekt, therm, sum1 (iesize), sum2 (iesize)
!                                                                       
!
!                                  
      iedim_local = 1
!
      vpar  (:,:)   = 0.0_rprec 
      eavg  (:,:,:) = 0.0_rprec
      eflux (:,:,:) = 0.0_rprec
!
      loop_j: DO j = 1, jsize
      loop_i: DO i = 1, isize
      IF (i <imin_j(j)) CYCLE
!
!           For each grid point, clear sum1 and sum2:
!
            sum1  = 0.0_rprec
            sum2  = 0.0_rprec
!
!
!           Now for each grid point, consider all species
!           present at that grid point, and compute sum1 and
!           sum2 for positive and negative particles separately:
!
            GRID_BASED: DO kc = 1, kcsize
   IF (fudgec(kc) == 0.0) CYCLE
             IF ( ABS(alamc(kc))*vm(i,j) > 500.0_rprec) THEN
               IF (alamc (kc) < 0.0_rprec) THEN
                  ie = 1 
               ELSE 
                  ie = 2 
!                 STOP 'BALGN4: ie is 2'
               END IF
               sum1(ie) = sum1(ie) + eeta(i,j,kc)*fudgec(kc)
               sum2(ie) = sum2(ie) + eeta(i,j,kc)*fudgec(kc)*ABS(alamc(kc))
             END IF
            END DO GRID_BASED 
!
!           For positive and negative particles separately,
!           compute precipitating number flux, average energy,
!           and parallel potential drop:
!
            DO ie = 1, iedim_local 
!                                                                       
               IF (sum1 (ie) > 0.0_rprec) THEN
!
!                compute thermal electron current, field-aligned
!                potential drop, electron energy flux,
!                and average electron energy at (i,j):          
!
                  en    = sum1 (ie) * vm (i, j)**1.5 / 6.38E+21
                  ekt   = (2.0_rprec/3.0_rprec) * sum2 (ie) * vm (i,j) / sum1 (ie)
                  therm = 0.02675 * en * SQRT(ekt*xmass(1)/xmass(ie))
!
                  IF (therm < 1.E-30) THEN 
                     therm      = 0.0_rprec
                     vpar (i,j) = 0.0_rprec
                     eflux(i,j,ie) = 0.0_rprec
                     eavg(i,j,ie) = 1000.
                  ELSE 
                     IF (- birk (i, j) / therm > 1.0_rprec ) THEN
                        vpar (i,j) = ekt * (- birk (i,j) / therm - 1.0_rprec)
                     ELSE 
                        vpar (i,j) = 1.0_rprec
                     END IF
                     vpar(i,j) = MIN (vpar (i, j), 10000.0_rprec)
!
!    !!!!!!!      ALERT: VPAR(I,J) IS SET TO 0 !!!!!!!!!!!!!!!!!
!
                  vpar (i, j) = 0.0_rprec
!
!
                  eflux(i,j,ie) = 0.002 * therm * &
                                 ( ekt + vpar(i,j) + 0.5*vpar(i,j)**2/ekt)
                  eavg(i,j,ie) = 2.0_rprec*(ekt + &
                                 vpar(i,j)+0.5*vpar(i,j)**2 /ekt) / &
                                 (1.0_rprec + vpar (i, j) / ekt)
                  END IF 
               ELSE 
!                                                                       
!                 Case fudge=0: we want eflux=0 and eavg=0 for no precipitation.
!
                  eflux (i, j, ie) = 0.0_rprec
                  eavg  (i, j, ie) = 0.0_rprec
!
               END IF 
!                                                                       
            END DO
!
      END DO loop_i
      END DO loop_j 
!                                                                       
!
      CALL Wrap_around_ghostcells (vpar, isize, jsize, n_gc)
      CALL Wrap_around_ghostcells (eflux(:,:,ie_ele), isize, jsize, n_gc)
      CALL Wrap_around_ghostcells (eavg (:,:,ie_ele), isize, jsize, n_gc)
!
      RETURN
      END SUBROUTINE Get_vparallel
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Floor_for_eflux ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: ivalue_max, i, j
      REAL (rprec)    :: eflux_max
      DO j = 1, jsize
         eflux_max = eflux (isize, j, ie_ele)
         ivalue_max = isize
         DO i = isize-1, imin_j(j), -1
            IF (eflux(i,j,ie_ele) > eflux(i+1,j,ie_ele)) THEN
               eflux_max  = eflux(i,j,ie_ele)
               ivalue_max = i
            END IF
         END DO
         DO i = imin_j(j), ivalue_max - 1
            eflux(i,j,ie_ele) = MAX (0.5*eflux_max, eflux(i,j,ie_ele))
         END DO
      END DO
      CALL Wrap_around_ghostcells (eflux(:,:,ie_ele),isize,jsize,n_gc)
      RETURN
      END SUBROUTINE Floor_for_eflux
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Correct_eflux ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: i, j, ikp_low, ikp_high
      REAL    (rprec) :: hardy_eflux_int, eflux_int, value_l, value_h, value,&
                         factor
!
!     IF (kp > 5.99) STOP ' kp too large in correct_eflux'
      ikp_low = INT (kp)
      ikp_high = ikp_low + 1
      IF (ikp_high > 6) THEN
         ikp_high = 6
         ikp_low  = 5
      END IF
      factor   = (kp - ikp_low)/(ikp_high-ikp_low)
      DO j = 1, jsize
!
!           1. Compute latitudinal integral of hardy's EFLUX:
!
            hardy_eflux_int = 0.0_rprec
            DO i = 2, isize-1
               CALL Elemod (1, ikp_low, &
                            90.0-colat(i,j)*RTD, &
                            MODULO (12.0+aloct(i,j)*RTH,24.0), &
                            value_l)
               CALL Elemod (1, ikp_high, &
                            90.0-colat(i,j)*RTD, &
                            MODULO (12.0+aloct(i,j)*RTH,24.0), &
                            value_h)
               value = value_l*(1.0_rprec-factor)+value_h*factor
               value = (10.0**value) *1.6E-09 * pi
               hardy_eflux_int = hardy_eflux_int + value*(colat(i,j)-colat(i-1,j))
            END DO
!
!
!           2. Compute latitudinal integral of uncorrected RCM's EFLUX:
!
            eflux_int  = 0.0_rprec
            DO i = imin_j(j), isize - 1
               eflux_int = eflux_int + eflux(i,j,1)*(colat(i,j)-colat(i-1,j))
            END DO
!
!
!           3. Make correction:
!
            IF (eflux_int > 0.0) THEN
               DO i = imin_j(j), isize
                  eflux (i,j,1) = eflux(i,j,1)*(0.5-0.3*SIN(aloct(1,j)))*&
                                  hardy_eflux_int / eflux_int
               END DO
            ELSE IF (eflux_int == 0.0) THEN
               IF (ANY(eflux(imin_j(j):isize,j,1) /=0.0)) THEN
                  call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
                      'EFLUX_INT = ZERO BUT EFLUX IS NOT')
               ELSE
                  eflux (imin_j(j):isize,j,1) = eflux (imin_j(j):isize,j,1)
               END IF
            ELSE
               call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
                   'EFLUX_INT IS NEGATIVE')
            END IF
!
      END DO
!
      DO i = 1, isize
      DO j = 1, jsize
         IF (eavg(i,j,ie_ele) < 1.0E-30) eflux(i,j,ie_ele) = 0.0
      END DO
      END DO
!
      CALL Wrap_around_ghostcells (eflux(:,:,ie_ele),isize,jsize,n_gc)
      CALL Wrap_around_ghostcells (eavg(:,:,ie_ele), isize,jsize,n_gc)
      RETURN
      END SUBROUTINE Correct_eflux
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_hardy_cond ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: i,j, kp_low, kp_high
      REAL    (rprec) :: value, value_low, value_high, factor 
!
      kp_low = INT (kp)
      kp_high = kp_low + 1
      factor = (kp-kp_low)/(kp_high-kp_low)
      DO i = 1, isize
      DO j = 1, jsize
         CALL Elemod (3, kp_low, 90.0-colat(i,j)*RTD, &
                      MODULO(12.+aloct(i,j)*RTH,24.0), value_low)
         CALL Elemod (3, kp_high, 90.0-colat(i,j)*RTD, &
                      MODULO(12.+aloct(i,j)*RTH,24.0), value_high)
         value = value_low*(1.0_rprec-factor)+value_high*factor
         hall (i,j) = qthall(i,j) + 2.0_rprec*value / sini(i,j)
         CALL Elemod (4, kp_low, 90.0-colat(i,j)*RTD, &
                      MODULO(12.+aloct(i,j)*RTH,24.0), value_low)
         CALL Elemod (4, kp_high, 90.0-colat(i,j)*RTD, &
                      MODULO(12.+aloct(i,j)*RTH,24.0), value_high)
         value = value_low*(1.0_rprec-factor)+value_high*factor
         pedpsi(i,j) = qtped(i,j) + 2.0_rprec*value
         pedlam(i,j) = qtplam(i,j) + 2.0_rprec*value/sini(i,j)**2
      END DO
      END DO
      CALL Wrap_around_ghostcells (pedpsi, isize, jsize, n_gc)
      CALL Wrap_around_ghostcells (pedlam, isize, jsize, n_gc)
      CALL Wrap_around_ghostcells (hall  , isize, jsize, n_gc)
      RETURN
      END SUBROUTINE Get_hardy_cond
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_active_cond ( )
      USE Rcm_variables
      IMPLICIT NONE
!
!______________________________________________________________________________
!
!  This subroutine calculates conductance enhancement
!  due to auroral electron precipitation and adds this 
!  to the quiet time conductances read in subroutine 
!  qtcond. This subroutine contains changes made for 
!  tjfr run and corrected formulas for conductances 
!  as put forth by robinson         
!  last update: 11-06-86                                                
!               02-07-96 frt - min_j array added                        
!                                                                       
!  Stanislav, april 14 1999: added arrays for precipitation
!             conductances so that they are smoothed and 
!             quiet-time conductances are not modified.
!             This was also accompanied by a change in BALGN4
!             that now runs from i = 1 not min_j.
!
!______________________________________________________________________________
!
      INTEGER (iprec) :: i, j
      REAL    (rprec) :: ezero, sigp, sigh
!
!
      IF (ifloor)   CALL Floor_for_eflux ()
      IF (icorrect) CALL Correct_eflux   ()
!
      pedpsi = 0.0_rprec
      pedlam = 0.0_rprec
      hall   = 0.0_rprec
!
!
      DO j = 1, jsize
!     DO i = Get_imin_for_grid(j)-1, isize
      DO i = 1, isize
      IF ( i < imin_j(j)) CYCLE
        IF (eflux(i,j,ie_ele) > 1.0E-6 .AND. eavg(i,j,ie_ele) < 1.E-5) THEN
           WRITE (*,*) 'stopping in cond, see the code'
           WRITE (*,*) i,j,eflux(i,j,1),eavg(i,j,1), imin_j(j)
           call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90')
        END IF 
        ezero = eavg (i, j, ie_ele) / 1.0E3
        sigp  = SQRT(eflux(i,j,ie_ele)) * 40.0 * ezero / (16.0 + ezero**2)
        sigh  = 0.45 * sigp * ezero**(0.85)
        pedpsi (i, j) = 2.0_rprec * sigp
        pedlam (i, j) = 2.0_rprec * sigp / (sini(i,j)**2)
        hall   (i, j) = 2.0_rprec * sigh / sini (i, j)
      END DO
      END DO
!
     CALL Wrap_around_ghostcells (pedlam, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (pedpsi, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (hall, isize, jsize, n_gc)
!
      CALL Smooth_j (pedpsi)
      CALL Smooth_j (pedlam)
      CALL Smooth_j (hall)
!                                                                       
      CALL Smooth_i (pedpsi)
      CALL Smooth_i (pedlam)
      CALL Smooth_i (hall)
     CALL Wrap_around_ghostcells (pedlam, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (pedpsi, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (hall, isize, jsize, n_gc)
!
      pedpsi (:,:) = pedpsi (:,:) + qtped (:,:)
      pedlam (:,:) = pedlam (:,:)+ qtplam (:,:)
      hall   (:,:) = hall   (:,:)+ qthall (:,:)
!
      RETURN
      CONTAINS
!
              SUBROUTINE Smooth_i (array)
              IMPLICIT NONE
              REAL (rprec), INTENT (IN OUT) :: array (:,:)
        !
              INTEGER (iprec) :: i,j,n, idim, jdim
              REAL (rprec), DIMENSION (SIZE(array,1),SIZE(array,2)) :: work
              idim = SIZE (array, DIM = 1)
              jdim = SIZE (array, DIM = 2)
        !
              DO n = 1, nsmthi
        !
                DO j = 1, jdim
!               DO  i = Get_imin_for_grid(j)+2, idim - 1
                DO  i = imin_j(j)+1, idim - 1
                   work(i, j) = (array(i-1,j) + 4.0 * array(i,j)+array(i+1,j))/6.0
                END DO
                work (imin_j(j), j) = array(imin_j(j), j)
                work (idim, j) = array (idim, j)
                END DO
        !
                DO j = 1, jdim
!               DO i = Get_imin_for_grid(j), idim
                DO i = imin_j(j), idim
                  array (i, j) = work (i, j)
                END DO
                END DO

        !
              END DO
        !
              CALL Wrap_around_ghostcells (array, isize, jsize, n_gc)
        !
              RETURN
              END SUBROUTINE Smooth_i
        !
        !
        !
              SUBROUTINE Smooth_j (array)
              IMPLICIT NONE
              REAL (rprec), INTENT (IN OUT) :: array (:,:)
        !
              INTEGER (iprec) :: i,j,n, idim, jdim
              REAL (rprec) :: work (SIZE(array,1), SIZE (array,2))
        !
              idim = SIZE (array, DIM = 1)
              jdim = SIZE (array, DIM = 2)
        !
              DO n = 1, nsmthj
        !
                DO i = 1, idim
                DO j = j1, j2
                  work(i,j)=(array(i,j-1)+4.0*array(i,j)+array(i,j+1))/6.0
                END DO
                END DO
        !
                CALL Wrap_around_ghostcells (work, isize, jsize, n_gc)
        !
                DO i = 1, idim
                DO j = 1, jdim
                  array (i, j) = work (i, j)
                END DO
                END DO
        !
              END DO
        !
              RETURN
              END SUBROUTINE Smooth_j
!
!
       
      END SUBROUTINE Get_active_cond
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_Coeff 
      USE Rcm_variables
      IMPLICIT NONE
!______________________________________________________________________________
!
! code based on subroutine coeff in spiro.agu83.fort         
! last update 06-25-85 by rws                                
!              02-07-96 frt                                 
!                - min_j replaced imin and min_j+1 replaces i1
!
! This subroutine computes the coefficients of the 
! discretized PDE of MI coupling, except for the  
! inhomogenious term. Formulas are from Jaggi and Wolf 1973
! JGR paper. The discretized PDE for the electrostatic 
! potential looks like:
!
! V(i,j) = C1(i,j)*V(i+1,j) + C2(i,j)*V(i-1,j)
!        + C3(i,j)*V(i,j+1) + C4(i,j)*V(i,j-1) + C5(i,j)
!         
! This subroutine computes the coefficients for grid
! points that are inside the high-latitude boundary
! as defined by MIN_J array; that is, for all points (i,j)
! such that i <= min_j(j). If the boundary is non-integer,
! then some coefficients will have to be modified. c1-c5
! will be recomputed in subroutine CASE3; however, that 
! subroutine requires that arrays A, B, and D be computed
! here and be available in CASE3 for all I values.
!
! STANISLAV: if the boundary coincides with an integer
!            I-value, then this subroutine computes 
!            coefficients for i = imin+1 to i=imax.
!            It does not compute the coefficients at
!            i = imin because we don't use the difference
!            equations there, but rather the Dirichlet
!            boundary condtion (value of V). 
!            On the row just inside of the boundary 
!            (i = imin + 1), we approximate first deri-
!            vatives in I by a 2-point forward difference
!            rather than 3-pt. This leads to a O(dlam)
!            accuracy approximation (compared to 
!            (O(dlam**2)) for other points, but this is
!            due to the conductivities changing sharply
!            (edge of auroral zone!), so the 2-pt diff
!            may simply be not representative of the 
!            derivative in question.
!
!-----------------------------------------------------------
!                                                             
!                                                            
!   this subroutine computes the coefficients c1,c2,c3 & c4.
!   these are coefficients of the elliptic magnetosphere-
!   ionosphere coupling  equation that is solved in potent.
!   computed values  of the coeffecients are stored in array c.  
!                                                             
!   this subroutine called from subroutine comput           
!
!______________________________________________________________________________
!
      INTEGER (iprec) :: i,j,k
      REAL (rprec) :: aa, bb, cc, dd, ee, ff, b_min, bc, hmin, hc, &
                      a (isize,jsize), b (isize,jsize), d (isize,jsize)
!                                                         
      c_pde (:,:,:) = 0.0_rprec
!                                                          
      DO j = 1, jsize, jint
      DO i = 1, isize, iint
         a (i, j) = alpha (i, j) * pedpsi (i, j) / beta  (i, j) 
         b (i, j) = beta  (i, j) * pedlam (i, j) / alpha (i, j) 
         d (i, j) = 2.0_rprec * ( b(i, j) / dlam**2 + a(i, j) / dpsi**2 )
      END DO
      END DO
!                                                        
      loop_30: DO  j = j1, j2, jint
         Loop_20: DO  i = imin_j(j), isize, iint
!
            IF (i < CEILING(bndloc(j)) ) THEN
!
!              Definitely outside modeling region, skip point:
!
               CYCLE Loop_20
!
            ELSE

               IF (i < isize .AND. i > imin_j(j) + 1) THEN
!
!                 Strictly inside the modeling region,
!                 Use central differences for I-derivatives:
!
                  bb = b (i + iint, j) - b (i - iint, j)
                  ee = hall(i+iint,j)-hall(i-iint,j)
                  ee = signbe * ee
!
               ELSE IF (i == isize) THEN
!
!                 On the equatorial boundary,
!                 Use backward 3-pt difference for I-derivatives:
!
                  bb = 3.0 * b (i,j) - 4.0 * b (i-1,j) + b (i - 2, j)
                  ee = 3.0*hall (i,j) - 4.0 * hall (i-1,j) + hall (i-2, j)
                  ee = signbe * ee
!
               ELSE
!
!                 On the second row of modeling region,
!                 Use forward 2-pt differences for I-derivatives:
!
                  b_min = 2.0_rprec * b (i, j) - b (i + 1, j)
                  bc = 0.5 * b (i, j)
                  IF (b_min < bc) b_min = bc
                  bb = b (i + 1, j) - b_min
                  hmin = 2.0_rprec * hall (i, j) - hall (i + 1, j)
                  hc = 0.5 * hall (i, j)
                  IF (ABS (hmin)  < ABS (hc) ) hmin = hc
                  ee = (hall (i + 1, j) - hmin)*signbe

               END IF
!
            END IF 
!                                                         
            cc = hall (i, j + jint) - hall (i, j - jint)
            cc = cc * signbe 
            dd = (bb - cc * dlam / dpsi) * 0.25_rprec
            aa = a (i, j + jint) - a (i, j - jint) 
            ff = (aa + ee * dpsi / dlam) * 0.25_rprec
            c_pde (1, i, j) = (b (i, j) + dd) / (d (i, j) * dlam**2) 
            c_pde (2, i, j) = (b (i, j) - dd) / (d (i, j) * dlam**2) 
            c_pde (3, i, j) = (a (i, j) + ff) / (d (i, j) * dpsi**2) 
            c_pde (4, i, j) = (a (i, j) - ff) / (d (i, j) * dpsi**2) 
!
         END DO loop_20
      END DO loop_30

      DO i = 1, isize
         DO k = 1, ncoeff
            c_pde (k, i, j1 - 2) = c_pde (k, i, j2 - 1) 
            c_pde (k, i, j1 - 1) = c_pde (k, i, j2) 
            c_pde (k, i, j2 + 1) = c_pde (k, i, j1) 
         END DO 
      END DO
!                                                        
!                                                           
      RETURN 
      END SUBROUTINE Comput_Coeff
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_C5_wind
      USE Rcm_variables
      IMPLICIT NONE
!                                                                       
      INTEGER (iprec) :: i, j
      REAL (rprec) :: th, th0, th1, th2, dth, dph, denom, dr0, dr1, dr2, &
                      djre, dp1, dp2, denom1, denom2, djph, dw (isize,jsize)
      REAL (rprec), PARAMETER :: bnorm = 1.0E-6_rprec

!
!
!     1. If iwind=0, c5w(i,j) is set equal to zero and                  
!        subroutine returns to calling program.       
!                                                                       
!
      IF (iwind == 0) THEN
         c5w (:,:) = 0.0_rprec
         RETURN
      ELSE
         call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
             'wind is not implemented yet')
      END IF
!
!                                                                       
!     2. If iwind.ne.0, calculation of c5w(i,j) in volts.               
!        Discretization of the rhs of div*(s*e)=birk-div*jw 
!        yields a term dw(i,j)*v(i,j). Then 
!        c5w(i,j)=-div*jw/dw(i,j).    
!        jw=part of ionospheric current due to action of
!        thermospheric winds.       
!                                                                       
!     2.1 calculation of dw(i,j)                                        
!                                                 
      DO j = j1, j2, jint
      DO i = CEILING(bndloc(j))+1, isize
         th = colat (i, j)
         IF (i == isize) THEN
            th0 = pi_by_two - colat (i, j)
            th1 = pi_by_two - colat (i - 1, j)
            th2 = pi_by_two - colat (i - 2, j)
            dth = - 0.5*(1.0_rprec/COS(th2)**2-4.0/COS(th1)**2+3.0/COS(th0)**2)
         ELSE
            th1 = colat (i - 1, j)
            th2 = colat (i + 1, j)
            dth = 0.5 * (1.0_rprec/SIN(th1)**2 - 1.0_rprec /SIN(th2)**2)
         END IF
         dph = 0.5*(aloct(i,j+1) - aloct(i,j-1))
         IF (dph < 0.0_rprec) dph = dph + pi
!                                                                       
         dw (i, j) = 2.0_rprec /ri*( &
            2.0_rprec*COS(th)*pedlam(i,j) / SIN(th)**2 / dth**2 + &
            SIN(th)**2*pedpsi(i,j) / (2.0_rprec*COS(th) ) / dph**2)
      END DO
      END DO
!                                                                       
!     2.2  calculation of -div*jw. meridional component.                
!        div*jw is multiplied by 1.e-6 to express bir in teslas.        
!
      DO j = j1, j2, jint
      DO i = CEILING(bndloc(j))+1, isize
         IF (i == isize) THEN
!                                                                       
!        2.2.1 meridional part at i=imax. derivative is approximated
!              by a 3-point forward difference formula.
!                                                                       
            th0 = pi_by_two - colat (i, j)
            th1 = pi_by_two - colat (i - 1, j)
            th2 = pi_by_two - colat (i - 2, j)
            dth = - 0.5*(1.0_rprec / COS(th2)**2 - 4.0/COS(th1)**2 + &
                         3.0/COS(th0)**2)
            denom = 2.0_rprec * dth
!
!           comment by yong on 7/26/90. introduce "signbe" for the
!           following lines of rcm where there is a hall:
            dr0 = cos (th0) * bnorm * bir (i, j) * ( &
                  pedlam (i,j) * pwe (i,j) + &
                  signbe * hall (i, j) * hwn (i, j) )
            dr1 = cos (th1) * bnorm * bir (i - 1, j) * ( &
                        pedlam (i-1,j) * pwe (i-1,j) + &
                        signbe * hall (i-1,j)*hwn(i-1,j) )
            dr2 = cos (th2) * bnorm * bir (i - 2, j) * ( &
                       pedlam(i-2,j)*pwe(i-2,j) + &
                       signbe * hall (i-2,j) * hwn (i-2,j))
!                                                                       
            djre = - (dr2 - 4.0 * dr1 + 3. * dr0) / denom
!                                                                       
         ELSE
!                                                                       
!           2.2.2 meridional part at i.lt.imax. derivative is
!                 approximated by central differences.
!                                                                       
            th1    = colat (i - 1, j)
            th2    = colat (i + 1, j)
            dth    = 0.5*(1.0_rprec / sin (th1)**2 - 1.0_rprec / sin (th2)**2)
            denom1 = 2.0_rprec * dth
!                                                                       
            dr2    = bnorm * bir(i-1,j) * sin(th1) * &
                     ( pedlam (i-1,j) * pwe(i-1,j) + &
                       hall(i-1,j) * hwn(i-1,j) * signbe )
!                                                                       
            dr1 = bnorm * bir (i + 1, j) * sin (th2) * &
                  ( pedlam (i + 1, j) * pwe (i + 1, j) + &
                    signbe * hall (i + 1, j) * hwn (i + 1, j) )
!                                                                       
            djre = (dr2 - dr1) / denom1
         END IF
!
!        2.2.3 zonal part.derivative is approximated by
!              central differences.
!                                                                       
         th1 = colat (i, j - 1)
         th2 = colat (i, j + 1)
         dph = 0.5 * (aloct (i, j + 1) - aloct (i, j - 1) )
         IF (dph < 0.0_rprec) dph = dph + pi
         denom2 = 2.0_rprec * dph
         dp2 = SIN(th2)**3 / (2.0_rprec * cos (th2) ) * &
               bnorm * bir (i, j + 1) * signbe * &
               (hall(i,j+1) * hwe(i,j+1) - pedpsi(i,j+1) * pwn(i,j+1))
         dp1 = (SIN(th1))**3 / (2.0_rprec * COS(th1) ) * &
               bnorm * bir (i,j-1) * signbe * &
               (hall(i,j-1) * hwe(i,j-1) - pedpsi(i,j-1) * pwn(i,j-1))
!
!        end of change for inserting "signbe" on 7/26/90
!                                                                       
         djph = (dp2 - dp1) / denom2
         c5w (i, j) = - (djre+djph) / dw (i, j)
      END DO
      END DO
!
      CALL Wrap_around_ghostcells (c5w, isize, jsize, n_gc)
!
      RETURN 
      END SUBROUTINE Comput_C5_wind
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_C5_total
      USE Rcm_variables
      IMPLICIT NONE
!
!
      INTEGER (iprec) :: i,j
      REAL (rprec)    :: d
!                                                                       
      DO j = 1, jsize
      DO i = imin_j(j), isize
         d = 2.0_rprec * ( beta(i,j)  * pedlam(i,j) / (alpha(i,j) * dlam**2)  &
                   + alpha(i,j) * pedpsi(i,j) / (beta(i,j) *  dpsi**2) )
         IF (d <= 1.0e-30_rprec) THEN
            c_pde (5, i, j) = 0.0_rprec
         ELSE 
            c_pde (5, i, j) = alpha(i,j) * beta(i,j) * (ri**2) * birk(i,j) / d + &
                          c5w(i,j)
         END IF
      END DO
      END DO
!
      RETURN 
      END SUBROUTINE Comput_c5_total
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_lowlat_boundary
      USE Rcm_variables
      IMPLICIT NONE
!                                                                       
!______________________________________________________________________________
!                                                                       
!  last update| 08-27-86        written by| g.a.mantjoukis     
!                                                                       
!  subroutine to compute equatorial bndy condition.based on             
!   mantjou.eqbndy.text                                                 
!                                                                       
!  current conservation requires that at i=imax,                        
!  v(imax,j)=c(1,imax,j)*v(imax+1,j)
!           +c(2,imax,j)*v(imax-1,j)            
!           +c(3,imax,j)*v(imax,j+1)
!           +c(4,imax,j)*v(imax,j-1)
!           +c(5,imax,j)
!  where v(imax+1,j) is outside the modeling region.                    
!                                                                       
!  the equatorial bndy condition gives an expression for
!  v(imax+1,j) in terms of quantities inside the modeling
!  region and is of the form    
!  v(imax+1,j)=ceq1*v(imax,j-1)+ceq2*v(imax,j)
!             +ceq3*v(imax,j+1)+ceq4*v(imax-1,j)+ceq5    
!  where ceq1 through ceq5 are calculated below.                        
!                                                                       
!  ss(j) is a cowling-type conductance (see mantjou.eqbndy.text)        
!       integrated over the cross-section of the equatorial band,at     
!       any given local time.                                           
!                                                                       
!  sw(j) is a wind-dependent quantity that contributes to ceq5          
!                                                                       
!                                                                       
!  to set bnd cond to no current across imax 
!  (ie., no eq electrojet) explicityly zero ss(j) for all j  
!                                                                       
!______________________________________________________________________________
!                                                                       
      INTEGER (iprec) :: i, j, n
      REAL (rprec) :: cf, ceq1, ceq2, ceq3, ceq4, ceq5, den
      REAL (rprec), PARAMETER :: bnorm = 1.0E-6_rprec
!
      i = isize
      DO j = 1, jsize
         cf = alpha (i,j) * dlam / beta (i,j) / dpsi / pedlam (i,j)
         ceq1 = cf * (signbe*hall (i, j) -  &
                0.5 * (ss (j + 1) - 4.0 * ss (j) - ss (j - 1) ) / dpsi)
         ceq2 = - 4.0 * cf * ss (j) / dpsi
         ceq3 = cf * ( - signbe * hall (i, j) +  &
                0.5 * (ss (j + 1) + 4.0 * ss (j) - ss (j - 1) ) / dpsi)
         ceq4 = 1.0_rprec
         ceq5 = - 2.0_rprec * ri * alpha (i, j) * dlam * bnorm * bir (i, j) &
                * (pwe(i,j) + signbe * hall(i,j) / pedlam (i,j) * hwn (i,j))
         ceq5 = ceq5 - cf * (sw (j + 1) - sw (j - 1) ) 
         den  = 1.0_rprec - ceq2 * c_pde (1, i, j)
!                                                                       
         c_pde (5, i, j) = (c_pde (5, i, j) + ceq5 * c_pde (1, i, j) ) / den
         c_pde (4, i, j) = (c_pde (4, i, j) + ceq1 * c_pde (1, i, j) ) / den
         c_pde (3, i, j) = (c_pde (3, i, j) + ceq3 * c_pde (1, i, j) ) / den
         c_pde (2, i, j) = (c_pde (2, i, j) + ceq4 * c_pde (1, i, j) ) / den
         c_pde (1, i, j) = 0.0_rprec
      END DO 
!                                                                       
      DO n = 1, ncoeff 
         CALL Wrap_around_ghostcells (c_pde (n,:,:), isize, jsize, n_gc)
      END DO 
!                                                                       
      RETURN 
      END SUBROUTINE Comput_lowlat_boundary
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE New_coeff
      USE Rcm_variables
      IMPLICIT NONE
!
!     ..Local variables:
!
      LOGICAL :: pt_loc (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc)
      LOGICAL :: a2_loc, a3_loc, a4_loc
!
      REAL(rprec) :: alp1, bet1, gamma, lam_f, lam_b, psi_f, psi_b
!
      REAL(rprec) :: a (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                    b (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc)
      REAL(rprec) :: bb, ee, aa, cc, denom
      REAL(rprec) :: g_lam, g_psi, f_lam, f_psi
      REAL(rprec) :: cl1_a1, cl1_a2, cl1_p
      REAL(rprec) :: cl2_a1, cl2_a2, cl2_p
      REAL(rprec) :: cp1_a3, cp1_a4, cp1_p
      REAL(rprec) :: cp2_a3, cp2_a4, cp2_p
      REAL(rprec) :: v_right, v_left, dis_r_1, dis_r_2
      REAL(rprec) :: dis_L_1, dis_L_2, b_min, bc, hmin, hc
!
      INTEGER (iprec), SAVE :: N_calls = 0
      INTEGER (iprec):: i, j, kindex
!
!         j-1    j     j+1        
!
!   i-1..  x     A2     x
!                |
!                |              
!   i ...  A4----P----A3
!                |
!                |
!   i+1..  x     A1     x
!
!     In this subroutine we presume that high-lat. bndy is
!     adequately specified before calling it, that is, 
!     ain(j)+1 > min_j(j) >= ain(j) and |min(j)-ain(j)| > 1E-6
!     and min_j(j) >=2 (check all this outside?)
!
!     .. Executable statements:
!
      N_calls = N_calls + 1
!
      IF (isize < 3) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
          'idim must be larger than 2 in NEW_COEFF')
!
!
!     1. Run over ALL grid points and flag them as being 
!     inside the modeling region (includes integer boundary
!     crossings) with pt_loc = 1 or outside the modeling
!     region with pt_loc = 0:
!
      DO j = 1, jsize
         DO i = 1, isize
            IF (i >= bndloc(j)) THEN
               pt_loc(i,j) = .TRUE.
            ELSE 
               pt_loc(i,j) = .FALSE.
            END IF
         END DO
      END DO
      CALL Wrap_around_ghostcells_log (pt_loc, isize, jsize, n_gc)
! 
!
!-------------------------------------------------------
!     Preliminary calculations:
!
      DO i = 1, isize
         DO j = 1, jsize
!
            a (i,j) = alpha(i,j)*pedpsi(i,j)/beta(i,j)
            b (i,j) = beta(i,j)*pedlam(i,j)/alpha(i,j)
!
         END DO
      END DO
      CALL Wrap_around_ghostcells (a, isize, jsize, n_gc)
      CALL Wrap_around_ghostcells (b, isize, jsize, n_gc)
!------------------------------------------------------- 
!
!     2. Run over grid points inside boundaries. For each
!     point in the modeling region, generate coefficients
!     of the difference equation. Our PDE has the form:
!
!     -b*D2V/Dlambda^2 -a*D2V/Dpsi^2 + 
!     + (-Db/Dlambda + Dhall/Dpsi)*DV/Dlambda +
!     + (-Da/Dpsi - Dhall/Dlambda)*DV/Dpsi = RHS
!     or
!     g_lambda*D2V/Dlambda^2 + g_psi*D2V/Dpsi^2 +
!     f_lambda*DV/Dlambda    + f_psi*DV/Dpsi     = RHS
!
!     We need: to approximate derivatives of V by finite
!     differences, and approximate "coefficients" f_lambda,
!     f_psi, g_lambda, g_psi by finite differences. For 
!     approximating coefficients, we can only use neighbors
!     that are integer grid points since conductances are 
!     evaluated on the grid. For approximating derivatives of
!     V, we will use values of V on non-integer grid neighbors
!     computed from the boundary condition on V.
!
      Loop_j: DO j = 1, jsize
!        Loop_i:  DO i = imin_j(j), isize
         Loop_i:  DO i = 1, isize
!
            IF (i < CEILING(bndloc(j))) CYCLE
!
!            IF (.NOT.pt_loc(i,j)) CYCLE
!
!           For each grid point (i,j)=P in the modeling
!           region, we will need four neighbors.
!           Determine how many of the grid neighbors are
!           inside the modeling region, and flag them. 
!
            a2_loc = pt_loc(MAX (i - 1, 1),j)
            a3_loc = pt_loc(i,j + 1)
            a4_loc = pt_loc(i,j - 1)
!
!
!           Determine distances from (i,j) to its 3 neighbors
!           (A_1 is always one grid unit away), in grid units:
!
            IF (a2_loc) THEN      ! A_2 is also grid point
               alp1 = 1.0_rprec
            ELSE                  ! A_2 is a boundary point
               alp1 = REAL(i,rprec) - bndloc(j)
            END IF
!
            IF (alp1 < 1E-6) THEN 
!
!               This is a special case when the point (i,j) is
!               on the boundary (integer boundary crossing).
!               Handle as Dirichlet boundary condition:
!
                c_pde(1:4,i,j) = 0.0_rprec
                c_pde(5,i,j)   = vbnd(j)
                CYCLE Loop_i
            END IF
!
            IF (a3_loc) THEN      ! A_3 is a grid point
               bet1 = 1.0_rprec
            ELSE                  ! A_3 is on the boundary
               bet1    = (i-bndloc(j)) / (bndloc(j+1)-bndloc(j))
               dis_r_1 = SQRT( (i-bndloc(j))**2 + bet1**2)
               dis_r_2 = SQRT( (bndloc(j+1)-i)**2 + (1.0_rprec-bet1)**2)
               v_right = (dis_r_1*vbnd(j+1)+dis_r_2*vbnd(j)) &
                         / (dis_r_1+dis_r_2)
            END IF
!
            IF (a4_loc) THEN      ! A_4 is a grid point
               gamma = 1.0_rprec
            ELSE                  ! A_4 is on the boundary
               gamma = (i-bndloc(j))/(bndloc(j-1)-bndloc(j))
               dis_L_1 = SQRT( (i-bndloc(j))**2 + gamma**2)
               dis_L_2 = SQRT( (bndloc(j-1)-i)**2 + (1.0_rprec-gamma)**2)
               v_left  = (dis_L_1*vbnd(j-1)+dis_L_2*vbnd(j)) &
                         / (dis_L_1+dis_L_2)
            END IF
!
!
!           Approximate coefficients with lambda-derivatives:
!
            IF (ABS(alp1-1.0_rprec) < 10.0*EPSILON(1.0_rprec)) THEN

               IF (i < isize) THEN
!           
!                 (i,j) is an interior grid point, and we
!                 can use central differences formula for
!                 lambda-derivatives:
! 
                  bb = b(i+1,j) - b(i-1,j)
                  ee = (hall(i+1,j) - hall(i-1,j))*signbe

               ELSE
!
!                 (i,j) in on low-latitude boundary, need
!                 to use backward differences for deriv:
!
                  bb = 3.0 * b(i,j) - 4.0 * b(i-1,j)+b(i-2,j)
                  ee = 3.0 * hall(i,j)-4.0 * hall(i-1,j)+hall(i-2,j)
                  ee = ee * signbe

               END IF
!
            ELSE
!         
!              alp1 < 1, so "i-1,j" grid point is outside
!              the boundary, and we need forward difference:
!
!!!!           bb = b(i+1,j) - &
!!!!                MAX (0.5*b(i,j),2.*b(i,j)-b(i+1,j))
!!!!           hmin = 2.*hall(i,j) - hall(i+1,j)
!!!!           hc   = 0.5*hall(i,j)
!!!!           IF (ABS(hmin) < ABS(hc)) hmin = hc
!!!!!!         ee   = signbe*(hall(i+1,j)-hmin)
!>                bb = -3.0*b(i,j)+ 4.0*b(i+1,j)-b(i+2,j)
!>                ee = -3.0*hall(i,j)+4.0*hall(i+1,j)-hall(i+2,j)
!>                ee = ee * signbe
              b_min = 2. * b (i, j) - b (i + 1, j)
              bc = 0.5 * b (i, j)
              IF (b_min < bc) b_min = bc
              bb = b (i + 1, j) - b_min
              hmin = 2. * hall (i, j) - hall (i + 1, j)
              hc = 0.5 * hall (i, j)
              IF (ABS (hmin)  < ABS (hc) ) hmin = hc
              ee = (hall (i + 1, j) - hmin)*signbe

!
            END IF
!
!
!           Approximate coefficients with psi-derivatives:
!
            IF (ABS(bet1-1.0_rprec) < EPSILON(1.0_rprec) .AND. &
                ABS(gamma-1.0_rprec) < EPSILON(1.0_rprec))            THEN
!
!               (i,j) is an inner point, can use central
!               differences:
!
                cc = signbe*(hall(i,j+1)-hall(i,j-1))
                aa = a(i,j+1) - a(i,j-1)
!
            ELSE IF (bet1 < 1.0_rprec .AND. ABS(gamma-1.0_rprec) < EPSILON(1.0_rprec)) THEN
!
!               use backward difference, mult. by 2:
!
                cc = 2.0_rprec*signbe*(hall(i,j)-hall(i,j-1))
                aa = 2.0_rprec*(a(i,j)-a(i,j-1))
!
            ELSE IF (ABS(bet1-1.0_rprec) < EPSILON(1.0_rprec) .AND. gamma < 1.0_rprec) THEN
!
!               use forward difference, mult. by 2:
!
                cc = 2.0_rprec * signbe*(hall(i,j+1)-hall(i,j))
                aa = 2.0_rprec * (a(i,j+1)-a(i,j))
!
            ELSE
!
!               gamma and bet1 are < 1, set derivs to zero:
!
                aa = 0.0_rprec
                cc = 0.0_rprec
!
            END IF
!               
!
            g_lam       = - b(i,j)
            g_psi       = - a(i,j)
            f_lam       = - bb/dlam/2. + cc/dpsi/2. 
            f_psi       = - ee/dlam/2. - aa/dpsi/2. 
!
!
!           Approximate partial derivatives of V in the PDE
!           DV/Dlambda, DV/Dpsi, D2V/Dlambda^2, D2V/Dpsi^2. Lambda
!           derivatives will be linear combinations of V(P), V(A_1),
!           and V(A_2); psi derivatives will be linear combinations
!           of V(P), V(A_3), and V(A_4); here P=(i,j). We use
!           notation
!           DV/Dlambda =    cl1_a1*V(A_1)+cl1_a2*V(A_2)+cl1_p*V(P)
!           D2V/Dlambda^2 = cl2_a1*V(A_1)+cl2_a2*V(A_2)+cl2_p*V(P)
!           DV/Dpsi =       cp1_a3*V(A_3)+cp1_a4*V(A_4)+cp1_p*V(P)
!           D2V/Dpsi^2 =    cp2_a3*V(A_3)+cp2_a4*V(A_4)+cp2_p*V(P)
!
!
!           Compute the distances to the 4 neighbors:
!
            lam_f = dlam
            lam_b = dlam * alp1
            psi_f = dpsi * bet1
            psi_b = dpsi * gamma
!
            cl1_a1 = + lam_b / lam_f / (lam_f+lam_b)
            cl1_a2 = - lam_f / lam_b / (lam_f+lam_b)
            cl1_p  = + (lam_f-lam_b) / (lam_f*lam_b)
!
            cl2_a1 = + 2. / lam_f / (lam_f+lam_b)
            cl2_a2 = + 2. / lam_b / (lam_f+lam_b)
            cl2_p  = - 2. / (lam_f*lam_b)
!
            cp1_a3 = + psi_b / psi_f / (psi_f+psi_b)
            cp1_a4 = - psi_f / psi_b / (psi_f+psi_b)
            cp1_p  = + (psi_f-psi_b) / (psi_f*psi_b)
!
            cp2_a3 = + 2. / psi_f / (psi_f+psi_b)
            cp2_a4 = + 2. / psi_b / (psi_f+psi_b)
            cp2_p  = - 2. / (psi_f*psi_b)
!
            denom     = g_lam*cl2_p  + g_psi*cp2_p &
                      + f_lam*cl1_p  + f_psi*cp1_p
!
            c_pde (1,i,j) = g_lam*cl2_a1 + f_lam*cl1_a1
            c_pde (2,i,j) = g_lam*cl2_a2 + f_lam*cl1_a2
            c_pde (3,i,j) = g_psi*cp2_a3 + f_psi*cp1_a3
            c_pde (4,i,j) = g_psi*cp2_a4 + f_psi*cp1_a4
!
            c_pde (1,i,j) = - c_pde(1,i,j)/denom
            c_pde (2,i,j) = - c_pde(2,i,j)/denom
            c_pde (3,i,j) = - c_pde(3,i,j)/denom
            c_pde (4,i,j) = - c_pde(4,i,j)/denom
            c_pde (5,i,j) = + c_pde(5,i,j)/denom
!
            IF (.NOT.a2_loc) THEN
               c_pde (5,i,j) = c_pde(5,i,j) + c_pde(2,i,j)*vbnd(j)
               c_pde (2,i,j) = 0.
            END IF
!
            IF (.NOT.a3_loc) THEN
               c_pde (5,i,j) = c_pde(5,i,j) + c_pde(3,i,j)*v_right
               c_pde (3,i,j) = 0.
            END IF
!
            IF (.NOT.a4_loc) THEN
               c_pde (5,i,j) = c_pde(5,i,j) + c_pde(4,i,j)*v_left
               c_pde (4,i,j) = 0.
            END IF
!
         END DO Loop_i
!        don't put anything between these two ENDDOs
      END DO Loop_j
!
!
      DO kindex = 1, ncoeff
         CALL Wrap_around_ghostcells (c_pde(kindex,:,:), isize, jsize, n_gc)
      END DO
!
      RETURN
      END SUBROUTINE New_coeff
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE New_cfive 
      USE Rcm_variables
      IMPLICIT NONE
!
      INTEGER (iprec) :: i,j
!
      DO j = 1, jsize
      DO i = 1, isize
!     DO i = imin_j(j), isize
!        c_pde (5,i,j) = alpha(i,j)*beta(i,j)*(Ri**2)*birk(i,j) + c5w(i,j)
         c_pde (5,i,j) = alpha(i,j)*beta(i,j)*(Ri**2)*birk(i,j) + c5w(i,j)
      END DO
      END DO
      CALL Wrap_around_ghostcells (c_pde(5,:,:), isize, jsize, n_gc)
!
      RETURN 
      END SUBROUTINE New_cfive                          
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_V_Potnt3 (imax,jmax,n_coef, imin_j, c, u)
      USE Rcm_variables, ONLY : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: imax,jmax,n_coef, imin_j(jmax)
      REAL (rprec) , INTENT (IN)   :: c (imax,jmax,n_coef)
      REAL (rprec), INTENT(IN OUT) :: u (imax,jmax)
!
      INTEGER (iprec), PARAMETER :: maxits = 50000
      REAL (rprec) :: cut,rjac, anorm, resid, omega, pi
      INTEGER  (iprec) :: n,i,j
!
!
      pi = 4.0_rprec * ATAN (1.0_rprec)
      cut = 10.0_rprec        ! cut for sum of residuals in volts
      rjac = 1.0_rprec - pi**2 /1.0_rprec / (jmax**2) ! what is it for a non-square grid?
!                                                                       
      omega = 1.0_rprec
!                                                                       
!
!
      iterate_loop: DO n = 1, maxits 
!
         anorm = 0.0_rprec
!                                                                       
!        Inner boundary using coefficients   
!
         DO j = 2, jmax - 1 
            i = imin_j (j) 
            u (i, j) = c (5, i, j) + c (1, i, j) * u (i + 1, j) &
                                   + c (2, i, j) * u (i - 1, j) &
                                   + c (3, i, j) * u (i, j + 1) &
                                   + c (4, i, j) * u (i, j - 1)
         END DO 
!                                                                       
!        Outer boundary using coefficients             
!        took out c1 because imax+1 is out of bounds  
!
         DO j = 2, jmax - 1 
            i = imax 
            u (i, j) = c (5, i, j) + c (2, i, j) * u (i - 1, j) &
                                   + c (3, i, j) * u (i, j + 1) &
                                   + c (4, i, j) * u (i, j - 1) 
         END DO 
!
!        Use periodicity to get other points:            
!                                                                       
         u (imin_j (1), 1) = u (imin_j (jmax - 2), jmax - 2) 
         u (imin_j (jmax), jmax) = u (imin_j (3), 3) 
!                                                                       
         u (imax, 1) = u (imax, jmax - 2) 
         u (imax, jmax) = u (imax, 3) 
!                                                                       
         DO j = 2, jmax-2
         DO i = imin_j(j)+1, imax - 1 
!
           u (i, 1)        = u (i,jmax-2)
           u (i, jmax - 1) = u (i, 2)
           u (i,jmax)      = u (i,3)
!
           IF (MOD (i + j, 2) == MOD (n, 2) ) THEN 
              resid = - c (1, i, j) * u (i + 1, j) &
                      - c (2, i, j) * u (i - 1, j) &
                      - c (3, i, j) * u (i, j + 1) &
                      - c (4, i, j) * u (i, j - 1) &
                      + u (i, j) - c (5, i, j) 
!
              anorm = anorm + ABS (resid) 
              u (i, j) = u (i, j) - omega * resid 
!
           END IF 
!
         END DO 
         END DO 
!
WRITE (*,*) 'anorm',n,anorm
         IF (n == 1) THEN 
            omega = 1.0_rprec / (1.0_rprec-0.5_rprec * rjac**2) 
         ELSE 
            omega = 1.0_rprec / (1.0_rprec-0.25_rprec * rjac**2 * omega) 
         END IF 
!                                                                       
         IF (n >  1 .AND. anorm < cut) THEN
             WRITE (*,*) 'in potnt3',n,anorm
             EXIT iterate_loop
         END IF
!                                                                       
      END DO iterate_loop
!
      IF (n >= maxits) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'//&
          'maxits exceeded')
!
      RETURN
      END SUBROUTINE Comput_V_Potnt3
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Thet (aa, bb, xc, yc, phi)
      USE Rcm_variables
!
!     this function gives solution to equation for ellipse in
!     in flat polar coordinates.  it specifies colatitude as a
!     function of local time.
!
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: aa, bb, xc, yc, phi
      REAL (rprec) :: thet
!
      REAL (rprec) :: cphi, sphi, ca, cb, cc
!
!
         cphi = COS(phi)
         sphi = SIN(phi)
         ca = (cphi/aa)**2 + (sphi/bb)**2
         cb = -xc*cphi/aa**2 - yc*sphi/bb**2
         cc = (xc/aa)**2 + (yc/bb)**2 - 1.
         thet = (-cb+SQRT(cb**2-ca*cc))/ca
      RETURN
      END FUNCTION Thet
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      SUBROUTINE Elemod (icase, ikp, glat, amlt, value)
      USE Rcm_variables, ONLY : iprec, rprec, lun, NameRcmDir
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: icase, ikp
      REAL (rprec),    INTENT (IN) :: glat, amlt
      REAL (rprec),    INTENT (OUT):: value
!
!______________________________________________________________________________
!
!   Stanislav Tue May 25 20:40:19 MDT 1999: adopt the hardy
!             code to use with RCM. This subroutine communicates
!             with RCM through the arguments, but also needs
!             module rcm_mod_hardy to store coeffs.
!
!   SUBROUTINE TO EVALUATE THE HARDY AVERAGE AURORAL MODEL
!     DESCRIBED IN:  HARDY, ET AL., "J. GEOPHYS. RES.",
!     VOL. 92, PAGE 12,275-12,294 (1987).
!
!
!   INPUTS:
!
!     ICASE=1   ENERGY FLUX
!           2   NUMBER FLUX
!           3   HALL CONDUCTIVITY
!           4   PEDERSON CONDUCTIVITY
!
!     IKP       KP DIVISION 0 TO 6
!
!     GLAT      GEOMAGNETIC LATITUDE
!
!     AMLT      MAGNETIC LOCAL TIME
!
!   OUTPUTS:
!
!     VALUE     LOG10 ENERGY FLUX IN KEV/CM**2-S-SR (ICASE=1)
!               LOG10 NUMBER FLUX IN ELECTRONS/CM**2-S-SR (ICASE=2)
!               CONDUCTIVITY IN MHOS (ICASE=3,4)
!  
!   INTERNAL VARIABLES
!
!     CRD       COEFFICIENTS FOR MAXIMUM FUNCTION VALUE
!     CHAT      COEFFICIENTS FOR LATITUDE OF MAXIMUM VALUE
!     CS1       COEFFICIENTS FOR UP-SLOPE
!     CS2       COEFFICIENTS FOR DOWN-SLOPE
!     CUTL      LOW LATITUDE CUT-OFF VALUE
!     CUTH      HIGH LATITUDE CUT-OFF VALUE
!
!   FILES:
!
!     THE FILE ELECOEF.DAT MUST BE PRESENT IN THE DEFAULT
!     DIRECTORY.
!
!   NOTES:
!
!     THIS VERSION OPERATES ON VAX/VMS OR IBM-PC MACHINES.
!
!   21 JUNE 1993 -- WJM
!______________________________________________________________________________
!
      REAL (rprec), DIMENSION(4), PARAMETER :: &
                    cutl=(/6.,6.,0.,0./), cuth=(/7.,7.,.55,.55/)
!      CHARACTER(LEN=80) :: aline
      INTEGER (iprec) :: j, jcase, jco, jkp, kp, ipc, ips
      REAL (rprec) :: xarg, rd, hat, s1, s2, xa, c, s
      INTEGER (iprec), SAVE :: iread_hardy_first = 0
      REAL (rprec), SAVE :: crd(13,7,4), chat(13,7,4), cs1(13,7,4), cs2(13,7,4)
!
      value = 0.0
!
      IF (ikp > 6) call CON_stop('ERROR in IM/RCM2/src/rcm_comput.f90:'// &
          'ELEMOD: KP > 6 !')  ! stanislav
!
!
      IF (iread_hardy_first == 0) THEN
        iread_hardy_first = 1
!
        OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'input/elecoef.dat', STATUS = 'OLD')
        READ (LUN,'(A80)')  ! aline
        DO jcase = 1, 4
           DO jkp = 1, 7
              DO jco = 1, 13
                 READ (LUN, '(26X,4F12.7)') crd (jco,jkp,jcase), &
                                            chat (jco,jkp,jcase),&
                                            cs1 (jco,jkp,jcase), &
                                            cs2 (jco,jkp,jcase)
              END DO
           END DO
        END DO
        CLOSE (LUN)
      END IF
!
!
      IF (glat < 50.0) RETURN
!
      kp   = ikp + 1
      xarg = amlt*3.14159265/12.
      rd   = crd(1,kp,icase)
      hat  = chat(1,kp,icase)
      s1   = cs1(1,kp,icase)
      s2   = cs2(1,kp,icase)
!
      DO j = 1, 6
         xa  = j*xarg
         c   = COS(xa)
         s   = SIN(xa)
         ipc = j+1
         ips = j+7
!
         rd = rd+c*crd(ipc,kp,icase)+s*crd(ips,kp,icase)
         hat= hat+c*chat(ipc,kp,icase)+s*chat(ips,kp,icase)
         s1 = s1+c*cs1(ipc,kp,icase)+s*cs1(ips,kp,icase)
         s2 = s2+c*cs2(ipc,kp,icase)+s*cs2(ips,kp,icase)
!
      END DO
!
      value = Epst (glat,rd,hat,s1,s2,cutl(icase),cuth(icase))
!
      RETURN
!
      CONTAINS
!
          FUNCTION Epst (clat, rd, hat, s1, s2, xmin, xmax)
          IMPLICIT NONE
          REAL (rprec), INTENT (IN) :: clat, rd, hat, s1, s2, xmin, xmax
          REAL (rprec):: Epst
    !
    !
          REAL (rprec) :: d, ex, xl, ep
    !
          d  = clat-hat
          ex = EXP (d)
          xl = (1.0 - s1 / s2 * ex) / (1.0 - s1 / s2)
          xl = LOG (xl)
          ep = rd+s1*d+(s2-s1)*xl
    !
          IF (clat < hat .AND. ep < xmin) ep=xmin
          IF (clat > hat .AND. ep < xmax) ep=xmax
    !
          epst = ep
          RETURN
          END FUNCTION epst
!
      END SUBROUTINE Elemod
