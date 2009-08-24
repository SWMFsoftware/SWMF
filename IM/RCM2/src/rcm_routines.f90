    SUBROUTINE Compute_vel_edges (dt)
    USE Rcm_variables
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
!_____________________________________________________________________________
!
!   Last update: 07-22-85                   by:rws
!                01-29-96 replace imin with min_j - frt
!                10-05-98 delete alast stuff - stanislav
!                                                                       
!   This subroutine computes dbidt and dbjdt, the
!   velocities of the string of bi,bj pts due to
!   simple drift motion.
!                                                                       
!   Dependency:  GNTRP
!_____________________________________________________________________________
!
    INTEGER (iprec) :: k, m, mbeg, mend, kdim, jdim, idim
    REAL (rprec) :: veff (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                    dvefdi (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                    dvefdj (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc),&
                    bik, bjk, bik1, bjk1, bik2, bjk2, bik3,&
                    bjk3, bik4, bjk4
!
    LOGICAL :: DoWarn = .true.
!
!
    kdim = SIZE (alam)
    jdim = SIZE (v, DIM = 2)
    idim = SIZE (v, DIM = 1)

    DO k = 1, kdim

!      1. Compute effective potential on the grid and bndy:

       veff   = v + vcorot - vpar + alam(k) * vm
       CALL Wrap_around_ghostcells (veff, isize, jsize, n_gc)
!
!      2. Differentiate it with respect to I and J:

       CALL Deriv_i (veff, isize, jsize, imin_j, dvefdi)
       CALL Deriv_j (veff, isize, jsize, imin_j, 1.0e+25, dvefdj)
       WHERE (ABS(dvefdj) > 1.0E+24)
          dvefdi = 0.0
          dvefdj = 0.0
       END WHERE
!
!
!      3. Estimate velocities at the locations of particles:
!
       mbeg = mpoint (k)
       mend = mbeg + npoint (k) - 1
!                                                                       
       DO m = mbeg, mend
!
          IF (ivoptn == 1) THEN
!
!           Interpolate dvefdi,dvefdj,bir,alpha,and beta for location
!           bi(m),bj(m):
!
            bik = bi (m)
            bjk = bj (m)
            CALL Check_ij_limits  (bik, bjk)
            dbidt (m) = +Fcn1 (dvefdj, bik, bjk )
            dbjdt (m) = -Fcn1 (dvefdi, bik, bjk )
!
          ELSE IF (ivoptn == 2) THEN
!
!           Compute dbidt and dbjdt using 4th order Runge-Kutta scheme:
!
            bik = bi (m)
            bjk = bj (m)
            CALL Check_ij_limits (bik, bjk)
            bik1 = + dt * Fcn1 (dvefdj, bik, bjk)
            bjk1 = - dt * Fcn1 (dvefdi, bik, bjk)
!
            bjk = bj (m) + 0.5_rprec * bjk1
            bik = bi (m) + 0.5_rprec * bik1
            CALL Check_ij_limits (bik, bjk)
            bik2 = + dt * Fcn1 (dvefdj, bik, bjk)
            bjk2 = - dt * Fcn1 (dvefdi, bik, bjk)
!
            bjk = bj(m) + 0.5_rprec * bjk2
            bik = bi(m) + 0.5_rprec * bik2
            CALL Check_ij_limits (bik, bjk)
            bik3 = + dt * Fcn1 (dvefdj, bik, bjk)
            bjk3 = - dt * Fcn1 (dvefdi, bik, bjk)
!
            bjk = bj(m) + bjk3
            bik = bi(m) + bik3
            CALL Check_ij_limits (bik, bjk)
            bik4 = + dt * Fcn1 (dvefdj, bik, bjk)
            bjk4 = - dt * Fcn1 (dvefdi, bik, bjk)
!
            dbidt (m) = (bik1 + 2.0_rprec*bik2 + 2.0_rprec*bik3 + bik4) / &
                        (6.0_rprec * dt)
            dbjdt (m) = (bjk1 + 2.0_rprec*bjk2 + 2.0_rprec*bjk3 + bjk4) / &
                        (6.0_rprec * dt)
!
          ELSE
            call CON_STOP('IVOPTN value is illegal in VCALC')
          END IF
!
       END DO
    END DO
    RETURN
!
    CONTAINS
       FUNCTION Fcn1 (array, bi_local, bj_local)
       IMPLICIT NONE
       REAL (rprec)  :: fcn1
       REAL (rprec), INTENT (IN) :: array (:,:), bi_local, bj_local
!
       REAL (rprec) :: az, bir2, alph2, beta2, Gntrp_2d
!
       az    = Gntrp_2d (array, isize, jsize, bi_local, bj_local)
       bir2  = Gntrp_2d (bir,   isize, jsize, bi_local, bj_local)
       alph2 = Gntrp_2d (alpha, isize, jsize, bi_local, bj_local)
       beta2 = Gntrp_2d (beta,  isize, jsize, bi_local, bj_local)
       fcn1  = signbe * az *1.0E3_rprec / &
            (bir2 * alph2 * beta2 * dlam * dpsi * ri**2)
       RETURN
       END FUNCTION Fcn1

       SUBROUTINE Check_ij_limits (bi, bj)
       IMPLICIT NONE
       REAL (KIND=rprec), INTENT (IN) :: bi, bj

       IF (bi > REAL(idim,rprec) .OR. bi < 1.0_rprec) THEN
           call CON_STOP('BI is out of range, stopping')
       ELSE IF ((bj > REAL(jdim,rprec) .OR. bj < 0.99_rprec) .AND. DoWarn) THEN
           write(*,*)'IM_WARNING Check_ij_limits: BJ is out of range, ',&
                'bi,bj=',bi,bj,' !!!'
           DoWarn = .false.
           !call CON_STOP('BJ is out of range, stopping')
       END IF
       RETURN
       END SUBROUTINE Check_ij_limits

    END SUBROUTINE Compute_vel_edges
!
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!9. Facilities for moving particles (both edges and grid) in one time step.
!
!
    SUBROUTINE Move_plasma ( dt )
      USE Rcm_variables
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
!_____________________________________________________________________________
!
!  Time step subroutine to do simple euler time step                    
!                                                                       
!  Last update:
!   8-29-86                                                 
!   1-29-96 frt added boundary arrays and calls to bndy     
!   3-19-97 rws ibtime and nbf added as calling parameters  
!   10-02-98 sts fudge is sized as kcdim for electrons on grid
!   may 99 sts removed hardy coeffs--they are in module
!_____________________________________________________________________________
!
!
!sys CALL Move_plasma_edges (dt)
    IF (i_advect == 1) THEN
       CALL Move_plasma_grid  (dt, 1)
       CALL Move_plasma_grid  (dt, 2)
       CALL Move_plasma_grid  (dt, 3)
    ELSE IF (i_advect == 2) THEN
       CALL Move_plasma_grid_new (n_gc, isize,jsize,kcsize, iesize, i1, i2, &
                                  j1, j2, imin_j, dt, vcorot, vpar, v, alamc,&
                                  etac, eeta, vm, xmass, fudgec, ikflavc, sini,&
                                  alpha,beta, bir, &
                                  dlam, dpsi, signbe, Ri, 3)
       CALL Move_plasma_grid_new (n_gc, isize,jsize,kcsize, iesize, i1, i2, &
                                  j1, j2, imin_j, dt, vcorot, vpar, v, alamc,&
                                  etac, eeta, vm, xmass, fudgec, ikflavc, sini,&
                                  alpha,beta, bir, &
                                  dlam, dpsi, signbe, Ri, 2)
       CALL Move_plasma_grid (dt, 1)
    ELSE IF (i_advect == 3) THEN
       CALL Move_plasma_grid_new (n_gc,isize,jsize,kcsize,iesize, i1, i2, &
                                  j1, j2, imin_j, dt, vcorot, vpar, v, alamc,&
                                  etac, eeta, vm, xmass, fudgec, ikflavc, sini, &
                                  alpha,beta, bir, &
                                  dlam, dpsi, signbe, Ri, 3)
       CALL Move_plasma_grid_new (n_gc,isize,jsize,kcsize,iesize, i1, i2, &
                                  j1, j2, imin_j, dt, vcorot, vpar, v, alamc,&
                                  etac, eeta, vm, xmass, fudgec, ikflavc, sini, &
                                  alpha,beta, bir, &
                                  dlam, dpsi, signbe, Ri, 2)
       CALL Move_plasma_grid_new (n_gc,isize,jsize,kcsize,iesize, i1, i2, &
                                  j1, j2, imin_j, dt, vcorot, vpar, v, alamc,&
                                  etac, eeta, vm, xmass, fudgec, ikflavc, sini, &
                                  alpha,beta, bir, &
                                  dlam, dpsi, signbe, Ri, 1)
    ELSE
       call CON_STOP('ILLEGAL I_ADVECT IN MOVING PLASMA')
    END IF
!
    RETURN
    END SUBROUTINE Move_plasma
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Move_plasma_edges (dt)
      USE Rcm_variables
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
!_____________________________________________________________________________
!   Move plasma edges one time step. Algorithm is as follows:
!   Loop over all edges. For each edges:
!     1. Compute particle velocities for the current time step--VCALC
!     2. Move ALL particles (even those outside boundary)
!        except for held ones
!     3. Correct motion of particles outside boundary--EDGFIX
!   CALLS:       VCALC, DCODE, EDGFIX
!   CALLED FROM: TSTEP1
!_____________________________________________________________________________
!
!
    INTEGER (iprec) :: k, mbeg, mend, m, is
    INTEGER (iprec) :: jjm, jjp, mlast, mnext, Bjmod_int
    REAL    (rprec) :: db_dt_n, Bndy, Bjmod_real
!
!
    CALL Compute_vel_edges (dt)
!
!
    DO k = 1, ksize
       mbeg = mpoint (k)
       mend = mbeg + npoint (k) - 1
       DO m = mbeg, mend

          CALL Dcode (itrack, SIZE(itrack), m, is)
!         IF (is < 0) CYCLE
          IF (is < 0) call CON_STOP('DECODE')

          mnext = m + 1
          IF (mnext > mend) mnext = mbeg
          mlast = m - 1
          IF (mlast < mbeg) mlast = mend
!
          IF ( bi(m) > Bndy(bndloc, jsize, bj(m))) THEN
!
             bi (m) = MIN ( bi(m) + dbidt(m) * dt, REAL(isize,rprec) )
             bj (m) = Bjmod_real ( bj(m) + dbjdt (m) * dt, jsize )
!
          ELSE
!
             jjm = Bjmod_int (INT (bj(m)),  jsize)
             jjp = jjm + 1
             db_dt_n = dbidt(m) - dbjdt(m)*(bndloc(jjp)-bndloc(jjm))
!
             IF (db_dt_n > 0.0_rprec) THEN
!
                bi (m) = MIN ( bi(m) + dbidt(m) * dt, REAL(isize,rprec) )
                bj (m) = Bjmod_real ( bj(m) + dbjdt (m) * dt, jsize )
!
             ELSE IF (bi(mlast) < Bndy(bndloc, jsize, bj(mlast)) .AND. &
                      bi (m)    < Bndy(bndloc, jsize, bj(m))     .AND. &
                      bi(mnext) < Bndy(bndloc, jsize, bj(mnext)) )       THEN
                bi (m) = Bndy(bndloc,jsize, bj(m)) - epslon_edge
!
             ELSE
!
!               DO NOTHING

             END IF

          END IF
!
!!!          bi (m) = MIN ( bi(m) + dbidt(m) * dt, REAL(isize,rprec) )
!!!          bj (m) = Bjmod_real ( bj(m) + dbjdt (m) * dt, jsize )
!

       END DO
       CALL Edgfix (k, dt, epslon_edge )
    END DO
!
    RETURN
    END SUBROUTINE Move_plasma_edges
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Edgfix (k, delta_t, epslon)
      USE Rcm_variables
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: delta_t, epslon
    INTEGER (KIND=iprec), INTENT (IN) :: k
!
!
!   Subroutine to correct motion of inner edge pts at the bndy of
!   modeling region.  Based on raw mac doc dated nov 5,1986
!   This subroutine called by a tstep routine.
!   Note: this routine needs to be rethought if tstep1 is not used
!   as time stepper
!
!   Last update:
!     11-05-86
!     01-29-96 frt added calls  to bndy and array ain(jdim) fr
!     may 99, sts added call to dcode
!     09/27/99 sts: this subroutine ASSUMES that ALL THE POINTS have been
!     moved in TSTEP1, and corrects that motion for some points. You must
!     first move all the points (except held points, of course).
!
!   Dependency: BNDY, DCODE, ADJUST_BJ
!
    INTEGER (iprec) :: mbeg, mend, m, mnext, mlast, is
    REAL (rprec) :: bj_mnext, bj_mlast, bj_m, bj_new_1, bj_new_2, &
                    bnd_pnt, bnd_next, bnd_last, bi_pnt, bi_next, bi_last,&
                    bj_pnt, bj_next, bj_last, Bndy
!
!
    IF (k > SIZE (alam )) call CON_STOP('K OUT OF RANGE IN EDGFIX')
!
    mbeg = mpoint (k)
    mend = mbeg + npoint (k) - 1
!
    Loop_k: DO m = mbeg, mend
!
!      We correct motion of points that are not held and are outside or on
!      the boundary; all others are ignored here:
!
       CALL Dcode (itrack, SIZE(itrack), m, is)
       IF (is < 0) CYCLE Loop_k

       IF (bi(m) < Bndy (bndloc, jsize, bj(m))) THEN
          bi (m) = Bndy (bndloc, jsize, bj(m)) - epslon/delta_t*delta_t
       END IF

!!!       bnd_pnt = Bndy (bndloc, jsize, bj(m) )
!!!       bi_pnt  = bi(m)
!!!       bj_pnt  = bj(m)
!!!       IF (bi_pnt > bnd_pnt) CYCLE loop_k
!
!!!       mnext = m + 1
!!!       IF (mnext > mend) mnext = mbeg
!!!       bnd_next = Bndy (bndloc, jsize, bj(mnext))
!!!       bi_next  = bi(mnext)
!!!       bj_next  = bj(mnext)
!
!!!       mlast = m - 1
!!!       IF (mlast < mbeg) mlast = mend
!!!       bnd_last = Bndy (bndloc, jsize, bj(mlast))
!!!       bi_last  = bi(mlast)
!!!       bj_last  = bj(mlast)
!
!!!       IF (bi_last <= bnd_last .AND. bi_next <= bnd_next ) THEN
!
!         This point should not have been moved, so let's put it back
!         to where it was before moving:
!
!!!          bi(m) = bi_pnt - dbidt(m)*delta_t
!!!          bj(m) = bj_pnt - dbjdt(m)*delta_t
!
!!!       ELSE IF (bi_last <= Bnd_last .AND. bi_next > bnd_next) THEN
!
!         Move the point back along the segment toward "mnext" and put it
!         just outside the boundary:
!
!!!          IF (ABS(bi_next-bi_pnt) < 1.E-5) THEN
!!!              bi_next = bi_next + 1.E-5
!!!          END IF
!!!          CALL Adjust_bj_2 (bj_pnt, bj_next, bj_m, bj_mnext)
!!!          bj(m) = bj_pnt - (bi_pnt-Bnd_pnt+epslon) * &
!!!                           (bj_mnext-bj_m) / (bi_next-bi_pnt)
!!!          bi(m) = Bnd_pnt - epslon
!
!!!       ELSE IF (bi_last >  Bnd_last .AND. bi_next <= bnd_next) THEN
!
!         Move the point back along the segment toward "last" and put it
!         just outside the boundary:
!
!!!          CALL Adjust_bj_2 (bj_pnt, bj_last, bj_m, bj_mlast)
!!!          bj(m) = bj_pnt - &
!!!                 (bi_pnt-Bnd_pnt+epslon) * (bj_mlast-bj_m)/(bi_last-bi_pnt)
!!!          bi(m) = Bnd_pnt - epslon
!
!!!       ELSE
!
!!!          CALL Adjust_bj_2 (bj_pnt, bj_mlast, bj_m, bj_mlast)
!!!          bj_new_1 = bj_pnt - &
!!!              (bi_pnt-bnd_pnt+epslon)*(bj_mlast-bj_m)/(bi_last-bi_pnt )
!
!!!          CALL Adjust_bj_2 (bj_pnt, bj_next, bj_m, bj_mnext)
!!!          bj_new_2 = bj(m) - &
!!!              (bi_pnt-Bnd_pnt+epslon)*(bj_mnext-bj_m)/(bi_next-bi_pnt)
!
!!!          bj(m) = 0.5_rprec * (bj_new_1 + bj_new_2)
!!!          bi(m) = Bnd_pnt - epslon
!!!       END IF
!!!       bj(m) = Bjmod_real ( bj(m), jsize )
    END DO Loop_k
    RETURN
    END SUBROUTINE Edgfix
!
!
!
!
    SUBROUTINE Move_plasma_grid (dt, ie_ask)
    USE Rcm_variables
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
    INTEGER (iprec), INTENT (IN) :: ie_ask
!_____________________________________________________________________________
!   Subroutine to advance eta distribution for a time step
!   by a lifetime-based algorithm (raw doc dated 5/12/87)
!                                                                       
!   Last update: 05-11-88
!                01-29-96 ain ,min_j and calls to bndy added - frt
!   rws          06-05-97 etamov changed to reflect new use of
!                        eeta array in rcm697 version
!
!   CALLS:       BNDY, GNTRP, RATEFN, BJMOD
!   CALLED FROM: TSTEP1
!_____________________________________________________________________________
!
!
    REAL (rprec), DIMENSION (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc) :: &
                    eeta2, veff, &
                    dvefdi, dvefdj
    REAL (rprec) :: &
                    didt, djdt, biold, bjold, rate, mass_factor, r_dist,&
                    Bjmod_real, Intrp_2d_grid, Gntrp_2d, Bndy
    REAL (rprec) :: Plasmasphere_refill_rate, cexrat
    INTEGER (iprec) :: i, j, kc, ie
!
!
!
    DO kc = 1, kcsize

       ie = ikflavc (kc)
       IF (ie /= ie_ask) CYCLE

       veff = v + vcorot -vpar + alamc(kc)*vm
       mass_factor = SQRT (xmass(1) / xmass(ie))
!
!
!      Compute spatial derivatives of effective potential for given energy:
!      (only inside the modeling region, boundaries possibly included).
!
       CALL Deriv_i (veff, isize, jsize, imin_j, dvefdi)
       CALL Deriv_j (veff, isize, jsize, imin_j, 1.0E+25, dvefdj)
       WHERE (ABS(dvefdj) > 1.0E+24) 
          dvefdj = 0.0_rprec
          dvefdi = 0.0_rprec
       END WHERE
!
!
!      3. Find the position from which particles of given energy have come
!         from to reach a given grid point (i,j) in one time step, and
!         interpolate EETA at that position. Then estimate that EETA at
!         (i,j) point is the same as that value, corrected for loss.
!
!         Notice: we are moving plasma only inside the modeling region:
!         (therefore, derivatives dvefdi and dvefdj are used only inside).
!
       eeta2 (:,:) = eeta (:,:,kc)
       DO j = 1, jsize
       DO i = 1, isize-1

          IF (i <= imin_j(j) ) CYCLE

          didt   =   dvefdj (i,j) / fac (i,j)
          djdt   = - dvefdi (i,j) / fac (i,j)
!         biold  = MAX (REAL(i,rprec) - didt * dt , REAL(imin_j(j)))
          biold  = REAL(i,rprec) - didt * dt 
          bjold  = Bjmod_real (REAL(j,rprec) - djdt * dt, jsize )

!         rate   = Ratefn (fudgec(kc), alamc(kc), sini (i,j), bir (i,j), &
!                          vm (i,j), mass_factor)
!         eeta (i,j,kc) = Intrp_2d_grid (eeta2, biold, bjold, isize, jsize, jwrap, imin_j)*&

          r_dist = SQRT(xmin(i,j)**2+ymin(i,j)**2)
          rate   = 0.0

          IF (ie == 1) THEN
              IF (kc ==1 ) THEN 
!                Cold plasmaspheric electrons. Dont' do pitch angle scattering,
!                but compute refilling rates for the plasmasphere
                 rate = - Plasmasphere_refill_rate (r_dist, doy, &
                                      sunspot_number, eeta(i,j,kc), vm(i,j))
                 eeta(i,j,kc) = eeta(i,j,kc)-rate*dt
                 rate = 0.0
               ELSE
!                 Plasmasheet electrons, pitch-angle scattering loss rate:
                  rate = Ratefn (fudgec(kc), alamc(kc), sini(i,j),&
                                         bir (i,j), vm(i,j), mass_factor)
               END IF
          ELSE 
!            Positive ions, compute charge-exchange rate is it is on:
             IF (L_dktime) THEN
              rate = Cexrat (ie, ABS(alamc(kc))*vm(i,j), &
                                         R_dist, &
                                         sunspot_number, dktime, &
                                         irdk,inrgdk,isodk,iondk)
              END IF
          END IF

          eeta(i,j,kc) = Gntrp_2d (eeta2, isize, jsize, biold, bjold)*EXP(-rate*dt)
       END DO
       END DO
!
       CALL Wrap_around_ghostcells (eeta(:,:,kc), isize, jsize, n_gc)
!
    END DO
!
    RETURN
!
    CONTAINS
!
    FUNCTION Ratefn (fudgx, alamx, sinix, birx, vmx, xmfact)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
    REAL (rprec)              :: Ratefn
!                                                                       
!   Function subprogram to compute precipitation rate
!   Last update:  04-04-88
!
    Ratefn = 0.0466_rprec*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
    Ratefn = xmfact * ratefn
    RETURN
    END FUNCTION Ratefn
    END SUBROUTINE Move_plasma_grid
!
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
      SUBROUTINE Read_bfield ()
      USE Rcm_variables
      USE Rcm_io
      IMPLICIT NONE
!
!_____________________________________________________________________________
!
!     Subroutine to read all of the labels of the input
!     bfield models to determine the event times for each
!     of the records.  This information is placed in the
!     1-d array ibtime. These mark times are set in the
!     program creating the B-field arrays.
!     rws    3/19/97
!_____________________________________________________________________________
!
      INTEGER (iprec) :: n, nbf
      LOGICAL :: error_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
IF (itype_bf == 2 .OR. itype_bf == 3) THEN
  nbf = 1
  ALLOCATE (ibtime (nbf))
  ibtime = -999
ELSE IF(itype_bf == 1) THEN 
      n = 1
      DO
         CALL Read_array ('input/rcmxmin_inp', n, label, ARRAY_2D = xmin, &
                           ERROR_FLAG = error_flag, ASCI = asci_flag)
         IF (error_flag) EXIT
         n = n + 1
      END DO
      nbf = n - 1

      ALLOCATE (ibtime (nbf))

      DO n = 1, nbf
        CALL Read_array ('input/rcmxmin_inp', n, label, ARRAY_2D = xmin, ASCI=asci_flag)
        ibtime (n) = label%intg (6)
      END DO 
ELSE
   call CON_STOP('ILLEGAL ITYPE_BF IN READ_BFIELD')
END IF
called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_bfield
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
      SUBROUTINE Read_eta_on_bndy ()
      USE Rcm_variables
      USE Rcm_io
      IMPLICIT NONE
!_____________________________________________________________________________
!
!_____________________________________________________________________________
!
      INTEGER (iprec) :: n, n_t
      LOGICAL :: error_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
IF (i_eta_bc == 1) THEN
      n = 1
      DO
         CALL Read_array ('input/rcmetac_inp', n, label, ARRAY_1D = etac, &
                           ERROR_FLAG = error_flag, ASCI = asci_flag)
         IF (error_flag) EXIT
         n = n + 1
      END DO
      n_t = n - 1

      ALLOCATE (itime_etac (n_t), etac_inp(kcsize,n_t))

      DO n = 1, n_t
        CALL Read_array ('input/rcmetac_inp', n, label, ARRAY_1D = etac, ASCI=asci_flag)
        itime_etac (n) = label%intg (6)
        etac_inp (:,n) = etac
      END DO 
ELSE IF (i_eta_bc == 2 .OR. i_eta_bc == 3) THEN
      n_t = 1
      ALLOCATE (itime_etac (n_t), etac_inp(kcsize,n_t))
      itime_etac = 0
      etac_inp = 0.0
ELSE 
      call CON_STOP('ILLEGAL VALUE OF I_eta_bd ON INPUT')
END IF
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_eta_on_bndy
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
    FUNCTION Eta_lambda_vgamma ( kbeg, kend, kcbeg, kcend, gamma)  RESULT (out)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: kbeg, kend, kcbeg, kcend
    REAL (rprec), INTENT (IN) :: gamma
    REAL (rprec) :: out (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,iesize)
!______________________________________________________________________________
!
!   Subroutine computes quantity ETA*ABS(ALAM)*V**GAMMA at each grid point
!   for electrons and ions separately. KBEG, KEND, KCBEG, KCEND can be used
!   to restrict species to be included in the sum (use 1, ksize, 1, kcsize for
!   no restrictions). GAMMA is an input parameter:
!   ** if GAMMA = 0, then the computed sum is the adiabatic parameter PVGAMMA
!   ** if GAMMA = -5/3, compute energy density (or pressure without the 2/3 factor)
!   ** if GAMMA = -2/3, compute total energy of particles
!______________________________________________________________________________
!
    INTEGER (iprec) :: nbi, k, m, mbeg, mend, ipmax, ncount, i,j,n, ie,kc
    REAL (rprec) :: q, bimax, bicrss (100), charge, Bndy
!
!
    out = 0.0_rprec
!
    DO ie = 1, iesize
!
        IF (ie == 1) THEN
           charge = - 1.0_rprec
        ELSE
           charge = + 1.0_rprec
        END IF
!
!       I. Compute sum for plasma on inner edges, electrons and ions separately :
!
!       1. find bimax, the innermost limit of relevant particles:
!
        bimax = 0.0_rprec
        DO  k = kbeg, kend
           q = alam(k) / charge
           IF (q > 0.0_rprec) THEN
              mbeg = mpoint(k)
              mend = mbeg + npoint(k) - 1
              DO m = mbeg, mend
                 bimax = MAX (bimax, bi(m), Bndy (bndloc, jsize, bj(m)) )
              END DO
           END IF
        END DO
        ipmax = INT (bimax)
!
!
!       2. Compute sum at each grid pt between boundary and ipmax:
!
        DO  j = 1, jsize
        DO  i = CEILING (bndloc(j)), ipmax
           loop_120: DO k = kbeg, kend
              q = alam (k) / charge
              IF (q < 0.0_rprec) CYCLE loop_120
              IF (i == CEILING(bndloc(j)) ) THEN
!
!                New j=const grid line. Find its all crossings by particles of
!                species k and store in bicrss. First, initialize bicrss array:
!
                 bicrss = 0.0_rprec
                 CALL Crschk (j, k, bicrss, SIZE(bicrss), nbi)
!
!                nbi = # of times k inner edge crosses j grid line
!                bicrss(1:nbi) are i locations of crossings.
!
              END IF
!
!             Check if species k is present at grid pt(i,j):
!
              ncount = 0
              IF (nbi /= 0) THEN
                 DO n = 1, nbi
                    IF (bicrss (n) > REAL(i,rprec) ) ncount = ncount + 1
                 END DO
              END IF
!
!
!             If ncount is odd, species k is present at (i,j), so add contibution:
!
              IF (MOD (ncount, 2) == 1) THEN
                 out (i,j,ie) = out (i,j,ie) + eta (k) * ABS (alam(k))
              END IF
!
           END DO loop_120
        END DO
        END DO
!
!
!   II. Compute the sum for grid_based electrons or ions:
!
        DO j = 1, jsize
        DO i = 1, isize
           IF (REAL(i,rprec) < Bndy(bndloc, jsize, REAL(j)) ) CYCLE
           DO kc = kcbeg, kcend
              q = alamc(kc) / charge
              IF (q > 0.0_rprec) THEN
                 out (i,j,ie) = out (i,j,ie) + ABS (alamc(kc) * eeta(i,j,kc))
              END IF
           END DO
        END DO
        END DO
!
    END DO
!
!
    DO ie = 1, iesize
       CALL Wrap_around_ghostcells (out(:,:,ie), isize, jsize, n_gc)
    END DO
!
!
    DO ie = 1, iesize
       out (:,:,ie) = out(:,:,ie) * (vm (:,:)**((-3.0_rprec/2.0_rprec)*gamma))
    END DO
!
    RETURN
    END FUNCTION Eta_lambda_vgamma
!
!
!
!
    SUBROUTINE Crschk (j_line, k_edge, bicrss, n_bicrss, nbc)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec),  INTENT (IN) :: j_line, k_edge, n_bicrss
    INTEGER (iprec),  INTENT (OUT):: nbc
    REAL(rprec),      INTENT (OUT):: bicrss (n_bicrss)
!
!   last update: 09-17-86
!               01-29-96 frt - added ain array                          
!                                                                       
!   subroutine to find all crossings of grid line j_line by species
!   k_edge inner edge.
!   nbc = # of crossings found
!   bicrss (n,kk)= array containing the bi values of the crossings.
!   If a segment is totally outside the modeling region, it is not considered.
!
!   Dependency: bjmod, bndy
!
!
    INTEGER (iprec) :: mbeg, mend, n, m, mnext, nn, ntest, kdim, ndim
    REAL (rprec)    :: bj_line, b1, b2, b3, bi_temp, bi_test, &
                       Bjmod_real, Bndy
!
    kdim = SIZE (alam)
    ndim = SIZE (bicrss)
    IF (k_edge < 1 .OR. k_edge > kdim) THEN
       WRITE (*,*) 'CRSCHK: k_edge=', k_edge
       call CON_STOP('CRSCHK: k_edge IS OUT OF RANGE')
    END IF
    n = 0
    bj_line = REAL (j_line,rprec)
!                                                                       
    mbeg = mpoint (k_edge)
    mend = mbeg + npoint (k_edge) - 1
!                                                                       
    loop_10: DO m = mbeg, mend
       mnext = m + 1
       IF (m == mend) mnext = mbeg
       b1 = Bjmod_real (bj_line, jsize)
       b2 = Bjmod_real (bj(m), jsize  )
       b3 = Bjmod_real (bj(mnext), jsize)
       IF (ABS(b1-b2) < 1.E-5_rprec ) b2 = b2 + 1.E-5_rprec
       IF (ABS(b1-b3) < 1.E-5_rprec ) b3 = b3 + 1.E-5_rprec
!
!
!      check if b1 is between b2 and b3:
!
       IF ((b1 <= b2 .AND. b1 > b3) .OR. (b1 >= b2 .AND. b1 < b3)) THEN
!
!         we have a crossing. Is the segment inside the region?
!
          IF (bi(m) <= Bndy (bndloc, jsize, bj(m)) .AND. &
                             bi(mnext) <= Bndy (bndloc, jsize, bj(mnext)) )&
              CYCLE loop_10
!
!         Yes. Linearly interpolate to find crossing pt:
!
          n = n + 1
!
          IF (n > ndim) THEN
             WRITE (*,*) 'CRSCHK: k, j are:',k_edge, j_line
             call CON_STOP('CRSCHK: array bicrss is too small')
          ELSE
             bicrss (n) = bi(m) + (bi(mnext)-bi(m)) * (b1-b2)/(b3-b2)
          END IF
!
       END IF
    END DO loop_10
!                                                                       
    nbc = n
    IF (nbc > 1) THEN
!
!      order pts out from earth in bicrss array:
!
       DO n = 1, nbc
          bi_temp = bicrss (n)
          bi_test = 0.0_rprec
          DO nn = n, nbc
             IF (bicrss (nn) >= bi_test) THEN
                bi_test = bicrss (nn)
                ntest = nn
             END IF
          END DO
          bicrss (n)     = bi_test
          bicrss (ntest) = bi_temp
       END DO
    END IF
!                                                                       
    RETURN
    END SUBROUTINE Crschk
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_vdrop ()
      USE Rcm_variables
      IMPLICIT NONE
!_____________________________________________________________________________
!
!     Subroutine to read cross polar cap potential drops
!     and place them in vinput array.  Times are stored in
!     ivtime array.  nvmax is actual number of potential drop
!     values.  These results are used in subroutine getv to
!     interpolate in time to get potential at any time of
!     interest.
!     rws  03-20-97
!_____________________________________________________________________________
!
      INTEGER (iprec) :: nv, nvmax, istat
      LOGICAL         :: logical_flag
      LOGICAL, SAVE   :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
      INQUIRE (FILE = trim(NameRcmDir)//'input/rcmpcp_inp', EXIST = logical_flag)
      IF (.NOT.logical_flag ) call CON_STOP('READV: RCMPCP_INP not found')
      INQUIRE (UNIT = LUN, OPENED = logical_flag)
      IF (logical_flag) call CON_STOP('READV: LUN is already open')
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = trim(NameRcmDir)//'input/rcmpcp_inp', IOSTAT=istat)
      IF (istat /=0) call CON_STOP('ERROR OPENING rcmpcp_inp')
      nvmax = 0
      DO
         READ (LUN,*, END = 19 )
         nvmax = nvmax + 1
      END DO
  19  CLOSE (UNIT = LUN)
!
      ALLOCATE (ivtime (nvmax), vinput (nvmax), vinput_phase(nvmax) )
!
      OPEN (UNIT = LUN, STATUS ='OLD', FILE = trim(NameRcmDir)//'input/rcmpcp_inp') 
      DO nv = 1, nvmax
         READ (LUN, *) ivtime (nv), vinput (nv), vinput_phase(nv)
      END DO
      CLOSE (UNIT = LUN)
      called_already = .TRUE.
!
      RETURN 
      END SUBROUTINE Read_vdrop
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_kp ()
      USE Rcm_variables
      IMPLICIT NONE
!_____________________________________________________________________________
!
!     Subroutine to read Kp index values
!     and place them in kpinput array.  Times are stored in
!     ikptime array.  nkpmax is actual number of Kp
!     values.  These results are used in subroutine getkp to
!     interpolate in time to get Kp at any time of interest.
!     rws  03-20-97 stanislav 05-28-99
!_____________________________________________________________________________
!
      INTEGER (iprec) :: nkp, nkpmax, istat
      LOGICAL :: logical_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!                                                                       
      IF (called_already) RETURN
!
      INQUIRE (FILE = trim(NameRcmDir)//'input/rcmkp_inp', EXIST = logical_flag)
      IF (.NOT.logical_flag ) call CON_STOP('READKP: RCMKP_INP not found')
      INQUIRE (UNIT = LUN, OPENED = logical_flag)
      IF (logical_flag) call CON_STOP('READKP: LUN is already open')
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = trim(NameRcmDir)//'input/rcmkp_inp', IOSTAT=istat)
      IF (istat /= 0) call CON_STOP('ERROR OPENING rcmkp_inp')
      nkpmax = 0
      DO
         READ (LUN, *, END = 19)
         nkpmax = nkpmax + 1
      END DO
   19 CLOSE (LUN)
!
      ALLOCATE (ikptime (nkpmax), kpinput (nkpmax) )
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = trim(NameRcmDir)//'input/rcmkp_inp')
      DO nkp = 1, nkpmax
         READ (LUN, *) ikptime (nkp), kpinput (nkp)
      END DO
      CLOSE (UNIT = LUN)
      called_already = .TRUE.
!
      RETURN 
      END SUBROUTINE Read_kp
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_winds
      USE Rcm_variables, junk_v => v
      IMPLICIT NONE
!
!_____________________________________________________________________________
!                                                                       
!     last update 08-27-86        written by g.a.mantjoukis
!                                                                       
!     subroutine to compute pedersen and hall winds on
!     the rcm-specified grid
!                                                                       
!     iwind=0    no wind                                                
!     iwind=1    tarpley-type wind (calculatoin)                        
!     iwind=2    roble sq winds                                         
!                                                                       
!-------------------------------------------------------------
!                                                                       
      INTEGER (iprec) :: i, j, is, n, istat
      CHARACTER (LEN=80) :: form_string
      REAL    (rprec) :: s, sm, vnorm, sv, su, phi, v, u, ath, c8, c6, c4, c2, c1,&
                         v8, v6, v4, v2, v0, u6, u4, u2, u0
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
sw = 0.0
pwe = 0.0
hwn = 0.0
hwe = 0.0
pwn =0.0
called_already = .TRUE.
return
      IF (iwind == 0) THEN
         sw   = 0.0_rprec
         pwe  = 0.0_rprec
         hwn  = 0.0_rprec
         hwe  = 0.0_rprec
         pwn  = 0.0_rprec
!          OPEN (LUN, FILE = trim(NameRcmDir)//'rcmwind', STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT=istat)
!          IF (istat /= 0) call CON_STOP('ERROR OPENING rcmwind')
!
!           READ (LUN,'(I10.10)') n
!           IF (n /= isize*jsize) call CON_STOP('size mismatch in read_winds')
!           READ (LUN,'(A80)') form_string
!           DO j = 1, jsize
!           DO i = 1, isize
!              READ (LUN, form_string) pwe(i,j), hwn(i,j), hwe(i,j), pwn(i,j)
!           END DO
!           END DO
!
!           READ (LUN,'(I10.10)') n
!           IF (n /= jsize) call CON_STOP('size mismatch in read_winds')
!           READ (LUN,'(A80)') form_string
!           DO j = 1, jsize
!              READ (LUN, form_string) sw(j)
!           END DO
!
!        CLOSE (LUN)
      ELSE IF (iwind == 1) THEN
!
          pwn = 0.0_rprec
          pwe = 0.0_rprec
!
          u0 =   0.4225615_rprec
          u2 = - 1.3858640_rprec
          u4 = - 1.1390120_rprec
          u6 = - 0.4121196_rprec
!                                                                       
          v0 = - 0.4549289_rprec
          v2 = + 3.0388490_rprec
          v4 = - 3.6561400_rprec
          v6 = - 4.5478570_rprec
          v8 = - 1.9232250_rprec
    !
          DO j = 1, jsize
          DO i = 1, isize
             IF (colat (i, j) > 1.029744 .AND. &
                 colat (i, j) < 1.064651         ) CYCLE
             c1 = COS (colat (i, j) )
             c2 = c1**2
             c4 = c1**4
             c6 = c1**6
             c8 = c1**8
             ath = 130.0 / (0.25_rprec - c2)
    !
             u = ath * 3.0_rprec * c1 * (u0 + u2 * c2 + u4 * c4 + u6 * c6)
             v = ath * (v0 + v2 * c2 + v4 * c4 + v6 * c6 + v8 * c8)
    !
             phi = aloct (i, j) + pi
             su = sin (phi + 250.0 * pi / 180.0)
             sv = sin (phi + 340.0 * pi / 180.0)
    !
    !        there is a minus sign in unorm since sin(phi)=-sin(phi+pi)
             vnorm = - 13.97155_rprec
             pwn (i, j) = - u * su / vnorm
             pwe (i, j) =   v * sv / vnorm
          END DO
          END DO
    !
    !     fixing singularity at colat=60.if for any j,more than
    !     one i's have 59<colat(i,j)<61 the fixing is no good.
    !     in latter case the program will stop.
    !
          DO_60: DO j = 1, jsize
             is = 0
             DO_50: DO i = 1, isize
             IF (pwn (i, j) /= 0.0_rprec) CYCLE DO_50
             is = is + 1
             IF (is /= 1) call CON_STOP('singularity in winds needs fixing')
             sm = colat (i + 1, j) - colat (i - 1, j)
             s = colat (i, j) - (colat (i + 1, j) + colat (i - 1, j) )*0.5_rprec
             pwn (i, j) = (0.5_rprec + s / sm) * pwn (i + 1, j) &
                        + (0.5_rprec - s / sm) * pwn (i - 1, j)
             pwe (i, j) = (0.5_rprec + s / sm) * pwe (i + 1, j) &
                        + (0.5_rprec - s / sm) * pwe (i - 1, j)
             END DO DO_50
          END DO DO_60
!
          hwn = pwn
          hwe = pwe
!
      END IF
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_winds
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Add_point (k, m, addi, addj, is )
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: k, m, is
    REAL (rprec),    INTENT (IN) :: addi, addj
!
!_____________________________________________________________________________
!   last update:  03-31-88             by:rws
!                                                                       
!   this subroutine adds a new pt at m+1 and moves
!   subsequent pts
!   by 1.  bi(m+1)=addi;  bj(m+1)=addj.
!                                                                       
!   if(is.lt.0) then the pt will be encoded as a held pt.
!   CALLS:       NCODE
!   CALLED FROM: EXPAND_EDGE
!_____________________________________________________________________________
!
!
    INTEGER (iprec) :: mend, mp, mm, mx, kk, mpp, newitk
    REAL (rprec) :: Bjmod_real
!
!
    mend = mpoint (ksize) + npoint (ksize) - 1
    mp = m + 1
    IF (m < mend) THEN
!       if adding not after last point, move pnts from m+1 to last right:
       DO mm = m, mend
          mx          = mend - (mm - mp)
          itrack (mx) = itrack (mx - 1)
          bi (mx)     = bi (mx - 1)
          bj (mx)     = bj (mx - 1)
          etab (mx)   = etab (mx - 1)
       END DO
    END IF
!
!
!   New point is going to be at m = mp:
!
    bi (mp)     = addi                          ! set bi for new point
    bj (mp)     = Bjmod_real (addj, jsize) ! set bj for new point
    CALL Ncode (itrack, SIZE(itrack), is, newitk)
    itrack (mp) = newitk                        ! set itrack for new point
    npoint (k) = npoint (k) + 1                 ! k-edge is longer now
!
    DO kk = 1, ksize   ! for edges to right of new point, shift beginnings:
       IF (mpoint (kk) >= mp) mpoint (kk) = mpoint (kk) + 1
    END DO
    mpoint (ksize + 1) = mpoint (ksize) + npoint (ksize)
!                                                                       
    mpp = mp + 1
    IF (mpp > mpoint (k) + npoint (k) - 1)     mpp = mpoint (k)
    etab (mp) = 0.5_rprec * (etab (m) + etab (mpp) )  ! set etab for new point
!                                                                       
    IF (mpoint (ksize) + npoint (ksize) - 1 >= nptmax) THEN
 print*,ksize,nptmax
print*,mpoint
print*,npoint
         call CON_STOP('ADDPT: # of pts exceeded nptmax-1')
    END IF
!
    RETURN
    END SUBROUTINE Add_point
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Expand_edge ( k )
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: k
!
!_____________________________________________________________________________
!                                                                       
!  last update: 03-31-88              by: rws                           
!               01-25-96                  frt - added ain               
!                                                                       
!  this subroutine marches around inner edge k, deciding 
!  whether to add new inner edge tracer pts.  a new point
!  is inserted  between pt m and pt m+1 if              
!                                                                       
!      (bi(m+1)-bi(m))**2+(bj(m+1)-bj(m))**2 > dstmax**2                
!                                                                       
!  Actual point addition procedure occurs in routine ADDAPT     
!                                                                       
!> Skeleton of this subprogram (algorithm):
!
!>    Find 1st point of the edge and assign it to m.
!>    Enter:
!>    Loop do:
!>       Find mnext
!>       Check distance between m and mnext: 
!>       IF (distance is not too large) THEN
!>          m = m + 1
!>          CYCLE loop with new m
!>       ELSE
!>          I.  Estimate new location of point
!>          II. Add new point between m and m+1
!>          CYCLE Loop for the same m to see if more
!>              points need to be added.
!>       END IF
!>    END loop DO
!   CALLS:       ADD_POINT, CORECT, GNTRP, DCODE, ADJUST_BJ, BJMOD
!   CALLED FROM: MAIN
!_____________________________________________________________________________
!
!
    INTEGER (iprec) :: m, mbeg, mend, mnext, ism, isp, is
    REAL (rprec)    :: sqmax, a3, a4, distsq, biz, bjz, biq, &
                       bjq, bi1, bj1, vm1, rho, vmbar, di, dj, &
                       bi2, bj2, vm2, vmdfmx, &
                       Bjmod_real, Gntrp_2d
!
!
    sqmax = dstmax**2
    mbeg = mpoint (k)
    m = mbeg
!
    loop_10: DO
       mend = mpoint (k) + npoint (k) - 1
       IF (m > mend) EXIT loop_10             ! Normal exit from loop
       mnext = m + 1
       IF (m == mend) mnext = mbeg
!
       CALL Adjust_bj_2 (bj (mnext), bj (m), a3, a4 )
       distsq = (bi(mnext)-bi(m))**2 + (a3-a4)**2
!
       IF (distsq <= sqmax) THEN
!
!         Points m and m+1 are close enough, proceed to check the next pair:
!
          m = m + 1
!
       ELSE
!
!         Points m and m+1 are too far apart, need to add a new point.
!         First estimate of new pt location:
!
          biz = 0.5_rprec * (bi (m) + bi (mnext) )
          bjz = Bjmod_real (0.5_rprec * (a3 + a4), jsize)
!
!
!         Calculate correction factor rho:
!
          biq = biz
          bjq = bjz
          CALL Corect (biq, bjq, m, mnext, a3, a4, di, dj, rho, vmbar)
!
!         Estimate corrected position only if correction is
!         not too large; otherwise initial estimate of biz and
!         bjz will be used to add a point:
!
          IF ( ABS(rho) <= rhomax) THEN
!                                                                       
!            Compute new estimate of point location:
!
             bi1 = biz + rho * dj
             bj1 = bjz - rho * di
             vm1 = Gntrp_2d (vm, isize, jsize, bi1, bj1)
             vmdfmx = vmfact * vmbar
             IF (ABS(vm1-vmbar) < vmdfmx) THEN
!
!               VM1 is close enough to vmbar, will use BI1 and BJ1
!
                biz = bi1
                bjz = bj1
!
             ELSE
!
!               Go through correction procedure again:
!
                biq = bi1
                bjq = bj1
                CALL Corect (biq, bjq, m, mnext, a3, a4, di, dj, rho, vmbar)
                IF (ABS(rho) <= rhomax) THEN
                   bi2 = bi1 + rho * dj
                   bj2 = bj1 - rho * di
                   vm2 = Gntrp_2d (vm, isize, jsize, bi2, bj2)
                   vmdfmx = vmfact * vmbar
                   IF (ABS(vm2-vmbar) < vmdfmx) then
!
!                     will use bi2 and bj2 for new location
!
                      biz = bi2
                      bjz = bj2
                   END IF
                END IF
             END IF
          END IF
          bjz = Bjmod_real (bjz, jsize )
!                                                                       
!         II. Adding a new point:
!
!         If both points m and mnext are held, new point should be held too:
!                                                                       
          is = 1
          CALL Dcode (itrack, SIZE(itrack), m,     ism)
          CALL Dcode (itrack, SIZE(itrack), mnext, isp)
          IF (ism < 0 .AND. isp < 0)  is = -1
!                                                                       
          CALL Add_point (k, m, biz, bjz, is )
       END IF
!
    END DO loop_10
!
    RETURN
    END SUBROUTINE Expand_edge
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Corect (biq, bjq, m, mnext, a3, a4, di, dj, rho, vmbar)
    USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: m, mnext
    REAL (rprec), INTENT (IN)  :: biq, bjq, a3, a4
    REAL (rprec), INTENT (OUT) :: rho, vmbar, di, dj
!
!_____________________________________________________________________________
!   last update: 05-16-85                by:rws
!                                                                       
!   this subroutine computes a correction factor rho
!   to be used in adding new inner edge tracer pts.
!   see 'rawolf.docmir2.text' for explanation
!   CALLS:       GNTRP, BJMOD
!   CALLED FROM: ADD_POINT
!_____________________________________________________________________________
!
!
    INTEGER (iprec) :: i, j, jnext, &
                       Bjmod_int
    REAL (rprec) ::  Gntrp_2d, Bjmod_real
    REAL (rprec)    :: a, d, b1, b2, b3, b, vvmp, vvmm, vmb, c
!
!
    i = MIN (INT (biq), SIZE(vm, DIM=1)-1)
    j = Bjmod_int ( INT(bjq), jsize )
!                                                                       
    di = bi (mnext) - bi (m)
    dj = a3 - a4
!                                                                       
    jnext = Bjmod_int (j + 1, jsize)
    d     = vm (i,j) + vm (i+1,jnext) - vm (i,jnext) - vm (i+1,j)
    a     = d * di * dj
    b1    = - dj * (vm (i + 1, j) - vm (i, j) )
    b2    = di * (vm (i, jnext) - vm (i, j) )
    b3    = d * ( (biq - i) * di - (bjq - INT(bjq) ) * dj)
    b     = b1 + b2 + b3
    vvmp  = Gntrp_2d (vm, isize, jsize, bi (mnext), bj (mnext))
    vvmm  = Gntrp_2d (vm, isize, jsize, bi (m), bj (m))
    vmbar = 0.5_rprec * (vvmp + vvmm)
!                                                                       
!
    vmb   = Gntrp_2d (vm, isize, jsize, biq, Bjmod_real (bjq,jsize))
    c     = vmbar - vmb
!                                                                       
    IF (ABS(a) > 1.E-5 .AND. b**2 > 1.E-5 .AND. 4.0*a*c/b**2 <= 1.0_rprec ) THEN
       rho = (0.5_rprec*b/a) * (- 1.0_rprec + SQRT (1.0_rprec - (4.0_rprec*a*c/b**2) ) )
    ELSE
      rho = 1.0E12_rprec
    END IF
    RETURN
    END SUBROUTINE Corect
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Efield (ipatt, pcp, deqdt, a, b, dx, dy, colat, aloct, &
                         mode, icntrl, itop, v, vbndry)
      USE Rcm_variables, ONLY : iprec, rprec, dtr, rtd, pi, pi_two
!
!     this subroutine returns values of the potential v at all grid
!     points, specified by colat(i,j) and aloct(i,j).
!
!
!    general call parameters for subroutine efield
!
!       ipatt = heppner-maynard pattern no.  must be specified as an
!               integer between 1 and 7.
!               ipatt=1 means pattern a, for bz south, by=0.
!               ipatt=2 means pattern bc, for bz south, bynorth>0
!               ipatt=3 means pattern de, for bz south,bynorth<0
!               ipatt=4 means pattern bcp,twisted bc, for bz>0
!               ipatt=5 means pattern bcpp,twisted bc,bz strong>0
!               ipatt=6 means patt. dep, twisted de,bz weak >0
!               ipatt=7 means patt. depp, twisted de,bz strong>0
!
!       pcp = polar-cap potential drop in volts.  must be specified
!             unless using unscaled heppner-maynard potential(icntrl=-1)
!
!       deqdt = estimated d/dt of latitude of equatorward edge of
!               auroral zone at local midnight, in degrees/hour.  this
!               is a parameter in low-latitude e-field model.  not
!               needed if icntrl.le.0.
!
!      vectors a(3), b(3), dx(3), dy(3) describe the ellipses that form
!      the boundaries between regions 1, 2, and 3.  boundary 1 is the
!      equatorward edge of the polar cap.  boundary 2 is the equatorward
!      edge of region 1 or the equatorward edge of the field-reversal
!      region.  boundary 3 is the shielding layer.
!       a(l) = radius of ellipse measured in x(sunward) direction.
!       b(l) = radius of ellipse measured in y(duskward) direction.
!       dx(l) = sunward displacement of coord.system center from pole.
!       dy(l) = duskward displ. of coord. system center from pole.
!      to use efield with rcm, choose a(2) =.5*(colat(imin,noon)+
!        colat(imin,midnt)), b(2)=.5*(colat(imin,dusk)+colat(imin,dawn))
!      dx(2) = -offset*180./pi, dy(2) = 0.  ellipse parameters for
!      boundaries 1 and 3 need not be specified. All a's, b's, dx's and
!      dy's are in degrees.
!
!      colat and aloct are the usual magnetic colatitude and magnetic
!      local time angles in radians.
!
!      mode is a dummy parameter, not currently used.
!
!      imax, jmax, jwrap have their usual meaning, same as in rcm.
!
!      itmdim = number of time labels in the v-matrix.  for use with
!               rcm, it should be set equal to 1.
!
!      itmcur = time-label number of v-matrix currently being computed.
!               for use with rcm, it should be set equal to 1.
!
!      icntrl = 1 means different formulas are used in regions 1,2,3.
!             = 0 means heppner-maynard formula is used for all
!                 latitudes, but scaled to externally specified pcp and
!                 ellipse parameters.
!             = -1 means heppner-maynard formula is used for all
!                  latitudes, unscaled.
!             = -2 is the mode used for use with rcm.  heppner-maynard
!                  is used poleward of ellipse 2, and on ellipse 2,
!                  scaled to fit specified potential drop and dimensions
!                  of ellipse 2.
!
!      itop = 1 means that a(1), b(1), dx(1), dy(1) are computed
!               internally in efield.  otherwise they must be passed to
!               efield from outside.  itop is normally = 1 for use in
!
!               msm.  itop should be set to 1 for use with rcm.
!      v = potential matrix computed in efield.  it is dimensioned at
!          latdim x ltdim x itmdim, where latdim and ltdim are specified
!          in parameter statement in efield and itmdim is passed.
!
!
!
!
!     subroutine efield calls subroutines epot(heppner-maynard),
!     low(low-latitudes), aurora(region 2), reg1(region 1).
!
!
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: icntrl, itop, ipatt, mode
      REAL (rprec), INTENT (IN) :: pcp, deqdt, colat (:,:), aloct (:,:)
      REAL (rprec), INTENT (IN OUT) :: a(3), b(3), dx(3), dy(3), v (:,:)
      REAL (rprec), INTENT (OUT) :: vbndry (:)
!
!
      REAL (rprec):: thetaa (SIZE(vbndry)),thetab (SIZE(vbndry)), &
!             thetac (SIZE(vbndry)),&
              coef (18,18), xco (18,18), ahm (2,7), &
              bhm (2,7), dxhm(2,7), dyhm(2,7), vmin(7), vmax(7),&
              fg, al, gmlat, glong, vtemp, &
              Thet
      INTEGER (iprec) :: ieff, l, j, i, nnmax, idim, jdim
!
!     parameter values
!
!       ntape = file number from which maynard-rich coefficients are
!               read.  dataset name is 'rawolf.efcoef.data'.
!
!
!     Stanislav: since now the grid is such that aloct(i,j) is
!     the same for all i and given J, ieff parameter is not
!     relevant. But watch out if this changes later. Also, the
!     boundary is assumed to be an ellipse (since we are
!     calling this subroutine!), but it does not have to 
!     coincide with the integer grid line. Location of bndy
!     is given fully by the parameters of the ellipse.
!
      ieff   = 1
      idim   = SIZE (v, DIM = 1)
      jdim   = SIZE (v, DIM = 2)
!
      IF (icntrl /= -2.OR.itop /= 1) call CON_STOP('EFIELD IS FOR RCM ONLY')
!
      IF (ABS(deqdt) > 0.0 .OR. mode /= 0) call CON_STOP('EFIELD ERROR')
      nnmax = 16
!
!
      CALL Input (coef, ipatt, xco, ahm, bhm, dxhm, dyhm, vmin, vmax)
!
!
!   fg = scaling factor that converts hm potential to
!         present situation
!
      IF (icntrl /= -1) THEN
         fg = pcp/(vmax(ipatt)-vmin(ipatt))
      ELSE
         fg = 1.0_rprec
      END IF
      IF (ABS(fg - 1.0_rprec) > EPSILON(1.0_rprec)) &
           call CON_STOP('ERROR IN EFIELD')
!
!
!     compute a(1),b(1),dx(1),dy(1) if itop = 1.
!
      IF (itop == 1) THEN
          a(1) = a(2) * ahm(1,ipatt) / ahm(2,ipatt)
          b(1) = b(2) * bhm(1,ipatt) / bhm(2,ipatt)
         dx(1) = dx(2) + a(2) * (dxhm(1,ipatt)-dxhm(2,ipatt)) /&
                                ahm(2,ipatt)
         dy(1) = dy(2) + b(2) * (dyhm(1,ipatt)-dyhm(2,ipatt)) /&
                                bhm(2,ipatt)
      END IF
!
!!    vmaxx = vmax(ipatt)*fg
!!    vminn = vmin(ipatt)*fg
!!    sa=sin((a(3)-dx(3))*pi/180.)
!!    vpenet = -13112.*g*deqdt*sa
!!    vbar = 0.6*vmaxx + 0.4*vminn -0.22*vpenet
!
!
!    computation of boundary locations as fcns of j.
!
!   guard against irrelevant a's or b's being zero or negative.
!
      DO  l=1,3
         if(a(l).lt.0.001) a(l)=0.001
         if(b(l).lt.0.001) b(l) = 0.001
      END DO
!
!
      DO j = 1, jdim
         al        = aloct (ieff,j)
         thetaa(j) = thet(a(1),b(1),dx(1),dy(1),al)*DTR
         thetab(j) = thet(a(2),b(2),dx(2),dy(2),al)*DTR
!        thetac(j) = thet(a(3),b(3),dx(3),dy(3),al)*DTR
      END DO
!
!     main do loop over grid
!
!
!
      DO j = 1, jdim
      DO i = 1, idim
!
         IF (colat(i,j) <= thetaa(j)) THEN
!
!          Region zero: polar cap:
!
            glong = MOD (aloct(i,j)+pi, pi_two)
            gmlat = 90.0 - colat(i,j)*RTD
            CALL Epot (coef, gmlat, glong, vtemp, &
                       ipatt, nnmax, a, b, dx, dy, xco, &
                       ahm, bhm, dxhm, dyhm, vmin, vmax, pcp, icntrl)
            v (i,j) = vtemp
!
         ELSE IF (colat(i,j) < thetab(j)+5.0E-4 .AND. &
                  colat(i,j) > thetaa(j) .AND. &
                  icntrl == -2) THEN
!
!          use heppner-maynard directly in region 1 for 
!          use with rcm.
!
            glong = MOD (aloct(i,j) + pi, pi_two)
            gmlat = 90.0 - colat(i,j)*RTD
            CALL Epot (coef, gmlat, glong, vtemp, &
                       ipatt, nnmax, a, b, dx, dy, xco, &
                       ahm, bhm, dxhm, dyhm, vmin, vmax, pcp, icntrl)
            v ( i,j ) = vtemp
         END IF
!
!
      END DO
      END DO
!
!
!   segment for defining boundary potentials and e-fields
!
!     potential at b
!
      DO j = 1, jdim
         glong = MOD (aloct(ieff,j)+pi, pi_two)
         gmlat = 90.0 - thetab(j)*RTD
         CALL Epot (coef, gmlat, glong, vbndry(j), ipatt,&
                    nnmax, a, b, dx, dy, xco, ahm, bhm, &
                    dxhm, dyhm, vmin, vmax, pcp, icntrl)
      END DO
!
!
      RETURN
      END SUBROUTINE Efield
!
!
!
      SUBROUTINE Input (coef, ipatt, xco, ahm, bhm, dxhm, dyhm, vmin, vmax)
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: ipatt
      REAL (rprec), INTENT (OUT) :: coef (18,18), dxhm(2,7), &
                            dyhm(2,7), xco(18,18), &
                            ahm(2,7), bhm(2,7), &
                            vmin(7), vmax(7)
!
!
      INTEGER (iprec) :: ip, n, i, j, ia, mmm, nnn
      REAL (rprec):: ccf
!     CHARACTER (LEN=50) :: mlbl
!
!  Subroutine that provides coefficients for Heppner-Maynard
!  polar cap convection model for use with RCM. yco arrays
!  holds Legendre coefficients (same for all HM models),
!  these come from the HM model.
! 
!  AAHM, BBHM, DDXHM, DDYHM arrays hold parameters of two
!  ellipses for each HM model. First ellipse (with index 1) 
!  is the poleward boundary of field-reversal region, and
!  second ellipse is the equatorward boundary of that region.
!  Second ellipse corresponds to the main RCM bounday.
!  VVMAX and VVMIN hold maximun (dawn) and minimun (dusk) 
!  potentials for each HM model.  All these arrays (as I
!  understand) were obtained by Bob Spiro.
!
!
      REAL (rprec), PARAMETER :: ycoTransposed (18,18) = RESHAPE ( (/ &
 .282095e+00, .488603e+00, .109255e+01, .228523e+01, .468333e+01, .951188e+01, .192265e+02, .387523e+02, .779645e+02, &
 .156658e+03, .314501e+03, .630964e+03, .126523e+04, .253611e+04, .508196e+04, .101809e+05, .203918e+05, .408366e+05,&
!
 .488603e+00, .488603e+00, .546274e+00, .144531e+01, .331161e+01, .719031e+01, .151999e+02, .316411e+02, .652298e+02,&
 .133599e+03, .272366e+03, .553392e+03, .112151e+04, .226837e+04, .458082e+04, .923904e+04, .186151e+05, .374743e+05,&
!
 .946175e+00, .109255e+01, .546274e+00, .590044e+00, .177013e+01, .440314e+01, .101333e+02, .223736e+02, .481754e+02,&
 .102038e+03, .213661e+03, .443701e+03, .915709e+03, .188083e+04, .384866e+04, .785168e+04, .159791e+05, .324537e+05,&
!
 .186588e+01, .228523e+01, .144531e+01, .590044e+00, .625836e+00, .207566e+01, .555021e+01, .134918e+02, .310971e+02,&
 .693209e+02, .151081e+03, .324033e+03, .686782e+03, .144253e+04, .300864e+04, .623988e+04, .128827e+05, .264983e+05,&
!
 .370249e+01, .468333e+01, .331161e+01, .177013e+01, .625836e+00, .656382e+00, .236662e+01, .674590e+01, .172496e+02,&
 .414272e+02, .955522e+02, .214328e+03, .471128e+03, .102002e+04, .218269e+04, .462762e+04, .973844e+04, .203694e+05, &
!
 .736787e+01, .951188e+01, .719031e+01, .440314e+01, .207566e+01, .656382e+00, .683184e+00, .264596e+01, .798499e+01,&
 .213929e+02, .534153e+02, .127330e+03, .293800e+03, .661878e+03, .146420e+04, .319336e+04, .688612e+04, .147131e+05,&
!
 .146845e+02, .192265e+02, .151999e+02, .101333e+02, .555021e+01, .236662e+01, .683184e+00, .707163e+00, .291571e+01,&
 .926339e+01, .259102e+02, .671087e+02, .165101e+03, .391572e+03, .903721e+03, .204248e+04, .454057e+04, .996084e+04,&
!
 .292940e+02, .387523e+02, .316411e+02, .223736e+02, .134918e+02, .674590e+01, .264596e+01, .707163e+00, .728927e+00,&
 .317732e+01, .105778e+02, .307916e+02, .825507e+02, .209304e+03, .509767e+03, .120459e+04, .278052e+04, .629979e+04,&
!
 .584734e+02, .779645e+02, .652298e+02, .481754e+02, .310971e+02, .172496e+02, .798499e+01, .291571e+01, .728927e+00,&
!  in line above, I changed the first constant from .30... to 
!    .31... since Heppner-Maynard official release has it
!    this way, in constract to Bob Spiro's verion. Stanislav
!    March 8 1999.
!
 .748901e+00, .343190e+01, .119255e+02, .360281e+02, .997819e+02, .260366e+03, .650553e+03, .157290e+04, .370647e+04,  &
!
 .116766e+03, .156658e+03, .133599e+03, .102038e+03, .693209e+02, .414272e+02, .213929e+02, .926339e+01, .317732e+01,&
 .748901e+00, .767395e+00, .368030e+01, .133043e+02, .416119e+02, .118840e+03, .318704e+03, .816138e+03,  .201755e+04, &
!
 .233240e+03, .314501e+03, .272366e+03, .213661e+03, .151081e+03, .955522e+02, .534153e+02, .259102e+02, .105778e+02, &
.343190e+01, .767395e+00, .784642e+00, .392321e+01, .147120e+02, .475361e+02, .139761e+03, .384731e+03,   .100877e+04, &
!
 .465998e+03, .630964e+03, .553392e+03, .443701e+03, .324033e+03, .214328e+03, .127330e+03, .671087e+02, .307916e+02, &
 .119255e+02, .368030e+01, .784642e+00, .800822e+00, .416119e+01, .161472e+02, .537941e+02, .162579e+03, .458849e+03, &
!
 .931187e+03, .126523e+04, .112151e+04, .915709e+03, .686782e+03, .471128e+03, .293800e+03, .165101e+03, .825507e+02, &
 .360281e+02, .133043e+02, .392321e+01, .800822e+00, .816077e+00, .439471e+01, .176082e+02, .603802e+02, .187325e+03, &
!
 .186100e+04, .253611e+04, .226837e+04, .188083e+04, .144253e+04, .102002e+04, .661878e+03, .391572e+03, .209304e+03,&
 .997819e+02, .416119e+02, .147120e+02, .416119e+01, .816077e+00, .830522e+00, .462415e+01, .190939e+02, .672889e+02, &
!
 .371962e+04, .508196e+04, .458082e+04, .384866e+04, .300864e+04, .218269e+04, .146420e+04, .903721e+03, .509767e+03,&
 .260366e+03, .118840e+03, .475361e+02, .161472e+02, .439471e+01, .830522e+00, .844251e+00, .484985e+01, .206029e+02,&
!
 .743510e+04, .101809e+05, .923904e+04, .785168e+04, .623988e+04, .462762e+04, .319336e+04, .204248e+04, .120459e+04, &
 .650553e+03, .318704e+03, .139761e+03, .537941e+02, .176082e+02, .462415e+01, .844251e+00, .857341e+00, .507210e+01, &
!
 .148629e+05, .203918e+05, .186151e+05, .159791e+05, .128827e+05, .973844e+04, .688612e+04, .454057e+04, .278052e+04,&
 .157290e+04, .816138e+03, .384731e+03, .162579e+03, .603802e+02, .190939e+02, .484985e+01, .857341e+00, .869857e+00,&
!
 .297130e+05, .408366e+05, .374743e+05, .324537e+05, .264983e+05, .203694e+05, .147131e+05, .996084e+04, .629979e+04,&
 .370647e+04, .201755e+04, .100877e+04, .458849e+03, .187325e+03, .672889e+02, .206029e+02, .507210e+01,   .869857e+00&
  /), (/18,18/)  )
!
      REAL, PARAMETER :: aahmTransposed(7,2) = RESHAPE ( (/ &
            17.45,12.13,16.06,13.72,14.79,12.82,15., &
            20.00,16.17,18.88,19.31,18.3,16.6,17.93 /), &
            (/7,2/) )
      REAL, PARAMETER :: bbhmTransposed(7,2) = RESHAPE ( (/ &
            14.26,13.78,14.31,13.78,14.73,11.70,12.07, &
            16.97,16.7,17.07,17.23,17.93,17.5,16.6 /), &
            (/7,2/) )
      REAL, PARAMETER :: ddxhmTransposed(7,2) = RESHAPE ( (/ &
            -2.66,-5.53,-3.19,-.11,-2.45,-1.65,-2.23, &
            -3.30,-7.45,-4.1,-1.76,-2.77,-2.98,-3.14/), &
            (/7,2/) )
      REAL, PARAMETER :: ddyhmTransposed(7,2) = RESHAPE ( (/ &
            1.60,.8,1.12,.05,.48,3.4,1.22, &
            1.33,.32,1.54,.85,.69,1.86,.21 /), &
            (/7,2/) )
      REAL, PARAMETER :: vvmax(7)= (/ &
            34007.,55354.,14390.,11287.,9329.,13249., 13221. /)
      REAL, PARAMETER :: vvmin(7)= (/ &
            -42280.,-16003.,-60935.,-16250.,-12947., -14428.,-13460./)
!
!
      DO ip = 1,7
         DO  n = 1, 2
            ahm(n,ip)  = aahmTransposed(ip,n)    
            bhm(n,ip)  = bbhmTransposed(ip,n)   
            dxhm(n,ip) = ddxhmTransposed(ip,n)
            dyhm(n,ip) = ddyhmTransposed(ip,n)
         END DO
         vmax(ip) = vvmax(ip)
         vmin(ip) = vvmin(ip)
       END DO
!
!
!   Initialize coef from unit ntape:
!
      DO i = 1, 18
         DO j = 1, 18
            coef (i,j) = 0.
            xco  (i,j) = ycoTransposed (j,i)
         END DO
      END DO
!
      OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'HMRCOEF', STATUS = 'OLD')
      DO ia = 1, ipatt
         READ (LUN, FMT = '(A)')
         DO
            READ (LUN,*) nnn,mmm,i,j,ccf
            IF (nnn == -1) EXIT
            IF (i > 18 .OR. j > 18) call CON_STOP('out of bounds in INPUT')
            IF (mmm < 0 ) call CON_STOP('MMM IN EPOT')
            coef (i,j) = ccf
         END DO
      END DO
      CLOSE (UNIT = LUN)
!
      RETURN
      END SUBROUTINE Input
!
!
!
      SUBROUTINE Epot (coef, tlat, tlon, value,  ipatt, nnmax, a, b, &
                       dx, dy, xco, ahm, bhm, dxhm, dyhm, vmin, vmax, &
                       pcp, icntrl)
      USE Rcm_variables, ONLY : iprec, rprec, pi
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: ahm(2,7), bhm(2,7), dxhm(2,7), dyhm(2,7),&
                                   vmin(7), vmax(7), a(3), b(3), dx(3), dy(3),&
                                   coef(18,18), xco (18,18)
      REAL (rprec), INTENT (IN) :: tlat, tlon, pcp
      REAL (rprec), INTENT (OUT) :: value
      INTEGER (iprec), INTENT (IN) :: icntrl, ipatt, nnmax
!
      REAL (rprec) :: dp(18,18), sp(18), cp(18),  &
!                     fn(18), fm(18),  &
                      p (18,18) = 0.0_rprec, const(18,18)= 0.0_rprec
!
!
      INTEGER (iprec) :: nmax, n, m
      REAL (rprec)   :: cph, sph, st, ct, beta, alpha, xlon, xl, xlat, &
                        xcol, yy, xx, tlong, tcol, pol
      REAL (rprec), PARAMETER :: cmin=50.0_rprec, cmax = 90.0_rprec
!
!
      nmax = nnmax
!
      value = -1.0E-9_rprec
      IF (tlat <= cmin) THEN
         WRITE (*,*) 'latitude ',tlat,' out of model range'
         RETURN
      END IF
!
!
!   segment for scaling grid-pt. location to fit h-m ellipse.
!
      tcol = 90.0_rprec - tlat
      tlong = tlon - pi
      IF (icntrl /= -1) THEN
         xx = dxhm(1,ipatt)+ahm(1,ipatt)*(tcol*COS(tlong)-dx(1))/a(1)
         yy = dyhm(1,ipatt)+bhm(1,ipatt)*(tcol*SIN(tlong)-dy(1))/b(1)
      ELSE
         xx = tcol*COS(tlong)
         yy = tcol*SIN(tlong)
      END IF
      xcol  = SQRT(xx**2 + yy**2)
      xlat  = 90.0_rprec - xcol
      xl    = ATAN2(yy,xx)
      xlon  = xl+pi
      alpha = 2.0_rprec / (cmax-cmin)
      beta  = 1.0_rprec - alpha*cmax
      ct    = xlat*alpha+beta
      st    = SQRT(1.0_rprec - ct*ct)
      sph   = SIN(xlon)
      cph   = COS(xlon)
!
      IF (ABS(p(1,1)) < TINY (1.0_rprec)) THEN
         p(1,1)  = 1.0_rprec 
         dp(1,1) = 0.0_rprec
         sp(1)   = 0.0_rprec
         cp(1)   = 1.0_rprec
         DO n = 2, 18
!            fn(n) =n
            DO m = 1, n
!               fm(m) = m-1
               const(n,m) = REAL((n-2)**2-(m-1)**2, rprec)/&
                            REAL((2*n-3)*(2*n-5), rprec)
            END DO
         END DO
      END IF
!
      sp(2)  = sph
      p(1,1) = 1.0_rprec
      cp(2)  = cph
      do m = 3, nmax
         sp(m)=sp(2)*cp(m-1)+cp(2)*sp(m-1)
         cp(m)=cp(2)*cp(m-1)-sp(2)*sp(m-1)
      END DO
!
      value=coef(1,1)
!
      do n = 2, nmax
      do m = 1, n
         IF (n == m) THEN
            p(n,n)=st*p(n-1,n-1)
           dp(n,n)=st*dp(n-1,n-1)+ct*p(n-1,n-1)
         ELSE
            IF (n /= 2) p(n,m)=ct*p(n-1,m)-const(n,m)*p(n-2,m)
            IF (n == 2) p(n,m)=ct*p(n-1,m)
            dp(n,m)=ct*dp(n-1,m)-st*p(n-1,m)-const(n,m)*dp(n-2,m)
         END IF
      END DO
      END DO
!
!
!
      p(1,1) = p(1,1)*xco(1,1)
!
      value = coef(1,1)*p(1,1)
      DO n = 2, nmax
      DO m = 1, n
         IF (m /= 1) THEN
            pol      = p(n,m)*xco(n,m)
            p(m-1,n) = cp(m)*pol
            p(n,m)   = sp(m)*pol
            value    = value+p(m-1,n)*coef(m-1,n)+p(n,m)*coef(n,m)
         ELSE
            p(n,m) = p(n,m)*xco(n,m)
            value  = value+p(n,m)*coef(n,m)
         END IF
      END DO
      END DO
      IF (icntrl /= -1) THEN
         value = (pcp/(vmax(ipatt)-vmin(ipatt)))*value*1000.
      ELSE
         value = 1000.0_rprec*value
      END IF
!
!
      RETURN
      END SUBROUTINE Epot
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
!   6. A set of facilities to interpolate 1-D, 2-D arrays and 2-D sections
!      of 3-D arrays, real numbers:


    FUNCTION Gntrp_1d (array, i_max, bi )
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i_max
    REAL (rprec), INTENT (IN)    :: array (i_max), bi
    REAL (rprec)                 :: Gntrp_1d
!                                                                       
!   This function subprogram interpolates 1-dim. ARRAY to return
!   the value of ARRAY at the non-integer point BI.
!
!   Stanislav: if Bi < 1, then the array is extrapolated
!              linearly based on the values A(1,:) and A(2,:).
!              If Bi > imax, then array is linearly
!              extrapolated on the values A(imax-1,:) and A(imax,:).
!
    INTEGER (iprec) :: ii
    REAL (rprec)    :: fi
!
    ii = MAX (1, MIN (INT (bi), i_max-1))
    fi = REAL (ii,rprec)
    Gntrp_1d = (1.0_rprec - (bi-fi) ) * array (ii) + (bi-fi) * array (ii+1)
    RETURN
    END FUNCTION Gntrp_1d
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Gntrp_2d (array, i_max, j_max, bi, bbj)
    USE Rcm_variables, ONLY : iprec, rprec, n_gc
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i_max, j_max
    REAL (rprec), INTENT (IN)    :: array (1-n_gc:i_max+n_gc,1-n_gc:j_max+n_gc), bi, bbj
    REAL (rprec)                 :: Gntrp_2d
!
    INTEGER (iprec) :: ii, jn, jj, jp1
    REAL (rprec)    :: fi,fj,a1,a2,bj, &
                       Bjmod_real
!
!
    ii = MAX (1, MIN (INT (bi), i_max-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
    bj    = Bjmod_real ( bbj, j_max)
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp.
!
      Gntrp_2d = (1.0_rprec-(bi-fi)) * array(ii,jn) + &
                     (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interp.
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (1.0_rprec-(bi-fi))*array(ii,jj) + (bi-fi)*array(ii+1,jj)
       a2 = (1.0_rprec-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
       Gntrp_2d = (1.0_rprec - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Gntrp_2d
!
!
!
    FUNCTION Gntrp_2d_ang (array, i_max, j_max, jwrap, bi, bbj)
    USE Rcm_variables, ONLY : iprec, rprec, pi_two
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i_max, j_max, jwrap
    REAL (rprec), INTENT (IN)    :: array (i_max,j_max), bi, bbj
    REAL (rprec)                 :: Gntrp_2d_ang
!
    INTEGER (iprec) :: ii, jn, jj, jp1
    REAL (rprec)    :: fi,fj,a1,a2,bj, &
                       Bjmod_real
!
!
    ii = MAX (1, MIN (INT (bi), i_max-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
    bj    = Bjmod_real ( bbj,j_max)
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp.
!
      Gntrp_2d_ang = (1.0_rprec-(bi-fi)) * array(ii,jn) + &
                     (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interp.
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (1.0_rprec-(bi-fi))*array(ii,jj) + (bi-fi)*array(ii+1,jj)
       a2 = (1.0_rprec-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
          IF (jp1 == j_max) a2 = a2 + pi_two
          IF (jj  == jwrap)      a1 = 0.0_rprec
          IF (jp1 == jwrap)      a2 = a2 + pi_two      ! sts, feb 22
!
       Gntrp_2d_ang = (1.0_rprec - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Gntrp_2d_ang
!
!
!
    FUNCTION Intrp_2d_grid (array, bi, bbj, i_max, j_max, jwrap, imin_j)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: jwrap, i_max, j_max,imin_j(j_max)
    REAL (rprec), INTENT (IN)    :: array (i_max,j_max), bi, bbj
    REAL (rprec)                 :: Intrp_2d_grid
!
    INTEGER (iprec) :: ii,  jj
    REAL (rprec)    :: bj, ca, cb, cc, cd, di, dj, &
                       Bjmod_real
!
!
    ii = MAX (1, MIN (INT (bi), i_max-1))
    bj    = Bjmod_real ( bbj,j_max)
    jj  = INT (bj)
!
!   (i,j)---------------------(i,j+1)
!   | A          |               B  |
!   |            |                  |
!   |            |di                |
!   |   dj       |       1-dj       |
!   |---------(bi,bj)---------------|
!   |            |                  |
!   |            |                  |
!   |            |1-di              |
!   |            |                  |
!   | C          |              D   |
!   (i+1,j)-----------------(i+1,j+1)
!
    di = bi - ii
    dj = bj - jj
    IF (ii < imin_j(jj)) THEN
       ca = 0.0
    ELSE
       ca = (1.0-di)*(1.0-dj) 
    END IF
    IF (ii+1 < imin_j(jj)) THEN
       cc = 0.0
    ELSE
       cc = di * (1.0-dj)
    END IF
    IF (ii < imin_j(jj+1)) THEN
       cb = 0.0
    ELSE
       cb = (1.0-di)*dj
    END IF
    IF (ii+1 < imin_j(jj+1)) THEN
       cd = 0.0
    ELSE
       cd = di*dj
    END IF
!                                                                       
    Intrp_2d_grid = ca*array(ii,jj)+cb*array(ii,jj+1)+ &
                      cc*array(ii+1,jj)+cd*array(ii+1,jj+1)
    Intrp_2d_grid = Intrp_2d_grid / (ca+cb+cc+cd)
!
    RETURN
    END FUNCTION Intrp_2d_grid
!
!
!
!
!
!   6.2 A set of facilities to interpolate 1-D, 2-D arrays and 2-D sections
!      of 3-D arrays, real numbers, but these arrays are not circularized in J:


    FUNCTION Interp_1d (array, bi )
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:), bi
    REAL (rprec)                 :: Interp_1d
!                                                                       
!   This function subprogram interpolates 1-dim. ARRAY to return
!   the value of ARRAY at the non-integer point BI.
!
!   Stanislav: if Bi < 1, then the array is extrapolated
!              linearly based on the values A(1,:) and A(2,:).
!              If Bi > imax, then array is linearly
!              extrapolated on the values A(imax-1,:) and A(imax,:).
!
    INTEGER (iprec) :: ii, imax_array
    REAL (rprec)    :: fi
!
    imax_array = SIZE (array)
!
    IF (bi < 1.0_rprec .OR. bi > imax_array) &
         call CON_STOP('OUT OF BOUNDS IN INTERP_1D')
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
    Interp_1d = (1.0_rprec - (bi-fi) ) * array (ii) + (bi-fi) * array (ii+1)
    RETURN
    END FUNCTION Interp_1d



    FUNCTION Interp_2d (array, bi, bj)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bj
    REAL (rprec)                 :: Interp_2d
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
!
!
    IF (bi < 1.0_rprec .OR. bi > imax_array .OR. bj < 1 .OR. bj > jmax_array) &
        call CON_STOP('OUT OF BOUNDS IN INTERP_2D')
!
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d = (1.0_rprec-(bi-fi)) * array(ii,jn) + (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (1.0_rprec-(bi-fi))*array(ii,jj)  + (bi-fi)*array(ii+1,jj)
       a2 = (1.0_rprec-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
       Interp_2d = (1.0_rprec - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d




    FUNCTION Interp_2d_of3d (array, bi, bj, index_3)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: index_3
    REAL (rprec), INTENT (IN)    :: array (:,:,:), bi, bj
    REAL (rprec)                 :: Interp_2d_of3d
!                                                                       
!   This is the same as Gntrp_2d but for a 3-dim array, see comments for
!   Gntrp_2d. A separate function is needed since if Gntrp_2d were used,
!   then we would need to pass array sections (the other option is Fortran
!   77 style of passing an offset array, but that should be avoided for
!   compiler checking and parallelization reasons).
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array, kmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
    kmax_array = SIZE (array, DIM = 3)
!
!
    IF (bi < 1.0_rprec .OR. bi > imax_array .OR. bj < 1.0_rprec .OR. bj > jmax_array) THEN
       WRITE (*,*) 'OUT OF BOUNDS IN INTERP_2D_OF_3D'
       call CON_stop('ERROR in IM/RCM2/src/rcm_routines.f90')
    END IF
!
!
    IF (index_3 > kmax_array .OR. index_3 < 1) &
         call CON_STOP('GNTRP_2D_OF3D: index_3 OUT OF RANGE')
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d_of3d = (1.0_rprec-(bi-fi)) * array(ii,jn,index_3) + &
                  (bi-fi)*array(ii+1,jn,index_3)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (1.0_rprec-(bi-fi))*array(ii,jj,index_3)  + (bi-fi)*array(ii+1,jj,index_3)
       a2 = (1.0_rprec-(bi-fi))*array(ii,jp1,index_3) + (bi-fi)*array(ii+1,jp1,index_3)
!                                                                       
       Interp_2d_of3d = (1.0_rprec - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d_of3d
!
!
!
!   7. A few routines needed for circulariation and normalization in J:
!
    FUNCTION Bjmod_real (bj, jsize)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    REAL (rprec),    INTENT (IN) :: bj
    INTEGER(iprec),  INTENT (IN) :: jsize
    REAL (rprec)                 :: Bjmod_real
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_real = bj
!                                                                       
    do_1: DO
       IF (Bjmod_real > REAL (jsize,rprec)) THEN
          Bjmod_real = Bjmod_real - REAL (jsize,rprec)
       ELSE
          EXIT do_1
       END IF
    END DO do_1
!
    do_2: DO
       IF (Bjmod_real < 1 ) THEN
           Bjmod_real = Bjmod_real + REAL(jsize,rprec)
       ELSE
           EXIT do_2
       END IF
    END DO do_2
!
    RETURN
    END FUNCTION Bjmod_real
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Bjmod_int (bj, jsize)
    USE Rcm_variables, ONLY : iprec, rprec, n_gc
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: bj
    INTEGER (iprec), INTENT (IN) :: jsize
    INTEGER (iprec)              :: Bjmod_int
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_int = bj
!                                                                       
    do_1: DO
       IF (Bjmod_int > jsize) THEN
          Bjmod_int = Bjmod_int - jsize 
       ELSE
          EXIT do_1
       END IF
    END DO do_1
!
    do_2: DO
       IF (Bjmod_int < 1) THEN
           Bjmod_int = Bjmod_int + jsize
       ELSE
           EXIT do_2
       END IF
    END DO do_2
!
    RETURN
    END FUNCTION Bjmod_int
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
    FUNCTION Bndy (bndloc, jsize, bj)
    USE Rcm_variables, ONLY : iprec, rprec, n_gc
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: jsize
    REAL (rprec), INTENT (IN) :: bndloc(1-n_gc:jsize+n_gc), bj
    REAL (rprec)              :: Bndy
!
!   Function to determine whether you are outside RCM
!   boundary as defined by array AIN; returns boundary 
!   location (non-integer I-value) for given bj
!   Written 1/25/96 frt                                   
!                                                                       
!   Dependency: BJMOD
!
    REAL (rprec)    :: bip, bim, bj_point, bjp, bjm, &
                       Bjmod_real
    INTEGER (iprec) :: jm, jp
!                                                                       
!   Call to BJMOD returns bj_point in [jwrap,jsize), then jm is in
!   [jwrap,jsize) and jp1 is in [jwrap+1,jsize], so all three indices are
!   within the required range.
!
    bj_point = Bjmod_real (bj, jsize)
    jm       = INT (bj_point)
    jp       = jm + 1
    bjp = REAL (jp,rprec)
    bjm = REAL (jm,rprec)
!                                                                       
    bim = bndloc (jm)
    bip = bndloc (jp)
!                                                                       
    bndy = bim * (bjp - bj_point) + bip * (bj_point - bjm)
!                                                                       
    RETURN
    END FUNCTION Bndy
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
    SUBROUTINE Adjust_bj_2 (b1_in, b2_in, b1_out, b2_out )
    USE Rcm_variables, ONLY : iprec, rprec, jsize
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)  :: b1_in,  b2_in
    REAL (rprec), INTENT (OUT) :: b1_out, b2_out
!                                                                       
!   Author: r.w. spiro
!   Last update: 1-3-84 by rws
!   This subroutine adjusts a pair of input bj values (b1 and b2)
!   to the same modulus.  b3 and b4 are the returned values with
!   b3 corresponding to b1, b4 to b2.
!   Here we assume that both B1 and B2 are normalized by BJMOD so that
!   they are in the interval [jwrap, jmax-1)
!   Dependency:  none
!
    REAL (rprec) :: dist_half, dist_full
!
!
    dist_full = REAL (jsize,rprec)
    dist_half = REAL (jsize,rprec) * 0.5_rprec
    IF (ABS (b1_in - b2_in) <= dist_half ) THEN
      b1_out = b1_in
      b2_out = b2_in
    ELSE IF (b1_in < b2_in) THEN
      b1_out = b1_in + dist_full
      b2_out = b2_in
    ELSE
      b1_out = b1_in
      b2_out = b2_in + dist_full
    END IF
!
    RETURN
    END SUBROUTINE Adjust_bj_2
!
!
!
!
    SUBROUTINE Adjust_bj_3 (b1_in, b2_in, b3_in, b1_out, b2_out, b3_out)
    USE Rcm_variables, ONLY : iprec, rprec, jsize
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: b1_in, b2_in, b3_in
    REAL (rprec), INTENT (OUT):: b1_out, b2_out, b3_out
!
!   Last update: 11-29-84          by:rws
!                                                                       
!   This subroutine adjusts a triad of input bj values (bj1, bj2, bj3)
!   to the same modulus.  b1, b2, and b3 are the returned values.
!
!   Dependency:  none
!
    LOGICAL :: ic1, ic2, ic3
    REAL (rprec)    :: diftst, dist_full
!
    diftst    = 0.5_rprec * REAL (jsize,rprec)
    dist_full = REAL (jsize,rprec)
!
    ic1 = (ABS (b1_in - b2_in) < diftst)
    ic2 = (ABS (b1_in - b3_in) < diftst)
    ic3 = (ABS (b2_in - b3_in) < diftst)
!                                                                       
    IF (ic1 .AND. ic2 .AND. ic3 ) THEN
       b1_out = b1_in
       b2_out = b2_in
       b3_out = b3_in
    ELSE IF (ic3) THEN
       b2_out = b2_in
       b3_out = b3_in
       IF (b1_in < b2_in) THEN
          b1_out = b1_in + dist_full
       ELSE
          b1_out = b1_in - dist_full
       END IF
    ELSE IF (ic2) THEN
       b1_out = b1_in
       b3_out = b3_in
       IF (b2_in < b1_in) then
          b2_out = b2_in + dist_full
       ELSE
          b2_out = b2_in - dist_full
       END IF
    ELSE
       b1_out = b1_in
       b2_out = b2_in
       IF (b3_in < b1_in) THEN
          b3_out = b3_in + dist_full
       ELSE
          b3_out = b3_in - dist_full
       END IF
    END IF
    RETURN
    END SUBROUTINE Adjust_bj_3
!
!
!   10. Facilities for inner edges cleanup (addition and deletion of points).
!       There are only 3 "front-end" routines, the rest are called from these
!       three. The main routines are to delete extra points (zap_edge),
!       routine to add needed points (expand_edge),and routine to chop off long
!       thing tails (farmers_wife).
!
!
    SUBROUTINE Farmers_wife (k)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (KIND=iprec), INTENT (IN) :: k
!
!
!   last update: 03-31-88              by: rws
!               09-19-97                  rws                           
!                                     changed algorithm used            
!   Sat May 22 11:09:53 MDT 1999  corrected a bug (y2 vs ym)
!                                                                       
!   this subroutine marches around inner edge k searching 
!   for long thin tails to eliminate.  control%fmrwif_dlim is a control
!   parameter that gives the minimum allowed tail thickness (RE).
!
!   Algorithm: go around the edge, consider each point "m" and its two
!   neighbors (mnext and mlast). There are three possible logical paths:
!   (1) if the three points do not form a "tail" (no extremum in BI or BJ
!   for the middle point "m"), then go to the next point and repeat the
!   test. (2) if there is an extremum but the 3-point tail is not too thin
!   (thickness is computed in the equatorial plane in the units of Re, see
!   the code for details), skip to the next point too. In both cases, we
!   count (mcount=mcount+1) how many points (triads) we tested, and if we
!   checked all points, exit the routine. (3) If the tail is too thin, then
!   we delete the middle point. In this case, we go back one point, to
!   "mlast", because its "next" neighbor is different now, and restart testing
!   from this point. We also reset "mcount" in this case, which allows us to
!   test the whole edge again. In other words, the looping over particles is
!   exited only after we have gone around the edge completely without elimi-
!   nating anything.
!
!   In the code, rprev, rthis, and rnext are two-dim vectors holding the
!   (x,y) locations of "mlast", "m", and "mnext" points respectively in the
!   equatorial GSM plane.
!
!   I don't know what to do if a held point must be deleted, so put a STOP
!   statement. For the future, I guess just don't delete it?
!   Second issue: doesn't consider if points are outside bndy; but hopefully
!   tails don't form there.
!
!   Dependency:  bmod3, compar, gntrp, elim
!
    INTEGER (iprec) :: mbeg, mend, m, mcount, mnext, mlast, is
    REAL (rprec) :: b1, b2, b3, d_to_last, d_to_next, f, rlast (2),&
                        rnext (2), rthis (2), &
                    Gntrp_2d
    LOGICAL :: Extremum
!
    mbeg   = mpoint (k)
    m      = mbeg
    mcount = 0
!                                                                       
    DO
       mend = mbeg + npoint (k) - 1
!
       IF (mcount >= npoint (k) ) EXIT   !   normal exit from loop
!
       mnext = m + 1
       mlast = m - 1
       IF (m == mend) mnext = mbeg
       IF (m == mbeg) mlast = mend
!                                                                       
       CALL Adjust_bj_3 (bj(mlast), bj(m), bj(mnext), b1, b2, b3 )
!                                                                       
       IF ( Extremum(b1,b2,b3) .OR. Extremum(bi(mlast), bi(m), bi(mnext))) THEN
!
!         calculate distance from m to m-1 and from m to m+1
!                                                                       
          rthis(1) = Gntrp_2d (xmin, isize, jsize, bi(m), bj(m))
          rthis(2) = Gntrp_2d (ymin, isize, jsize, bi(m), bj(m))
!                                                                       
          rlast(1) = Gntrp_2d (xmin, isize, jsize, bi(mlast), bj(mlast))
          rlast(2) = Gntrp_2d (ymin, isize, jsize, bi(mlast), bj(mlast))
!                                                                       
          rnext(1) = Gntrp_2d (xmin, isize, jsize, bi(mnext), bj(mnext))
          rnext(2) = Gntrp_2d (ymin, isize, jsize, bi(mnext), bj(mnext))
!
          d_to_last = SQRT (DOT_PRODUCT (rthis-rlast, rthis-rlast ) )
          d_to_next = SQRT (DOT_PRODUCT (rthis-rnext, rthis-rnext ) )
!
          IF (d_to_last <= d_to_next) THEN
!
!            Adjust position of mnext to be along line from "m" to "mnext",
!            but same distance as from "m" to "mlast":
!
             f = d_to_last / d_to_next
             rnext = f * rnext + (1.0_rprec - f) * rthis
          ELSE
             f = d_to_next / d_to_last
             rlast = f * rlast + (1.0_rprec - f) * rthis
          END IF
!                                                                       
!
!         Check if the middle point is held:
          CALL Dcode (itrack, SIZE(itrack), m, is)
!
          IF (SQRT (DOT_PRODUCT(rnext-rlast,rnext-rlast)) < fmrwif_dlim) THEN
!
!            "mnext" is too close to "mlast", i.e, tail is too thin;
!            therefore, eliminate the leading point m:
!
             IF (is < 0) call CON_STOP( &
                  'FARMERS WIFE MUST DELETE A POINT, BUT IT IS HELD')
!
             CALL Delete_point (k, m)
!                                                                       
!            must consider previous pt for elimination
!
             m = m - 1
             IF (m < mbeg) m = mbeg + npoint (k) - 1
             mcount = 0
          ELSE
             mcount = mcount + 1
             m      = m + 1
             IF (m > mbeg + npoint(k) - 1)  m = mbeg
          END IF
       ELSE
          mcount = mcount + 1
          m      = m + 1
          IF (m > mbeg + npoint(k) - 1)  m = mbeg
       END IF
    END DO
!
    RETURN
    END SUBROUTINE Farmers_wife
!
!
!
!
    SUBROUTINE Zap_edge (k)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER(KIND=iprec), INTENT (IN)  :: k
!
!   last update: 03-31-88              by: rws
!                01-29-96   added calls to bndy and ain arrays - frt
!                                                                       
!   this subroutine marches around an inner edge deciding whether to    
!   delete(zap) pts that are too close together.
!   point m+1 is zapped if
!   (bi(m+1)-bi(m))**2 + (bj(m+1)-bj(m))**2 < dstmin**2
!                                                                       
!     where dstmin= min spacing parameter and neither point m+2
!     nor point m are being held at the boundary, and point m+1
!     is not an extremum.
!                                                                       
!   Elimination stategy: if m+1 is too close to m, then
!   m+1 will be eliminated only if:
!       1. Neither m nor m+2 (two neighbors) are held,
!       2. m+1 has no extremum in either J or I (except when m+1 is
!         (outside the modeling region),
!       3. if m+1 is outside the modeling region, then m and m+2
!          (both neighbors) must be outside too.
!   To determine whether a point is being held, DCODE is called to
!   check for negative sign of itrack(m).
!
!   Stanislav: add checking whether "m+1" itself held; if yes, don't
!   eliminate no matter what. Oct 15 1999.
!
!   Dependency:  bndy, dcode, bmod, bmod3, compar, elim
!
    REAL (KIND=rprec) :: b1, b2, b3, &
                         Bndy
    INTEGER (KIND=iprec) :: m, mp1, mp2, mend, is
    LOGICAL :: Extremum
!
!
    IF (dstmin > 0.5_rprec*dstmax ) &
         call CON_STOP('trouble with sqmin or sqmax - stopping')
!                                                                       
    IF (npoint(k) <= 2) &
         call CON_STOP('ZAP_EDGE: THE EDGE IS TWO POINTS OR LESS IN LENGTH')
!
    m = mpoint(k)
    Loop_10: DO
       mend = mpoint (k) + npoint (k) - 1
       IF (m > mend) EXIT    ! normal exit from loop
!
       mp1 = m + 1
       IF (m == mend) mp1 = mpoint (k)
       mp2 = mp1 + 1
       IF (mp1 == mend) mp2 = mpoint (k)
!
!      Point "mp1" will be considered for elimination now. It will be
!      eliminated if it is too close to "m"; but some cases must be
!      excluded. The criteria depend on whether the point is with respect
!      to the boundary.
!
!      First, in any case it cannot be eliminated if any one of the two
!      neighbors is held:
!
       CALL Dcode (itrack, SIZE(itrack), mp1, is)
       IF (is < 0) THEN
          m = m + 1
          CYCLE Loop_10
       END IF
       CALL Dcode (itrack, SIZE(itrack), m, is)
       IF (is < 0) THEN
          m = m + 1
          CYCLE Loop_10
       END IF
       CALL Dcode (itrack, SIZE(itrack), mp2, is)
       IF (is < 0) THEN
          m = m + 1
          CYCLE Loop_10
       END IF
!
!
!      Second, is the point inside or outside the boundary?
!
       IF (bj(mp1) <= Bndy (bndloc, jsize, bi(mp1))) THEN
!
!         Point "mp1" is outside. It can only be considered for elimination
!         if both neighbors are outside too (don't care about extrema since
!         edges outside the boundary don't generate any currents anyway).
!         Make sure this is the case; if not, skip considering "mp1" altogether:
!
          IF (bj(m) > Bndy (bndloc,  jsize, bi(m)) .OR. &
                            bj(mp2) > Bndy (bndloc, jsize, bi(mp2))) THEN
             m = m + 1
             CYCLE Loop_10
          END IF
       ELSE
!
!         Point "mp1" is inside. Before considering it for elimination, we
!         must make sure that there is no extremum in either I or J;
!         otherwise, just skip to consider the next point:
!
          CALL Adjust_bj_3 (bj(m), bj(mp1), bj(mp2), b1, b2, b3 )
          IF (Extremum(bi(m), bi(mp1), bi(mp2)) .OR. Extremum (b1, b2, b3)) THEN
             m = m + 1
             CYCLE Loop_10
          END IF
       END IF
!
!
!      Now it is safe to consider elimination:
!
       CALL Adjust_bj_2 (bj(mp1), bj(m), b1, b2)
       IF ( SQRT( (bi(mp1)-bi(m))**2 + (b1-b2)**2) < dstmin) THEN
          CALL Delete_point (k, m)
          IF (mp1 < m) THEN
!            We have just deleted the first point of the edge; therefore,
!            the whole edge was shifted left by 1, so the current value
!            of m is now m-1:
             m = m - 1
          END IF
          CYCLE loop_10
       ELSE
          m = m + 1
          CYCLE Loop_10
       END IF
    END DO loop_10
!
    RETURN
    END SUBROUTINE Zap_edge
!
!
!
!
    SUBROUTINE Delete_point (k, mzap)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (KIND=iprec), INTENT (IN) :: k, mzap
!
!   last update: 03-31-88            by:rws
!   this subroutine eliminates pt at mzap by shifting all points with
!   numbers higher than "mzap" by one to the left, squeezing out mzap.
!
!   This routine does not check if a point is being held, you must do this
!   in the calling routine explicitely.
!
!   Dependency:  none
!
    INTEGER (KIND=iprec) :: mend, m, kk, kdim
!
    kdim = SIZE (alam)
    mend = mpoint (kdim) + npoint (kdim) - 1
!
!   1. Shift all points after "mzap" to the "left" by one, squeezing out mzap:
    DO  m = mzap, mend
       bi (m) = bi (m+1)
       bj (m) = bj (m+1)
       etab (m) = etab (m+1)
       itrack(m) = itrack (m+1)
    END DO
!
!   2. The number of points on edge "k" is less by one:
    npoint (k) = npoint (k) - 1
!
!   3. For all edges starting after mzap, starting points shifted:
    DO kk = 1, kdim
       IF (npoint (kk) /= 0 .AND. mpoint(kk) > mzap) mpoint(kk) = mpoint(kk)-1
    END DO
!                                                                       
    mpoint (kdim + 1) = mpoint (kdim + 1) - 1
!                                                                       
    RETURN
    END SUBROUTINE Delete_point
!
!
!
!
!   Miscelaneous routines needed by others, such as computing derivatives,
!   encoding/decoding particles' status, checking for extrema, etc.
!
    SUBROUTINE Dcode (itrack, nptmax, m, is)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: nptmax, itrack (nptmax), m
    INTEGER (iprec), INTENT (OUT) :: is
!
!   Last update:  05-16-85                by: rws
!   This subroutine decodes the identification information
!   contained in itrack(m).
!   If a point "m" is held, itrack (m) is negative, else positive.
!
!   Dependency:  none
!
    IF (itrack (m) < 0) THEN
       is = - 1
    ELSE
       is = 1
    END IF
    RETURN
    END SUBROUTINE Dcode
!
!
!
!
    SUBROUTINE Ncode (itrack, nptmax, is, newitk)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (KIND=iprec), INTENT (IN)     :: is, nptmax
    INTEGER (KIND=iprec), INTENT (IN OUT) :: itrack(nptmax)
    INTEGER (KIND=iprec), INTENT (OUT)    :: newitk
!
!   last update: 05-24-85             by: rws
!                                                                       
!   function subprogram to assign a unique identifier to an
!   inner edge tracer pt.
!                                                                       
!   ncode.lt.0 means pt is being held
!   June 3, 1999--Stanislav: changed function to subroutine.
!   what used to be 'NCODE' is now 'NEWITK'
!
!   Dependency:  none
!
    newitk = is * itrack (nptmax)
    itrack (nptmax) = itrack (nptmax) + 1
    RETURN
    END SUBROUTINE Ncode
!
!
!
    FUNCTION Extremum (b1, b2, b3)
    USE Rcm_variables, ONLY : iprec, rprec
    IMPLICIT NONE
    REAL (KIND=rprec), INTENT (IN) :: b1, b2, b3
    LOGICAL :: extremum
!
!   last update:  05-16-85           by: rws
!                                                                       
!   this subroutine tests the triad b1,b2,b3 for an
!   extremum at b2.
!   Dependency:  none
!
    IF (b2 > b1 .AND. b2 > b3) THEN  ! B2 IS A MAXIMUM
       extremum = .TRUE.
    ELSE IF (b2 < b1 .AND. b2 < b3) THEN     !  b2 is a minimum
       extremum = .TRUE.
    ELSE
       extremum = .FALSE.
    END IF
    RETURN
    END FUNCTION Extremum
!
!
    SUBROUTINE Outbij (k_init, k_final, NCOL, ntp)
      USE Rcm_variables
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: k_init, k_final, ncol, ntp
!
!
    INTEGER (iprec) :: ikind, k, mbeg, mend, mfinal, mk, m, kdim
    REAL(rprec) :: vmtrp (5), &
                   Gntrp_2d
!                                                                       
    ikind = 0
    kdim  = SIZE (alam)
    IF (k_init < 1 .OR. k_final > kdim) THEN
       WRITE (*,'(T2,A)') 'OUTBIJ: K VALUES ARE OUT OF RANGE, STOPPING'
       call CON_stop('ERROR in IM/RCM2/src/rcm_routines.f90')
    END IF
!
    DO k = k_init, k_final
!
       WRITE (ntp, 810) k, alam(k), eta(k), mpoint(k), npoint(k)
!
       mbeg = mpoint (k)
       mend = mbeg + npoint (k) - 1
       DO
          mfinal = ncol/40 - 1 + mbeg
          IF (mfinal > mend) mfinal = mend
          mk = 0
          DO m = mbeg, mfinal
             mk = mk + 1
             vmtrp (mk) = Gntrp_2d (vm, isize, jsize, bi(m), bj(m))
          END DO
          WRITE (ntp, 820) (itrack(m), bi(m), bj(m), etab(m), &
                            vmtrp (m-mbeg+1), m = mbeg, mfinal )
          IF (mfinal == mend) EXIT
          mbeg = mfinal + 1
       END DO
!
    END DO
!
 810  FORMAT (/ T02, 'K=', I3.3, TR2, 'ALAM=', ES9.2, TR2, 'ETA=', ES10.2,&
                TR2, 'MPOINT=', I5, TR2, 'NPOINT=', I5)
 820  FORMAT (3(TR2, I6.6, '(', TR1, F5.2, ',', TR1, F5.2, ',', &
              E10.2, ',', F5.2,')') )
!
    RETURN
    END SUBROUTINE Outbij
!
     FUNCTION Lt_from_aloct (phi)
      USE Rcm_variables
     IMPLICIT NONE
     REAL (rprec), INTENT (IN) :: phi
     REAL (rprec) :: Lt_from_aloct
!
!    Convert an RCM phi angle (aloct in ionosphere) to MLT
!    Output (result) is:   0.0  <=  result < 24.00
!
     IF (phi < 0.0_rprec .OR. phi > pi_two) THEN
        WRITE (*,*) 'IN LT_FROM_ALOCT, PHI IS OUT OF BOUNDS'
        call CON_stop('ERROR in IM/RCM2/src/rcm_routines.f90')
     ELSE
        Lt_from_aloct = MODULO((phi-pi)*RTH, 24.0_rprec)
        Lt_from_aloct = MODULO(Lt_from_aloct, 24.0_rprec)
     END IF
     RETURN
     END FUNCTION Lt_from_aloct
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Dipole_Bfield (theta, phi, arg)
      USE Rcm_variables, ONLY : iprec, rprec, besu
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: theta, phi
      INTEGER (iprec), INTENT (IN) :: arg
      REAL (rprec) :: Dipole_Bfield
!
!       THETA and PHI are in radians, THETA measured from n.pole down,
!       PHI measured from noon to dusk to midnight etc.
!       Distances RMIN, XMIN, and YMIN are in units of RE
!       Since in the RCM, BESU is in [nT], and the factor of RE is ommited
!       from the formula, VM has units of (RE/nT)**(-2/3)
!
      REAL (rprec) :: rmin, xmin, ymin, bmin, vm
!
      rmin = 1.0_rprec / SIN(theta)**2
      xmin = rmin * COS (phi)
      ymin = rmin * SIN (phi)
      bmin = besu * (1.0_rprec/rmin**3)
      vm   = (32.0_rprec/35.0_rprec* rmin**4 / Besu *  &
              SQRT(1.0_rprec-1.0_rprec/rmin)* &
              (1.0_rprec+0.5_rprec/rmin + &
               3.0_rprec/8.0_rprec/rmin**2 + &
               5.0_rprec/8.0_rprec/2.0_rprec/rmin**3) &
              ) ** (-2.0_rprec/3.0_rprec)
!
      IF (arg == 1) THEN
         Dipole_Bfield = xmin
      ELSE IF (arg == 2) THEN
         Dipole_Bfield = ymin
      ELSE IF (arg == 3) THEN
         Dipole_Bfield = bmin
      ELSE IF (arg == 4) THEN
         Dipole_Bfield = vm
      ELSE
         call CON_STOP('ILLEGAL ARGUMENT FOR DIPOLE_BFIELD')
      END IF
!
      RETURN
      END FUNCTION Dipole_Bfield
!
!
!
!
!
      SUBROUTINE Deriv_i (a, i_size, j_size, imin_j, d_di)
      USE Rcm_variables, ONLY : iprec, rprec, n_gc
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: i_size, j_size, imin_j(1-n_gc:j_size+n_gc)
      REAL (rprec), INTENT (IN) :: a    (1-n_gc:i_size+n_gc,1-n_gc:j_size+n_gc)
      REAL (rprec), INTENT (OUT) ::d_di (1-n_gc:i_size+n_gc,1-n_gc:j_size+n_gc)
!
!_________________________________________________________________________
!     The idea is to use central differences for a second-order accuracy
!     of the first-order derivative, if we can. If one of the neighboring
!     points is outside the boundary, use forward or back difference
!     (first-order accuracy) if we can. If both points are outside, set
!     derivative to zero.
!_________________________________________________________________________
!
      INTEGER (iprec) :: i, j
!
!
      DO j = 1, j_size
!
         d_di (1:imin_j(j)-1,j) = 0.0_rprec
!
         i = imin_j(j)
         d_di (i,j) = -1.5*a(i,j) + 2.0*a(i+1,j) - 0.5*a(i+2,j)
!
         DO i = imin_j(j)+1, i_size - 1
            d_di (i,j) = 0.5_rprec*(a(i+1,j) - a(i-1,j))
         END DO
!
         i = i_size
         d_di (i,j) = 1.5*a(i,j) - 2.0*a(i-1,j) + 0.5*a(i-2,j)
!
      END DO
      CALL Wrap_around_ghostcells (d_di, i_size, j_size, n_gc)
!
      RETURN
      END SUBROUTINE Deriv_i
!
!
!
      SUBROUTINE Deriv_j (a, i_size, j_size, imin_j, error_value, d_dj)
      USE Rcm_variables, ONLY : iprec, rprec, n_gc
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: i_size, j_size, imin_j(1-n_gc:j_size+n_gc)
      REAL (rprec), INTENT (IN) :: a (1-n_gc:i_size+n_gc,1-n_gc:j_size+n_gc), error_value
      REAL (rprec), INTENT (OUT) :: d_dj (1-n_gc:i_size+n_gc,1-n_gc:j_size+n_gc)
!___________________________________________________________________________
!
!     Take derivative of a with respect to J. The idea is to use
!     grid points that are inside the modeling region. Stanislav 10/27/2000.
!     3/6/2001: modified to be the same as Frank's version in rcm296.
!     Notice that pt_jpp is defined differently (>) than pt_jp(>=), and ditto
!     for pt_jmm and pt_jm. Took this from Frank's version.
!___________________________________________________________________________
!
      INTEGER (iprec) :: i_bnd, i, j, jm, jmm, jp, jpp
      LOGICAL :: pt_jp, pt_jpp, pt_jm, pt_jmm
!
!
      i_bnd = MAXVAL (imin_j)
!
      DO j = 1, j_size
!
         jp = j + 1
         jpp= j + 2
         jm = j - 1
         jmm= j - 2
!
         d_dj (1:imin_j(j)-1,j) = 0.0
!
         DO i = imin_j(j), i_bnd-1
!
            pt_jp  = (i >= imin_j(jp))
            pt_jm  = (i >= imin_j(jm))
            pt_jpp = (i > imin_j(jpp))
            pt_jmm = (i > imin_j(jmm))
!
            IF (pt_jp .AND. pt_jm) THEN
               d_dj (i,j) = 0.5_rprec*(a (i,jp) - a (i,jm))
            ELSE IF (pt_jp) THEN !j-1 is outside
               IF (pt_jpp) THEN
                  d_dj (i,j)= -1.5*a(i,j)+2.0*a(i,jp)-0.5*a(i,jpp)
               ELSE
                  d_dj (i,j) = a(i,jp) - a(i,j)
               END IF
            ELSE IF (pt_jm) THEN  ! j+1 is outside
               IF (pt_jmm) THEN
                  d_dj (i,j) = +1.5*a(i,j)-2.0*a(i,jm)+0.5*a(i,jmm)
               ELSE
                  d_dj (i,j) = (a (i,j) - a (i,jm))
               END IF
            ELSE
               d_dj (i,j) = error_value
!!$print*,i,j,'deriv_j'
            END IF 
!
         END DO
!
         d_dj (i_bnd:,j) = 0.5_rprec*( a(i_bnd:,jp) - a(i_bnd:,jm))
      END DO
!
      CALL Wrap_around_ghostcells (d_dj, i_size, j_size, n_gc)
!
      RETURN
      END SUBROUTINE Deriv_j
!
      SUBROUTINE Wrap_around_ghostcells (array, isize, jsize, n_gc)
      USE rcm_variables, only : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n_gc, isize, jsize
      REAL (rprec), INTENT (IN OUT) :: array (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc)
!
      array (1:isize, -1) = array (1:isize, jsize - 1)
      array (1:isize, 0)  = array (1:isize, jsize )
      array (1:isize, jsize+1) = array (1:isize, 1)
      array (1:isize, jsize+2) = array (1:isize, 2)
!
      array (-1, :) = array (1, :)
      array (0,  :) = array (1, :)
      array (isize+1, :) = array (isize,:)
      array (isize+2, :) = array (isize,:)
!
      RETURN
      END SUBROUTINE Wrap_around_ghostcells
!
!
!
      SUBROUTINE Wrap_around_ghostcells_log (array, isize, jsize, n_gc)
      USE rcm_variables, only : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n_gc, isize, jsize
      LOGICAL, INTENT (IN OUT) :: array (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc)
!
      array (1:isize, -1) = array (1:isize, jsize - 1)
      array (1:isize, 0)  = array (1:isize, jsize )
      array (1:isize, jsize+1) = array (1:isize, 1)
      array (1:isize, jsize+2) = array (1:isize, 2)
!
      array (-1, :) = array (1, :)
      array (0,  :) = array (1, :)
      array (isize+1, :) = array (isize,:)
      array (isize+2, :) = array (isize,:)
!
      RETURN
      END SUBROUTINE Wrap_around_ghostcells_log
!
!
SUBROUTINE Move_plasma_grid_NEW (n_gc,isize, jsize, kcsize, iesize, i1, i2, j1, j2, imin_j, &
                                 dt, vcorot, vpar, v, &
                                 alamc, etac, eeta, vm, xmass, fudgec, ikflavc, sini,&
                                 alpha, beta, bir, dlam, dpsi, signbe, Ri, ie_ask )
  USE Rcm_variables, ONLY : rprec, sunspot_number, doy, rmin, L_dktime, xmin, ymin, &
                            irdk,inrgdk,isodk,iondk, dktime, iTimeT1
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: isize, jsize, kcsize, i1, i2, j1, j2, iesize,n_gc
  INTEGER, INTENT (IN) :: imin_j(1-n_gc:jsize+n_gc), ikflavc (kcsize)
  REAL(rprec), INTENT (IN), DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc) :: &
                       vcorot, vpar, v, vm, alpha, beta, bir,sini
  REAL(rprec), INTENT (IN) :: dlam, dpsi, signbe, ri, alamc (kcsize), &
                       etac (kcsize), xmass (iesize), fudgec(kcsize)
  REAL(rprec), INTENT (IN OUT) :: eeta (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,kcsize)  
  REAL(rprec), INTENT (IN) :: dt
  INTEGER, INTENT (IN) :: ie_ask

  REAL (rprec) :: Plasmasphere_refill_rate, Cexrat

!_____________________________________________________________________________
!   Subroutine to advance eta distribution for a time step
!   by using new CLAWPACK advection routines
!                                                                       
!   Created:     12-05-00
!_____________________________________________________________________________
!
!
  REAL(rprec) ::         veff     (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  dvefdi   (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  dvefdj   (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  loc_didt (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  loc_djdt (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  loc_Eta  (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  loc_rate (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  fac      (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), &
                  didt, djdt, rate, mass_factor, r_dist
  INTEGER :: i, j, kc, ie
!
  INTEGER :: CLAWiter, icut
  DOUBLE PRECISION :: xlower,xupper,ylower,yupper, T1,T2
  LOGICAL, SAVE :: FirstTime=.true.
! 
! 
! 
  T1=iTimeT1
!   IF (FirstTime) THEN
!      T1=0.
!   ELSE
!      T1=T2
!   END IF
  T2=T1+dt
!
!
!
  DO kc = 1, kcsize
!
do i = 1, isize
do j = 1, jsize
  if (eeta(i,j,kc) < 0.) then
  !  print*,'eeta error_i: ',i,imin_j(j),j,kc,eeta(i,j,kc)
!!$     print*,'kc=',kc,' i=',i,' eeta(i,:,kc)='
!!$     print*,eeta(i,:,kc)
!!$     print*,'kc=',kc,' j=',j,' eeta(:,j,kc)='
!!$     print*,eeta(:,j,kc)
!!$     call CON_STOP('abort')
!!$     eeta(i,j,kc) = 0.0
     eeta(i,j,kc) = 0.5*(eeta(max(1,i-1),j,kc)+eeta(min(isize,i+1),j,kc))
  end if
end do
end do
!
     ie = ikflavc(kc)

     IF (ie /= ie_ask) CYCLE 
     mass_factor = SQRT (xmass(1)/xmass(ikflavc(kc)))
!
     veff = v + vcorot -vpar + vm*alamc(kc)
     CALL Wrap_around_ghostcells (veff, isize, jsize, n_gc)
!
!    2. Differentiate Veff with respect to I and J:
!
     CALL Deriv_i (veff, isize, jsize, imin_j, dvefdi)
     CALL Deriv_j (veff, isize, jsize, imin_j, 1.0E+32, dvefdj)
     

     WHERE (ABS(dvefdj) > 1.0E+30)
        dvefdi = 0.0
        dvefdj = 0.0
     END WHERE
!
!
!   3. New advection scheme.
     loc_Eta  = 0.0
     loc_didt = 0.0
     loc_djdt = 0.0
     loc_rate = 0.0
!
     icut = 0
     DO j = 1, jsize
        icut = MAX (icut, imin_j(j))
        DO i = 1, isize
           IF (i < imin_j(j)) CYCLE
           IF (eeta(i,j,kc) > 1.0) icut = MAX (icut, i)
        END DO
     END DO
     icut = icut + 5
!
     fac = 1.0E-3*signbe*bir*alpha*beta*dlam*dpsi*ri**2

     DO j = 1, jsize
        DO i = 1, isize-1

           IF (i <= imin_j(j)) CYCLE


           ! ss, 04/19/2005--comment this out
!$         IF (dvefdj(i-1,j) > 0.0 .AND. dvefdj(i,j) > 0.0) THEN
!$            loc_didt (i,j) = dvefdj (i-1,j) / fac
!$         ELSE IF (dvefdj(i-1,j) < 0.0 .AND. dvefdj(i,j) < 0.0) THEN
!$            loc_didt (i,j) = dvefdj (i,j) / fac
!$         ELSE
!$            loc_didt (i,j) = 0.5*(dvefdj(i-1,j)+dvefdj(i,j))/fac
!$         END IF
!$         IF (dvefdi(i,j-1) < 0.0 .AND. dvefdi(i,j) < 0.0) THEN
!$            loc_djdt (i,j) = - dvefdi (i,j-1) / fac
!$         ELSE IF (dvefdi(i,j-1) > 0.0 .AND. dvefdi(i,j) > 0.0) THEN
!$            loc_djdt (i,j) = - dvefdi (i,j) / fac
!$         ELSE
!$            loc_djdt (i,j) = - 0.5*(dvefdi(i,j-1)+dvefdi(i,j))/fac
!$         END IF
           ! ss, 04/19/2005--end of comment this out

           loc_didt (i,j) = dvefdj (i-1,j) / fac(i-1,j)
           loc_djdt (i,j) = - dvefdi (i,j-1) / fac(i-1,j)


           IF (i > icut) THEN   
              loc_didt(i,j) = 0.0
              loc_djdt(i,j) = 0.0
           END IF

!$         loc_rate (i,j) = Ratefn (fudgec(kc), alamc(kc), sini (i,j),&
!$                                  bir (i,j), vm (i,j), mass_factor)
!          Determine the rate of loss for a given species:
!
           r_dist = SQRT(xmin(i,j)**2+ymin(i,j)**2)

           IF (ie == 1) THEN
              IF (kc == 1) THEN
!                Cold Electrons (plasmasphere). No pitch-angle scattering,
!                but account for the refilling from the ionosphere 
!                rate is > 0, so don't use is as a loss term!!!
                 if (i >= imin_j(j)) then
                 loc_rate(i,j) = &
                    Plasmasphere_refill_rate (rmin(i,j), doy, sunspot_number,&
                                              eeta(i,j,kc), vm(i,j))
!                eeta (i,j,kc) = eeta(i,j,kc)+(loc_rate(i,j)/eta_scale)*dt
                 eeta (i,j,kc) = eeta(i,j,kc)+(loc_rate(i,j))*dt
                 end if
                 loc_rate(i,j) = 0.0
              ELSE
!                Plasmasheet electrons, strong pitch-angle scattering loss:
                 loc_rate(i,j) = Ratefn (fudgec(kc), alamc(kc), sini(i,j),&
                                           bir (i,j), vm(i,j), mass_factor)
              END IF
           ELSE 
!             Positive ions, compute charge-exchange rate if it is on:
              IF (L_dktime) THEN
              loc_rate(i,j) = Cexrat (ie, ABS(alamc(kc))*vm(i,j), &
                                           R_dist, &
                                           sunspot_number, dktime, &
                                           irdk,inrgdk,isodk,iondk)
              END IF
           END IF

        END DO
!       eeta (1:imin_j(j)-1,j,kc) = etac (kc)
        loc_didt(isize,j) = loc_didt(isize-1,j)
        loc_djdt(isize,j) = loc_djdt(isize-1,j)
     END DO
!
!Copy to local variables
     loc_Eta (1:isize, 1:jsize) = eeta (1:isize, 1:jsize, kc)
!   
! 
!Call clawpack
     xlower = 0.0
     xupper = isize
     ylower = 0.0
     yupper = jsize
     FirstTime=.true.
     CALL Claw2ez (FirstTime, T1,T2, xlower,xupper, ylower,yupper, &
                   CLAWiter, 2,isize,jsize, &
                   loc_Eta, loc_didt, loc_djdt, loc_rate)
     FirstTime=.false.
 
!!$ i = 20
!!$ j = 1
!!$ print*,'eeta before:',i,j,kc,eeta(i,j,kc)
!!$ print*,'kc=',kc,' i=',i,' eeta(i,:,kc)='
!!$ print*,eeta(i,:,kc)
!!$ print*,'kc=',kc,' j=',j,' eeta(:,j,kc)='
!!$ print*,eeta(:,j,kc)

!Copy out
     eeta (:,:,kc) = loc_eta
!    DO j = 1, jsize
!       eeta (1:imin_j(j)-1, j, kc) = etac (kc)
!    END DO
!

do i = 1, isize
do j = 1, jsize
  if (eeta(i,j,kc) < 0.) then
     eeta(i,j,kc) = 0.0_rprec
!     print*,'eeta error: ',i,imin_j(j),j,kc,eeta(i,j,kc)
!!$     print*,'kc=',kc,' i=',i,' eeta(i,:,kc)='
!!$     print*,eeta(i,:,kc)
!!$     print*,'kc=',kc,' j=',j,' eeta(:,j,kc)='
!!$     print*,eeta(:,j,kc)
!!$     call CON_STOP('abort')
!!$     eeta(i,j,kc) = 0.0
!     eeta(i,j,kc) = 0.5*(eeta(max(1,i-1),j,kc)+eeta(min(isize,i+1),j,kc))
  end if
end do
end do
     CALL Wrap_around_ghostcells (eeta( :,:,kc), isize, jsize, n_gc)
!
  END DO
!  
  RETURN
!
  CONTAINS
!
    FUNCTION Ratefn (fudgx, alamx, sinix, birx, vmx, xmfact)
    IMPLICIT NONE
    REAL, INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
    REAL             :: Ratefn
!
!   Function subprogram to compute precipitation rate
!   Last update:  04-04-88
!
    Ratefn = 0.0466*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
    Ratefn = xmfact * ratefn
    RETURN
    END FUNCTION Ratefn
! 
END SUBROUTINE Move_plasma_grid_NEW
!
!
      SUBROUTINE Rcm_plasma_bc (iflag, i_where)
      USE Rcm_variables, ONLY : n_gc, isize, jsize, kcsize, iesize, &
                                alamc, etac, ikflavc, eeta, xmass, &
                                vm, Re, pi, imin_j, density, temperature,&
                                densityHp, densityOp, temperatureHp, temperatureOp, &
                                kmin, kmax, &
                                iprec, rprec, &
                                x_h, x_o,DoMultiFluidGMCoupling
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: iflag, i_where
!
!------------------------------------------------------------------------
! Given density and temperature on the rcm grid, ALAM channels, and the 
! magnetic field on the rcm grid, this subroutine will compute ion and 
! electron ETAs on the first grid points inside the boundary, by 
! converting N and T into kappa-distributions.
! Before calling, you must have values of:
!  * Re, pi
!  * alamc, xmass, ikflavc
!  * vm, imin_j
!  * density, temperature
! After calling, you will have values of EETA on boundary points.
! Above the boundary, EETA is extended using boundary values along J-lines.
! If I_WHERE is 1, do conversion only along the boundary (boundary condition),
! If I_WHERE is 2, do conversion everywhere in the region (initial condition).
!
! At the moment (10/31/2005), the MHD moments passed into this subroutine are:
! density--mass density in amu (rho/m_p)/m3
! temperature--P/(rho/m_p) in [eV], so temperature is essentially 
!   energy density without 1.5 factor
!------------------------------------------------------------------------
!
!
      INTEGER (iprec) :: kappa_e, kappa_i
      REAL (rprec) :: N_e, N_i, Eavg_e, Eavg_i, E0_e, E0_i, &
              E_i, E_e, eta_bnd, press_pt_old, press_pt_new
      INTEGER (iprec) :: k, k_beg, k_end, k_local, j, i, i_start, i_stop, ie, kc
      REAL (rprec) :: delta_e (kcsize), energy (kcsize)

      REAL (rprec) :: n_species (isize,jsize,iesize), temp_species (isize,jsize,iesize)
!
      REAL (rprec), PARAMETER :: a_conv = 6.371E+6/1.0E-9/1.0E-6
      REAL (rprec) :: a_factor, b_factor, s1, s2, s3, denom, pressure_rcm, density_rcm, &
           densityHp_rcm,pressureHp_rcm,densityOp_rcm,pressureOp_rcm
      INTEGER (iprec) :: unit_debug=59
      LOGICAL :: Flag_found_kuse, Flag_correct_ok (isize,jsize,iesize)=.TRUE.
!

      IF (iflag /= 2 .AND. iflag /= 3) call CON_STOP('IFLAG IN PLASMA_BC')
      IF (i_where /= 1 .AND. i_where /= 2) call CON_STOP('I_WHERE IN PLASMA_BC')
      kappa_e   = 6
      kappa_i   = 6
!
!     
      ! Using "density" from MHD, set partial number densities for species
      ! and their "temperatures" (average energies):

      DO j = 1, jsize

         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop
            if(.not. DoMultiFluidGMCoupling)then
               n_species (i,j,2) = density(i,j) * &
                    xmass(2) * x_h / (xmass(2)*x_h+xmass(3)*x_o)
               n_species (i,j,3) = density(i,j) * &
                    xmass(2) * x_o / (xmass(2)*x_h+xmass(3)*x_o)
               n_species (i,j,1) = n_species(i,j,2) + n_species(i,j,3)
               
               temp_species (i,j,2) = temperature (i,j) * &
                    (xmass(2)*x_h+xmass(3)*x_o) / xmass(2) / &
                    (1.0+1.0/7.8)
               temp_species (i,j,3) = temp_species (i,j,2)
               temp_species (i,j,1) = temp_species (i,j,2) / 7.8
            else
               n_species(i,j,2) = densityHp(i,j)
               n_species(i,j,3) = densityOp(i,j)
               n_species (i,j,1) = n_species(i,j,2) + n_species(i,j,3)
               
               temp_species (i,j,2) = temperatureHp (i,j)
               temp_species (i,j,3) = temperatureOp (i,j)
               temp_species (i,j,1) = temp_species (i,j,2)/ 7.8
               
        end if

         END DO
      END DO


!
!     H+
      ie = 2; k_beg = kmin(ie); k_end = kmax (ie)
!
      DO j = 1, jsize

         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop
         energy = ABS(alamc(:)) * vm(i,j)
         CALL Get_delta_e (kcsize, k_beg, k_end, energy, delta_e)
         Eavg_i    = 1.5 * temp_species(i,j,ie) ! in eV
         E0_i      = Convert_Eavg_to_E0 (Eavg_i, kappa_i) ! in eV
         N_i       = n_species (i,j,ie)
         DO k = k_beg, k_end
            E_i       = energy (k)
            eeta (i,j,k) = Jflux_kappa (N_i, xmass(ie), E0_i, E_i, kappa_i) * &
                                       4.0 * pi * SQRT(xmass(ie) / 2.0 / &
                                       (E_i*1.6E-19)) * &
                                       (delta_e(k)*1.6E-19)    ! number density in SI units
!           Convert to the ETA in program units ([ETA]=1/Weber=1/([B][L][L])):
            eeta (i,j,k) = eeta (i,j,k) * (vm(i,j)**(-1.5) / 1.0E-9 * Re * 1.0E+3)
         END DO
         END DO
      END DO
      DO k = k_beg, k_end
         DO j = 1, jsize
            eeta (1:imin_j(j)-1,j,k) = eeta(imin_j(j),j,k)
         END DO
      END DO
!

!     O+
      ie = 3; k_beg = kmin(ie); k_end = kmax (ie)

      DO j = 1, jsize

         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop
         energy = ABS(alamc(:)) * vm(i,j)
         CALL Get_delta_e (kcsize, k_beg, k_end, energy, delta_e)
         DO k = k_beg, k_end
            Eavg_i    = 1.5 * temp_species(i,j,ie) ! in eV
            E0_i      = Convert_Eavg_to_E0 (Eavg_i, kappa_i) ! in eV
            N_i       = n_species (i,j,ie)
            E_i       = energy (k)
            eeta (i,j,k) = Jflux_kappa (N_i, xmass(ie), E0_i, E_i, kappa_i) * &
                                       4.0 * pi * SQRT(xmass(ie) / 2.0 / &
                                       (E_i*1.6E-19)) * &
                                       (delta_e(k)*1.6E-19)    ! number density in SI units
!           Convert to the ETA in program units ([ETA]=1/Weber=1/([B][L][L])):
            eeta (i,j,k) = eeta (i,j,k) * (vm(i,j)**(-1.5) / 1.0E-9 * Re * 1.0E+3)
            END DO
         END DO
      END DO
      DO k = k_beg, k_end
         DO j = 1, jsize
            eeta (1:imin_j(j)-1,j,k) = eeta(imin_j(j),j,k)
         END DO
      END DO

!
!     Electrons:
!
      ie = 1; k_beg = kmin(ie)+1; k_end = kmax (ie)
      ! SS: "+1" in the above line is to avoid counting k=1 (lowest) 
      !     electron energy channel, which is now reserved for the plasmasphere
      !     population.
!
!
      DO j = 1, jsize

         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop
         energy = ABS(alamc(:)) * vm(i,j)
         CALL Get_delta_e (kcsize, k_beg, k_end, energy, delta_e)
         DO k = k_beg, k_end
            Eavg_e    = 1.5 * temp_species(i,j,ie) ! in eV
            E0_e      = Convert_Eavg_to_E0 (Eavg_e, kappa_e) ! in eV
            N_e       = n_species (i,j,ie)
            E_e       = energy (k)
            eeta (i,j,k) = Jflux_kappa (N_e, xmass(ie), E0_e, E_e, kappa_e) * &
                                       4.0 * pi * SQRT(xmass(ie) / 2.0 / &
                                       (E_e*1.6E-19)) * &
                                       (delta_e(k)*1.6E-19)    ! number density in SI units
!           Convert to the ETA in program units ([ETA]=1/Weber=1/([B][L][L])):
            eeta (i,j,k) = eeta (i,j,k) * (vm(i,j)**(-1.5) / 1.0E-9 * Re * 1.0E+3)
         END DO
      END DO
      END DO
      DO k = k_beg, k_end
      DO j = 1, jsize
         eeta (1:imin_j(j)-1,j,k) = eeta(imin_j(j),j,k)
      END DO
      END DO

      eeta (:,:,1) = 1.E+10  ! hack to keep plasmasphere in the code--will need to change later

      ! Check for negative eta's:
      DO j=1,jsize
      DO i=1,isize
      DO k=1,kcsize
         IF (eeta(i,j,k) < 0.0) THEN
             OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
             WRITE (UNIT_DEBUG,'(/////)')
             WRITE (UNIT_DEBUG,*) 'RCM, RCM_PLASMA_BC, ERROR:'
             WRITE (UNIT_DEBUG,*) 'NEGATIVE ETA VALUE AFTER FIRST PASS TO ASSIGN ETAS'
             WRITE (UNIT_DEBUG,*) 'I=',i,' J=',j,' K=',k,' IKFLAVC=',IKFLAVC(k),&
                                  ' ETA=', eeta(i,j,k)
             WRITE (UNIT_DEBUG,*) imin_j
             WRITE (UNIT_DEBUG,*) imin_j
             WRITE (UNIT_DEBUG,*) imin_j
             CLOSE (UNIT_DEBUG)
             CALL CON_STOP ('RCM_ERROR, SEE FILE RCM_DEBUG')
         END IF
      END DO
      END DO
      END DO
!
! Debugging code, not needed for regular runs:
   !  DO j = 1, jsize
   !  CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)
   !  DO i = i_start, i_stop
   !  density_rcm = 0.0
   !  pressure_rcm=0.0
   !  DO k = 2, kcsize
   !      density_rcm = density_rcm + (xmass(ikflavc(k))/xmass(2)) * (eeta(i,j,k)/6.37E+21) * vm(i,j)**1.5
   !      pressure_rcm= pressure_rcm+ ((ABS(alamc(k))*eeta(i,j,k))*1.67E-20) * ((vm(i,j)**2.5)*1.0E-6) ! nPa
   !  END DO
   !        IF (temperature(i,j) /= 0.0) THEN
   !            IF (ABS((temperature(i,j)-pressure_rcm/density(i,j)/1.6E-4)/temperature(i,j)) > 0.01) THEN
   !             WRITE (UNIT_DEBUG+1,'(/////)')
   !             WRITE (UNIT_DEBUG+1,*) 'RCM_D: IN RCM_PLASMA_BC THERE IS A MISMATCH IN ENERGY DENSITY MOMENTS:'
   !             WRITE (UNIT_DEBUG+1,*) 'MHD SENT: ',' N=',density(i,j),' T=',temperature(i,j)
   !             WRITE (UNIT_DEBUG+1,*) 'RCM  HAS: ',' N=',density_rcm ,' T=',pressure_rcm/density(i,j)/1.6E-4
   !             WRITE (UNIT_DEBUG+1,*) 'P RCM vs MHD: ',  pressure_rcm, temperature(i,j)*1.6E-4*density(i,j)
   !             WRITE (UNIT_DEBUG+1,*) 'LOCATION: ','I=',i,' J=',j
   !             DO k = 1, kcsize
   !                WRITE (UNIT_DEBUG+1,*) k,eeta(i,j,k) 
   !             END DO
   !            END IF
   !        END IF
   !  END DO
   !  END DO
!
!---------------------------------------------------------------------------
!   Cut off high-energy particle pressure (run this and the correction below
!   for any value of I_WHERE):
!
   !  DO j = 1, jsize
   !     DO i = imin_j(j), isize
   !        DO k = 2, kcsize
   !           IF (ABS(alamc(k))*vm(i,j) > 900000.0) eeta (i,j,k) = 0.0
   !        END DO
   !     END DO
   !  END DO
      !
      !
      ! Now adjust particle content so that the final distribution of EETA
      ! at each grid point yields the original moments (p and T) that came
      ! out of MHD. This will also correct for possible errors in the nume-
      ! rical procedure that splits the distribution function into discrete
      ! energy channels. For the second reason, this procedure should be run
      ! over all grid points, even if this routine is used for boundary 
      ! conditions only.

      Flag_correct_ok = .TRUE.

      DO j = 1, jsize

         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop

         Loop_check_species: DO ie = 1, iesize
           IF (ie == 1) THEN
              k_beg = kmin(ie)+1
           ELSE
              k_beg = kmin(ie)
           END IF

           ! Before we start with this, a quick check if the original MHD
           ! density that RCM got is zero. If this is the case, temperature
           ! is irrelevant and there should be nothing to correct, so skip 
           ! the rest of the loop:
           if(.not. DoMultiFluidGMCoupling)then
              IF (density(i,j) == 0.) THEN
                 CYCLE Loop_check_species
              END IF
           else
              if(densityOp(i,j) == 0. .or. densityHp(i,j) == 0.)then
                 CYCLE Loop_check_species
              end if
           end if
           
           ! Here start iterative process. We would like to apply the correction
           ! procedure to all energy channels. However, because the correction 
           ! factor is linear in energy, it is possible for it to become negative at high energies.
           ! It seems to happen for channels with flux-tube content of negligible (1.0E-9) fraction
           ! of the peak value. It is safe to zero out those energy channels, but to do so and 
           ! preserve the algorithm, we need to iterate. We start with the whole distribution function
           ! (k_end = kmax) and if that includes problematic channels, we descend in energy while
           ! repeating procedure, until it works.

           Loop_find_kuse: DO k_end = kmax(ie), k_beg+1, -1

              Flag_found_kuse = .TRUE.

              s1 = SUM(eeta(i,j,k_beg:k_end))/ a_conv
              s2 = SUM(eeta(i,j,k_beg:k_end)*ABS(alamc(k_beg:k_end))) / a_conv
              s3 = SUM(eeta(i,j,k_beg:k_end)*ABS(alamc(k_beg:k_end))*ABS(alamc(k_beg:k_end)))/a_conv
              denom = s1*s3 - s2*s2

              a_factor = ((n_species(i,j,ie)/vm(i,j)**1.5)*s3 - &
                (1.5*n_species(i,j,ie)*temp_species(i,j,ie) /vm(i,j)**2.5)*s2) / denom
              b_factor = (-(n_species(i,j,ie)/vm(i,j)**1.5)*s2 + &
                (1.5*n_species(i,j,ie)*temp_species(i,j,ie) /vm(i,j)**2.5)*s1) / denom

              DO k_local = k_beg, k_end
                 IF (a_factor + b_factor*ABS(alamc(k_local)) < 0.0) Flag_found_kuse = .FALSE.
              END DO

              IF (Flag_found_kuse) EXIT Loop_find_kuse

           END DO Loop_find_kuse

           ! It is possible that even for k_end=k_beg+1, the correction factor
           ! is still negative. Thus, add a check for it:
           IF (.NOT.Flag_found_kuse) THEN
               OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
               WRITE (UNIT_DEBUG,'(/////)')
               WRITE (UNIT_DEBUG,*) 'RCM, RCM_PLASMA_BC, ERROR OCCURED:'
               WRITE (UNIT_DEBUG,*) 'CORRECTION FACTOR IS NEGATIVE AT SOME ENERGIES'
               WRITE (UNIT_DEBUG,*) 'CORRECTION WILL NOT BE MADE FOR THIS SPECIES' 
               WRITE (UNIT_DEBUG,*) 'THIS HAPPENED AT I=',i,' J=',j,'IE=',ie
               WRITE (UNIT_DEBUG,*) 'S1*S3-S2**2=',s1*s3-s2*s2
               WRITE (UNIT_DEBUG,*) 'a_factor=', a_factor
               WRITE (UNIT_DEBUG,*) 'B_FACTOR=', b_factor
               WRITE(UNIT_DEBUG,*)' '
               if(.not. DoMultiFluidGMCoupling)then
                  WRITE(UNIT_DEBUG,*)'MHD PASSED TO RCM: ','N=',density(i,j), ' T=',temperature(i,j)
               else
                  WRITE(UNIT_DEBUG,*)&
                       'MHD PASSED TO RCM: ',&
                       'N_Hp=',densityHp(i,j), ' T_Hp=',temperatureHp(i,j),&
                       'N_Op=',densityOp(i,j), ' T_Op=',temperatureOp(i,j)
               end if
               WRITE(UNIT_DEBUG,*)'MHD ADAPTED FOR THIS SPECIES: ', 'N=',n_species(i,j,ie), ' T=',temp_species(i,j,ie)
               write(unit_debug,*)'all species: n=',n_species(i,j,:)
               WRITE(UNIT_DEBUG,*)''
               WRITE (UNIT_DEBUG,*) imin_j
               WRITE (UNIT_DEBUG,*) imin_j
               CLOSE (UNIT_DEBUG)
!              CALL CON_STOP ('RCM: ERROR IN RCM_PLASMA_BC, SEE FILE RCM_DEBUG')
               ! Force no correction (what else?)

               a_factor=1.0
               b_factor=0.0
               k_end = kmax(ie)
               Flag_correct_ok (i,j,ie) = .FALSE.
           END IF


           IF (s1*s3-s2*s2 < 0.01*s1*s3 .OR. &
               ABS(a_factor-1.0) > 0.2  .OR. &
               (s2 > 0.0 .AND. ABS(b_factor) > 0.2*s1/s2)) THEN

               ! Write out error message:
               OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
               WRITE (UNIT_DEBUG,'(/////)')
               WRITE (UNIT_DEBUG,*) 'RCM, RCM_PLASMA_BC, POSSIBLE PROBLEM OCCURED:'
               WRITE (UNIT_DEBUG,*) 'S1*S3-S2**2 IS CLOSE TO BEING NEGATIVE?'
               WRITE (UNIT_DEBUG,*) 'S1*S3-S2**2=',s1*s3-s2*s2
               WRITE (UNIT_DEBUG,*) 'A_FACTOR IS TOO BIG?'
               WRITE (UNIT_DEBUG,*) 'a_factor=', a_factor
               WRITE (UNIT_DEBUG,*) 'b_factor is too big?'
               WRITE (UNIT_DEBUG,*) 'B_FACTOR=', b_factor
               WRITE (UNIT_DEBUG,*) 'CORRECTION WILL NOT BE MADE FOR THIS SPECIES' 
               WRITE (UNIT_DEBUG,*) 'THIS HAPPENED AT I=',i,' J=',j,'IE=',ie
               WRITE(UNIT_DEBUG,*)'FOR THIS SPECIES, ENERGY RANGES ARE:',alamc(k_beg)*vm(i,j), alamc(kmax(ie))*vm(i,j),' eV'
               WRITE(UNIT_DEBUG,*)' '
               if(.not. DoMultiFluidGMCoupling)then
                  WRITE(UNIT_DEBUG,*)'MHD PASSED TO RCM: ','N=',density(i,j), ' T=',temperature(i,j)
               else
                  WRITE(UNIT_DEBUG,*)&
                       'MHD PASSED TO RCM: ',&
                       'N_Hp=',densityHp(i,j), ' T_Hp=',temperatureHp(i,j),&
                      'N_Op=',densityOp(i,j), ' T_Op=',temperatureOp(i,j)
               end if
               WRITE(UNIT_DEBUG,*)'MHD ADAPTED FOR THIS SPECIES: ', 'N=',n_species(i,j,ie), ' T=',temp_species(i,j,ie)
               write(unit_debug,*)'all species: n=',n_species(i,j,:)
               WRITE(UNIT_DEBUG,*)''
               WRITE(UNIT_DEBUG,*)'RESULTS OF CORRECTION PROCEDURE:'
               WRITE(UNIT_DEBUG,*)'     S1=',s1
               WRITE(UNIT_DEBUG,*)'     S2=',s2
               WRITE(UNIT_DEBUG,*)'     S3=',s3
               WRITE(UNIT_DEBUG,*)'     A_FACTOR=',a_factor
               WRITE(UNIT_DEBUG,*)'     B_FACTOR=',b_factor
               WRITE (UNIT_DEBUG,'(/////)')
               WRITE (UNIT_DEBUG,*) 'HERE IS INFORMATION FOR ALL CHANNELS OF THIS SPECIES:'
                   DO kc=k_beg,kmax(ie)
                      WRITE(UNIT_DEBUG,*) &
                           'K=',kc,' ALAM(K)=',alamc(kc), &
                           ' EETA=',eeta(i,j,kc)/(a_factor+b_factor*ABS(ALAMC(kc))),&
                           ' EETA_NEW=',eeta(i,j,kc)
                   END DO
               WRITE (UNIT_DEBUG,*) imin_j
               WRITE (UNIT_DEBUG,*) imin_j
               CLOSE (UNIT_DEBUG)

               IF (s1*s3-s2*s2 < 0.01*s1*s3) THEN 
                  ! Force no correction:
                  a_factor = 1.0
                  b_factor = 0.0
                  Flag_correct_ok (i,j,ie) = .FALSE.
               END IF

!              CALL CON_STOP ('RCM: ERROR IN RCM_PLASMA_BC, SEE FILE RCM_DEBUG')
           END IF

           DO k = k_beg, k_end
              eeta (i,j,k) = eeta (i,j,k)*(a_factor+b_factor*ABS(alamc(k)))
           END DO
           DO k = k_end+1, kmax(ie)
              eeta (i,j,k) = 0.0
           END DO

           ! another check for negative etas:
           DO k = k_beg, k_end
              if (eeta(i,j,k) < 0.0) then
                 OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
                 WRITE (UNIT_DEBUG,'(/////)')
                 WRITE(UNIT_DEBUG,*)'NEGATIVE EETA VALUE FOUND IN RCM_PLASMA_BC'
                 WRITE(UNIT_DEBUG,*)'AFTER THE CORRECTION PROCEDURE'
                 WRITE(UNIT_DEBUG,*)'THE FOLLOWING DETAILS SHOULD HELP:'
                 WRITE(UNIT_DEBUG,*)'SPATIAL GRID LOCATION: ', 'I=',I, ' J=',J, ' VM=',vm(i,j),' BND_I=',imin_j(j)
                 WRITE(UNIT_DEBUG,*)' '
                 WRITE(UNIT_DEBUG,*)'PARTICLE DETAILS: ','K=',K, ' IKFLAV=', Ie,' ALAM=',alamc(k)
                 WRITE(UNIT_DEBUG,*)'FOR THIS SPECIES, ENERGY RANGES ARE:',alamc(k_beg)*vm(i,j), alamc(k_end)*vm(i,j),' eV'
                 WRITE(UNIT_DEBUG,*)' '
                 if(.not. DoMultiFluidGMCoupling)then
                    WRITE(UNIT_DEBUG,*)'MHD PASSED TO RCM: ','N=',density(i,j), ' T=',temperature(i,j)
                 else
                    WRITE(UNIT_DEBUG,*)'MHD PASSED TO RCM: ','N_Hp=',densityHp(i,j), ' T_Hp=',temperatureHp(i,j),&
                         'N_Op=',densityOp(i,j), ' T_Op=',temperatureOp(i,j)
                 end if
                 WRITE(UNIT_DEBUG,*)'MHD ADAPTED FOR THIS SPECIES: ', 'N=',n_species(i,j,ie), ' T=',temp_species(i,j,ie)
                 write(unit_debug,*)'all species: n=',n_species(i,j,:)
                 WRITE(UNIT_DEBUG,*)'EETA VALUE DERIVED INITIALLY:', eeta(i,j,k)/(a_factor+b_factor*ABS(alamc(k)))
                 WRITE(UNIT_DEBUG,*)''
                 WRITE(UNIT_DEBUG,*)'RESULTS OF CORRECTION PROCEDURE:'
                 WRITE(UNIT_DEBUG,*)'     S1=',s1
                 WRITE(UNIT_DEBUG,*)'     S2=',s2
                 WRITE(UNIT_DEBUG,*)'     S3=',s3
                 WRITE(UNIT_DEBUG,*)'     A_FACTOR=',a_factor
                 WRITE(UNIT_DEBUG,*)'     B_FACTOR=',b_factor
                 WRITE(UNIT_DEBUG,*)'     NEW EETA=',eeta(i,j,k)
                 WRITE (UNIT_DEBUG,'(/////)')
                 WRITE (UNIT_DEBUG,*) 'HERE IS INFORMATION FOR ALL CHANNELS OF THIS SPECIES:'
                     DO kc=k_beg,k_end
                        WRITE(UNIT_DEBUG,*) &
                             'K=',kc,' ALAM(K)=',alamc(kc), &
                             ' EETA=',eeta(i,j,kc)/(a_factor+b_factor*ABS(ALAMC(kc))),&
                             ' EETA_NEW=',eeta(i,j,kc)
                     END DO
                 CLOSE (UNIT_DEBUG)
                 CALL CON_STOP ('RCM DEBUG: STOPPING')
              END IF
           END DO

         END DO Loop_check_species

      END DO
      END DO
!
!---------------------------------------------------------------------------
!

      ! Final check: the RCM moments must match the MHD-supplied moments
      ! exactly within given accuracy since we distorted the shape of the
      ! initial distribution function. Verify this:

      DO j = 1, jsize
      
         CALL Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)

         DO i = i_start, i_stop

            IF (ANY (.NOT.Flag_correct_ok (i,j,:))) CYCLE

            if(.not. DoMultiFluidGMCoupling)then
            density_rcm = 0.0
            pressure_rcm = 0.0

            DO k = 2, kcsize
            density_rcm = density_rcm + (xmass(ikflavc(k))/xmass(2)) * (eeta(i,j,k)/6.37E+21) * vm(i,j)**1.5
            pressure_rcm= pressure_rcm+ ((ABS(alamc(k))*eeta(i,j,k))*1.67E-20) * ((vm(i,j)**2.5)*1.0E-6) ! nPa
            END DO

            IF (density(i,j) /= 0.0 ) THEN
                IF (ABS(density_rcm-density(i,j))/density(i,j) > 0.01) THEN
                 OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
                 WRITE (UNIT_DEBUG,'(/////)')
                 WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN MASS DENSITY MOMENTS:'
                 WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N=',density(i,j),' T=', temperature(i,j)
                 WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N=',density_rcm, ' T=', pressure_rcm/density(i,j)/1.6E-4
                 WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j, ' IMIN_J(J)=',imin_j(J)
                 WRITE (UNIT_DEBUG,*) ' '
                 DO k = 1, kcsize
                    WRITE (UNIT_DEBUG,*) k,eeta(i,j,k) 
                 END DO
                 CLOSE (UNIT_DEBUG)
                 CALL CON_STOP ('RCM DEBUG: STOPPING')
                END IF
            END IF

            IF (temperature(i,j) /= 0.0 .and. density(i,j) /= 0.0) THEN
                IF (ABS((temperature(i,j)-pressure_rcm/density(i,j)/1.6E-4)/temperature(i,j)) > 0.01) THEN
                 OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG')
                 WRITE (UNIT_DEBUG,'(/////)')
                 WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN ENERGY DENSITY MOMENTS:'
                 WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N=',density(i,j),' T=',temperature(i,j)
                 WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N=',density_rcm ,' T=',pressure_rcm/density(i,j)/1.6E-4
                 WRITE (UNIT_DEBUG,*) 'P RCM vs MHD: ',  pressure_rcm, temperature(i,j)*1.6E-4*density(i,j)
                 WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j
                 DO k = 1, kcsize
                       WRITE (UNIT_DEBUG,*) k, alamc(k), ikflavc(k), eeta(i,j,k)
                 END DO
                 CLOSE (UNIT_DEBUG)
                 CALL CON_STOP ('RCM DEBUG: STOPPING')
                END IF
            END IF
         else
            !MultiFluid                                                                                                                         
            densityHp_rcm = 0.0
            pressureHp_rcm = 0.
            densityOp_rcm = 0.0
            pressureOp_rcm = 0.0
            
            do k=kmin(2),kmax(2)
               densityHp_rcm = densityHp_rcm + (xmass(ikflavc(k))/xmass(2)) * (eeta(i,j,k)/6.37E+21) * vm(i,j)**1.5
               pressureHp_rcm= pressureHp_rcm+ ((ABS(alamc(k))*eeta(i,j,k))*1.67E-20) * ((vm(i,j)**2.5)*1.0E-6) ! nPa                           
            end do
            do k=kmin(3),kmax(3)
               densityOp_rcm = densityOp_rcm + (xmass(ikflavc(k))/xmass(3)) * (eeta(i,j,k)/6.37E+21) * vm(i,j)**1.5
               pressureOp_rcm= pressureOp_rcm+ ((ABS(alamc(k))*eeta(i,j,k))*1.67E-20) * ((vm(i,j)**2.5)*1.0E-6) ! nPa                           
            end do
            IF (densityHp(i,j) /= 0.0 ) THEN
               IF (ABS(densityHp_rcm-densityHp(i,j))/densityHp(i,j) > 0.01) THEN
                  OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG1')
                  WRITE (UNIT_DEBUG,'(/////)')
                  WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN MASS DENSITY MOMENTS:'
                  WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N_Hp=',densityHp(i,j),' T_Hp=', temperatureHp(i,j)
                  WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N_Hp=',densityHp_rcm, ' T_Hp=', pressureHp_rcm/densityHp_rcm/1.6E-4
                  WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j, ' IMIN_J(J)=',imin_j(J)
                  WRITE (UNIT_DEBUG,*) ' '
                  DO k = kmin(2),kmax(2)
                     WRITE (UNIT_DEBUG,*) k,eeta(i,j,k)
                  END DO
                  CLOSE (UNIT_DEBUG)
                  CALL CON_STOP ('RCM DEBUG densityHp: STOPPING')

               END IF
            END IF

            IF (densityOp(i,j) /= 0.0 ) THEN
               IF (ABS(densityOp_rcm-densityOp(i,j))/densityOp(i,j) > 0.01) THEN
                  OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG2')
                  WRITE (UNIT_DEBUG,'(/////)')
                  WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN MASS DENSITY MOMENTS:'
                  WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N_Op=',densityOp(i,j),' T_Op=', temperatureOp(i,j)
                  WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N_Op=',densityOp_rcm, ' T_Op=', pressureHp_rcm/densityOp_rcm/1.6E-4
                  WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j, ' IMIN_J(J)=',imin_j(J)
                  WRITE (UNIT_DEBUG,*) ' '
                  DO k = kmin(3),kmax(3)
                     WRITE (UNIT_DEBUG,*) k,eeta(i,j,k)
                  END DO
                  CLOSE (UNIT_DEBUG)
                  CALL CON_STOP ('RCM DEBUG densityOp: STOPPING')
               END IF
            END IF

            IF (temperatureHp(i,j) /= 0.0 .and. densityHp(i,j) /= 0.0) THEN
               IF (ABS((temperatureHp(i,j)-pressureHp_rcm/densityHp_rcm/1.6E-4)/temperatureHp(i,j)) > 0.01) THEN
                  OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG3')
                  WRITE (UNIT_DEBUG,'(/////)')
                  WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN ENERGY DENSITY MOMENTS:'
                  WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N_Hp=',densityHp(i,j),' T_Hp=',temperatureHp(i,j)
                  WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N_Hp=',densityHp_rcm ,' T_Hp=',pressureHp_rcm/densityHp_rcm/1.6E-4
                  WRITE (UNIT_DEBUG,*) 'P_Hp RCM vs MHD: ',  pressureHp_rcm, temperatureHp(i,j)*1.6E-4*densityHp(i,j)
                  WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j
                  DO k = kmin(2), kmax(2)
                     WRITE (UNIT_DEBUG,*) k, alamc(k), ikflavc(k), eeta(i,j,k)
                  END DO
                  CLOSE (UNIT_DEBUG)
                  CALL CON_STOP ('RCM DEBUG THp : STOPPING')
               END IF
            END IF
            
            IF (temperatureOp(i,j) /= 0.0 .and. densityOp(i,j) /= 0.0) THEN
               IF (ABS((temperatureOp(i,j)-pressureOp_rcm/densityOp_rcm/1.6E-4)/temperatureOp(i,j)) > 0.01) THEN
                  OPEN (unit=UNIT_DEBUG,FILE='RCM_DEBUG4')
                  WRITE (UNIT_DEBUG,'(/////)')
                  WRITE (UNIT_DEBUG,*) 'RCM: IN RCM_PLASMA_BC THERE IS A MISMATCH IN ENERGY DENSITY MOMENTS:'
                  WRITE (UNIT_DEBUG,*) 'MHD SENT: ',' N_Op=',densityOp(i,j),' T_Op=',temperatureOp(i,j)
                  WRITE (UNIT_DEBUG,*) 'RCM  HAS: ',' N_Op=',densityOp_rcm ,' T_Op=',pressureOp_rcm/densityOp_rcm/1.6E-4
                  WRITE (UNIT_DEBUG,*) 'P_Op RCM vs MHD: ',  pressureOp_rcm, temperatureOp(i,j)*1.6E-4*densityOp(i,j)
                  WRITE (UNIT_DEBUG,*) 'LOCATION: ','I=',i,' J=',j
                  DO k = kmin(2), kmax(2)
                     WRITE (UNIT_DEBUG,*) k, alamc(k), ikflavc(k), eeta(i,j,k)
                  END DO
                  CLOSE (UNIT_DEBUG)
                  CALL CON_STOP ('RCM DEBUG TOp: STOPPING')
               END IF
            END IF

         end if

         END DO

      END DO


      DO k = 1, kcsize
         CALL Wrap_around_ghostcells (eeta(:,:,k), isize, jsize, n_gc)
      END DO
!
      IF (iflag == 2) THEN
         RETURN
      ELSE IF (iflag == 3) THEN
         DO k = 1, kcsize
            eta_bnd = eeta(imin_j(jsize/2),jsize/2,k)
            DO j = 1-n_gc, jsize+n_gc
               eeta(imin_j(j),j,k) = eta_bnd
               eeta (1:imin_j(j)-1,j,k) = eta_bnd
            END DO
         END DO
      DO k = 1, kcsize
         CALL Wrap_around_ghostcells (eeta(:,:,k), isize, jsize, n_gc)
      END DO
!
      ELSE
         call CON_STOP('ILLEGAL CONTROL FLAG IN BOUNDARY PLASMA ROUTINE')
      END IF
!
      RETURN
      CONTAINS


      SUBROUTINE Rcm_plasma_bc_get_i_range (i_where, j, i_start, i_stop)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: i_where, j
      INTEGER(iprec), INTENT (OUT) :: i_start, i_stop

      i_start = imin_j(j)
      IF (i_where == 1) THEN
         i_stop = i_start
      ELSE
         i_stop = isize
      END IF
      RETURN
      END SUBROUTINE Rcm_plasma_bc_get_i_range



      FUNCTION Fdist_kappa_dflux (N, mass, E0, E, kappa)
!
!     Generic kappa distribution function that returns differential
!     (directional) particle flux (particles/cm2/sec/srad) in SI units
!     for particles of given MASS, density N, characteristic energy E0,
!     spectral index KAPPA, and at the energy E.
!
!     Units: N      in m-3
!            mass   in kg
!            E0     in Joules
!            E      in Joules
!            OUTPUT in particles/m2/sec/srad
!
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: kappa
      REAL(rprec), INTENT (IN) :: N, mass, E0, E
      REAL(rprec)              :: Fdist_kappa_dflux
      DOUBLE PRECISION ::  DGamma
      EXTERNAL             DGamma
!
      DOUBLE PRECISION :: c1, c2, c3, c4
      
      IF (N == 0.0) THEN
         Fdist_kappa_dflux = 0.0
         RETURN
      ELSE
         c1 = N / SQRT (2.0*mass) / SQRT (E0)    ! to avoid under/overflows
         c2 = DGamma (DBLE(kappa+1)) / DGamma (DBLE(kappa-0.5)) / (pi*kappa)**1.5
         c3 = E/E0
         c4 = (1.0 + E / E0 / kappa)**(-kappa-1)
!        c4 = EXP(LOG(1.0 + E / E0 / kappa)*(-kappa-1))
         Fdist_kappa_dflux = c1 * c2 * c3 * c4
         RETURN
      END IF
      END FUNCTION Fdist_kappa_dflux
!
!
!
!
!
      FUNCTION Jflux_kappa (N, mass, E0, E, kappa)
      IMPLICIT NONE
      REAL(rprec), INTENT (IN) :: N, mass, E0, E
      INTEGER (iprec), INTENT(IN) ::  kappa
      REAL(rprec) :: Jflux_kappa
!-------------------------------------------------------------------------------------
!     Kappa energy distribution function. Energy is in eV, mass in kg, density in cm-3
!     Returns differential (directional) flux is in SI units.
!-------------------------------------------------------------------------------------

      IF (E0 == 0.0) CALL CON_STOP ('E0 IS ZERO IN RCM[JFLUX_KAPPA]')

      Jflux_kappa = Fdist_kappa_dflux (N*1.0E+6, mass, E0*1.6E-19, E*1.6E-19, kappa)
      RETURN
      END FUNCTION Jflux_kappa
!
!
!
!
!
!     FUNCTION Convert_E0_to_Eavg (E0, kappa)
!     IMPLICIT NONE
!     REAL, INTENT (IN) :: E0, kappa
!     REAL              :: Convert_E0_to_Eavg
!
!     Compute Eavg (average energy) from E0 for a kappa distribution
!
!     Convert_E0_to_Eavg = 1.5*E0*kappa/(kappa-1.5)
!     RETURN
!     END FUNCTION Convert_E0_to_Eavg
!
!
!
!
!
      FUNCTION Convert_Eavg_to_E0 (Eavg, kappa)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: kappa
      REAL(rprec), INTENT (IN) :: Eavg
      REAL(rprec)              :: Convert_Eavg_to_E0
!
!     Compute E0 from Eavg for a kappa distribution
!
      Convert_Eavg_to_E0 = Eavg/1.5 *(kappa-1.5)/kappa
      RETURN
      END FUNCTION Convert_Eavg_to_E0
!
      END SUBROUTINE Rcm_plasma_bc
!
!
!
      SUBROUTINE Get_delta_E (kcsize, k_beg, k_end, energy, delta_E)
      USE Rcm_variables, ONLY : rprec, iprec
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: kcsize, k_beg, k_end
      REAL(rprec), INTENT (IN)    :: energy (kcsize)
      REAL(rprec), INTENT (IN OUT)   :: delta_e (kcsize)
!
!  Compute widths of energy channels numbered by index K running
!  from KBEG to KEND (KEND > KBEG).
!  We assume that ENERGY(kbeg) is the lowest-energy channel.
!  Then the width of the k-th channel is the distance
!  between centers of intervals [energy(k), energy(k+1)] and
!  [energy(k-1),energy(k)], with the exception of the boundary cases, namely:
!  if k=kbeg, then the width is from energy=0 to the midpoint of
!     [energy(kbeg), energy(kbeg+1)],
!  if k=kend, then width of the last channel is extended to the right of
!     energy(kend) by (energy(kend)-energy(kend-1))/2.
!
!
      INTEGER :: k
!
!
      IF (k_beg >= k_end) THEN
         call CON_STOP('k_end is not larger than k_beg, cannot proceed')
      END IF
!
      DO k = k_beg, k_end
         IF (k == k_beg) THEN
            delta_E(k) = (energy(k+1) + energy(k) ) / 2.0
         ELSE IF (k == k_end) THEN
            delta_E(k) = energy(k) - energy(k-1)
         ELSE
            delta_E(k) = ( energy(k+1) - energy(k-1) ) / 2.0
         END IF
      END DO
      RETURN
      END SUBROUTINE Get_delta_E




      FUNCTION Plasmasphere_refill_rate (r, doy, sunspot_number, eta_pt, vm_pt)
      USE Rcm_variables, ONLY : rprec, iprec, trf
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: r, doy, sunspot_number, eta_pt, vm_pt
      REAL (rprec) :: Plasmasphere_refill_rate
      REAL (rprec) :: Plasmasphere_den_CA92_sat
!
!     Compute Lambour et al refill rate for a given point at radial
!     distance R that has flux tube content ETA_PT. The rate is
!     d(eta)/dt = (ETA_sat - Eta_pt) / tau, where ETA_sat is the 
!     saturated plasmasphere density, Eta_pt is current density, and 
!     tau is the refilling time constant interpolated to given radial
!     distance for given conditions. Rate is positive, but in the rcm
!     advection scheme must be negative since this is treated as a loss term.
!
      INTEGER (iprec) :: irf
      REAL (rprec) :: taurf, eta_sat
!
      DO irf = 1, 18
         IF (r >= trf(irf,1) .AND. r <= trf(irf+1,1)) THEN
             taurf = trf(irf,5)+ ((trf(irf+1,5)-trf(irf,5))/&
                     (trf(irf+1,1)-trf(irf,1))*(r-trf(irf,1)))
         END IF
      END DO
      IF (r > trf(19,1)) THEN
          taurf = trf(19,5)+((trf(19,5)-trf(18,5))/(trf(19,1)-&
                  trf(18,1)))*(r-trf(19,1))
      END IF
      IF (r < trf(1,1)) taurf = trf (1,5)
      taurf = taurf * 86400.0 ! convert from days to seconds
      eta_sat = 6.38E+21*vm_pt**(-1.5)* &
                Plasmasphere_den_CA92_sat(r,doy,sunspot_number)
      Plasmasphere_refill_rate = (eta_sat - eta_pt)/taurf
      RETURN
      END FUNCTION Plasmasphere_refill_rate



!
      FUNCTION Cexrat (isp,enrg,rloc,ssn,dktime,irdk,inrgdk,isoldk, &
                       iondk)
      USE Rcm_variables, ONLY : rprec
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: isp, irdk, inrgdk, isoldk, iondk
      REAL (rprec), INTENT (IN) :: enrg, rloc, ssn, dktime (irdk,inrgdk,isoldk,iondk)
      REAL (rprec) :: Cexrat
!
!-------------------------------------------------------------------------
!  copyright rice university, 1993
!
!  version 1.00                                 05.09.90
!          2.00                                 02.04.90
!                                       msm delivery version
!          2.10                                 06.11.93
!               error output routed to unit 9
!
!  programmer: r. w. spiro
!
!  purpose:  function subprogram to return charge exchange loss rate
!          (sec**(-1)) for ions of species isp, energy enrg (ev) at
!          l=rloc (re) for sunspot number ssn.  this routine is based
!          on a table generated by james bishop of u. of michigan.
!
!  calling parameters
!        isp       species identifier
!                    isp=2 for h+ ions
!                    isp=3 for o+ ions
!        enrg      energy in ev
!        rloc      radial location (re)
!        ssn       sunspot number
!        dktime    table of ion decay times
!        irdk      radial dimension of dktime array
!        inrgdk    energy dimension of dktime array
!        isoldk    sunspot number dimension of dktime array
!        iondk     number of ion species in dktime array
!--------------------------------------------------------------------------------
!
      INTEGER, PARAMETER ::irsiz=18,inrgsz=13,isolsz=2,ionsiz=2
      REAL ::  elgvec(inrgsz), rvec(irsiz),ssnvec(2), &
               enrglg, br, bnrg, ssnuse, bssn, decayt, G3ntrp
      INTEGER :: ispndx, ir, inrg
!
      DATA elgvec /2.50,2.75,3.00,3.25,3.50,3.75,4.00, &
                   4.25,4.50,4.75,5.00,5.25,5.50/
!
      DATA rvec /1.50,2.00,2.50,3.00, &
                 3.50,4.00,4.50,5.00, &
                 5.50,6.00,6.50,7.00, &
                 7.50,8.00,8.50,9.00, &
                 9.50,10.00/
!
      DATA ssnvec /0.0,100./
!
!
      IF (irsiz /= irdk .OR. inrgsz /= inrgdk .OR. &
          ionsiz /= iondk .OR. isolsz /= isoldk) THEN
         write(*,*) 'dimension error in function cexrat'
         write(*,*) 'irdk,inrgdk,iondk,isoldk',irdk,inrgdk,iondk,isoldk
         write(*,*) 'irsiz,inrgsz,ionsiz,isolsz',irsiz,inrgsz,ionsiz,isolsz
         write(*,*) 'stopping program in cexrat'
         call CON_stop('ERROR in IM/RCM2/src/rcm_routines.f90')
      END IF
!
      enrglg = LOG10(enrg) !  work with log10 of particle energy
      ispndx=isp-1
!
      if_1: IF (rloc <= rvec(1)) THEN !  find br for interpolation
         br=1.0
      ELSE IF (rloc > rvec(irdk)) THEN
         br=irdk
      ELSE
         do_1: DO ir=1,irdk-1
            IF (rloc <= rvec(ir+1)) THEN
               br=ir+(rloc-rvec(ir))/(rvec(ir+1)-rvec(ir))
               EXIT do_1
            END IF
         END DO do_1
      END IF if_1
!
      if_2: IF (enrglg.le.elgvec(1)) THEN !  find bnrg for interpolation
         bnrg = 1.0
      ELSE IF (enrglg > elgvec(inrgdk)) THEN
         bnrg = inrgdk
      ELSE
         do_2: DO inrg=1,inrgdk-1
            IF (enrglg <= elgvec(inrg+1)) THEN
               bnrg=inrg+(enrglg-elgvec(inrg))/(elgvec(inrg+1)-elgvec(inrg))
               EXIT do_2
            END IF
         END DO do_2
      END IF if_2
!
!**********  change 9/30/91  *****************************************
!  if ssn.gt.ssnvec(2), then use ssnvec(2) for ssn
      ssnuse=ssn
      IF (ssnuse > ssnvec(2)) ssnuse=ssnvec(2)
!
!*********  end change  9/30/91  ************************************
!
!  find bssn for interpolation
      bssn=1.0+(ssnuse-ssnvec(1))/(ssnvec(2)-ssnvec(1))
!
!  decayt is decay time in seconds
      decayt = G3ntrp (dktime(1,1,1,ispndx),irdk,inrgdk,isoldk,br,bnrg,bssn)
!
      IF (ABS(decayt) < 1.0E-20) THEN
         write(*,*) 'decayt is less than 1.e-20 sec in cexrat'
         write(*,*) 'decayt=',decayt,'  br=',br,'  bnrg=',bnrg,'bssn=',bssn
         write(*,*) 'isp=',isp,'  enrg=',enrg,' rloc=',rloc,' ssn=',ssn
         write(*,*) 'ssnuse=',ssnuse
      END IF
!
!  to get charge exchange rate (sec**9-1)) cexrat, invert decayt
!
      cexrat=1.0/decayt
      RETURN
      END FUNCTION Cexrat
!
!
      FUNCTION G3ntrp (a,imax,jmax,kmax,bi,bj,bk)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: imax, jmax, kmax
      REAL, INTENT (IN) :: a(imax,jmax,kmax), bi, bj, bk
      REAL :: G3ntrp
!
!---------------------------------------------------------------------------
!  copyright Rice University, 1993
!
!  VERSION 1.00                            DATE: 01.11.88
!          1.01A                                 02.02.89
!          2.00  MSM DELIVERY VERSION            01.28.93
!
!  PURPOSE: FUNCTION SUBPROGRAM TO PERFORM A GENERAL 3-D LINEAR
!           INTERPOLATION OF ARRAY A(I,J,K) AT PT(BV(1),BV(2),BV(3))
!
!  INPUT:
!       A          3-D ARRAY TO BE INTERPOLATED
!       IMAX       I DIMENSION OF ARRAY A
!       JMAX       J DIMENSION OF ARRAY A
!       KMAX       K DIMENSION OF ARRAY A
!       BI         FLOATING POINT VALUE TO INTERPOLATE IN I DIMENSION
!       BJ         FLOATING POINT VALUE TO INTERPOLATE IN J DIMENSION
!       BK         FLOATING POINT VALUE TO INTERPOLATE IN K DIMENSION
!
!
!  OUTPUT:
!       G3NTRP     INTERPOLATED VALUES OF ARRAY A
!----------------------------------------------------------------------
!
!
      INTEGER ::  ndx(3),ndim(3), kstop, jstop, L, i, j, k
      REAL ::  BV(3),COEF(3,2), fndx
!
      NDIM(1)=IMAX
      NDIM(2)=JMAX
      NDIM(3)=KMAX
      BV(1)=BI
      BV(2)=BJ
      BV(3)=BK
      DO L=1,3
         NDX(L)=BV(L)
         IF(NDX(L).LT.1) NDX(L)=1
         IF(NDX(L).GT.NDIM(L)-1) NDX(L)=NDIM(L)-1
         IF(NDX(L).LE.0) NDX(L)=1
         FNDX=REAL(NDX(L))
         COEF(L,1)=1.-BV(L)+FNDX
         COEF(L,2)=BV(L)-FNDX
      END DO
!
      G3NTRP=0.
      kstop = MIN(KMAX,2)
      jstop = MIN(JMAX,2)
      DO I=1,2
      DO J=1,jstop
      DO K=1,kstop
         G3ntrp=G3ntrp+ &
          coef(1,i)*coef(2,j)*coef(3,k)*a(ndx(1)+i-1,ndx(2)+j-1,ndx(3)+k-1)
      END DO
      END DO
      END DO
!
      RETURN
      END FUNCTION G3ntrp
!
!
!
      FUNCTION Plasmasphere_den_CA92_sat (r, doy, sunspot_number)
      USE Rcm_variables, ONLY : rprec, pi
!
!     Carpenter_and_Anderson_JGR_1992 saturated plasmaspheric density in cm-3
!     as a function of radial distance, day of year, and mean monthly sunspot
!     number. Reduced by 5% following Lambour (accounting for ISEE-1 off-equat
!     location?) 
!     Valid earthward of the plasmaspheric trough.
!
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: r, doy, sunspot_number
      REAL (rprec) :: Plasmasphere_den_CA92_sat
!
      Plasmasphere_den_CA92_sat = 10.0_rprec**( &
        (-0.3145*r+3.9043) + ( &
           0.15*COS(2*pi*(doy+9)/365.) - 0.075*COS(4*pi*(doy+9)/365.) + &
           0.00127*sunspot_number-0.0635) *EXP((2.-r)/1.5) &
                                               )
      Plasmasphere_den_CA92_sat = Plasmasphere_den_CA92_sat*0.95
      RETURN
      END FUNCTION Plasmasphere_den_CA92_sat
