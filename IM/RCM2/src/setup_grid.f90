    PROGRAM Rcm_setup_grid
    USE Rcm_variables, ONLY : iprec, rprec, pi, pi_two, DTR, besu, HTR,&
                             dlam, dpsi, Re, Ri, RTH, RTD, &
                             LUN, label, n_gc, &
                             colat, aloct, alpha, beta, vcorot, sini, bir
    USE Rcm_io
    IMPLICIT NONE
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   MSM grid info:
    INTEGER (iprec), PARAMETER :: isize_msm = 62, &   ! MSM grid params
                                  jsize_msm = 51, &
                                  jwrap_msm = 3,  &
                                  delta_i = 2,    &
                                  i_pc = 6,       &
                                  i_pp = 55
    REAL (rprec), PARAMETER    :: dtheta_di_max = 2.0_rprec,  &
                                  dtheta_di_min = 0.4_rprec
    REAL(rprec) ::  colat_msm (isize_msm, jsize_msm), &
                    aloct_msm (isize_msm, jsize_msm)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!   GRID #2: MSM + more points at lower latitudes down to 10 degrees:
    INTEGER (iprec), PARAMETER :: isize_2 = 78,  &
                                  jsize_2 = 48,  &
                                  jwrap_2 = 3
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
    INTEGER (iprec), PARAMETER :: isize_3 = 155, &    ! high-resolution grid
                                  jsize_3 = 99,  &
                                  jwrap_3 = 3
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
    INTEGER (iprec) :: isize, jsize, jwrap, ncol
!   REAL (rprec), ALLOCATABLE :: colat(:,:), aloct(:,:), alpha(:,:), &
!                                beta(:,:), vcorot(:,:), sini(:,:), &
!                                bir (:,:)
    INTEGER(iprec), ALLOCATABLE :: msm_i_index (:,:), &
                                   msm_j_index (:,:)
    INTEGER(iprec)       :: grid_type, i, j, istat, i_msm, j_msm
    REAL (rprec)         :: dthetadi
    CHARACTER (LEN=80)   :: filename, form_string
    CHARACTER (LEN=1)    :: answer_1char
    LOGICAL, ALLOCATABLE :: msm_mask (:,:)
    LOGICAL              :: logic_flag_1, logic_flag_2
!
!
     IF (.NOT.Check_logical_units () ) THEN
        STOP 'LOGICAL UNITS ARE NOT AVAILABLE, STOP'
     END IF
!
    INQUIRE (FILE = trim(NameRcmDir)//'input/rcmcrd11', EXIST = logic_flag_1)
    INQUIRE (FILE = trim(NameRcmDir)//'input/rcmcrd21', EXIST = logic_flag_2)
    IF (logic_flag_1 .OR. logic_flag_2) THEN
       WRITE (*,'(//,T2,A)') 'AT LEAST ONE OF THE 2 GRID FILES EXISTS'
       WRITE (*,'(T2,A)')    'IF YOU PROCEED, THESE FILES WILL BE OVERWRITTEN'
       WRITE (*,'(T2,A)',ADVANCE='NO') 'PROCEED WITH GRID SETUP? (Y/N)----->'
       READ (*,*) answer_1char
       IF (answer_1char /= 'Y' .AND. answer_1char /= 'y') STOP
    END IF
    DO
      WRITE (*,'(T2,/)')
      WRITE (*,'(T2,A)') '   OPTIONS TO CREATE IONOSPHERIC RCM GRID:'
      WRITE (*,'(T2,A)') '   ======================================='
      WRITE (*,'(T2,A)') ' 1. USE MSM GRID (62 X 51)'
      WRITE (*,'(T2,A)') ' 2. USE MSM+MORE POINTS IN LATITUDE (78 X 51)'
      WRITE (*,'(T2,A)') ' 3. USE OPTION 2 WITH TWICE RESOLUTION (155 X 99)'
      WRITE (*,'(T2,A)') ' 4. USE UNIFORM GRID 80 X 51  (DEBUG!)'
      WRITE (*,'(T2,A)') ' 5. USE UNIFORM GRID 160 X 99 (DEBUG!)'
      WRITE (*,'(T2,A)',ADVANCE='NO') 'ENTER OPTION (1, 2, 3, 4, OR 5)-->:'
      READ (*,*) grid_type
      IF (grid_type > 0 .AND. grid_type <= 5 ) EXIT
    END DO
!
!
    CALL Make_msm_grid ()
!
    SELECT CASE (grid_type)
!
    CASE (1) 
!
       stop 'disabled'
       isize = isize_msm
       jsize = jsize_msm
       jwrap = jwrap_msm

!      CALL Make_grid_arrays (isize,jsize)

       DO i = 1, isize
       DO j = 1, jsize
          colat(i,j) = colat_msm (i,j)
          aloct(i,j) = aloct_msm (i,j)
       END DO
       END DO
!
    CASE (2)
!
       isize = isize_2
       jsize = jsize_2
       jwrap = jwrap_2

!      CALL Make_grid_arrays (isize,jsize)

       DO i = 1, isize
       DO j = 1, jsize
          IF (i<= isize_msm) THEN
             colat (i,j) = colat_msm (i,j+2)
             aloct (i,j) = aloct_msm (i,j+2)
          ELSE
             colat (i,j) = colat_msm (1,j+2) + &
                   Theta_of_i (i, dtheta_di_max, dtheta_di_min, delta_i,&
                               i_pc, i_pp) * DTR
             aloct (i,j ) = aloct_msm (1,j+2)
          END IF
       END DO
       END DO
!
       CALL Wrap_around_ghostcells (colat, isize, jsize, n_gc)
       CALL Wrap_around_ghostcells (aloct, isize, jsize, n_gc)
!
!
    CASE (3)
!
       stop 'disabled'
       isize = isize_3
       jsize = jsize_3
       jwrap = jwrap_3

!      CALL Make_grid_arrays (isize,jsize)

       DO i = 1, isize
       DO j = jwrap, jsize-1
          colat (i,j) = colat_msm (1,1) + &
                Theta_of_i (i, 0.5*dtheta_di_max, 0.5*dtheta_di_min, &
                            2*delta_i, 2*i_pc-1, 2*i_pp-1) * DTR
          aloct(i,j) = 2*pi*REAL(j-jwrap)/REAL(jsize-jwrap)
          IF (aloct(i,j) < 0.) aloct(i,j) = aloct(i,j) + 2*pi
       END DO
       END DO
!
       aloct (:,jsize) = aloct(:,jwrap)
       aloct (:, 1   ) = aloct(:,jsize-jwrap+1)
       aloct (:, 2   ) = aloct(:,jsize-jwrap+2)
!
       colat (:,jsize) = colat(:,jwrap)
       colat (:, 1   ) = colat(:,jsize-jwrap+1)
       colat (:, 2   ) = colat(:,jsize-jwrap+2)
!
    CASE (4)
!
       isize = 80
       jsize = 51
       jwrap = 3
!      CALL Make_grid_arrays (isize, jsize)
!
       DO i = 1, isize
          colat (i,:) = REAL(i,rprec)*DTR
       END DO
!
       DO i = 1, isize
       DO j = jwrap, jsize-1
          aloct (i,j) = 0.5_rprec * REAL (j-jwrap, rprec)*HTR
       END DO
       END DO
!
       aloct (:,jsize) = aloct(:,jwrap)
       aloct (:, 1   ) = aloct(:,jsize-jwrap+1)
       aloct (:, 2   ) = aloct(:,jsize-jwrap+2)
!
    CASE (5)
!
       isize = 80*2
       jsize = 51
       jwrap = 3
!      CALL Make_grid_arrays (isize, jsize)
!
       DO i = 1, isize
          colat (i,:) = 0.5_rprec*REAL(i,rprec)*DTR
       END DO
!
       DO i = 1, isize
       DO j = jwrap, jsize-1
          aloct (i,j) = 0.5_rprec * REAL (j-jwrap, rprec)*HTR
       END DO
       END DO
!
       aloct (:,jsize) = aloct(:,jwrap)
       aloct (:, 1   ) = aloct(:,jsize-jwrap+1)
       aloct (:, 2   ) = aloct(:,jsize-jwrap+2)
!
    CASE DEFAULT
       STOP 'ILLEGAL TYPE OF GRID'
    END SELECT
!
!
    Re   = 6375.0_rprec        ! Earth's radius in km
    Ri   = Re + 100.0_rprec     ! radius of RCM's ionosphere in km
    dlam = 1.0_rprec / REAL(isize-1,rprec)
    dpsi = pi_two / REAL(jsize-jwrap,rprec)
!
    DO i = 1, isize
    DO j = 1, jsize
       IF (i == 1) THEN
          dthetadi = colat (i+1,j) - colat (i,j)
       ELSE IF (i == isize) THEN
          dthetadi = colat (i,j) - colat (i-1,j)
       ELSE
          dthetadi = 0.5_rprec * (colat (i+1,j) - colat (i-1,j))
       END IF
       alpha (i,j) = dthetadi / dlam
    END DO
    END DO
!
    CALL Wrap_around_ghostcells (alpha, isize, jsize, n_gc)
!
    vcorot = -92400.0_rprec*(Re / Ri)*SIN(colat)**2
    beta   = SIN(colat)
    sini   = 2.0_rprec*COS(colat)/SQRT(1.0_rprec+3.0_rprec*COS(colat)**2)
    bir    = 2.0_rprec*(Re / Ri)**3*besu*COS(colat)
!
!
!
!   Identify which grid points belong to the original MSM grid:
!
    ALLOCATE (msm_mask(isize,jsize), msm_i_index(isize,jsize), &
                                      msm_j_index(isize,jsize))
!
    DO i = 1, isize
    DO j = 1, jsize
!
       j_msm =  1
       logic_flag_1 = .FALSE.
       one_loop: DO i_msm = 1, isize_msm
          IF ( ABS(colat(i,j)-colat_msm(i_msm,j_msm)) <= &
               10.0_rprec*EPSILON(1.0_rprec)*colat(i,j) ) THEN
               logic_flag_1 = .TRUE.
               EXIT one_loop
          END IF
       END DO one_loop

       IF ( logic_flag_1 ) THEN
          logic_flag_2 = .FALSE.
          two_loop: DO j_msm = jwrap_msm, jsize_msm-1
             IF ( ABS(aloct(i,j)-aloct(i_msm,j_msm)) <= &
                  10.0_rprec*EPSILON(1.0_rprec)*aloct(i,j) ) THEN
                  logic_flag_2 = .TRUE.
                  EXIT two_loop
             END IF
          END DO two_loop
       END IF
!
       IF (logic_flag_1 .AND. logic_flag_2) THEN
          msm_mask (i,j) = .TRUE.
          msm_i_index (i,j) = i_msm
          msm_j_index (i,j) = j_msm
       ELSE
          msm_mask (i,j) = .FALSE.
          msm_i_index (i,j) = -1
          msm_j_index (i,j) = -1
       END IF
!
    END DO
    END DO
!
!
    CALL Write_grid (isize, jsize)
!
!
!   Write info about relation of the grid to the MSM grid:
!
    filename = 'rcmcrd_msm'
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//TRIM(filename), STATUS = 'REPLACE', &
          FORM = 'UNFORMATTED', IOSTAT = istat)
      IF (istat /= 0) then
         write(*,*)'ERROR OPENING ',trim(NameRcmDir),TRIM(filename)
         STOP
      end IF
      WRITE (LUN) isize_msm, jsize_msm, jwrap_msm
      WRITE (LUN) colat_msm
      WRITE (LUN) aloct_msm
    CLOSE (LUN)
!
    filename = 'rcmcrd_mask'
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//TRIM(filename), STATUS = 'REPLACE', &
          FORM = 'UNFORMATTED')
      WRITE (LUN) msm_mask
      WRITE (LUN) msm_i_index
      WRITE (LUN) msm_j_index
    CLOSE (LUN)
!
!   Make a formatted printout of the variables written to unformatted files:
!                                                        
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'setup_grid.list', STATUS = 'REPLACE')
!
    WRITE (LUN,*) 'ISIZE     ', isize
    WRITE (LUN,*) 'JSIZE     ', jsize
    WRITE (LUN,*) 'JWRAP     ', jwrap
    WRITE (LUN,*) 'DLAM      ', dlam
    WRITE (LUN,*) 'DPSI      ', dpsi
    WRITE (LUN,*) 'Ri        ', Ri
    WRITE (LUN,*) 'Re        ', Re
!
    ncol = 132
    CALL Outp (R = alpha, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'alpha', NTP = LUN, NCOL = ncol)
    CALL Outp (R = beta , ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'beta ', NTP = LUN, NCOL = ncol)
    CALL Outp (R = colat*RTD, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'colat in degrees', NTP=LUN, NCOL=ncol)
    CALL Outp (R = aloct, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'aloct', NTP = LUN, NCOL = ncol)
    CALL Outp (R = vcorot, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'vcorot', NTP = LUN, NCOL = ncol)
    CALL Outp (R = sini, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'sini', NTP = LUN,  NCOL = ncol)
    CALL Outp (R = bir, ISIZE=isize, JSIZE= jsize,IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'bir ', NTP = LUN, NCOL = ncol)
    CLOSE (UNIT = LUN)
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'rcmcrd_mask.list', STATUS = 'REPLACE')
    CALL Outp (R = msm_i_index, IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 1.0,       &
               ILABEL = label%intg, TITLE = 'msm_i_index ', NTP = LUN, NCOL = ncol)
!
    CALL Outp (R = msm_j_index, IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 1.0,       &
               ILABEL = label%intg, TITLE = 'msm_j_index ', NTP = LUN, NCOL = ncol)
!
    CALL Outp (R = msm_mask, IBEG = 1, IEND = isize, IINC = 1, JBEG = 1, &
               JEND = jsize, JINC = 1, XSCALE = 0.0,       &
               ILABEL = label%intg, TITLE = 'msm_mask ', NTP = LUN, NCOL = ncol)
!
    CLOSE (UNIT = LUN)
!
!
    STOP
!
    CONTAINS
!
    SUBROUTINE Make_msm_grid ()
    IMPLICIT NONE
    INTEGER (iprec) :: i, j
!
    colat_msm (1,:) = 1.000344_rprec * DTR
    DO i = 2, isize_msm
       colat_msm (i,:) = colat_msm(1,:) + &
                         Theta_of_i (i, dtheta_di_max, dtheta_di_min, &
                                     delta_i, i_pc, i_pp) * DTR
    END DO
!
    DO j = jwrap_msm, jsize_msm - 1
    DO i = 1, isize_msm                 ! teta is colatitude in degrees:
      aloct_msm (i,j) = 2*pi * REAL(j-jwrap_msm) / REAL(jsize_msm - jwrap_msm)
      IF (aloct_msm(i,j) < 0.) aloct_msm(i,j) = aloct_msm(i,j) + 2*pi
    END DO
    END DO
!
!
    aloct_msm (:,jsize_msm) = aloct_msm(:,jwrap_msm)
    aloct_msm (:, 1   ) = aloct_msm(:,jsize_msm-jwrap_msm+1)
    aloct_msm (:, 2   ) = aloct_msm(:,jsize_msm-jwrap_msm+2)
!
    RETURN
    END SUBROUTINE Make_msm_grid
!
!
!=========================================
!
    FUNCTION Theta_of_i (i, dtheta_di_max, dtheta_di_min, delta_i, i_pc, i_pp)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i, delta_i, i_pc, i_pp
    REAL (rprec), INTENT (IN) :: dtheta_di_max, dtheta_di_min
    REAL (rprec) :: Theta_of_i
!
!   MSM colatitudes in degrees as a function of I. Params set in calling sub.
!
    Theta_of_i =  (i-1) * dtheta_di_max + &
                  (delta_i / pi) * (dtheta_di_max - dtheta_di_min) * &
                  (-F_wolf ( REAL(i_pc-i) / REAL(delta_i) )+    &
                    F_wolf ( REAL(i-i_pp) / REAL(delta_i) )+    &
                    F_wolf ( REAL(i_pc-1) / REAL(delta_i) )-    &
                    F_wolf ( REAL(1-i_pp) / REAL(delta_i) )     &
                  )
    RETURN
    END FUNCTION Theta_of_i
!
!
!==========================================
!
    FUNCTION DTheta_di_of_i (i,dtheta_di_max, dtheta_di_min, delta_i, &
                             i_pc, i_pp)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: i, delta_i, i_pc, i_pp
    REAL (rprec), INTENT (IN) :: dtheta_di_max, dtheta_di_min
    REAL (rprec) :: DTheta_di_of_i
!
!   MSM's latitudinal grid spacing.
!
    DTheta_di_of_i = dtheta_di_max + (dtheta_di_max-dtheta_di_min) / pi * &
                     (ATAN (REAL(i_pc-i)/REAL(delta_i)) + &
                      ATAN (REAL(i-i_pp)/REAL(delta_i)) )
    RETURN
    END FUNCTION DTheta_di_of_i
!
!
!===========================================
!
    FUNCTION F_wolf (x)
       IMPLICIT NONE
       REAL (rprec), INTENT (IN) :: x
       REAL (rprec) :: F_wolf
       F_wolf = x*ATAN(x)-0.5*LOG(x**2+1)
       RETURN
    END FUNCTION F_wolf
!
!
!
!   SUBROUTINE Make_grid_arrays (isize,jsize)
!   IMPLICIT NONE
!   INTEGER (iprec), INTENT (IN) :: isize,jsize
!   ALLOCATE ( colat(1-n_gc:isize,1-n_gc:jsize), aloct (1-n_gc:isize,1-n_gc:jsize), &
!              alpha(1-n_gc:isize,1-n_gc:jsize), beta  (1-n_gc:isize,1-n_gc:jsize), &
!              sini (1-n_gc:isize,1-n_gc:jsize), bir   (1-n_gc:isize,1-n_gc:jsize), &
!              vcorot(1-n_gc:isize,1-n_gc:jsize),msm_mask (1-n_gc:isize,1-n_gc:jsize),&
!              msm_i_index(1-n_gc:isize,jsize), msm_j_index(1-n_gc:isize,1-n_gc:jsize),&
!              STAT = istat)
!   IF (istat /= 0) STOP 'FAILURE TO ALLOCATE'
!   RETURN
!   END SUBROUTINE Make_grid_arrays
!
    END PROGRAM Rcm_setup_grid
