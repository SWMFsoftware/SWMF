  PROGRAM Rcm_plasma_init_cond
  USE Rcm_variables
  USE Rcm_io
  IMPLICIT NONE
  INTEGER (iprec) :: itimei, itimef, irdr, irdw, idt, &
                         idt1, idt2, idt3, idebug

!   Plasma parameters:
  REAL :: N_i, kT_i, kT_e,&
          U_e, U_i, P_e, P_i, &
          eta_i_numeric, eta_e_numeric, &
          kT_e_numeric, kT_i_numeric,&
          N_e_numeric, N_i_numeric, &
          p_i_numeric, p_e_numeric, &
          U_i_numeric, U_e_numeric, &
          Eavg_i_numeric, Eavg_e_numeric, &
          pvgamma_edges, pvgamma_grid
!
!
  REAL :: p(kcsize), pc(kcsize)
!
!
!  Local variables:
      INTEGER :: Lrec, i, j, k, n_pt, nn, kc, &
                 option_eeta, option_edges, i_deepest, plasma_type, k_in
      CHARACTER (LEN=1) :: answer
      CHARACTER (LEN=14), DIMENSION(ksize)  :: remark_edg
      CHARACTER (LEN=14), DIMENSION(kcsize) :: remark_grd
      CHARACTER (LEN=12) :: char_date, char_time
      CHARACTER (LEN=50) :: char_label, title_string
      CHARACTER (LEN=80) :: form_string
      REAL :: vm_13, bi_13, x_pt, y_pt, x_1, y_1, x_2, y_2, &
              bi_initial_edges (1-n_gc:jsize+n_gc,ksize), bi_initial_grid (1-n_gc:jsize+n_gc,kcsize),&
              Gntrp_2d, eta_tmp, eta_i, eta_e, fraction
!
!
      OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'rcm.params', STATUS = 'OLD', &
            ACTION = 'READ', FORM = 'FORMATTED')
!
         READ (LUN,*) itimei   ! 1.  start time                          ! local
         READ (LUN,*) itimef   ! 2.  end time                            ! local
         READ (LUN,*) irdr     ! 3.  record # to read in                 ! local
         READ (LUN,*) irdw     ! 4.  record # to write out               ! local
         READ (LUN, '(a80)') label%char! 5.  text label                  ! in mod
         READ (LUN,*) idebug   ! 6.  0 <=> do disk printout              ! local
!
         READ (LUN,*) idt   !  1.  basic time step in program            ! local
         READ (LUN,*) idt1  !  2.  t-step for changing disk write rcds   ! local
         READ (LUN,*) idt2  !  3.  t-step for writing formatted output   ! local
         READ (LUN,*) idt3  !  4.  t-step for invoking ADD & ZAP         ! local
         READ (LUN,*) imin  !  5.  i-value of poleward bndy              ! in mod
         READ (LUN,*) ipot  !  6.  which potential solver to use         ! in mod
         READ (LUN,*) iwind !  9.  0 is no neutral winds                 ! in mod
         READ (LUN,*) ivoptn! 10.  1 is euler time derivs                ! in mod
!                                  in vcalc, 2--runge-kutta
         READ (LUN,*) ibnd_type  ! 14.  type of bndy (1-eq.p, 2-iono)
         READ (LUN,*) ipcp_type  ! 14.  type of bndy (1-eq.p, 2-iono)    ! in mod
         READ (LUN,*) nsmthi! 15.  How much to smooth cond in I          ! in mod
         READ (LUN,*) nsmthj! 16.  How much to smooth cond in J          ! in mod
         READ (LUN,*) icond ! 17.                                        ! in mod
         READ (LUN,*) ifloor! 18.                                        ! in mod
         READ (LUN,*) icorrect! 19.                                      ! in mod
!
         READ (LUN,*) dstmin ! min allowed dist. between adjacent pnts   ! in mod
         READ (LUN,*) dstmax ! max allowed dist. between adjacent pnts   ! in mod
         READ (LUN,*) rhomax                                             ! in mod
         READ (LUN,*) vmfact                                             ! in mod
         READ (LUN,*) epslon_edge  !  in MODULE rcm_mod_edgpt            ! in mod
         READ (LUN,*) fmrwif_dlim    ! min tail thickness for new FMRWIF ! in mod
         READ (LUN,*) cmax    ! in rcm_mod_balgn                         ! in mod
         READ (LUN,*) eeta_cutoff ! as a fraction
         READ (LUN,*) tol_gmres
         READ (LUN,*) itype_bf
         READ (LUN,*) i_advect
         READ (LUN,*) i_eta_bc
!
      CLOSE (UNIT = LUN)
!
      CALL Read_grid ()
!
!
      OPEN (UNIT=LUN, FILE=trim(NameRcmDir)//'alam_values.dat', &
           STATUS='OLD', ACTION='READ')
      DO k = 1, kcsize
         READ (LUN,*) k_in, alamc(k), ikflavc(k)
         IF (k_in /= k) STOP 'mismatch in k'
         IF (ikflavc(k) == 1 .AND. alamc(k) > 0.0) STOP 'ALAM>0 FOR ELECTRONS'
      END DO
      CLOSE (LUN)
!
      alam   = alamc
      ikflav = ikflavc
!
!
!
!     CALL Rcm_plasma_bc (i_eta_bc)
!
      DO k = 1, kcsize
         etac (k) = 0.0
      END DO
      eta = etac
!
!
      WRITE (*,'(A)',ADVANCE='NO') 'ENTER 1 FOR ION EDGES, 2 FOR GRID IONS: '
      READ (*,*) plasma_type
      IF (plasma_type == 1) THEN
         DO k = 1, kcsize
            IF (alamc(k) > 0.0) etac (k) = 0.0
         END DO
      ELSE IF (plasma_type == 2) THEN
         DO k = 1, kcsize
            IF (alamc(k) > 0.0) eta (k) = 0.0
         END DO
      ELSE
         STOP 'bad choice'
      END IF
      WRITE (*,'(a)',ADVANCE='NO') 'ELECTRONS WILL BE ON THE GRID (WITH EMPTY EDGES)'
      READ (*,*)
      DO k = 1, kcsize
         IF (alamc(k) < 0.0) eta (k) = 0.0
      END DO
!
      DO kc = 1, kcsize
           WRITE (remark_grd(kc),'(A,F6.3)') 'FUDGE=',fudgec(kc)
      END DO
      DO kc = 1, ksize
          IF (alam(kc) < 0.0) CYCLE
          IF (eta(kc) == 0.0) THEN
             WRITE (remark_edg(kc),'(A12)') 'TRACER ONLY'
          ELSE
             WRITE (remark_edg(kc),'(A12)') 'REAL EDGE'
          END IF
      END DO
      DO kc = 1, kcsize
           IF (alam(kc) > 0.0) CYCLE
           WRITE (remark_edg(kc),'(A12)') 'TRACER ONLY'
      END DO
!
      DO k = 1, kcsize
         CALL Wrap_around_ghostcells (eeta(:,:,k), isize, jsize, n_gc)
      END DO
!
!
      vm_13 = vm (imin_j(jsize/2),jsize/2)
!
!    7. Compute pressure for energy channels:
!        p = (2/3)*alam*eta*vm**(5/2)*2.51*e-34 [dynes/cm**2]
!
      DO k = 1, ksize
         IF (eta(k) == 0) THEN
            eta_tmp = etac(k)
         ELSE IF (etac(k) == 0.0) THEN
            eta_tmp = eta (k)
         ELSE
            WRITE (*,*) 'FOR K=',K, 'YOU HAVE NON-EMPTY EDGE AND GRID-BASED CHANNEL'
            STOP
         END IF
         p(k) = ABS(alam(k))*eta_tmp*vm_13**2.5*(2./3.)
         p(k) = p(k) * 2.51E-34
      END DO
!
      DO k = 1, kcsize
         pc(k) = ABS(alamc(k))*etac(k)*vm_13**2.5*(2./3.)
         pc(k) = pc(k) * 2.51E-34
      END DO
!
      pvgamma_edges = SUM ( ABS(alam)*eta   )
      pvgamma_grid  = SUM ( ABS(alamc)*etac )
!
!
GOTO 111
!
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    Compute plasma properties for the ions edges and electrons on grid:
!
!    ** U is energy density in erg/cm-3
!    ** Eavg is average particle density in eV: = U/density
!    ** kT is (2/3)*Eavg, in eV
!    ** pressure p = 2/3*U, in erg/cm-3
!
      eta_i_numeric = 0.0
      eta_e_numeric = 0.0
      U_i_numeric   = 0.0
      U_e_numeric   = 0.0
      Kt_i_numeric  = 0.0
      Kt_e_numeric  = 0.0
      DO k = 1, ksize
          IF (alam(k) > 0.0) THEN
             eta_i_numeric = eta_i_numeric + MAX(eta(k), etac(k))
             U_i_numeric   = U_i_numeric + 2.51E-34*(vm_13**2.5)* &
                             MAX(eta(k),etac(k))*ABS(alam(k))
             kT_i_numeric  = (2./3.)*MAX(eta(k),etac(k))*ABS(alam(k))*vm_13 + kT_i_numeric
          ELSE
             eta_e_numeric = eta_e_numeric + MAX(eta(k), etac(k))
             U_e_numeric   = U_e_numeric + 2.51E-34*(vm_13**2.5)* &
                             MAX(eta(k),etac(k))*ABS(alam(k))
             kT_e_numeric  = (2./3.)*MAX(eta(k),etac(k))*ABS(alam(k))*vm_13 + kT_e_numeric
          END IF
      END DO
      N_i_numeric   = eta_i_numeric * (vm_13**1.5)*1.E-9/Re/1.E+3/1.E+6
      N_e_numeric   = eta_e_numeric * (vm_13**1.5)*1.E-9/Re/1.E+3/1.E+6
      kT_i_numeric  = kT_i_numeric / eta_i_numeric ! eV
      kT_e_numeric  = kT_e_numeric / eta_e_numeric ! eV
!
      Eavg_i_numeric = 1.5 * kT_i_numeric
      Eavg_e_numeric = 1.5 * kT_e_numeric
      P_i_numeric    = (2./3.)*U_i_numeric
      P_e_numeric    = (2./3.)*U_e_numeric
!
      eta_i   = N_i / ((vm_13**1.5)*1.E-9/Re/1.E+3/1.E+6)
      eta_e   = N_i / ((vm_13**1.5)*1.E-9/Re/1.E+3/1.E+6)
!     Print summary plasma info (all info fits one 25-line screen)
!
      WRITE (6,'(T2,A)')  ' ION SETUP INFORMATION:'
      WRITE (6,'(T2,A)')  '-----------------------------------------'
      WRITE (6,'(T2,A)')  '-----------------------------------------'
      WRITE (6,'(T2,A)')  ' PARAM       |   ANALYTIC   |   NUMERIC  ' 
      WRITE (6,'(T2,A)')  '-----------------------------------------'
      WRITE (6, 801)      'N,  cm-3',  '|', N_i,      '|', N_i_numeric
      WRITE (6, 801)      'kT, eV',    '|', kT_i,     '|', kT_i_numeric
      WRITE (6, 801)      '<E>,eV',    '|', kT_i*1.5, '|', Eavg_i_numeric
      WRITE (6, 801)      'P, erg/cm3','|', P_i,      '|', P_i_numeric
      WRITE (6, 801)      'ETA, Wb-1', '|', eta_i,    '|', eta_i_numeric
      WRITE (6, 801)      'Pv**gamma ','|', 0.0 ,     '|', Pvgamma_edges
      WRITE (6,'(T2,A)')  '-----------------------------------------'
      WRITE (6,'(T2,A)')  '-----------------------------------------'
!
!
      WRITE (6,'(/,T2,A)') ' ELECTRON SETUP INFORMATION: '
      WRITE (6,'(T2,A)')   '-----------------------------------------'
      WRITE (6,'(T2,A)')   '-----------------------------------------'
      WRITE (6,'(T2,A)')   ' PARAM       |   ANALYTIC   |   NUMERIC  ' 
      WRITE (6,'(T2,A)')   '-----------------------------------------'
      WRITE (6, 801)       'N, cm-3',   '|', N_i,      '|', N_e_numeric
      WRITE (6, 801)       'kT, eV',    '|', kT_e,     '|', kT_e_numeric
      WRITE (6, 801)       '<E>, eV',   '|', kT_e*1.5, '|', Eavg_e_numeric
      WRITE (6, 801)       'P, erg/cm3','|', P_e,      '|', P_e_numeric
      WRITE (6, 801)       'ETA, Wb-1', '|', eta_e,    '|', eta_e_numeric
      WRITE (6, 801)       'Pv**gamma ','|', 0.0 ,     '|', Pvgamma_grid
      WRITE (6,'(T2,A)')     '-----------------------------------------'
      WRITE (6,'(T2,A)')     '-----------------------------------------'
 801  FORMAT (T2,A, T15,A1, ES9.2, T30,A1, ES9.2)
!
111 CONTINUE
!
!  6. Set up inner edges:
!
      DO k = 1, ksize  ! Initial location of edges:
         bi_initial_edges (:,k) = bndloc (:) + 1.0
      END DO
!
      n_pt = 0  ! point counter
      DO k = 1, ksize
        mpoint(k) = n_pt + 1
        Nn=0
        DO j = 1, jsize, 2
          n_pt         = n_pt + 1
          Nn           = Nn + 1
          bi(n_pt)     = bi_initial_edges(j,k)
          bj(n_pt)     = REAL(j,rprec)
          itrack(n_pt) = n_pt
          etab(n_pt)   = eta(k)
!
!         x_pt = Gntrp_2d (xmin, isize, jsize, bi(n_pt), bj(n_pt))
!         y_pt = Gntrp_2d (ymin, isize, jsize, bi(n_pt), bj(n_pt))
!         write (*,'(2(TR1,i4),2(TR2,F5.2),2(TR2,F8.2))') &
!               k, n_pt, bi(n_pt), bj(n_pt), x_pt, y_pt
        
        END DO
        npoint(k)      = Nn
      END DO
      mpoint(ksize+1) = n_pt+1
!
!
!
      WRITE (*,'(T2,A)',ADVANCE='NO') 'PRESS ENTER TO WRITE TO FILES, CTRL-C TO QUIT:'
      READ (*,*)
!
!
!
      CALL Write_plasma ()
!
!
!  Write inner edges information:
!
      label%intg = 0
      label%real = 0.0
      label%char = ''
      Lrec = 1
      label%intg (1) = Lrec
      CALL Write_array ('rcmbi_inp',     Lrec, Label, ARRAY_1D=bi,     SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmbj_inp',     Lrec, Label, ARRAY_1D=bj,     SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmetab_inp',   Lrec, Label, ARRAY_1D=etab,   SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmitrack_inp', Lrec, Label, ARRAY_1D=itrack, SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmmpoint_inp', Lrec, Label, ARRAY_1D=mpoint, SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmnpoint_inp', Lrec, Label, ARRAY_1D=npoint, SETUP = .TRUE., ASCI=asci_flag)
      CALL Write_array ('rcmeeta_inp',   Lrec, Label, ARRAY_3D=eeta,   SETUP = .TRUE., ASCI=asci_flag)
    CALL Write_array ('rcmetac_inp',     Lrec, label, ARRAY_1D=etac,   SETUP = .TRUE., asci=asci_flag)
!
STOP 'STOP HERE'
!
!=============================================================
!===  Print out a summary of the plasma setup: ===============
!
      OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//'setup_plasma.list', STATUS = 'REPLACE')
      CALL date_and_time (char_date, char_time)
!
!     A. Write a header:
!
      WRITE (LUN,901) 
      WRITE (LUN,901) 
      WRITE (LUN,'(/,T16,A)') 'SUMMARY OF RCM INITIAL PLASMA SETUP'
      WRITE (LUN,'(T16,A6,A2,A1,A2,A1,A2,TR5,A6,A2,A1,A2,A1,A2,/)') &
             'DATE: ', char_date(3:4),'/',char_date(5:6),'/',char_date(7:8),&
             'TIME: ', char_time(1:2),':',char_time(3:4),':',char_time(5:6)
!
!
!     B. Write out inner edges info:
!
      WRITE (LUN,'(A)') ' 1. PLASMA DESCRIBED by INNER EDGES:' 
      WRITE (LUN,901)
      WRITE (LUN, 903, ADVANCE = 'NO')
903   FORMAT (T3,'K', T10,'ALAM', T19,'W @-13Re', &
              T33,'ETA', T42,'P, %', T48,'Eta, %', T60,'Remark')
      WRITE (LUN, 902)
902   FORMAT (T6,'|',T17,'|',T28,'|',T39,'|',T47,'|', T55,'|', T71,'|')
      WRITE (LUN,901) 
!
      DO k = 1, ksize
         WRITE (LUN, 904, ADVANCE='NO') &
            k, alam(k), alam(k)*vm_13, eta(k), &
            ABS(alam(k))*eta(k)*vm_13**2.5*(2./3.)/SUM(p)*100.0, &
            eta(k)/eta_i_numeric*100.0, remark_edg(k)
         WRITE (LUN, 902)
      END DO
      WRITE (LUN, '(T28,ES10.2)', ADVANCE='NO') SUM(eta)
      WRITE (LUN, 902)
904   FORMAT (I3, T8,F8.1, T19,F8.1, T28,ES10.2, T40,F6.1, T49,F6.1, T58,A14)
      WRITE (LUN, 901) 
!
!
!     C. Write out grid plasma info:
!
      WRITE (LUN,'(A)') ' 2. PLASMA DESCRIBED ON THE GRID:' 
      WRITE (LUN,901)
      WRITE (LUN, 903, ADVANCE = 'NO')
      WRITE (LUN, 902)
      WRITE (LUN,901) 
!
      DO kc = 1, kcsize
         WRITE (LUN, 904, ADVANCE= 'NO') &
            kc, alamc (kc), alamc (kc)*vm_13, etac (kc), &
            pc(kc) / SUM(p)*100.,&
            etac (kc) / eta_e_numeric*100.,  remark_grd (kc)
         WRITE (LUN, 902)
      END DO
      WRITE (LUN, 901)
!
!
      WRITE (LUN,'(T2,A)')  ' ION SETUP INFORMATION'
      WRITE (LUN,'(T2,A)')  '-----------------------------------------'
      WRITE (LUN,'(T2,A)')  '-----------------------------------------'
      WRITE (LUN,'(T2,A)')  ' PARAM       |   ANALYTIC   |   NUMERIC  '
      WRITE (LUN,'(T2,A)')  '-----------------------------------------'
      WRITE (LUN, 801)      'N,  cm-3',  '|', N_i,      '|', N_i_numeric
      WRITE (LUN, 801)      'kT, eV',    '|', kT_i,     '|', kT_i_numeric
      WRITE (LUN, 801)      '<E>,eV',    '|', kT_i*1.5, '|', Eavg_i_numeric
      WRITE (LUN, 801)      'P, erg/cm3','|', P_i,      '|', P_i_numeric
      WRITE (LUN, 801)      'ETA, Wb-1', '|', eta_i,    '|', eta_i_numeric
      WRITE (LUN, 801)      'Pv**gamma ','|', 0.0 ,     '|', Pvgamma_edges
      WRITE (LUN,'(T2,A)')  '-----------------------------------------'
      WRITE (LUN,'(T2,A)')  '-----------------------------------------'
!
!
      WRITE (LUN,'(/,T2,A)') ' ELECTRON SETUP INFORMATION: '
      WRITE (LUN,'(T2,A)')   '-----------------------------------------'
      WRITE (LUN,'(T2,A)')   '-----------------------------------------'
      WRITE (LUN,'(T2,A)')   ' PARAM       |   ANALYTIC   |   NUMERIC  '
      WRITE (LUN,'(T2,A)')   '-----------------------------------------'
      WRITE (LUN, 801)       'N, cm-3',   '|', N_i,      '|', N_e_numeric
      WRITE (LUN, 801)       'kT, eV',    '|', kT_e,     '|', kT_e_numeric
      WRITE (LUN, 801)       '<E>, eV',   '|', kT_e*1.5, '|', Eavg_e_numeric
      WRITE (LUN, 801)       'P, erg/cm3','|', P_e,      '|', P_e_numeric
      WRITE (LUN, 801)       'ETA, Wb-1', '|', eta_e,    '|', eta_e_numeric
      WRITE (LUN, 801)       'Pv**gamma ','|', 0.0 ,     '|', Pvgamma_grid
      WRITE (LUN,'(T2,A)')     '-----------------------------------------'
      WRITE (LUN,'(T2,A)')     '-----------------------------------------'
      WRITE (LUN,*)
      WRITE (LUN,*)
      WRITE (LUN,'(T10,A//)') 'END OF RCM INTIAL PLASMA SETUP SUMMARY'
!
      WRITE (LUN,901)
      WRITE (LUN,901)
901   FORMAT ('------------------------------------------------------------')
!
!
!
!
      CLOSE (UNIT = LUN)
!
      STOP 
  END PROGRAM Rcm_plasma_init_cond
