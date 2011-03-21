MODULE Rcm_variables

  ! SWMF: use the temporary unit number provided by ModIoUnit
  use ModIoUnit, ONLY : UNITTMP_, STDOUT_

    IMPLICIT NONE
    SAVE
!
    LOGICAL :: IsPartofFramework = .true.
!
    INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND (9)
    INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND (12,100)
!
    ! SWMF: changed the value of LUN from 11 to UNITTMP_
    INTEGER, PARAMETER :: LUN = UNITTMP_

    ! SWMF: parameter attribute removed for LUN_2 and LUN_3
    ! SWMF: logical units are now obtained from ModIoUnit::io_unit_new()
    INTEGER :: LUN_2 = 12, LUN_3 = 13, LUN_4=14
!
    ! SWMF: added variable directory name for flexible IO
    character(len=100) :: NameRcmDir='IM/'
!
    ! SWMF: added restart logical, time step, and plot frequency
    logical :: DoRestart=.false.
    integer :: iDtRcm=5, nFilesPlot=0, iDnPlot(9)=-1, iDtPlot(9)=-1
    character(len=3) :: plot_area(9), plot_var(9), plot_format(9)
!
    ! SWMF: added unit number and prefix for output
    integer          :: iUnitOut=STDOUT_
    character(len=3) :: StringPrefix='IM:'

    ! SWMF: added logical to check for initialized grid
    logical :: IsUninitialized = .true., DoneGmCoupling = .false.
    ! GM-IM coupling is multifluid passing
     logical :: DoMultiFluidGMCoupling = .false.

    ! SWMF: added values to control satellite tracing.
    logical :: DoWriteSats  = .false.  !!!DTW 2007
    logical :: IsFirstWrite = .true.
    integer :: nImSats = 0, iStartIter = 0
    real,               allocatable :: SatLoc_3I(:,:,:)
    character(len=100), allocatable :: NameSat_I(:)

    ! SWMF: set initial values for composition and charge exchange information
    real(rprec) :: x_h=0.70, x_o=1.00-0.70
    logical     :: L_dktime = .TRUE.
    real(rprec) :: sunspot_number=125., f107=169., doy=90.
!          comment: allow doy to be floating in case will want to use fractionals
!                   later (stanislav, 6/15/2003).

    ! SWMF: set exponential decay of ring current
    logical     :: UseDecay = .FALSE.
    real(rprec) :: DecayTimescale = 36000.

    REAL (rprec), PARAMETER ::     &
         pi           = 3.141592654_rprec, &
         pi_two       = 2.0_rprec * pi, &
         pi_by_two    = pi / 2.0_rprec, &
         rtd          = 180.0_rprec/pi, &
         dtr          = pi/180.0_rprec, &
         rth          = 12.0_rprec / pi,&
         htr          = 1.0_rprec / rth      
!
    INTEGER (iprec), PARAMETER :: &
       n_gc  = 2,  &
       isize =  78, &
       jsize =  48, &
       iesize =   3, &
       kmin (iesize) = (/1,31, 116/), &
       kmax (iesize) = (/30,115, 200/), &
       ksize = 200, kcsize = ksize, &
       nptmax = 50000, &
       ncoeff =   5
    LOGICAL :: asci_flag = .TRUE.
!
    REAL (rprec), PARAMETER :: &
         xmass (iesize) = (/ 9.1E-31_rprec, 1.67E-27_rprec, 16*1.67E-27_rprec /), &
         besu           = 3.0584E+4_rprec, &
         signbe         = 1.0_rprec, &
         romeca         = 0.0_rprec, &
         charge_e       = 1.6E-19_rprec, &
         sgn (ksize)    = 1.0_rprec
!
    INTEGER (iprec), PARAMETER :: ie_ele = 1, ie_prt = 2, ie_oxg = 3
    CHARACTER (LEN=2), DIMENSION(3), PARAMETER :: species_char = &
                          (/ 'e-', 'H+', 'O+' /)
!
    TYPE :: label_def    !   Definition of the label structure, for I/O:
       INTEGER (iprec)   :: intg (20)
       REAL (rprec)      :: real (20)
       CHARACTER(LEN=80) :: char
    END TYPE label_def
!
    TYPE (label_def) :: label
!
!
    TYPE :: ellipse_def   ! Definition of an ellipse:
       REAL(rprec) :: aa, bb, xx, yy
    END TYPE ellipse_def
!
!
!   Grid info:
    REAL (rprec) :: dlam, dpsi, Ri, Re, &
                    alpha (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    beta  (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    colat (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    aloct (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    bir   (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    sini  (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    vcorot(1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    fac   (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc)
    INTEGER (iprec) :: i1, i2, iint, j1, j2, jint, imin, &
           imin_j(1-n_gc:jsize+n_gc)=isize, ibnd_type=4
    TYPE (ellipse_def), DIMENSION (2) :: boundary
!
!
!   Plasma on inner edges:
    INTEGER (iprec) :: ikflav (ksize), mpoint (ksize+1), npoint(ksize), &
                       itrack(nptmax)=-999, ivoptn
    REAL (rprec) :: bi (nptmax), bj (nptmax), etab (nptmax), alam (ksize), &
                    eta (ksize), dbidt (nptmax), dbjdt (nptmax), &
                    fudge (ksize), &
                    fmrwif_dlim, dstmin, dstmax, epslon_edge, rhomax, vmfact
!
!
!   Plasma on grid:
    REAL (rprec) :: alamc (kcsize)=0.0, etac (kcsize), fudgec (kcsize), &
                    eeta_cutoff, cmax, &
                    precipitation_tau(iesize)=(/0.3,0.0,0.0/), &
                    density     (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    pressure    (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    temperature (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    densityHp   (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    densityOp   (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    pressureHp  (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    pressureOp  (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    temperatureHp(1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), &
                    temperatureOp(1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc)


! SWMF: allocatable array (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,kcsize)
    REAL (rprec), allocatable:: eeta(:,:,:) !, eeta_avg(:,:,:) never used


    INTEGER (iprec) :: ikflavc (kcsize), i_advect, i_eta_bc
    INTEGER (iprec), PARAMETER :: irdk=18, inrgdk=13, isodk=2, iondk=2
    REAL (rprec) :: dktime (irdk, inrgdk, isodk, iondk)
    REAL (rprec), DIMENSION (19,5) :: trf !plasmasphere refilling rates, cm-3/day
!
!
!   Magnetic field:
    REAL (rprec), DIMENSION (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc) :: &
                    xmin, ymin, &
                    bmin, vm=0.0, &
                    rmin=0.0, pmin,&
                    xmin_1, xmin_2, &
                    ymin_1, ymin_2, &
                    bmin_1, bmin_2, &
                    vm_1, vm_2, &
                    bndloc (1-n_gc:jsize+n_gc)=isize
    INTEGER (iprec), ALLOCATABLE :: ibtime (:)
    REAL    (rprec) :: fstoff=0.0, fclps=0.0, fdst=0.0, fmeb=0.0, ftilt=0.0
    INTEGER (iprec) :: itype_bf
!
!
!   Ionospheric quantities:
    REAL (rprec), DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc) :: &
                    qtped, pedpsi, &
                    qtplam, pedlam, &
                    qthall, hall, &
                    ss (1-n_gc:jsize+n_gc), &
                    pwe, pwn, &
                    hwe, hwn , &
                    sw  (1-n_gc:jsize+n_gc), &
                    eflux (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,iesize), &
                    eavg (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,iesize),&
                    c_pde (ncoeff, 1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc), c5w 
!
!
!------------------------------------------------------------------------------
!   Potential solver variables to hold the linear system:
!
    INTEGER (iprec), PARAMETER :: ni_pde = isize , &
                                  nj_pde = jsize , &
                                  nij_pde = ni_pde * nj_pde, &
                                  nzmax_pde = ni_pde * nj_pde * ncoeff, &
                                  m_gmres   = 300
    REAL (rprec) ::               a_mtrx_pde (nzmax_pde), &
                                  b_mtrx_pde (nij_pde), &
                                  pivots_pde (nij_pde), &
                                  x0_pde     (nij_pde)
    INTEGER (iprec) ::            i_column_pde (nzmax_pde), &
                                  row_ptr_pde (nij_pde+1), &
                                  diag_ptr_pde (nij_pde)
    REAL (rprec) :: tol_gmres, & ! Potential solver GMRESM tolerance:
                    tol_bicgstab
    INTEGER (iprec), PARAMETER :: iter_max_bicgstab = 10000, &
                                  iter_max_gmres   = 300
    INTEGER (iprec) :: iter_bicgstab, iter_gmresm
    INTEGER (iprec) :: iTimeT1=0
!------------------------------------------------------------------------------
!
!
    INTEGER (iprec) :: icond, nsmthi, nsmthj, iwind
    LOGICAL :: ifloor, icorrect
!
!
!   Magnetospheric quantities:
    REAL (rprec), DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc) :: &
                    v, vpar, vbnd (1-n_gc:jsize+n_gc), &
                    birk, pvgamma (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,iesize), &
                    pressrcm, &
                    v_avg, birk_avg, birk_mhd, sigmaH_mhd,sigmaP_mhd
    INTEGER (iprec) :: ipcp_type=11, ipot=6
!
!
!   Input PCP drop and its current value:
    INTEGER (iprec), ALLOCATABLE :: ivtime (:)
    REAL    (rprec), ALLOCATABLE :: vinput (:), vinput_phase(:)
    REAL    (rprec)              :: vdrop,      vdrop_phase
!
!
!   Input Kp values and its current value:
    INTEGER (iprec), ALLOCATABLE :: ikptime (:)
    REAL    (rprec), ALLOCATABLE :: Kpinput (:)
    REAL    (rprec)              :: Kp
!
!
!   Input ETAC values:
    INTEGER (iprec), ALLOCATABLE :: itime_etac (:)
    REAL (rprec),    ALLOCATABLE :: etac_inp(:,:)

! Timing variables
    REAL (rprec) :: time0,time1

! - New variables for a future parallelized RCM:
! 
!   integer, parameter :: n_blk = 1, &
!   REAL (rprec), DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc,n_blk) :: &
!                 v_blk, birk_blk, &
!                 alpha_blk, beta_blk, aloct_blk, colat_blk, vcorot_blk, &
!                 sini_blk, bir_blk, &
!                 qtpedlam_blk, qtped_blk, qthall_blk, pedlam_blk, pedpsi_blk,&
!                 hall_blk, c5w_blk 
! 
! SWMF: allocatable arrays(1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc,kcsize)
!    REAL (rprec), allocatable :: eeta_blk(:,:,:)
!
END MODULE Rcm_variables
