Module ModMain
  !\
  ! Get values of nBLK, nI, nJ, nK, and gcn from ModSize
  !/
  use ModSize

  character (len=10) :: CodeVersion='6.07'

  !\
  ! Block parameters.
  !/
  integer :: nMultiBlkLevels, globalBLK

  !\
  ! Problem definition.
  !/
  integer :: problem_type
  integer, parameter :: problem_uniform    =1, &
                        problem_shocktube  =2, &
                        problem_heliosphere=3, &
                        problem_comet      =5, &
                        problem_rotation   =6, &
                        problem_diffusion  =7, &
                        problem_earth      =11,&
                        problem_saturn     =12,&
                        problem_jupiter    =13,&
                        problem_venus      =14,&
                        problem_cylinder   =21,&
                        problem_sphere     =22

  integer :: fluxfcn_type
  integer, parameter :: fluxfcn_roe=1, fluxfcn_rusanov=2, fluxfcn_linde=3

  ! Radius where we switch between conservative E and non-conservative P update
  real :: Rconservative

  ! Index parameters for conservative variables
  integer, parameter:: rho_=1, rhoU_=1, rhoUx_=2, rhoUy_=3, rhoUz_=4, &
                       B_=4, Bx_=5, By_=6, Bz_=7, E_=8, P_=9

  ! Primitive variable names
  integer, parameter :: U_=rhoU_, Ux_=rhoUx_, Uy_=rhoUy_, Uz_=rhoUz_

  ! Primitive variables in local r,phi,theta coordinate system
  integer, parameter :: Ur_=Ux_, Uph_=Uy_, Uth_=Uz_, &
                        Br_=Bx_, Bph_=By_, Bth_=Bz_

  !\
  ! Face array parameters.
  !/
  integer, parameter :: nFacesIx=nI+1, &
                        nFacesJx=nJ,   &
                        nFacesKx=nK
  integer, parameter :: nFacesIy=nI,   &
                        nFacesJy=nJ+1, &
                        nFacesKy=nK
  integer, parameter :: nFacesIz=nI,   &
                        nFacesJz=nJ,   &
                        nFacesKz=nK+1
  
  !\
  ! Convergence parameters.
  !/
  real, parameter :: rTOL=1.e-10
  real, parameter :: eps=.00001 ! (a small number)
  
  !\
  ! Geometry parameters.
  !/
  real  ::    x1, x2, y1, y2, z1, z2

  logical :: unused

  real :: dxyz(3),xyzend(4,3)
  
  !\
  ! Time stepping parameters and values.
  !/
  integer :: n_step, nITER, nORDER, nSTAGE, iteration_number=0
  real :: t_max=-1, dt, cfl, dt_BLK(nBLK), cputime_max
  logical :: time_accurate, point_implicit, boris_correction, &
             time_loop, check_stopfile

  !\
  ! Model Coupling variables
  !/

  logical :: UseIonosphere
  logical :: UseUAM = .false.
  logical :: UAMtoMHD =.false.
  logical :: MHDtoUAM = .false.
  logical :: UseIMM = .false.
  logical :: UseAMIE = .false.
  logical :: UseRAM = .false.
  logical :: UseFakeRegion2 = .false.

  integer :: dn_couple_ionosphere
  real    :: dt_couple_ionosphere
  integer, parameter :: iono_init=1, iono_fac=2, iono_read=3, iono_save=4, &
                        iono_save_restart=5, iono_solve=6

  integer :: dn_raytrace=100
  logical :: check_rayloop=.false.

  logical :: SetDipoleTilt
  real    :: ThetaTiltDeg
  real    :: dt_UpdateB0

  !\
  ! Refinement parameters.
  !/
  integer :: initial_refine_levels
  integer :: dn_refine
  integer :: min_block_level, max_block_level
  real :: min_cell_dx, max_cell_dx
  logical :: automatic_refinement, fix_body_level

  !\
  ! Load balance parameters
  !/
  integer, parameter :: maxloadbalancelevels = 32
  integer, dimension(maxloadbalancelevels) :: loadbalance
 
  !\
  ! Block cell-centered MHD solution array definitions.
  !/
  real,  dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK) :: &
    x_BLK,y_BLK,z_BLK,R_BLK,R2_BLK, &                !Block cell coordinates
    rho_BLK,rhoUx_BLK,rhoUy_BLK,rhoUz_BLK, &         !Block cell state
    Bx_BLK,By_BLK,Bz_BLK,E_BLK,p_BLK, &           
    B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, &           !Block intrinsic field
    time_BLK, &                                      !Block time 
    tmp1_BLK, tmp2_BLK                               !Temporary block variables
  common /MHDSolnBlock1/ &
    x_BLK,y_BLK,z_BLK,R_BLK,R2_BLK, &
    rho_BLK,rhoUx_BLK,rhoUy_BLK,rhoUz_BLK, & 
    Bx_BLK,By_BLK,Bz_BLK,E_BLK,p_BLK, &           
    B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, & 
    time_BLK, &
    tmp1_BLK, tmp2_BLK

  !Block cell old state, body forces & heat sources
  real,  dimension(1:nI, 1:nJ, 1:nK,nBLK) :: &
    rho_o_BLK,rhoUx_o_BLK,rhoUy_o_BLK,rhoUz_o_BLK, &
    Bx_o_BLK,By_o_BLK,Bz_o_BLK,E_o_BLK,p_o_BLK, &
    fbody_x_BLK,fbody_y_BLK,fbody_z_BLK,qheat_BLK
 
  common /MHDSolnBlock2/ &
    rho_o_BLK,rhoUx_o_BLK,rhoUy_o_BLK,rhoUz_o_BLK, & 
    Bx_o_BLK,By_o_BLK,Bz_o_BLK,E_o_BLK,p_o_BLK, &
    fbody_x_BLK,fbody_y_BLK,fbody_z_BLK,qheat_BLK  

  ! Shocktube initial state values
  real, dimension(rho_:P_) :: shock_Lstate, shock_Rstate
  logical :: shock_rotated

  ! Local cell source terms and divB.
  real, dimension(1:nI,1:nJ,1:nK) :: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,SE,SP,Theat0
  real, dimension(0:nI+1,0:nJ+1,0:nK+1):: Sdivb

  real, dimension( 0:nI+1, 0:nJ+1, 0:nK+1) :: &
      gradX_rho,gradX_Ux,gradX_Uy,gradX_Uz,gradX_Bx,gradX_By,gradX_Bz,gradX_p,&
      gradX_VAR,&
      gradY_rho,gradY_Ux,gradY_Uy,gradY_Uz,gradY_Bx,gradY_By,gradY_Bz,gradY_p,&
      gradY_VAR,&
      gradZ_rho,gradZ_Ux,gradZ_Uy,gradZ_Uz,gradZ_Bx,gradZ_By,gradZ_Bz,gradZ_p,&
      gradZ_VAR
  
  !\
  ! Block face-centered intrinsic magnetic field array definitions.
  !/
  real,  dimension(0:nFacesIx+1,0:nFacesJx+1,0:nFacesKx+1,nBLK) :: &
    B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK 
  
  real,  dimension(0:nFacesIy+1,0:nFacesJy+1,0:nFacesKy+1,nBLK) :: &
    B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK
  
  real,  dimension(0:nFacesIz+1,0:nFacesJz+1,0:nFacesKz+1,nBLK) :: &
    B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK

  !\
  ! Other block solution and geometry parameters.
  !/
  real, dimension(nBLK) :: dx_BLK, dy_BLK, dz_BLK, Rmin_BLK, Rmin2_BLK, &
                           fAx_BLK, fAy_BLK, fAz_BLK, cV_BLK, cVinv_BLK

  logical, dimension(nBLK) :: unusedBLK

  integer, dimension(nBLK) :: global_block_number

  !\
  ! X Face local MHD solution array definitions.
  !/
  real, dimension(0:nFacesIx+1,0:nFacesJx+1,0:nFacesKx+1) ::     &
    rhoFaceL_x,UxFaceL_x,UyFaceL_x,UzFaceL_x,              & !\
    BxFaceL_x,ByFaceL_x,BzFaceL_x,pFaceL_x,                & ! Face Left X
    rhoFaceR_x,UxFaceR_x,UyFaceR_x,UzFaceR_x,              & !\
    BxFaceR_x,ByFaceR_x,BzFaceR_x,pFaceR_x,                & ! Face Right X
    rhoFaceF_x,rhoUxFaceF_x,rhoUyFaceF_x,rhoUzFaceF_x,     & !\
    BxFaceF_x,ByFaceF_x,BzFaceF_x,EFaceF_x,PFaceF_x,       & ! Face Flux X
    VdtFace_x                                                !V/dt Face X
  common /MHDSolnFaceX/ &
    rhoFaceL_x,UxFaceL_x,UyFaceL_x,UzFaceL_x,              &
    BxFaceL_x,ByFaceL_x,BzFaceL_x,pFaceL_x,                &
    rhoFaceR_x,UxFaceR_x,UyFaceR_x,UzFaceR_x,              &
    BxFaceR_x,ByFaceR_x,BzFaceR_x,pFaceR_x,                &
    rhoFaceF_x,rhoUxFaceF_x,rhoUyFaceF_x,rhoUzFaceF_x,     &
    BxFaceF_x,ByFaceF_x,BzFaceF_x,EFaceF_x,PFaceF_x,       &
    VdtFace_x

  !\
  ! X Face conservative or corrected flux.
  !/
  real, dimension(1:nFacesJx, &
                  1:nFacesKx, &
                  1:2,nBLK) :: &
    rhoFaceFC_x_BLK, &
    rhoUxFaceFC_x_BLK,rhoUyFaceFC_x_BLK,rhoUzFaceFC_x_BLK, &
    BxFaceFC_x_BLK,ByFaceFC_x_BLK,BzFaceFC_x_BLK, &
    EFaceFC_x_BLK,BxFaceAC_x_BLK

  !\
  ! Y Face local MHD solution array definitions.
  !/
  real, dimension(0:nFacesIy+1,0:nFacesJy+1,0:nFacesKy+1) ::     &
    rhoFaceL_y,UxFaceL_y,UyFaceL_y,UzFaceL_y,              & !\
    BxFaceL_y,ByFaceL_y,BzFaceL_y,pFaceL_y,                & ! Face Left Y
    rhoFaceR_y,UxFaceR_y,UyFaceR_y,UzFaceR_y,              & !\
    BxFaceR_y,ByFaceR_y,BzFaceR_y,pFaceR_y,                & ! Face Right Y
    rhoFaceF_y,rhoUxFaceF_y,rhoUyFaceF_y,rhoUzFaceF_y,     & !\
    BxFaceF_y,ByFaceF_y,BzFaceF_y,EFaceF_y,PFaceF_y,       & ! Face Flux Y
    VdtFace_y                                                !V/dt Face Y
  common /MHDSolnFaceY/ &
    rhoFaceL_y,UxFaceL_y,UyFaceL_y,UzFaceL_y,              &
    BxFaceL_y,ByFaceL_y,BzFaceL_y,pFaceL_y,                &
    rhoFaceR_y,UxFaceR_y,UyFaceR_y,UzFaceR_y,              &
    BxFaceR_y,ByFaceR_y,BzFaceR_y,pFaceR_y,                &
    rhoFaceF_y,rhoUxFaceF_y,rhoUyFaceF_y,rhoUzFaceF_y,     &
    BxFaceF_y,ByFaceF_y,BzFaceF_y,EFaceF_y,PFaceF_y,       &
    VdtFace_y                                               

  !\
  ! Y Face conservative or corrected flux.
  !/
  real, dimension(1:nFacesIy, &
                  1:nFacesKy, &
                  1:2,nBLK) :: &
    rhoFaceFC_y_BLK, &
    rhoUxFaceFC_y_BLK,rhoUyFaceFC_y_BLK,rhoUzFaceFC_y_BLK, &
    BxFaceFC_y_BLK,ByFaceFC_y_BLK,BzFaceFC_y_BLK, &
    EFaceFC_y_BLK,ByFaceAC_y_BLK

  !\
  ! Z Face local MHD solution array definitions.
  !/
  real, dimension(0:nFacesIz+1,0:nFacesJz+1,0:nFacesKz+1) ::     &
    rhoFaceL_z,UxFaceL_z,UyFaceL_z,UzFaceL_z,              & !\
    BxFaceL_z,ByFaceL_z,BzFaceL_z,pFaceL_z,                & ! Face Left Z
    rhoFaceR_z,UxFaceR_z,UyFaceR_z,UzFaceR_z,              & !\
    BxFaceR_z,ByFaceR_z,BzFaceR_z,pFaceR_z,                & ! Face Right Z
    rhoFaceF_z,rhoUxFaceF_z,rhoUyFaceF_z,rhoUzFaceF_z,     & !\
    BxFaceF_z,ByFaceF_z,BzFaceF_z,EFaceF_z,PFaceF_z,       & ! Face Flux Z
    VdtFace_z                                                !V/dt Face Z
  common /MHDSolnFaceZ/ &
    rhoFaceL_z,UxFaceL_z,UyFaceL_z,UzFaceL_z,              &
    BxFaceL_z,ByFaceL_z,BzFaceL_z,pFaceL_z,                &
    rhoFaceR_z,UxFaceR_z,UyFaceR_z,UzFaceR_z,              &
    BxFaceR_z,ByFaceR_z,BzFaceR_z,pFaceR_z,                &
    rhoFaceF_z,rhoUxFaceF_z,rhoUyFaceF_z,rhoUzFaceF_z,     &
    BxFaceF_z,ByFaceF_z,BzFaceF_z,EFaceF_z,PFaceF_z,       &
    VdtFace_z

  !\
  ! Z Face conservative or corrected flux.
  !/
  real, dimension(1:nFacesIz, &
                  1:nFacesJz, &
                  1:2,nBLK) :: &
    rhoFaceFC_z_BLK, &
    rhoUxFaceFC_z_BLK,rhoUyFaceFC_z_BLK,rhoUzFaceFC_z_BLK, &
    BxFaceFC_z_BLK,ByFaceFC_z_BLK,BzFaceFC_z_BLK, &
    EFaceFC_z_BLK,BzFaceAC_z_BLK

  ! Parameters to hold direction names

  integer, parameter:: east_=1, west_=2, south_=3, north_=4, bot_=5, top_=6

  ! Inner and Outer boundary conditions

  character (len=20) :: innerBCtype, outerBCtype(east_:top_)

  !\
  ! Parallel AMR:
  ! Neighbor solution block refinement levels
  ! ( 0=neighbors at same level, 
  !  -1=neighbors at lower level,
  !  +1=neighbors at higher level,
  !  NOBLK=no neighbors).
  !/
  integer, parameter :: NOBLK=-100

  integer, dimension(nBLK) :: &
       neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth

  integer, dimension(east_:top_,nBLK):: neiLEV
  
  !\
  ! Parallel AMR:
  ! Neighbor processor and block numbers (a value of NOBLK
  ! means not used).  As only one level change is permitted
  ! between neighboring solution blocks, there are either 1 or 4 
  ! neighboring blocks in each of the six directions.
  !/
  integer, dimension(4,nBLK) :: &
       neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth, &
       neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth

  integer, dimension(4,east_:top_,nBLK) :: neiPE, neiBLK
  
  integer, dimension( -1:1, -1:1, -1:1, 4, nBLK) :: &
     BLKneighborPE, BLKneighborBLK, BLKneighborCHILD
  integer, dimension( -1:1, -1:1, -1:1, nBLK) :: BLKneighborLEV
  
  !\
  ! Variables for nonblocking sends used in fix_refine and fix_coarsen.
  !/
  integer, parameter    :: NO_BLK=-1, &
                           NUM_NB_BUFFS=12
  
  integer :: nb_buffs(3,NUM_NB_BUFFS,0:nBLK)
  integer :: curr_nb_buff(0:nBLK)
  integer :: nb_buff_req(NUM_NB_BUFFS,0:nBLK)
  
  ! Parameters for block location among eight subcubes.
  integer, parameter    :: LOC_TSE=1, &
                           LOC_TSW=2, &
                           LOC_BSW=3, &
                           LOC_BSE=4, &
                           LOC_BNE=5, &
                           LOC_BNW=6, &
                           LOC_TNW=7, &
                           LOC_TNE=8
  
  ! Use corners/edges in gradients and message passing
  logical :: UseCorners

  ! Use raytracing
  logical :: UseRaytrace=.false.

  ! How to deal with div B = 0
  logical :: UseDivbSource, UseDivbDiffusion, UseProjection, UseConstrainB

  ! Update check parameters
  logical :: UseUpdateCheck
  real :: percent_max_rho(2), percent_max_p(2)

  ! Default value for diffusion coefficient. Should be stable for CFL<1.
  ! Set 11/7/2000 - criteria defined/choosen by Gabor Toth
  real :: divb_diffcoeff=1.0/6.0

  ! Choice of limiter
  character*6 :: limiter_type='minmod'

  ! Prolongation order (1 or 2) and type ('central','lr','minmod','central2'..)
  integer :: prolong_order=1
  character*10 :: prolong_type='lr'

  ! Variables describing cells inside boundaries

  logical :: body_BLK(nBLK)     ! true when the block (including ghost cells)
                                !   is intersected by any body 
  logical :: true_BLK(nBLK)     ! true when all cells in block 
                                !   (not including ghost cells) are true_cells 
  logical :: true_cell(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK)
                                ! true cells are cells that are not ghost 
                                !   cells out side or inside a body

  logical :: far_field_BCs_BLK(nBLK)

  ! Logical for bodies
  logical :: body1=.false., body2=.false.

  ! Logical for corotation
  logical :: UseCorotation=.false.
  logical :: AlignDipoleCorotation=.true.
  logical :: SetCorotationTilt=.false.
  real    :: CorotationTilt=0.0, CorotationTiltDeg=0.0
  real    :: CorotationLon=0.0, CorotationLonDeg=0.0
  real    :: Magnetic_Pole_Colat, Magnetic_Pole_East_Lon
  real    :: Max_DoY_Tilt
  integer :: CoordinateSystem

  ! Logical for mass loading
  logical :: UseMassLoading=.false.
  logical :: AccelerateMassLoading=.false.

  ! String variable for debugging. A space separated list of words,
  ! typically names of subroutines. Listing "QUIET" suppresses progress reports

  character (len=79) :: test_string

  ! Location for test

  integer :: Itest=1, Jtest=1, Ktest=1, BLKtest=1, PROCtest=0, ITERtest=-1
  integer :: VARtest=1, DIMtest=1
  real    :: Xtest, Ytest, Ztest, Ttest
  real    :: Xtest_mod, Ytest_mod, Ztest_mod
  logical :: coord_test=.false.

  !\
  ! Debug parameters. 
  !/
  logical, parameter :: debug1=.false., debug2=.false.

  logical :: okdebug=.false., ShowGhostCells=.true.

  ! My processor number and number of processors

  integer :: me_world, numprocs

  ! error codes, so we do not need to define them in every subroutine

  integer :: ira, erno

  ! Time variables
  real :: Time_Simulation=-1.
  integer, dimension(7) :: Start_Time_Array, Time_Array
  real*8                :: Real_Time_of_Year

  ! Optimization switches
  logical :: optimize_faceflux=.false., optimize_conservativeflux=.true.
  character*10 :: optimize_message_pass='dir'

  ! Solar Wind Input Parameters

  logical :: UseUpstreamInputFile

  integer, parameter :: Max_Upstream_Npts = 5000
  integer :: Upstream_Npts

  real, dimension(Max_Upstream_Npts, 8)   :: Upstream_Data
  real*8, dimension(Max_Upstream_Npts)    :: Upstream_Time

  ! Solar Wind Input Variables
  integer, parameter :: GSM = 1
  integer, parameter :: GSE = 2
  integer, parameter :: GEO = 3
  integer :: Input_Coor_System
  real :: Propagation_Plane_XY, Propagation_Plane_XZ
  real :: Satellite_Y_Pos, Satellite_Z_Pos

  ! Timing variables

  logical:: UseTiming
  integer:: dn_timing
  integer:: niter_timing

  real :: CodeSpeed

  real*8 :: time_total, time_start, time_end, time_tmp
  real*8 :: time_advance_expl = 0.0
  real*8 :: time_projection=0.0
  real*8 :: time_set_BCs = 0.0
  real*8 :: time_calc_facevalues_bfo = 0.0
  real*8 :: time_calc_facefluxes_bfo = 0.0
  real*8 :: time_send = 0.0
  real*8 :: time_calc_facevalues = 0.0
  real*8 :: time_calc_facefluxes = 0.0
  real*8 :: time_calc_sources = 0.0
  real*8 :: time_update = 0.0
  real*8 :: time_update_check = 0.0
  real*8 :: time_exchange = 0.0
  real*8 :: time_barrier = 0.0
  real*8 :: time_other = 0.0

end module ModMain
