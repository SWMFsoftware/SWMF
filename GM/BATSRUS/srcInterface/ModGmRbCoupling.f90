!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
module ModGmRbCoupling

  use ModMpi
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use CON_coupler, ONLY: Grid_C, RB_, ncells_decomposition_d

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace, ONLY : ray,rayface
  use ModPhysics, ONLY: No2Si_V, Si2No_V, UnitP_, UnitRho_, UnitTemperature_, &
       UnitB_, Bdp, rCurrents, rBody
  implicit none

  character (len=*), parameter :: NameMod='ModGmRbCoupling'

  ! RB Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the RB grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RB_lat, RB_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_SUM_vol, &
       MHD_SUM_rho, &
       MHD_SUM_p, &
       MHD_Beq, &
       MHD_Xeq, &
       MHD_Yeq
  real, parameter :: noValue=-99999.
  real :: qb(3),eqB,xL,yL,zL
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s6,s8,factor1,factor2

  integer, parameter :: maxMessages=10
  integer :: itag, NewMsg, iError
  integer :: nRECVrequests, RECVrequests(maxMessages), &
       nSENDrequests, SENDrequests(maxMessages), &
       MESGstatus(MPI_STATUS_SIZE, maxMessages)  

  logical :: dbg0=.false.

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

contains

  subroutine allocate_gm_rb(iSizeIn,jSizeIn)
    use CON_comp_param, ONLY: RB_

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_rb'

    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(RB_lat))then
       nCells_D=ncells_decomposition_d(RB_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RB_lat(iSize), RB_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       RB_lat = Grid_C(RB_) % Coord1_I
       RB_lon = Grid_C(RB_) % Coord2_I
    end if

    if(.not.allocated(MHD_SUM_vol))then
       allocate( MHD_SUM_vol(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_vol")
       MHD_SUM_vol = 0.
    end if

    if(.not.allocated(MHD_SUM_rho))then
       allocate( MHD_SUM_rho(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_rho")
       MHD_SUM_rho = 0.
    end if

    if(.not.allocated(MHD_SUM_p))then
       allocate( MHD_SUM_p(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_p")
       MHD_SUM_p = 0.
    end if

    if(.not.allocated(MHD_Beq))then
       allocate( MHD_Beq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Beq")
       MHD_Beq = 0.
    end if

    if(.not.allocated(MHD_Xeq))then
       allocate( MHD_Xeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Xeq")
       MHD_Xeq = 0.
    end if

    if(.not.allocated(MHD_Yeq))then
       allocate( MHD_Yeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Yeq")
       MHD_Yeq = 0.
    end if

    if(.not.allocated(MHD_lat_boundary))then
       allocate( MHD_lat_boundary(jsize), stat=iError )
       call alloc_check(iError,"MHD_lat_boundary")
       MHD_lat_boundary = 0
    end if

  end subroutine allocate_gm_rb

  !============================================================================
  subroutine write_integrated_data_tec
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: j2, nCall=0
    real :: tmpT, tmpV1,tmpV2, lonShift
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
         ', "Xeq", "Yeq"', &
         ', "Volume", "Volume**(-2/3)"', &
         ', "MHD `r", "MHD p", "MHD T", "Beq"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)*Si2No_V(UnitP_)) / &
               (MHD_SUM_rho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          write(UNITTMP_,'(2i4,12G14.6)') j2,i,RB_lon(j)+lonShift,RB_lat(i), &
               MHD_lat_boundary(j), &
               MHD_Xeq(i,j),MHD_Yeq(i,j), &
               tmpV1,tmpV2, &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_tec

  !============================================================================
  subroutine write_integrated_data_idl
    use ModPhysics, ONLY: No2Si_V, UnitN_, UnitU_, UnitB_, UnitP_,UnitRho_
    use ModIoUnit, ONLY: UNITTMP_
    use ModMain,   ONLY: time_simulation
    CHARACTER (LEN=100) :: filename
    integer :: nCall = 0
    !-------------------------------------------------------------------------

    !write values to plot file
    nCall = nCall+1
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,6
    write(UNITTMP_,'(3i4)')            jSize,iSize
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)')            'Lon Lat Xeq Yeq vol rho p Beq nothing'

    do i=isize,1,-1
       do j=1,jsize
          write(UNITTMP_,'(100(1pe18.10))') &
               modulo(RB_lon(j)-180.0,360.0),RB_lat(i), &
               MHD_Xeq(i,j),MHD_Yeq(i,j),&
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Beq(i,j) * No2Si_V(UnitB_)
       end do
    end do
    CLOSE(UNITTMP_)
    
  end subroutine write_integrated_data_idl


end module ModGmRbCoupling
