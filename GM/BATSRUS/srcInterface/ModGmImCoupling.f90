!^CFG COPYRIGHT UM
!^CMP FILE IM
module ModGmImCoupling

  use ModMpi
  use ModNumConst
  use CON_coupler

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace, ONLY : ray,rayface
  use ModPhysics, ONLY: unitSI_p, unitSI_rho, unitSI_temperature, unitSI_b, &
       Bdp_dim, rCurrents, rBody
  implicit none

  character (len=*), parameter :: NameMod='ModGmImCoupling'

  ! IM Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the IM grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RCM_lat, RCM_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_PE_vol,  MHD_SUM_vol, MHD_tmp, &
       MHD_PE_rho,  MHD_SUM_rho, &
       MHD_PE_p,    MHD_SUM_p, &
       MHD_PE_Beq,  MHD_Beq, &
       MHD_PE_Xeq,  MHD_Xeq, &
       MHD_PE_Yeq,  MHD_Yeq, &
       MHD_Fluxerror
  real, parameter :: noValue=-99999.
  real :: qb(3),eqB,xL,yL,zL
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s8,factor1,factor2

  integer, parameter :: maxMessages=10
  integer :: itag, NewMsg, iError
  integer :: nRECVrequests, RECVrequests(maxMessages), &
       nSENDrequests, SENDrequests(maxMessages), &
       MESGstatus(MPI_STATUS_SIZE, maxMessages)  

  logical :: dbg0=.false.

  integer, save :: nCalls=0

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

contains

  subroutine allocate_gm_im(iSizeIn,jSizeIn)

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_im'

    if(allocated(MHD_lat_boundary)) deallocate(MHD_lat_boundary)
    if(allocated(MHD_PE_vol))       deallocate(MHD_PE_vol)
    if(allocated(MHD_SUM_vol))      deallocate(MHD_SUM_vol)
    if(allocated(MHD_PE_rho))       deallocate(MHD_PE_rho)
    if(allocated(MHD_SUM_rho))      deallocate(MHD_SUM_rho)
    if(allocated(MHD_PE_p))         deallocate(MHD_PE_p)
    if(allocated(MHD_SUM_p))        deallocate(MHD_SUM_p)
    if(allocated(MHD_tmp))          deallocate(MHD_tmp)
    if(allocated(MHD_PE_Beq))       deallocate(MHD_PE_Beq)
    if(allocated(MHD_PE_Xeq))       deallocate(MHD_PE_Xeq)
    if(allocated(MHD_PE_Yeq))       deallocate(MHD_PE_Yeq)
    if(allocated(MHD_Beq))          deallocate(MHD_Beq)
    if(allocated(MHD_Xeq))          deallocate(MHD_Xeq)
    if(allocated(MHD_Yeq))          deallocate(MHD_Yeq)
    if(allocated(MHD_Fluxerror))    deallocate(MHD_Fluxerror)

    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(RCM_lat))then
       nCells_D=ncells_decomposition_d(IM_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RCM_lat(iSize), RCM_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       RCM_lat = (cHalfPi - Grid_C(IM_) % Coord1_I) * cRadToDeg
       RCM_lon =  Grid_C(IM_) % Coord2_I             * cRadToDeg
    end if

    ! Arrays needed for the field line integrals
    allocate( MHD_SUM_vol(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_vol")
    MHD_SUM_vol = 0.

    allocate( MHD_SUM_rho(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_rho")
    MHD_SUM_rho = 0.

    allocate( MHD_SUM_p(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_p")
    MHD_SUM_p = 0.

    allocate( MHD_Beq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Beq")
    MHD_Beq = 0.

    allocate( MHD_Xeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Xeq")
    MHD_Xeq = 0.

    allocate( MHD_Yeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Yeq")
    MHD_Yeq = 0.

    allocate( MHD_PE_vol(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_vol")
    MHD_PE_vol = 0.

    allocate( MHD_PE_rho(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_rho")
    MHD_PE_rho = 0.

    allocate( MHD_PE_p(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_p")
    MHD_PE_p = 0.

    allocate( MHD_PE_Beq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Beq")
    MHD_PE_Beq = noValue

    allocate( MHD_PE_Xeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Xeq")
    MHD_PE_Xeq = noValue

    allocate( MHD_PE_Yeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Yeq")
    MHD_PE_Yeq = noValue

    allocate( MHD_tmp(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_tmp")
    MHD_tmp = 0.

    allocate( MHD_Fluxerror(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Fluxerror")
    MHD_Fluxerror = 0.

    allocate( MHD_lat_boundary(jsize), stat=iError )
    call alloc_check(iError,"MHD_lat_boundary")
    MHD_lat_boundary = 0

  end subroutine allocate_gm_im

  !============================================================================
  subroutine write_integrated_data
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: j2
    real :: tmpT, tmpV1,tmpV2, lonShift
    !-------------------------------------------------------------------------

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCalls,".dat"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
         ', "Xeq", "Yeq"', &
         ', "Volume", "Volume**(-2/3)"', &
         ', "MHD `r", "MHD p", "MHD T", "Beq"', &
         ', "FluxError"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)/unitSI_p)/(MHD_SUM_rho(i,j)/unitSI_rho)) &
               * unitSI_temperature
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          write(UNITTMP_,'(2i4,12G14.6)') j2,i,RCM_lon(j)+lonShift,RCM_lat(i), &
               MHD_lat_boundary(j), &
               MHD_Xeq(i,j),MHD_Yeq(i,j), &
               tmpV1,tmpV2, &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j), &
               MHD_Fluxerror(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data

  !============================================================================
  subroutine write_integrated_data_idl

    use ModIoUnit, ONLY: UNITTMP_
    use ModMain,   ONLY: time_simulation
    integer :: nStepLast = -1, nSubStep = 1
    CHARACTER (LEN=100) :: filename
    !-------------------------------------------------------------------------

    !write values to plot file

    if(nStepLast/=n_step)then
       nSubStep = 1
       write(filename,'(a,i6.6,a)')"rayValues_n=",n_step,".out"
    else
       nSubStep = nSubStep + 1
       write(filename,'(a,i6.6,a,i1,a)')&
            "rayValues_n=",n_step,"_",nSubStep,".out"
    end if
    nStepLast = n_step

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
               RCM_lon(j),RCM_lat(i), &
               MHD_Xeq(i,j),MHD_Yeq(i,j),&
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Beq(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_idl

  !============================================================================
  subroutine MHD_compute_from_raytrace
    integer :: iray
    real :: rLat,rLon, minLat,maxLat,minLon,maxLon
    real, dimension(3,2,1:nI+1,1:nJ+1,1:nK+1):: rayface_tmp
    real, dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d !assumes cube!!!!!

    logical :: dbg,found, shifted
    real :: fvol,rvol,pvol
    !-----------------------------------------------------------------------
    call timing_start('compute_ray')
 
    dbg=.false.

    MHD_PE_vol = 0.
    MHD_PE_rho = 0.
    MHD_PE_p   = 0.

    do iBLK=1,nBlockMax
       if(unusedBLK(iBLK)) CYCLE

       !copy to rayface_tmp
       rayface_tmp=rayface(:,:,:,:,:,iBLK)

       !convert to LatLonStatus variables
       call convFaces2LatLon(rayface_tmp)

       !put x,y,z into rayface_tmp, overwrite southern hemisphere values
       rayface_tmp(1,2,:,:,:)=x_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dx_BLK(iBLK)
       rayface_tmp(2,2,:,:,:)=y_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dy_BLK(iBLK)
       rayface_tmp(3,2,:,:,:)=z_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dz_BLK(iBLK)

       !\
       ! Find closed fieldlines on any plane of a block and integrate for volume
       !/
       if(any(rayface_tmp(3,1,1:nI+1,1:nJ+1,1:nK+1)==3.))then

          !zmin face
          zL=z_BLK(1,1,1,iBLK)-0.5*dz_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,:,:,1)
          rayface_2d(1,2,:,:)=rayface_tmp(1,2,:,:,1)
          rayface_2d(2,2,:,:)=rayface_tmp(2,2,:,:,1)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
                   if(abs(xL)>.01 .or. abs(yL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

          !zmax face
          zL=z_BLK(1,1,nK+1,iBLK)-0.5*dz_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,:,:,nK+1)
          rayface_2d(1,2,:,:)=rayface_tmp(1,2,:,:,nK+1)
          rayface_2d(2,2,:,:)=rayface_tmp(2,2,:,:,nK+1)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
                   if(abs(xL)>.01 .or. abs(yL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

          !ymin face
          yL=y_BLK(1,1,1,iBLK)-0.5*dy_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,:,1,:)
          rayface_2d(1,2,:,:)=rayface_tmp(1,2,:,1,:)
          rayface_2d(2,2,:,:)=rayface_tmp(3,2,:,1,:)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,xL,zL,rayface_2d)
                   if(abs(xL)>.01 .or. abs(zL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

          !ymax face
          yL=y_BLK(1,nJ+1,1,iBLK)-0.5*dy_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,:,nJ+1,:)
          rayface_2d(1,2,:,:)=rayface_tmp(1,2,:,nJ+1,:)
          rayface_2d(2,2,:,:)=rayface_tmp(3,2,:,nJ+1,:)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,xL,zL,rayface_2d)
                   if(abs(xL)>.01 .or. abs(zL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

          !xmin face
          xL=x_BLK(1,1,1,iBLK)-0.5*dx_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,1,:,:)
          rayface_2d(1,2,:,:)=rayface_tmp(2,2,1,:,:)
          rayface_2d(2,2,:,:)=rayface_tmp(3,2,1,:,:)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,yL,zL,rayface_2d)
                   if(abs(yL)>.01 .or. abs(zL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

          !xmax face
          xL=x_BLK(nI+1,1,1,iBLK)-0.5*dx_BLK(iBLK)
          rayface_2d=0.
          rayface_2d(:,1,:,:)=rayface_tmp(:,1,nI+1,:,:)
          rayface_2d(1,2,:,:)=rayface_tmp(2,2,nI+1,:,:)
          rayface_2d(2,2,:,:)=rayface_tmp(3,2,nI+1,:,:)
          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then
             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)
             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE
                   call find_rayface_location(dbg,iBLK,rLat,rLon,yL,zL,rayface_2d)
                   if(abs(yL)>.01 .or. abs(zL)>.01)then
                      !Trace value at xL,yL,zL
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
                      MHD_PE_vol(i1,i2) = MHD_PE_vol(i1,i2)+fvol
                      MHD_PE_rho(i1,i2) = MHD_PE_rho(i1,i2)+rvol
                      MHD_PE_p(i1,i2)   = MHD_PE_p(i1,i2)+pvol
                   end if
                end do
             end do
          end if

       end if

       !\
       ! Find closed fieldlines on Z=0 plane
       !/
       if( abs(z_BLK(1,1,1,iBLK)-0.5*dz_BLK(iBLK)) < 0.01 )then

          found=.false.
          !copy to rayface_tmp
          rayface_tmp=rayface(:,:,:,:,:,iBLK)

          !convert to LatLonStatus variables
          call convFaces2LatLon(rayface_tmp)

          rayface_2d(:,1,:,:)=rayface_tmp(:,1,:,:,1)

          where(rayface_2d(1,1,:,:)<0.) &
               rayface_2d(1,1,:,:)=-rayface_2d(1,1,:,:)

          !put x,y,z into rayface_tmp, overwrite southern hemisphere values
          rayface_tmp(1,2,:,:,:)=x_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dx_BLK(iBLK)
          rayface_tmp(2,2,:,:,:)=y_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dy_BLK(iBLK)
          rayface_tmp(3,2,:,:,:)=z_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dz_BLK(iBLK)

          rayface_2d(:,2,:,:)=rayface_tmp(:,2,:,:,1)

          if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then

             call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)

             do i1=1,isize
                rLat=RCM_lat(i1)
                if(rLat<minLat .or. rLat>maxLat) CYCLE
                do i2=1,jsize
                   rLon=RCM_lon(i2)
                   if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                   if(rLon<minLon .or. rLon>maxLon) CYCLE

                   call find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
                   if(abs(xL)>.01 .or. abs(yL)>.01)then

                      MHD_PE_Xeq(i1,i2)=xL
                      MHD_PE_Yeq(i1,i2)=yL
                      zL=1.E-6
                      call find_eqB(dbg,iBLK,xL,yL,zL,eqB)
                      MHD_PE_Beq(i1,i2)=eqB

                   end if
                end do
             end do

          end if

       end if

    end do


    !collect the computed values to iProc0
    call MPI_REDUCE(MHD_PE_vol(1,1),      MHD_SUM_vol(1,1),      &
         isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
    call MPI_REDUCE(MHD_PE_rho(1,1),      MHD_SUM_rho(1,1),      &
         isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
    call MPI_REDUCE(MHD_PE_p(1,1),        MHD_SUM_p(1,1),        &
         isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
    call MPI_REDUCE(MHD_PE_Beq(1,1),      MHD_Beq(1,1),          &
         isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)
    call MPI_REDUCE(MHD_PE_Xeq(1,1),      MHD_Xeq(1,1),          &
         isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)
    call MPI_REDUCE(MHD_PE_Yeq(1,1),      MHD_Yeq(1,1),          &
         isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)


    call timing_stop('compute_ray')

  end subroutine MHD_compute_from_raytrace

  !===========================================================================
  subroutine process_integrated_data

    if(dbg0)then
       if(iProc==0)write(*,*)' -D-'
    end if

    where(MHD_SUM_vol>0.) &
         MHD_SUM_p=MHD_SUM_p/MHD_SUM_vol

    where(MHD_SUM_vol>0.) &
         MHD_SUM_rho=MHD_SUM_rho/MHD_SUM_vol

    !Set volume floor
    MHD_SUM_vol = max(1.E-8,MHD_SUM_vol)

    if(dbg0)then
       if(iProc==0)write(*,*)' -E-'
    end if

    ! If the field-line tracer returned a good value, we may need to
    ! add contribution to V from the part of the field line that
    ! goes inside the body. If the field-line tracer did not
    ! return a good value, we will compute total V assuming dipole.
    Ri=(6378.+100.)/6378.
    factor=(unitSI_B*1.E+9) * 2. * (Ri**4) / 31100.
    do i=1,isize
       colat = cPi/2. - RCM_lat(i)*(cPi/180.)
       s2=(sin(colat))**2
       s8=(sin(colat))**8
       Ci=cos(colat)

!!$       if( s2 < Ri/Rbody .and. any(MHD_SUM_vol(i,:)>0.) )then
!!$          !Fieldline goes beyond Rbody, add piece of fieldline volume
!!$          Cs=sqrt(1.-(Rbody/Ri)*s2)
!!$          FCiCs=(Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) - (1./7.)*(Ci**7-Cs**7)
!!$          Vol=factor*FCiCs/s8
!!$
!!$          where(MHD_SUM_vol(i,:)>1.1E-8)
!!$             MHD_SUM_vol(i,:)=MHD_SUM_vol(i,:) + Vol
!!$          end where
!!$       end if

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values

          !Compute the full volume
          Cs=0.
          FCiCs=(Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) - (1./7.)*(Ci**7-Cs**7)
          Vol=factor*FCiCs/s8

          if( s2 > Ri/Rbody )then
             !Fieldline stays inside of Rbody, recompute exact volume
             MHD_SUM_vol(i,:)=Vol

             !Fix the grid inside Rbody
             MHD_Xeq(i,:) = (Ri/s2)*cos(RCM_lon(:)*(cPi/180.))
             MHD_Yeq(i,:) = (Ri/s2)*sin(RCM_lon(:)*(cPi/180.))

             !Fix the equatorial B value
             xL=MHD_Xeq(i,1)
             yL=MHD_Yeq(i,1)
             zL=0.
             call get_b0(xL,yL,zL,qb)
             eqB=sqrt(qb(1)**2+qb(2)**2+qb(3)**2)
             MHD_Beq(i,:)=eqB

          else
             !Fieldline stays inside of Rcurrents but outside Rbody
             factor1=( 1. / ((Ri/Rbody)-(Ri/Rcurrents)) )*((Ri/Rbody)-s2)
             factor2=( 1. / ((Ri/Rbody)-(Ri/Rcurrents)) )*(s2-(Ri/Rcurrents))

             !Blend volume with exact volume
             where(MHD_SUM_vol(i,:)>1.1E-8)
                MHD_SUM_vol(i,:)= factor1*MHD_SUM_vol(i,:) + factor2*Vol
             end where
             where(MHD_SUM_vol(i,:)<1.1E-8)
                MHD_SUM_vol(i,:)=Vol
             end where

          end if
       end if
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -F-'
    end if

    !set open fieldline values
    do j=1,jsize
       i0=isize
       do i=isize,1,-1
          if(MHD_SUM_vol(i,j)<1.1E-8 .OR. &
               MHD_Xeq(i,j)<=noValue .OR. &
               MHD_Yeq(i,j)<=noValue ) then
             i0=i+1
             exit
          end if
       end do

       i=i0
       MHD_lat_boundary(j)=i

       MHD_SUM_vol(1:i-1,j)=MHD_SUM_vol(i,j)
       MHD_Beq    (1:i-1,j)=MHD_Beq    (i,j)
       MHD_Xeq    (1:i-1,j)=MHD_Xeq    (i,j)
       MHD_Yeq    (1:i-1,j)=MHD_Yeq    (i,j)
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -G-'
    end if

    !error checking
    do j=1,jsize; do i=isize,1,-1
       IF (iProc ==0) then
          if (mhd_sum_vol(i,j) > 1.0E-8 .and. &
               (MHD_Xeq(i,j) < -999. .or. MHD_Yeq(i,j) < -999.)) then
             print*,'999999:',i,j,MHD_Xeq(i,j),MHD_Yeq(i,j), mhd_sum_vol(i,j)
             call stop_mpi('stanislav')
          end if
          if (mhd_lat_boundary(j) < i) then
             if (mhd_sum_vol(i,j) < 1.1E-8) then
                print*,'7777',i,j,mhd_lat_boundary(j),mhd_sum_vol(i,j), &
                     MHD_Xeq(i,j),MHD_Yeq(i,j)
                CALL STOP_MPI ('HOLE ')
             end if
          end if
       end if
    end do; enddo

    if(dbg0)then
       if(iProc==0)write(*,*)' -H-'
    end if

    !set values for open fieldlines
    do j=1,jsize
       MHD_SUM_vol(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_rho(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_p  (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Beq    (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Xeq    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Xeq(int(MHD_lat_boundary(j)),j)
       MHD_Yeq    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Yeq(int(MHD_lat_boundary(j)),j)
    end do

    !dimensionalize values
    where(MHD_SUM_vol>0.)
       MHD_SUM_vol = MHD_SUM_vol  / unitSI_B
       MHD_SUM_rho = MHD_SUM_rho * unitSI_rho
       MHD_SUM_p   = MHD_SUM_p * unitSI_p
    elsewhere
       MHD_SUM_vol = -1.
       MHD_SUM_rho = -1.
       MHD_SUM_p   = -1.
    end where

    where(MHD_Beq>0.)
       MHD_Beq = MHD_Beq * unitSI_B
    elsewhere
       MHD_Beq = -1.
    end where

    if(dbg0)then
       if(iProc==0)write(*,*)' -I-'
    end if

    do j=1,jsize-1; do i=1,isize-1
       if ( MHD_SUM_vol(i  ,j) > 0.0 .AND. &
            MHD_SUM_vol(i+1,j) > 0.0 .AND. &
            MHD_SUM_vol(i,j+1) > 0.0 .AND. &
            MHD_SUM_vol(i+1,j+1)> 0.0) THEN
          MHD_Fluxerror(i,j) = &
               0.25E+9*(MHD_Beq(i,j  ) + MHD_Beq(i+1,j  ) + &
               MHD_Beq(i,j+1) + MHD_Beq(i+1,j+1)) * &
               0.5*(ABS((MHD_Xeq(i,j+1)-MHD_Xeq(i,j))*(MHD_Yeq(i+1,j)-MHD_Yeq(i,j)) - &
               (MHD_Yeq(i,j+1)-MHD_Yeq(i,j))*(MHD_Xeq(i+1,j)-MHD_Xeq(i,j))) +&
               ABS((MHD_Xeq(i+1,j)-MHD_Xeq(i+1,j+1))*(MHD_Yeq(i,j+1)-MHD_Yeq(i+1,j+1)) -&
               (MHD_Yeq(i+1,j)-MHD_Yeq(i+1,j+1))*(MHD_Xeq(i,j+1)-MHD_Xeq(i+1,j+1))))/&
               (ABS(Bdp_dim)*(SIN(Rcm_lat(i)/57.296)**2-SIN(Rcm_lat(i+1)/57.296)**2)* &
               (RCM_lon(j+1)-RCM_lon(j))/57.296 )- 1.0
       ELSE
          MHD_Fluxerror (i,j) = 0.0
       END IF
    end do; end do
    MHD_fluxerror(:,jsize) = MHD_Fluxerror (:,1)
    MHD_fluxerror(isize,:) = MHD_FLuxerror(isize-1,:)

    if(dbg0)then
       if(iProc==0)write(*,*)' -J-'
    end if

  end subroutine process_integrated_data

end module ModGmImCoupling
