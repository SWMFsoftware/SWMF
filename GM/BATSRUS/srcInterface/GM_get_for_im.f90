!^CFG COPYRIGHT UM
!^CMP FILE IM
!==========================================================================
!==========================================================================
!==========================================================================
subroutine GM_get_for_im(Buffer_IIV,iSize,jSize,nVar,NameVar)

  !call stop_mpi('RCM is OFF') !^CFG UNCOMMENT IF NOT RCM
  !^CFG IF RCM BEGIN
  use ModMpi
  use ModNumConst
  use CON_coupler

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics, nVarGm => nVar
  use ModRaytrace, ONLY : ray,rayface
  implicit none

  integer::nCells_D(2)
  character (len=*), parameter :: NameSub='GM_get_for_im'

  integer, intent(in) :: iSize,jSize,nVar
  real, intent(out), dimension(iSize,jSize,nVar) :: Buffer_IIV
  character (len=*), intent(in) :: NameVar

  ! Information about the IM grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RCM_lat, RCM_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK


  integer, save, dimension(:,:), allocatable :: &
       MHD_PE_Cnumber,  MHD_SUM_Cnumber, &
       MHD_PE_Fnumber,  MHD_SUM_Fnumber, &
       MHD_PE_Z0number, MHD_SUM_Z0number
  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_PE_vol,  MHD_SUM_vol, MHD_SUM_vol0, MHD_tmp, &
       MHD_PE_rho,  MHD_SUM_rho, &
       MHD_PE_p,    MHD_SUM_p, &
       MHD_PE_Z0T,  MHD_Z0T, &
       MHD_PE_Bmin, MHD_Bmin, &
       MHD_PE_Z0x,  MHD_Z0x, &
       MHD_PE_Z0y,  MHD_Z0y

  real, parameter :: noValue=-99999.
  real :: tmpValue(10),tmpCount1,tmpCount2,tmpCount3
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s8,factor1,factor2

  integer, parameter :: maxMessages=10
  integer :: itag, NewMsg, iError
  integer :: nRECVrequests, RECVrequests(maxMessages), &
       nSENDrequests, SENDrequests(maxMessages), &
       MESGstatus(MPI_STATUS_SIZE, maxMessages)  

  real*8 :: time_this

  logical :: dbg0=.false.

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6
  !--------------------------------------------------------------------------
  if(NameVar /= 'vol:z0x:z0y:bmin:rho:p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  time_this=MPI_WTIME()

  ! Make sure that we have ray tracing
  call ray_trace

  !\
  ! First, deallocate all arrays
  !/
  if(allocated(MHD_PE_Cnumber))   deallocate(MHD_PE_Cnumber)
  if(allocated(MHD_PE_Fnumber))   deallocate(MHD_PE_Fnumber)
  if(allocated(MHD_PE_Z0number))  deallocate(MHD_PE_Z0number)
  if(allocated(MHD_SUM_Cnumber))  deallocate(MHD_SUM_Cnumber)
  if(allocated(MHD_SUM_Fnumber))  deallocate(MHD_SUM_Fnumber)
  if(allocated(MHD_SUM_Z0number)) deallocate(MHD_SUM_Z0number)
  if(allocated(MHD_lat_boundary)) deallocate(MHD_lat_boundary)
  if(allocated(MHD_PE_vol))       deallocate(MHD_PE_vol)
  if(allocated(MHD_SUM_vol))      deallocate(MHD_SUM_vol)
  if(allocated(MHD_PE_rho))       deallocate(MHD_PE_rho)
  if(allocated(MHD_SUM_rho))      deallocate(MHD_SUM_rho)
  if(allocated(MHD_PE_p))         deallocate(MHD_PE_p)
  if(allocated(MHD_SUM_p))        deallocate(MHD_SUM_p)
  if(allocated(MHD_tmp))          deallocate(MHD_tmp)
  if(allocated(MHD_SUM_vol0))     deallocate(MHD_SUM_vol0)
  if(allocated(MHD_PE_Z0T))       deallocate(MHD_PE_Z0T)
  if(allocated(MHD_PE_Bmin))      deallocate(MHD_PE_Bmin)
  if(allocated(MHD_PE_Z0x))       deallocate(MHD_PE_Z0x)
  if(allocated(MHD_PE_Z0y))       deallocate(MHD_PE_Z0y)
  if(allocated(MHD_Z0T))          deallocate(MHD_Z0T)
  if(allocated(MHD_Bmin))         deallocate(MHD_Bmin)
  if(allocated(MHD_Z0x))          deallocate(MHD_Z0x)
  if(allocated(MHD_Z0y))          deallocate(MHD_Z0y)

  if(.not.allocated(RCM_lat))then
     nCells_D=ncells_decomposition_d(IM_)
     if(  iSize /= nCells_D(1) .or. &
          jSize /= nCells_D(2) ) then
        write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
             iSize,jSize, nCells_D(1:2)
        call CON_stop(NameSub)
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

  allocate( MHD_Bmin(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_Bmin")
  MHD_Bmin = 0.

  allocate( MHD_Z0x(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_Z0x")
  MHD_Z0x = 0.

  allocate( MHD_Z0y(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_Z0y")
  MHD_Z0y = 0.

  allocate( MHD_PE_Cnumber(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Cnumber")
  MHD_PE_Cnumber = 0

  allocate( MHD_PE_Fnumber(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Fnumber")
  MHD_PE_Fnumber = 0

  allocate( MHD_PE_Z0number(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Z0number")
  MHD_PE_Z0number = 0

  allocate( MHD_PE_vol(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_vol")
  MHD_PE_vol = 0.

  allocate( MHD_PE_rho(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_rho")
  MHD_PE_rho = 0.

  allocate( MHD_PE_p(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_p")
  MHD_PE_p = 0.

  allocate( MHD_PE_Z0T(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Z0T")
  MHD_PE_Z0T = 0.

  allocate( MHD_PE_Bmin(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Bmin")
  MHD_PE_Bmin = noValue

  allocate( MHD_PE_Z0x(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Z0x")
  MHD_PE_Z0x = noValue

  allocate( MHD_PE_Z0y(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_PE_Z0y")
  MHD_PE_Z0y = noValue

  allocate( MHD_SUM_Cnumber(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_SUM_Cnumber")
  MHD_SUM_Cnumber = 0

  allocate( MHD_SUM_Fnumber(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_SUM_Fnumber")
  MHD_SUM_Fnumber = 0

  allocate( MHD_SUM_Z0number(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_SUM_Z0number")
  MHD_SUM_Z0number = 0

  allocate( MHD_tmp(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_tmp")
  MHD_tmp = 0.

  allocate( MHD_SUM_vol0(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_SUM_vol0")
  MHD_SUM_vol0 = 0.

  allocate( MHD_Z0T(isize,jsize), stat=iError )
  call alloc_check(iError,"MHD_Z0T")
  MHD_Z0T = 0.

  allocate( MHD_lat_boundary(jsize), stat=iError )
  call alloc_check(iError,"MHD_lat_boundary")
  MHD_lat_boundary = 0

  if(dbg0)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
     if(iProc==0)write(*,*)' -A-'
  end if

  !compute various values
  call MHD_compute_from_raytrace

  if(dbg0)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
     if(iProc==0)write(*,*)' -B-'
  end if

  !collect the computed values to iProc0
  call MPI_REDUCE(MHD_PE_Cnumber(1,1),  MHD_SUM_Cnumber(1,1),  &
       isize*jsize, MPI_INTEGER, MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Fnumber(1,1),  MHD_SUM_Fnumber(1,1),  &
       isize*jsize, MPI_INTEGER, MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Z0number(1,1), MHD_SUM_Z0number(1,1), &
       isize*jsize, MPI_INTEGER, MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_vol(1,1),      MHD_SUM_vol(1,1),      &
       isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_rho(1,1),      MHD_SUM_rho(1,1),      &
       isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_p(1,1),        MHD_SUM_p(1,1),        &
       isize*jsize, MPI_REAL,    MPI_SUM, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Z0T(1,1),      MHD_Z0T(1,1),          &
       isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Bmin(1,1),     MHD_Bmin(1,1),         &
       isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Z0x(1,1),      MHD_Z0x(1,1),          &
       isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)
  call MPI_REDUCE(MHD_PE_Z0y(1,1),      MHD_Z0y(1,1),          &
       isize*jsize, MPI_REAL,    MPI_MAX, 0, iComm, iError)

  if(dbg0)then
     call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
     if(iProc==0)write(*,*)' -C-'
  end if

  if (iProc == 0) then
     call process_integrated_data
     if(dbg0) call write_integrated_data

     ! Put results into output buffer
     Buffer_IIV(:,:,vol_) = MHD_SUM_vol
     Buffer_IIV(:,:,z0x_) = MHD_Z0x
     Buffer_IIV(:,:,z0y_) = MHD_Z0y
     Buffer_IIV(:,:,bmin_)= MHD_Bmin
     Buffer_IIV(:,:,rho_) = MHD_SUM_rho
     Buffer_IIV(:,:,p_)   = MHD_SUM_p

  end if

contains

  !============================================================================
  subroutine write_integrated_data

    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    !-------------------------------------------------------------------------

    !write values to plot file
    write(filename,'(a,i6.6,a)')"rayValues_n=",n_step,".dat"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary"', &
         ', "Z=0 X", "Z=0 Y", "Z=0 Count"', &
         ', "Cell Count", "Face Count"', &
         ', "Orig Volume", "Volume"', &
         ', "MHD `r", "MHD p", "Equatorial T", "Equatorial B"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j=1,jsize
          write(UNITTMP_,'(2i4,5G14.6,3i7,6G14.6)') &
               j,i,RCM_lon(j),RCM_lat(i), &
               MHD_lat_boundary(j), &
               MHD_Z0x(i,j),MHD_Z0y(i,j),MHD_SUM_Z0number(i,j),&
               MHD_SUM_Cnumber(i,j),MHD_SUM_Fnumber(i,j), &
               MHD_SUM_vol0(i,j),MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Z0T(i,j),MHD_Bmin(i,j)
       end do
       j=1
       write(UNITTMP_,'(2i4,5G14.6,3i7,6G14.6)') &
            j+jsize,i,RCM_lon(j)+360.,RCM_lat(i), &
            MHD_lat_boundary(j), &
            MHD_Z0x(i,j),MHD_Z0y(i,j),MHD_SUM_Z0number(i,j), &
            MHD_SUM_Cnumber(i,j),MHD_SUM_Fnumber(i,j), &
            MHD_SUM_vol0(i,j),MHD_SUM_vol(i,j), &
            MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Z0T(i,j),MHD_Bmin(i,j)
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data

  subroutine MHD_compute_from_raytrace
    integer :: iray
    real :: rLat,rLon, minLat,maxLat,minLon,maxLon, xL,yL,zL, eqT,eqB
    real, dimension(3,2,1:nI+1,1:nJ+1,1:nK+1):: rayface_tmp
    real, dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d !assumes cube!!!!!

    logical :: dbg,found, shifted
    real :: fvol,rvol,pvol

    dbg=.false.

    if(dbg0)then
       write(*,*)'MHD_compute_from_raytrace starting, iProc,isize, jsize=',&
            iProc,isize,jsize
       !write(*,*)'MHD_compute_from_raytrace iProc,RCM_lat=',iProc,RCM_lat
       !write(*,*)'MHD_compute_from_raytrace iProc,RCM_lon=',iProc,RCM_lon
    end if

    MHD_PE_vol = 0.
    MHD_PE_rho = 0.
    MHD_PE_p   = 0.

    do iBLK=1,nBlockMax
       if(unusedBLK(iBLK)) CYCLE

       !\
       ! Check to see if cell centers are on closed fieldline
       !/
       if(any(ray(3,1,1:nI,1:nJ,1:nK,iBLK)==3.))then
          do i=1,nI; do j=1,nJ; do k=1,nK
             if(ray(3,1,i,j,k,iBLK)==3.) then
                rLat=ray(1,1,i,j,k,iBLK)
                rLon=ray(2,1,i,j,k,iBLK)

                !NOTE: RCM_lat in decending order
                i1=1
                do n=2,isize
                   if(rLat<0.5*(RCM_lat(n-1)+RCM_lat(n)))then
                      i1=n
                   end if
                end do

                !NOTE: RCM_lon in accending order
                i2=1
                do n=2,jsize
                   if(rLon>0.5*(RCM_lon(n-1)+RCM_lon(n)))then
                      i2=n
                   end if
                end do
                if(rLon>0.5*(RCM_lon(jsize)+RCM_lon(1)+360.)) i2=1

                MHD_PE_Cnumber(i1,i2) = MHD_PE_Cnumber(i1,i2)+1
             end if
          end do; end do; end do
       end if

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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. (i2==19 .and. iBLK==14)) dbg=.true.
!!$ if(i1==45 .and. (i2==31 .and. iBLK==32)) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing zmin face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. i2==19) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing zmax face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. i2==19) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing ymin face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. i2==19) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing ymax face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. i2==19) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing xmin face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
                      MHD_PE_Fnumber(i1,i2) = MHD_PE_Fnumber(i1,i2)+1
!!$ if(i1==45 .and. i2==19) dbg=.true.
                      call integrate_ray(dbg,iBLK,xL,yL,zL,fvol,rvol,pvol)
!!$ dbg=.false.
!!$ if(i1==45 .and. (i2==19 .or. i2==31))then
!!$    write(*,*)iProc,iBLK,'  i=',i1,'j=',i2,' Tracing xmax face at ',xL,yL,zL,'  returning ',fvol,rvol,pvol
!!$ end if
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
          do iray=1,2
             !copy to rayface_tmp
             rayface_tmp=rayface(:,:,:,:,:,iBLK)

             !convert to LatLonStatus variables
             call convFaces2LatLon(rayface_tmp)

             rayface_2d(:,1,:,:)=rayface_tmp(:,iray,:,:,1)

             where(rayface_2d(1,1,:,:)<0.)
                rayface_2d(1,1,:,:)=-rayface_2d(1,1,:,:)
             end where

             !put x,y,z into rayface_tmp, overwrite southern hemisphere values
             rayface_tmp(1,2,:,:,:)=x_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dx_BLK(iBLK)
             rayface_tmp(2,2,:,:,:)=y_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dy_BLK(iBLK)
             rayface_tmp(3,2,:,:,:)=z_BLK(1:nI+1,1:nJ+1,1:nK+1,iBLK)-0.5*dz_BLK(iBLK)

             rayface_2d(:,2,:,:)=rayface_tmp(:,2,:,:,1)

             if(any(rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))then

!!$             !Write Z=0 plane if it contains any closed fieldlines
!!$             write(filename,'(a,i2.2,a,i4.4,a)')"Z=0_pe",iProc,"_blk",iBLK,".dat"
!!$             OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
!!$             write(UNITTMP_,'(a)') 'TITLE="Debug Test"'
!!$             write(UNITTMP_,'(a)') 'VARIABLES="X", "Y", "Lat", "Lon", "Status"'
!!$             write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
!!$                  'ZONE T="PE=',iProc,'", I=',nI+1,', J=',nJ+1,', K=1, F=POINT'
!!$             do j=1,nJ+1; do i=1,nI+1
!!$                write(UNITTMP_,'(5G14.6)') &
!!$                     rayface_2d(1,2,i,j),rayface_2d(2,2,i,j), &
!!$                     rayface_2d(1,1,i,j),rayface_2d(2,1,i,j),rayface_2d(3,1,i,j)
!!$             end do; end do
!!$             CLOSE(UNITTMP_)

                call getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)

!!$                if(iProc==1 .and. iBLK==1063)then
!!$                   i1=41; i2=47
!!$                   write(*,*)' '
!!$                   write(*,'(a)') 'VARIABLES="X", "Y", "Lat", "Lon", "Status"'
!!$                   do j=1,nJ+1; do i=1,nI+1
!!$                      write(*,'(5G14.6)') &
!!$                           rayface_2d(1,2,i,j),rayface_2d(2,2,i,j), &
!!$                           rayface_2d(1,1,i,j),rayface_2d(2,1,i,j),rayface_2d(3,1,i,j)
!!$                   end do; end do
!!$                   write(*,*)' '
!!$                   write(*,*)'Trying PE=',iProc,' BLK=',iBLK,' iray=',iray, &
!!$                        ' rLat=',RCM_lat(i1),' rLon=',RCM_lon(i2),'    ',i1,i2
!!$                   write(*,*)' '
!!$                   write(*,*)' minLat=',minLat,' maxLat=',maxLat
!!$                   write(*,*)' minLon=',minLon,' maxLon=',maxLon
!!$                   write(*,*)' '
!!$                end if

                do i1=1,isize
                   rLat=RCM_lat(i1)
                   if(rLat<minLat .or. rLat>maxLat) CYCLE
                   do i2=1,jsize
                      rLon=RCM_lon(i2)
                      if(shifted .and. RCM_lon(i2)<180.) rLon=RCM_lon(i2)+360.
                      if(rLon<minLon .or. rLon>maxLon) CYCLE

!!$                      if(iProc==1 .and. iBLK==1063 .and. i1==41 .and. i2 == 47)then
!!$                         dbg=.true.
!!$                      else
!!$                         dbg=.false.
!!$                      end if

                      call find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
                      if(abs(xL)>.01 .or. abs(yL)>.01)then

                         MHD_PE_Z0x(i1,i2)=xL
                         MHD_PE_Z0y(i1,i2)=yL
                         MHD_PE_Z0number(i1,i2) = MHD_PE_Z0number(i1,i2)+1
                         zL=1.E-4
                         call find_eqT(iBLK,xL,yL,zL,eqT,eqB)
                         MHD_PE_Z0T(i1,i2)=eqT
                         MHD_PE_Bmin(i1,i2)=eqB

                      end if
                   end do
                end do

             end if
          end do

       end if

    end do

  end subroutine MHD_compute_from_raytrace

  !===========================================================================
  subroutine process_integrated_data

    where(MHD_SUM_vol>0.) &
         MHD_SUM_p=MHD_SUM_p/MHD_SUM_vol

    where(MHD_SUM_vol>0.) &
         MHD_SUM_rho=MHD_SUM_rho/MHD_SUM_vol

    !Set volume floor
    MHD_SUM_vol = max(1.E-8,MHD_SUM_vol)

    !Copy original values
    MHD_SUM_vol0 = MHD_SUM_vol

    !Overwrite computed fieldline volumes inside body and
    ! add portion of fieldline inside body to others
    Ri=(6378.+100.)/6378.
    factor=(unitSI_B*1.E+9) * 2. * (Ri**4) / 31100.
    do i=1,isize
       colat = cPi/2. - RCM_lat(i)*(cPi/180.)
       s2=(sin(colat))**2
       s8=(sin(colat))**8
       Ci=cos(colat)

       if( s2 < Ri/Rbody .and. any(MHD_SUM_vol(i,:)>0.) )then
          !Fieldline goes beyond Rbody, add piece of fieldline volume
          Cs=sqrt(1.-(Rbody/Ri)*s2)
          FCiCs=(Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) - (1./7.)*(Ci**7-Cs**7)
          Vol=factor*FCiCs/s8

          where(MHD_SUM_vol(i,:)>1.1E-8)
             MHD_SUM_vol(i,:)=MHD_SUM_vol(i,:) + Vol
          end where
       end if

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values

          !Fix the grid inside Rbody
          MHD_Z0x(i,:) = (Ri/s2)*cos(RCM_lon(:)*(cPi/180.))
          MHD_Z0y(i,:) = (Ri/s2)*sin(RCM_lon(:)*(cPi/180.))

          !Now fix the volumes
          Cs=0.
          FCiCs=(Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) - (1./7.)*(Ci**7-Cs**7)
          Vol=factor*FCiCs/s8

          if( s2 > Ri/Rbody )then
             !Fieldline stays inside of Rbody, recompute exact volume
             MHD_SUM_vol(i,:)=Vol

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

    !set open fieldline values
    do j=1,jsize
       i0=1
       do i=2,isize
          if(MHD_SUM_vol(i,j)>1.1E-8 .and. MHD_SUM_vol(i,j)<MHD_SUM_vol(i-1,j))then
             i0=i-1
             exit
          end if
       end do

       !Fix grid overlaps
       i1=isize
       i2=isize
       do i=isize-1,1,-1
          if(i1==isize)then
             i1=i2
             i2=i
             CYCLE
          end if
          if(MHD_Z0x(i,j)<noValue+1. .or. MHD_Z0y(i,j)<noValue+1.)then
             MHD_Z0x(i,j) = 1.02*MHD_Z0x(i2,j)
             MHD_Z0y(i,j) = 1.02*MHD_Z0y(i2,j)
             i1=i2
             i2=i
             CYCLE
          end if
          i1=i2
          i2=i
       end do

       i=i0
       MHD_lat_boundary(j)=i

       MHD_SUM_vol(1:i-1,j)=MHD_SUM_vol(i,j)
       MHD_Z0T    (1:i-1,j)=MHD_Z0T    (i,j)
       MHD_Bmin   (1:i-1,j)=MHD_Bmin   (i,j)
       MHD_Z0x    (1:i-1,j)=MHD_Z0x    (i,j)
       MHD_Z0y    (1:i-1,j)=MHD_Z0y    (i,j)
    end do

    !Fill missing values
    do n=1,9
       do j=1,jsize; do i=isize,1,-1
          tmpValue=0.; tmpCount1=0.; tmpCount2=0.; tmpCount3=0.
          do i0=i-1,i+1; do j0=j-1,j+1
             if(i0/=i+1)then
                i1=i0; i2=i0+1; j1=j0; j2=j0
                i1=min(i1,isize); i1=max(i1,1)
                i2=min(i2,isize); i2=max(i1,1)
                if(j1<1) j1=j1+jsize; if(j1>jsize) j1=j1-jsize
                if(j2<1) j2=j2+jsize; if(j2>jsize) j2=j2-jsize

                if(MHD_SUM_vol(i1,j1)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i1,j1)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i1,j1)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i1,j1)
                   tmpCount1=tmpCount1+1.
                end if
                if(MHD_SUM_vol(i2,j2)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i2,j2)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i2,j2)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i2,j2)
                   tmpCount1=tmpCount1+1.
                end if

                if(MHD_Z0x(i1,j1)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i1,j1)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i1,j1)
                   tmpCount2=tmpCount2+1.
                end if
                if(MHD_Z0x(i2,j2)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i2,j2)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i2,j2)
                   tmpCount2=tmpCount2+1.
                end if

                if(MHD_Z0T(i1,j1)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i1,j1)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i1,j1)
                   tmpCount3=tmpCount3+1.
                end if
                if(MHD_Z0T(i2,j2)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i2,j2)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i2,j2)
                   tmpCount3=tmpCount3+1.
                end if

             end if
             if(j0/=j+1)then
                i1=i0; i2=i0; j1=j0; j2=j0+1
                i1=min(i1,isize); i1=max(i1,1)
                i2=min(i2,isize); i2=max(i1,1)
                if(j1<1) j1=j1+jsize; if(j1>jsize) j1=j1-jsize
                if(j2<1) j2=j2+jsize; if(j2>jsize) j2=j2-jsize

                if(MHD_SUM_vol(i1,j1)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i1,j1)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i1,j1)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i1,j1)
                   tmpCount1=tmpCount1+1.
                end if
                if(MHD_SUM_vol(i2,j2)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i2,j2)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i2,j2)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i2,j2)
                   tmpCount1=tmpCount1+1.
                end if

                if(MHD_Z0x(i1,j1)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i1,j1)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i1,j1)
                   tmpCount2=tmpCount2+1.
                end if
                if(MHD_Z0x(i2,j2)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i2,j2)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i2,j2)
                   tmpCount2=tmpCount2+1.
                end if

                if(MHD_Z0T(i1,j1)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i1,j1)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i1,j1)
                   tmpCount3=tmpCount3+1.
                end if
                if(MHD_Z0T(i2,j2)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i2,j2)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i2,j2)
                   tmpCount3=tmpCount3+1.
                end if

             end if
          end do; end do
          if(tmpCount1>0.) tmpValue(1:3)=tmpValue(1:3)/tmpCount1
          if(tmpCount2>0.) tmpValue(4:5)=tmpValue(4:5)/tmpCount2
          if(tmpCount3>0.) tmpValue(6:7)=tmpValue(6:7)/tmpCount3

          if(MHD_SUM_vol(i,j)<1.1E-8)MHD_SUM_vol(i,j)=tmpValue(1)
          if(MHD_SUM_rho(i,j)<1.1E-8)MHD_SUM_rho(i,j)=tmpValue(2)
          if(MHD_SUM_p  (i,j)<1.1E-8)MHD_SUM_p  (i,j)=tmpValue(3)
          if(MHD_Z0x    (i,j)<1.1E-8)MHD_Z0x    (i,j)=tmpValue(4)
          if(MHD_Z0y    (i,j)<1.1E-8)MHD_Z0y    (i,j)=tmpValue(5)
          if(MHD_Z0T    (i,j)<1.1E-8)MHD_Z0T    (i,j)=tmpValue(6)
          if(MHD_Bmin   (i,j)<1.1E-8)MHD_Bmin   (i,j)=tmpValue(7)
       end do; end do
    end do

    !Average values
    do n=1,1
       do j=1,jsize; do i=isize,1,-1
          tmpValue=0.; tmpCount1=0.; tmpCount2=0.; tmpCount3=0.
          do i0=i-1,i+1; do j0=j-1,j+1
             if(i0/=i+1)then
                i1=i0; i2=i0+1; j1=j0; j2=j0
                i1=min(i1,isize); i1=max(i1,1)
                i2=min(i2,isize); i2=max(i1,1)
                if(j1<1) j1=j1+jsize; if(j1>jsize) j1=j1-jsize
                if(j2<1) j2=j2+jsize; if(j2>jsize) j2=j2-jsize

                if(MHD_SUM_vol(i1,j1)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i1,j1)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i1,j1)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i1,j1)
                   tmpCount1=tmpCount1+1.
                end if
                if(MHD_SUM_vol(i2,j2)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i2,j2)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i2,j2)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i2,j2)
                   tmpCount1=tmpCount1+1.
                end if

                if(MHD_Z0x(i1,j1)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i1,j1)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i1,j1)
                   tmpCount2=tmpCount2+1.
                end if
                if(MHD_Z0x(i2,j2)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i2,j2)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i2,j2)
                   tmpCount2=tmpCount2+1.
                end if

                if(MHD_Z0T(i1,j1)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i1,j1)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i1,j1)
                   tmpCount3=tmpCount3+1.
                end if
                if(MHD_Z0T(i2,j2)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i2,j2)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i2,j2)
                   tmpCount3=tmpCount3+1.
                end if

             end if
             if(j0/=j+1)then
                i1=i0; i2=i0; j1=j0; j2=j0+1
                i1=min(i1,isize); i1=max(i1,1)
                i2=min(i2,isize); i2=max(i1,1)
                if(j1<1) j1=j1+jsize; if(j1>jsize) j1=j1-jsize
                if(j2<1) j2=j2+jsize; if(j2>jsize) j2=j2-jsize

                if(MHD_SUM_vol(i1,j1)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i1,j1)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i1,j1)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i1,j1)
                   tmpCount1=tmpCount1+1.
                end if
                if(MHD_SUM_vol(i2,j2)>1.1E-8)then
                   tmpValue(1)=tmpValue(1)+MHD_SUM_vol(i2,j2)
                   tmpValue(2)=tmpValue(2)+MHD_SUM_rho(i2,j2)
                   tmpValue(3)=tmpValue(3)+MHD_SUM_p  (i2,j2)
                   tmpCount1=tmpCount1+1.
                end if

                if(MHD_Z0x(i1,j1)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i1,j1)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i1,j1)
                   tmpCount2=tmpCount2+1.
                end if
                if(MHD_Z0x(i2,j2)>-200.)then
                   tmpValue(4)=tmpValue(4)+MHD_Z0x(i2,j2)
                   tmpValue(5)=tmpValue(5)+MHD_Z0y(i2,j2)
                   tmpCount2=tmpCount2+1.
                end if

                if(MHD_Z0T(i1,j1)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i1,j1)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i1,j1)
                   tmpCount3=tmpCount3+1.
                end if
                if(MHD_Z0T(i2,j2)>1.1E-8)then
                   tmpValue(6)=tmpValue(6)+MHD_Z0T(i2,j2)
                   tmpValue(7)=tmpValue(7)+MHD_Bmin(i2,j2)
                   tmpCount3=tmpCount3+1.
                end if

             end if
          end do; end do
          if(tmpCount1>0.) tmpValue(1:3)=tmpValue(1:3)/tmpCount1
          if(tmpCount2>0.) tmpValue(4:5)=tmpValue(4:5)/tmpCount2
          if(tmpCount3>0.) tmpValue(6:7)=tmpValue(6:7)/tmpCount3

          MHD_SUM_vol(i,j)=tmpValue(1)
          MHD_SUM_rho(i,j)=tmpValue(2)
          MHD_SUM_p  (i,j)=tmpValue(3)
          MHD_Z0x    (i,j)=tmpValue(4)
          MHD_Z0y    (i,j)=tmpValue(5)
          MHD_Z0T    (i,j)=tmpValue(6)
          MHD_Bmin   (i,j)=tmpValue(7)
       end do; end do
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -D-'
    end if

    !set values for open fieldlines
    do j=1,jsize
       MHD_SUM_vol(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_rho(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_p  (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Z0T    (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Bmin   (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Z0x    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Z0x(int(MHD_lat_boundary(j)),j)
       MHD_Z0y    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Z0y(int(MHD_lat_boundary(j)),j)
    end do

    !Overwrite computed fieldline volumes inside body
    Ri=(6378.+100.)/6378.
    factor=(unitSI_B*1.E+9) * 2. * (Ri**4) / 31100.
    do i=1,isize
       colat = cPi/2. - RCM_lat(i)*(cPi/180.)
       s2=(sin(colat))**2
       s8=(sin(colat))**8
       Ci=cos(colat)

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values
          Cs=0.
          FCiCs=(Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) - (1./7.)*(Ci**7-Cs**7)
          Vol=factor*FCiCs/s8

          if( s2 > Ri/Rbody )then
             !Fieldline stays inside of Rbody, recompute exact volume
             MHD_SUM_vol(i,:)=Vol
          else
             !Fieldline stays inside of Rcurrents but outside Rbody
             factor1=( 1. / ((Ri/Rbody)-(Ri/Rcurrents)) )*((Ri/Rbody)-s2)
             factor2=( 1. / ((Ri/Rbody)-(Ri/Rcurrents)) )*(s2-(Ri/Rcurrents))

             !Blend volume with exact volume
             MHD_SUM_vol(i,:)= factor1*MHD_SUM_vol(i,:) + factor2*Vol
          end if
       end if
    end do

    !dimensionalize values
    where(MHD_SUM_vol0>1.1E-8)
       MHD_SUM_vol0 = MHD_SUM_vol0 / unitSI_B
    elsewhere
       MHD_SUM_vol0 = -1.
    end where

    where(MHD_SUM_vol>1.1E-8)
       MHD_SUM_vol  = MHD_SUM_vol  / unitSI_B
       MHD_SUM_rho  = MHD_SUM_rho * unitSI_rho
       MHD_SUM_p    = MHD_SUM_p * unitSI_p
    elsewhere
       MHD_SUM_vol  = -1.
       MHD_SUM_rho  = -1.
       MHD_SUM_p    = -1.
    end where

    where(MHD_Z0T>1.1E-8)
       MHD_Z0T      = MHD_Z0T * unitSI_temperature
       MHD_Bmin     = MHD_Bmin * unitSI_B
    elsewhere
       MHD_Z0T      = -1.
       MHD_Bmin     = -1.
    end where

  end subroutine process_integrated_data

!^CFG END RCM
end subroutine GM_get_for_im

!==========================================================================
!==========================================================================
!==========================================================================
!^CFG IF RCM BEGIN
subroutine getMinMaxLatLon(minLat,maxLat, minLon,maxLon, rayface_2d, shifted)


  use ModMain, ONLY : nI,nJ
  implicit none

  real, intent(inout) :: minLat,maxLat, minLon,maxLon
  real, intent(inout), dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d
  logical, intent(inout) :: shifted

  shifted=.false.
  minLat=minval(rayface_2d(1,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  maxLat=maxval(rayface_2d(1,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  minLon=minval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  maxLon=maxval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
       MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  if(maxLon-minLon > 180.)then
     shifted=.true.
     where(rayface_2d(2,1,1:nI+1,1:nJ+1)<180.) &
          rayface_2d(2,1,1:nI+1,1:nJ+1)=rayface_2d(2,1,1:nI+1,1:nJ+1)+360.
     minLon=minval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
          MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
     maxLon=maxval(rayface_2d(2,1,1:nI+1,1:nJ+1), &
          MASK=   (rayface_2d(3,1,1:nI+1,1:nJ+1)==3.))
  end if
end subroutine getMinMaxLatLon

!==========================================================================
!==========================================================================
!==========================================================================
subroutine find_eqT(iBLK,xL,yL,zL,eqT,eqB)
  use ModMain, ONLY : nI,nJ,nK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModAdvance, ONLY : State_VGB, rho_, p_, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in) :: xL,yL,zL
  real, intent(out) :: eqT,eqB

  integer :: i,j,k

  do i=1,nI; do j=1,nJ; do k=1,nK
     if(  abs(x_BLK(i,j,k,iBLK)-xL) <= 0.5*dx_BLK(iBLK) .and. &
          abs(y_BLK(i,j,k,iBLK)-yL) <= 0.5*dy_BLK(iBLK) .and. &
          abs(z_BLK(i,j,k,iBLK)-zL) <= 0.5*dz_BLK(iBLK) ) then
        eqT=State_VGB(p_,i,j,k,iBLK)/State_VGB(rho_,i,j,k,iBLK)
        eqB=sqrt( &
             (State_VGB(Bx_,i,j,k,iBLK)+B0xCell_BLK(i,j,k,iBLK))**2 + &
             (State_VGB(By_,i,j,k,iBLK)+B0yCell_BLK(i,j,k,iBLK))**2 + &
             (State_VGB(Bz_,i,j,k,iBLK)+B0zCell_BLK(i,j,k,iBLK))**2 )
     end if
  end do; end do; end do

end subroutine find_eqT
!==========================================================================
subroutine find_rayface_location(dbg,iBLK,rLat,rLon,xL,yL,rayface_2d)
  use ModProcMH
  use ModMain, ONLY : nI,nJ
  use ModRaytrace
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in) :: rLat,rLon
  real, intent(inout) :: xL,yL
  real, intent(in), dimension(3,2,1:nI+1,1:nJ+1):: rayface_2d
  logical, intent(in) :: dbg

  integer, parameter :: k=1
  integer :: i,j, iNear,jNear,nNear,nNear2, idir,jdir
  real :: gradLatI,gradLatJ,gradLatR2, gradLonI,gradLonJ,gradLonR2
  real :: xp1,yp1,Lat1,Lon1, xp2,yp2,Lat2,Lon2

  logical, dimension(nI+1,nJ+1,2) :: Lnear
  logical, dimension(nI+1,nJ+1) :: Tnear

  Lnear = .false.
  Tnear = .false.

  xL=0.; yL=0.

  do i=1,nI+1; do j=1,nJ+1
     if(rayface_2d(3,1,i,j)==3.) then
        if(abs(rayface_2d(1,1,i,j)-rLat)<1.E-4)then
           Lnear(i,j,1)=.true.
           if(i<nI+1)then
              if(rayface_2d(3,1,i+1,j)==3.) Lnear(i+1,j,1)=.true.
           end if
           if(j<nJ+1)then
              if(rayface_2d(3,1,i,j+1)==3.) Lnear(i,j+1,1)=.true.
           end if
        else
           if(i<nI+1)then
              if(rayface_2d(3,1,i+1,j)==3.)then
                 if(  (rayface_2d(1,1,i+1,j)-rLat)/ &
                      (rayface_2d(1,1,i  ,j)-rLat) <= 0. )then
                    Lnear(i:i+1,j,1)=.true.
                 end if
              end if
           end if
           if(j<nJ+1)then
              if(rayface_2d(3,1,i,j+1)==3.)then
                 if(  (rayface_2d(1,1,i,j+1)-rLat)/ &
                      (rayface_2d(1,1,i,j  )-rLat) <= 0. )then
                    Lnear(i,j:j+1,1)=.true.
                 end if
              end if
           end if
        end if

        if(abs(rayface_2d(2,1,i,j)-rLon)<1.E-4)then
           Lnear(i,j,2)=.true.
           if(i<nI+1)then
              if(rayface_2d(3,1,i+1,j)==3.) Lnear(i+1,j,2)=.true.
           end if
           if(j<nJ+1)then
              if(rayface_2d(3,1,i,j+1)==3.) Lnear(i,j+1,2)=.true.
           end if
        else
           if(i<nI+1)then
              if(rayface_2d(3,1,i+1,j)==3.)then
                 if(  (rayface_2d(2,1,i+1,j)-rLon)/ &
                      (rayface_2d(2,1,i  ,j)-rLon) <= 0. )then
                    Lnear(i:i+1,j,2)=.true.
                 end if
              end if
           end if
           if(j<nJ+1)then
              if(rayface_2d(3,1,i,j+1)==3.)then
                 if(  (rayface_2d(2,1,i,j+1)-rLon)/ &
                      (rayface_2d(2,1,i,j  )-rLon) <= 0. )then
                    Lnear(i,j:j+1,2)=.true.
                 end if
              end if
           end if
        end if
     end if
  end do; end do

  Tnear = (Lnear(:,:,1) .and. Lnear(:,:,2))

  if(dbg)then
     write(*,*)' '
     write(*,*)'PE=',iProc,' BLK=',iBLK
     write(*,*)'  Searching for value at Lat=',rLat,' Lon=',rLon
     write(*,*)' '
     write(*,*)' Lnear-1'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1-5 K=1  Lnear(i,j,1)=',Lnear(i,:,1)
     end do
     write(*,*)' Lnear-2'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1-5 K=1  Lnear(i,j,2)=',Lnear(i,:,2)
     end do
     write(*,*)' Tnear'
     do i=1,nI+1
        write(*,*) 'I=',i,' J=1-5 K=1  Tnear(i,j)=  ',Tnear(i,:)
     end do
     write(*,*)' X'
     do i=1,nI+1
        write(*,*) &
             'I=',i,' J=1-5 K=1  rayface_2d(1,2,i,j)=',rayface_2d(1,2,i,:)
     end do
     write(*,*)' Y'
     do i=1,nI+1
        write(*,*) &
             'I=',i,' J=1-5 K=1  rayface_2d(2,2,i,j)=',rayface_2d(2,2,i,:)
     end do
     write(*,*)' Lat'
     do i=1,nI+1
        write(*,*) &
             'I=',i,' J=1-5 K=1  rayface_2d(1,1,i,j)=',rayface_2d(1,1,i,:)
     end do
     write(*,*)' Lon'
     do i=1,nI+1
        write(*,*) &
             'I=',i,' J=1-5 K=1  rayface_2d(2,1,i,j)=',rayface_2d(2,1,i,:)
     end do
     write(*,*)' Status'
     do i=1,nI+1
        write(*,*) &
             'I=',i,' J=1-5 K=1  rayface_2d(3,1,i,j)=',rayface_2d(3,1,i,:)
     end do
     write(*,*)' '
  end if

  !Return if no points found
  if(.not.any(Tnear)) return

  !Find nearest from Tnear list
  iNear=0
  jNear=0
  nNear=0
  do i=1,nI+1; do j=1,nJ+1
     if(Tnear(i,j))then
        nNear2=1
        if(i>1)then
           if(Tnear(i-1,j))   nNear2=nNear2+8
           if(Lnear(i-1,j,1)) nNear2=nNear2+1
           if(Lnear(i-1,j,2)) nNear2=nNear2+1
        else
           nNear2=nNear2+3
        end if
        if(i<nI+1)then
           if(Tnear(i+1,j))   nNear2=nNear2+8
           if(Lnear(i+1,j,1)) nNear2=nNear2+1
           if(Lnear(i+1,j,2)) nNear2=nNear2+1
        else
           nNear2=nNear2+3
        end if
        if(j>1)then
           if(Tnear(i,j-1))   nNear2=nNear2+8
           if(Lnear(i,j-1,1)) nNear2=nNear2+1
           if(Lnear(i,j-1,2)) nNear2=nNear2+1
        else
           nNear2=nNear2+3
        end if
        if(j<nJ+1)then
           if(Tnear(i,j+1))   nNear2=nNear2+8
           if(Lnear(i,j+1,1)) nNear2=nNear2+1
           if(Lnear(i,j+1,2)) nNear2=nNear2+1
        else
           nNear2=nNear2+3
        end if
        if(nNear2>nNear)then
           nNear=nNear2
           iNear=i
           jNear=j
        end if
     end if
  end do; end do

  if(iNear==1)then
     if(rayface_2d(3,1,iNear+1,jNear)/=3.) return
     idir= 1
  elseif(iNear==nI+1)then
     if(rayface_2d(3,1,iNear-1,jNear)/=3.) return
     idir=-1
  else
     if    (Tnear(iNear+1,jNear))then
        idir= 1
     elseif(Tnear(iNear-1,jNear))then
        idir=-1
     elseif(any(Lnear(iNear+1,jNear,:)))then
        idir= 1
     elseif(any(Lnear(iNear-1,jNear,:)))then
        idir=-1
     elseif(rayface_2d(3,1,iNear+1,jNear)==3.)then
        idir= 1
     elseif(rayface_2d(3,1,iNear-1,jNear)==3.)then
        idir=-1
     else
        return
     end if
  end if

  if(jNear==1)then
     if(rayface_2d(3,1,iNear,jNear+1)/=3.) return
     jdir= 1
  elseif(jNear==nJ+1)then
     if(rayface_2d(3,1,iNear,jNear-1)/=3.) return
     jdir=-1
  else
     if    (Tnear(iNear,jNear+1))then
        jdir= 1
     elseif(Tnear(iNear,jNear-1))then
        jdir=-1
     elseif(any(Lnear(iNear,jNear+1,:)))then
        jdir= 1
     elseif(any(Lnear(iNear,jNear-1,:)))then
        jdir=-1
     elseif(rayface_2d(3,1,iNear,jNear+1)==3.)then
        jdir= 1
     elseif(rayface_2d(3,1,iNear,jNear-1)==3.)then
        jdir=-1
     else
        return
     end if
  end if

  if(dbg)then
     write(*,*)'Nearest Point:',iNear,jNear
     write(*,*)'dirs: ',idir,jdir
  end if

  !compute gradients
  if(abs(rayface_2d(1,1,iNear+idir,jNear)-rayface_2d(1,1,iNear,jNear))<0.01)then
     gradLatI=(rLat-rayface_2d(1,1,iNear,jNear))/0.01
  else
     gradLatI=(rLat-rayface_2d(1,1,iNear,jNear))/ &
          (rayface_2d(1,1,iNear+idir,jNear)-rayface_2d(1,1,iNear,jNear))
  end if
  if(abs(rayface_2d(1,1,iNear,jNear+jdir)-rayface_2d(1,1,iNear,jNear))<0.01)then
     gradLatJ=(rLat-rayface_2d(1,1,iNear,jNear))/0.01
  else
     gradLatJ=(rLat-rayface_2d(1,1,iNear,jNear))/ &
          (rayface_2d(1,1,iNear,jNear+jdir)-rayface_2d(1,1,iNear,jNear))
  end if
  gradLatR2=gradLatI**2 + gradLatJ**2

  if(abs(rayface_2d(2,1,iNear+idir,jNear)-rayface_2d(2,1,iNear,jNear))<0.01)then
     gradLonI=(rLon-rayface_2d(2,1,iNear,jNear))/0.01
  else
     gradLonI=(rLon-rayface_2d(2,1,iNear,jNear))/ &
          (rayface_2d(2,1,iNear+idir,jNear)-rayface_2d(2,1,iNear,jNear))
  end if
  if(abs(rayface_2d(2,1,iNear,jNear+jdir)-rayface_2d(2,1,iNear,jNear))<0.01)then
     gradLonJ=(rLon-rayface_2d(2,1,iNear,jNear))/0.01
  else
     gradLonJ=(rLon-rayface_2d(2,1,iNear,jNear))/ &
          (rayface_2d(2,1,iNear,jNear+jdir)-rayface_2d(2,1,iNear,jNear))
  end if
  gradLonR2=gradLonI**2 + gradLonJ**2

  !assign final value
  if(gradLatR2>gradLonR2)then
     xp1=rayface_2d(1,2,iNear,jNear) + &
          gradLatI*(rayface_2d(1,2,iNear+idir,jNear)-rayface_2d(1,2,iNear,jNear))
     yp1=rayface_2d(2,2,iNear,jNear) + &
          gradLatI*(rayface_2d(2,2,iNear+idir,jNear)-rayface_2d(2,2,iNear,jNear))
     Lon1=rayface_2d(2,1,iNear,jNear) + &
          gradLatI*(rayface_2d(2,1,iNear+idir,jNear)-rayface_2d(2,1,iNear,jNear))

     xp2=rayface_2d(1,2,iNear,jNear) + &
          gradLatJ*(rayface_2d(1,2,iNear,jNear+jdir)-rayface_2d(1,2,iNear,jNear))
     yp2=rayface_2d(2,2,iNear,jNear) + &
          gradLatJ*(rayface_2d(2,2,iNear,jNear+jdir)-rayface_2d(2,2,iNear,jNear))
     Lon2=rayface_2d(2,1,iNear,jNear) + &
          gradLatJ*(rayface_2d(2,1,iNear,jNear+jdir)-rayface_2d(2,1,iNear,jNear))

     if(abs(Lon2-Lon1)<01.E-4)then
!!$        write(*,*)'WARNING: Small difference found (Lon):',Lon2-Lon1
!!$        call debugging_and_stop
        return
     end if
     xL=xp1+((rLon-Lon1)/(Lon2-Lon1))*(xp2-xp1)
     yL=yp1+((rLon-Lon1)/(Lon2-Lon1))*(yp2-yp1)
  else
     xp1=rayface_2d(1,2,iNear,jNear) + &
          gradLonI*(rayface_2d(1,2,iNear+idir,jNear)-rayface_2d(1,2,iNear,jNear))
     yp1=rayface_2d(2,2,iNear,jNear) + &
          gradLonI*(rayface_2d(2,2,iNear+idir,jNear)-rayface_2d(2,2,iNear,jNear))
     Lat1=rayface_2d(1,1,iNear,jNear) + &
          gradLonI*(rayface_2d(1,1,iNear+idir,jNear)-rayface_2d(1,1,iNear,jNear))

     xp2=rayface_2d(1,2,iNear,jNear) + &
          gradLonJ*(rayface_2d(1,2,iNear,jNear+jdir)-rayface_2d(1,2,iNear,jNear))
     yp2=rayface_2d(2,2,iNear,jNear) + &
          gradLonJ*(rayface_2d(2,2,iNear,jNear+jdir)-rayface_2d(2,2,iNear,jNear))
     Lat2=rayface_2d(1,1,iNear,jNear) + &
          gradLonJ*(rayface_2d(1,1,iNear,jNear+jdir)-rayface_2d(1,1,iNear,jNear))

     if(abs(Lat2-Lat1)<01.E-4)then
!!$        write(*,*)'WARNING: Small difference found (Lat):',Lat2-Lat1
!!$        call debugging_and_stop
        return
     end if
     xL=xp1+((rLat-Lat1)/(Lat2-Lat1))*(xp2-xp1)
     yL=yp1+((rLat-Lat1)/(Lat2-Lat1))*(yp2-yp1)
  end if

  if(dbg)then
     write(*,*)'gradLat:',gradLatI,gradLatJ
     write(*,*)'gradLon:',gradLonI,gradLonJ
     write(*,*)' point1:',xp1,yp1,Lat1,Lon1
     write(*,*)' point2:',xp2,yp2,Lat2,Lon2
     write(*,*)' tentative location: ',xL,yL
  end if

  !nudge to block edge
  if(abs(xL-rayface_2d(1,2,   1,1))<1.E-4) xL=rayface_2d(1,2,   1,1)
  if(abs(xL-rayface_2d(1,2,nI+1,1))<1.E-4) xL=rayface_2d(1,2,nI+1,1)
  if(abs(yL-rayface_2d(2,2,1,   1))<1.E-4) yL=rayface_2d(2,2,1,   1)
  if(abs(yL-rayface_2d(2,2,1,nJ+1))<1.E-4) yL=rayface_2d(2,2,1,nJ+1)

  !check that value is in block
  if(  xL<=rayface_2d(1,2,   1,1) .or. &
       xL> rayface_2d(1,2,nI+1,1) .or. &
       yL<=rayface_2d(2,2,1,   1) .or. &
       yL> rayface_2d(2,2,1,nJ+1) )then
     if(dbg)write(*,*)'Not in block, returning:', &
          '  x=',xL,rayface_2d(1,2,1,1),rayface_2d(1,2,nI+1,nJ+1), &
          '  y=',yL,rayface_2d(2,2,1,1),rayface_2d(2,2,nI+1,nJ+1)
     xL=0.
     yL=0.
     return
  end if

  if(dbg)then
     write(*,*)' '
     write(*,*)' Final Location: ',xL,yL
     write(*,*)' '
  end if

contains

  subroutine debugging_and_stop
    write(*,*)' '
    write(*,*)'PE=',iProc,' BLK=',iBLK
    write(*,*)'  Searching for value at Lat=',rLat,' Lon=',rLon
    write(*,*)' '
    write(*,*)' Lnear-1'
    do i=1,nI+1
       write(*,*) 'I=',i,' J=1-5 K=1  Lnear(i,j,1)=',Lnear(i,:,1)
    end do
    write(*,*)' Lnear-2'
    do i=1,nI+1
       write(*,*) 'I=',i,' J=1-5 K=1  Lnear(i,j,2)=',Lnear(i,:,2)
    end do
    write(*,*)' Tnear'
    do i=1,nI+1
       write(*,*) 'I=',i,' J=1-5 K=1  Tnear(i,j)=  ',Tnear(i,:)
    end do
    write(*,*)' X'
    do i=1,nI+1
       write(*,'(a,i1,a,5f10.3)') &
            'I=',i,' J=1-5 K=1  rayface_2d(1,2,i,j)=',rayface_2d(1,2,i,:)
    end do
    write(*,*)' Y'
    do i=1,nI+1
       write(*,'(a,i1,a,5f10.3)') &
            'I=',i,' J=1-5 K=1  rayface_2d(2,2,i,j)=',rayface_2d(2,2,i,:)
    end do
    write(*,*)' Lat'
    do i=1,nI+1
       write(*,'(a,i1,a,5f10.3)') &
            'I=',i,' J=1-5 K=1  rayface_2d(1,1,i,j)=',rayface_2d(1,1,i,:)
    end do
    write(*,*)' Lon'
    do i=1,nI+1
       write(*,'(a,i1,a,5f10.3)') &
            'I=',i,' J=1-5 K=1  rayface_2d(2,1,i,j)=',rayface_2d(2,1,i,:)
    end do
    write(*,*)' Status'
    do i=1,nI+1
       write(*,'(a,i1,a,5f10.1)') &
            'I=',i,' J=1-5 K=1  rayface_2d(3,1,i,j)=',rayface_2d(3,1,i,:)
    end do
    write(*,*)' '
    write(*,*)'Nearest Point:',iNear,jNear
    write(*,*)'dirs: ',idir,jdir
    write(*,*)'gradLat:',gradLatI,gradLatJ
    write(*,*)'gradLon:',gradLonI,gradLonJ
    write(*,*)' point1:',xp1,yp1,Lat1,Lon1
    write(*,*)' point2:',xp2,yp2,Lat2,Lon2
    write(*,*)' tentative location: ',xL,yL

    call stop_mpi("debug-stop")
  end subroutine debugging_and_stop

end subroutine find_rayface_location
!==========================================================================
!^CFG END RCM
