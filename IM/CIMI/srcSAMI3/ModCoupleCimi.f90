Module ModCoupleCimi
    implicit none 
    private !except

    logical, public :: DoCoupleCimi = .false.
    
    integer, public ::  iStartTime_I(7)  =(/1976,6,28,0,0,0,0/)
 
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal
    
    ! cimi Latitude (radians) and lon (radians) grid
    real, allocatable :: LatCimi_C(:),LonCimi_C(:)
    integer :: nLatCIMI, nLonCIMI
    
    ! cimi potential on cimi grid [V]
    real, allocatable :: PotCIMI_C(:,:)
    
    !cimi potential on sami grid [statV]
    real, public, allocatable :: PotCimiOnSamiGrid_C(:,:)
    
    ! sami grid parameters for this module
    integer :: nLatSami, nLonSami
    real,allocatable    :: bLatSami(:,:), bLonSami(:,:)

    ! Electron density at magnetic equator from SAMI3
    real,allocatable :: PlasEq_C(:,:)
    ! H+ density at magnetic equator from SAMI3
    real,allocatable :: PlasHpEq_C(:,:)
    ! He+ density at magnetic equator from SAMI3
    real,allocatable :: PlasHepEq_C(:,:)
    ! O+ density at magnetic equator from SAMI3
    real,allocatable :: PlasOpEq_C(:,:)
    
    
    !public methods
    public :: sami_set_global_mpi
    public :: sami_put_init_from_cimi
    public :: sami_get_from_cimi
    public :: set_sami_grid_for_mod
    public :: sami_send_to_cimi
    public :: sami_get_init_for_cimi
  contains
    !==========================================================================
    ! call by all procs
    subroutine sami_set_global_mpi(iProc0CimiIN,iProc0SamiIN, iCommGlobalIN)
      integer, intent(in) ::iProc0CimiIN,iProc0SamiIN, iCommGlobalIN
      !------------------------------------------------------------------------
      iProc0CIMI=iProc0CimiIN
      iProc0Sami=iProc0SamiIN
      iCommGlobal=iCommGlobalIN

      ! Set the coupling to be true
      DoCoupleCimi=.true.
    end subroutine sami_set_global_mpi
    !========================================================================
    subroutine sami_put_init_from_cimi
      use ModSAMI, ONLY: iComm,iProc
      use ModMpi
      use ModTimeConvert, ONLY: time_int_to_real
      use CON_axes,         ONLY: init_axes
      integer ::  iStartTimeCIMI_I(7),iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      ! real version of start time for axes
      real :: StartTime
      !----------------------------------------------------------------------

      ! recieve and then bcast the starttime from cimi
      if (iProc == 0) call MPI_recv(iStartTime_I,7,MPI_INTEGER,iProc0CIMI,&
           1,iCommGlobal,iStatus_I,iError)
      call MPI_bcast(iStartTime_I,7,MPI_LOGICAL,0,iComm,iError)

      !set startime and axes
      call time_int_to_real(iStartTime_I,StartTime)
      !\
      ! Set axes for coord transform 
      !/
      call init_axes(StartTime)
      
      ! recieve grid size
      if (iProc == 0) call MPI_recv(nLatCIMI,1,MPI_INTEGER,iProc0CIMI,&
           2,iCommGlobal,iStatus_I,iError)

      if (iProc == 0) call MPI_recv(nLonCIMI,1,MPI_INTEGER,iProc0CIMI,&
           3,iCommGlobal,iStatus_I,iError)
      
      ! bcast grid size
      call MPI_bcast(nLatCIMI,1,MPI_LOGICAL,0,iComm,iError)
      call MPI_bcast(nLonCIMI,1,MPI_LOGICAL,0,iComm,iError)
      
      ! allocate CIMI grid
      allocate(LatCimi_C(nLatCIMI))
      allocate(LonCimi_C(nLonCIMI))

      if (iProc == 0) call MPI_recv(LatCimi_C,nLatCIMI,MPI_REAL,iProc0CIMI,&
           4,iCommGlobal,iStatus_I,iError)

      if (iProc == 0) call MPI_recv(LonCimi_C,nLonCIMI,MPI_REAL,iProc0CIMI,&
           5,iCommGlobal,iStatus_I,iError)

      ! allocate arrays that hold data from cimi
      allocate(PotCIMI_C(nLatCIMI,nLonCIMI))
      
    end subroutine sami_put_init_from_cimi
    
    !==========================================================================
    ! 
    subroutine sami_get_from_cimi(TimeSimulation)
      use ModMPI
      use ModSAMI, ONLY: iComm,iProc
      real, intent(in) :: TimeSimulation
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      !if(iProc ==0) write(*,*) 'sami_get_from_cimi',nLatCIMI,nLonCIMI
      ! Send the starttime from cimi to sami
      if(iProc ==0) then
         call MPI_recv(PotCIMI_C,nLatCIMI*nLonCIMI,MPI_REAL,&
              iProc0CIMI,1,iCommGlobal,iStatus_I,iError)
         !KLUDGE temporarily overwrite CIMI potential with 0 to test rest of 
         ! coupling
         !PotCIMI_C(:,:)=0.0
      endif
      ! write out max and min of pot
      if(iProc ==0) write(*,*) 'Max and min of CIMI Potential [V]: ', &
           maxval(PotCIMI_C), minval(PotCIMI_C)
      
      ! interpolate cimi potential to sami grid
      if(iProc ==0) call interpolate_cimi_to_sami(TimeSimulation)
      
      if(iProc ==0) &
           write(*,*) 'Max and min of CIMI Potential on Sami Grid [statV]: ', &
           maxval(PotCimiOnSamiGrid_C), minval(PotCimiOnSamiGrid_C)

      ! now bcast cimi potential on sami grid to all procs
      call MPI_bcast(PotCimiOnSamiGrid_C,nLatSami*nLonSami, &
           MPI_REAL,0,iComm,iError)
      
    end subroutine sami_get_from_cimi

    !==========================================================================
    subroutine interpolate_cimi_to_sami(TimeSimulation)
      use ModInterpolate, ONLY: bilinear
      use ModNumConst,    ONLY: cDegToRad
      real, intent(in) :: TimeSimulation
      integer :: iLat, iLon
      real :: LatLon_D(2)
      real,parameter :: cVToStatV = 0.00333564095198
      do iLon = 1, nLonSami-1
         do iLat = 1, nLatSami
            !write(*,*) 'iLat,iLon,cDegToRad',iLat,iLon,cDegToRad
            !LatLon_D(1) = bLatSami(iLat,iLon) * cDegToRad
            !LatLon_D(2) = bLonSami(iLat,iLon) * cDegToRad
            
            !convert point on SAMI grid to CIMI coords for interpolation 
!            call convert_sami_to_cimi_lat_lon(TimeSimulation, &
!                 bLatSami(iLat,iLon),bLonSami(iLat,iLon), &
!                 LatLon_D(1), LatLon_D(2))

            ! Now assume SAMI in SM coords but azimuthal angle 0 defined at 
            ! midnight. For CIMI 0 is at noon so 180 degree rotation.
            LatLon_D(1) = bLatSami(iLat,iLon)     * cDegToRad
            LatLon_D(2) = modulo(bLonSami(iLat,iLon)+180.0,360.0) * cDegToRad
            
            !deal with coord transformation here

            if (LatLon_D(1) > maxval(LatCimi_C) .or. &
                 LatLon_D(1) < minval(LatCimi_C) ) then
               PotCimiOnSamiGrid_C(iLon,iLat) = 0.0
            else
               PotCimiOnSamiGrid_C(iLon,iLat) = &
                    bilinear(PotCIMI_C,1,nLatCIMI,1,nLonCIMI,LatLon_D, &
                    LatCimi_C,LonCimi_C,DoExtrapolate=.true.)
               !convert from V to stat volts
               PotCimiOnSamiGrid_C(iLon,iLat) = &
                    PotCimiOnSamiGrid_C(iLon,iLat)*cVToStatV
            endif
         enddo
      enddo
      
      ! add the ghost cell for potential
      PotCimiOnSamiGrid_C(nLonSami,:) = PotCimiOnSamiGrid_C(1,:) 
      !call plot_coupled_values
      

    end subroutine interpolate_cimi_to_sami

    !==========================================================================
    subroutine set_sami_grid_for_mod(nzp1,nfp1,nlt,nnx,nny,blatpt,blonpt)
      integer, intent(in) :: nzp1, nfp1, nlt, nnx, nny
      real,    intent(in) :: blatpt(nzp1,nfp1,nlt),blonpt(nzp1,nfp1,nlt)


      allocate(bLatSami(nny,nnx-1),bLonSami(nny,nnx-1))
      
      nLonSami = nnx
      nLatSami = nny
      bLatSami = blatpt(nzp1,1:nny,1:nnx-1)
      bLonSami = blonpt(nzp1,1:nny,1:nnx-1)

      allocate(PotCimiOnSamiGrid_C(nLonSami,nLatSami))
    end subroutine set_sami_grid_for_mod

    !==========================================================================
    ! subroutine that takes SAMI lat lon and returns CIMI lat lon
    subroutine convert_sami_to_cimi_lat_lon(TimeSimulation, SamiLat,SamiLon, &
         CimiLat, CimiLon)
      use CON_axes, ONLY: transform_matrix
      use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg
      
      real, intent(in) :: TimeSimulation, SamiLat, SamiLon
      real, intent(out):: CimiLat, CimiLon
      real             :: SamiToCimi_DD(3,3), theta, phi 
      real             :: XyzCimi_D(3), XyzSami_D(3)
      character(len=*),  parameter :: NameCimiCoord='SMG'
      character(len=*),  parameter :: NameSamiCoord='MAG'
      !------------------------------------------------------------------------
      
      ! Get polar angle (theta) and azimuthal angle (phi)
      theta = 0.5*cPi - SamiLat*cDegToRad
      phi   = SamiLon*cDegToRad
      
      ! Get xyzCimi_D from CimiLat and CimiLon
      xyzSami_D(1) = sin(theta)*cos(phi)
      xyzSami_D(2) = sin(theta)*sin(phi)
      xyzSami_D(3) = cos(theta)
      
      !\
      ! get equivalent geographic coords Hemisphere 1
      !/
      
      ! Get transform matrix 
      SamiToCimi_DD = &
           transform_matrix(TimeSimulation, NameSamiCoord, NameCimiCoord)
      
      ! Transform xyzCimi_D to XyzGg_DI
      XyzCimi_D = matmul( SamiToCimi_DD, XyzSami_D)
      
      ! Calculate CimiLat and CimiLon 
      CimiLon = modulo(atan2(XyzCimi_D(2), XyzCimi_D(1)), cTwoPi) * cRadToDeg
      CimiLat = 90.0 - (acos(max(-1.0,min(1.0, XyzCimi_D(3))))*cRadToDeg)
      
      
    end subroutine convert_sami_to_cimi_lat_lon

    !==========================================================================
    ! only call by proc 0
    subroutine sami_get_init_for_cimi
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      ! Send the grid parameters from sami to cimi
      call MPI_send(nLatSami,1,MPI_INTEGER,iProc0CIMI,&
           1,iCommGlobal,iError)

      call MPI_send(nLonSami,1,MPI_INTEGER,iProc0CIMI,&
           2,iCommGlobal,iError)
      
      call MPI_send(bLatSami(1:nLatSami,1),nLatSami,MPI_REAL,iProc0CIMI,&
           3,iCommGlobal,iError)

      call MPI_send(bLonSami(1,1:nLonSami-1),nLonSami-1,MPI_REAL,iProc0CIMI,&
           4,iCommGlobal,iError)

    end subroutine sami_get_init_for_cimi
    !==========================================================================
    ! 
    subroutine sami_send_to_cimi
      use ModSAMI, ONLY: iProc
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      !gather plasmaspheric densities for CIMI
      if (iProc == 0 .and. .not.allocated(PlasEq_C)) &
           allocate(PlasEq_C(nLatSami,nLonSami))
      if (iProc == 0 .and. .not.allocated(PlasHpEq_C)) &
           allocate(PlasHpEq_C(nLatSami,nLonSami))
      if (iProc == 0 .and. .not.allocated(PlasHepEq_C)) &
           allocate(PlasHepEq_C(nLatSami,nLonSami))
      if (iProc == 0 .and. .not.allocated(PlasOpEq_C)) &
           allocate(PlasOpEq_C(nLatSami,nLonSami))
      if(iProc ==0) write(*,*) 'iProc',iProc,'calling gather_dene'      
      call gather_dene
      !if(iProc ==0) PlasEq_C(:,:) =0.0
      if(iProc ==0) write(*,*) 'iProc',iProc,'finished gather_dene'      
      !if(iProc==0)  write(*,*) 'cimi_send_to_sami',nLat,nLon
      ! Send the starttime from cimi to sami
      if(iProc ==0) then
         call MPI_send(PlasEq_C,nLatSami*nLonSami,MPI_REAL, &
              iProc0CIMI,1,iCommGlobal,iError)
         call MPI_send(PlasHpEq_C,nLatSami*nLonSami,MPI_REAL, &
              iProc0CIMI,1,iCommGlobal,iError)
         call MPI_send(PlasHepEq_C,nLatSami*nLonSami,MPI_REAL, &
              iProc0CIMI,1,iCommGlobal,iError)
         call MPI_send(PlasOpEq_C,nLatSami*nLonSami,MPI_REAL, &
              iProc0CIMI,1,iCommGlobal,iError)
      end if
      if(iProc ==0) write(*,*) 'SAMI3 iProc=',iProc,'sending PlasEq_C to cimi'
    end subroutine sami_send_to_cimi

    !=========================================================================  
    subroutine gather_dene
      use ModSAMI, ONLY: iComm,iProc,nProc
      use ModMPI
      include 'param3_mpi-1.98.inc'
!      include 'com3_mpi-1.98.inc' 
      integer :: ierr
      real denitmp(nz,nf,nl)
      integer :: nntmp, itmp, jtmp,ktmp, iwrk,nn,k,kk,i,j,ni
      integer :: status(MPI_STATUS_SIZE)
      real, allocatable :: denet(:,:,:),denit(:,:,:,:)
      integer,parameter :: nion1=1, nion2=7
      ! variables we need for the commonblocks
      real deni(nz,nf,nl,nion),ne(nz,nf,nl),denn(nz,nf,nl,nneut)
      real denni(nz,nf,nl,nneut),dennf(nz,nf,nl,nneut)
      real vsi(nz,nf,nl,nion),vsid(nz,nf,nl,nion)
      real sumvsi(nz,nf,nl,nion),vsic(nz,nf,nl,nion)
      real te(nz,nf,nl),ti(nz,nf,nl,nion),tn(nz,nf,nl)
      real tni(nz,nf,nl),tnf(nz,nf,nl)
      real u(nz,nf,nl),v(nz,nf,nl),vpi(nz,nf,nl)
      real ui(nz,nf,nl),vi(nz,nf,nl)
      real uf(nz,nf,nl),vf(nz,nf,nl)
      real vnq(nz,nf,nl),vnp(nz,nf,nl),vnphi(nz,nf,nl)
      real vhi(nz,nf,nl),nuen(nz,nf,nl), hruti,hrutf,dt,ftc
      
      real :: deniCombined(nz,nf,nl,nion)

      common / var / deni,denn,ne,te,ti,tn,u,v,vsi,vsid,denni,dennf&
           ,sumvsi,vpi,vsic,ftc,dt,vhi,nuen,vnq,vnp,vnphi,&
           tni,tnf,ui,vi,uf,vf,hruti,hrutf
      

!            namelist / go / fmtout,maxstep,hrmax,dthr,hrpr,dt0,
!     .                grad_in,glat_in,glon_in,
!     .                fejer,
!     .                rmin,rmax,
!     .                altmin,
!     .                fbar,f10p7,ap,
!     .                year,day,mmass,
!     .                nion1,nion2,hrinit,tvn0,tvexb0,ver,veh,vw,
!     .                 gams1,gams1m,gamp1,nz1,
!     .                 gams2,gams2m,gamp2,nz2,
!     .                 gams3,gams3m,gamp3,nz3,
!     .                 gams4,gams4m,gamp4,nz4,
!     .                 r_min1,r_max1,
!     .                 r_max2,alt_crit_avg,
!     .                 blat_max3,blat_max4,
!     .                snn,stn,denmin,alt_crit,cqe,plat,plon,
!     .                dellon,psmooth,hall,restart,
!     .                storm_ti,storm_tf,vexb_max,
!     .                lmadala,lcr,lvs,lweimer,decay_time,pcrit,
!     .                lhwm93,lhwm14,
!     .                 vsi0,delta_vsi0,anu_drag0

      ! pass the ion densities 
      if (iProc .ne. 0) then
         !find local portion of dene
         do nntmp = 1,nion
            do itmp = 1,nz
               do jtmp = 1,nf
                  do ktmp = 1,nl
                     denitmp(itmp,jtmp,ktmp) & 
                          = deni(itmp,jtmp,ktmp,nntmp)
                  enddo
               enddo
            enddo
            !send local portion of deni to master
            call mpi_send(denitmp, nz*nf*nl, MPI_REAL, 0, 0, &
                 iComm, ierr)
!            write(*,*) 'Sending from worker', iProc,'nntmp',nntmp
         end do
      else
         if (.not.allocated(denet)) allocate(denet(nz,nf,nlt))
         if (.not.allocated(denit)) allocate(denit(nz,nf,nlt,nion))
         
         WORKER: do iwrk=1,nProc-1
            do nntmp = 1,nion
               call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                    iComm, status, ierr)
!               write(*,*) 'recieving from worker', iwrk,'nntmp',nntmp
               !now combine portions of ion density array
               do itmp = 1,nz
                  do jtmp = 1,nf
                     do ktmp = 1,nl
                        deniCombined(itmp,jtmp,ktmp,nntmp) &
                             =  denitmp(itmp,jtmp,ktmp)
                     enddo
                  enddo
               enddo
            enddo
            
            
            ! Put the submatrices into the correct matrix
            
            do nn = 1,nion
               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +k -1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        denit(i,j,kk,nn) = deniCombined(i,j,k,nn)
                     enddo
                  enddo
               enddo
            end do
         end do WORKER
         ! calculate the electron density from the ion densities
         do k = 1,nlt
            do j = 1,nf
               do i = 1,nz
                  denet(i,j,k) = 0
                  do ni = nion1,nion2
                     denet(i,j,k) = denet(i,j,k) + denit(i,j,k,ni)
                  enddo
               enddo
            enddo
         enddo
         ! H+(ni =1), O+(ni=2), He+ (ni=5) see for instance ptop in param3 file 
         ! get2D plasmasphere density map on SAMI grid
         PlasEq_C(1:nLatSami,1:nLonSami-1) = denet(nz/2,1:nLatSami,1:nLonSami-1)
         PlasHpEq_C(1:nLatSami,1:nLonSami-1) = &
              denit(nz/2,1:nLatSami,1:nLonSami-1,pthp)
         PlasHepEq_C(1:nLatSami,1:nLonSami-1) = &
              denit(nz/2,1:nLatSami,1:nLonSami-1,pthep)
         PlasOpEq_C(1:nLatSami,1:nLonSami-1) = &
              denit(nz/2,1:nLatSami,1:nLonSami-1,ptop)
         !fill ghost cells
         PlasEq_C(1:nLatSami,nLonSami) = denet(nz/2,1:nLatSami,1)
         PlasHpEq_C(1:nLatSami,nLonSami)  = PlasHpEq_C(1:nLatSami,1)
         PlasHepEq_C(1:nLatSami,nLonSami) = PlasHepEq_C(1:nLatSami,1)
         PlasOpEq_C(1:nLatSami,nLonSami)  = PlasOpEq_C(1:nLatSami,1)
         ! deallocate 3D arrays to save memory
         deallocate(denit,denet)
         !call plot_2dplas
      endif
      
      
      
    end subroutine gather_dene
    !===========================================================================
    subroutine plot_coupled_values
      use ModIoUnit,    ONLY: UnitTmp_
      use ModNumConst, ONLY: cDegToRad, cRadToDeg
      integer :: iLat, iLon
      
      open(UnitTmp_,FILE='SamiGridPot.dat')
      write(UnitTmp_,'(a)') &
           'VARIABLES = "Lat", "Lon", "Phi"'
      write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLatSami, &
           ', J=', nLonSami-1,', DATAPACKING=POINT'
      
      do iLon = 1, nLonSami-1
         do iLat = 1, nLatSami
            write(UnitTmp_,"(100es18.10)") bLatSami(iLat,iLon),&
                 bLonSami(iLat,iLon),PotCimiOnSamiGrid_C(iLon,iLat)
         enddo
      enddo
      close(UnitTmp_)
      
      open(UnitTmp_,FILE='CimiGridPot.dat')
      write(UnitTmp_,'(a)') &
           'VARIABLES = "Lat", "Lon", "Phi"'
      write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLatCimi, &
           ', J=', nLonCimi,', DATAPACKING=POINT'
      
      do iLon = 1, nLonCimi
         do iLat = 1, nLatCimi
            write(UnitTmp_,"(100es18.10)") LatCimi_C(iLat)*cRadToDeg,&
                 LonCimi_C(iLon)*cRadToDeg,PotCimi_C(iLat,iLon)
         enddo
      enddo
      close(UnitTmp_)
      
    end subroutine plot_coupled_values


    !===========================================================================
    subroutine plot_2dplas
      use ModIoUnit,    ONLY: UnitTmp_
      use ModNumConst, ONLY: cDegToRad, cRadToDeg
      integer :: iLat, iLon
      
      open(UnitTmp_,FILE='Sami3PlasIn2D.dat')
      write(UnitTmp_,'(a)') &
           'VARIABLES = "Lat", "Lon", "edens[/cc]" "Hpdens[/cc]" "Hepdens[/cc]" "Opdens[/cc]"'
      write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLatSami, &
           ', J=', nLonSami-1,', DATAPACKING=POINT'
      
      do iLon = 1, nLonSami-1
         do iLat = 1, nLatSami
            write(UnitTmp_,"(100es18.10)") bLatSami(iLat,iLon),&
                 bLonSami(iLat,iLon),PlasEq_C(iLat,iLon),PlasHpEq_C(iLat,iLon),&
                 PlasHepEq_C(iLat,iLon),PlasOpEq_C(iLat,iLon)
         enddo
      enddo
      close(UnitTmp_)
      
    end subroutine plot_2dplas


  
end Module ModCoupleCimi
  
