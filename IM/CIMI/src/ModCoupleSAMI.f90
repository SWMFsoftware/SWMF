Module ModCoupleSami
    implicit none 
    private !except

    logical, public :: DoCoupleSami = .false.
    
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal

    ! SAMI Latitude (radians) and lon (radians) grid
    real, allocatable :: LatSami_C(:),LonSami_C(:)
    integer :: nLatSAMI, nLonSAMI
    
    ! sami plasmasphere on sami grid
    real, allocatable :: PlasSAMI_C(:,:),PlasHpSAMI_C(:,:),PlasHepSAMI_C(:,:),&
         PlasOpSAMI_C(:,:)
    
    ! sami plasmasphere on cimi grid
    real, public, allocatable :: PlasSamiOnCimiGrid_C(:,:)
    real, public, allocatable :: PlasHpSamiOnCimiGrid_C(:,:)
    real, public, allocatable :: PlasHepSamiOnCimiGrid_C(:,:)
    real, public, allocatable :: PlasOpSamiOnCimiGrid_C(:,:)

    !public methods
    public :: cimi_set_global_mpi
    public :: cimi_get_init_for_sami
    public :: cimi_send_to_sami
    public :: cimi_put_init_from_sami
    public :: cimi_get_from_sami
  contains
    !==========================================================================
    subroutine cimi_set_global_mpi(iProc0CimiIN,iProc0SamiIN, iCommGlobalIN)
      integer, intent(in) ::iProc0CimiIN,iProc0SamiIN, iCommGlobalIN
      !------------------------------------------------------------------------
      iProc0CIMI=iProc0CimiIN
      iProc0Sami=iProc0SamiIN
      iCommGlobal=iCommGlobalIN

      ! Set the coupling to be true
      DoCoupleSami=.true.
    end subroutine cimi_set_global_mpi
    !==========================================================================
    ! only call by proc 0
    subroutine cimi_get_init_for_sami
      use ModCimiGrid, ONLY: xlatr, nLat=>np, phi, nLon=>nt
      use ModImTime,   ONLY: iStartTime_I
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      ! Send the starttime from cimi to sami
      call MPI_send(iStartTime_I,7,MPI_INTEGER,iProc0SAMI,&
           1,iCommGlobal,iError)

      ! Send the grid parameters from cimi to sami
      call MPI_send(nLat,1,MPI_INTEGER,iProc0SAMI,&
           2,iCommGlobal,iError)
      
      call MPI_send(nLon,1,MPI_INTEGER,iProc0SAMI,&
           3,iCommGlobal,iError)

      call MPI_send(xlatr,nLat,MPI_REAL,iProc0SAMI,&
           4,iCommGlobal,iError)

      call MPI_send(phi,nLon,MPI_REAL,iProc0SAMI,&
           5,iCommGlobal,iError)
      


    end subroutine cimi_get_init_for_sami
    !==========================================================================
    ! 
    subroutine cimi_send_to_sami
      use ModCimiGrid, ONLY: nLat=>np, phi, nLon=>nt, iProc
      use ModIeCimi,   ONLY: pot
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      if(iProc==0) write(*,*) 'cimi_send_to_sami',nLat,nLon
      ! Send the starttime from cimi to sami
      if(iProc ==0) &
           call MPI_send(pot,nLat*nLon,MPI_REAL,iProc0SAMI,1,iCommGlobal,iError)
      if(iProc==0) write(*,*) 'Max and Min of pot: ', MAXVAL(pot), MINVAL(pot)
      if(iProc==0) write(*,*) 'cimi_send_to_sami done'
    end subroutine cimi_send_to_sami
    
    !========================================================================
    subroutine cimi_put_init_from_sami
      use ModCimiGrid, ONLY: iProc,iComm, nLat=>np, nLon=>nt
      use ModMpi
      use ModNumConst,    ONLY: cDegToRad
      integer :: iStatus_I(MPI_STATUS_SIZE)
      integer :: iError
      !----------------------------------------------------------------------

      ! recieve SAMI grid size
      if (iProc == 0) call MPI_recv(nLatSAMI,1,MPI_INTEGER,iProc0SAMI,&
           1,iCommGlobal,iStatus_I,iError)


      if (iProc == 0) call MPI_recv(nLonSAMI,1,MPI_INTEGER,iProc0SAMI,&
           2,iCommGlobal,iStatus_I,iError)

      ! bcast grid size
      call MPI_bcast(nLatSAMI,1,MPI_LOGICAL,0,iComm,iError)
      call MPI_bcast(nLonSAMI,1,MPI_LOGICAL,0,iComm,iError)

      ! allocate SAMI grid
      allocate(LatSami_C(nLatSAMI))
      allocate(LonSami_C(nLonSAMI))

      if (iProc == 0) then 
         call MPI_recv(LatSami_C,nLatSAMI,MPI_REAL,iProc0SAMI,&
              3,iCommGlobal,iStatus_I,iError)
         call MPI_recv(LonSami_C(1:nLonSAMI-1),nLonSAMI-1,MPI_REAL,iProc0SAMI,&
              4,iCommGlobal,iStatus_I,iError)
         ! last lon point in sami is 0, but make it equivilently 360 for 
         ! interpolation purposes
         LonSami_C(nLonSAMI) = 360.0

         ! convert sami lat and lon into radians from degrees
         LatSami_C = LatSami_C*cDegToRad
         LonSami_C = LonSami_C*cDegToRad
      endif

      ! allocate arrays that hold data from sami
      allocate(PlasSAMI_C(nLatSAMI,nLonSAMI))
      allocate(PlasHpSAMI_C(nLatSAMI,nLonSAMI))
      allocate(PlasHepSAMI_C(nLatSAMI,nLonSAMI))
      allocate(PlasOpSAMI_C(nLatSAMI,nLonSAMI))
      allocate(PlasSamiOnCimiGrid_C(nLat,nLon))      
      allocate(PlasHpSamiOnCimiGrid_C(nLat,nLon))      
      allocate(PlasHepSamiOnCimiGrid_C(nLat,nLon))      
      allocate(PlasOpSamiOnCimiGrid_C(nLat,nLon))      
    end subroutine cimi_put_init_from_sami
    
    !==========================================================================
    ! 
    subroutine cimi_get_from_sami(TimeSimulation)
      use ModMPI
      use ModCimiGrid, ONLY: iProc, nProc, iComm, nLat=>np, nLon=>nt      
      use DensityTemp,    ONLY: density
      
      real, intent(in) :: TimeSimulation
      integer :: iwrk, iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      real, parameter :: cCm3ToM3 = 1.0e6
      !------------------------------------------------------------------------
      write(*,*) 'Starting cimi_get_from_sami on iProc=', iProc
      !if(iProc ==0) write(*,*) 'sami_get_from_cimi',nLatCIMI,nLonCIMI
      ! Send the starttime from cimi to sami
      if(iProc ==0) then
         call MPI_recv(PlasSAMI_C,nLatSAMI*nLonSAMI,MPI_REAL,&
              iProc0SAMI,1,iCommGlobal,iStatus_I,iError)
         call MPI_recv(PlasHpSAMI_C,nLatSAMI*nLonSAMI,MPI_REAL,&
              iProc0SAMI,1,iCommGlobal,iStatus_I,iError)
         call MPI_recv(PlasHepSAMI_C,nLatSAMI*nLonSAMI,MPI_REAL,&
              iProc0SAMI,1,iCommGlobal,iStatus_I,iError)
         call MPI_recv(PlasOpSAMI_C,nLatSAMI*nLonSAMI,MPI_REAL,&
              iProc0SAMI,1,iCommGlobal,iStatus_I,iError)
      endif
      ! write out max and min of plas
      if(iProc ==0) then
         write(*,*) 'Max and min of SAMI Plas recv in CIMI: ', &
              maxval(PlasSAMI_C), minval(PlasSAMI_C)
         write(*,*) 'Max and min of SAMI H+ Plas recv in CIMI: ', &
              maxval(PlasHpSAMI_C), minval(PlasHpSAMI_C)
         write(*,*) 'Max and min of SAMI He+ Plas recv in CIMI: ', &
              maxval(PlasHepSAMI_C), minval(PlasHepSAMI_C)
         write(*,*) 'Max and min of SAMI O+ Plas recv in CIMI: ', &
              maxval(PlasOpSAMI_C), minval(PlasOpSAMI_C)
      endif
      ! interpolate sami plas to cimi grid
      if(iProc ==0) call interpolate_sami_to_cimi(TimeSimulation)
      
      if(iProc ==0) then
         write(*,*) 'Max and min of SAMI Plas on cimi Grid: ',&
              maxval(PlasSamiOnCimiGrid_C), minval(PlasSamiOnCimiGrid_C)
         write(*,*) 'Max and min of SAMI H+ Plas on cimi Grid: ',&
              maxval(PlasHpSamiOnCimiGrid_C), minval(PlasHpSamiOnCimiGrid_C)
         write(*,*) 'Max and min of SAMI He+ Plas on cimi Grid: ',&
              maxval(PlasHepSamiOnCimiGrid_C), minval(PlasHepSamiOnCimiGrid_C)
         write(*,*) 'Max and min of SAMI O+ Plas on cimi Grid: ',&
              maxval(PlasOpSamiOnCimiGrid_C), minval(PlasOpSamiOnCimiGrid_C)
      endif
!      write(*,*) iProc,'started send/recv of PlasSamiOnCimiGrid_C'
!      call MPI_BARRIER(iComm,iError)
!     if(iProc==0) then
!         do iwrk=1,nProc
!            call MPI_send(PlasSamiOnCimiGrid_C,nLat*nLon,MPI_REAL,iwrk,&
!                 1,iComm,iError)
!         enddo
!      else
!         call MPI_recv(PlasSamiOnCimiGrid_C,nLat*nLon,MPI_REAL,&
!              0,1,iComm,iStatus_I,iError)
!      endif
!
!      ! now bcast sami plas on cimi grid to all procs
      call MPI_bcast(PlasSamiOnCimiGrid_C,nLat*nLon, &
           MPI_REAL,0,iComm,iError)
      call MPI_bcast(PlasHpSamiOnCimiGrid_C,nLat*nLon, &
           MPI_REAL,0,iComm,iError)
      call MPI_bcast(PlasHepSamiOnCimiGrid_C,nLat*nLon, &
           MPI_REAL,0,iComm,iError)
      call MPI_bcast(PlasOpSamiOnCimiGrid_C,nLat*nLon, &
           MPI_REAL,0,iComm,iError)
      
      write(*,*) iProc,'finished send/recv of PlasSamiOnCimiGrid_C'

      ! Set plas in module for waves
      density=PlasSamiOnCimiGrid_C*cCm3ToM3

    end subroutine cimi_get_from_sami

    !==========================================================================
    subroutine interpolate_sami_to_cimi(TimeSimulation)
      use ModCimiGrid, ONLY: nLat=>np, nLon=>nt,Lat_C=>xlatr,Lon_C=>phi      
      use ModInterpolate, ONLY: bilinear
      use ModNumConst,    ONLY: cDegToRad,cPi,cTwoPi
      
      real, intent(in) :: TimeSimulation
      integer :: iLat, iLon
      real :: LatLon_D(2)
      real, parameter :: DensityMin=0.1 !/cc

      do iLon = 1, nLon
         do iLat = 1, nLat
            !write(*,*) 'iLat,iLon,cDegToRad',iLat,iLon,cDegToRad
            !LatLon_D(1) = bLatSami(iLat,iLon) * cDegToRad
            !LatLon_D(2) = bLonSami(iLat,iLon) * cDegToRad
            
            !convert point on CIMI grid to SAMI coords for interpolation 
            !call convert_cimi_to_sami_lat_lon(TimeSimulation, &
            !     Lat_C(iLat),Lon_C(iLon), &
            !     LatLon_D(1), LatLon_D(2))
            LatLon_D(1) = Lat_C(iLat)
            LatLon_D(2) = modulo(Lon_C(iLon)+cPi,cTwoPi)

            if (LatLon_D(1) > maxval(LatSami_C) .or. &
                 LatLon_D(1) < minval(LatSami_C) ) then
!!$               PlasSamiOnCimiGrid_C(iLat,iLon)    = DensityMin
!!$               PlasHpSamiOnCimiGrid_C(iLat,iLon)  = 1E-6
!!$               PlasHepSamiOnCimiGrid_C(iLat,iLon) = 1E-7
!!$               PlasOpSamiOnCimiGrid_C(iLat,iLon)  = 1E-6
               PlasSamiOnCimiGrid_C(iLat,iLon)    = 0.0
               PlasHpSamiOnCimiGrid_C(iLat,iLon)  = 0.0
               PlasHepSamiOnCimiGrid_C(iLat,iLon) = 0.0
               PlasOpSamiOnCimiGrid_C(iLat,iLon)  = 0.0
            else
!!$               PlasSamiOnCimiGrid_C(iLat,iLon) = &
!!$                    max(bilinear(PlasSAMI_C,1,nLatSAMI,1,nLonSAMI,LatLon_D, &
!!$                    LatSami_C,LonSami_C,DoExtrapolate=.true.),DensityMin)
               PlasSamiOnCimiGrid_C(iLat,iLon) = &
                    bilinear(PlasSAMI_C,1,nLatSAMI,1,nLonSAMI,LatLon_D, &
                    LatSami_C,LonSami_C,DoExtrapolate=.true.)
               PlasHpSamiOnCimiGrid_C(iLat,iLon) = &
                    bilinear(PlasHpSAMI_C,1,nLatSAMI,1,nLonSAMI,LatLon_D, &
                    LatSami_C,LonSami_C,DoExtrapolate=.true.)
               PlasHepSamiOnCimiGrid_C(iLat,iLon) = &
                    bilinear(PlasHepSAMI_C,1,nLatSAMI,1,nLonSAMI,LatLon_D, &
                    LatSami_C,LonSami_C,DoExtrapolate=.true.)
               PlasOpSamiOnCimiGrid_C(iLat,iLon) = &
                    bilinear(PlasOpSAMI_C,1,nLatSAMI,1,nLonSAMI,LatLon_D, &
                    LatSami_C,LonSami_C,DoExtrapolate=.true.)
            endif
         enddo
      enddo
      
      !      ! add the ghost cell for potential
!      PotCimiOnSamiGrid_C(:,nLonSami) = PotCimiOnSamiGrid_C(:,1) 
      !call plot_coupled_values
      

    end subroutine interpolate_sami_to_cimi
    !==========================================================================
    ! subroutine that takes CIMI lat lon and returns SAMI lat lon
    subroutine convert_cimi_to_sami_lat_lon(TimeSimulation, CimiLat,CimiLon, &
         SamiLat, SamiLon)
      use CON_axes, ONLY: transform_matrix
      use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg
      
      real, intent(in) :: TimeSimulation, CimiLat, CimiLon
      real, intent(out):: SamiLat, SamiLon
      real             :: CimiToSami_DD(3,3), theta, phi 
      real             :: XyzCimi_D(3), XyzSami_D(3)
      character(len=*),  parameter :: NameCimiCoord='SMG'
      character(len=*),  parameter :: NameSamiCoord='MAG'
      !------------------------------------------------------------------------
      
      ! Get polar angle (theta) and azimuthal angle (phi)
      theta = 0.5*cPi - CimiLat
      phi   = CimiLon
      
      ! Get xyzCimi_D from CimiLat and CimiLon
      xyzCimi_D(1) = sin(theta)*cos(phi)
      xyzCimi_D(2) = sin(theta)*sin(phi)
      xyzCimi_D(3) = cos(theta)
      
      !\
      ! get equivalent geographic coords Hemisphere 1
      !/
      
      ! Get transform matrix 
      CimiToSami_DD = &
           transform_matrix(TimeSimulation, NameCimiCoord, NameSamiCoord)
      
      ! Transform xyzCimi_D to XyzSami_DI
      XyzSami_D = matmul( CimiToSami_DD, XyzCimi_D)
      
      ! Calculate SamiLat and SamiLon 
      SamiLon = modulo(atan2(XyzSami_D(2), XyzSami_D(1)), cTwoPi) * cRadToDeg
      SamiLat = 90.0 - (acos(max(-1.0,min(1.0, XyzSami_D(3))))*cRadToDeg)
      
      
    end subroutine convert_cimi_to_sami_lat_lon
    !===========================================================================
    subroutine plot_coupled_values
      use ModCimiGrid, ONLY: nLat=>np, nLon=>nt,Lat_C=>xlatr,Lon_C=>phi      
      use ModIoUnit,    ONLY: UnitTmp_
      use ModNumConst, ONLY: cDegToRad, cRadToDeg
      integer :: iLat, iLon
      
      open(UnitTmp_,FILE='SamiGridPlas.dat')
      write(UnitTmp_,'(a)') &
           'VARIABLES = "Lat", "Lon", "ePlas", "HpPlas", "HepPlas", "OpPlas"'
      write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLatSami, &
           ', J=', nLonSami-1,', DATAPACKING=POINT'
      
      do iLon = 1, nLonSami-1
         do iLat = 1, nLatSami
            write(UnitTmp_,"(100es18.10)") LatSami_C(iLat),&
                 LonSami_C(iLon),PlasSAMI_C(iLat,iLon),PlasHpSAMI_C(iLat,iLon),&
                 PlasHepSAMI_C(iLat,iLon),PlasOpSAMI_C(iLat,iLon)
         enddo
      enddo
      close(UnitTmp_)
      
      open(UnitTmp_,FILE='CimiGridPlas.dat')
      write(UnitTmp_,'(a)') &
           'VARIABLES = "Lat", "Lon", "ePlas", "HpPlas", "HepPlas", "OpPlas"'
      write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLat, &
           ', J=', nLon,', DATAPACKING=POINT'
      
      do iLon = 1, nLon
         do iLat = 1, nLat
            write(UnitTmp_,"(100es18.10)") Lat_C(iLat),&
                 Lon_C(iLon),PlasSamiOnCimiGrid_C(iLat,iLon),&
                 PlasHpSamiOnCimiGrid_C(iLat,iLon),&
                 PlasHepSamiOnCimiGrid_C(iLat,iLon),&
                 PlasOpSamiOnCimiGrid_C(iLat,iLon)
         enddo
      enddo
      close(UnitTmp_)
      
    end subroutine plot_coupled_values

    
  end Module ModCoupleSami
  
