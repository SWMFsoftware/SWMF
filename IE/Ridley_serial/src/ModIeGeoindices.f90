!==============================================================================
module ModIeGeoindices
  
  ! ModIeGeoindices is a module for calculating geomagnetic indices through
  ! the use of virtual ground based magnetometers.  At present, only Kp is
  ! instituted.  
  ! Because indices require the combination of iono- and magneto-spheric
  ! contributions, most of the work performed in this module prepares
  ! values to hand to GM such that the actual index may be calculated on the
  ! GM side of the coupling.  GM also writes the final index file.

  implicit none
  save
  
  logical, private :: IsInitialized=.false.

  ! Specify what indices should be calculated and handed to GM:
  logical :: DoCalcIndices=.false., DoCalcKp=.false., DoCalcDst=.false.

  ! Each index has a set number of mags.  Total number passed to GM
  ! is nIndexMag = sum of all mags used.
  integer           :: nIndexMag=0
  integer, parameter:: nKpMag=24, nDstMag=4 ! Number of mags per index.
  real              :: XyzKp_DI(3, nKpMag)  ! Locations of kp mags, SMG.
  real, parameter   :: fakeplat = 50.0      ! Latitude of Kp stations.

contains

  !===========================================================================
  subroutine init_geoindices
    ! Initialize all arrays, etc.
    use ModNumConst,    ONLY: cDegToRad, cTwoPi
    use ModPlanetConst, ONLY: rPlanet_I, Earth_

    integer :: i
    real    :: radXY, phi

    character(len=*), parameter :: NameSub='init_geoindices'
    logical :: DoTest, DoTestMe   
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Tally up total number of index-related magnetometers.
    if(DoCalcKp ) nIndexMag=nIndexMag+nKpMag
    if(DoCalcDst) nIndexMag=nIndexMag+nDstMag

    if(DoTest)write(*,*)'IE: '//NameSub//' Total number of mags = ', nIndexMag

    ! Initialize grid and arrays.  FaKe_p uses stations fixed in SMG coords.
    XyzKp_DI(3,:) = rPlanet_I(Earth_)*sin(fakepLat * cDegToRad)
    radXY         = rPlanet_I(Earth_)*cos(fakepLat * cDegToRad)
    do i=1, nKpMag
       phi = cTwoPi * (i-1)/24.0
       XyzKp_DI(1,i) = radXY * cos(phi)
       XyzKp_DI(2,i) = radXY * sin(phi)
       !write(*,'(a, 3(1x, e13.3))')'Coords = ', xyzkp_di(:,1)
    end do

    IsInitialized=.true.

  end subroutine init_geoindices

  !===========================================================================
  subroutine get_index_mags(MagOut_DI)
    ! Obtain index-related magnetometer data to pass to BATS-R-US.
    use ModIonoMagPerturb, ONLY: iono_mag_perturb
    use ModProcIE,         ONLY: iProc, nProc, iComm
    use ModMpi

    real, intent(out) :: MagOut_DI(3,nIndexMag) ! Mag B in NED coordinates.

    integer:: i, nTmpMag, iError
    real, dimension(3,nIndexMag) :: &
         XyzSmg_DI, MagPerturb_Jh_DI, MagPerturb_Jp_DI, &
         MagPertTotal_DI, MagSum_dI

    character(len=*), parameter :: NameSub='get_index_mags'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Fill XyzSmg depending on what indices we are calculating.
    i = 1
    if(DoCalcKp) then
       XyzSmg_DI(:,i:i+nKpMag-1) = XyzKp_DI
       i = i + nKpMag
    end if
    if(DoCalcDst) then
       !XyzSmzg_DI(:,i:i+nKpMag-1) = XyzDst_DI ! NOT IMPLEMENTED YET.
       i = i + nDstMag
    end if

    ! Check to ensure indexing of array was correct.
    if( nIndexMag .ne. (i-1) ) call CON_stop(&
         NameSub//' Indexing error!  Not all magnetometers accounted for.')

    ! Collect Hall and Peterson B pertubations, sum them.
    call iono_mag_perturb(nIndexMag, XyzSmg_DI, MagPerturb_Jh_DI, MagPerturb_Jp_DI)
    MagPertTotal_DI = MagPerturb_Jh_DI + MagPerturb_Jp_DI

    ! Collect the variables from all the PEs
    MagOut_DI=0.0
    MagSum_DI=0.0
    if(nProc>1)then 
       call MPI_reduce(MagPertTotal_DI, MagSum_DI, 3*nIndexMag, &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       if(iProc==0) MagOut_DI = MagSum_DI
       call MPI_bcast(MagOut_DI, 3*nIndexMag, MPI_REAL,0,iComm,iError)
    else
       MagOut_DI=MagPertTotal_DI
    end if

  end subroutine get_index_mags

  !===========================================================================

end module ModIeGeoindices
!==============================================================================
