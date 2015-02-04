!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
  
  logical :: IsInitialized=.false.

  ! Specify what indices should be calculated and handed to GM:
  logical :: DoCalcIndices=.false., DoCalcKp=.false., DoCalcDst=.false.

  ! Each index has a set number of mags.  Total number passed to GM
  ! is nIndexMag = sum of all mags used.
  integer           :: nIndexMag=0
  integer, parameter:: nKpMag=24, nDstMag=4 ! Number of mags per index.
  real              :: XyzKp_DI(3, nKpMag)  ! Locations of kp mags, SMG.
  real, parameter   :: FakePLat = 60.0      ! Latitude of Kp stations.

contains

  !===========================================================================
  subroutine init_geoindices
    ! Initialize all arrays, etc.
    use ModNumConst,    ONLY: cDegToRad, cTwoPi
    use ModPlanetConst, ONLY: rPlanet_I, Earth_

    integer :: i
    real    :: RadXY, Phi

    character(len=*), parameter :: NameSub='init_geoindices'
    logical :: DoTest, DoTestMe   
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Tally up total number of index-related magnetometers.
    nIndexMag=0
    if(DoCalcKp ) nIndexMag=nIndexMag+nKpMag
    if(DoCalcDst) nIndexMag=nIndexMag+nDstMag

    if(DoTest)write(*,*)'IE: '//NameSub//' Total number of mags = ', nIndexMag

    ! Initialize grid and arrays.  FaKe_p uses stations fixed in SMG coords.
    XyzKp_DI(3,:) = rPlanet_I(Earth_)*sin(FakePLat * cDegToRad)
    RadXY         = rPlanet_I(Earth_)*cos(FakePLat * cDegToRad)
    do i=1, nKpMag
       Phi = cTwoPi * (i-1)/24.0
       XyzKp_DI(1,i) = RadXY * cos(Phi)
       XyzKp_DI(2,i) = RadXY * sin(Phi)
       !write(*,'(a, 3(1x, e13.3))')'Coords = ', xyzkp_di(:,1)
    end do

    IsInitialized=.true.

  end subroutine init_geoindices

  !===========================================================================
  subroutine get_index_mags(MagJhOut_DI, MagJpOut_DI)
    ! Obtain index-related magnetometer perturbations from Hall (MagJhOut_DI)
    ! Pederson (MagJpOut_DI) currents to pass to the GM module.

    use ModIonoMagPerturb, ONLY: iono_mag_perturb
    use ModProcIE,         ONLY: nProc, iComm
    use ModMpi

    ! Mag B in NED coordinates.
    real, intent(out), dimension(3,nIndexMag) :: MagJhOut_DI, MagJpOut_DI 

    integer:: i, iError
    real, dimension(3,nIndexMag) :: &
         XyzSmg_DI, MagPerturbJh_DI, MagPerturbJp_DI, &
         MagPertTotalJh_DI, MagPertTotalJp_DI


    character(len=*), parameter :: NameSub='get_index_mags'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe) write(*,*)'IE: '//NameSub//' Collecting geomagnetic indices.'
    if(.not. IsInitialized) call CON_stop(NameSub// &
         ': Error: Geomagnetic indices not initialized!')

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

    ! Collect Hall and Pedersen B perturbations.
    call iono_mag_perturb(nIndexMag, XyzSmg_DI, MagPerturbJh_DI, MagPerturbJp_DI)

    ! Collect the variables from all the PEs
    if(nProc>1)then 
       call MPI_allreduce(MagPertTotalJh_DI, MagJhOut_DI, 3*nIndexMag, &
            MPI_REAL, MPI_SUM, iComm, iError)
       call MPI_allreduce(MagPertTotalJp_DI, MagJpOut_DI, 3*nIndexMag, &
            MPI_REAL, MPI_SUM, iComm, iError)
    else
       MagJhOut_DI = MagPertTotalJh_DI
       MagJpOut_DI = MagPertTotalJp_DI
    end if

  end subroutine get_index_mags

  !===========================================================================

end module ModIeGeoindices
!==============================================================================
