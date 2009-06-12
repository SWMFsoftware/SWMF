! Wrapper for Internal Magnetosphere (IM) component
!=============================================================================

subroutine IM_set_param(CompInfo,TypeAction)

  use CON_comp_info
  use ModProcIM
  use ModHeidiMain
  use ModReadParam, only: i_session_read
  use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case
  use ModHeidiIO, ONLY : IsFramework, StringPrefix
  use ModIoUnit, only: STDOUT_

  implicit none
  character (len=*), parameter :: NameSub='IM_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  !LOCAL VARIABLES:
  character (len=100) :: NameCommand, StringPlot
  logical             :: DoEcho=.false.
  logical             :: UseStrict=.true.  
  integer             :: iUnitOut
  !---------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                         &
          Use=.true.,                           &
          NameVersion='HEIDI (Liemohn)', &
          Version=1.1)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     if(nProc>4)call CON_stop( NameSub // &
          ' IM_ERROR this version can run on 4 PE !')

     IsFramework = .true.
  case('READ')
     call heidi_read
  case('CHECK')
     !We should check and correct parameters here
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
     call heidi_check
  case('STDOUT')
     iUnitOut = STDOUT_
     if(nProc==1)then
        StringPrefix='IM:'
     else
        write(StringPrefix,'(a,i3.3,a)')'IM',iProc,':'
     end if
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call IM_set_grid
  case default
     call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine IM_set_param

!============================================================================

subroutine IM_set_grid

  use ModIonoHeidi
  use ModHeidiSize
  use ModHeidiMain
  use ModProcIM
  use ModNumConst
  use CON_coupler, ONLY: set_grid_descriptor, is_proc, IM_

  implicit none
  character (len=*), parameter :: NameSub='IM_set_grid'
  real :: Radius_I(1)
  real :: Colat_I(2*IONO_nTheta-1)
  logical :: IsInitialized = .false.
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)'IM_set_grid_descriptor called, IsInitialized=',&
       IsInitialized
  if(IsInitialized) RETURN
  IsInitialized=.true.

  Radius_I(1) = IONO_Radius ! radial size of the ionosphere in meters

  ! Total hack, because I don't actually know the real grid we will use in the 
  ! actual code.

  ! The colatitudes for both hemispheres
  Colat_I(            1:  IONO_nTheta) = IONO_NORTH_Theta(:,1)
  Colat_I(IONO_nTheta:2*IONO_nTheta-1) = IONO_SOUTH_Theta(:,1)

  call set_grid_descriptor(                        &
       IM_,                          &! component index
       nDim=2,                       &! dimensionality
       nRootBlock_D=(/1,1/),         &! north+south hemispheres
       nCell_D =(/2*IONO_nTheta-1,IONO_nPsi/), &! size of node based grid
       XyzMin_D=(/cOne, cOne/),      &! min colat and longitude indexes
       XyzMax_D=(/real(2*IONO_nTheta-1),&
       real(IONO_nPsi)/),            &! max colat and longitude indexes
       TypeCoord='SMG',                            &! solar magnetic coord.
       Coord1_I=Colat_I,             &! colatitudes
       Coord2_I=IONO_NORTH_Psi(1,:),               &! longitudes
       Coord3_I=(/Radius_I/),                          &! radial size in meters
       IsPeriodic_D=(/.false.,.true./))

!!!  ! IM grid size in generalized coordinates
!!!  call set_grid_descriptor( IM_,                         & ! component index
!!!       nDim=4,                                           & ! dimensionality
!!!       nRootBlock_D=(/1,1/),                             & ! single block
!!!       nCell_D=(/nT,nR/),                          & ! size of cell based grid
!!!       XyzMin_D=(/cHalf, cHalf/),                        & ! min gen.coords for cells
!!!       XyzMax_D=(/nT+cHalf,nR+cHalf/),             & ! max gen.coords for cells
!!!       TypeCoord='SMG',                                   & ! solar magnetic coord
!!!       Coord1_I=phi,  & ! magnetic local times
!!!       Coord2_I=z,    & ! l-shell?
!!!!       Coord1_I=real(colat(1:NT,1)),                     & ! colatitudes
!!!!       Coord2_I=real(aloct(1,1:NR)),                     & ! longitudes
!!!!       Coord3_I=Radius_I,                                & ! radial size in meters
!!!!       Coord4_I=NPA(0:90)                                & ! Grid in pitch angle
!!!!       COORD5_I=NE                                       & ! Grid in energy
!!!       IsPeriodic_D=(/.true.,.false./))    ! periodic in longitude
!!!!       IsPeriodic_D=(/.true.,.false.,.false.,.false./))    ! periodic in longitude

end subroutine IM_set_grid
!==============================================================================
subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide current for IE
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V(??)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModIonoHeidi, ONLY: IONO_NORTH_RCM_JR,IONO_SOUTH_RCM_JR, IONO_nTheta, IONO_nPsi

  implicit none
  character(len=*), parameter :: NameSub='IM_get_for_ie'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iLat, iLon, iBlock, iPoint
  real    :: w

  !---------------------------------------------------------------------------
  Buff_V = 0.0

  do iPoint = iPointStart, iPointStart + nPoint - 1

     iLat   = Index % iCB_II(1,iPoint)
     iLon   = Index % iCB_II(2,iPoint)
     iBlock = Index % iCB_II(3,iPoint)
     w      = Weight % Weight_I(iPoint)

     if(iBlock/=1)then
        write(*,*)NameSub,': iPoint,Index % iCB_II=',&
             iPoint,Index%iCB_II(:,iPoint)
        call CON_stop(NameSub//&
             ' SWMF_ERROR iBlock should be 1=North in IM-IE coupling')
     end if

     if(iLat<1 .or. iLat>IONO_nTheta*2 .or. iLon<1 .or. iLon>IONO_nPsi+1)then
        write(*,*)'iLat,iLon=',iLat, IONO_nTheta*2, iLon, IONO_nPsi
        call CON_stop(NameSub//' SWMF_ERROR index out of range')
     end if

     ! Only worry about the northern hemisphere....  IE can fix the southern hemisphere.
     if (iLat <= IONO_nTheta .and. iLon <= IONO_nPsi) &
          Buff_V(1) = Buff_V(1) + w * IONO_NORTH_RCM_JR(iLat,iLon)

     if (iLat > IONO_nTheta .and. iLon <= IONO_nPsi) &
          Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_RCM_JR(2*IONO_nTheta-iLat+1,iLon)

  end do

end subroutine IM_get_for_ie

!============================================================================
subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)

  use ModPlotFile, ONLY: save_plot_file
  implicit none
  integer, intent(in):: nTheta, nPhi
  real,    intent(in):: Potential_II(nTheta, nPhi)

  !call save_plot_file(...)

end subroutine IM_put_from_ie_mpi

!==============================================================================
subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModIonoHeidi, ONLY: IONO_NORTH_PHI, IONO_SOUTH_PHI, IONO_nTheta, IONO_nPsi

  implicit none
  character(len=*), parameter   :: NameSub='IM_put_from_ie'
  integer,intent(in)            :: nPoint, iPointStart, nVar
  real, intent(in)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight
  logical,intent(in)            :: DoAdd
  integer :: iBlock,i,j
  !---------------------------------------------------------------------------
  if(nPoint>1)then
     write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
          nPoint,iPointStart,Weight % Weight_I
     call CON_stop(NameSub//': should be called with 1 point')
  end if
  if(DoAdd)then
     write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
          nPoint,iPointStart,Weight % Weight_I
     write(*,*)NameSub,': WARNING DoAdd is true'
  end if

  i = Index % iCB_II(1,iPointStart)
  j = Index % iCB_II(2,iPointStart)

  if(i<1.or.i>2*IONO_nTheta-1.or.j<1.or.j>IONO_nPsi+1)then
     write(*,*)'i,j,DoAdd=',i,2*IONO_nTheta-1,j,IONO_nPsi+1,DoAdd
     call CON_stop('IM_put_from_ie (in IM_wrapper): index out of range')
  end if

  if (i <= IONO_nTheta .and. j <= IONO_nPsi) then
     if(DoAdd)then
        IONO_NORTH_PHI(i,j)        = IONO_NORTH_PHI(i,j)        + Buff_V(1)
     else
        IONO_NORTH_PHI(i,j)        = Buff_V(1)
     end if
  endif

  if (i > IONO_nTheta .and. j <= IONO_nPsi) then
     if(DoAdd)then
        IONO_SOUTH_PHI(i-IONO_nTheta,j) = &
             IONO_SOUTH_PHI(i-IONO_nTheta,j) + Buff_V(1)
     else
        IONO_SOUTH_PHI(i-IONO_nTheta,j) = Buff_V(1)
     end if
  endif

end subroutine IM_put_from_ie
!==============================================================================
subroutine IM_put_from_ie_complete

  write(*,*) "Don't know what this is really supposed to do.  I think that it is"
  write(*,*) "Supposed to be applying periodic boundaries...?"

end subroutine IM_put_from_ie_complete

!==============================================================================

subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

  ! This should be similar to RBE coupling

  use ModIonoHeidi
  use ModConst

  implicit none

  character (len=*),parameter :: NameSub='IM_put_from_gm'

  integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
  real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
  character (len=*),intent(in)       :: NameVar

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' starting with NameVar=',NameVar

  IonoGmVolume   = Buffer_IIV(:,:,vol_)
  IonoGmXPoint   = Buffer_IIV(:,:,z0x_)
  IonoGmYPoint   = Buffer_IIV(:,:,z0y_)
  IonoGmBField   = Buffer_IIV(:,:,bmin_)
  ! I think that this is mass density in SI units.  Change to number density
  ! in #/cc.  Then get rid of -1 values.
  IonoGmDensity  = Buffer_IIV(:,:,rho_)/cProtonMass/1.0e6
  where (IonoGmDensity < 0.0) IonoGmDensity = 0.0

  ! This is in Pascals
  IonoGmPressure = Buffer_IIV(:,:,p_)
  where (IonoGmPressure < 0.0) IonoGmPressure = 0.0

  IonoGmTemperature = 0.0
  where (IonoGmDensity > 0) &
       IonoGmTemperature = IonoGmPressure/(IonoGmDensity*1.0e6*cBoltzmann)/&
       11604.0 ! k -> eV

  !  write(*,*) 'This is not working'

end subroutine IM_put_from_gm

!==============================================================================

subroutine IM_put_from_gm_line(BufferLine_VI, nVarLineIn, nPointLineIn, NameVar)

  implicit none
  character (len=*),parameter :: NameSub='IM_put_from_gm_line'

  integer, intent(in) :: nVarLineIn, nPointLineIn
  real,    intent(in) :: BufferLine_VI(nVarLineIn,nPointLineIn)
  character(len=*), intent(in) :: NameVar

  ! This should be used eventually

end subroutine IM_put_from_gm_line

!==============================================================================

subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
  ! Puts satellite locations and names from GM into IM variables.
!!!DTW 2007

  use ModHeidiSatellites
  use ModNumConst,   ONLY: cDegToRad

  implicit none
  character (len=*),parameter :: NameSub='IM_put_sat_from_gm'

  ! Arguments
  integer, intent(in)            :: nSats
  real, intent(in)               :: Buffer_III(3,2,nSats)
  character(len=100), intent(in) :: Buffer_I(nSats)

  ! Internal variables
  integer :: iError, iSat, l1, l2

  DoWriteSats = .true.
  nImSats = nSats

  if (nImSats > nMaxSatellites) then
     write(*,*) "nImSats > nMaxSatellites"
     call CON_stop("Stoping in routine " // NameSub)
  endif

  ! Assign incoming values, remove path and extension from name.
  SatLoc_3I = Buffer_III
  do iSat=1, nSats
     l1 = index(Buffer_I(iSat), '/', back=.true.) + 1
     l2 = index(Buffer_I(iSat), '.') - 1
     if (l1-1<=0) l1=1
     if (l2+1<=0) l2=len_trim(Buffer_I(iSat))
     NameSat_I(iSat) = Buffer_I(iSat)(l1:l2)
  end do

  ! Change to correct units (degrees to radians)
  SatLoc_3I(1,2,:) = (90. - SatLoc_3I(1,2,:)) * cDegToRad
  SatLoc_3I(2,2,:) =        SatLoc_3I(2,2,:)  * cDegToRad

end subroutine IM_put_sat_from_gm

!==============================================================================

subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  use CON_time, ONLY : get_time
  use ModNumConst, ONLY: cPi, cDegToRad
  use ModConst, ONLY: cProtonMass
  use ModIonoHeidi
  use ModHeidiSize
  use ModHeidiCurrents
  implicit none
  character (len=*),parameter :: NameSub='IM_get_for_gm'

  integer, intent(in)                                :: iSizeIn,jSizeIn,nVar
  real, dimension(iSizeIn,jSizeIn,nVar), intent(out) :: Buffer_IIV
  character (len=*),intent(in)                       :: NameVar

  integer, parameter :: pres_=1, dens_=2

  integer :: iLat, iLon, l, k
  real :: T, P, latsHeidi(NR), mltsHeidi(NT)

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if (DoTestMe) &
       write(*,*)NameSub,' starting with iSizeIn,jSizeIn,nVar,NameVar=',&
       iSizeIn,jSizeIn,nVar,NameVar

  if(NameVar /= 'p:rho') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(iSizeIn /= IONO_nTheta*2-1 .or. jSizeIn /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  Buffer_IIV = -1.0

  ! eden and rnht are defined on a nr,nt grid
  ! where do I get latitude and mlt on nr,nt grid?

  ! the ionosphere and magnetosphere grid are shifted by 1, such that the
  ! ionosphere grid has an extra point at the lower end (and 2 at the upper)

  do iLon=1,jo
     mltsHeidi(iLon) = LonFac(iLon) * cPi / 12.0 
  enddo
  mltsHeidi(jo+1) = mltsHeidi(1) + 2.0 * cPi

  do iLat=1,io
     latsHeidi(iLat) = Latfac(iLat+1) * cDegToRad
  enddo

  do iLat = 1, IONO_nTheta
     do iLon = 1, IONO_nPsi

        T = cPi/2.0 - IONO_NORTH_Theta(iLat,iLon)
        P = mod(IONO_NORTH_Psi(iLat,iLon) + cPi, cPi*2)

        if ((T < latsHeidi(1)).or.(T > latsHeidi(io))) then
           Buffer_IIV(iLat,iLon,:) = -1.0
        else 

           k = 1
           do while (T > latsHeidi(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mltsHeidi(l))
              l = l + 1
           enddo

           ! This takes the nearest cell, and does not do linear interpolation

           ! Add together pressures from H+ (2) and O+ (4)
           ! Convert from keV/cc to Pa
           Buffer_IIV(iLat,iLon,pres_) = &
                eden(k,l,2)*0.1602*1.0e-9 + &
                eden(k,l,4)*0.1602*1.0e-9

           ! Add together density from H+ (2) and O+ (4)
           ! Convert from #/cc to kg/m3
           Buffer_IIV(iLat,iLon,dens_) = &
                rnht(k,l,2)*1.0e6*cProtonMass + &
                rnht(k,l,4)*1.0e6*cProtonMass*16.0

        endif

     enddo

  enddo

  do iLat = 1, IONO_nTheta
     do iLon = 1, IONO_nPsi

        T = IONO_SOUTH_Theta(iLat,iLon) - cPi/2
        P = mod(IONO_SOUTH_Psi(iLat,iLon) + cPi, cPi*2)

        if ((T < latsHeidi(1)).or.(T > latsHeidi(io))) then
           Buffer_IIV(iLat,iLon,:) = -1.0
        else 

           k = 1
           do while (T > latsHeidi(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mltsHeidi(l))
              l = l + 1
           enddo

           if (l > 1) l = l - 1

           ! This takes the nearest cell, and does not do linear interpolation

           ! Add together pressures from H+ (2) and O+ (4)
           ! Convert from keV/cc to Pa
           Buffer_IIV(iLat,iLon,pres_) = &
                eden(k,l,2)*0.1602*1.0e-9 + &
                eden(k,l,4)*0.1602*1.0e-9

           ! Add together density from H+ (2) and O+ (4)
           ! Convert from #/cc to kg/m3
           Buffer_IIV(iLat,iLon,dens_) = &
                rnht(k,l,2)*1.0e6*cProtonMass + &
                rnht(k,l,4)*1.0e6*cProtonMass*16.0

        endif

     enddo

  enddo

  ! species = e, H, he, o

!!! RNHT(colat,mlt,species) = density in #/cc
!!! EDEN("                ) = equatorial pressure (keV/cc) (*0.1602 = nPa)

end subroutine IM_get_for_gm

!==============================================================================

subroutine IM_init_session(iSession, TimeSimulation)
  use ModHeidiIO, ONLY: time
  implicit none
  
  integer,  intent(in) :: iSession       ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  logical :: IsUninitialized = .true.
  !--------------------------------------------------------------------------

  Time = TimeSimulation 

  if(IsUninitialized)then
     call heidi_init

     IsUninitialized = .false.
  end if

end subroutine IM_init_session

!==============================================================================
subroutine IM_finalize(TimeSimulation)

  use ModProcIM
  use ModInit, ONLY:nS
  use ModHeidiIO, ONLY :iUnitSw1,iUnitSw2,&
       iUnitMpa,iUnitSopa,iUnitPot,iUnitSal
  
  implicit none
  
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  !--------------------------------------------------------------------------

  close(iUnitSal)           ! Closes continuous output file
  close(iUnitSw1)           ! Closes sw1 input file
  close(iUnitSw2)           ! Closes sw2 input file
  close(iUnitMpa)           ! Closes MPA input file
  close(iUnitSopa)          ! Closes SOPA input file
  close(iUnitPot)           ! Closes FPOT input file

end subroutine IM_finalize

!=============================================================================

subroutine IM_run(SimTime, SimTimeLimit)

  use ModHeidiSize, only: dt, dtMax

  implicit none

  real, intent(inout) :: SimTime   ! current time of component
  real, intent(in)    :: SimTimeLimit ! simulation time not to be exceeded

  !--------------------------------------------------------------------------
  Dt = min(DtMax, (SimTimeLimit - SimTime)/2 )

  call heidi_run 

  SimTime = SimTime + dt*2

end subroutine IM_run

!===========================================================================

subroutine IM_save_restart(TimeSimulation)
  implicit none

  real,     intent(in)        :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='IM_save_restart'
  !-------------------------------------------------------------------------
!!! call heidi_save_restart

end subroutine IM_save_restart

!===========================================================================


