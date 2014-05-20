!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IM_wrapper

  ! Wrapper for HEIDI Internal Magnetosphere (IM) component

  implicit none

  private ! except

  public:: IM_set_param
  public:: IM_init_session
  public:: IM_run
  public:: IM_save_restart
  public:: IM_finalize

  ! Coupling with IE
  public:: IM_get_for_ie
  public:: IM_put_from_ie_mpi
  public:: IM_put_from_ie
  public:: IM_put_from_ie_complete

  ! Coupling with GM
  public:: IM_get_for_gm
  public:: IM_put_from_gm
  public:: IM_put_from_gm_line
  public:: IM_put_from_gm_crcm
  public:: IM_put_sat_from_gm

contains

  !============================================================================

  subroutine IM_set_param(CompInfo,TypeAction)

    use CON_comp_info
    use ModProcHeidi
    use ModHeidiMain
    use ModHeidiSize, ONLY: tMax
    use ModReadParam, ONLY: i_session_read
    use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case
    use ModHeidiIO!,   ONLY: IsFramework, StringPrefix
    use ModIoUnit,    ONLY: STDOUT_


    character (len=*), parameter      :: NameSub='IM_set_param'
    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*),  intent(in)   :: TypeAction ! What to do

    !Local Variables:
    character (len=100)               :: NameCommand, StringPlot
    logical                           :: DoEcho=.false.
    logical                           :: UseStrict=.true.  
    integer                           :: iUnitOut
    !--------------------------------------------------------------------------

    IsBFieldNew = .true.
    write(*,*) 'HEIDI set_param Action=', TypeAction

    select case(TypeAction)
    case('VERSION')
       write(*,*) 'VERSION'

       call put(CompInfo,                         &
            Use=.true.,                           &
            NameVersion='RAM_HEIDI',              &
            Version=1.1)

    case('MPI')
       write(*,*) 'MPI'

       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       if(nProc>4)call CON_stop( NameSub // &
            ' IM_ERROR this version can run on 4 PE !')
       IsFramework = .true.


    case('CHECK')
       write(*,*) 'CHECK'

       !We should check and correct parameters here
       if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
       call heidi_check


    case('GRID')
       write(*,*) 'GRID'

       call IM_set_grid  

    case('READ')
       write(*,*) 'READ'

       call heidi_read

       write(*,*) 'READ--------->DT,TMAX,TINT,TIME',DT,TMAX,TINT,TIME

    case('STDOUT')
       write(*,*) 'STDOUT'

       iUnitOut = STDOUT_
       if(nProc==1)then
          StringPrefix='IM:'
       else
          write(StringPrefix,'(a,i3.3,a)')'IM',iProc,':'
       end if

    case('FILEOUT')
       write(*,*) 'FILEOUT'

       call get(CompInfo,iUnitOut=iUnitOut)
       StringPrefix=''

    case default
       call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)

    end select

  end subroutine IM_set_param

  !============================================================================
  subroutine IM_set_grid

    use ModNumConst,  ONLY: cTwoPi
    use CON_coupler,  ONLY: set_grid_descriptor, is_proc, IM_
    use ModHeidiSize, ONLY: RadiusMin, RadiusMax,NT,NR
    use ModHeidiMain, ONLY: LZ, Z, DL1, DPHI,PHI

    character (len=*), parameter :: NameSub='IM_set_grid'
    logical                      :: IsInitialized=.false.
    logical                      :: DoTest, DoTestMe
    integer                      :: i, j  
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)'IM_set_grid called, IsInitialized=', &
         IsInitialized
    if(IsInitialized) return

    IsInitialized = .true.

    !\
    ! IM grid: the equatorial grid is described by Coord1_I and Coord2_I
    ! Occasional +0.0 is used to convert from single to double precision
    !/

    call set_grid_descriptor( IM_,           & ! component index
         nDim     = 2,                       & ! dimensionality
         nRootBlock_D = (/1,1/),             & ! number of blocks
         nCell_D =(/nR, nT-1/),              & ! size of equatorial grid
         XyzMin_D=(/RadiusMin+0.0,0.0/),     & ! min coordinates
         XyzMax_D=(/RadiusMax+0.0,cTwoPi/),  & ! max coordinates
         Coord1_I = LZ(1:nR)+0.0,            & ! radial coordinates
         Coord2_I = Phi(1:nT-1)+0.0,         & ! longitudinal coordinates
         TypeCoord= 'SMG',                   & ! solar magnetic coord
         nVar = 2,                           & ! number of "fluid" vars                         
         NameVar = "p rho")                    ! names of "fluid" vars

    if(DoTest)then
       write(*,*)NameSub,' NR = ', NR
       write(*,*)NameSub,' NT = ', NT
    end if


  end subroutine IM_set_grid
  !============================================================================
  subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)
    !\
    ! Provide current for IE
    ! The value should be interpolated from nPoints with
    ! indexes stored in Index and weights stored in Weight
    ! The variables should be put into Buff_V(??)
    !/
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    use ModIonoHeidi, ONLY: &
         IONO_NORTH_RCM_JR,IONO_SOUTH_RCM_JR, IONO_nTheta, IONO_nPsi

    character(len=*), parameter    :: NameSub='IM_get_for_ie'
    integer,intent(in)             :: nPoint, iPointStart, nVar
    real,intent(out)               :: Buff_V(nVar)
    type(IndexPtrType),intent(in)  :: Index
    type(WeightPtrType),intent(in) :: Weight
    integer                        :: iLat, iLon, iBlock, iPoint
    real                           :: w
    !--------------------------------------------------------------------------
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

       ! Only worry about the northern hemisphere....  
       ! IE can fix the southern hemisphere.
       if (iLat <= IONO_nTheta .and. iLon <= IONO_nPsi) &
            Buff_V(1) = Buff_V(1) + w * IONO_NORTH_RCM_JR(iLat,iLon)

       if (iLat > IONO_nTheta .and. iLon <= IONO_nPsi) &
            Buff_V(1) = Buff_V(1) &
            + w * IONO_SOUTH_RCM_JR(2*IONO_nTheta-iLat+1,iLon)

    end do

  end subroutine IM_get_for_ie

  !============================================================================
  subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)

    use ModHeidiIO,  ONLY: time
    use ModPlotFile, ONLY: save_plot_file

    integer, intent(in):: nTheta, nPhi
    real,    intent(in):: Potential_II(nTheta, nPhi, 1)
    character(len=100) :: NameFile
    !-------------------------------------------------------------------------
    write(NameFile,'(a,i5.5,a)') &
         "IM/plots/potential_t",nint(Time),".out"

    call save_plot_file(NameFile, &
         StringHeaderIn = 'Ionospheric potential', &
         TimeIn         = time+0.0, &
         NameVarIn      = 'Theta Phi Pot', &
         CoordMinIn_D   = (/0.0, 0.0/), &
         CoordMaxIn_D   = (/180.0,360.0/), &
         VarIn_IIV = Potential_II)

  end subroutine IM_put_from_ie_mpi

  !============================================================================

  subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    use ModIonoHeidi, ONLY: &
         IONO_NORTH_PHI, IONO_SOUTH_PHI, IONO_nTheta, IONO_nPsi

    character(len=*), parameter   :: NameSub='IM_put_from_ie'
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd
    integer :: iBlock,i,j
    !--------------------------------------------------------------------------
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

  !============================================================================

  subroutine IM_put_from_ie_complete

    !--------------------------------------------------------------------------

  end subroutine IM_put_from_ie_complete

  !============================================================================

  subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

    ! This should be similar to RBE coupling

    use ModIonoHeidi
    use ModConst

    character (len=*),parameter :: NameSub='IM_put_from_gm'

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
    real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
    character (len=*),intent(in)       :: NameVar

    integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
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
  !============================================================================
  subroutine IM_put_from_gm_line(nRadiusIn, nLonIn, Map_DSII, &
       nVarLineIn, nPointLineIn, BufferLine_VI, NameVar)

    use ModHeidiMain,      ONLY: &
         nR, nT, LZ, BHeidi_III, SHeidi_III, RHeidi_III,&
         bGradB1xHeidi_III,bGradB1yHeidi_III, bGradB1zHeidi_III,&
         BxHeidi_III, ByHeidi_III, BzHeidi_III,Xyz_VIII, pHeidi_III, &
         RhoHeidi_III, MhdEqPressure_I, MhdEqDensity_I
    use ModHeidiMain,      ONLY: Phi, IsBFieldNew
    use ModHeidiIO,        ONLY: Time
    use ModHeidiSize,      ONLY: RadiusMin, RadiusMax, iPointBMin_II, iPointEq
    use ModIoUnit,         ONLY: UnitTmp_
    use ModPlotFile,       ONLY: save_plot_file
    use ModHeidiBField,    ONLY: dipole_length
    use ModCoordTransform, ONLY: xyz_to_sph
    use ModNumConst,       ONLY: cPi
    use ModPlanetConst,    ONLY: DipoleStrengthPlanet_I, rPlanet_I, Earth_

    integer, intent(in)          :: nRadiusIn, nLonIn
    real,    intent(in)          :: Map_DSII(3,2,nRadiusIn,nLonIn)
    integer, intent(in)          :: nVarLineIn, nPointLineIn
    real,    intent(in)          :: BufferLine_VI(nVarLineIn,nPointLineIn)
    character(len=*), intent(in) :: NameVar
    character(len=*), parameter  :: NameSub='IM_put_from_gm_line'
    character(len=100)           :: NameFile
    logical                      :: IsFirstCall = .true.
    logical                      :: DoTest, DoTestMe
    !\
    ! These variables should either be in a module, OR
    ! there is no need for them, and BufferLine_VI should be put 
    ! into HEIDI variables right here. 
    ! Note that this routine is only called on the root processor !!!
    !/8

    integer                :: nVarLine   = 0    ! number of vars per line point
    integer                :: nPointLine = 0    ! number of points in all lines
    real, save, allocatable:: StateLine_VI(:,:)       ! state along all lines
    integer,save           :: iRiTiDIr_DI(3,2*nR*nT)  ! line index 
    character(LEN=500)     :: StringVarName, StringHeader
    character(len=20)      :: TypePosition
    character(len=20)      :: TypeFile = 'ascii'
    !\
    ! Local Variables
    !/
    integer, parameter:: &
         I_=1, S_=2, X_=3, Y_=4, Z_=5, rho_= 6, ux_=7, uy_=8, uz_=9, p_=13
    integer, parameter:: Bx_=10, By_=11, Bz_=12, gx_=14, gy_=15, gz_=16 
    integer, parameter:: nStepInside = 10, nStepInterp = 40
    integer, parameter:: nStep = 2*(nStepInside + nStepInterp)+1
    real,    parameter:: rBoundary = 3.0

    real, dimension(3,nStepInside+1)     :: bDipoleS_VI,bDipoleN_VI,XyzDipoleN_VI,XyzDipoleS_VI
    real, dimension(nStepInside+1)       :: sDipoleS_I, sDipoleN_I,rDipoleS_I,rDipoleN_I
    real, dimension(nStepInside+1)       :: bDipoleMagnS_I, bDipoleMagnN_I
    real, dimension(nStepInside+1,nR,nT) :: BDipoleN_III, BxDipoleN_III, ByDipoleN_III, BzDipoleN_III
    real, dimension(3,nStep)             :: XyzDipole_VI, bDipole_VI
    real, dimension(nStep)               :: bDipoleMagn_I, sDipole_I, rDipole_I
    real, dimension(3,nStepInside,nR,nT) :: XyzDipoleN_VIII,XyzDipoleS_VIII
    real, dimension(nStepInside,nR,nT)   :: BDipoleS_III, BxDipoleS_III, ByDipoleS_III, BzDipoleS_III
    real, dimension(nStep,nR,nT)         :: BDipoleMagn_III, STemp
    real, dimension(3,nStep,nR,nT)       :: bDipole_VIII
    real, dimension(nStepInside,nR)      :: rDipoleN_II, sDipoleN_II, rDipoleS_II, sDipoleS_II
    real, dimension(nStep,nR)            :: rDipole_II, sDipole_II 
    real, dimension(nStep,nR)            :: bGradB1x_II, bGradB1y_II
    real, dimension(nStepInterp)         :: LengthHeidi_I,BHeidi_I,RHeidi_I,LengthHeidinew_I
    real, dimension(nStepInterp)         :: XHeidi_I,YHeidi_I,ZHeidi_I,XHeidinew_I,XHeidi1new_I
    real, dimension(nStepInterp)         :: LengthHeidi1new_I
    real, dimension(nStepInterp)         :: BxHeidi_I,ByHeidi_I,BzHeidi_I
    real, dimension(nStepInterp)         :: pHeidi_I, rhoHeidi_I
    real, dimension(nStepInterp)         :: bGradB1xHeidi_I, bGradB1yHeidi_I, bGradB1zHeidi_I
    real, allocatable                    :: B_I(:), Length_I(:),RadialDist_I(:)
    real, allocatable                    :: bGradB1x_I(:), bGradB1y_I(:), bGradB1z_I(:)
    real, allocatable                    :: Bx_I(:), By_I(:), Bz_I(:)
    real, allocatable                    :: X_I(:),Y_I(:),Z_I(:)
    real, allocatable                    :: p_I(:), rho_I(:)
    real                                 :: LatBoundaryN, LatBoundaryS
    real                                 :: LatMax, LatMin, Lat, dLat,x,y,z,a,dLength
    real                                 :: Tr,Ttheta, r,gradB0R1,gradB0R2, gradB0Theta1,gradB0Theta2
    integer                              :: iStep,ns,np,k
    integer                              :: iR, iT, iDir, n
    integer                              :: iPoint,ip, iPhi
    integer                              :: iMax, i, iLineLast,iLine,iLineFirst,j
    real                                 :: sMax
    real                                 :: LatDipole(nStep,nR),LatDipoleN(nStepInside,nR),LatDipoleS(nStepInside,nR)
    real                                 :: xS, yS, zS, xN, yN, zN
    real                                 :: LengthExS(nR,nT,2), LengthExN(nR,nT,2)
    real                                 :: Re, DipoleFactor
    !--------------------------------------------------------------------------
    Re = rPlanet_I(Earth_)
    DipoleFactor = DipoleStrengthPlanet_I(Earth_)*(Re)**3

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    !\
    ! Save total number of points along all field lines
    !/
    nPointLine = nPointLineIn
    nVarLine   = nVarLineIn

    !\
    ! Allocate buffer
    !/
    if (allocated(StateLine_VI)) deallocate(StateLine_VI)
    allocate(StateLine_VI(nVarLine,nPointLine))
    !\
    ! Copy into local variables
    !/
    StateLine_VI = BufferLine_VI

    if(DoTest)then
       write(*,*)NameSub,' nVarLine,nPointLine=',nVarLine,nPointLine

       ! Set the file name
       write(NameFile,'(a,i5.5,a)') &
            "IM/plots/ray_data_t",nint(Time),".out"
       open(UnitTmp_, FILE=NameFile, STATUS="replace")
       ! Same format as in GM/BATSRUS/src/ray_trace_new.f90
       write(UnitTmp_, *) 'nRadius, nLon, nPoint=',nR, nT, nPointLine
       write(UnitTmp_, *) 'iLine l x y z rho ux uy uz bx by bz p bgradb1x bgradb1y bgradb1z'
       do iPoint = 1, nPointLine
          write(UnitTmp_, *) StateLine_VI(:, iPoint)
       end do
       close(UnitTmp_)
       !\
       ! Now save the mapping files (+0.0 for real precision)
       !/
       write(NameFile,'(a,i5.5,a)') &
            "IM/plots/map_north_t",nint(Time),".out"

       call save_plot_file( &
            NameFile, &
            StringHeaderIn = 'Mapping to northern ionosphere', &
            TimeIn       = Time+0.0, &
            NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
            CoordMinIn_D = (/RadiusMin+0.0,   0.0/), &
            CoordMaxIn_D = (/RadiusMax+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,1,:,:))

       write(NameFile,'(a,i5.5,a)') &
            "IM/plots/map_south_t",nint(Time),".out"
       call save_plot_file( &
            NameFile, &
            StringHeaderIn = 'Mapping to southern ionosphere', &
            TimeIn       = Time+0.0, &
            NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
            CoordMinIn_D = (/RadiusMin+0.0,   0.0/), &
            CoordMaxIn_D = (/RadiusMax+0.0, 360.0/), &
            VarIn_VII    = Map_DSII(:,2,:,:))
    end if
    !\
    ! Convert Units here. Input is in SI !!!
    ! Check Map_DSII for open-closed field lines, also use it for mapping
    ! to the ionosphere for electric potential.
    !/

    !Create index array that converts radial and local time index to line index
    iPoint = 0
    do iT = 1, nT
       do iR = 1, nR
          do iDir = 1, 2
             iPoint =iPoint +1
             iRiTiDir_DI(:,iPoint) = (/iR, iT, iDir/)
          end do
       end do
    end do

    !\
    ! Count the maximum size of the field line array (iMax)
    !/
    iMax = 0 ; iLineLast = -1; i = 1
    do iPoint =1 ,nPointLine
       if (StateLine_VI(1,iPoint) == iLineLast) then
          i = i + 1
          iMax = max(iMax,i)
       else
          i =1
       end if
       iLineLast = StateLine_VI(1,iPoint)
    end do

    allocate(X_I(iMax), Y_I(iMax), Z_I(iMax))
    allocate(B_I(iMax), Length_I(iMax), RadialDist_I(iMax))
    allocate(bGradB1x_I(iMax), bGradB1y_I(iMax), bGradB1z_I(iMax));
    allocate(Bx_I(iMax), By_I(iMax), Bz_I(iMax));
    allocate(p_I(iMax), rho_I(iMax));

    iLineFirst = StateLine_VI(1,1)
    iLineLast = -1
    i = 1
    j = 1

    LengthExN=0.0
    LengthExS=0.0
    sTemp = 0.0 

    do iPoint = 1, nPointLine !+ 1
       ! Check if this is the point after the last one, or if it is a new line segment

       !     if(iPoint > nPointLine .or. &
       !          (iPoint > 1 .and. StateLine_VI(1,min(nPointLine,iPoint)) /= iLineLast)) then

       if(iPoint > 1 .and. StateLine_VI(1,iPoint) /= iLineLast) then

          np = i-1  ! np = number of points on each half field line

          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(1:np-1), X_I(2:np), XHeidinew_I,LengthHeidinew_I)
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   X_I(2:np), XHeidi_I, LengthHeidi_I)
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   Y_I(2:np), YHeidi_I, LengthHeidi_I)
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   Z_I(2:np), ZHeidi_I, LengthHeidi_I)
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   B_I(2:np), BHeidi_I, LengthHeidi_I)
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   RadialDist_I(2:np),RHeidi_I,LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   bGradB1x_I(2:np), bGradB1xHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   bGradB1y_I(2:np), bGradB1yHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   bGradB1z_I(2:np), bGradB1zHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   Bx_I(2:np),  BxHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   By_I(2:np),  ByHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(3,(np-1),nStepInterp,Length_I(2:np),   Bz_I(2:np),  BzHeidi_I, LengthHeidi_I) 

          call interpolate_mhd(2,(np-1),nStepInterp,Length_I(2:np),   p_I(2:np),   PHeidi_I, LengthHeidi_I) 
          call interpolate_mhd(2,(np-1),nStepInterp,Length_I(2:np),   Rho_I(2:np), RhoHeidi_I, LengthHeidi_I) 

          dLength =  Length_I(np-1)/nStepInterp

          iLine = StateLine_VI(1,iPoint-1)
          iR   = iRiTiDir_DI(1,iLine)
          iT   = iRiTiDir_DI(2,iLine)
          iDir = iRiTiDir_DI(3,iLine)

          ! if (iLine ==9) write(*,*) 'iLine, B_I =', iLine, B_I
          ! if (iLine ==90) write(*,*) 'iLine, B_I =', iLine, B_I
          ! if (iLine ==30) write(*,*) 'iLine, B_I =', iLine, B_I
          ! if (iR==7 .and. iT==1 )  write(*,*) 'iR, iT, iDir, B in MHD', iR, iT, iDir, B_I
          ! if (iR==18 .and. iT==1 )  write(*,*) 'iR, iT, iDir, B in MHD', iR, iT, iDir, B_I


          if (iDir ==1) then  ! Northern hemisphere
             dLength =  Length_I(np-1)/(nStepInterp-1)
             do k = 1, nStepInterp
                LengthHeidi_I(k) = k * dLength
             end do
             Xyz_VIII(1,(iPointEq+ 1):(nStep - nStepInside),iR,iT) = XHeidi_I
             Xyz_VIII(2,(iPointEq+ 1):(nStep - nStepInside),iR,iT) = YHeidi_I
             Xyz_VIII(3,(iPointEq+ 1):(nStep - nStepInside),iR,iT) = ZHeidi_I
             BHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT) = BHeidi_I
             STemp((iPointEq+ 1):(nStep - nStepInside),iR,iT)      = LengthHeidi_I
             RHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT) = RHeidi_I
             bGradB1xHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT) = bGradB1xHeidi_I
             bGradB1yHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT) = bGradB1yHeidi_I
             bGradB1zHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT) = bGradB1zHeidi_I
             BxHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT)       = BxHeidi_I(:)
             ByHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT)       = ByHeidi_I(:)
             BzHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT)       = BzHeidi_I(:)

             pHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT)      = pHeidi_I(:)
             rhoHeidi_III((iPointEq+ 1):(nStep - nStepInside),iR,iT)    = RhoHeidi_I(:)
             LengthExN(iR,iT,iDir) = Length_I(2)
          end if

          if (iDir ==2) then     ! Southern hemisphere
             LengthExS(iR,iT,iDir) = Length_I(2)
             dLength =  Length_I(np-1)/(nStepInterp-1)
             do k = 1, nStepInterp
                LengthHeidi_I(k) = k * dLength
             end do
             Xyz_VIII(1,(nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = XHeidi_I(nStepInterp:1:-1)
             Xyz_VIII(2,(nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = YHeidi_I(nStepInterp:1:-1)
             Xyz_VIII(3,(nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = ZHeidi_I(nStepInterp:1:-1)
             BHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = BHeidi_I(nStepInterp:1:-1)
             STemp((nStepInside + 1):(nStepInside + nStepInterp),iR,iT)  = LengthHeidi_I(:) 
             RHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = RHeidi_I(nStepInterp:1:-1) 
             bGradB1xHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = bGradB1xHeidi_I(nStepInterp:1:-1) 
             bGradB1yHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = bGradB1yHeidi_I(nStepInterp:1:-1) 
             bGradB1zHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = bGradB1zHeidi_I(nStepInterp:1:-1) 
             BxHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = BxHeidi_I(nStepInterp:1:-1) 
             ByHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = ByHeidi_I(nStepInterp:1:-1) 
             BzHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = BzHeidi_I(nStepInterp:1:-1) 

             pHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT)   = pHeidi_I(nStepInterp:1:-1) 
             rhoHeidi_III((nStepInside + 1):(nStepInside + nStepInterp),iR,iT) = rhoHeidi_I(nStepInterp:1:-1)

          end if

          Xyz_VIII(1,iPointEq,iR,iT) = X_I(1)
          Xyz_VIII(2,iPointEq,iR,iT) = Y_I(1)
          Xyz_VIII(3,iPointEq,iR,iT) = Z_I(1)
          BHeidi_III(iPointEq,iR,iT) = B_I(1)
          STemp(iPointEq,iR,iT) = 0.5*(LengthExS(iR,iT,2) + LengthExN(iR,iT,1))
          RHeidi_III(iPointEq,iR,iT) = RadialDist_I(1)
          bGradB1xHeidi_III(iPointEq,iR,iT) = bGradB1x_I(1)
          bGradB1yHeidi_III(iPointEq,iR,iT) = bGradB1y_I(1)
          bGradB1zHeidi_III(iPointEq,iR,iT) = bGradB1z_I(1)
          BxHeidi_III(iPointEq,iR,iT) = Bx_I(1)
          ByHeidi_III(iPointEq,iR,iT) = By_I(1)
          BzHeidi_III(iPointEq,iR,iT) = Bz_I(1)

          pHeidi_III(iPointEq,iR,iT)   = p_I(1)
          rhoHeidi_III(iPointEq,iR,iT) = rho_I(1)

          i = 1
          iLineLast = StateLine_VI(1,iPoint)

          !if (iPoint > nPointLine) EXIT

       end if

!!$     write(*,*) '~~~~~~~~~~~~~~~~~~~~'
!!$     write(*,*)  'StateLine_VI(:,iPoint) = ',  StateLine_VI(:,iPoint)
!!$     write(*,*) '~~~~~~~~~~~~~~~~~~~~'

       X_I(i) = StateLine_VI(X_,iPoint)
       Y_I(i) = StateLine_VI(Y_,iPoint)
       Z_I(i) = StateLine_VI(Z_,iPoint)
       B_I(i) = sqrt(StateLine_VI(BX_,iPoint)**2+ &
            StateLine_VI(BY_,iPoint)**2 + StateLine_VI(BZ_,iPoint)**2)
       Length_I(i) = StateLine_VI(S_,iPoint)
       RadialDist_I(i) = sqrt(StateLine_VI(X_,iPoint)**2 + &
            StateLine_VI(Y_,iPoint)**2 + StateLine_VI(Z_,iPoint)**2 )
       bGradB1x_I(i) = StateLine_VI(gx_,iPoint)
       bGradB1y_I(i) = StateLine_VI(gy_,iPoint)
       bGradB1z_I(i) = StateLine_VI(gz_,iPoint)
       Bx_I(i)       = StateLine_VI(BX_,iPoint)
       By_I(i)       = StateLine_VI(BY_,iPoint)
       Bz_I(i)       = StateLine_VI(BZ_,iPoint)
       P_I(i)        = StateLine_VI(p_,iPoint)
       Rho_I(i)      = StateLine_VI(rho_, iPoint)

       iLineLast = StateLine_VI(1,iPoint)
       i = i + 1


    end do

    deallocate(X_I);         deallocate(Y_I);         deallocate(Z_I);
    deallocate(B_I);         deallocate(Length_I);    deallocate(RadialDist_I);
    deallocate(bGradB1x_I);  deallocate(bGradB1y_I);  deallocate(bGradB1z_I);
    deallocate(Bx_I);        deallocate(By_I);        deallocate(Bz_I);
    deallocate(p_I) ;        deallocate(rho_I);

    !\
    ! Convert the MHD variables (B field) to Heidi grid
    !/

    ns = nStepInside

    do iT = 1,nT-1
       do iR = 1, nR
          if (LZ(iR) > rBoundary) then
             xS = Xyz_VIII(1, nStepInside+1, iR, iT)
             yS = Xyz_VIII(2, nStepInside+1, iR, iT)
             zS = Xyz_VIII(3, nStepInside+1, iR, iT)

             xN = Xyz_VIII(1, nStep - nStepInside, iR, iT)
             yN = Xyz_VIII(2, nStep - nStepInside, iR, iT)
             zN = Xyz_VIII(3, nStep - nStepInside, iR, iT)

             if (xS == 0.0) xS =  xN
             if (yS == 0.0) yS =  yN
             if (zS == 0.0) xS = -zN

             LatBoundaryS = -atan(zS/(sqrt(xS**2 + yS**2)))
             LatBoundaryN =  atan(zN/(sqrt(xN**2 + yN**2)))

             call fill_dipole_north(nStepInside+1, LZ(iR), Phi(iT), LatBoundaryN, XyzDipoleN_VI, bDipoleN_VI,&
                  bDipoleMagnN_I, sDipoleN_I, rDipoleN_I)
             call fill_dipole_south(nStepInside+1, LZ(iR), Phi(iT), LatBoundaryS, XyzDipoleS_VI, bDipoleS_VI,&
                  bDipoleMagnS_I,sDipoleS_I, rDipoleS_I)

             Xyz_VIII(:,1:nStepInside,iR,iT)                 = XyzDipoleS_VI(:,1:nStepInside)
             Xyz_VIII(:,(nStep-nStepInside +1):nStep,iR,iT)  = XyzDipoleN_VI(:,2:nStepInside+1)

             BHeidi_III(1:nStepInside,iR,iT)                 = bDipoleMagnS_I(1:nStepInside)
             BHeidi_III((nStep-nStepInside +1):nStep,iR,iT)  = bDipoleMagnN_I(2:nStepInside+1)

             RHeidi_III(1:nStepInside,iR,iT)                 = rDipoleS_I(1:nStepInside)   
             RHeidi_III((nStep-nStepInside +1):nStep,iR,iT)  = rDipoleN_I(2:nStepInside+1)

             SHeidi_III(1:nStepInside,iR,iT)                 = sDipoleS_I(1:nStepInside)   

             do i = nStepInside+1, nStepInside + nStepInterp
                SHeidi_III(i,iR,iT) = sDipoleS_I(nStepInside) + STemp(i,iR,iT)
             end do

             do i = iPointEq, iPointEq
                SHeidi_III(i,iR,iT) = sDipoleS_I(nStepInside) + STemp(iPointEq-1,iR,iT)+ STemp(iPointEq,iR,iT)
             end do

             do i = nStepInside + nStepInterp + 2, nStep - nStepInside
                SHeidi_III(i,iR,iT) = SHeidi_III(nStepInside + nStepInterp+1, iR, iT) + STemp(i,iR,iT) 
             end do

             do i = nStep-nStepInside+1 , nStep
                j = i - (nStep-nStepInside) 
                SHeidi_III(i,iR,iT) = SHeidi_III(nStep-nStepInside,iR,iT) + sDipoleN_I(j)
             end do

             !\
             ! Fill in dipole values for mhd field lines beyond rBoundary: LatBoundary = atan(z/sqrt(x^2+y^2))
             !/

             BxHeidi_III(1:nStepInside,iR,iT)                = bDipoleS_VI(1,1:nStepInside) 
             BxHeidi_III((nStep-nStepInside +1):nStep,iR,iT) = bDipoleN_VI(1,2:nStepInside+1)

             ByHeidi_III(1:nStepInside,iR,iT)                = bDipoleS_VI(2,1:nStepInside) 
             ByHeidi_III((nStep-nStepInside +1):nStep,iR,iT) = bDipoleN_VI(2,2:nStepInside+1)

             BzHeidi_III(1:nStepInside,iR,iT)                = bDipoleS_VI(3,1:nStepInside) 
             BzHeidi_III((nStep-nStepInside +1):nStep,iR,iT) = bDipoleN_VI(3,2:nStepInside+1)

             pHeidi_III(1:nStepInside,iR,iT)                = pHeidi_III(nStepInside+1,iR,iT)  
             pHeidi_III((nStep-nStepInside +1):nStep,iR,iT) = pHeidi_III(nStep-nStepInside,iR,iT)

             RhoHeidi_III(1:nStepInside,iR,iT)                = RhoHeidi_III(nStepInside+1,iR,iT)  
             RhoHeidi_III((nStep-nStepInside +1):nStep,iR,iT) = RhoHeidi_III(nStep-nStepInside,iR,iT)


          end if

          if (LZ(iR) <= rBoundary) then
             call fill_dipole(nStep, LZ(iR), Phi(iT), XyzDipole_VI, bDipole_VI, bDipoleMagn_I, sDipole_I, rDipole_I)
             Xyz_VIII(:,:,iR,iT)  = XyzDipole_VI(:,:)
             BHeidi_III(:,iR,iT)  = bDipoleMagn_I(:)
             RHeidi_III(:,iR,iT)  = rDipole_I(:)
             SHeidi_III(:,iR,iT)  = sDipole_I(:)
             BxHeidi_III(:,iR,iT) = bDipole_VI(1,:)
             ByHeidi_III(:,iR,iT) = bDipole_VI(2,:)
             BzHeidi_III(:,iR,iT) = bDipole_VI(3,:)
          end if

       end do  ! L loop
    end do     ! Phi loop

    ! Fill in periodic boundary conditions
    Xyz_VIII(:,:,:,nT)        = Xyz_VIII(:,:,:,1)
    BHeidi_III(:,:,nT)        = BHeidi_III(:,:,1)
    RHeidi_III(:,:,nT)        = RHeidi_III(:,:,1)
    SHeidi_III(:,:,nT)        = SHeidi_III(:,:,1)
    BxHeidi_III(:,:,nT)       = BxHeidi_III(:,:,1)
    ByHeidi_III(:,:,nT)       = ByHeidi_III(:,:,1)
    BzHeidi_III(:,:,nT)       = BzHeidi_III(:,:,1)
    bGradB1xHeidi_III(:,:,nT) = bGradB1xHeidi_III(:,:,1)
    bGradB1yHeidi_III(:,:,nT) = bGradB1yHeidi_III(:,:,1)  
    bGradB1zHeidi_III(:,:,nT) = bGradB1zHeidi_III(:,:,1)
    pHeidi_III(:,:,nT)        = pHeidi_III(:,:,1)
    rhoHeidi_III(:,:,nT)      = rhoHeidi_III(:,:,1)

    ! This should be gone !!!
    Xyz_VIII(:,:,nR,:)        = Xyz_VIII(:,:,nR-1,:)
    BHeidi_III(:,nR,:)        = BHeidi_III(:,nR-1,:)
    RHeidi_III(:,nR,:)        = RHeidi_III(:,nR-1,:)
    SHeidi_III(:,nR,:)        = SHeidi_III(:,nR-1,:)
    BxHeidi_III(:,nR,:)       = BxHeidi_III(:,nR-1,:)
    ByHeidi_III(:,nR,:)       = ByHeidi_III(:,nR-1,:)
    BzHeidi_III(:,nR,:)       = BzHeidi_III(:,nR-1,:)
    bGradB1xHeidi_III(:,nR,:) = bGradB1xHeidi_III(:,nR-1,:)
    bGradB1yHeidi_III(:,nR,:) = bGradB1yHeidi_III(:,nR-1,:)  
    bGradB1zHeidi_III(:,nR,:) = bGradB1zHeidi_III(:,nR-1,:)
    rhoHeidi_III(:,nR,:)      = rhoHeidi_III(:,nR-1,:)
    pHeidi_III(:,nR,:)        = pHeidi_III(:,nR-1,:)


    !Find the location of minimum B
    do iR = 1, nR
       do iT = 1, nT
          iPointBmin_II(iR,iT) = minloc(BHeidi_III(:,iR,iT), 1) 
       end do
    end do


    !Get the pressure, density and velocity at the HEIDI boundary in the equatorial plane

    do iT = 1, nT
       MhdEqPressure_I(iT) = pHeidi_III(iPointBmin_II(nR,iT),nR,iT)
       MhdEqDensity_I(iT)    = rhoHeidi_III(iPointBmin_II(nR,iT),nR,iT)
       ! MhdEqVelocity = sqrt(StateLine_VI(ux_,iPointBmin)**2 + &
       !      StateLine_VI(uy_,iPointBmin)**2 + StateLine_VI(uz_,iPointBmin)**2)


       write(*,*) 'pressure and rho from MHD: iT, P, rho = ', iT,  MhdEqPressure_I(iT),  MhdEqDensity_I(iT)


    end do

    !STOP

    !~~~~~~~~~~~~~~~~~~~~~~~ Write out files for testing ~~~~~~~~~~~~~~~~~~~~~~

    NameFile      = 'BFieldMagn.out'
    StringHeader  = 'Magnetic field in the equatorial plane'
    StringVarName = 'R MLT B'
    TypePosition  = 'rewind'

    call save_plot_file(NameFile, & 
         TypePositionIn = TypePosition,&
         TypeFileIn     = TypeFile,&
         StringHeaderIn = StringHeader, &
         nStepIn = 0, &
         TimeIn = 0.0, &
         NameVarIn = StringVarName, &
         nDimIn = 2, & 
         CoordMinIn_D = (/1.75, 0.0/),&
         CoordMaxIn_D = (/6.5, 24.0/),&
         VarIn_VII = BHeidi_III(iPointEq:iPointEq,:,:))
    TypePosition = 'rewind' 


    ! Use UnitTmp_ (or 9)
    open(unit=9,file='interface.out')
    write(9,*) 'from interface'
    write(9,*)'iPoint  iR  b  bx  by  bz  x  y  z  s Lat'

    do iPhi =1, 1
       do iR = 1,20
          do iPoint =1, nStep
             write(9,'(2I5,  9E11.3)') iPoint, iR, BHeidi_III(iPoint,iR,iPhi),&
                  BxHeidi_III(iPoint,iR,iPhi), ByHeidi_III(iPoint,iR,iPhi),&
                  BzHeidi_III(iPoint,iR,iPhi), Xyz_VIII(1,iPoint,iR,iPhi),Xyz_VIII(2,iPoint,iR,iPhi),&
                  Xyz_VIII(3,iPoint,iR,iPhi),SHeidi_III(iPoint,iR,iPhi),&
                  atan(Xyz_VIII(3,iPoint,iR,iPhi)/(sqrt(Xyz_VIII(1,iPoint,iR,iPhi)**2+ Xyz_VIII(2,iPoint,iR,iPhi)**2)))
          end do
       end do
    end do
    close(9)

    !  call stop_mpi('DEBUG')



  end subroutine IM_put_from_gm_line

  !============================================================================

  subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
    ! Puts satellite locations and names from GM into IM variables.

    use ModHeidiSatellites
    use ModNumConst,   ONLY: cDegToRad

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

  !============================================================================

  subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

    use CON_time, ONLY : get_time
    use ModNumConst, ONLY: cPi, cDegToRad
    use ModConst, ONLY: cProtonMass
    use ModHeidiCurrents, ONLY:eden, rnht
    use ModHeidiSize,  ONLY:  nR, nT

    character (len=*),parameter :: NameSub='IM_get_for_gm'

    integer, intent(in)                                :: iSizeIn,jSizeIn,nVar
    real, dimension(iSizeIn,jSizeIn,nVar), intent(out) :: Buffer_IIV
    character (len=*),intent(in)                       :: NameVar
    integer, parameter :: pres_=1, dens_=2
    integer :: iLat, iLon, l, k
    real :: T, P, latsHeidi(NR), mltsHeidi(NT)
    logical :: DoTest, DoTestMe
    real    :: Pmin
    integer :: i,j, iSize, jSize
    !--------------------------------------------------------------------------
    iSize = nR
    jSize = nT-1


    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if (DoTestMe) &
         write(*,*)NameSub,' starting with iSizeIn, jSizeIn, nVar, NameVar=',&
         iSizeIn,jSizeIn,nVar,NameVar

    if(NameVar /= 'p:rho') &
         call CON_stop(NameSub//' invalid NameVar='//NameVar)

    if(iSizeIn /= iSize .or. jSizeIn /= jSize)then
       write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    Buffer_IIV = 0.0

    !Fill pressure and density
    do i=1,iSize
       do j=1,jSize

          !write(*,*) 'i,j', i, j
          !write(*,*) 'iSize,jSize = ', iSize, jSize
          !write(*,*) 'i, j, eden(i,j,:)=', i, j, eden(i,j,1), eden(i,j,2), rnht(i,j,1), rnht(i,j,2)
          !STOP

          ! make sure pressure passed to GM is not lower than Pmin [nPa]
          ! to avoid too low GM pressure 
          Pmin = minval(eden(i,j,:))
          Buffer_IIV(i,j,pres_) = max(sum(eden(i,j,:)), Pmin)*1e-9

          ! Add together density from H+, He+ and O+
          ! Convert from #/cc to kg/m3
          Buffer_IIV(i,j,dens_) = &
               rnht(i,j,1)*1.0e6*cProtonMass * 4 + & ! He+
               rnht(i,j,2)*1.0e6*cProtonMass + &     ! H+
               rnht(i,j,4)*1.0e6*cProtonMass*16.0    ! O+

          ! Only a not-a-number can be less than zero and larger than one
          if(  .not. Buffer_IIV(i,j,pres_) > 0 .and. &
               .not. Buffer_IIV(i,j,pres_) < 1) then
             write(*,*)NameSub,': ERROR IN PRESSURE'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,pres_)
             call CON_stop(NameSub // ' ERROR: Not a number found in IM pressure !')
          end if
          if(  .not. Buffer_IIV(i,j,dens_) > 0 .and. &
               .not. Buffer_IIV(i,j,dens_) < 1) then
             write(*,*)NameSub,': ERROR IN DENSITY'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,dens_)
             call CON_stop(NameSub // ' ERROR: Not a number found in IM density !')
          end if
       end do
    end do

    if(DoTestMe)write(*,*) NameSub,' finished'


    !STOP
    ! species = e, H, he, o
    ! RNHT(colat,mlt,species) = density in #/cc
    ! EDEN("                ) = equatorial pressure (keV/cc) (*0.1602 = nPa)

  end subroutine IM_get_for_gm

  !============================================================================

  subroutine IM_init_session(iSession, TimeSimulation)
    use ModHeidiIO, ONLY: time

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

  !============================================================================

  subroutine IM_finalize(TimeSimulation)

    use ModProcHeidi
    use ModInit, ONLY:nS
    use ModHeidiIO, ONLY :iUnitSw1,iUnitSw2,&
         iUnitMpa,iUnitSopa,iUnitPot,iUnitSal

    real,     intent(in) :: TimeSimulation   ! seconds from start time
    !--------------------------------------------------------------------------

    close(iUnitSal)           ! Closes continuous output file
    close(iUnitSw1)           ! Closes sw1 input file
    close(iUnitSw2)           ! Closes sw2 input file
    close(iUnitMpa)           ! Closes MPA input file
    close(iUnitSopa)          ! Closes SOPA input file
    close(iUnitPot)           ! Closes FPOT input file

  end subroutine IM_finalize

  !============================================================================

  subroutine IM_run(SimTime, SimTimeLimit)

    use ModHeidiSize,  ONLY: dt, dtMax, tmax
    use ModHeidiMain,  ONLY: t, IsBFieldNew
    use ModHeidiIO,    ONLY: time
    use ModInit,       ONLY: i3

    real, intent(inout) :: SimTime      ! current time of component
    real, intent(in)    :: SimTimeLimit ! simulation time not to be exceeded

    !--------------------------------------------------------------------------

    dT = min(DtMax, (SimTimeLimit - SimTime)/2 )

    write(*,*) 'WRAPPER TIMES====>Dt, DtMax, SimTimeLimit, SimTime:', &
         Dt, DtMax, SimTimeLimit, SimTime

    write(*,*)'IM_run before heidi_run t,time,i3=',t,time,i3

    call heidi_run 


    i3 = i3 + 1

    write(*,*)'IM_run after  heidi_run t,time,i3=',t,time,i3

    SimTime = SimTime + dT*2

    write(*,*) 'WRAPPER: SimTime=', SimTime

  end subroutine IM_run

  !===========================================================================

  subroutine IM_save_restart(TimeSimulation)

    real,     intent(in)        :: TimeSimulation   ! seconds from start time
    character(len=*), parameter :: NameSub='IM_save_restart'
    !-------------------------------------------------------------------------
!!! call heidi_save_restart

  end subroutine IM_save_restart

  !===========================================================================
  subroutine IM_put_from_gm_crcm(Integral_IIV,iSizeIn,jSizeIn,nIntegralIn,&
       BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)

    integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
    real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
    integer, intent(in) :: nVarLine, nPointLine
    real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)
    real,    intent(in) :: tSimulation
    character (len=*), intent(in) :: NameVar

    character (len=*), parameter :: NameSub='IM_put_from_gm_crcm'

    call CON_stop(NameSub//': CRCM version cannot be used for HEIDI!')

  end subroutine IM_put_from_gm_crcm


  !===========================================================================

  subroutine interpolate_linear_b( &
       nP, Length_I, B_I, nPoint, LengthHeidi_I, BHeidi_I)


    integer, intent(in)  :: nP ! number of points along the MHD field line
    real,    intent(in)  :: Length_I(nP)
    real,    intent(in)  :: B_I(nP)
    integer, intent(in)  :: nPoint ! number of points on new grid
    real,    intent(out) :: LengthHeidi_I(nPoint)
    real,    intent(out) :: BHeidi_I(nPoint)

    !Local variables
    real    :: dLength, LengthMax,LengthMin
    integer :: iP, iPoint, i
    !--------------------------------------------------------------------------

    LengthMax = Length_I(nP)
    LengthMin = Length_I(1)
    dLength = (LengthMax - LengthMin)/(nPoint - 1)

    ! Linear Interpolation
    do iPoint = 1, nPoint
       LengthHeidi_I(iPoint) = LengthMin + (iPoint - 1) * dLength

       do iP = 1, nP
          if (Length_I(iP) > LengthHeidi_I(iPoint)) then
             i = iP - 1
             BHeidi_I(iPoint) = B_I(i) + (LengthHeidi_I(iPoint) - Length_I(i))&
                  *(B_I(i+1) - B_I(i))/(Length_I(i+1) - Length_I(i) )
             EXIT
          end if
       end do
    end do

  end subroutine interpolate_linear_b

  !============================================================================
  real function lagrange(lHeidi, lMhd_I, bMhd_I, nStepMHD, nOrder)

    real    :: lHeidi ! find the value of B at this point along the fiedl line
    integer :: nStepMHD         ! number of points along the MHD field line
    real    :: lMhd_I(nStepMHD) ! field line length values from MHD
    real    :: bMhd_I(nStepMHD) ! magnetic fiedl values from MHD
    integer :: nOrder              ! order of interpolation
    real    :: func_I(nStepMHD)
    integer :: i, j, k, l, m
    real    :: y

    !-------------------------------------------------------------------------
    ! Check if the size of the array is larger than the order of interpolation 
    !if (nOrder > nStepMHD)  nOrder = nStepMHD

    ! Check if lHeidi is outside the lMhd(1)-lMhd(nStepMHD) interval. 
    ! If yes set a boundary value.
    if (lHeidi <= lMhd_I(1)) then
       lagrange = bMhd_I(1)
       return
    end if
    if (lHeidi >= lMhd_I(nStepMHD)) then
       lagrange = bMhd_I(nStepMHD)
       return
    end if

    ! Search to find i so that lMhd(i) < lHeidi < lMhd(i+1)
    i = 1
    j = nStepMHD
    do while (j > i+1)
       k = (i+j)/2
       if (lHeidi < lMhd_I(k)) then
          j = k
       else
          i = k
       end if
    end do

    ! shift i so that will correspond to n-th order of interpolation
    ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
    i = i + 1 - nOrder/2

    if (i < 1) i=1
    if (i + nOrder > nStepMHD) i = nStepMHD - nOrder + 1

    ! Lagrange interpolation
    y = 0.0
    do m = i, i + nOrder - 1
       func_I(m)=1.0
       do l = i, i + nOrder -1
          if(l /= m) func_I(m) = func_I(m) &
               *(lHeidi - lMhd_I(l))/(lMhd_I(m) - lMhd_I(l))
       end do
       y = y + bMhd_I(m)*func_I(m)
    end do
    lagrange = y

  end function lagrange

  !============================================================================

  subroutine interpolate_mhd( &
       nOrder, nStepMhd, nStep, lMhd_I, bMhd_I, bHeidi_I, lHeidi_I)

    integer, intent(in)  :: nOrder    ! order of interpolation
    integer, intent(in)  :: nStepMHD  ! number of points along MHD field line
    integer, intent(in)  :: nStep     ! number of points to interpolate onto
    real,    intent(in)  :: lMhd_I(nStepMHD) ! field line length from MHD
    real,    intent(in)  :: bMhd_I(nStepMHD) ! magnetic field values from MHD
    real,    intent(out) :: bHeidi_I(nStep)  ! magnetic field interpolated 
    real,    intent(out) :: lHeidi_I(nStep)  ! field line length values 
    ! Local Variables
    real    :: LengthMax, LengthMin, dLength
    integer :: iStep
    !-------------------------------------------------------------------------

    LengthMax = lMhd_I(nStepMhd)
    LengthMin = lMhd_I(1)

    dLength = (LengthMax - LengthMin)/(nStep - 1)

    do iStep = 1, nStep
       lHeidi_I(iStep) = LengthMin + (iStep - 1) * dLength
       bHeidi_I(iStep) = &
            lagrange(lHeidi_I(iStep),lMhd_I,bMhd_I,nStepMhd, nOrder+1)
    end do

  end subroutine interpolate_mhd

  !============================================================================

  subroutine fill_dipole(nStep, L, Phi, XyzDipole_VI, bDipole_VI, &
       bDipoleMagn_I, sDipole_I, rDipole_I)

    use ModHeidiBField,    ONLY: dipole_length
    use ModPlanetConst,    ONLY: DipoleStrengthPlanet_I, rPlanet_I, Earth_

    integer, intent(in)  :: nStep
    real,    intent(in)  :: L,Phi
    real,    intent(out) :: XyzDipole_VI(3,nStep)
    real,    intent(out) :: bDipole_VI(3,nStep), bDipoleMagn_I(nStep)
    real,    intent(out) :: sDipole_I(nStep),rDipole_I(nStep)
    real                 :: LatDipole_I(nStep)
    real                 :: LatMax, LatMin, Lat, dLat
    real                 :: r, x, y, z, a
    integer              :: iStep
    real                 :: Re, DipoleFactor
    !--------------------------------------------------------------------------

    Re = rPlanet_I(Earth_)
    DipoleFactor = DipoleStrengthPlanet_I(Earth_)*(Re)**3


    LatMax =  acos(sqrt(1./L))
    LatMin = -LatMax
    dLat   = (LatMax-LatMin)/(nStep-1)
    Lat = LatMin
    do iStep = 1, nStep
       r = Re * L * (cos(Lat))**2
       x = r * cos(Lat) * cos(Phi)
       y = r * cos(Lat) * sin(Phi)
       z = r * sin(Lat)
       a = (sqrt(x**2 + y**2 +z**2))**5

       XyzDipole_VI(1,iStep) = x 
       XyzDipole_VI(2,iStep) = y 
       XyzDipole_VI(3,iStep) = z

       bDipole_VI(1,iStep) = DipoleFactor * (3. * z * x)/a
       bDipole_VI(2,iStep) = DipoleFactor * (3. * z * y)/a
       bDipole_VI(3,iStep) = DipoleFactor * (2. * z**2 - x**2 - y**2)/a

       bDipoleMagn_I(iStep)  = sqrt((bDipole_VI(1,iStep))**2+&
            (bDipole_VI(2,iStep))**2 + (bDipole_VI(3,iStep))**2)
       sDipole_I(iStep) = dipole_length(Re * L ,LatMin,Lat) 
       rDipole_I(iStep) = r
       LatDipole_I(iStep) = Lat

       Lat = Lat + dLat
    end do

  end subroutine fill_dipole

  !============================================================================

  subroutine fill_dipole_south(&
       nStep, L, Phi, LatBoundary, XyzDipoleS_VI, bDipoleS_VI,&
       bDipoleMagnS_I, sDipoleS_I, rDipoleS_I)

    use ModHeidiBField,    ONLY: dipole_length
    use ModPlanetConst,    ONLY: DipoleStrengthPlanet_I, rPlanet_I, Earth_

    integer, intent(in)  :: nStep
    real,    intent(in)  :: L, LatBoundary, Phi
    real,    intent(out) :: XyzDipoleS_VI(3,nStep)
    real,    intent(out) :: bDipoleS_VI(3,nStep)
    real,    intent(out) :: sDipoleS_I(nStep)
    real,    intent(out) :: rDipoleS_I(nStep)
    real,    intent(out) :: bDipoleMagnS_I(nStep)
    real                 :: LatDipoleS_I(nStep)
    real                 :: LatMax, LatMin, Lat, dLat
    real                 :: x, y, z, r, a
    integer              :: iStep
    real                 :: Re, DipoleFactor
    !--------------------------------------------------------------------------

    Re = rPlanet_I(Earth_)
    DipoleFactor = DipoleStrengthPlanet_I(Earth_)*(Re)**3

    LatMin = LatBoundary
    LatMax =  acos(sqrt(1./L))
    dLat   = (LatMax-LatMin)/(nStep-1)         
    Lat = -LatMax

    do iStep = 1, nStep
       r = Re * L * (cos(Lat))**2
       x = r * cos(Lat) * cos(Phi)
       y = r * cos(Lat) * sin(Phi)
       z = r * sin(Lat)
       a = (sqrt(x**2 + y**2 +z**2))**5

       XyzDipoleS_VI(1,iStep) = x 
       XyzDipoleS_VI(2,iStep) = y 
       XyzDipoleS_VI(3,iStep) = z

       bDipoleS_VI(1,iStep)   = DipoleFactor * (3. * z * x)/a
       bDipoleS_VI(2,iStep)   = DipoleFactor * (3. * z * y)/a
       bDipoleS_VI(3,iStep)   = DipoleFactor * (2. * z**2 - x**2 - y**2)/a
       bDipoleMagnS_I(iStep)  = sqrt((bDipoleS_VI(1,iStep))**2+&
            (bDipoleS_VI(2,iStep))**2 + (bDipoleS_VI(3,iStep))**2)

       sDipoleS_I(iStep)      = dipole_length(Re*L, -LatMax, Lat) 
       rDipoleS_I(iStep)      = r
       LatDipoleS_I(iStep)    = Lat
       Lat = Lat + dLat

    end do


  end subroutine fill_dipole_south

  !============================================================================

  subroutine fill_dipole_north(&
       nStep, L, Phi, LatBoundary, XyzDipoleN_VI, bDipoleN_VI,&
       bDipoleMagnN_I, sDipoleN_I, rDipoleN_I)

    use ModHeidiBField,    ONLY: dipole_length
    use ModPlanetConst,    ONLY: DipoleStrengthPlanet_I, rPlanet_I, Earth_

    integer, intent(in)  :: nStep
    real,    intent(in)  :: L, LatBoundary, Phi
    real,    intent(out) :: XyzDipoleN_VI(3,nStep), bDipoleN_VI(3,nStep)
    real,    intent(out) :: sDipoleN_I(nStep)
    real,    intent(out) :: rDipoleN_I(nStep)
    real,    intent(out) :: bDipoleMagnN_I(nStep)
    real                 :: LatDipoleN_I(nStep)
    real                 :: LatMax, LatMin, Lat, dLat
    real                 :: x, y, z, r, a
    integer              :: iStep
    real                 :: Re, DipoleFactor
    !--------------------------------------------------------------------------

    Re = rPlanet_I(Earth_)
    DipoleFactor = DipoleStrengthPlanet_I(Earth_)*(Re)**3


    LatMin = LatBoundary
    LatMax =  acos(sqrt(1./L))
    dLat   = (LatMax-LatMin)/(nStep-1)         
    Lat = LatMin

    do iStep = 1, nStep
       r = Re * L * (cos(Lat))**2
       x = r * cos(Lat) * cos(Phi)
       y = r * cos(Lat) * sin(Phi)
       z = r * sin(Lat)
       a = (sqrt(x**2 + y**2 +z**2))**5

       XyzDipoleN_VI(1,iStep) = x 
       XyzDipoleN_VI(2,iStep) = y 
       XyzDipoleN_VI(3,iStep) = z

       bDipoleN_VI(1,iStep)   = DipoleFactor * (3. * z * x)/a
       bDipoleN_VI(2,iStep)   = DipoleFactor * (3. * z * y)/a
       bDipoleN_VI(3,iStep)   = DipoleFactor * (2. * z**2 - x**2 - y**2)/a
       bDipoleMagnN_I(iStep)  = sqrt((bDipoleN_VI(1,iStep))**2+&
            (bDipoleN_VI(2,iStep))**2 + (bDipoleN_VI(3,iStep))**2)

       sDipoleN_I(iStep)      = dipole_length(Re*L, LatMin-dLat, Lat) 
       rDipoleN_I(iStep)      = r
       LatDipoleN_I(iStep)    = Lat
       Lat = Lat + dLat

    end do

  end subroutine fill_dipole_north
  !============================================================================

end module IM_wrapper

