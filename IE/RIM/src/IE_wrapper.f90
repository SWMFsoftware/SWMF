! Wrapper for RIM
!==============================================================================
subroutine IE_set_param(CompInfo, TypeAction)

  use ModProcIE
  use ModRIM
  use ModIoRIM
  use ModParamRIM

  use ModIoUnit
  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='IE_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  real :: IMFBx, IMFBy, IMFBz, SWVx, HemisphericPower

  integer :: iError

  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use=.true.,                                    &
          NameVersion='Ridley Ionosphere Model (RIM)',&
          Version=0.1)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

     if( NameOutputDir(1:3) /= 'IE/' ) NameOutputDir = 'IE/'//NameOutputDir
  case('READ','CHECK')
     call read_param
  case('STDOUT')
     iUnitOut=STDOUT_
     if(nProc==1)then
        StringPrefix='IE:'
     else
        write(StringPrefix,'(a,i1,a)')'IE',iProc,':'
     end if
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call IE_set_grid
  case default
     call CON_stop(NameSub//' IE_ERROR: invalid TypeAction='//TypeAction)
  end select

contains

  subroutine read_param

    use ModReadParam
    use ModFiles
    use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case

    ! The name of the command
    character (len=100) :: NameCommand

    ! Read parameters
    logical :: DoEcho=.false., UseStrict=.true.

    ! Plot file parameters
    integer :: iFile, i, iError, iDebugProc
    character (len=50) :: plot_string

    character (len=100), dimension(100) :: cTempLines

    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('CHECK')
       ! We should check and correct parameters here
       if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()

       RETURN
    case('READ')
       if(iProc==0)write(*,*) NameSub,': READ iSession =',i_session_read(),&
            ' iLine=',i_line_read(),' nLine =',n_line_read()
    end select

    ! Read input data from text via ModReadParam
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#STRICT")
          call read_var('UseStrict',UseStrict)
       case("#OUTPUTDIR")
          call read_var("NameOutputDir",NameOutputDir)
          call fix_dir_name(NameOutputDir)
          if (iProc==0) call check_dir(NameOutputDir)
       case("#SAVEPLOT", "#IE_SAVEPLOT")
          call read_var('nPlotFile',nFile)
          if (nFile > MaxFile)call CON_stop(NameSub//&
               ' IE_ERROR number of ouput files is too large in #IE_SAVEPLOT:'&
               //' nFile>MaxFile')
          if (nFile>0.and.iProc==0) call check_dir(NameOutputDir)
          do iFile=1,nFile

             call read_var('StringPlot',plot_string)
             call lower_case(plot_string)

             ! Check to see if the ionosphere directory exists...
             if(iProc==0)call check_dir(NameOutputDir)

             ! Plotting frequency
             call read_var('DnSavePlot',dn_output(iFile))
             call read_var('DtSavePlot',dt_output(iFile))

             ! Plot file format
             if(index(plot_string,'idl')>0)then
                plot_form(iFile)='idl'
             elseif(index(plot_string,'tec')>0)then 
                plot_form(iFile)='tec'
             else
                call CON_stop(NameSub//&
                     ' IE_ERROR format (idl,tec) missing from plot_string='&
                     //plot_string)
             end if
             if(index(plot_string,'min')>0)then
                plot_vars(iFile)='min'
             elseif(index(plot_string,'max')>0)then
                plot_vars(iFile)='max'
             elseif(index(plot_string,'aur')>0)then
                plot_vars(iFile)='aur'
             elseif(index(plot_string,'uam')>0)then
                plot_vars(iFile)='uam'
             else
                call CON_stop(NameSub//&
                     ' IE_ERROR variable definition missing in #IE_SAVEPLOT'//&
                     ' from plot_string='//plot_string)
             end if
          end do

       case ("#SOLARWIND")
          call read_var('IMFBx',IMFBx)
          call read_var('IMFBy',IMFBy)
          call read_var('IMFBz',IMFBz)
          call read_var('SWVx',SWVx)
          call IO_set_imf_by_single(IMFby)
          call IO_set_imf_bz_single(IMFbz)
          call IO_set_sw_v_single(abs(SWvx))
          UseStaticIMF = .true.

        case ("#HPI")
           call read_var('HemisphericPower',HemisphericPower)
           call IO_set_hpi_single(HemisphericPower)

        case ("#MHD_INDICES")
           cTempLines(1) = NameCommand
           if(read_line()) then
              cTempLines(2) = NameCommand
              cTempLines(3) = " "
              cTempLines(4) = "#END"
              call IO_set_inputs(cTempLines)
              call read_MHDIMF_Indices(iError)
           else 
              iError = 1
           endif

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful"
              EXIT
           else
              UseStaticIMF = .false.
           endif

       case("#SOLVE")
          call read_var('DoSolve', DoSolve)
          call read_var('HighLatBoundary', HighLatBoundary)
          call read_var('LowLatBoundary', LowLatBoundary)
          call read_var('DoFold', DoFold)
          HighLatBoundary = HighLatBoundary * cDegToRad
          LowLatBoundary  = LowLatBoundary  * cDegToRad

       case("#IONOSPHERE")
          call read_var('iConductanceModel',iConductanceModel)
          call read_var('F10.7 Flux',f107flux)
          call read_var('StarLightPedConductance',StarLightPedConductance)
          call read_var('PolarCapPedConductance',PolarCapPedConductance)

       case("#IM")
          call read_var('TypeImCouple',TypeImCouple)
          call lower_case(TypeImCouple)

       case("#KRYLOV")
          call read_var('UsePreconditioner',UsePreconditioner)
          call read_var('UseInitialGuess',UseInitialGuess)
          call read_var('Tolerance',Tolerance)
          call read_var('MaxIteration',MaxIteration)

       case("#DEBUG")
          call read_var('iDebugLevel',iDebugLevel)
          call read_var('iDebugProc',iDebugProc)
          if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
             iDebugLevel = -1
          endif

       case("#AMIEFILES")
          call read_var('NameAmieFileNorth',AMIEFileNorth)
          call read_var('NameAmieFileSouth',AMIEFileSouth)
          NameEFieldModel = "amie"

       case("#BACKGROUND")
          call read_var('NameEFieldModel',NameEFieldModel)
          call read_var('NameAuroralModel',NameAuroralModel)
          call read_var('NameSolarModel',NameSolarModel)

          if (index(NameAuroralModel,'IHP') > 0) &
               NameAuroralModel = 'ihp'
          if (index(NameAuroralModel,'PEM') > 0) &
               NameAuroralModel = 'pem'

          if (index(NameEFieldModel,'AMIE') > 0) &
               NameEFieldModel = 'amie'

          if (index(NameEFieldModel,'weimer01') > 0) &
               NameEFieldModel = 'weimer01'
          if (index(NameEFieldModel,'Weimer01') > 0) &
               NameEFieldModel = 'weimer01'
          if (index(NameEFieldModel,'WEIMER01') > 0) &
               NameEFieldModel = 'weimer01'

          if (index(NameEFieldModel,'weimer') > 0 .and. &
               index(NameEFieldModel,'01') == 0) &
               NameEFieldModel = 'weimer96'
          if (index(NameEFieldModel,'Weimer') > 0 .and. &
               index(NameEFieldModel,'01') == 0) &
               NameEFieldModel = 'weimer96'
          if (index(NameEFieldModel,'WEIMER') > 0 .and. &
               index(NameEFieldModel,'01') == 0) &
               NameEFieldModel = 'weimer96'

          if (index(NameEFieldModel,'weimer96') > 0) &
               NameEFieldModel = 'weimer96'
          if (index(NameEFieldModel,'Weimer96') > 0) &
               NameEFieldModel = 'weimer96'
          if (index(NameEFieldModel,'WEIMER96') > 0) &
               NameEFieldModel = 'weimer96'

          if (index(NameEFieldModel,'SAMIE') > 0) &
               NameEFieldModel = 'samie'

       case("#SAVELOGFILE")
          call read_var('DoSaveLogfile',DoSaveLogfile)
          if(DoSaveLogfile)then
             if(iProc==0)call check_dir(NameOutputDir)
          endif

       case default
          if(iProc==0) then
             write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                  ' invalid command '//trim(NameCommand)
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if
       end select
    end do

    if (.not.DoSolve) then
       HighLatBoundary = 0.0
       LowLatBoundary  = 0.0
    endif

  end subroutine read_param

end subroutine IE_set_param
!=============================================================================
subroutine IE_set_grid

  ! Set the grid descriptor for IE
  ! Since IE has a static grid the descriptor has to be set once.
  ! There can be many couplers that attempt to set the descriptor,
  ! so we must check IsInitialized.
  use ModProcIE
  use ModRIM
  use CON_coupler
  use ModPlanetConst

  implicit none
  character (len=*), parameter :: NameSub='IE_set_grid'
  logical :: IsInitialized=.false.
  integer, allocatable :: iProc_A(:)
  integer :: i

  logical :: DoTest, DoTestMe

  !------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  allocate(iProc_A(nProc))
  do i=1,nProc
     iProc_A(i)=i-1
  end do

  !\
  ! When coupling, all models expect the ionosphere solution to go from
  ! the north pole to the south pole, and the "longitudes" to start from
  ! 12 MLT instead of 00 MLT (the other models interpret 0 as 12 MLT...)
  ! Further, we need to pad the solution at the north and south poles
  ! so we have a solution at exactly +/- 90.
  !/

  call set_grid_descriptor(                  &
       IE_,                                  &! component index
       nDim=2,                               &! dimensionality
       nRootBlock_D=(/1,nProc/),             &! north+south hemispheres
       nCell_D =(/nLons*nProc+1,nLats+2/),   &! size of node based grid
       XyzMin_D=(/cOne, cOne/),              &! min colat and longitude indexes
       XyzMax_D=(/real(nLons*nProc+1),       &
                  real(nLats+2)/),           &! max colat and longitude indexes
       TypeCoord='SMG',                      &! solar magnetic coord.
       Coord1_I=LatitudeAll(:,1),           &! colatitudes
       Coord2_I=LongitudeAll(1,:),            &! longitudes
       Coord3_I=(/rPlanet_I(Planet_) +       &
                  IonoHeightPlanet_I(Planet_)/),  &! radial size in meters
       iProc_A = iProc_A)                          ! processor assigment

end subroutine IE_set_grid

!==============================================================================

subroutine IE_get_for_gm(Buffer_II,iSize,jSize,tSimulation)

  use ModProcIE, only:nProc
  use ModSizeRIM
  use ModRIM

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_gm'

  integer, intent(in)           :: iSize,jSize
  real, intent(out)             :: Buffer_II(iSize,jSize)
  real (Real8_),    intent(in)  :: tSimulation

  integer :: i,j,k
  real    :: tSimulationTmp
  !--------------------------------------------------------------------------
  if(iSize /= nLats+2 .or. jSize /= nLons*nProc+1)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' nLats+2,nLons*nProc+1=',nLats+2,nLons*nProc+1
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call IE_run(tSimulationTmp,tSimulation)

  Buffer_II = PotentialAll

end subroutine IE_get_for_gm

!==============================================================================

subroutine IE_put_from_gm(Buffer_IIV,iSize,jSize,nVar)

  use ModRIM
  use ModProcIE
  use ModMpi

  implicit none
  character (len=*), parameter :: NameSub = 'IE_put_from_gm'
  integer,          intent(in) :: iSize, jSize, nVar
  real                         :: Buffer_IIV(iSize, jSize, nVar)

  integer :: iError, iMessageSize
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' starting'

  IsNewInput = .true.

  if (nProc > 1) then
     iMessageSize = iSize*jSize*nVar
     call MPI_Bcast(Buffer_IIV,iMessageSize,MPI_Real,0,iComm,iError)
  endif

  OuterMagJrAll = Buffer_IIV(:,:,1)
  if(nVar>1)then
     OuterMagInvBAll = Buffer_IIV(:,:,2)
     OuterMagRhoAll  = Buffer_IIV(:,:,3)
     OuterMagPAll    = Buffer_IIV(:,:,4)
  else
     OuterMagInvBAll = -1.0e32
     OuterMagRhoAll  = -1.0e32
     OuterMagPAll    = -1.0e32
  endif

  !\
  ! This seems like a total hack, but the latitude boundary
  ! is stored in this region, so we have to zero it out....
  !/
  LatBoundaryGm = Buffer_IIV(nLats/2+1,1,1)
  OuterMagJrAll(nLats/2+1:nLats/2+2,1) = 0.0

  if(DoTest)write(*,*)NameSub,' finished'

end subroutine IE_put_from_gm

!==============================================================================

subroutine IE_get_for_rb(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
     tSimulation)

  use ModSizeRIM
  use ModRIM
!  use ModProcIE
!  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_rb'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
  character (len=*),intent(in)  :: NameHem
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar
  real    :: tSimulationTmp

  return

!!!     !--------------------------------------------------------------------------
!!!     if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
!!!        write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
!!!             ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
!!!        call CON_stop(NameSub//' SWMF_ERROR')
!!!     end if
!!!   
!!!     ! Make sure that the most recent result is provided
!!!     tSimulationTmp = tSimulation
!!!     call IE_run(tSimulationTmp,tSimulation)
!!!   
!!!     select case(NameHem)
!!!   
!!!     case('North')
!!!   
!!!        if(iProc /= 0) RETURN
!!!        do iVar = 1, nVar
!!!           select case(Name_V(iVar))
!!!           case('Pot')
!!!              Buffer_IIV(:,:,iVar) = IONO_NORTH_Phi
!!!           case default
!!!              call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
!!!           end select
!!!        end do
!!!   
!!!     case('South')
!!!   
!!!        if(iProc /= nProc - 1) RETURN
!!!        do iVar = 1, nVar
!!!           select case(Name_V(iVar))
!!!           case('Pot')
!!!              Buffer_IIV(:,:,iVar) = IONO_SOUTH_Phi
!!!           case default
!!!              call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
!!!           end select
!!!        end do
!!!   
!!!     case default
!!!   
!!!        call CON_stop(NameSub//' invalid NameHem='//NameHem)
!!!   
!!!     end select

end subroutine IE_get_for_rb

!==============================================================================

subroutine IE_get_for_pw(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
     tSimulation)

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_pw'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
  character (len=*),intent(in)  :: NameHem
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_pw


!==============================================================================

subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide potential and current for IM
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V

  use CON_coupler,   ONLY: IndexPtrType, WeightPtrType
  use ModRIM, ONLY: nLats, nLons, &
       Potential, Jr, SigmaH, SigmaP, &
       cpcpn, cpcps
  use ModParamRIM,    ONLY: TypeImCouple

  implicit none
  character(len=*), parameter :: NameSub='IE_get_for_im'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iBlock, iLat, iLon, iLatSouth, iPoint
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
             ' SWMF_ERROR iBlock should be 1=North in IE-IM coupling')
     end if

     if(iLat<1 .or. iLat>nLats .or. iLon<1 .or. iLon>nLons+1)then
        write(*,*)'iLat,iLon=',iLat,iLon
        call CON_stop(NameSub//' SWMF_ERROR index out of range')
     end if

     ! Index for the same latitude on the southern hemisphere
     iLatSouth = nLats + 1 - iLat

     select case(TypeImCouple)
     case('north')
        Buff_V(1) = Buff_V(1) + w * Potential(iLon,iLat)
        Buff_V(2) = Buff_V(2) + w * Jr(iLon,iLat)
        Buff_V(3) = Buff_V(3) + w * SigmaH(iLon,iLat)
        Buff_V(4) = Buff_V(4) + w * SigmaP(iLon,iLat)
     case('south')
        Buff_V(1) = Buff_V(1) + w * Potential(iLon,iLatSouth)
        Buff_V(2) = Buff_V(2) + w * Jr(iLon,iLatSouth)
        Buff_V(3) = Buff_V(3) + w * SigmaH(iLon,iLatSouth)
        Buff_V(4) = Buff_V(4) + w * SigmaP(iLon,iLatSouth)
     case('cpcpmin')
        if(cpcpn < cpcps)then
           Buff_V(1) = Buff_V(1) + w * Potential(iLon,iLat)
           Buff_V(2) = Buff_V(2) + w * Jr(iLon,iLat)
           Buff_V(3) = Buff_V(3) + w * SigmaH(iLon,iLat)
           Buff_V(4) = Buff_V(4) + w * SigmaP(iLon,iLat)
        else
           Buff_V(1) = Buff_V(1) + w * Potential(iLon,iLatSouth)
           Buff_V(2) = Buff_V(2) + w * Jr(iLon,iLatSouth)
           Buff_V(3) = Buff_V(3) + w * SigmaH(iLon,iLatSouth)
           Buff_V(4) = Buff_V(4) + w * SigmaP(iLon,iLatSouth)
        end if
     case('average')
        Buff_V(1) = Buff_V(1) + w * &
             0.5*(Potential(iLon,iLat)+Potential(iLon,iLatSouth))
        Buff_V(2) = Buff_V(2) + w * &
             0.5*(Jr(iLon,iLat)+Jr(iLon,iLatSouth))
        Buff_V(3) = Buff_V(3) + w * 0.5*( &
             SigmaH(iLon,iLat)  + &
             SigmaH(iLon,iLatSouth))
        Buff_V(4) = Buff_V(4) + w * 0.5*( &
             SigmaP(iLon,iLat)  + &
             SigmaP(iLon,iLatSouth))
     case default
        call CON_stop(NameSub//' ERROR: Unknown value for TypeImCouple='// &
             TypeImCouple)
     end select
  end do

end subroutine IE_get_for_im

!==============================================================================

subroutine IE_get_for_ps(Buffer_IIV, iSize, jSize, nVar)

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_ps'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  !NOTE: The Buffer variables must be collected to i_proc0(IE_) before return.

  write(*,*) NameSub,' -- called but not yet implemented.'


end subroutine IE_get_for_ps
!==============================================================================

subroutine initialize_ie_ua_buffers(iOutputError)

  implicit none

  integer :: iOutputError

  character (len=*),parameter :: NameSub='initialize_ie_ua_buffers'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine initialize_ie_ua_buffers

!==============================================================================

subroutine IE_put_from_UA(Buffer_III, iBlock, nMLTs, nLats, nVarsToPass)

  implicit none

  integer, intent(in) :: nMlts, nLats, iBlock, nVarsToPass
  real, dimension(nMlts, nLats, nVarsToPass), intent(in) :: Buffer_III

  character (len=*),parameter :: NameSub='IE_put_from_UA'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_put_from_UA

!==============================================================================

subroutine IE_get_for_ua(Buffer_II,iSize,jSize,NameVar,NameHem,tSimulation)

  implicit none

  integer,          intent(in)  :: iSize,jSize
  real,             intent(out) :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  character (len=*),intent(in)  :: NameHem
  real,             intent(in)  :: tSimulation

  character (len=*),parameter :: NameSub='IE_get_for_ua'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_get_for_ua

!==============================================================================

subroutine SPS_put_into_ie(Buffer_II, iSize, jSize, NameVar, iBlock)

  implicit none

  integer, intent(in)           :: iSize,jSize
  real, intent(in)              :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  integer,intent(in)            :: iBlock

  character (len=*), parameter :: NameSub='SPS_put_into_ie'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine SPS_put_into_ie


!==============================================================================

subroutine IE_init_session(iSession, tSimulation)

  ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

  use CON_physics, ONLY: get_time, get_planet, get_axes
  use ModRIM,      ONLY: IsTimeAccurate, ThetaTilt, DipoleStrength, StartTime
  use ModIoRIM,    ONLY: dt_output, t_output_last
  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: tSimulation   ! seconds from start time

  !DESCRIPTION:
  ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

  character(len=*), parameter :: NameSub='IE_init_session'

  logical :: IsUninitialized=.true.

  logical :: DoTest,DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(IsUninitialized)then
     call init_RIM
     IsUninitialized = .false.
  end if

  call get_time(  DoTimeAccurateOut = IsTimeAccurate)
  call get_planet(DipoleStrengthOut = DipoleStrength)
  call get_axes(tSimulation, MagAxisTiltGsmOut = ThetaTilt)

  call get_time(tCurrentOut = StartTime)

  DipoleStrength = DipoleStrength*1.0e9 ! Tesla -> nT

  write(*,*)NameSub,': DipoleStrength, ThetaTilt =',DipoleStrength,ThetaTilt,&
       StartTime

  ! Reset t_output_last in case the plotting frequency has changed
  if(IsTimeAccurate)then
     where(dt_output>0.) &
          t_output_last=int(tSimulation/dt_output)
  end if

end subroutine IE_init_session

!==============================================================================
subroutine IE_finalize(tSimulation)

  use ModProcIE
  use ModRIM, ONLY: TimeArray, nSolve
  use ModIoRIM, ONLY: nFile
  use CON_physics, ONLY: get_time
  use ModTimeConvert, ONLY: time_real_to_int
  use ModKind, ONLY: Real8_
  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_finalize'

  integer :: iFile
  real(Real8_) :: tCurrent
  !---------------------------------------------------------------------------
  call get_time(tCurrentOut = tCurrent)
  call time_real_to_int(tCurrent, TimeArray)

  if(nSolve>0)then
     do iFile=1,nFile
        call write_output_RIM(iFile)
     end do
  end if

end subroutine IE_finalize

!==============================================================================

subroutine IE_save_restart(tSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_save_restart'

  RETURN

end subroutine IE_save_restart

!==============================================================================

subroutine IE_run(tSimulation,tSimulationLimit)

  use ModProcIE
  use ModRIM
  use ModParamRIM, only: iDebugLevel
  use CON_physics, ONLY: get_time, get_axes, time_real_to_int
  use ModKind
  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: tSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  integer      :: nStep

  character(len=*), parameter :: NameSub='IE_run'

  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if (iDebugLevel > 2) DoTest = .true.

  if(DoTest)write(*,*)NameSub,': iProc,tSimulation,tSimulationLimit=',&
       iProc,tSimulation,tSimulationLimit

  if(DoTest)write(*,*)NameSub,': iProc,IsNewInput=',iProc,IsNewInput

  ! Do not solve if there is no new input from GM or UA
  if(.not.IsNewInput) RETURN

  CurrentTime = StartTime + tSimulation
  call time_real_to_int(CurrentTime, TimeArray)

  if (iDebugLevel >= 0) write(*,*) "Current Time : ",TimeArray

  ! Since IE is not a time dependent component, it may advance to the 
  ! next coupling time in a time accurate run
  if (IsTimeAccurate) tSimulation = tSimulationLimit

  ! Obtain the position of the magnetix axis
  call get_axes(tSimulation,MagAxisTiltGsmOut = ThetaTilt)

  call advance_RIM

  IsNewInput = .false.

!   ! Solve for the ionosphere potential
!   call IE_solve

!   ! Save logfile if required
!   call IE_save_logfile

end subroutine IE_run

!=================================================================


