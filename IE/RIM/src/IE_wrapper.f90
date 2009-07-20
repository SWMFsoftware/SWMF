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
    character (len=100) :: imffilename
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

       case ("#TEST")
          call read_var('UseTests',UseTests)
          call read_var('TestName',TestName)

       case ("#HPI")
          call read_var('HemisphericPower',HemisphericPower)
          call IO_set_hpi_single(HemisphericPower)

       case ("#MHD_INDICES")
          cTempLines(1) = NameCommand
          call read_var('imffilename', imffilename)
          cTempLines(2) = imffilename
          cTempLines(3) = " "
          cTempLines(4) = "#END"
          call IO_set_inputs(cTempLines)
          call read_MHDIMF_Indices(iError)

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
          Tolerance = Tolerance*1000.0
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
          UseAmie = .true.

          cTempLines(1) = "#AMIEFILES"
          cTempLines(2) = AMIEFileNorth
          cTempLines(3) = AMIEFileSouth
          cTempLines(4) = ""
          cTempLines(5) = ""
          cTempLines(6) = ""
          cTempLines(7) = "#DEBUG"
          cTempLines(8) = "0"
          cTempLines(9) = "0"
          cTempLines(10) = ""
          cTempLines(11) = "#END"

          call EIE_set_inputs(cTempLines)

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

          cTempLines(1) = "#BACKGROUND"
          cTempLines(2) = "EIE/"
          cTempLines(3) = NameEFieldModel
          cTempLines(4) = NameAuroralModel
          cTempLines(5) = "idontknow"
          cTempLines(6) = ""
          cTempLines(7) = "#DEBUG"
          cTempLines(8) = "0"
          cTempLines(9) = "0"
          cTempLines(10) = ""
          cTempLines(11) = "#END"

          call EIE_set_inputs(cTempLines)

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
  use ModNumConst, only:cPi

  implicit none
  character (len=*), parameter :: NameSub='IE_set_grid'
  logical :: IsInitialized=.false.
  integer :: iProc_A(1)
  integer :: i

  logical :: DoTest, DoTestMe

  !------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  iProc_A(1)=0

  if (nLonsAll <= 0) then 
     allocate( &
          LatitudeAll(nLats+2,nLons+1), &
          LongitudeAll(nLats+2,nLons+1))
     LatitudeAll = 0.0
     LongitudeAll = 0.0
     nLonsAll = nLons
  endif

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
       nRootBlock_D=(/1,1/),             &! north+south hemispheres
       nCell_D =(/nLats+2,nLonsAll+1/),   &! size of node based grid
       XyzMin_D=(/cOne, cOne/),              &! min colat and longitude indexes
       XyzMax_D=(/real(nLats+2),             &! max colat and longitude indexes
                  real(nLonsAll+1)/),     &
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

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModRIM, ONLY: nLats, nLonsAll, &
       PotentialAll, OuterMagJrAll, SigmaHAll, SigmaPAll, &
       cpcpn, cpcps
  use ModParamRIM,    ONLY: TypeImCouple

  implicit none
  character(len=*), parameter :: NameSub='IE_get_for_im'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iBlock, iLat, iLon, iLatSouth, iPoint
  real    :: w, tSimulationTmp
  !---------------------------------------------------------------------------
  Buff_V = 0.0

  do iPoint = iPointStart, iPointStart + nPoint - 1

     iLat   = Index % iCB_II(1,iPoint)
     iLon   = Index % iCB_II(2,iPoint)
     iBlock = Index % iCB_II(3,iPoint)
     w      = Weight % Weight_I(iPoint)

!     write(*,*) "ilat, ilon : ",iLat, iLon, w

     if(iBlock/=1)then
        write(*,*)NameSub,': iPoint,Index % iCB_II=',&
             iPoint,Index%iCB_II(:,iPoint)
        call CON_stop(NameSub//&
             ' SWMF_ERROR iBlock should be 1=North in IE-IM coupling')
     end if

     if(iLat<1 .or. iLat>nLats+2 .or. iLon<1 .or. iLon>nLonsAll+1)then
        iLon = mod(iLon,nLonsAll+1)
        if(iLat<1 .or. iLat>nLats+2 .or. iLon<1 .or. iLon>nLonsAll+1)then
           call CON_stop(NameSub//' SWMF_ERROR index out of range')
        endif
     end if

     ! Index for the same latitude on the southern hemisphere
!     iLatSouth = nLats + 1 - iLat
     iLatSouth = iLat

     select case(TypeImCouple)
     case('north')
        Buff_V(1) = Buff_V(1) + w * PotentialAll(iLat,iLon)
        Buff_V(2) = Buff_V(2) + w * OuterMagJrAll(iLat,iLon)
        Buff_V(3) = Buff_V(3) + w * SigmaHAll(iLat,iLon)
        Buff_V(4) = Buff_V(4) + w * SigmaPAll(iLat,iLon)
     case('south')
        Buff_V(1) = Buff_V(1) + w * PotentialAll(iLatSouth,iLon)
        Buff_V(2) = Buff_V(2) + w * OuterMagJrAll(iLatSouth,iLon)
        Buff_V(3) = Buff_V(3) + w * SigmaHAll(iLatSouth,iLon)
        Buff_V(4) = Buff_V(4) + w * SigmaPAll(iLatSouth,iLon)
     case('cpcpmin')
        if(cpcpn < cpcps)then
           Buff_V(1) = Buff_V(1) + w * PotentialAll(iLat,iLon)
           Buff_V(2) = Buff_V(2) + w * OuterMagJrAll(iLat,iLon)
           Buff_V(3) = Buff_V(3) + w * SigmaHAll(iLat,iLon)
           Buff_V(4) = Buff_V(4) + w * SigmaPAll(iLat,iLon)
        else
           Buff_V(1) = Buff_V(1) + w * PotentialAll(iLatSouth,iLon)
           Buff_V(2) = Buff_V(2) + w * OuterMagJrAll(iLatSouth,iLon)
           Buff_V(3) = Buff_V(3) + w * SigmaHAll(iLatSouth,iLon)
           Buff_V(4) = Buff_V(4) + w * SigmaPAll(iLatSouth,iLon)
        end if
     case('average')
        Buff_V(1) = Buff_V(1) + w * &
             0.5*(PotentialAll(iLat,iLon)+PotentialAll(iLatSouth,iLon))
        Buff_V(2) = Buff_V(2) + w * &
             0.5*(OuterMagJrAll(iLat,iLon)+OuterMagJrAll(iLatSouth,iLon))
        Buff_V(3) = Buff_V(3) + w * 0.5*( &
             SigmaHAll(iLat,iLon)  + &
             SigmaHAll(iLatSouth,iLon))
        Buff_V(4) = Buff_V(4) + w * 0.5*( &
             SigmaPAll(iLat,iLon)  + &
             SigmaPAll(iLatSouth,iLon))
     case default
        call CON_stop(NameSub//' ERROR: Unknown value for TypeImCouple='// &
             TypeImCouple)
     end select
  end do

end subroutine IE_get_for_im

!==============================================================================
subroutine IE_put_from_im(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModRIM

  implicit none
  character(len=*), parameter   :: NameSub='IE_put_from_im'
  integer,intent(in)            :: nPoint, iPointStart, nVar
  real, intent(in)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight
  logical,intent(in)            :: DoAdd
  integer :: iBlock,iLat,iLon, iLM
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

  iLat = Index % iCB_II(1,iPointStart)
  iLon = Index % iCB_II(2,iPointStart)

!  if(iLat<1.or.iLat>nLats+2.or.iLon<0.or.iLon>nLonsAll+1)then
!     write(*,*)'iLat,iLon,DoAdd=',iLat,nLats,iLon,nLonsAll+1,DoAdd
!     call CON_stop('IE_put_from_im: index out of range')
!  end if

  if ( iLat >= 1 .and. iLat <= nLats+2 .and. &
       iLon >= 0 .and. iLon <=nLonsAll+1) then
     iLM = nLats+2 - iLat + 1
     if(DoAdd)then
        InnerMagJrAll(iLat,iLon)    = InnerMagJrAll(iLat,iLon)    + Buff_V(1)
        InnerMagEFluxAll(iLat,iLon) = InnerMagEFluxAll(iLat,iLon) + Buff_V(2)
        InnerMagAveEAll(iLat,iLon)  = InnerMagAveEAll(iLat,iLon)  + Buff_V(3)
!        if (iLat < nLats/2) then
!           InnerMagJrAll(iLM,iLon)    = InnerMagJrAll(iLM,iLon)    + Buff_V(1)
!           InnerMagEFluxAll(iLM,iLon) = InnerMagEFluxAll(iLM,iLon) + Buff_V(2)
!           InnerMagAveEAll(iLM,iLon)  = InnerMagAveEAll(iLM,iLon)  + Buff_V(3)
!        endif
     else
        InnerMagJrAll(iLat,iLon)    = Buff_V(1)
        InnerMagEFluxAll(iLat,iLon) = Buff_V(2)
        InnerMagAveEAll(iLat,iLon)  = Buff_V(3)
!        if (iLat < nLats/2) then
!           InnerMagJrAll(iLM,iLon)    = Buff_V(1)
!           InnerMagEFluxAll(iLM,iLon) = Buff_V(2)
!           InnerMagAveEAll(iLM,iLon)  = Buff_V(3)
!        endif
     end if
  endif

  IsNewInput = .true.

end subroutine IE_put_from_im

!==============================================================================

subroutine IE_put_from_im_complete

  write(*,*)"Don't know what IE_put_from_im_complete is really supposed to do."
  write(*,*)"I think that it is"
  write(*,*)"Supposed to be applying periodic boundaries...?"

end subroutine IE_put_from_im_complete

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

  ! Set the starttime, by getting the current time and subtracting off
  ! the simulation time (the subtraction is needed for restarts).
  call get_time(tCurrentOut = StartTime)
  StartTime = StartTime - tSimulation

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
  use ModParamRIM, only: iDebugLevel, DoSolve, UseTests
  use CON_physics, ONLY: get_time, get_axes, time_real_to_int
  use ModKind
  use ModMpi, only: mpi_wtime
  use ModKind, ONLY: Real8_

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: tSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  integer      :: nStep

  character(len=*), parameter :: NameSub='IE_run'

  logical :: DoTest,DoTestMe

  real (Real8_) :: TimingStart, TimingEnd

  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if (iDebugLevel > 2) DoTest = .true.

  if(DoTest)write(*,*)NameSub,': iProc,tSimulation,tSimulationLimit=',&
       iProc,tSimulation,tSimulationLimit

  if(DoTest)write(*,*)NameSub,': iProc,IsNewInput=',iProc,IsNewInput

  if (tSimulationLimit < 1.0e30) then
     tSimulation = tSimulationLimit
  endif

  ! Do not solve if there is no new input from GM or UA
  if (DoSolve .and. (.not.IsNewInput .and. .not.UseTests)) RETURN

  call timing_start('IE_run')
  TimingStart = mpi_wtime()

  CurrentTime = StartTime + tSimulation
  if (tSimulation == 0) OldTime = CurrentTime
  call time_real_to_int(CurrentTime, TimeArray)

  if (iDebugLevel >= 0) &
       write(*,"(a,i4,5i3,i4,a,i8,a)") &
       "IE Current Time (nSolve): ",TimeArray, " (",nSolve,")"

  ! Since IE is not a time dependent component, it may advance to the 
  ! next coupling time in a time accurate run
  if (IsTimeAccurate) tSimulation = tSimulationLimit

  ! Obtain the position of the magnetix axis
  call get_axes(tSimulation,MagAxisTiltGsmOut = ThetaTilt)

  call advance_RIM

  IsNewInput = .false.

  OldTime = CurrentTime

  call timing_stop('IE_run')
  TimingEnd = mpi_wtime()

  if (iDebugLevel > 1) write(*,*) "RIM==> Timing : ",TimingEnd-TimingStart

end subroutine IE_run

!=================================================================


