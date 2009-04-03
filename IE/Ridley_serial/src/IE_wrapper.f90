! Wrapper for Ridley's ionosphere
!==============================================================================
subroutine IE_set_param(CompInfo, TypeAction)

  use ModProcIE
  use ModIonosphere
  use IE_ModIo
  use IE_ModMain

  use ModIoUnit
  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='IE_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: iError

  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use=.true.,                                    &
          NameVersion='Serial Potential Solver (Ridley)',&
          Version=1.1)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

     if( nProc>2 )call CON_stop(NameSub//' IE_ERROR '//&
          'this version can run on 1 or 2 PE-s only!')
     if( NameIonoDir(1:3) /= 'IE/' ) NameIonoDir = 'IE/'//NameIonoDir
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
    use ModIE_Interface
    use ModFiles
    use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case

    ! The name of the command
    character (len=100) :: NameCommand

    ! Read parameters
    logical :: DoEcho=.false., UseStrict=.true., IsUninitialized=.true.

    ! Plot file parameters
    integer :: iFile, i, iError, iDebugProc
    character (len=50) :: plot_string

    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('CHECK')
       if(IsUninitialized)call set_defaults
       IsUninitialized=.false.

       ! We should check and correct parameters here
       if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()

       RETURN
    case('READ')
       if(iProc==0)write(*,*) NameSub,': READ iSession =',i_session_read(),&
            ' iLine=',i_line_read(),' nLine =',n_line_read()

       if(IsUninitialized)call set_defaults
       IsUninitialized=.false.
    end select

    ! Read input data from text via ModReadParam
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#STRICT")
          call read_var('UseStrict',UseStrict)
       case("#IONODIR")
          call read_var("NameIonoDir",NameIonoDir)
          call fix_dir_name(NameIonoDIr)
          if (iProc==0) call check_dir(NameIonoDir)
       case("#SAVEPLOT", "#IE_SAVEPLOT")
          call read_var('nPlotFile',nFile)
          if (nFile > MaxFile)call CON_stop(NameSub//&
               ' IE_ERROR number of ouput files is too large in #IE_SAVEPLOT:'&
               //' nFile>MaxFile')
          if (nFile>0.and.iProc==0) call check_dir(NameIonoDir)
          do iFile=1,nFile

             call read_var('StringPlot',plot_string)
             call lower_case(plot_string)

             ! Check to see if the ionosphere directory exists...
             if(iProc==0)call check_dir(NameIonoDir)

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
                plot_vars(iFile)='minimum'
             elseif(index(plot_string,'max')>0)then
                plot_vars(iFile)='maximum'
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
       case("#IONOSPHERE")
          call read_var('iConductanceModel',conductance_model)
          call read_var('UseFullCurrent' ,UseFullCurrent)
          call read_var('UseFakeRegion2' ,UseFakeRegion2)
          call read_var('F10.7 Flux',f107_flux)
          call read_var('StarLightPedConductance',StarLightPedConductance)
          call read_var('PolarCapPedConductance',PolarCapPedConductance)

       case("#IM")
          call read_var('TypeImCouple',TypeImCouple)
          call lower_case(TypeImCouple)
       case("#BOUNDARY")
          call read_var('LatBoundary',LatBoundary)
          LatBoundary = LatBoundary * cDegToRad
       case("#UA")
          call read_var('DoCoupleUaCurrent',DoCoupleUaCurrent)
          if(DoCoupleUaCurrent)then
             call read_var('LatBoundary',LatBoundary)
             LatBoundary = LatBoundary * cDegToRad
          endif
       case("#SPS")
          call read_var('UseSPS',UseSPS)
          IE_NameOfEFieldModel = "SPS"
          UseGridBasedIE = .true.
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
          IE_NameOfEFieldModel = "amie"
          UseGridBasedIE = .true.
          UseAMIE = .true.

       case("#BACKGROUND")

          call read_var('NameOfModelDir',IE_NameOfModelDir)
          call read_var('NameOfEFieldModel',IE_NameOfEFieldModel)
          call read_var('NameOfAuroralModel',IE_NameOfAuroralModel)
          call read_var('NameOfSolarModel',IE_NameOfSolarModel)

          if (index(IE_NameOfAuroralModel,'IHP') > 0) &
               IE_NameOfAuroralModel = 'ihp'
          if (index(IE_NameOfAuroralModel,'PEM') > 0) &
               IE_NameOfAuroralModel = 'pem'

          if (index(IE_NameOfEFieldModel,'AMIE') > 0) &
               IE_NameOfEFieldModel = 'amie'

          if (index(IE_NameOfEFieldModel,'weimer01') > 0) &
               IE_NameOfEFieldModel = 'weimer01'
          if (index(IE_NameOfEFieldModel,'Weimer01') > 0) &
               IE_NameOfEFieldModel = 'weimer01'
          if (index(IE_NameOfEFieldModel,'WEIMER01') > 0) &
               IE_NameOfEFieldModel = 'weimer01'

          if (index(IE_NameOfEFieldModel,'weimer') > 0 .and. &
               index(IE_NameOfEFieldModel,'01') == 0) &
               IE_NameOfEFieldModel = 'weimer96'
          if (index(IE_NameOfEFieldModel,'Weimer') > 0 .and. &
               index(IE_NameOfEFieldModel,'01') == 0) &
               IE_NameOfEFieldModel = 'weimer96'
          if (index(IE_NameOfEFieldModel,'WEIMER') > 0 .and. &
               index(IE_NameOfEFieldModel,'01') == 0) &
               IE_NameOfEFieldModel = 'weimer96'

          if (index(IE_NameOfEFieldModel,'weimer96') > 0) &
               IE_NameOfEFieldModel = 'weimer96'
          if (index(IE_NameOfEFieldModel,'Weimer96') > 0) &
               IE_NameOfEFieldModel = 'weimer96'
          if (index(IE_NameOfEFieldModel,'WEIMER96') > 0) &
               IE_NameOfEFieldModel = 'weimer96'

          if (index(IE_NameOfEFieldModel,'SAMIE') > 0) &
               IE_NameOfEFieldModel = 'samie'

          UseGridBasedIE = .false.

       case("#SAVELOGFILE")
          call read_var('DoSaveLogfile',DoSaveLogfile)
          if(DoSaveLogfile)then
             if(iProc==0)call check_dir(NameIonoDir)
          endif

       case("#CONDUCTANCE")
          call read_var('OvalWidthFactor',OvalWidthFactor)
          call read_var('OvalStrengthFactor',OvalStrengthFactor)
          if (conductance_model/=4) then
             write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                  ' command '//trim(NameCommand)// &
                  ' can only be used with conductance model 4'
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if

       case default
          if(iProc==0) then
             write(*,'(a,i4,a)')NameSub//' IE_ERROR at line ',i_line_read(),&
                  ' invalid command '//trim(NameCommand)
             if(UseStrict)call CON_stop('Correct PARAM.in!')
          end if
       end select
    end do

  end subroutine read_param
  !===========================================================================
  subroutine set_defaults

    conductance_model       = 5
    UseFullCurrent          = .false.
    UseFakeRegion2          = .false.
    StarLightPedConductance = 0.25
    PolarCapPedConductance  = 0.25
    f107_flux               = 150.0

  end subroutine set_defaults

end subroutine IE_set_param
!=============================================================================
subroutine IE_set_grid

  ! Set the grid descriptor for IE
  ! Since IE has a static grid the descriptor has to be set once.
  ! There can be many couplers that attempt to set the descriptor,
  ! so we must check IsInitialized.
  use ModProcIE
  use ModIonosphere
  use IE_ModIo
  use IE_ModMain
  use CON_coupler
  use ModNumConst

  implicit none
  character (len=*), parameter :: NameSub='IE_set_grid'
  logical :: IsInitialized=.false.
  integer :: iProc_A(2)

  real :: Colat_I(2*IONO_nTheta-1)

  logical :: DoTest, DoTestMe

  !------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' IsInitialized=',IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  ! IE runs on 1 or 2 PE-s so processor array is (/0,0/) or (/0,1/)
  iProc_A(1)=0
  iProc_A(2)=nProc-1

  if (nProc<0)then
     IONO_NORTH_Theta=0.
     IONO_SOUTH_Theta=0.
     IONO_NORTH_Psi=0.
     IONO_Radius=1.
     IONO_Height=1.
  end if

  ! The colatitudes for both hemispheres
  Colat_I(            1:  IONO_nTheta) = IONO_NORTH_Theta(:,1)
  Colat_I(IONO_nTheta:2*IONO_nTheta-1) = IONO_SOUTH_Theta(:,1)

  call set_grid_descriptor(                        &
       IE_,                          &! component index
       nDim=2,                       &! dimensionality
       nRootBlock_D=(/2,1/),         &! north+south hemispheres
       nCell_D =(/IONO_nTheta - 1,IONO_nPsi - 1/), &! size of node based grid
       XyzMin_D=(/cOne, cOne/),      &! min colat and longitude indexes
       XyzMax_D=(/real(2*IONO_nTheta-1),&
       real(IONO_nPsi)/),            &! max colat and longitude indexes
       TypeCoord='SMG',                            &! solar magnetic coord.
       Coord1_I=Colat_I,                           &! colatitudes
       Coord2_I=IONO_NORTH_Psi(1,:),               &! longitudes
       Coord3_I=(/IONO_Radius + IONO_Height/),     &! radial size in meters
       iProc_A = iProc_A)                           ! processor assigment

end subroutine IE_set_grid

!==============================================================================

subroutine IE_get_for_gm(Buffer_II,iSize,jSize,tSimulation)

  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_gm'

  integer, intent(in)           :: iSize,jSize
  real, intent(out)             :: Buffer_II(iSize,jSize)
  real,             intent(in)  :: tSimulation

  integer :: i,j,k
  real    :: tSimulationTmp
  !--------------------------------------------------------------------------
  if(iSize /= IONO_nTheta*2-1 .or. jSize /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call IE_run(tSimulationTmp,tSimulation)

  Buffer_II = IONO_Phi

end subroutine IE_get_for_gm
!==============================================================================

subroutine IE_get_for_pw(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
     tSimulation)

  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_pw'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
  character (len=*),intent(in)  :: NameHem
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar
  real    :: tSimulationTmp
  !--------------------------------------------------------------------------
  if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call IE_run(tSimulationTmp,tSimulation)

  select case(NameHem)

  case('North')

     if(iProc /= 0) RETURN
     do iVar = 1, nVar
        select case(Name_V(iVar))
        case('Pot')
           Buffer_IIV(:,:,iVar) = IONO_NORTH_Phi
        case('Jr')
           Buffer_IIV(:,:,iVar) = IONO_NORTH_Jr
        case default
           call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
        end select
     end do

  case('South')

     if(iProc /= nProc - 1) RETURN
     do iVar = 1, nVar
        select case(Name_V(iVar))
        case('Pot')
           Buffer_IIV(:,:,iVar) = IONO_SOUTH_Phi
        case('Jr')
           Buffer_IIV(:,:,iVar) = IONO_SOUTH_Jr
        case default
           call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
        end select
     end do

  case default

     call CON_stop(NameSub//' invalid NameHem='//NameHem)

  end select

end subroutine IE_get_for_pw
!==============================================================================

subroutine IE_get_for_rb(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
     tSimulation)

  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_rb'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
  character (len=*),intent(in)  :: NameHem
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar
  real    :: tSimulationTmp
  !--------------------------------------------------------------------------
  if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call IE_run(tSimulationTmp,tSimulation)

  select case(NameHem)

  case('North')

     if(iProc /= 0) RETURN
     do iVar = 1, nVar
        select case(Name_V(iVar))
        case('Pot')
           Buffer_IIV(:,:,iVar) = IONO_NORTH_Phi
        case default
           call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
        end select
     end do

  case('South')

     if(iProc /= nProc - 1) RETURN
     do iVar = 1, nVar
        select case(Name_V(iVar))
        case('Pot')
           Buffer_IIV(:,:,iVar) = IONO_SOUTH_Phi
        case default
           call CON_stop(NameSub//' invalid NameVar='//Name_V(iVar))
        end select
     end do

  case default

     call CON_stop(NameSub//' invalid NameHem='//NameHem)

  end select

end subroutine IE_get_for_rb
!==============================================================================

subroutine IE_get_for_ua(Buffer_II,iSize,jSize,NameVar,NameHem,tSimulation)

  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_ua'

  integer,          intent(in)  :: iSize,jSize
  real,             intent(out) :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
  character (len=*),intent(in)  :: NameHem
  real,             intent(in)  :: tSimulation

  integer :: i,j,k
  real    :: tSimulationTmp
  !--------------------------------------------------------------------------
  if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call IE_run(tSimulationTmp,tSimulation)

  select case(NameHem)

  case('North')

     if(iProc /= 0) RETURN

     select case(NameVar)

     case('Pot')
        Buffer_II = IONO_NORTH_Phi
     case('Ave')
        Buffer_II = IONO_NORTH_Ave_E
     case('Tot')
        Buffer_II = IONO_NORTH_EFlux
     case default
        call CON_stop(NameSub//' invalid NameVar='//NameVar)

     end select

  case('South')

     if(iProc /= nProc - 1) RETURN

     select case(NameVar)

     case('Pot')
        Buffer_II = IONO_SOUTH_Phi
     case('Ave')
        Buffer_II = IONO_SOUTH_Ave_E
     case('Tot')
        Buffer_II = IONO_SOUTH_EFlux
     case default
        call CON_stop(NameSub//' invalid NameVar='//NameVar)

     end select

  case default

     call CON_stop(NameSub//' invalid NameHem='//NameHem)

  end select

end subroutine IE_get_for_ua

!==============================================================================
subroutine IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

  use IE_ModMain, ONLY: IsNewInput, LatBoundaryGm
  use ModProcIE
  use ModIonosphere
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

  if (iProc == 0) then
     LatBoundaryGm = Buffer_IIV(IONO_nTheta,1,1)
     if(DoTest)write(*,*) "LatBoundary : ",LatBoundaryGm*180.0/3.1415926
     Iono_North_Jr = Buffer_IIV(1:IONO_nTheta,:,1)
     Iono_North_Jr(IONO_nTheta-1:IONO_nTheta,1) = 0.0
     if(nVar>1)then
        Iono_North_invB           = Buffer_IIV(1:IONO_nTheta,:,2)
        Iono_North_rho            = Buffer_IIV(1:IONO_nTheta,:,3)
        Iono_North_p              = Buffer_IIV(1:IONO_nTheta,:,4)
        if(DoTest) call write_dataN
     end if
  endif
  if (iProc == nProc-1) then
     LatBoundaryGm = Buffer_IIV(IONO_nTheta,1,1)
     if(DoTest)write(*,*) "LatBoundary2 : ",LatBoundaryGm*180.0/3.1415926
     Iono_South_Jr = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,1)
     Iono_South_Jr(1:2,1) = 0.0
     if(nVar>1)then
        Iono_South_invB           = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,2)
        Iono_South_rho            = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,3)
        Iono_South_p              = Buffer_IIV(IONO_nTheta:2*IONO_nTheta-1,:,4)
        if(DoTest) call write_dataS
     end if
  endif

  if(DoTest)write(*,*)NameSub,' finished'

  contains

    !=========================================================================
    !write values to North plot file
    subroutine write_dataN
      use ModIoUnit, ONLY: UNITTMP_
      CHARACTER (LEN=80) :: filename
      integer :: i,j
      integer, save :: nCall=0
      !------------------------------------------------------------------------

      nCall=nCall+1
      write(filename,'(a,i5.5,a)')"gm2ie_debugN_",nCall,".dat"
      OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
      write(UNITTMP_,'(a)') 'TITLE="gm2ie debugN values"'
      write(UNITTMP_,'(a)') &
           'VARIABLES="J", "I", "Theta", "Psi", "JR", "1/B", "rho", "p"'
      write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
           'ZONE T="PE=',iProc,'", I=',jsize,', J=',isize,', K=1, F=POINT'
      do i=1,isize; do j=1,jsize
         write(UNITTMP_,'(2i4,6G14.6)') j,i, &
              Iono_North_Theta(i,j),Iono_North_Psi(i,j),Iono_North_Jr(i,j), &
              Iono_North_invB(i,j),Iono_North_rho(i,j),Iono_North_p(i,j)
      end do; end do
      CLOSE(UNITTMP_)
    end subroutine write_dataN

    !write values to South plot file
    subroutine write_dataS
      use ModIoUnit, ONLY: UNITTMP_
      CHARACTER (LEN=80) :: filename
      integer :: i,j
      integer, save :: nCall=0
      !-------------------------------------------------------------------------
      
      nCall=nCall+1
      write(filename,'(a,i5.5,a)')"gm2ie_debugS_",nCall,".dat"
      OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
      write(UNITTMP_,'(a)') 'TITLE="gm2ie debugS values"'
      write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Theta", "Psi", "JR", "1/B", "rho", "p"'
      write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
           'ZONE T="PE=',iProc,'", I=',jsize,', J=',isize,', K=1, F=POINT'
      do i=1,isize; do j=1,jsize
         write(UNITTMP_,'(2i4,6G14.6)') j,i, &
              Iono_South_Theta(i,j),Iono_South_Psi(i,j),Iono_South_Jr(i,j), &
              Iono_South_invB(i,j),Iono_South_rho(i,j),Iono_South_p(i,j)
      end do; end do
      CLOSE(UNITTMP_)
    end subroutine write_dataS

end subroutine IE_put_from_gm

!==============================================================================

subroutine IE_put_from_ua(Buffer_III, iBlock, nMLTs, nLats, nVarsToPass)

  use IE_ModMain, ONLY: IsNewInput, DoCoupleUaCurrent, StarLightPedConductance
  use ModIonosphere
  use ModConst
  use ModUtilities, ONLY: check_allocate

  implicit none

  save

  integer, intent(in) :: nMlts, nLats, iBlock, nVarsToPass
  real, dimension(nMlts, nLats, nVarsToPass), intent(in) :: Buffer_III

  !\
  ! UA_Lats and UA_Mlts are the latitudes and magnetic local times of the
  !    UA magnetic grid.
  ! iLat and iMlt are the indices of where to find the points in the
  !    UA magnetic grid.
  ! rLat and rMlt are the multiplication factors to get stuff from the
  !    UA magnetic grid.
  !/

  real,    dimension(:,:,:), allocatable :: UA_Lats, UA_Mlts
  integer, dimension(Iono_nTheta,2) :: iLat
  integer, dimension(Iono_nPsi,2)   :: iMlt
  real, dimension(Iono_nTheta,2)    :: rLat
  real, dimension(Iono_nPsi,2)      :: rMlt

  integer :: iError, i, j, ii, jj
  real    :: t, p

  integer, parameter :: Fac_ = 1
  integer, parameter :: Ped_ = 2
  integer, parameter :: Hal_ = 3
  integer, parameter :: Lat_ = 4
  integer, parameter :: Mlt_ = 5

  character(len=*), parameter :: NameSub='IE_put_from_ua'

  !--------------------------------------------------------------------------

  IsNewInput=.true.

  if (nVarsToPass == 5) then
     if (.not.allocated(UA_Lats)) then
        allocate(UA_Lats(nMLTs, nLats,2), &
             UA_Mlts(nMLTs, nLats,2),     &
             stat=iError)
        call check_allocate(iError,NameSub//'UA_Lats,UA_Mlts')
     endif
     UA_Lats(:,:,iBlock) = Buffer_III(:,:,Lat_)
     UA_Mlts(:,:,iBlock) = Buffer_III(:,:,Mlt_)

     !\
     ! In this instance, t = theta
     !/

     do i = 1, IONO_nTheta
        if (iBlock == 1) t = 90.0 - Iono_North_Theta(i,1)*cRadToDeg
        if (iBlock == 2) t = 90.0 - Iono_South_Theta(i,1)*cRadToDeg

        if (t >= maxval(UA_Lats(1,:,iBlock))) then
           ii = nLats-1
           iLat(i,iBlock) = ii
           rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
        else if (t < minval(UA_Lats(1,:,iBlock))) then
           ii = 1
           iLat(i,iBlock) = ii
           rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
        else
           ii = 1
           do while (ii < nLats)
              if ((t >= UA_Lats(1,ii,iBlock) .and. &
                   t <  UA_Lats(1,ii+1,iBlock)) .or.  &
                   (t <= UA_Lats(1,ii,iBlock) .and. &
                   t >  UA_Lats(1,ii+1,iBlock))) then
                 iLat(i,iBlock) = ii
                 rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                      (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
                 ii = nLats
              endif
              ii = ii+1
           enddo

        endif

     enddo

     !\
     ! In this instance, p = psi
     !/

     do j = 1, IONO_nPsi
        if (iBlock == 1) p = mod(Iono_North_Psi(1,j)*12.0/cPi + 12.0,24.0)
        if (iBlock == 2) p = mod(Iono_South_Psi(1,j)*12.0/cPi + 12.0,24.0)

        jj = 1
        do while (jj < nMlts)
           if ((p >= UA_Mlts(jj,1,iBlock) .and. &
                p <  UA_Mlts(jj+1,1,iBlock)) .or.  &
                (p <= UA_Mlts(jj,1,iBlock) .and. &
                p >  UA_Mlts(jj+1,1,iBlock))) then
              iMlt(j,iBlock) = jj
              rMlt(j,iBlock) = 1.0 - (p - UA_Mlts(jj,1,iBlock)) / &
                   (UA_Mlts(jj+1,1,iBlock) - UA_Mlts(jj,1,iBlock))
              jj = nMlts
           end if
           jj = jj+1
        enddo

     enddo

     !write(*,*)NameSub,' iBlock=',iBlock,&
     !     ' iLat(1:5)=',iLat(1:5,iBlock), &
     !     ' rLat(1:5)=',rLat(1:5,iBlock)
     !
     !write(*,*)NameSub,' iBlock=',iBlock,&
     !     ' iLat(Iono_nTheta-5:Iono_nTheta)=',&
     !     iLat(Iono_nTheta-5:Iono_nTheta,iBlock), &
     !     ' rLat(Iono_nTheta-5:Iono_nTheta)=',&
     !     rLat(Iono_nTheta-5:Iono_nTheta,iBlock)

  end if


  if (iBlock == 1) then ! iBlock == 1: Northern Hemisphere

     do i = 1, Iono_nTheta

        !\
        ! Now t = 0.0 - 1.0, and is the interpolation coefficient for theta
        !/

        ii = iLat(i,iBlock)
        t  = rLat(i,iBlock)

        do j = 1, Iono_nPsi

           !\
           ! Now p = 0.0 - 1.0, and is the interpolation coefficient for psi
           !/

           jj = iMlt(j,iBlock)
           p  = rMlt(j,iBlock)

           IONO_NORTH_SigmaH(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Hal_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Hal_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Hal_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Hal_)

           IONO_NORTH_SigmaP(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Ped_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Ped_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Ped_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Ped_)

           if(DoCoupleUaCurrent)then
              IONO_NORTH_TGCM_JR(i,j) =                          &
                   (    t)*(    p)*Buffer_III(jj  ,ii  ,Fac_) + &
                   (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Fac_) + &
                   (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Fac_) + &
                   (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Fac_)
           else
              IONO_NORTH_TGCM_JR(i,j) = 0.0
           end if

        enddo
     enddo
     ! Limit the conductance with the StarLightConductance
     IONO_NORTH_SigmaP = max(IONO_NORTH_SigmaP,   StarLightPedConductance)
     IONO_NORTH_SigmaH = max(IONO_NORTH_SigmaH, 2*StarLightPedConductance)

  else ! iBlock == 2: Southern Hemisphere

     do i = 1, Iono_nTheta

        !\
        ! Now t = 0.0 - 1.0, and is the interpolation coefficient for theta
        !/

        ii = iLat(i,iBlock)
        t  = rLat(i,iBlock)

        do j = 1, Iono_nPsi

           !\
           ! Now p = 0.0 - 1.0, and is the interpolation coefficient for psi
           !/

           jj = iMlt(j,iBlock)
           p  = rMlt(j,iBlock)

           IONO_SOUTH_SigmaH(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Hal_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Hal_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Hal_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Hal_)

           IONO_SOUTH_SigmaP(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Ped_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Ped_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Ped_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Ped_)

           if(DoCoupleUaCurrent)then
              IONO_SOUTH_TGCM_JR(i,j) =                          &
                   (    t)*(    p)*Buffer_III(jj  ,ii  ,Fac_) + &
                   (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Fac_) + &
                   (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Fac_) + &
                   (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Fac_)
           else
              IONO_SOUTH_TGCM_JR(i,j) = 0.0
           end if
        enddo
     enddo
     ! Limit the conductance with the StarLightConductance
     IONO_SOUTH_SigmaP = max(IONO_SOUTH_SigmaP,   StarLightPedConductance)
     IONO_SOUTH_SigmaH = max(IONO_SOUTH_SigmaH, 2*StarLightPedConductance)

  end if

end subroutine IE_put_from_ua

!==============================================================================

subroutine IE_put_from_im(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  use IE_ModMain, only:IsNewInput
  use ModIonosphere
  use ModProcIE

  use CON_router,   ONLY: IndexPtrType, WeightPtrType

  implicit none
  character(len=*), parameter   :: NameSub='IE_put_from_im'
  integer,intent(in)            :: nPoint, iPointStart, nVar
  real, intent(in)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight
  logical,intent(in)            :: DoAdd
  integer :: iBlock,iLat,iLon
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

  if ( iLat >= 1 .and. iLat <= iono_nTheta .and. &
       iLon >= 1 .and. iLon <= iono_nPsi) then
     if (.not.IsFilledWithIm(iLat,iLon)) then
        iono_north_im_jr(iLat,iLon)    = buff_v(1)
        iono_north_im_eflux(iLat,iLon) = buff_v(2)
        iono_north_im_avee(iLat,iLon)  = buff_v(3)
        IsFilledWithIm(iLat,iLon) = .true.
     endif
  endif

  IsNewInput = .true.

end subroutine IE_put_from_im

!==============================================================================

subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide potential and current for IM
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi, &
       IONO_NORTH_PHI, IONO_NORTH_JR, IONO_SOUTH_PHI, IONO_SOUTH_JR, &
       IONO_NORTH_SigmaH, IONO_NORTH_SigmaP, &
       IONO_SOUTH_SigmaH, IONO_SOUTH_SigmaP, &
       cpcp_north, cpcp_south
  use IE_ModMain,    ONLY: TypeImCouple

  implicit none
  character(len=*), parameter :: NameSub='IE_get_for_im'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iBlock, i, j, iSouth, iPoint
  real    :: w

  !---------------------------------------------------------------------------
  Buff_V = 0.0

  do iPoint = iPointStart, iPointStart + nPoint - 1

     i      = Index % iCB_II(1,iPoint)
     j      = Index % iCB_II(2,iPoint)
     iBlock = Index % iCB_II(3,iPoint)
     w      = Weight % Weight_I(iPoint)

     if(iBlock/=1)then
        write(*,*)NameSub,': iPoint,Index % iCB_II=',&
             iPoint,Index%iCB_II(:,iPoint)
        call CON_stop(NameSub//&
             ' SWMF_ERROR iBlock should be 1=North in IE-IM coupling')
     end if

     if(i<1 .or. i>IONO_nTheta .or. j<1 .or. j>IONO_nPsi)then
        write(*,*)'i,j=',i,j
        call CON_stop(NameSub//' SWMF_ERROR index out of range')
     end if

     ! Index for the same latitude on the southern hemisphere
     iSouth = IONO_nTheta + 1 - i

     select case(TypeImCouple)
     case('north')
        Buff_V(1) = Buff_V(1) + w * IONO_NORTH_PHI(i,j)
        Buff_V(2) = Buff_V(2) + w * IONO_NORTH_JR(i,j)
        Buff_V(3) = Buff_V(3) + w * IONO_NORTH_SigmaH(i,j)
        Buff_V(4) = Buff_V(4) + w * IONO_NORTH_SigmaP(i,j)
     case('south')
        Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
        Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
        Buff_V(3) = Buff_V(3) + w * IONO_SOUTH_SigmaH(iSouth,j)
        Buff_V(4) = Buff_V(4) + w * IONO_SOUTH_SigmaP(iSouth,j)
     case('cpcpmin')
        if(cpcp_north < cpcp_south)then
           Buff_V(1) = Buff_V(1) + w * IONO_NORTH_PHI(i,j)
           Buff_V(2) = Buff_V(2) + w * IONO_NORTH_JR(i,j)
           Buff_V(3) = Buff_V(3) + w * IONO_NORTH_SigmaH(i,j)
           Buff_V(4) = Buff_V(4) + w * IONO_NORTH_SigmaP(i,j)
        else
           Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
           Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
           Buff_V(3) = Buff_V(3) + w * IONO_SOUTH_SigmaH(iSouth,j)
           Buff_V(4) = Buff_V(4) + w * IONO_SOUTH_SigmaP(iSouth,j)
        end if
     case('average')
        Buff_V(1) = Buff_V(1) + w * &
             0.5*(IONO_NORTH_PHI(i,j) + IONO_SOUTH_PHI(iSouth,j))
        Buff_V(2) = Buff_V(2) + w * &
             0.5*(IONO_NORTH_JR(i,j)  + IONO_SOUTH_JR(iSouth,j))
        Buff_V(3) = Buff_V(3) + w * 0.5*( &
             IONO_NORTH_SigmaH(i,j)  + &
             IONO_SOUTH_SigmaH(iSouth,j))
        Buff_V(4) = Buff_V(4) + w * 0.5*( &
             IONO_NORTH_SigmaP(i,j)  + &
             IONO_SOUTH_SigmaP(iSouth,j))
     case default
        call CON_stop(NameSub//' ERROR: Unknown value for TypeImCouple='// &
             TypeImCouple)
     end select
  end do

contains

  real function minmod(a,b)
    real, intent(in) :: a,b
    minmod = (sign(0.5, a) + sign(0.5, b)) * min(abs(a), abs(b))
  end function minmod

end subroutine IE_get_for_im

!==============================================================================

subroutine IE_put_from_im_complete

  use ModProcIE
  use ModIonosphere
  use ModMpi

  implicit none

  integer iError, i 

  iono_north_im_eflux(:,iono_npsi) = iono_north_im_eflux(:,1)
  iono_north_im_avee(:,iono_npsi)  = iono_north_im_avee(:,1)
  iono_north_im_jr(:,iono_npsi)  = iono_north_im_jr(:,1)

  if (nProc > 1) then
     iError = 0
     call MPI_Bcast(iono_north_im_eflux, iono_nTheta*iono_nPsi, &
          MPI_Real, 0, iComm, iError)
     call MPI_Bcast(iono_north_im_avee, iono_nTheta*iono_nPsi, &
          MPI_Real, 0, iComm, iError)
     call MPI_Bcast(iono_north_im_jr, iono_nTheta*iono_nPsi, &
          MPI_Real, 0, iComm, iError)
  endif

  do i = 1, IONO_nTheta
     iono_south_im_jr(i,:) = iono_north_im_jr(Iono_nTheta-i+1,:)
  enddo

  IsFilledWithIm = .false.

end subroutine IE_put_from_im_complete

!==============================================================================

subroutine IE_init_session(iSession, tSimulation)

  ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

  use CON_physics,   ONLY: get_time, get_planet, get_axes
  use ModIonosphere, ONLY: IONO_Bdp
  use IE_ModMain,    ONLY: time_accurate, time_simulation, ThetaTilt
  use IE_ModIo,      ONLY: dt_output, t_output_last
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
     call ionosphere_fine_grid
     call ionosphere_init
     IsUninitialized = .false.
  end if

  time_simulation = tSimulation

  call get_time(  DoTimeAccurateOut = time_accurate)
  call get_planet(DipoleStrengthOut = IONO_Bdp)
  call get_axes(tSimulation, MagAxisTiltGsmOut = ThetaTilt)

  IONO_Bdp = IONO_Bdp*1.0e9 ! Tesla -> nT

  if(DoTest)write(*,*)NameSub,': IONO_Bdp, ThetaTilt =',IONO_Bdp,ThetaTilt

  ! Reset t_output_last in case the plotting frequency has changed
  if(time_accurate)then
     where(dt_output>0.) &
          t_output_last=int(time_simulation/dt_output)
  end if

end subroutine IE_init_session

!==============================================================================
subroutine IE_finalize(tSimulation)

  use ModProcIE
  use IE_ModMain, ONLY: Time_Array, time_simulation, nSolve
  use IE_ModIo, ONLY: nFile
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
  call time_real_to_int(tCurrent, Time_Array)
  time_simulation = tSimulation

  if(nSolve>0)then
     do iFile=1,nFile
        if(iProc==0)      call ionosphere_write_output(iFile, 1)
        if(iProc==nProc-1)call ionosphere_write_output(iFile, 2)
     end do
  end if

end subroutine IE_finalize

!==============================================================================

subroutine IE_save_restart(tSimulation)

  use ModProcIE, ONLY:  nProc
  use IE_ModMain, ONLY: nSolve

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: tSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IE_save_restart'

  RETURN
  ! IE does not need a restart file at least not in the framework
  ! ionosphere_write_restart_file is still not parallel

  if(nProc==1)call ionosphere_write_restart_file(nSolve)

end subroutine IE_save_restart

!==============================================================================

subroutine IE_run(tSimulation,tSimulationLimit)

  use ModProcIE
  use IE_ModMain
  use CON_physics, ONLY: get_time, get_axes, time_real_to_int
  use ModKind
  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: tSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: tSimulationLimit ! simulation time not to be exceeded

  real(Real8_) :: tStart
  integer      :: nStep

  character(len=*), parameter :: NameSub='IE_run'

  logical :: DoTest,DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  if(DoTest)write(*,*)NameSub,': iProc,tSimulation,tSimulationLimit=',&
       iProc,tSimulation,tSimulationLimit

  ! Store the current time
  time_simulation = tSimulation

  ! Since IE is not a time dependent component, it may advance to the 
  ! next coupling time in a time accurate run
  if(time_accurate)tSimulation = tSimulationLimit

  if(DoTest)write(*,*)NameSub,': iProc,IsNewInput=',iProc,IsNewInput

  ! Do not solve if there is no new input from GM or UA
  if(.not.IsNewInput) RETURN

  ! Check if we can have a reasonable magnetic field already
  call get_time(nStepOut=nStep)

  if(DoTest)write(*,*)NameSub,': iProc,nStep = ',iProc,nStep

  ! After the solve this input can be considered old
  IsNewInput = .false.

  ! Obtain the position of the magnetix axis
  call get_axes(time_simulation,MagAxisTiltGsmOut = ThetaTilt)

  call get_time(tStartOut=tStart)
  call time_real_to_int(tStart + time_simulation, Time_Array)

  nSolve = nSolve + 1

  if(DoTest)write(*,*) 'solve'

  ! Solve for the ionosphere potential
  call IE_solve

  if(DoTest)write(*,*) 'done with solve'

  ! Save solution (plot files) into file if required
  call IE_output

  if(DoTest)write(*,*) 'done with output'

  call IE_gather

  if(DoTest)write(*,*) 'gather done'

  ! Save logfile if required
  call IE_save_logfile

  if(DoTest)write(*,*) 'done with IE_run'

end subroutine IE_run

!=================================================================
subroutine IE_get_for_ps(Buffer_IIV, iSize, jSize, nVar)

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_ps'

  integer, intent(in)           :: iSize, jSize, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)

  !NOTE: The Buffer variables must be collected to i_proc0(IE_) before return.

  write(*,*) NameSub,' -- called but not yet implemented.'

end subroutine IE_get_for_ps
!==============================================================================

subroutine IE_setnMlts(iComponent, nMLTsIn, iError)

  implicit none

  integer, intent(in)  :: iComponent, nMLTsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setnMlts'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setnMlts

!==============================================================================

subroutine IE_setnLats(iComponent, nLatsIn, iError)

  implicit none

  integer, intent(in)  :: iComponent, nLatsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setnLats'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setnLats

!==============================================================================

subroutine IE_setgrid(iComponent, MLTsIn, LatsIn, iError)

  integer, intent(in) :: iComponent
  real, intent(in) :: MLTsIn,LatsIn
  integer, intent(out) :: iError

  character (len=*), parameter :: NameSub='IE_setgrid'

  call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

end subroutine IE_setgrid

!==============================================================================

