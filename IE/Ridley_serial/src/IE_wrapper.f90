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
          call read_var('nFile',nFile)
          if (nFile > MaxFile)call CON_stop(NameSub//&
               ' IE_ERROR number of ouput files is too large in #IE_SAVEPLOT:'&
               //' nFile>MaxFile')
          if (nFile>0.and.iProc==0) call check_dir(NameIonoDir)
          do iFile=1,nFile

             call read_var('plot_string',plot_string)
             call lower_case(plot_string)

             ! Check to see if the ionosphere directory exists...
             if(iProc==0)call check_dir(NameIonoDir)

             ! Plotting frequency
             call read_var('dn_output',dn_output(iFile))
             call read_var('dt_output',dt_output(iFile))

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
          call read_var('conductance_model',conductance_model)
          call read_var('UseFullCurrent' ,UseFullCurrent)
          call read_var('UseFakeRegion2' ,UseFakeRegion2)
          call read_var('F10.7 Flux',f107_flux)
          call read_var('StarLightPedConductance',StarLightPedConductance)
          call read_var('PolarCapPedConductance',PolarCapPedConductance)

       case("#IM")
          call read_var('TypeImCouple',TypeImCouple)
          call lower_case(TypeImCouple)
       case("#SPS")
          call read_var('UseSPS',UseSPS)
          IE_NameOfEFieldModel = "SPS"
          UseGridBasedIE = .true.

       case("#DEBUG")
          call read_var('iDebugLevel',iDebugLevel)
          call read_var('iDebugProc',iDebugProc)
          if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
             iDebugLevel = -1
          endif

       case("#AMIEFILES")
          call read_var('AMIEFileNorth',AMIEFileNorth)
          call read_var('AMIEFileSouth',AMIEFileSouth)
          IE_NameOfEFieldModel = "amie"
          UseGridBasedIE = .true.
          UseAMIE = .true.

       case("#BACKGROUND")

          call read_var('IE_NameOfModelDir',IE_NameOfModelDir)
          call read_var('IE_NameOfEFieldModel',IE_NameOfEFieldModel)
          call read_var('IE_NameOfAuroralModel',IE_NameOfAuroralModel)
          call read_var('IE_NameOfSolarModel',IE_NameOfSolarModel)

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

  ! In the model Ridley_serial the grid is not exactly uniform and not
  ! exactly periodic.

end subroutine IE_set_grid

!==============================================================================

subroutine IE_get_for_gm(Buffer_II,iSize,jSize,NameVar,tSimulation)

  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*),parameter :: NameSub='IE_get_for_gm'

  integer, intent(in)           :: iSize,jSize
  real, intent(out)             :: Buffer_II(iSize,jSize)
  character (len=*),intent(in)  :: NameVar
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
  call IE_run(tSimulationTmp,tSimulationTmp)

  select case(NameVar)
  case('PotNorth')
     if(iProc /= 0) RETURN
     Buffer_II = IONO_NORTH_Phi
  case('PotSouth')
     if(iProc /= nProc - 1) RETURN
     Buffer_II= IONO_SOUTH_Phi
  case default
     call CON_stop(NameSub//' SWMF_ERROR invalid NameVar='//NameVar)
  end select

end subroutine IE_get_for_gm
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
  call IE_run(tSimulationTmp,tSimulationTmp)

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
subroutine IE_put_from_gm(Buffer_II, iSize, jSize, NameVar)

  use IE_ModMain, ONLY: IsNewInput
  use ModProcIE
  use ModIonosphere

  implicit none
  character (len=*), parameter :: NameSub = 'IE_put_from_gm'
  integer,          intent(in) :: iSize, jSize
  real,             intent(in) :: Buffer_II(iSize, jSize)
  character(len=*), intent(in) :: NameVar

  integer :: iError
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,' starting with NameVar=',NameVar

  IsNewInput = .true.

  select case(NameVar)

  case('JrNorth')
     if (iProc /= 0) RETURN
     Iono_North_Jr         = Buffer_II
  case('JrSouth')
     if (iProc /= nProc-1) RETURN
     Iono_South_Jr         = Buffer_II
  case default
     call CON_stop(NameSub//' SWMF_ERROR invalid NameVar='//NameVar)
  end select

  if(DoTest)write(*,*)NameSub,' finished'

end subroutine IE_put_from_gm

!==============================================================================

subroutine IE_put_from_UA(Buffer_III, iBlock, nMLTs, nLats, nVarsToPass)

  use IE_ModMain, ONLY: IsNewInput
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

  character(len=*), parameter :: NameSub='IE_put_for_UA'

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
        if (iBlock == 1) t = 90.0 - Iono_North_Theta(i,1)*180.0/cPi
        if (iBlock == 2) t = 90.0 - Iono_South_Theta(i,1)*180.0/cPi

        if (t > maxval(UA_Lats(1,:,iBlock))) then
           ii = nLats-1
           iLat(i,iBlock) = ii
           rLat(i,iBlock) = 1.0 - (t - UA_Lats(1,ii,iBlock)) / &
                (UA_Lats(1,ii+1,iBlock) - UA_Lats(1,ii,iBlock))
        else
           if (t < minval(UA_Lats(1,:,iBlock))) then
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

  end if

  if (iBlock == 1) then

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

           IONO_NORTH_TGCM_JR(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Fac_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Fac_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Fac_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Fac_)

        enddo
     enddo

  else

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

           IONO_SOUTH_TGCM_JR(i,j) =                          &
                (    t)*(    p)*Buffer_III(jj  ,ii  ,Fac_) + &
                (1.0-t)*(    p)*Buffer_III(jj  ,ii+1,Fac_) + &
                (    t)*(1.0-p)*Buffer_III(jj+1,ii  ,Fac_) + &
                (1.0-t)*(1.0-p)*Buffer_III(jj+1,ii+1,Fac_)

        enddo
     enddo

  end if

end subroutine IE_put_from_UA

!==============================================================================

!==============================================================================

subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide potential and current for IM
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V

  use CON_coupler,   ONLY: IndexPtrType, WeightPtrType
  use ModIonosphere, ONLY: IONO_nTheta, IONO_nPsi, &
       IONO_NORTH_PHI, IONO_NORTH_JR, IONO_SOUTH_PHI, IONO_SOUTH_JR, &
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
     case('south')
        Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
        Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
     case('cpcpmin')
        if(cpcp_north < cpcp_south)then
           Buff_V(1) = Buff_V(1) + w * IONO_NORTH_PHI(i,j)
           Buff_V(2) = Buff_V(2) + w * IONO_NORTH_JR(i,j)
        else
           Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_PHI(iSouth,j)
           Buff_V(2) = Buff_V(2) + w * IONO_SOUTH_JR(iSouth,j)
        end if
     case('average')
        Buff_V(1) = Buff_V(1) + w * &
             0.5*(IONO_NORTH_PHI(i,j) + IONO_SOUTH_PHI(iSouth,j))
        Buff_V(2) = Buff_V(2) + w * &
             0.5*(IONO_NORTH_JR(i,j)  + IONO_SOUTH_JR(iSouth,j))
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

subroutine IE_init_session(iSession, tSimulation)

  ! Initialize the Ionosphere Electrostatic (IE) module for session iSession

  use CON_physics, ONLY: get_time, get_planet, get_axes
  use ModIonosphere, ONLY: IONO_Bdp
  use IE_ModMain
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

  call get_time(  DoTimeAccurateOut = time_accurate)
  call get_planet(DipoleStrengthOut = IONO_Bdp)
  call get_axes(tSimulation, MagAxisTiltGsmOut = ThetaTilt)

  IONO_Bdp = IONO_Bdp*1.0e9 ! Tesla -> nT

  write(*,*)NameSub,': IONO_Bdp, ThetaTilt =',IONO_Bdp,ThetaTilt

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

  ! Solve for the ionosphere potential
  call IE_solve

  ! Save solution (plot files) into file if required
  call IE_output

end subroutine IE_run

!=================================================================
