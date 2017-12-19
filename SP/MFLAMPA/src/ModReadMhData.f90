!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModReadMhData

  ! This module contains methods for reading input MH data

  use SP_ModSize, ONLY: nDim, nParticleMax

  use SP_ModGrid, ONLY: get_node_indexes, &
       nVarRead, nVar, nBlock, iShock_IB, iNode_B, FootPoint_VB, &
       nParticle_B,  State_VIB, LagrID_, Z_, Shock_, NameVar_V

  use SP_ModAdvance, ONLY: TimeGlobal, iIterGlobal

  use ModPlotFile, ONLY: read_plot_file

  implicit none

  SAVE

  private ! except
 
  public:: &
       set_read_mh_data_param, init_read_mh_data, read_mh_data, &
       DoReadMhData


  !\
  !----------------------------------------------------------------------------
  ! Format of output files
  integer, parameter:: &
       Tec_ = 0, &
       Idl_ = 1
  !----------------------------------------------------------------------------
  ! the input directory
  character (len=100):: NameInputDir=""
  ! the input file name base
  character (len=100):: NameFileBase="MH_data"
  character (len=4)  :: NameFormat
  character (len=20) :: TypeFile

  ! buffer is larger than the data needed to be read in the case 
  ! the input file has additional data
  real:: Buffer_II(nVar,nParticleMax)
  !
  ! buffer for Lagrangian coordinate
  real:: Buffer_I(nParticleMax)

  real:: TimeRead, TimeReadStart, TimeReadMax, DtRead
  integer:: iIterRead, iIterReadStart, DnRead

  logical:: DoReadMhData = .false.
  !/
contains
  
  subroutine set_read_mh_data_param
    use ModReadParam, ONLY: read_var
    ! set parameters of output files: file format, kind of output etc.
    character(len=300):: StringPlot
    ! loop variables
    integer:: iFile, iNode
    character(len=*), parameter :: NameSub='SP:set_read_mh_data_param'
    !--------------------------------------------------------------------------
    !
    call read_var('DoReadMhData', DoReadMhData)
    if(.not. DoReadMhData)&
         RETURN
    ! the input directory
    call read_var('NameInputDir', NameInputDir)
    ! ADD "/" IF NOT PRESENT
    if(NameInputDir(len_trim(NameInputDir):len_trim(NameInputDir))/='/')&
         NameInputDir = trim(NameInputDir) // '/'
    !
    call read_var('TypeFile', TypeFile)

    ! time step
    call read_var('DtRead', DtRead)
    call read_var('DnRead', DnRead)

    ! the format of output file must be set
    select case(trim(TypeFile))
    case('tec')
       NameFormat='.dat'
    case('idl','ascii')
       TypeFile = 'ascii'
       NameFormat='.out'
    case('real4','real8')
       NameFormat='.out'
    case default
       call CON_stop(NameSub//': input format was not set in PARAM.in')
    end select
  end subroutine set_read_mh_data_param

  !============================================================================

  subroutine init_read_mh_data
    ! initialize by setting the time and interation index of input files
    !-------------------------------------------------------------------------
    TimeRead = TimeGlobal
    iIterRead= iIterGlobal
  end subroutine init_read_mh_data

  !============================================================================

  subroutine read_mh_data(TimeOut)
    real,    intent(out):: TimeOut
    ! read 1D MH data, which are produced by write_mh_1d n ModWrite
    ! separate file is read for each field line, name format is (usually)
    ! MH_data_<iLon>_<iLat>_t<ddhhmmss>_n<iIter>.{out/dat}
    !------------------------------------------------------------------------
    ! name of the input file
    character(len=100):: NameFile
    ! loop variables
    integer:: iBlock, iParticle, iVarPlot
    ! indexes of corresponding node, latitude and longitude
    integer:: iNode, iLat, iLon
    ! number of particles saved in the input file
    integer:: nParticleInput
    ! auxilary parameter index
    integer, parameter:: RShock_ = Z_ + 2
    ! additional parameters of lines
    real:: Param_I(LagrID_:RShock_)
    ! timestamp
    character(len=8):: StringTime
    !------------------------------------------------------------------------

    ! the simulation time corresponding to the input file
    TimeOut = TimeRead

    ! read the data
    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       call get_time_string(TimeRead, StringTime)
       write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
            trim(NameInputDir)//trim(NameFileBase)//'_',iLon,'_',iLat,&
            '_t'//StringTime//'_n',iIterRead, NameFormat

       ! read the header first
       call read_plot_file(&
            NameFile   = NameFile,&
            TypeFileIn = TypeFile,&
            n1out      = nParticleInput,&
            Coord1Out_I= Buffer_I,&
            VarOut_VI  = Buffer_II,&
            ParamOut_I = Param_I(LagrID_:RShock_)&
            )
       State_VIB(LagrID_   , 1:nParticleInput, iBlock) = &
            Buffer_I(             1:nParticleInput)
       State_VIB(1:nVarRead, 1:nParticleInput, iBlock) = &
            Buffer_II(1:nVarRead, 1:nParticleInput)

       nParticle_B(  iBlock) = nParticleInput

       !Parameters
       FootPoint_VB(LagrID_:Z_,iBlock) = Param_I(LagrID_:Z_)
       ! shock location
       iShock_IB(Shock_,iBlock) = nint(Param_I(RShock_-1))

    end do

    ! advance read time and iteration
    TimeRead  = TimeRead  + DtRead
    iIterRead = iIterRead + DnRead

  end subroutine read_mh_data

  !==========================================================================

  subroutine get_time_string(Time, StringTime)
    ! the subroutine converts real variable Time into a string,
    ! the structure of the string is 'ddhhmmss', 
    ! i.e shows number of days, hours, minutes and seconds 
    ! after the beginning of the simulation
    real,             intent(in) :: Time
    character(len=8), intent(out):: StringTime
    !--------------------------------------------------------------------------
    ! This is the value if the time is too large
    StringTime = '99999999'
    if(Time < 8640000) &
         write(StringTime,'(i2.2,i2.2,i2.2,i2.2)') &
         int(                  Time          /86400), & ! # days
         int((Time-(86400*int(Time/86400)))/ 3600), & ! # hours
         int((Time-( 3600*int(Time/ 3600)))/   60), & ! # minutes
         int( Time-(  60*int(Time/   60)))           ! # seconds
  end subroutine get_time_string

end module SP_ModReadMhData
