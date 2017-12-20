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

  use ModUtilities, ONLY: fix_dir_name

  implicit none

  SAVE

  private ! except
 
  public:: &
       set_read_mh_data_param, init_read_mh_data, read_mh_data, &
       DoReadMhData


  !\
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

  ! buffer for Lagrangian coordinate
  real:: Buffer_I(nParticleMax)

  ! number of input files
  integer:: nFileRead = 0
  ! index of a current input file
  integer:: iFileRead

  ! time/iteration stamps of input files
  character(len=17), allocatable:: NameFileStamp_I(:)

  logical:: DoReadMhData = .false.
  !/
contains
  
  subroutine set_read_mh_data_param(NameCommand)
    use ModReadParam, ONLY: read_var
    ! set parameters of input files with background data
    character (len=*), intent(in):: NameCommand ! From PARAM.in  
    ! loop variables
    integer:: iFile
    character(len=*), parameter :: NameSub='SP:set_read_mh_data_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#READMHDATA')
       !
       call read_var('DoReadMhData', DoReadMhData)
       if(.not. DoReadMhData)&
            RETURN
       ! the input directory
       call read_var('NameInputDir', NameInputDir)
       call fix_dir_name(NameInputDir) ! adds "/" if not present
       
    case('#MHDATA')
       ! type of data files
       call read_var('TypeFile', TypeFile)

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
       
       ! number of input files
       call read_var('nFileRead', nFileRead)

       ! prepare the container for file names
       allocate(NameFileStamp_I(nFileRead))
       
       ! list of files
       do iFile = 1, nFileRead
          call read_var('NameFile', NameFileStamp_I(iFile))
       end do
    end select
  end subroutine set_read_mh_data_param

  !============================================================================

  subroutine init_read_mh_data
    ! initialize by setting the time and interation index of input files
    character (len=*), parameter :: NameSub='SP:init_read_mh_data'
    !-------------------------------------------------------------------------
    if(.not.DoReadMhData)&
         RETURN
    iFileRead= 0
    ! check whether increments are properly set
    if(nFileRead <= 0)&
         call CON_stop(NameSub//&
         " invalid number of input files, change PARAM.in")
    ! read the first input file
    call read_mh_data(TimeGlobal)
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
    integer:: iBlock
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
    ! increase file counter
    iFileRead = iFileRead + 1

    ! read the data
    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameInputDir)//trim(NameFileBase)//'_',iLon,'_',iLat,&
            '_'//NameFileStamp_I(iFileRead)//NameFormat

       ! read the header first
       call read_plot_file(&
            NameFile   = NameFile,&
            TypeFileIn = TypeFile,&
            TimeOut    = TimeOut,&
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
