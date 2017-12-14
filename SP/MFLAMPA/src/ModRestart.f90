module SP_ModRestart

  ! This module contains methods for writing output files

  use SP_ModSize, ONLY: &
       nLon, nLat, nParticleMax, nMomentumBin

  use SP_ModGrid, ONLY: &
       get_node_indexes, &
       iProc, Z_,&
       nVarRead, nBlock, State_VIB, iGridLocal_IB, iNode_B, &
       RMin, RBufferMin, RBufferMax, RMax, &
       Distribution_IIB,  ParamLocal_IB, &
       End_, &
       nBlockParam, nBlockIndexes

   use SP_ModAdvance, ONLY: TimeGlobal, iIterGlobal

   use ModPlotFile, ONLY: save_plot_file, read_plot_file

  use ModUtilities, ONLY: open_file, close_file
  use ModIoUnit,    ONLY: UnitTmp_

  implicit none

  SAVE

  private ! except
 
  public:: set_restart_param, write_restart, read_restart
  public:: NameRestartInDir, NameRestartOutDir

  integer, parameter:: nBufferMax = &
       nBlockParam+nBlockIndexes+nParticleMax*(Z_+nMomentumBin)
  real, allocatable:: Buffer_I(:)


  !----------------------------------------------------------------------------
  ! the restart directory
  character (len=100) :: NameRestartOutDir="SP/restartOUT/"
  character (len=100) :: NameRestartInDir ="SP/restartIN/"
  ! name of the header file
  character (len=100) :: NameHeaderFile   ="restart.H"
  !----------------------------------------------------------------------------
  !/
contains
  
  subroutine set_restart_param
    use ModReadParam, ONLY: read_var
    ! set parameters of restart files
    character(len=*), parameter :: NameSub='SP:set_restart_param'
    !--------------------------------------------------------------------------
  end subroutine set_restart_param

  !============================================================================

  subroutine write_restart
    ! write the restart data
 
    ! name of the output file
    character(len=100):: NameFile
    ! loop variables
    integer:: iBlock, iParticle, i
    ! indexes of corresponding node, latitude and longitude
    integer:: iNode, iLat, iLon
    !
    integer:: nBuffer
    ! index of first/last particle on the field line
    integer:: iLast

    character(len=*), parameter:: NameSub = 'SP:write_restart'
    !--------------------------------------------------------------------------
    if(.not.allocated(Buffer_I))&
         allocate(Buffer_I(nBufferMax))

    call write_restart_header

    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartOutDir)//'data_',iLon,'_',iLat,&
            '.rst'

       ! fill the output buffer
       nBuffer = 0

       ! general parameters
       do i = 1, nBlockIndexes
          nBuffer = nBuffer + 1
          Buffer_I(nBuffer) = &
               real(iGridLocal_IB(i, iBlock))
       end do
       do i = 1, nBlockParam
          nBuffer = nBuffer + 1
          Buffer_I(nBuffer) = ParamLocal_IB(i, iBlock)
       end do

       ! get max particle indexes on this field line
       iLast  = iGridLocal_IB(End_,   iBlock)

       do iParticle = 1, iLast
          ! background plasma paramters
          do i = 0, Z_
             nBuffer = nBuffer + 1
             Buffer_I(nBuffer) = State_VIB(i, iParticle, iBlock)
          end do

          ! distribution
          do i = 1, nMomentumBin
             nBuffer = nBuffer + 1
             Buffer_I(nBuffer) = Distribution_IIB(i, iParticle, iBlock)
          end do
       end do

       ! print data to file
       call save_plot_file(&
            NameFile   = NameFile, &
            TypeFileIn = 'real8', &
            TimeIn     = TimeGlobal, &
            nStepIn    = iIterGlobal, &
            CoordMinIn_D = (/real(1)/), &
            CoordMaxIn_D = (/real(iLast)/), &
            VarIn_I    = Buffer_I(1:nBuffer)&
            )
    end do

  end subroutine write_restart

  !============================================================================

  subroutine read_restart
    ! read the restart data

    ! name of the input file
    character(len=100):: NameFile
    ! loop variables
    integer:: iBlock, iParticle, i
    ! indexes of corresponding node, latitude and longitude
    integer:: iNode, iLat, iLon
    !
    integer:: nBuffer
    ! index of first/last particle on the field line
    integer:: iLast

    character(len=*), parameter:: NameSub = 'SP:read_restart'
    !--------------------------------------------------------------------------
    if(.not.allocated(Buffer_I))&
         allocate(Buffer_I(nBufferMax))

    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartInDir)//'data_',iLon,'_',iLat,&
            '.rst'

       ! read data from file
       call read_plot_file(&
            NameFile   = NameFile, &
            TypeFileIn = 'real8', &
            VarOut_I    = Buffer_I(:)&
            )

       ! process buffer
       nBuffer = 0

       ! general parameters
       do i = 1, nBlockIndexes
          nBuffer = nBuffer + 1
          iGridLocal_IB(i, iBlock) = nint(Buffer_I(nBuffer))
       end do

       do i = 1, nBlockParam
          nBuffer = nBuffer + 1
          ParamLocal_IB(i, iBlock) = Buffer_I(nBuffer)
       end do

       ! get min and max particle indexes on this field line
       iLast  = iGridLocal_IB(End_,   iBlock)

       do iParticle = 1, iLast
          ! background plasma paramters
          do i = 0, Z_
             nBuffer = nBuffer + 1
             State_VIB(i, iParticle, iBlock) = Buffer_I(nBuffer)
          end do

          ! distribution
          do i = 1, nMomentumBin
             nBuffer = nBuffer + 1
             Distribution_IIB(i, iParticle, iBlock) = Buffer_I(nBuffer)
          end do
       end do
    end do

  end subroutine read_restart


  subroutine write_restart_header

    ! full name of the header file
    character(len=100):: NameFile

    character(len=*), parameter:: NameSub='write_restart_header'
    !--------------------------------------------------------------------------
    if (iProc/=0) RETURN
    NameFile = trim(NameRestartOutDir)//trim(NameHeaderFile)

    call open_file(file=NameFile, NameCaller=NameSub)
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#RESTART'
    write(UnitTmp_,'(a)')'T'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#CHECKGRIDSIZE'
    write(UnitTmp_,'(i8,a32)') nParticleMax,'nParticleMax'
    write(UnitTmp_,'(i8,a32)') nLon,     'nLon'
    write(UnitTmp_,'(i8,a32)') nLat,     'nLat'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#NSTEP'
    write(UnitTmp_,'(i8,a32)')iIterGlobal,'nStep'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#TIMESIMULATION'
    write(UnitTmp_,'(es22.15,a18)')TimeGlobal,'tSimulation'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#GRID'
    write(UnitTmp_,'(es22.15,a18)')RMin,      'RMin'
    write(UnitTmp_,'(es22.15,a18)')RBufferMin,'RBufferMin'
    write(UnitTmp_,'(es22.15,a18)')RBufferMax,'RBufferMax'
    write(UnitTmp_,'(es22.15,a18)')RMax,      'RMin'
    write(UnitTmp_,*)
    write(UnitTMP_,'(a)')'#END'
    write(UnitTmp_,*)
    call close_file

  end subroutine write_restart_header

end module SP_ModRestart
