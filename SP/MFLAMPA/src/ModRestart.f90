!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_ModRestart
  ! This module contains methods for writing output files
  use SP_ModSize,   ONLY: nLon, nLat, nParticleMax
  use SP_ModGrid,   ONLY: get_node_indexes, LagrID_, Z_,&
       nBlock, State_VIB, iShock_IB, iNode_B, &
       !RMin, RBufferMin, RBufferMax, RMax, &
       FootPoint_VB, nParticle_B, nShockParam
  use SP_ModDistribution, ONLY: Distribution_IIB
  use SP_ModTime,   ONLY: SPTime, iIter
  use ModPlotFile,  ONLY: save_plot_file, read_plot_file
  use ModUtilities, ONLY: open_file, close_file
  use ModIoUnit,    ONLY: UnitTmp_
  implicit none
  SAVE
  private ! except
  !Public members
  public:: save_restart, read_restart
  public:: NameRestartInDir, NameRestartOutDir

  !----------------------------------------------------------------
  ! the restart directory
  character (len=100) :: NameRestartOutDir="SP/restartOUT/"
  character (len=100) :: NameRestartInDir ="SP/restartIN/"
  ! name of the header file
  character (len=100) :: NameHeaderFile   ="restart.H"
  !----------------------------------------------------------------
  !/
contains
  !================================================================
  subroutine save_restart
    use ModIoUnit,     ONLY: UnitTmp_
    use ModUtilities,  ONLY: open_file, close_file
    ! write the restart data
 
    ! name of the output file
    character(len=100):: NameFile
    ! loop variable
    integer:: iBlock
    ! indexes of corresponding node, latitude and longitude
    integer:: iNode, iLat, iLon
    character(len=*), parameter:: NameSub = 'SP:save_restart'
    !--------------------------------------------------------------
 
    call write_restart_header

    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartOutDir)//'data_',iLon,'_',iLat,&
            '.rst'
       call open_file(file=NameFile, form='UNFORMATTED',&
            NameCaller=NameSub)
       write(UnitTmp_)real(nParticle_B(iBlock)),&
            real(iShock_IB(:, iBlock))
       write(UnitTmp_)&
            FootPoint_VB(:, iBlock),&
            State_VIB(LagrID_:Z_,1:nParticle_B(iBlock), iBlock),&
            Distribution_IIB(:,1:nParticle_B(iBlock), iBlock)
       call close_file
    end do
  end subroutine save_restart
  !================================================================
  subroutine read_restart
    ! read the restart data

    ! name of the input file
    character(len=100):: NameFile
    ! loop variables
    integer:: iBlock
    ! indexes of corresponding node, latitude and longitude
    integer:: iNode, iLat, iLon
    real   :: Aux, Aux_I(nShockParam) !For reading integers
    integer:: iError
    character(len=*), parameter:: NameSub = 'SP:read_restart'
    !--------------------------------------------------------------
   
    do iBlock = 1, nBlock
       iNode = iNode_B(iBlock)
       call get_node_indexes(iNode, iLon, iLat)

       ! set the file name
       write(NameFile,'(a,i3.3,a,i3.3,a)') &
            trim(NameRestartInDir)//'data_',iLon,'_',iLat,&
            '.rst'
       call open_file(file=NameFile,  status='old',&
            form='UNFORMATTED', NameCaller=NameSub)
       read(UnitTmp_,iostat = iError)Aux, Aux_I
       if(iError>0)then
          write(*,*)'Error in reading nPoint in Block=', iBlock
          call close_file
          call CON_stop('Run stops')
       end if
       ! process buffer
       nParticle_B(iBlock) = nint(Aux)
       ! general parameters
       iShock_IB(:, iBlock) = nint(Aux_I)
       read(UnitTmp_, iostat = iError) &
            FootPoint_VB(:, iBlock),&
            State_VIB(LagrID_:Z_,1:nParticle_B(iBlock), iBlock),&
            Distribution_IIB(:,1:nParticle_B(iBlock), iBlock)
       call close_file
    end do
  end subroutine read_restart
  !================================================================
  subroutine write_restart_header
    use SP_ModPlot, ONLY: nTag
    use SP_ModProc, ONLY: iProc
    use SP_ModDistribution, ONLY: nP, EnergyInjIo, EnergyMaxIo
    ! full name of the header file
    character(len=100):: NameFile

    character(len=*), parameter:: NameSub='write_restart_header'
    !--------------------------------------------------------------
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
    write(UnitTmp_,'(i8,a32)')iIter,'nStep'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#TIMESIMULATION'
    write(UnitTmp_,'(es22.15,a18)')SPTime,'tSimulation'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#NTAG'
    write(UnitTmp_,'(i8,a32)')nTag,'nTag'
    write(UnitTmp_,*)
    write(UnitTmp_,'(a)')'#MOMENTUMGRID'
    write(UnitTmp_,'(es22.15,a18)')EnergyInjIo,'EnergyMin'
    write(UnitTmp_,'(es22.15,a18)')EnergyMaxIo,'EnergyMax'
    write(UnitTmp_,'(i8,a32)')nP,'nP'
    write(UnitTmp_,*)
    write(UnitTMP_,'(a)')'#END'
    write(UnitTmp_,*)
    call close_file
  end subroutine write_restart_header
end module SP_ModRestart
