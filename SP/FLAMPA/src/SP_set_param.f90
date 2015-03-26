!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModSetParam
  use ModReadParam
  implicit none
contains
  subroutine SP_set_parameters(TypeAction)
    use SP_ModMain
    character (len=*), intent(in)     :: TypeAction ! What to do
    character (len=*), parameter :: NameSub='SP_set_parameters'
    ! The name of the command
    character (len=100) :: NameCommand
    integer::iVerbose=0
    !-------------
    write(iStdOut,*)NameSub//': CHECK iSession =',i_session_read()
    if(DoInit)then
       DoInit=.false.
    end if
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          DoRestart=.true.
       case('#LINE')
          call read_var('XyzLine_D(x_)',XyzLine_D(1))
          call read_var('XyzLine_D(y_)',XyzLine_D(2))
          call read_var('XyzLine_D(z_)',XyzLine_D(3))
          call read_var('RBoundSC'     ,RBoundSC    )!^CMP IF SC 
          call read_var('RBoundIH'     ,RBoundIH    )!^CMP IF IH
       case('#RTRANSIENT')
          !           call read_var('rTransient',rTransient)
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#SAVEMHDATA')
          call read_var('SaveMhData',SaveMhData)
       case('#DOREADMHDATA')
          call read_var('DoReadMhData',DoReadMhData)
       case('#NSTEP')
          call read_var('nStep',iDataSet)
       case('#TSIMULATION')
          call read_var('tSimulation',SP_Time)
       case('#PLOT')
          !           call read_var('DnPlot',kfriss)
       case('#VERBOSE')
          call read_var('iVerbose',iVerbose)
          if(iVerbose>0)DoWriteAll=.true.
          if(DoWriteAll.and.iProc==0)&
               write(*,*)prefix,' Verbose everything'
       case('#NSMOOTH')
          call read_var('nSmooth',nSmooth)
       case default
          call CON_stop(NameSub//&
               ': Unknown command '&
               //NameCommand)
       end select
    end do
  end subroutine SP_set_parameters
end module SP_ModSetParam
