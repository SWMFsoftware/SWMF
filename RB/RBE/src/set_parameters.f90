subroutine RB_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use rbe_time, ONLY: dt
  use rbe_constant
  use rbe_cread1
  use rbe_cread2
  use rbe_io_unit
  
  implicit none

  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=*), parameter  :: NameSub = 'RB_set_parameters'
  Logical :: IsStandAlone=.true.

  !\
  ! Description:
  ! This subroutine gets the inputs for RBE
  !/

  !---------------------------------------------------------------------------
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case('#STOP')
        if(IsStandAlone)then
           call read_var('Tmax',Tmax)
        else
           write(*,*)'RB WARNING: #STOP command is ignored in the framework'
        end if
     case('#SAVEPLOT')
        call read_var('DtSavePlot',tint)   ! output results every tint seconds
        call read_var('OutName',OutName)
     case('#RESTART')
        call read_var('iType',iType)       ! 1=initial run,  2=continuous run
     case('#TIMESTEP')
        call read_var('Dt',Dt)             ! time step in s. 
                                           ! Summer 2006: read dt from *.dat
     case('#STARTTIME')
        call read_var('tStart',tStart)
     case('#SPECIES')
        call read_var('js',js)             ! species: 1=RB e-, 2=RB H+
     case('#STARTUPTIME')
        call read_var('tStartup',trans)    ! startup time in sec when itype=1
     case('#BMODEL')
        call read_var('iModel',iMod)       ! 1=t96_01, 2=t0UnitTmp__s, 3=MHD
        call read_var('FixedB',ires)       ! 0=fixed B config or 
                                           ! 1=changing B config
     case('#CONVECT')
        call read_var('iConvect',iConvect) ! 1=Weimer, 2=MHD
     case('#STORM')
        call read_var('NameStorm',storm)
     case('#PLASMASPHERE')
        call read_var('PlasmaSphere',iplsp)! 0=no plasmasphere, 1=plasmasphere

     endselect


  enddo
  

end subroutine RB_set_parameters
