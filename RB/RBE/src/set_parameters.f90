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
  logical :: IsRestart, UseFixedB, UsePlasmaSphere, UseSmooth
  logical :: DoDiffusePA, DoDiffuseE
  character (len=100) :: NameSpecies, NameModel
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
        
     case('#PLOTELECTRODYNAMICS')
        call read_var('DoSaveIe',DoSaveIe)

     case('#RESTART')
        call read_var('IsRestart',IsRestart) !T:Continuous run
                                             !F:Initial run
        if (IsRestart) then
           iType = 2
        else
           iType = 1
        endif
     
     case('#TIMESTEP')
        call read_var('Dt',Dt)             ! time step in s. 
                                           ! Summer 2006: read dt from *.dat

     case('#SPLITING')
        call read_var('UseSplitting', UseSplitting)
        if (.not. UseSplitting) call read_var('UseCentralDiff', UseCentralDiff)
        
     case('#LIMITER')
        call read_var('UseMcLimiter', UseMcLimiter)
        if(UseMcLimiter) call read_var('BetaLimiter', BetaLimiter)

     case('#TIMESIMULATION')
        call read_var('TimeSimulation',tStart)
     
     case('#SPECIES')
        call read_var('NameSpecies',NameSpecies)
        if (NameSpecies == 'e') then       ! species: 1=RB e-, 2=RB H+
           js=1
        elseif(NameSpecies == 'H')then
           js =2
        else
           call con_stop('Error: Species not found')
        endif
     
     case('#STARTUPTIME')
        call read_var('tStartup',trans)    ! startup time in sec when itype=1
     
     case('#BMODEL')
        call read_var('NameModel',NameModel)!t96,t04,MHD
        call read_var('UseFixedB',UseFixedB)!T=fixed B config or 
                                            !F=changing B config
        if (NameModel == 't96') then
           iMod=1
        elseif(NameModel == 't04') then
           iMod=2
        elseif(NameModel == 'MHD')then
           iMod=3
        else
           call con_stop('Error: Model not found') 
        endif
        if (UseFixedB) then
           ires = 0
        else
           ires = 1
        endif
       
     case('#IEMODEL')
        call read_var('iConvect',iConvect) ! 1=Weimer, 2=MHD
     
     case('#INPUTDATA')
        call read_var('NameStorm',storm)
     
     case('#PLASMASPHERE')
        call read_var('UsePlasmaSphere',UsePlasmaSphere)
        if(UsePlasmaSphere)then
           iplsp=1
        else
           iplsp=0
        endif
     case('#SMOOTH')
        call read_var('UseSmooth',UseSmooth)
        if(UseSmooth)then
           ismo=1
        else
           ismo=0
        end if
     case('#WAVES')
        call read_var('DoDiffusePA',DoDiffusePa)
        if(DoDiffusePA)then
           idfa=1
        else
           idfa=0
        end if
        call read_var('DoDiffuseE',DoDiffuseE)
        if(DoDiffuseE)then
           idfe=1
        else
           idfe=0
        end if

     end select


  enddo
  

end subroutine RB_set_parameters
