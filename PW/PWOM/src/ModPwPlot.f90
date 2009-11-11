Module ModPwPlots
  
  implicit none
  
  private ! except 
  
  public :: PW_print_plot

  character(len=5),  public    :: TypePlot   = 'ascii'
  character(len=22), parameter :: NameHeader = 'Polarwind output_var11'  
contains
  !=============================================================================
  subroutine PW_print_plot
    use ModCommonVariables
    use ModIoUnit, ONLY: UnitTmp_
    use ModPWOM,   ONLY: iLine,NameGraphics
    use ModPlotFile,ONLY: save_plot_file
    
    real MO,MH,MHe,Me
    real, allocatable :: PlotState_IV(:,:)
    real, allocatable :: Coord_I(:)
    integer :: iAlt, iIon
    !---------------------------------------------------------------------------
    ! Allocate PlotState and Coord arrays
    if (.not.allocated(PlotState_IV)) allocate(PlotState_IV(0:nDim+1,nPlotVar))
    if (.not.allocated(Coord_I)) allocate(Coord_I(0:nDim+1))
    
    !\
    ! Fill PlotState_IV array 
    !/

    ! Set Lat Lon
    PlotState_IV (0:nDim+1, 1) = GmLat
    PlotState_IV (0:nDim+1, 2) = GmLong

    ! Set Velocity 
    do iIon=1,nIon
       PlotState_IV (0:nDim+1, iIon+2) = State_GV(0:nDim+1,iU_I(iIon))/1.E5
    end do
        
    ! Set ion densities
    do iIon=1,nIon
       do iAlt= 0,nDim+1
          PlotState_IV (iAlt, iIon+2+nIon) = &
               alog10_check(State_GV(iAlt,iRho_I(iIon))/Mass_I(iIon))
       end do
    end do
    
    ! Set Temperatures
    do iIon=1,nIon
       PlotState_IV (0:nDim+1, iIon+2+2*nIon)  = State_GV(0:nDim+1,iT_I(iIon))
    end do

    ! Set Mach number output
    do iIon=1,nIon
       PlotState_IV (0:nDim+1, iIon+2+3*nIon) = & 
         State_GV(0:nDim+1,iU_I(iIon))/sqrt(gamma*State_GV(0:nDim+1,iP_I(iIon))&
         / State_GV(0:nDim+1,iRho_I(iIon)))
    end do
    
    ! Set Efield output
    PlotState_IV (1:nDim, 3+4*nIon)   = Efield(1:nDim)
    PlotState_IV (nDim+1, 3+4*nIon)   = Efield(nDim)
    PlotState_IV (0, 3+4*nIon)        = Efield(1)
    
    ! Set Pe for output
    PlotState_IV (0:nDim+1, 4+4*nIon)   = State_GV(0:nDim+1,pE_)
    
    ! Set altitude for output
    Coord_I (1:nDim) = AltD(1:nDim)
    Coord_I (0)      = AltMin
    Coord_I (nDim+1) = AltMax

    !write plot
    call save_plot_file(NameGraphics(iLine), TypePositionIn='append',     &
         TypeFileIn=TypePlot,StringHeaderIn = NameHeader,                 & 
         NameVarIn = NamePlotVar, nStepIn= nint(time/dt),TimeIn=time,     &
         nDimIn=1,CoordIn_I = Coord_I, VarIn_IV = PlotState_IV,           &
         ParamIn_I = (/gamma/))
    
    ! Deallocate to save memory
    deallocate(PlotState_IV)
    deallocate(Coord_I)
    RETURN
  END SUBROUTINE PW_print_plot
  
  !========================================================================
  real function alog10_check(x)
    
    implicit none
    real, intent(in) :: x
    
    !---------------------------------------------------------------------- 
    
    if(x < 0.0) then
       write(*,*)'negative argument for alog10_check:',x
       call CON_stop('PWOM ERROR: negative argument for alog10')
    endif
    alog10_check=alog10(x)
    
  end function alog10_check
  
  
end Module ModPwPlots
