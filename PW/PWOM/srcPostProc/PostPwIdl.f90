
program PostPwIdl
  
  ! Read PW lines and create 2D altitude slices

  use ModPlotFile, ONLY: save_plot_file, read_plot_file
  use ModNumConst, ONLY: cDegToRad
  implicit none

  ! This is copied from ModKind, because PostIDL.exe may be compiled with
  ! different precision then the rest of the codes.
  integer, parameter :: Real4_=selected_real_kind(6,30)
  integer, parameter :: Real8_=selected_real_kind(12,100)
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0
  
  integer, parameter :: UnitTmp_=99, nDimIn = 1, nDimOut = 2
  integer, parameter :: z_=1, Lat_=1, Lon_=2
  character(len=20) ::  TypePlot='ascii'
!  character(len=20) ::  TypePlot='real8'
  integer :: nLine, nAltOut, nAlt, nTime, iLine, iTime, nVar, nParam, nStep
  character(len=100), allocatable :: NameFile_I(:)
  character(len=100) :: NameHeader, Type, NameOut, NamePlotVar
  real, allocatable :: PlotState_IV(:,:), PlotState_IIV(:,:,:), &
                       Coord1_I(:), Coord2_I(:), Coord_I(:),Param_I(:)
  real :: theta, phi, Time
  integer :: iVar
  
   ! Read information from STDIN
  read(*,'(a)') TypePlot
  read(*,'(i5)')nLine
  read(*,'(i5)')nAltOut
  read(*,'(i5)')nTime
  read(*,'(a)') NameOut

  ! Allocate and fill filename array
  allocate(NameFile_I(nLine))
  do iLine=1,nLine
     read(*,'(a)') NameFile_I(iLine)
  end do
  
  ! Read Header from the first file to get nVar and nAlt
  call read_plot_file(NameFile_I(1), n1Out = nAlt, nVarOut = nVar,    &
       TypeFileIn=TypePlot,StringHeaderOut = NameHeader,              & 
       NameVarOut = NamePlotVar, nStepOut= nStep,TimeOut=Time,        &
       nParamOut=nParam,iUnitIn = (100+1))  

  !Now allocate PlotState and coord arrays
  allocate (PlotState_IV(nAlt,nVar),PlotState_IIV(nLine,1,nVar))
  allocate (Coord1_I(nLine), Coord2_I(nLine), Param_I(nParam), Coord_I(nAlt))

  !open file for writing
 ! select case(TypePlot)
 ! case('formatted','ascii')
     open(UnitTmp_, file=NameOut, status='replace')
 ! case('real8','real4')
 !    open(UnitTmp_, file=NameOut, status='replace',form='unformatted')
 ! end select

  NamePlotVar(1:9) = 'Y X Z Lon'
  

  do iTime=1,nTime 
     do iLine = 1, nLine
        call read_plot_file(NameFile_I(iLine), ParamOut_I=Param_I,          &
             TypeFileIn=TypePlot, nStepOut= nStep,TimeOut=Time,             &
             CoordOut_I = Coord_I, VarOut_IV = PlotState_IV,                &
             iUnitIn = (100+iLine))
        
        ! Fill PlotState_IIV
        PlotState_IIV(iLine, 1, :) = PlotState_IV(nAltOut, :)
        
        ! Replace Lat,r, and Lon with x, y and z
        theta = (90.0-PlotState_IV(nAltOut, Lat_))*cDegToRad
        phi   = PlotState_IV(nAltOut, Lon_)       *cDegToRad
        PlotState_IIV(iLine, 1, z_) =  1.0 * cos(theta)
        !y goes into Coord1, x goes into Coord2
        Coord1_I(iLine) = -1.0 * sin(theta)*sin(phi)
        Coord2_I(iLine) =  1.0 * sin(theta)*cos(phi)
        !PlotState_IIV(iLine, 1, x_) = 1.0 * sin(theta)*cos(phi)
     enddo

     ! write out new plotfile
 !    select case(TypePlot)
 !    case('formatted', 'ascii')
        write(UnitTmp_, "(a)")             trim(NameHeader)
        write(UnitTmp_, "(i7,es13.5,3i3)") nStep, Time, -2, 1, nVar
        write(UnitTmp_, "(i8,i8)")         nLine,1
        write(UnitTmp_, "(100es13.5)")     Param_I
        write(UnitTmp_, "(a)")             trim(NamePlotVar)
        
        ! write out coordinates and variables line by line
        do iLine=1,nLine
           write(UnitTmp_, "(100es18.10)") &
                Coord1_I(iLine), Coord2_I(iLine),PlotState_IIV(iLine,1,:) 
        end do
        
 !    case('real8','real4')
 !       write(UnitTmp_)  trim(NameHeader)
 !       write(UnitTmp_)  nStep, Time, -2, 1, nVar
 !       write(UnitTmp_)  nLine,1
 !       write(UnitTmp_)  Param_I
 !       write(UnitTmp_)  trim(NamePlotVar)
 !       
 !       ! write out coordinates and variables line by line
 !       do iLine=1,nLine
 !          write(UnitTmp_) &
 !               Coord1_I(iLine), Coord2_I(iLine),PlotState_IIV(iLine,1,:) 
 !       end do
 !          
 !    end select
        
     
     

!
!   call save_plot_file(NameOut, TypePositionIn='rewind',     &
!             TypeFileIn=TypePlot,StringHeaderIn = NameHeader,               & 
!             NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=Time,           &
!             nDimIn=2,Coord1In_I = Coord1_I,Coord2In_I = Coord2_I,         &
!             VarIn_IIV = PlotState_IIV, ParamIn_I = Param_I,                &
!             IsCartesianIn=.false.)
!     else
!        call save_plot_file(NameOut, TypePositionIn='append',     &
!             TypeFileIn=TypePlot,StringHeaderIn = NameHeader,               & 
!             NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=Time,           &
!             nDimIn=2,Coord1In_I = Coord1_I,Coord2In_I = Coord2_I,         &
!             VarIn_IIV = PlotState_IIV, ParamIn_I = Param_I,                &
!             IsCartesianIn=.false.)
!     end if
  end do
  
  close(UnitTmp_)

end program PostPwIdl
!=============================================================================

subroutine CON_stop(String)

  ! This routine is needed for ModPlotFile

  implicit none

  character(len=*), intent(in):: String
  write(*,*) 'ERROR in PostIDL: '//String
  stop

end subroutine CON_stop
