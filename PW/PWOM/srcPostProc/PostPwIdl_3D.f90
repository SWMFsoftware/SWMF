
program PostPwIdl
  
  ! Read PW lines and create 2D altitude slices

  use ModPlotFile, ONLY: save_plot_file, read_plot_file
  use ModNumConst, ONLY: cDegToRad,cPi
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
  real, parameter   :: rEarth=6375.0
  logical, parameter:: UseDipole=.true.
  integer :: nLine, nAltOut, nAlt, nTime, iLine, iTime, nVar, nParam, nStep
  character(len=100), allocatable :: NameFile_I(:)
  character(len=100) :: NameHeader, Type, NameOut, NamePlotVar
  real, allocatable :: PlotState_IV(:,:), PlotState_IIIV(:,:,:,:), &
                       Coord1_III(:,:,:), Coord2_III(:,:,:), Coord3_III(:,:,:), Coord_I(:),Param_I(:)
  real :: theta0,theta, phi, Time,L,r
  integer :: iVar,iAlt
  
   ! Read information from STDIN
  read(*,'(a)') TypePlot
  read(*,'(i5)')nLine
  read(*,'(i5)')nTime

  ! Allocate and fill filename array
  allocate(NameFile_I(nLine))
  do iLine=1,nLine
     read(*,'(a)') NameFile_I(iLine)
  end do
  
  ! Read Header from the first file to get nVar and nAlt
  call read_plot_file(NameFile_I(1), n1Out = nAlt, nVarOut = nVar,    &
       TypeFileIn=TypePlot,StringHeaderOut = NameHeader,              & 
       NameVarOut = NamePlotVar, nStepOut= nStep,TimeOut=Time,        &
       nParamOut=nParam,iUnitIn = UnitTmp_)  

  !Now allocate PlotState and coord arrays
  allocate (PlotState_IV(nAlt,nVar),PlotState_IIIV(nTime,nLine,nAlt,nVar))
  allocate (Coord1_III(nTime,nAlt,nLine), Coord2_III(nTime,nAlt,nLine), &
       Coord3_III(nTime,nAlt,nLine), Param_I(nParam), Coord_I(nAlt))

  !open file for writing
 ! select case(TypePlot)
 ! case('formatted','ascii')
 ! case('real8','real4')
 !    open(UnitTmp_, file=NameOut, status='replace',form='unformatted')
 ! end select

  NamePlotVar(1:9) = 'Y X Z Lon'
  

  do iLine = 1, nLine
     do iTime=1,nTime 
        write(*,*) NameFile_I(iLine),itime
        call read_plot_file(NameFile_I(iLine), ParamOut_I=Param_I,          &
             TypeFileIn=TypePlot, nStepOut= nStep,TimeOut=Time,             &
             CoordOut_I = Coord_I, VarOut_IV = PlotState_IV,                &
             iUnitIn = UnitTmp_)
        
        ! Fill PlotState_IIV
        PlotState_IIIV(iTime,iLine, :, :) = PlotState_IV(:, :)
        
        ! Replace Lat,r, and Lon with x, y and z
        theta0 = (90.0-PlotState_IV(1, Lat_))*cDegToRad
        phi   = PlotState_IV(1, Lon_)       *cDegToRad
        L=(1.0/sin(theta0))**2.0
        !y goes into Coord1, x goes into Coord2
        if (UseDipole) then
           do iAlt=1,nAlt
              r=(rEarth+Coord_I(iAlt))/rEarth
              theta=cPi/2.0-acos(min(sqrt((r*(sin(theta0))**2.0)),1.0))

              Coord1_III(iTime,iAlt,iLine) = -r * sin(theta)*sin(phi)
              Coord2_III(iTime,iAlt,iLine) =  r * sin(theta)*cos(phi)
              Coord3_III(iTime,iAlt,iLine) =  r * cos(theta)
           enddo
        else
           Coord1_III(iTime,iAlt,iLine) = -1.0 * sin(theta0)*sin(phi)
           Coord2_III(iTime,iAlt,iLine) =  1.0 * sin(theta0)*cos(phi)
           Coord3_III(iTime,iAlt,iLine) = Coord_I(iAlt)/rEarth
        endif
        !PlotState_IIV(iLine, 1, x_) = 1.0 * sin(theta)*cos(phi)
     enddo
     close(UnitTmp_)
  enddo
  
  ! write out new plotfile
  
        ! write out coordinates and variables line by line
  
  do iTime=1,nTime
     write(NameOut,"(a,i8.8,a)") &
          'plots/3DPw',iTime,'.dat'
     open(UnitTmp_, file=NameOut, status='replace')
     write(UnitTmp_,'(a)') &
          'VARIABLES = "X", "Y", "Z", "Lat", "Lon", "uO", "uH", "uHe", "ue", "lgnO", "lgnH", '//&
          '"lgnHe", "lgne", "TO", "TH", "THe", "Te", "MO", "MH", "MHe", "Me", "Ef", "Pe"' 
     write(UnitTmp_,'(a,i3,a,i3,a,i9,a)') 'Zone I=', 1, ', J=', 1,', K=',nAlt*nLine,', DATAPACKING=POINT'
     do iLine=1,nLine
        do iAlt=1,nAlt
           write(UnitTmp_, "(100es18.10)") &
              Coord2_III(iTime,iAlt,iLine), Coord1_III(iTime,iAlt,iLine),&
              Coord3_III(iTime,iAlt,iLine),PlotState_IIIV(iTime,iLine,iAlt,:) 
        end do
     enddo
     close(UnitTmp_)
  enddo
  
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
  
  

end program PostPwIdl
!=============================================================================

subroutine CON_stop(String)

  ! This routine is needed for ModPlotFile

  implicit none

  character(len=*), intent(in):: String
  write(*,*) 'ERROR in PostIDL: '//String
  stop

end subroutine CON_stop
