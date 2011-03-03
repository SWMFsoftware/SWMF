module ModSolarwind

  use ModKind
  use ModVarIndexes

  implicit none

  SAVE

  private ! except
  public :: read_solar_wind_file
  public :: normalize_solar_wind_data
  public :: get_solar_wind_point

  logical, public :: UseSolarwindFile = .false.

  character(len=500), public :: NameSolarwindFile

  integer, parameter :: nTimeVar = 7
  integer, parameter :: MaxData  = 50000

  ! Number of data points in the input file
  integer :: nData = 0

  ! Number of variables in the input file
  integer :: nVarInput

  ! Name of the input variables
  character(len=10) :: NameInputVar_I(nVar)

  ! Variable index corresponding to the i-th input variable
  integer :: iVarInput_V(max(nVar,8)) 

  ! true if the variable is read from the data file
  logical :: IsInput_V(nVar)

  ! true if number density is read instead of density
  logical :: UseNumberDensity = .true.

  ! true if temperature is read instead of pressure
  logical :: UseTemperature = .true.

  ! Arrays of input data and time
  real,         allocatable :: Solarwind_VI(:, :)
  real(Real8_), allocatable :: Time_I(:)

  ! Coordinate system of input file
  character(len=3) :: NameInputCoord = 'GSM'

  ! Normal direction to the (tilted) plane of input data
  real :: Normal_D(3) = (/ 1.0, 0.0, 0.0 /)

  ! Position of the satellite
  real :: SatelliteXyz_D(3)=0.

  ! Shall we reread the file
  logical :: DoReadAgain = .false.

contains
  !============================================================================
  subroutine read_solar_wind_file

    use ModProcMH
    use ModMain, ONLY: lVerbose, StartTime, iStartTime_I, UseStrict, &
         TypeCoordSystem
    use ModPhysics, ONLY: SW_Bx_dim, SW_By_dim, SW_Bz_dim, &
         SW_Ux_dim, SW_Uy_dim, SW_Uz_dim, SW_n_dim, SW_T_dim

    use ModIoUnit, ONLY : UNITTMP_
    use ModNumConst, ONLY: cDegToRad, cHalfPi
    use ModConst, ONLY: cDay => cSecondPerDay
    use CON_geopack, ONLY: CON_recalc, HgiGse_DD, GsmGse_DD
    use ModIO, ONLY: iUnitOut, write_prefix
    use ModTimeConvert, ONLY: time_int_to_real
    use ModUtilities, ONLY: upper_case, lower_case, split_string

    implicit none

    character(len=500):: StringInputVar

    integer :: iError, i , iVar, jVar, iYear
    logical :: UseZeroBx

    ! One line of input
    character (len=100) :: line
    character (len=10) :: String, NameVar

    real    :: TmpData_V(nVar)
    integer :: iTime_I(nTimeVar)

    ! Tilt angle of the plane of the input data
    real :: PlaneAngleXY=0.0, PlaneAngleXZ=0.0

    real :: TimeDelay
    real(Real8_) :: Time, DtData1, DtData2
    real, save:: HgiGsm_DD(3,3)
    real    :: Solarwind_V(nVar)
    logical :: DoSetMatrix=.true.

    character (len=*), parameter:: NameSub='read_solar_wind_file'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(.not.allocated(Solarwind_VI)) &
         allocate(Solarwind_VI(nVar, MaxData), Time_I(MaxData))

    ! Set defaults
    UseZeroBx = .false.
    TimeDelay = 0.0
    PlaneAngleXY = 0.0
    PlaneAngleXZ = 0.0

    nVarInput = 8
    iVarInput_V(1:8) = (/ Bx_, By_, Bz_, Ux_, Uy_, Uz_, Rho_, p_ /)
    if(UseMultiSpecies) iVarInput_V(7) = SpeciesFirst_

    ! Read and convert INPUT data on processor 0 then broadcast to others
    open(UNITTMP_, file=NameSolarwindFile, status="old", iostat = iError)

    if (iError /= 0) call stop_mpi(NameSub// &
         ": Unable to open file "//trim(NameSolarwindFile))

    nData = 0
    if(lVerbose>0 .and. iProc==0)then
       call write_prefix
       write(iUnitOut,*) NameSub,' reading ',trim(NameSolarwindFile)
    end if

    ! Read header information
    do
       read(UNITTMP_,'(a)', iostat = iError ) line
       if (iError /= 0) call stop_mpi(NameSub// &
            ': could not find #START in '//trim(NameSolarwindFile))

       if(index(line,'#REREAD')>0) read(UNITTMP_,*) DoReadAgain

       if(index(line,'#COOR')>0)then
          read(UNITTMP_,'(a)') NameInputCoord
          call upper_case(NameInputCoord)
       endif

       if(index(line,'#PLANE')>0)then
          read(UNITTMP_,*) PlaneAngleXY
          read(UNITTMP_,*) PlaneAngleXZ
          PlaneAngleXY = PlaneAngleXY * cDegToRad
          PlaneAngleXZ = PlaneAngleXZ * cDegToRad

          ! Calculate normal vector
          if( abs(abs(PlaneAngleXY)-cHalfPi) < 1.0e-3 )then
             Normal_D=(/ 0.0, sign(1.0,PlaneAngleXY), 0.0 /)
          else if ( abs(abs(PlaneAngleXZ)-cHalfPi) < 1.0e-3 )then
             Normal_D=(/ 0.0, 0.0, sign(1.0, PlaneAngleXZ) /)
          else
             Normal_D(2) = tan(PlaneAngleXY)
             Normal_D(3) = tan(PlaneAngleXZ)
             Normal_D =  Normal_D / sqrt(sum(Normal_D**2))
          end if
          if(DoTest)write(*,*)'Normal propagation direction is',Normal_D
       endif

       if(index(line,'#VAR')>0)then
          read(UNITTMP_,'(a)', iostat = iError) StringInputVar
          call split_string(StringInputVar, nVar, NameInputVar_I, nVarInput)
          UseNumberDensity = .false.
          UseTemperature   = .false.
          iVarInput_V = 0
          if(DoTest)then
             write(*,*)'StringInputVar=', StringInputVar
             write(*,*)'NameInputVar_I=',NameInputVar_I
             write(*,*)'Max Number of input variables =',nVar
             write(*,*)'user input number of variables=',nVarInput
          end if

          do iVar= 1, nVarInput
             String = NameInputVar_I(iVar)
             call lower_case(String)
             select case (String)
             case ('n')
                iVarInput_V(iVar)  = rho_
                UseNumberDensity = .true.
             case ('ux')
                iVarInput_V(iVar) = Ux_
             case ('uy')
                iVarInput_V(iVar) = Uy_
             case ('uz')
                iVarInput_V(iVar) = Uz_
             case ('bx')
                iVarInput_V(iVar) = Bx_
             case ('by')
                iVarInput_V(iVar) = By_
             case ('bz')
                iVarInput_V(iVar) = Bz_
             case ('t')
                iVarInput_V(iVar) = p_
                UseTemperature  = .true.
             case  default
                do jVar=1, nVar
                   NameVar = NameVar_V(jVar)
                   call lower_case(NameVar)
                   if(NameVar /=String) CYCLE
                   iVarInput_V(iVar)=jVar
                end do
                if(iVarInput_V(iVar)==0.)then
                   write(*,*)NameSub,': unknown solarwind variable',&
                        iVar, String
                   call stop_mpi(NameSub//' wrong solarwind variables')
                end if
             end select
          end do
       end if

       if(index(line,'#POSITION')>0)then
          read(UNITTMP_,*) SatelliteXyz_D(2)
          read(UNITTMP_,*) SatelliteXyz_D(3)
       endif

       if(index(line,'#SATELLITEXYZ')>0)then
          read(UNITTMP_,*) SatelliteXyz_D(1)
          read(UNITTMP_,*) SatelliteXyz_D(2)
          read(UNITTMP_,*) SatelliteXyz_D(3)
       endif

       if(index(line,'#ZEROBX')>0)    read(UNITTMP_,*) UseZeroBx

       if(index(line,'#TIMEDELAY')>0) read(UNITTMP_,*) TimeDelay

       if(index(line,'#START')>0) EXIT
    end do

    ! Set logicals telling if a variable is read from the input file
    do iVar = 1, nVar
       IsInput_V(iVar) = any(iVarInput_V(1:nVarInput) == iVar)
    end do

    ! Initialize array so all elements are set
    Solarwind_VI = 0.0

    ! Read the data
    do
       read(UNITTMP_, *, iostat=iError) iTime_I, TmpData_V(1:nVarInput)
       if (iError /= 0) EXIT

       if (nData >= MaxData) then
          write(*,*) "=> The solar file ",trim(NameSolarwindFile)
          write(*,*) "   contains too many lines! (Max lines =", MaxData,")."
          write(*,*) "   Try reading (overlapping) files multiple times."
          EXIT
       end if

       nData = nData + 1

       do iVar =1, nVarInput
          jVar = iVarInput_V(iVar)
          Solarwind_VI(jVar, nData) = TmpData_V(iVar)              
       end do

       if(DoTest)then
          write(*,*)'nData=', nData
          write(*,*)'Solarwind_VI=', Solarwind_VI(:, nData)
       end if

       ! Fix 2 digit years
       iYear = iTime_I(1)
       if (iYear >= 65 .and. iYear < 100) iYear = iYear + 1900
       if (iYear < 65)                    iYear = iYear + 2000
       iTime_I(1) = iYear

       ! Convert integer time into double precision real time and apply delay
       call time_int_to_real(iTime_I, Time)
       Time_I(nData) = Time + TimeDelay

       if (UseZeroBx) Solarwind_VI(Bx_, nData) = 0.0

       if (NameInputCoord /= TypeCoordSystem) then 

          if (TypeCoordSystem == 'GSM' .and. NameInputCoord == 'GSE') then 

             call CON_recalc(&
                  iTime_I(1), &
                  iTime_I(2), &
                  iTime_I(3), &
                  iTime_I(4), &
                  iTime_I(5), &
                  iTime_I(6))

             Solarwind_VI(Bx_:Bz_, nData)=&
                  matmul(GsmGse_DD, Solarwind_VI(Bx_:Bz_, nData))
             Solarwind_VI(Ux_:Uz_, nData)=&
                  matmul(GsmGse_DD, Solarwind_VI(Ux_:Uz_, nData))

          elseif(TypeCoordSystem == 'HGI' .and. NameInputCoord == 'GSM') then

             if(DoSetMatrix)then
                DoSetMatrix=.false.
                call CON_recalc( &
                     iStartTime_I(1),& ! year
                     iStartTime_I(2),& ! month
                     iStartTime_I(3),& ! day
                     iStartTime_I(4),& ! hour
                     iStartTime_I(5),& ! minute
                     iStartTime_I(6))  ! second
                HgiGsm_DD = matmul(HgiGse_DD,transpose(GsmGse_DD))
             end if
             Solarwind_VI(Bx_:Bz_, nData) = &
                  matmul(HgiGsm_DD, Solarwind_VI(Bx_:Bz_, nData))
             Solarwind_VI(Ux_:Uz_, nData)=&
                  matmul(HgiGsm_DD, Solarwind_VI(Ux_:Uz_, nData))
          else
             write(*,*) 'Transformation between input ',&
                  ' coordinate system=',NameInputCoord,&
                  ' and code system=',TypeCoordSystem,&
                  ' is not implemented'
             call stop_mpi('GM_ERROR')
          endif
       endif
    enddo

    close(UNITTMP_)

    ! Check if the start time is within 1 day of the input data
    if (StartTime < Time_I(1)-cDay .or. StartTime > Time_I(nData)+cDay) then
       write(*,*) "**********************************************************"
       write(*,*) "*                                                        *"
       write(*,*) "*      Warning! Warning! Warning! Warning! Warning!      *"
       write(*,*) "*                                                        *"
       write(*,*) "*  Time dependent solar wind file disagrees with the     *"
       write(*,*) "*  starting time of the simulation by more than 24 hours.*"
       write(*,*) "*  This could cause results which you may not enjoy.     *"
       write(*,*) "*                                                        *"
       if (UseStrict) call stop_mpi('Correct PARAM.in or the INPUT data file')

       write(*,*) "*  I am assuming that you are smarter than I am.         *"
       write(*,*) "*  Continuing with simulation.                           *"
       write(*,*) "*                                                        *"
       write(*,*) "**********************************************************"
    endif

    if(lVerbose>0 .and. iProc==0)then
       call write_prefix; write(iUnitOut,*) NameSub, &
            ' read ',nData,' points from ',trim(NameSolarwindFile)
    end if

    ! This part is only needed for solar wind normalization based on
    ! the input file. This should be eliminated.
    if (SW_T_dim <= 0.0) then

       if (nData == 1) then
          Solarwind_V = Solarwind_VI(:, 1)

       else

          ! Find the index with time before start time
          i = 1

          do while ((i < nData).and.   &
               (Time_I(i) < StartTime))
             i = i + 1
          enddo

          if ((i == nData .and. &
               Time_I(i) <= StartTime).or.(i==1)) then

             Solarwind_V = Solarwind_VI(:, i)

          else

             DtData2 = (Time_I(i) - StartTime) / &
                  (Time_I(i) - Time_I(i-1) + 1.0e-6)
             DtData1 = 1.0 - DtData2

             Solarwind_V = DtData1*Solarwind_VI(:, i) &
                  +        DtData2*Solarwind_VI(:, i-1)
          endif

       endif

       ! add up species densities
       if (UseMultiSpecies) & 
            Solarwind_V(Rho_) = sum(Solarwind_V(SpeciesFirst_:SpeciesLast_))

       ! These scalars should be removed eventually
       SW_Bx_dim  = Solarwind_V(Bx_)
       SW_By_dim  = Solarwind_V(By_)
       SW_Bz_dim  = Solarwind_V(Bz_)
       SW_Ux_dim  = Solarwind_V(Ux_)
       SW_Uy_dim  = Solarwind_V(Uy_)
       SW_Uz_dim  = Solarwind_V(Uz_)
       SW_n_dim   = Solarwind_V(rho_)
       SW_T_dim   = Solarwind_V(p_)
    endif

  end subroutine read_solar_wind_file

  !============================================================================

  subroutine normalize_solar_wind_data

    use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure
    use ModPhysics, ONLY: &
         Io2No_V, UnitTemperature_, UnitN_, UnitRho_, UnitP_, UnitU_, UnitB_, &
         LowDensityRatio, ElectronPressureRatio
    use ModMultiFluid
    use ModConst

    implicit none

    integer:: iData
    integer:: T_= p_
    real :: Solarwind_V(nVar)

    character(len=*), parameter:: NameSub = 'normalize_solar_wind_data'
    logical:: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call set_oktest(NameSub,DoTest,DoTestMe)

    do iData=1, nData
       ! Put this point into a small array
       Solarwind_V = Solarwind_VI(:, iData)

       ! normalize B and U 
       Solarwind_V(Bx_:Bz_) = Solarwind_V(Bx_:Bz_)*Io2No_V(UnitB_)
       Solarwind_V(Ux_:Uz_) = Solarwind_V(Ux_:Uz_)*Io2No_V(UnitU_)

       if(UseMultiSpecies) then

          ! We assume that species densities are always given as number density
          where(IsInput_V(SpeciesFirst_:SpeciesLast_))
             ! Normalize species number densities
             Solarwind_V(SpeciesFirst_:SpeciesLast_) = &
                  Solarwind_V(SpeciesFirst_:SpeciesLast_)*Io2No_V(UnitN_)
          elsewhere
             ! Use very small value for densities not given
             Solarwind_V(SpeciesFirst_:SpeciesLast_) = &
                  LowDensityRatio * Solarwind_V(SpeciesFirst_)
          end where

          if(UseTemperature)then
             ! calculate normalized p = n*T 
             Solarwind_V(p_) = sum(Solarwind_V(SpeciesFirst_:SpeciesLast_)) &
                  *Solarwind_V(T_)*Io2No_V(UnitTemperature_)
          else
             ! normalize pressure
             Solarwind_V(p_) = Solarwind_V(p_)*Io2No_V(UnitP_)
          end if

          ! calculate mass density of each ion species 
          Solarwind_V(SpeciesFirst_:SpeciesLast_) = &
               Solarwind_V(SpeciesFirst_:SpeciesLast_)*MassSpecies_V
          ! add up species densities
          Solarwind_V(Rho_) = sum(Solarwind_V(SpeciesFirst_:SpeciesLast_))
       else
          ! Single species
          if(UseNumberDensity) then
             Solarwind_V(Rho_) = Solarwind_V(Rho_)*Io2No_V(UnitN_)*MassIon_I(1)
          else
             Solarwind_V(Rho_) = Solarwind_V(Rho_)*Io2No_V(UnitRho_)
          end if
          if(UseTemperature) then
             Solarwind_V(p_) = Solarwind_V(T_)*Io2No_V(UnitTemperature_)&
                  *Solarwind_V(Rho_)/MassIon_I(1)
          else
             Solarwind_V(p_) = Solarwind_V(p_)*Io2No_V(UnitP_)
          end if
       end if

       ! Modify total pressure with electron pressure
       if(UseTemperature .and. .not.UseElectronPressure) &
            Solarwind_V(p_) = Solarwind_V(p_) * (1.0 + ElectronPressureRatio)

       ! Set or normalize other fluids for multi-fluid equations
       do iFluid = 2, nFluid
          call select_fluid

          ! Set fluid density
          if(.not.IsInput_V(iRho)) then
             if(iFluid == IonFirst_)then
                ! By default the solar wind contains the first ion fluid
                Solarwind_V(iRho) = Solarwind_V(Rho_)* &
                     (1.0 - LowDensityRatio*(nFluid - IonFirst_))
             else
                ! Other ion fluids are set to a very small fraction
                Solarwind_V(iRho) = Solarwind_V(Rho_)*LowDensityRatio
             end if
          elseif(UseNumberDensity)then
             Solarwind_V(iRho) = Solarwind_V(iRho) &
                  * Io2No_V(UnitN_) * MassFluid_I(iFluid)
          else
             Solarwind_V(iRho) = Solarwind_V(iRho) * Io2No_V(UnitRho_)
          end if

          ! Set fluid pressure
          if(.not.IsInput_V(iP)) then
             !By default all fluids have the same temperature as the first ion
             Solarwind_V(iP) = Solarwind_V(p_)/Solarwind_V(Rho_) &
                  *Solarwind_V(iRho)*MassIon_I(1)/MassFluid_I(iFluid)
          elseif(UseTemperature)then
             ! Calculate normalized fluid pressure from temperature
             Solarwind_V(iP) = Solarwind_V(iP)*Io2No_V(UnitTemperature_) &
                  *Solarwind_V(iRho)/MassFluid_I(iFluid)
          else
             ! Normalize fluid pressure
             Solarwind_V(iP) = Solarwind_V(iP) * Io2No_V(UnitP_)
          end if

          ! Normalize or set fluid velocities
          if(IsInput_V(iUx)) then
             Solarwind_V(iUx) = Solarwind_V(iUx) * Io2No_V(UnitU_)
          else
             Solarwind_V(iUx) = Solarwind_V(Ux_)
          end if
          if(IsInput_V(iUy)) then
             Solarwind_V(iUy) = Solarwind_V(iUy) * Io2No_V(UnitU_)
          else
             Solarwind_V(iUy) = Solarwind_V(Uy_)
          end if
          if(IsInput_V(iUz)) then
             Solarwind_V(iUz) = Solarwind_V(iUz) * Io2No_V(UnitU_)
          else
             Solarwind_V(iUz) = Solarwind_V(Uz_)
          end if

       end do ! iFluid

       ! Fix total pressure if necessary
       if(IsMhd .and. UseMultiIon)then
          Solarwind_V(p_) = sum(Solarwind_V(iPIon_I))
          if(.not.UseElectronPressure) &
               Solarwind_V(p_) = (1+ElectronPressureRatio)*Solarwind_V(p_)
       end if

       if(UseAnisoPressure .and. .not. IsInput_V(Ppar_)) &
            Solarwind_V(pPar_) = Solarwind_V(p_)

       if(UseElectronPressure .and. .not. IsInput_V(Pe_))then
          if(ElectronPressureRatio > 0.0)then
             Solarwind_V(Pe_) = Solarwind_V(p_)*ElectronPressureRatio
          else
             Solarwind_V(Pe_) = Solarwind_V(p_)
          end if
       end if

       ! Put back results in big array
       Solarwind_VI(:,iData) = Solarwind_V

    end do ! iData

    if(DoTestMe)then
       write(*,*)'Io2No_V(UnitP_)',Io2No_V(UnitP_)
       write(*,*)'Io2No_V(UnitN_)*Io2No_V(UnitTemperature_)',&
            Io2No_V(UnitN_)*Io2No_V(UnitTemperature_)     
       write(*,*)'After normalization, Solarwind_VI(:,1)=',Solarwind_VI(:,1)
       write(*,*)'After normalization, Solarwind_VI(:,2)=',Solarwind_VI(:,2)
       write(*,*)'After normalization, Solarwind_VI(:,3)=',Solarwind_VI(:,3)
    end if

  end subroutine normalize_solar_wind_data

  !===========================================================================

  subroutine get_solar_wind_point(TimeSimulation, Xyz_D, SolarWind_V)

    use ModKind
    use ModMain
    use ModPhysics
    use ModVarIndexes
    use ModGeometry, ONLY: x1, x2
    use ModUtilities, ONLY: sleep

    implicit none
    real, intent(in)  :: TimeSimulation
    real, intent(in)  :: Xyz_D(3)
    real, intent(out) :: SolarWind_V(nVar) ! Varying solar wind parameters

    integer :: iData
    real(Real8_) :: DtData1, DtData2, Time
    real         :: SatDistance_D(3), U_D(3)

    logical, parameter:: DoTestCell = .false.
    character(len=*), parameter :: NameSub = 'get_solar_wind_point'
    !--------------------------------------------------------------------------

    if(nData <= 0 .or. .not.UseSolarwindFile)then
       ! Use fixed boundary conditon if there is no input data
       SolarWind_V = FaceState_VI(:,east_)
       RETURN
    end if

    if (nData == 1) then
       SolarWind_V = Solarwind_VI(:,1)
       RETURN
    end if

    ! DoTestCell = maxval(abs(Xyz_D - (/ xTest, yTest, zTest /))) < 0.1
    if(DoTestCell) write(*,*) NameSub, 'TimeSim, Xyz_D=',TimeSimulation, Xyz_D

    ! Set the X coordinate of the sattellite if not yet set.
    if(SatelliteXyz_D(1) == 0.0)then
       ! Assume that the satellite is at the east or west boundary
       if(Solarwind_VI(Ux_,1) > 0.0)then
          SatelliteXyz_D(1) = x1  
       else
          SatelliteXyz_D(1) = x2
       endif
    end if
    ! Calculate absolute time
    Time = StartTime + TimeSimulation

    ! Read the input file again if Time exceeds last time in the file
    if(DoReadAgain .and. Time > Time_I(nData))then
       do
          call read_solar_wind_file
          call normalize_solar_wind_data
          if(Time < Time_I(nData)) EXIT
          call sleep(1.0) ! sleep for 1 second CPU time
       end do
    end if

    ! Find data point preceeding this time
    do iData = 1, nData-1
       if(Time_I(iData) >= Time) EXIT
    end do

    if(DoTestCell)then
       write(*,*)NameSub,' Time, Time_I(iData), iData=',&
            Time, Time_I(iData), iData
    end if

    ! Take into account the propagation time from the satellite plane
    if(abs(Xyz_D(1)-SatelliteXyz_D(1)) > cTiny .or. Normal_D(1) /= 1.0)then

       SatDistance_D = Xyz_D - SatelliteXyz_D

       if((iData==nData .and. Time_I(iData)<=Time) .or. iData==1) then 
          U_D = Solarwind_VI(Ux_:Uz_,iData)
       else           
          DtData2 = (Time_I(iData) - Time) / &
               (Time_I(iData) - Time_I(iData-1) + 1.0e-6)
          DtData1 = 1.0 - DtData2

          U_D =  DtData1 * Solarwind_VI(Ux_:Uz_, iData)  &
               + DtData2 * Solarwind_VI(Ux_:Uz_, iData-1)

       endif

       if(DoTestCell) write(*,*)NameSub,' u_D=',u_D

       ! Shift Time with travel time from the satellite plane
       if(abs(dot_product(U_D, Normal_D))>1.0e-5)then  
          ! Absolute time is in SI units
          Time = Time - &
               dot_product(SatDistance_D, Normal_D)/ &
               dot_product(U_D, Normal_D) * No2Si_V(UnitT_)


          ! Find data point preceeding this time
          do iData = 1, nData-1
             if(Time_I(iData) >= Time) EXIT
          end do

          if(DoTestCell) write(*,*)NameSub,&
               ' corrected Time, Time_I(iData), iData=',&
               Time, Time_I(iData), iData

       end if

    end if

    if ((iData == nData .and. Time_I(iData) <= Time) .or. iData==1) then 
       ! Use end point value
       SolarWind_V = Solarwind_VI(:, iData)

    else
       ! Interpolate in time
       DtData2 = (Time_I(iData) - Time) / &
            (Time_I(iData) - Time_I(iData-1) + 1.0e-6)
       DtData1 = 1.0 - DtData2

       SolarWind_V = DtData1 * Solarwind_VI(:, iData) &
            +        DtData2 * Solarwind_VI(:, iData-1)

       if(DoTestCell)write(*,*)NameSub,' DtData1, DtData2=',DtData1, DtData2

    endif

    if(DoTestCell)write(*,*)NameSub,' SolarWind_V =',SolarWind_V


  end subroutine get_solar_wind_point

end module ModSolarwind

