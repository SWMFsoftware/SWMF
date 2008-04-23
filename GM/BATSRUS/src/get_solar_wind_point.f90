!Modified by Yingjuan for user defined IMF input variables

module ModUpstreamData
  use ModKind
  use ModVarIndexes
  use ModMultiFluid, ONLY: select_fluid, iFluid, &
       iRho, iRhoUx, iRhoUy, iRhoUz, iP, MassIon_I
  use ModPhysics, ONLY: LowDensityRatio, inv_g
  use ModNumConst, ONLY: cTiny8  
  implicit none

  integer, parameter :: MaxVarTime=7
  integer, parameter :: MaxVarImf =nVar
  integer, parameter :: Max_Upstream_Npts = 50000
  integer :: Upstream_Npts=0

  integer :: nVarImf
  integer :: iVarImf_V(MaxVarImf) !indexed by imf variable indexes 
                                  !and returns a normal variable index
  logical :: UseNumberDensity=.true., UseTemperature=.true.
  logical :: UseInputImfVars=.false.  !default old list of imf vars
  character (LEN=255):: NameImfVars
  character (LEN=10) :: NameImfVar_I(MaxVarImf)
  
  real         :: Upstream_Data(Max_Upstream_Npts, MaxVarImf)
  real(Real8_) :: Upstream_Time(Max_Upstream_Npts)

  
  ! Solar Wind Input Variables
  character(len=3) :: TypeInputCoordSystem = 'GSM'
  real :: Propagation_Plane_XY=0.0, Propagation_Plane_XZ=0.0
  real :: Normal_D(3)=(/1.0, 0.0, 0.0/), MagNormal

end Module ModUpstreamData

!============================================================================

subroutine read_upstream_input_file(upstreamfilename)
  use ModProcMH
  use ModMain
  use ModIoUnit, ONLY : UNITTMP_
  use ModPhysics
  use ModNumConst
  use ModUpstreamData
  use CON_geopack
  use ModIO, ONLY: iUnitOut, write_prefix
  use ModTimeConvert, ONLY: time_int_to_real
  use ModUtilities, ONLY: upper_case, lower_case
  use ModMpi
  implicit none

  character (len=100), intent(in) :: upstreamfilename

  integer :: iError, i , iVar, jVar
  logical :: UseZeroBx

  ! One line of input
  character (len=100) :: line
  character (len=10) :: string, NameVar


  real :: TmpData_V(MaxVarImf)
  integer :: Tmp_Time(MaxVarTime)
  integer :: Upstream_integer_time(Max_Upstream_Npts,MaxVarTime)


  real :: TimeDelay
  real(Real8_) :: DtData1, DtData2
  real, save:: HgiGsm_DD(3,3)
  logical   :: DoSetMatrix=.true.

  character (len=*), parameter:: NameSub='read_upstream_input_file'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  UseZeroBx = .false.
  TimeDelay = cZero

  ! Read and convert IMF data on processor 0 then broadcast to others
  if (iProc == 0) then

     open(UNITTMP_, file=upstreamfilename, status="old", iostat = iError)

     if (iError /= 0) call stop_mpi(NameSub// &
          ": Unable to open file "//trim(upstreamfilename))

     Upstream_Npts = 0
     Propagation_Plane_XY = 0.0
     Propagation_Plane_XZ = 0.0

     nVarImf = 8
     iVarImf_V(1)=Bx_
     iVarImf_V(2)=By_
     iVarImf_V(3)=Bz_
     iVarImf_V(4)=Ux_
     iVarImf_V(5)=Uy_
     iVarImf_V(6)=Uz_
     iVarImf_V(7)=Rho_
     iVarImf_V(8)=p_
        
     if(lVerbose>0)then
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' reading ',trim(upstreamfilename)
     end if

     ! Read header information
     do
        read(UNITTMP_,'(a)', iostat = iError ) line
        if (iError /= 0) call stop_mpi(NameSub// &
             ': could not find #START in '//trim(upstreamfilename))

        if(index(line,'#COOR')>0)then
           read(UNITTMP_,'(a)') TypeInputCoordSystem
           call upper_case(TypeInputCoordSystem)
        endif

        if(index(line,'#PLANE')>0)then
           read(UNITTMP_,*) Propagation_Plane_XY
           read(UNITTMP_,*) Propagation_Plane_XZ
           Propagation_Plane_XY = Propagation_Plane_XY * cDegToRad
           Propagation_Plane_XZ = Propagation_Plane_XZ * cDegToRad
           Normal_D(2) = tan(Propagation_Plane_XY)
           Normal_D(3) = tan(Propagation_Plane_XZ)
           MagNormal = sqrt(dot_product(Normal_D, Normal_D))
           if(MagNormal > 1.0e5 )call stop_mpi('too large angle with X plane')
           Normal_D =  Normal_D /MagNormal
        endif

        if(index(line,'#VAR')>0)then
           read(UNITTMP_,'(a)', iostat = iError )NameImfVars
           call split_str(NameImfVars, MaxVarImf, NameImfVar_I,nVarImf)
           UseNumberDensity=.false.
           UseTemperature =.false.
           iVarImf_V(:)=0
           if(DoTest)then
              write(*,*)'NameImfVars=', NameImfVars
              write(*,*)'NameImfVar_I=',NameImfVar_I
              write(*,*)'Max Number of input variables =',MaxVarImf
              write(*,*)'user input number of variables=',nVarImf
           end if
           Do iVar= 1, nVarImf
              String=NameImfVar_I(iVar)
              call lower_case(String)
              select case ( String)
              case ('n')
                 iVarImf_V(iVar)=rho_
                 UseNumberDensity=.true.
              case ('ux')
                 iVarImf_V(iVar)=Ux_
              case ('uy')
                 iVarImf_V(iVar)=Uy_
              case ('uz')
                 iVarImf_V(iVar)=Uz_
              case ('bx')
                 iVarImf_V(iVar)=Bx_
              case ('by')
                 iVarImf_V(iVar)=By_
              case ('bz')
                 iVarImf_V(iVar)=Bz_
              case ('t')
                 iVarImf_V(iVar)=p_
                 UseTemperature =.true.
              case  default
                 do jVar=1, nVar
                    NameVar=NameVar_V(jVar)
                    call lower_case(NameVar)
                    if(NameVar /=String) CYCLE
                    iVarImf_V(iVar)=jVar
                 end do
                 if(iVarImf_V(iVar)==0.)then
                    write(*,*)'unknown input variable', iVar, String
                    !call stop_mpi('wrong input variables',iVar, String)
                 end if
              end select
           end do
        end if

        if(index(line,'#POSITION')>0)then
           read(UNITTMP_,*) Satellite_Y_Pos
           read(UNITTMP_,*) Satellite_Z_Pos
        endif

        if(index(line,'#ZEROBX')>0)then
           read(UNITTMP_,*) UseZeroBx
        endif

        if(index(line,'#TIMEDELAY')>0)then
           read(UNITTMP_,*) TimeDelay
        endif

        if(index(line,'#START')>0) EXIT
     end do

     ! Initialize array so all elements are set to approriate values
     Upstream_Data = 0.0
     Upstream_Data(:, Rho_)   = 1.0
     Upstream_Data(:, P_)     = inv_g
     do iFluid = IonFirst_+1, nFluid
        call select_fluid 
        Upstream_Data(:, Rho_)   = LowDensityRatio
        Upstream_Data(:, iP)     = inv_g*LowDensityRatio &
             *MassIon_I(1)/MassFluid_I(iFluid)
     end do
     if(UseMultiSpecies)then
        Upstream_Data(:, SpeciesFirst_) = &
             1.0 - cTiny8*(SpeciesLast_-SpeciesFirst_)
        Upstream_Data(:, SpeciesFirst_+1:SpeciesLast_)= &
             cTiny8
     end if
     
     ! Read the data
     do
        read(UNITTMP_,*,iostat=iError) &
             Tmp_Time, TmpData_V(1:nVarImf)
        write(*,*)'Upstream_Npts=',Upstream_Npts
        write(*,*) "Tmp_Time=", Tmp_Time
        write(*,*) "TmpData_V=",TmpData_V
        if (iError /= 0) EXIT

        if (Upstream_Npts >= Max_Upstream_Npts) then
           write(*,*) "=> This upstream boundary condition file"
           write(*,*) "   contains too many lines! (Max lines ="
           write(*,*) Max_Upstream_Npts,"). One method for getting"
           write(*,*) "   around this limitation is to have multiple"
           write(*,*) "   files which are read in at different times"
           write(*,*) "   during the run.  These files can overlap"
           write(*,*) "   in time."
           EXIT
        end if

        Upstream_Npts = Upstream_Npts + 1
        
        do iVar =1, nVarImf
           jVar=iVarImf_V(iVar)
           Upstream_Data(Upstream_Npts,jVar) = TmpData_V(iVar)              
        end do

        if(DoTest)then
           write(*,*)'Upstream_Npts=', Upstream_Npts
           write(*,*)'Upstream_Data=', Upstream_Data(Upstream_Npts,:)
        end if

        Upstream_integer_time(Upstream_Npts,:) = Tmp_Time

        if (Upstream_integer_time(Upstream_Npts,1) >= 65 .and. &
             Upstream_integer_time(Upstream_Npts,1) < 100)      &
             Upstream_integer_time(Upstream_Npts,1) =            &
             Upstream_integer_time(Upstream_Npts,1) + 1900

        if (Upstream_integer_time(Upstream_Npts,1) < 65)      &
             Upstream_integer_time(Upstream_Npts,1) =           &
             Upstream_integer_time(Upstream_Npts,1) + 2000
        
        if (UseZeroBx) Upstream_Data(Upstream_Npts,1) = 0.0

        if (TypeInputCoordSystem /= TypeCoordSystem) then 

           if (TypeCoordSystem == 'GSM' .and. &
                TypeInputCoordSystem == 'GSE') then 

              call CON_recalc(&
                   Upstream_integer_time(Upstream_Npts,1),&
                   Upstream_integer_time(Upstream_Npts,2),&
                   Upstream_integer_time(Upstream_Npts,3),&
                   Upstream_integer_time(Upstream_Npts,4),&
                   Upstream_integer_time(Upstream_Npts,5),&
                   Upstream_integer_time(Upstream_Npts,6))
              
              Upstream_Data(Upstream_Npts,Bx_:Bz_)=&
                   matmul(GsmGse_DD,&
                   Upstream_Data(Upstream_Npts,Bx_:Bz_))
              Upstream_Data(Upstream_Npts,Ux_:Uz_)=&
                   matmul(GsmGse_DD,&
                   Upstream_Data(Upstream_Npts,Ux_:Uz_))
           elseif(TypeCoordSystem == 'HGI' .and. &
                TypeInputCoordSystem == 'GSM') then
              if(DoSetMatrix)then
                 DoSetMatrix=.false.
                 call CON_recalc( &
                      iStartTime_I(1),& ! year
                      iStartTime_I(2),& ! month
                      iStartTime_I(3),& ! day
                      iStartTime_I(4),& ! hour
                      iStartTime_I(5),& ! minute
                      iStartTime_I(6))  ! second
                 HgiGsm_DD=&
                      matmul(HgiGse_DD,transpose(GsmGse_DD))
              end if
              Upstream_Data(Upstream_Npts,Bx_:Bz_)=&
                   matmul(HgiGsm_DD,&
                   Upstream_Data(Upstream_Npts,Bx_:Bz_))
              Upstream_Data(Upstream_Npts,Ux_:Uz_)=&
                   matmul(HgiGsm_DD,&
                   Upstream_Data(Upstream_Npts,Ux_:Uz_))
           else
              write(*,*) 'Transformation between input ',&
                   ' coordinate system=',TypeInputCoordSystem,&
                   ' and code system=',TypeCoordSystem,&
                   ' is not implemented'
              call stop_mpi('GM_ERROR')
           endif
        endif
     enddo

     close(UNITTMP_)

     do i=1,Upstream_Npts
        call time_int_to_real(Upstream_integer_Time(i,:),Upstream_Time(i))
        Upstream_Time(i) = Upstream_Time(i) + TimeDelay
     end do

     ! Check if the IMF data is within 1 day of the start time.
     if (minval(abs(Upstream_Time(1:Upstream_Npts)-StartTime)) &
          > 24.0*3600.0) then
        write(*,*) "**********************************************************"
        write(*,*) "*                                                        *"
        write(*,*) "*         Warning! Warning! Warning! Warning! Warning!   *"
        write(*,*) "*                                                        *"
        write(*,*) "*  Time dependent solar wind file disagrees with the     *"
        write(*,*) "*  starting time of the                                  *"
        write(*,*) "*  simulation by more than 24 hours.  This could cause   *"
        write(*,*) "*  results which you                                     *"
        write(*,*) "*  may not enjoy.                                        *"
        write(*,*) "*                                                        *"
        if (UseStrict) call stop_mpi('Correct PARAM.in or the IMF data file')
        
        write(*,*) "*  I am assuming that you are smarter than I am.         *"
        write(*,*) "*  Continuing with simulation.                           *"
        write(*,*) "*                                                        *"
        write(*,*) "**********************************************************"
     endif

     if(lVerbose>0)then
        call write_prefix; write(iUnitOut,*) NameSub, &
             ' read ',Upstream_Npts,' points from ',trim(upstreamfilename)
     end if

  endif

  call MPI_Bcast(Upstream_Npts,1,MPI_Integer,0,iComm,iError)
  if(iError>0)call stop_mpi( NameSub//": Upstream_Npts could not be broadcast")

  call MPI_Bcast(Upstream_Time,Max_Upstream_Npts,MPI_DOUBLE_PRECISION, &
       0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub//": Upstream_Time could not be broadcast")

  call MPI_Bcast(Upstream_Data,Max_Upstream_Npts*MaxVarImf,MPI_Real, &
       0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub//": Upstream_Data could not be broadcast")

  call MPI_Bcast(TypeInputCoordSystem,3,MPI_CHARACTER,0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub// &
       ": TypeInputCoordSystem could not be broadcast")

  call MPI_Bcast(Propagation_Plane_XY,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub// &
       ": Propagation_Plane_XY could not be broadcast")

  call MPI_Bcast(Propagation_Plane_XZ,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub// &
       ": Propagation_Plane_XZ could not be broadcast")

  call MPI_Bcast(Satellite_Y_Pos,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub// &
       ": Satellite_Y_Pos could not be broadcast")

  call MPI_Bcast(Satellite_Z_Pos,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(NameSub// &
       ": Satellite_Z_Pos could not be broadcast")

  ! This part will be removed by KC.
  if (SW_T_dim <= 0.0) then

     if (Upstream_Npts == 1) then

        SW_Bx_dim  = Upstream_Data(1,Bx_)
        SW_By_dim  = Upstream_Data(1,By_)
        SW_Bz_dim  = Upstream_Data(1,Bz_)
        
        SW_Ux_dim  = Upstream_Data(1,Ux_)
        SW_Uy_dim  = Upstream_Data(1,Uy_)
        SW_Uz_dim  = Upstream_Data(1,Uz_)

        SW_n_dim   = Upstream_Data(1,rho_)
        SW_T_dim   = Upstream_Data(1,p_)

     else

        ! Find the index with time before start time
        i = 1

        do while ((i < Upstream_Npts).and.   &
             (Upstream_Time(i) < StartTime))
           i = i + 1
        enddo

        if ((i == Upstream_Npts .and. &
             Upstream_Time(i) <= StartTime).or.(i==1)) then 

           SW_Bx_dim  = Upstream_Data(i,Bx_)
           SW_By_dim  = Upstream_Data(i,By_)
           SW_Bz_dim  = Upstream_Data(i,Bz_)

           SW_Ux_dim  = Upstream_Data(i,Ux_)
           SW_Uy_dim  = Upstream_Data(i,Uy_)
           SW_Uz_dim  = Upstream_Data(i,Uz_)

           SW_n_dim = Upstream_Data(i,rho_)
           SW_T_dim   = Upstream_Data(i,p_)

        else

           DtData2 = (Upstream_Time(i) - StartTime) / &
                (Upstream_Time(i) - Upstream_Time(i-1) + 1.0e-6)
           DtData1 = 1.0 - DtData2

           SW_Bx_dim  =       DtData1  * Upstream_Data(i,Bx_) + &
                DtData2 * Upstream_Data(i-1,Bx_)
           SW_By_dim  =       DtData1  * Upstream_Data(i,By_) + &
                DtData2 * Upstream_Data(i-1,By_)
           SW_Bz_dim  =       DtData1  * Upstream_Data(i,Bz_) + &
                DtData2 * Upstream_Data(i-1,Bz_)

           SW_Ux_dim  =       DtData1  * Upstream_Data(i,Ux_) + &
                DtData2 * Upstream_Data(i-1,Ux_)
           SW_Uy_dim  =       DtData1  * Upstream_Data(i,Uy_) + &
                DtData2 * Upstream_Data(i-1,Uy_)
           SW_Uz_dim  =       DtData1  * Upstream_Data(i,Uz_) + &
                DtData2 * Upstream_Data(i-1,Uz_)

           SW_n_dim =       DtData1  * Upstream_Data(i,rho_) + &
                DtData2 * Upstream_Data(i-1,rho_)
           SW_T_dim   =       DtData1  * Upstream_Data(i,p_) + &
                DtData2 * Upstream_Data(i-1,p_)
        endif
     endif
  endif

end subroutine read_upstream_input_file

!============================================================================

subroutine normalize_upstream_data

  use ModPhysics, ONLY: &
       Io2No_V, UnitTemperature_, UnitN_, UnitP_, UnitU_, UnitB_,&
       AverageIonCharge, ElectronTemperatureRatio
  use ModUpstreamData
  use ModMultiFluid, ONLY: MassIon_I
  use ModConst

  implicit none

  integer:: iVar, jVar
  integer:: T_= p_

  character(len=*), parameter:: NameSub = 'normalize_upstream_data'
  logical:: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call set_oktest(NameSub,DoTest,DoTestMe)

  ! normalize B and U 
  Upstream_Data(:,Bx_:Bz_) = Upstream_Data(:, Bx_:Bz_)*Io2No_V(UnitB_)
  Upstream_Data(:,Ux_:Uz_) = Upstream_Data(:, Ux_:Uz_)*Io2No_V(UnitU_)
  
  do iFluid = IonFirst_+1, nFluid
     call select_fluid
     Upstream_Data(:, iRho  )=Upstream_Data(:, iRho)*Io2No_V(UnitN_)
     Upstream_Data(:, iRhoUx)=Upstream_Data(:, iRhoUx)*Io2No_V(UnitU_)
     Upstream_Data(:, iRhoUy)=Upstream_Data(:, iRhoUy)*Io2No_V(UnitU_)
     Upstream_Data(:, iRhoUz)=Upstream_Data(:, iRhoUz)*Io2No_V(UnitU_)

     !change number density to mass density
     if (UseNumberDensity) &
          Upstream_Data(:,irho)=Upstream_Data(:,irho)*MassFluid_I(iFluid)

     !change temperature to pressure
     if(UseTemperature)  Upstream_Data(:,p_)= &
          Upstream_Data(:,iP)*Io2No_V(UnitTemperature_)&
          *Upstream_Data(:,irho)/MassFluid_I(iFluid)
  end do

  if(UseMultiSpecies) then
     Upstream_Data(:,rho_)=0.0

     !calculate number density to get the right pressure p=n*T
     do jVar=SpeciesFirst_, SpeciesLast_
        Upstream_Data(:,jVar)=Upstream_Data(:,jVar)*Io2No_V(UnitN_) 
        Upstream_Data(:,rho_)=Upstream_Data(:,rho_)+Upstream_Data(:,jVar)
     end do

     !normalize p = n*T 
     if(UseTemperature)  Upstream_Data(:,p_)= &
          Upstream_Data(:,T_)*Io2No_V(UnitTemperature_)*Upstream_Data(:,rho_)

     !calculate mass density of each ion species and the total mass density
     Upstream_Data(:,rho_)=0.0
     do jVar=SpeciesFirst_, SpeciesLast_
        Upstream_Data(:,jVar)=Upstream_Data(:,jVar)*MassSpecies_V(jVar) 
        Upstream_Data(:,rho_)=Upstream_Data(:,rho_)+Upstream_Data(:,jVar)
     end do

  else
     Upstream_Data(:,rho_) = &
          Upstream_Data(:, rho_)*Io2No_V(UnitN_)
     
     if (UseNumberDensity) &
          Upstream_Data(:,rho_)=Upstream_Data(:,rho_)*MassIon_I(1)

     if(UseTemperature)  Upstream_Data(:,p_)= &
          Upstream_Data(:,T_)*Io2No_V(UnitTemperature_)&
          *Upstream_Data(:,rho_)/MassIon_I(1)
     
  end if

  ! normalize pressure or temperature
  if(UseTemperature)then
     Upstream_Data(:,p_) = Upstream_Data(:,p_) * &
          (1.0 + AverageIonCharge*ElectronTemperatureRatio)
  else
     Upstream_Data(:, p_) = Upstream_Data(:, p_)*Io2No_V(UnitP_)
  end if

  if(DoTestMe)then
     write(*,*)'Io2No_V(UnitP_)',Io2No_V(UnitP_)
     write(*,*)'Io2No_V(UnitN_)*Io2No_V(UnitTemperature_)',&
          Io2No_V(UnitN_)*Io2No_V(UnitTemperature_)     
     write(*,*)'After normalization, Upstream_Data(1,1:MaxVarImf)=',&
          Upstream_Data(1,1:MaxVarImf)
     write(*,*)'After normalization, Upstream_Data(2,1:MaxVarImf)=',&
          Upstream_Data(2,1:MaxVarImf)
     write(*,*)'After normalization, Upstream_Data(3,1:MaxVarImf)=',&
          Upstream_Data(3,1:MaxVarImf)
  end if

end subroutine normalize_upstream_data

!======================================================================

subroutine get_solar_wind_point(TimeSimulation, x, y, z, SolarWind_V)

  use ModKind
  use ModMain
  use ModUpstreamData
  use ModPhysics
  use ModVarIndexes
  use ModGeometry, ONLY: x1, x2

  implicit none
  real, intent(in) :: TimeSimulation,x,y,z
  real, intent(out) :: SolarWind_V(MaxVarImf) ! Varying solar wind parameters

  integer :: iData
  real(Real8_) :: DtData1, DtData2, time_now
  real    :: Satellite_X_Pos
  real    :: SatDistance_D(3), U_D(3)
  !--------------------------------------------------------------------------

  time_now = StartTime+TimeSimulation

  if (UseUpstreamInputFile.and.(Upstream_Npts > 0)) then
     
     if (Upstream_Npts == 1) then

        SolarWind_V = Upstream_Data(1,:)
        
     else
        if(Upstream_Data(1,Ux_)>0.0)then
           Satellite_X_Pos=x1  
        else
           Satellite_X_Pos=x2
        endif
        
        iData = 1
        do while ((iData < Upstream_Npts).and. &
             (Upstream_Time(iData) < time_now))
           iData = iData + 1
        enddo
        
        if ((abs(x-Satellite_X_Pos) > 1.0e-5) .or.   &
             (Propagation_Plane_XY /= 0.0)    .or.   &
             (Propagation_Plane_XZ /= 0.0)) then

           SatDistance_D = &
                (/x-Satellite_X_Pos, y-Satellite_Y_Pos, z-Satellite_Z_Pos /)

           if ((iData == Upstream_Npts .and. &
                Upstream_Time(iData) <= time_now).or.(iData==1)) then 
              U_D = Upstream_Data(iData,Ux_:Uz_)
           else           
              DtData2 = (Upstream_Time(iData) - time_now) / &
                   (Upstream_Time(iData) - Upstream_Time(iData-1) + 1.0e-6)
              DtData1 = 1.0 - DtData2
              
              U_D = DtData1 * Upstream_Data(iData,  Ux_:Uz_) + &
                   DtData2 * Upstream_Data(iData-1, Ux_:Uz_)

           endif
           
           ! Time is in SI units
           time_now = time_now - &
                dot_product(SatDistance_D, Normal_D)/ &
                dot_product(U_D, Normal_D) * No2Si_V(UnitT_)
           
           iData = 1
           do while ((iData < Upstream_Npts).and. &
                (Upstream_Time(iData) < time_now))
              iData = iData + 1
           enddo
        end if
        
        if ((iData == Upstream_Npts .and. &
             Upstream_Time(iData) <= time_now).or.(iData==1)) then 

           SolarWind_V(:)  = Upstream_Data(iData,:)

        else

           DtData2 = (Upstream_Time(iData) - time_now) / &
                (Upstream_Time(iData) - Upstream_Time(iData-1) + 1.0e-6)
           DtData1 = 1.0 - DtData2

           SolarWind_V(:)  = DtData1 * Upstream_Data(iData,:) + &
                DtData2 * Upstream_Data(iData-1,:)

        endif

     endif
  else
     SolarWind_V = FaceState_VI(:,east_)
  endif

end subroutine get_solar_wind_point
