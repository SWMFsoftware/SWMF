Module ModUpstreamData
  integer, parameter :: Max_Upstream_Npts = 5000
  integer :: Upstream_Npts

  real, dimension(Max_Upstream_Npts, 8)   :: Upstream_Data
  real*8, dimension(Max_Upstream_Npts)    :: Upstream_Time

  ! Solar Wind Input Variables
  character(len=3) :: TypeInputCoordSystem = 'GSM'
  real :: Propagation_Plane_XY, Propagation_Plane_XZ
end Module ModUpstreamData

subroutine read_upstream_input_file(upstreamfilename)
  use ModProcMH
  use ModMain
  use ModIoUnit, ONLY : UNITTMP_
  use ModPhysics
  use ModNumConst
  use ModUpstreamData
  use CON_geopack
  use ModIO, ONLY: iUnitOut, write_prefix
  use CON_physics, ONLY: time_int_to_real
  use ModUtilities, ONLY: upper_case
  use ModMpi
  implicit none

  character (len=100), intent(in) :: upstreamfilename

  integer :: iError, i 
  logical :: done, zerobx

  ! One line of input
  character (len=100) :: line

  integer, dimension(1:7) :: Tmp_Time
  integer, dimension(Max_Upstream_Npts,7) :: Upstream_integer_time

  real :: timedelay
  real :: dt_min
  real,save::HgiGsm_DD(3,3)
  logical,save::DoSetMatrix=.true.

  !--------------------------------------------------------------------------
  zerobx = .false.
  timedelay = cZero

  if (iProc == 0) then

     open(UNITTMP_, file=upstreamfilename, status="old", iostat = iError)

     if (iError /= 0) &
          call stop_mpi("Unable to open file "//trim(upstreamfilename))

     done = .false.
     Upstream_Npts = 0
     Propagation_Plane_XY = 0.0
     Propagation_Plane_XZ = 0.0

     if(lVerbose>0)then
        call write_prefix; write(iUnitOut,*) &
             "=> Reading Upstream Solar Wind Conditions File."
     end if

     do while (.not.done)

        read(UNITTMP_,'(a)', iostat = iError ) line

        if (iError.ne.0) done = .true.

        if(index(line,'#COOR')>0)then
           read(UNITTMP_,'(a)') TypeInputCoordSystem
           call upper_case(TypeInputCoordSystem)
        endif

        if(index(line,'#PLANE')>0)then
           read(UNITTMP_,*) Propagation_Plane_XY
           read(UNITTMP_,*) Propagation_Plane_XZ
           Propagation_Plane_XY = Propagation_Plane_XY * cPi / 180.0
           Propagation_Plane_XZ = Propagation_Plane_XZ * cPi / 180.0
        endif

        if(index(line,'#POSITION')>0)then
           read(UNITTMP_,*) Satellite_Y_Pos
           read(UNITTMP_,*) Satellite_Z_Pos
        endif

        if(index(line,'#ZEROBX')>0)then
           read(UNITTMP_,*) zerobx
        endif

        if(index(line,'#TIMEDELAY')>0)then
           read(UNITTMP_,*) timedelay
        endif


        if(index(line,'#START')>0)then

           do while (.not.done)

              Upstream_Npts = Upstream_Npts + 1

              read(UNITTMP_,*,iostat=iError) &
                   (Tmp_Time(i),i=1,7), &
                   (Upstream_Data(Upstream_Npts,i),i=1,8)

              if (iError.ne.0) then
                 done = .true.
                 Upstream_Npts = Upstream_Npts - 1
              else
                 if (Upstream_Npts >= Max_Upstream_Npts) then
                    done = .true.
                    write(*,*) "=> This upstream boundary condition file"
                    write(*,*) "   contains too many lines! (Max lines ="
                    write(*,*) Max_Upstream_Npts,"). One method for getting"
                    write(*,*) "   around this limitation is to have multiple"
                    write(*,*) "   files which are read in at different times"
                    write(*,*) "   during the run.  These files can overlap"
                    write(*,*) "   in time."
                 else

                    Upstream_integer_time(Upstream_Npts,:) = Tmp_Time(:)

                    if (Upstream_integer_time(Upstream_Npts,1) >= 65 .and. &
                         Upstream_integer_time(Upstream_Npts,1) < 100)      &
                         Upstream_integer_time(Upstream_Npts,1) =            &
                         Upstream_integer_time(Upstream_Npts,1) + 1900

                    if (Upstream_integer_time(Upstream_Npts,1) < 65)      &
                         Upstream_integer_time(Upstream_Npts,1) =           &
                         Upstream_integer_time(Upstream_Npts,1) + 2000

                    if (zerobx) Upstream_Data(Upstream_Npts,1) = 0.0

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

                          Upstream_Data(Upstream_Npts,1:3)=&
                               matmul(GsmGse_DD,&
                               Upstream_Data(Upstream_Npts,1:3))
                          Upstream_Data(Upstream_Npts,4:6)=&
                               matmul(GsmGse_DD,&
                               Upstream_Data(Upstream_Npts,4:6))
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
                          Upstream_Data(Upstream_Npts,1:3)=&
                               matmul(HgiGsm_DD,&
                               Upstream_Data(Upstream_Npts,1:3))
                          Upstream_Data(Upstream_Npts,4:6)=&
                               matmul(HgiGsm_DD,&
                               Upstream_Data(Upstream_Npts,4:6))
                       else
                          write(*,*) 'Transformation between input ',&
                               ' coordinate system=',TypeInputCoordSystem,&
                               ' and code system=',TypeCoordSystem,&
                               ' is not implemented'
                          call stop_mpi('GM_ERROR')
                       endif
                    endif
                 endif
              endif
           enddo
        endif
     enddo

     close(UNITTMP_)

  endif

  call MPI_Bcast(Upstream_Npts,1,MPI_Integer,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Upstream_Npts could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(Upstream_integer_time,Max_Upstream_Npts*7,MPI_Integer, &
       0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Upstream_integer_Time could not be broadcast by read_upstream_input_file")

  dt_min = 366.0*24.0*3600.0
  do i=1,Upstream_Npts
     call time_int_to_real(Upstream_integer_Time(i,:),Upstream_Time(i))
     if (dt_min > abs(Upstream_Time(i) - StartTime)) then
        dt_min = abs(Upstream_Time(i) - StartTime)
     endif
  end do

  if ((dt_min > 24.0*3600.0).and.(iProc == 0)) then
     write(*,*) "******************************************************************************"
     write(*,*) "*                                                                            *"
     write(*,*) "*  Warning! Warning! Warning! Warning! Warning! Warning! Warning! Warning!   *"
     write(*,*) "*                                                                            *"
     write(*,*) "*  Time dependent solar wind file disagrees with the starting time of the    *"
     write(*,*) "*  simulation by more than 24 hours.  This could cause results which you     *"
     write(*,*) "*  may not enjoy.                                                            *"
     write(*,*) "*                                                                            *"
     if (UseStrict) call stop_mpi('Correct PARAM.in or the IMF data file')

     write(*,*) "*  I am assuming that you are smarter than I am. Continuing with simulation. *"
     write(*,*) "*                                                                            *"
     write(*,*) "******************************************************************************"
  endif

  call MPI_Bcast(Upstream_Data,Max_Upstream_Npts*8,MPI_Real, &
       0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Upstream_Data could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(TypeInputCoordSystem,3,MPI_CHARACTER,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "TypeInputCoordSystem could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(Propagation_Plane_XY,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Propagation_Plane_XY could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(Propagation_Plane_XZ,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Propagation_Plane_XZ could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(Satellite_Y_Pos,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Satellite_Y_Pos could not be broadcast by read_upstream_input_file")

  call MPI_Bcast(Satellite_Z_Pos,1,MPI_Real,0,iComm,iError)
  if(iError>0)call stop_mpi(&
       "Satellite_Z_Pos could not be broadcast by read_upstream_input_file")

  !call MPI_Bcast(Time_Array,7,MPI_Integer, 0,iComm,iError)
  !if(iError>0)call stop_mpi(&
  !     "Time_Array could not be broadcast by read_upstream_input_file")

  if (SW_T_dim <= 0.0) then

     if (Upstream_Npts == 1) then

        SW_Bx_dim  = Upstream_Data(1,1)
        SW_By_dim  = Upstream_Data(1,2)
        SW_Bz_dim  = Upstream_Data(1,3)

        SW_Ux_dim  = Upstream_Data(1,4)
        SW_Uy_dim  = Upstream_Data(1,5)
        SW_Uz_dim  = Upstream_Data(1,6)

        SW_rho_dim = Upstream_Data(1,7)
        SW_T_dim   = Upstream_Data(1,8)

     else

        i = 1

        do while ((i < Upstream_Npts).and.   &
             (Upstream_Time(i) < StartTime))
           i = i + 1
        enddo

        if ((i == Upstream_Npts .and. &
             Upstream_Time(i) <= StartTime).or.(i==1)) then 

           SW_Bx_dim  = Upstream_Data(i,1)
           SW_By_dim  = Upstream_Data(i,2)
           SW_Bz_dim  = Upstream_Data(i,3)

           SW_Ux_dim  = Upstream_Data(i,4)
           SW_Uy_dim  = Upstream_Data(i,5)
           SW_Uz_dim  = Upstream_Data(i,6)

           SW_rho_dim = Upstream_Data(i,7)
           SW_T_dim   = Upstream_Data(i,8)

        else

           dt = 1.0 - (Upstream_Time(i) - StartTime) / &
                (Upstream_Time(i) - Upstream_Time(i-1) + 1.0e-6)

           SW_Bx_dim  =       dt  * Upstream_Data(i,1) + &
                (1.0 - dt) * Upstream_Data(i-1,1)
           SW_By_dim  =       dt  * Upstream_Data(i,2) + &
                (1.0 - dt) * Upstream_Data(i-1,2)
           SW_Bz_dim  =       dt  * Upstream_Data(i,3) + &
                (1.0 - dt) * Upstream_Data(i-1,3)

           SW_Ux_dim  =       dt  * Upstream_Data(i,4) + &
                (1.0 - dt) * Upstream_Data(i-1,4)
           SW_Uy_dim  =       dt  * Upstream_Data(i,5) + &
                (1.0 - dt) * Upstream_Data(i-1,5)
           SW_Uz_dim  =       dt  * Upstream_Data(i,6) + &
                (1.0 - dt) * Upstream_Data(i-1,6)

           SW_rho_dim =       dt  * Upstream_Data(i,7) + &
                (1.0 - dt) * Upstream_Data(i-1,7)
           SW_T_dim   =       dt  * Upstream_Data(i,8) + &
                (1.0 - dt) * Upstream_Data(i-1,8)
        endif
     endif
  endif

end subroutine read_upstream_input_file
!============================================================================
!BOP
!IROUTINE: reusable procedure to normalize the physicsl data file
!INTERFACE
subroutine normalize_upstream_data
  !USES:
  use ModPhysics,ONLY:unitSI_rho,unitSI_p,unitSI_U,unitSI_B
  use ModUpstreamData, ONLY: Upstream_Data
  use ModConst
!EOP
  implicit none
!BOP
!DESCRIPTION:
!     According to conventions, temperature in the IMF files is given
!     in K, magnetic field in nanoTesla's, velocity in km/sec:
  real,parameter:: cNanoTesla=cOne/cE9 ![Tesla]
  real,parameter:: cKmPerSec=cThousand ![m/s]
! Concentration is related per cubic centimeter, density is
! the concentration times the proton mass
  real,parameter::cPerCm3=cE6 ![per cubic meter]
  real,parameter::ConcentrationToDensitySI=cProtonMass*cPerCm3
!EOP
  ! Normalize B
  Upstream_Data(:,1:3) = &
       Upstream_Data(:,1:3)*cNanotesla/unitSI_B
  ! Normalize U
  Upstream_Data(:,4:6) = &
       Upstream_Data(:,4:6)*cKmPerSec/unitSI_U
  ! Convert T into Normalized P:p=n K_B T
  Upstream_Data(:,8) = &
       (Upstream_Data(:,7)*cPerCm3)*cBoltzmann*&
       Upstream_Data(:,8)/unitSI_p
  ! Normalize rho
  Upstream_Data(:,7) = Upstream_Data(:,7)*&
       ConcentrationToDensitySI/unitSI_rho
  
end subroutine normalize_upstream_data
!======================================================================
subroutine get_solar_wind_point(TimeSimulation,y,z,&
     current_SW_rho, current_SW_Ux, current_SW_Uy, current_SW_Uz, &
     current_SW_Bx, current_SW_By, current_SW_Bz, current_SW_p)

  use ModKind
  use ModMain
  use ModUpstreamData
  use ModPhysics

  implicit none
  real, intent(in) :: TimeSimulation,y,z
  real, intent(out) :: & ! Varying solar wind parameters
       current_SW_rho, current_SW_Ux, current_SW_Uy, current_SW_Uz, &
       current_SW_Bx, current_SW_By, current_SW_Bz,current_SW_p


  integer :: iData
  real(Real8_) :: dtData1, dtData2, time_now
  real    :: Ux

  !--------------------------------------------------------------------------
  time_now = StartTime+TimeSimulation

  if (UseUpstreamInputFile.and.(Upstream_Npts > 0)) then

     if (Upstream_Npts == 1) then

        current_SW_Bx  = Upstream_Data(1,1)
        current_SW_By  = Upstream_Data(1,2)
        current_SW_Bz  = Upstream_Data(1,3)

        current_SW_Ux  = Upstream_Data(1,4)
        current_SW_Uy  = Upstream_Data(1,5)
        current_SW_Uz  = Upstream_Data(1,6)

        current_SW_rho = Upstream_Data(1,7)
        current_SW_p   = Upstream_Data(1,8)

     else

        iData = 1

        do while ((iData < Upstream_Npts).and.   &
             (Upstream_Time(iData) < time_now))
           iData = iData + 1
        enddo

        if ((Propagation_Plane_XY /= 0.0).or.   &
             (Propagation_Plane_XZ /= 0.0)) then

           if ((iData == Upstream_Npts .and. &
                Upstream_Time(iData) <= time_now).or.(iData==1)) then 
              Ux = Upstream_Data(iData,4) * unitSI_U
           else

              dtData2 = (Upstream_Time(iData) - time_now) / &
                   (Upstream_Time(iData) - Upstream_Time(iData-1) + 1.0e-6)
              dtData1 = 1.0 - dtData2

              Ux = dtData1  * Upstream_Data(iData,4) + &
                   dtData2 * Upstream_Data(iData-1,4)

              Ux = Ux * unitSI_U

           endif

           time_now = time_now + &
                ((Y - Satellite_Y_Pos) * sin(Propagation_Plane_XY) + &
                (Z - Satellite_Z_Pos) * sin(Propagation_Plane_XZ)) &
                * unitSI_x / Ux

           iData = 1

           do while ((iData < Upstream_Npts).and.   &
                (Upstream_Time(iData) < time_now))
              iData = iData + 1
           enddo

        endif

        if ((iData == Upstream_Npts .and. &
             Upstream_Time(iData) <= time_now).or.(iData==1)) then 

           current_SW_Bx  = Upstream_Data(iData,1)
           current_SW_By  = Upstream_Data(iData,2)
           current_SW_Bz  = Upstream_Data(iData,3)

           current_SW_Ux  = Upstream_Data(iData,4)
           current_SW_Uy  = Upstream_Data(iData,5)
           current_SW_Uz  = Upstream_Data(iData,6)

           current_SW_rho = Upstream_Data(iData,7)
           current_SW_p   = Upstream_Data(iData,8)

        else

           dtData2 = (Upstream_Time(iData) - time_now) / &
                (Upstream_Time(iData) - Upstream_Time(iData-1) + 1.0e-6)
           dtData1 = 1.0 - dtData2

           current_SW_Bx  = dtData1 * Upstream_Data(iData,1) + &
                dtData2 * Upstream_Data(iData-1,1)
           current_SW_By  = dtData1 * Upstream_Data(iData,2) + &
                dtData2 * Upstream_Data(iData-1,2)
           current_SW_Bz  = dtData1 * Upstream_Data(iData,3) + &
                dtData2 * Upstream_Data(iData-1,3)

           current_SW_Ux  = dtData1 * Upstream_Data(iData,4) + &
                dtData2 * Upstream_Data(iData-1,4)
           current_SW_Uy  = dtData1 * Upstream_Data(iData,5) + &
                dtData2 * Upstream_Data(iData-1,5)
           current_SW_Uz  = dtData1 * Upstream_Data(iData,6) + &
                dtData2 * Upstream_Data(iData-1,6)

           current_SW_rho = dtData1 * Upstream_Data(iData,7) + &
                dtData2 * Upstream_Data(iData-1,7)
           current_SW_p   = dtData1 * Upstream_Data(iData,8) + &
                dtData2 * Upstream_Data(iData-1,8)

        endif

     endif

  else

     current_SW_rho = SW_rho
     current_SW_p   = SW_p
     current_SW_Ux  = SW_Ux
     current_SW_Uy  = SW_Uy
     current_SW_Uz  = SW_Uz
     current_SW_Bx  = SW_Bx
     current_SW_By  = SW_By
     current_SW_Bz  = SW_Bz

  endif

end subroutine get_solar_wind_point
