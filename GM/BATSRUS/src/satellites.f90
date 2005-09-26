!^CFG COPYRIGHT UM
!==============================================================================
subroutine read_satellite_input_files
  use ModProcMH
  use ModMain, ONLY : nDim, lVerbose, TypeCoordSystem, StartTime
  use ModIO
  use ModPhysics, ONLY : unitSI_x
  use CON_axes, ONLY: transform_matrix
  use ModTimeConvert, ONLY: time_int_to_real
  use ModMpi
  use ModKind
  implicit none

  integer :: iError, i, iSat , nPoint

  ! One line of input
  character (len=100) :: line

  integer      :: iTime_I(7)
  real         :: Xyz_D(nDim)
  real(Real8_) :: DateTime
  real         :: Time_I(Max_Satellite_Npts)
  real         :: Xyz_DI(nDim, Max_Satellite_Npts)

  character(len=*), parameter :: NameSub = 'read_satellite_input_files'

  logical :: DoTest, DoTestMe

  !---------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  SATELLITES: do iSat=1, nSatellite

     if(.not.UseSatelliteFile(iSat)) CYCLE SATELLITES

     ! Read file on the root processor
     if (iProc == 0) then

        filename = Satellite_Name(iSat)

        if(lVerbose>0)then
           call write_prefix; write(iUnitOut,*) NameSub, &
                " reading: ",trim(filename)
        end if

        open(unit_tmp, file=filename, status="old", iostat = iError)

        if (iError /= 0) call stop_mpi(NameSub // &
             ' ERROR: unable to open file ' // trim(filename))

        nPoint = 0

        ! Read the file: read #COOR TypeCoord, #START and points
        ! Default coordinate system is the one used by BATSRUS (or GSM?)
        TypeSatCoord_I(iSat) = TypeCoordSystem
        READFILE: do

           read(unit_tmp,'(a)', iostat = iError ) line

           if (iError /= 0) EXIT READFILE

           if(index(line,'#COOR')>0) &
                read(unit_tmp,'(a)') TypeSatCoord_I(iSat)
           
           if(index(line,'#START')>0)then

              READPOINTS: do

                 read(unit_tmp,*, iostat=iError) iTime_I, Xyz_D

                 if (iError /= 0) EXIT READFILE

                 if (nPoint >= Max_Satellite_Npts) then
                    call write_prefix;
                    write(*,*) NameSub,' WARNING: trajectory file: ',&
                         trim(filename),' contains too many lines! '
                    call write_prefix; write(*,*) NameSub, &
                         ': max number of lines =',Max_Satellite_Npts
                    EXIT READFILE
                 endif

                 ! Add new point
                 nPoint = nPoint + 1

                 ! Store coordinates
                 Xyz_DI(:,nPoint) = Xyz_D

                 ! Convert integer date/time to simulation time
                 call time_int_to_real(iTime_I, DateTime)
                 Time_I(nPoint) = DateTime - StartTime

              enddo READPOINTS

           endif

        enddo READFILE

        close(unit_tmp)

        if(DoTest)write(*,*) NameSub,': nPoint=',nPoint

        ! Convert the coordinates if necessary
        if(TypeSatCoord_I(iSat) /= TypeCoordSystem)then
           do i = 1, nPoint
              Xyz_DI(:,i) = matmul( &
                   transform_matrix( Time_I(i), &
                   TypeSatCoord_I(iSat), TypeCoordSystem), Xyz_DI(:,i) )
           end do
        end if

     end if

     ! Tell the number of points to the other processors
     call MPI_Bcast(nPoint, 1, MPI_INTEGER, 0, iComm, iError)
     Satellite_Npts(iSat) = nPoint

     ! Tell the other processors the satellite time
     call MPI_Bcast(Time_I, nPoint, MPI_REAL, 0, iComm, iError)

     ! Tell the other processors the coordinates
     call MPI_Bcast(Xyz_DI, nDim*nPoint, MPI_REAL, 0, iComm, iError)

     ! Store time and positions for satellite iSat on all PE-s

     Satellite_Time(iSat, 1:nPoint) = Time_I(1:nPoint)
     do i = 1, nPoint
        xSatellite_traj(iSat, i, :) = Xyz_DI(:, i)
     end do

     if(DoTest)then
        nPoint = min(10,nPoint)
        write(*,*) NameSub,': tSat=', Satellite_Time( iSat, 1:nPoint)
        write(*,*) NameSub,': xSat=', xSatellite_traj(iSat, 1:nPoint,1)
        write(*,*) NameSub,': ySat=', xSatellite_traj(iSat, 1:nPoint,2)
        write(*,*) NameSub,': zSat=', xSatellite_traj(iSat, 1:nPoint,3)
     end if

  end do SATELLITES

end subroutine read_satellite_input_files

!==========================================================================

subroutine set_satellite_flags
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBlockMax,PROCtest,unusedBLK
  use ModGeometry, ONLY : XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModGeometry, ONLY : TypeGeometry               !^CFG IF COVARIANT
  use ModIO, ONLY : iBLKsatellite,iPEsatellite,Xsatellite,&
       SatelliteInBLK,DoTrackSatellite_I,nSatellite
  use ModNumConst
  use ModMpi
  implicit none

  integer :: isat, iPE,iPEtmp, iBLK, iBLKtemp
  real    :: XSat,YSat,ZSat
  integer :: i,j,k, iError
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  call set_oktest('set_satellite_flags',oktest, oktest_me)

  if (oktest_me) &
       write(*,*)'Starting set_satellite_flags',&
       nSatellite,', call set_satellite_positions'

  call set_satellite_positions

  do iSat=1, nSatellite
     if(.not.DoTrackSatellite_I(iSat))CYCLE !Position is not defined
     select case(TypeGeometry)           !^CFG IF COVARIANT
     case('cartesian')                   !^CFG IF COVARIANT
        xSat=XSatellite(iSat,1)
        ySat=XSatellite(iSat,2)
        zSat=XSatellite(iSat,3)
     case('spherical')                   !^CFG IF COVARIANT BEGIN
        call xyz_to_spherical(XSatellite(iSat,1), XSatellite(iSat,2),&
             XSatellite(iSat,3), xSat, ySat, zSat)
     case('spherical_lnr')                   
        call xyz_to_spherical(XSatellite(iSat,1), XSatellite(iSat,2),&
             XSatellite(iSat,3), XSat, YSat, ZSat)
        xSat=log(max(xSat, cTiny))
     case default
        call stop_mpi('Unknown TypeGeometry='//TypeGeometry)
     end select                          !^CFG END COVARIANT

     iPE = -1
     iBLKtemp = -1

     do iBLK = 1, nBlockMax

        SatelliteInBLK(isat,iBLK) =.not.unusedBLK(iBLK).and.&
             XSat >  XyzStart_BLK(1,iBLK) - cHalf*dx_BLK(iBLK) .and. &
             XSat <= XyzStart_BLK(1,iBLK) + (nI-cHalf)*dx_BLK(iBLK) .and. &
             YSat >  XyzStart_BLK(2,iBLK) - cHalf*dy_BLK(iBLK) .and. &
             YSat <= XyzStart_BLK(2,iBLK) + (nJ-cHalf)*dy_BLK(iBLK) .and. &
             ZSat >  XyzStart_BLK(3,iBLK) - cHalf*dz_BLK(iBLK) .and. &
             ZSat <= XyzStart_BLK(3,iBLK) + (nK-cHalf)*dz_BLK(iBLK)

   	if(SatelliteInBLK(isat,iBLK))then 
           iPE = iProc
	   iBLKtemp = iBLK
        end if

     end do

     call MPI_ALLREDUCE(iPE, iPEtmp, 1, MPI_INTEGER, MPI_MAX, iComm, iError)
     iPEsatellite(isat) = iPEtmp

     if (iPEsatellite(isat) >= 0) &
          call MPI_Bcast(iBLKtemp ,1,MPI_INTEGER,iPEsatellite(isat),iComm,iError)
     iBLKsatellite(isat) = iBLKtemp

     if (iPEsatellite(isat) == -1) DoTrackSatellite_I(isat) = .false.

     if (oktest_me) write(*,*)'set_satellite_flags (Proc',PROCtest,')(isat=', &
          isat,'): iPE,iBLK,TrackSatellite:', &
          iPEsatellite(isat),iBLKsatellite(isat),DoTrackSatellite_I(isat) 

  end do

end subroutine set_satellite_flags


!=============================================================================

subroutine set_satellite_positions
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,StartTime,time_simulation
  use ModGeometry, ONLY : x1,x2,y1,y2,z1,z2,XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModGeometry, ONLY : TypeGeometry               !^CFG IF COVARIANT
  use ModNumConst
  use ModIO
  implicit none


  integer :: i, iSat, iBLK
  real    :: XSat,YSat,ZSat
  real    :: dtime

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  if (iProc==0) &
       call set_oktest('set_satellite_positions',oktest, oktest_me)

  do iSat = 1, nSatellite

     if (UseSatelliteFile(iSat)) then

        if (Satellite_Npts(iSat) > 0) then

           i = icurrent_satellite_position(iSat)

           do while ((i < Satellite_Npts(iSat)) .and.   &
                (Satellite_Time(iSat,i) < Time_Simulation))
              i = i + 1
           enddo

           icurrent_satellite_position(iSat) = i

           if ((i == Satellite_Npts(iSat).and. &
                Satellite_Time(iSat,i) <= Time_Simulation ).or.(i==1)) then 

              DoTrackSatellite_I(iSat) = .false.

           else

              DoTrackSatellite_I(iSat) = .true.

              dTime = 1.0 - (Satellite_Time(iSat,i) - Time_Simulation) / &
                   (Satellite_Time(iSat,i) - Satellite_Time(iSat,i-1) + 1.0e-6)

              xSatellite(iSat,:) = dTime * xSatellite_traj(iSat,i,:) + &
                   (1.0 - dTime) * xSatellite_traj(iSat,i-1,:) 

           endif

        endif

     else 

        call satellite_trajectory_formula(iSat)

     end if

  end do

end subroutine set_satellite_positions


!=============================================================================

subroutine satellite_trajectory_formula(iSat)
  use ModIO, ONLY : DoTrackSatellite_I,XSatellite,Satellite_name
  implicit none

  integer, intent(in) :: iSat
  character (len=100) :: name_string
  real :: Xvect(3)

  name_string = trim(Satellite_name(iSat))
  Xvect(:) = XSatellite(iSat,:)

  ! Case should be for a specific satellite.  The trajectories can depend
  ! on the 'real' time so that the satellite knows where it is at.  For
  ! example, Cassini could be if'd on the date so that the code knows 
  ! whether the run is for a time near Earth, Jupiter or Saturn. 

  ! This routine should always set the TrackSatellite(iSat) flag. When
  ! the satellite is at a useless position the time should return a
  ! do not track flag (DoTrackSatellite_I(iSat) = .false.).


  select case(name_string)
  case ('earth')
     Xvect(1) = 5.0
     Xvect(2) = 5.0
     Xvect(3) = 5.0
     DoTrackSatellite_I(iSat) = .true.
  case ('cassini')
     Xvect(1) = 5.0
     Xvect(2) = 5.0
     Xvect(3) = 5.0
     DoTrackSatellite_I(iSat) = .true.
  case default 
     Xvect(1) = 1.0
     Xvect(2) = 1.0
     Xvect(3) = 1.0
     DoTrackSatellite_I(iSat) = .false.
  end select

end subroutine satellite_trajectory_formula

!=============================================================================

subroutine open_satellite_output_files

  use ModMain,   ONLY: n_step
  use ModIoUnit, ONLY: io_unit_new
  use ModIO,     ONLY: nSatellite, Satellite_name, filename,&
       NamePlotDir, iUnitSat_I
  implicit none

  integer :: iSat, l1, l2
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  call set_oktest('open_satellite_output_files', oktest, oktest_me)

  do iSat = 1, nSatellite
     l1 = index(Satellite_name(iSat), '/', back=.true.) + 1
     l2 = index(Satellite_name(iSat), '.') - 1
     if (l1-1<=0) l1=1
     if (l2+1<=0) l2=len_trim(Satellite_name(iSat))

     write(filename,'(a,i6.6,a)')trim(NamePlotDir)//&
          'sat_'//Satellite_Name(iSat)(l1:l2)//'_n',n_step,'.sat'

     if(oktest) then
        write(*,*) 'open_satellite_output_files: satellitename:', &
          Satellite_name(iSat)
        write(*,*) 'iSat,l1,l2: ',iSat,l1,l2
        write(*,*) 'open_satellite_output_files: filename:', filename
     end if

     iUnitSat_I(iSat)=io_unit_new()
     open(iUnitSat_I(iSat), file=filename, status="replace")
  end do

end subroutine open_satellite_output_files

!=============================================================================

subroutine close_satellite_output_files

  use ModIO, ONLY : nSatellite,iUnitSat_I
  implicit none

  integer :: iSat

  do iSat = 1, nSatellite
     close(iUnitSat_I(iSat))
  end do

end subroutine close_satellite_output_files
!=============================================================================

