!^CFG COPYRIGHT UM
!==============================================================================
subroutine read_satellite_input_files
  use ModProcMH
  use ModMain, ONLY : lVerbose
  use ModIO
  use ModPhysics, ONLY : unitSI_x
  use CON_physics, ONLY: time_int_to_real
  use ModMpi
  implicit none

  integer :: ierror, i,isat , Npts_tmp
  logical :: done

  ! One line of input
  character (len=100) :: line

  integer, dimension(Max_Satellite_Npts,7) :: T_tmp
  real, dimension(Max_Satellite_Npts,3) :: X_tmp
!  real*8, external :: Time_of_Year

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('read_satellite_input_files',oktest, oktest_me)

  do iSat=1,nSatellite

     if (iProc == 0) then

        if (UseSatelliteFile(isat)) then

           filename = Satellite_Name(isat)

           if(lVerbose>0)then
              call write_prefix; write(iUnitOut,*) &
                   "=> Reading Satellite Trajectory File = ",filename
           end if

           open(unit_tmp, file=filename, status="old", iostat = ierror)

           if (ierror.ne.0) then
              call stop_mpi("Satellites: Unable to open file "//filename//". Stopping.")
           endif

           done = .false.
           Npts_tmp = 0

           do while (.not.done)

              read(unit_tmp,'(a)', iostat = ierror ) line

              if (ierror.ne.0) done = .true.

              if(index(line,'#COOR')>0) &
                   read(unit_tmp,'(a)') TypeSatCoord_I(iSat)

              if(index(line,'#START')>0)then

                 do while (.not.done)

                    Npts_tmp = Npts_tmp + 1

                    read(unit_tmp,*,iostat=ierror) &
                         (T_tmp(Npts_tmp,i),i=1,7), &
                         (X_tmp(Npts_tmp,i),i=1,3)

                    if (ierror.ne.0) then
                       done = .true.
                       Npts_tmp = Npts_tmp - 1
                    else
                       if (Npts_tmp >= Max_Satellite_Npts) then
                          done = .true.
                          call write_prefix;
                          write(*,*) '=> Satellite trajectory file: ',&
                               trim(Satellite_Name(isat)),&
                               ' contains too many lines! '
                          call write_prefix; write(*,*) &
                               '(Max lines =',Max_Satellite_Npts,').'
                       endif
                    endif

                 enddo

              endif

           enddo

           close(unit_tmp)

        end if
     end if

     call MPI_Bcast(Npts_tmp,1,MPI_Integer,0,iComm,iError)
     if(iError>0)call stop_mpi(&
  	  "Satellite_Npts could not be broadcast by read_satellite_input_files")
     Satellite_Npts(isat) = Npts_tmp         

     call MPI_Bcast(TypeSatCoord_I(iSat),3,MPI_CHARACTER,0,iComm,iError)
     if(iError>0)call stop_mpi(&
  	  "TypeSatCoord_I could not be broadcast by read_satellite_input_files")

     call MPI_Bcast(T_tmp,Max_Satellite_Npts*7,MPI_Integer, &
  	  0,iComm,iError)
     if(iError>0)call stop_mpi(&
  	  "Satellite_Time could not be broadcast by read_satellite_input_files")
     do i=1,Satellite_Npts(isat)
        call time_int_to_real(T_tmp(i,:),Satellite_Time(isat,i))
     end do

     call MPI_Bcast(X_tmp,Max_Satellite_Npts*3,MPI_Real, &
  	  0,iComm,iError)
     if(iError>0)call stop_mpi(&
  	  "XSatellite_traj could not be broadcast by read_satellite_input_files")
     XSatellite_traj(isat,:,:) = X_tmp

  end do

end subroutine read_satellite_input_files

!==========================================================================

subroutine set_satellite_flags
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBlockMax,PROCtest,unusedBLK
  use ModGeometry, ONLY : XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModGeometry, ONLY : TypeGeometry               !^CFG IF NOT CARTESIAN
  use ModIO, ONLY : iBLKsatellite,iPEsatellite,Xsatellite,&
       SatelliteInBLK,DoTrackSatellite_I,nsatellite
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
       nsatellite,', call set_satellite_positions'

  call set_satellite_positions

  do iSat=1, nSatellite
     if(.not.DoTrackSatellite_I(iSat))cycle !Position is not defined
     select case(TypeGeometry)           !^CFG IF NOT CARTESIAN
     case('cartesian')                   !^CFG IF NOT CARTESIAN
        XSat=XSatellite(iSat,1)
        YSat=XSatellite(iSat,2)
        ZSat=XSatellite(iSat,3)
     case('spherical')                   !^CFG IF NOT CARTESIAN BEGIN
        call xyz_to_spherical(XSatellite(iSat,1),XSatellite(iSat,2),&
             XSatellite(iSat,3), XSat,YSat,ZSat)
     case('spherical_lnr')                   
        call xyz_to_spherical(XSatellite(iSat,1),XSatellite(iSat,2),&
             XSatellite(iSat,3), XSat,YSat,ZSat)
        XSat=log(max(XSat,cTiny))
     case default
        call stop_mpi('Unknown TypeGeometry='//TypeGeometry)
     end select                          !^CFG END CARTESIAN

     iPE = -1
     iBLKtemp = -1

     do iBLK = 1,nBlockMax

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
  use ModGeometry, ONLY : TypeGeometry               !^CFG IF NOT CARTESIAN
  use ModNumConst
  use ModIO
  implicit none


  integer :: i, isat, iBLK
  real    :: XSat,YSat,ZSat
  real*8  :: dtime, time_now

  !  real*8, external :: Time_Of_Year
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  if (iProc==0) &
       call set_oktest('set_satellite_positions',oktest, oktest_me)

  time_now = StartTime + Time_Simulation

  do isat = 1, nsatellite

     if (UseSatelliteFile(isat)) then

        if (Satellite_Npts(isat) > 0) then

           i = icurrent_satellite_position(isat)

           do while ((i < Satellite_Npts(isat)) .and.   &
                (Satellite_Time(isat,i) < time_now))
              i = i + 1
           enddo

           icurrent_satellite_position(isat) = i

           if ((i == Satellite_Npts(isat).and. &
                Satellite_Time(isat,i) <= time_now ).or.(i==1)) then 

              DoTrackSatellite_I(isat) = .false.

           else

              DoTrackSatellite_I(isat) = .true.

              dtime = 1.0 - (Satellite_Time(isat,i) - time_now) / &
                   (Satellite_Time(isat,i) - Satellite_Time(isat,i-1) + 1.0e-6)

              XSatellite(isat,:) = dtime * XSatellite_traj(isat,i,:) + &
                   (1.0 - dtime) * XSatellite_traj(isat,i-1,:) 

           endif

        endif

     else 

        call satellite_trajectory_formula(isat)

     end if

  end do

end subroutine set_satellite_positions


!=============================================================================

subroutine satellite_trajectory_formula(isat)
  use ModIO, ONLY : DoTrackSatellite_I,XSatellite,Satellite_name
  implicit none

  integer, intent(in) :: isat
  character (len=100) :: name_string
  real :: Xvect(3)

  name_string = trim(Satellite_name(isat))
  Xvect(:) = XSatellite(isat,:)

  ! Case should be for a specific satellite.  The trajectories can depend
  ! on the 'real' time so that the satellite knows where it is at.  For
  ! example, Cassini could be if'd on the date so that the code knows 
  ! whether the run is for a time near Earth, Jupiter or Saturn. 

  ! This routine should always set the TrackSatellite(isat) flag. When
  ! the satellite is at a useless position the time should return a
  ! do not track flag (DoTrackSatellite_I(isat) = .false.).


  select case(name_string)
  case ('earth')
     Xvect(1) = 5.0
     Xvect(2) = 5.0
     Xvect(3) = 5.0
     DoTrackSatellite_I(isat) = .true.
  case ('cassini')
     Xvect(1) = 5.0
     Xvect(2) = 5.0
     Xvect(3) = 5.0
     DoTrackSatellite_I(isat) = .true.
  case default 
     Xvect(1) = 1.0
     Xvect(2) = 1.0
     Xvect(3) = 1.0
     DoTrackSatellite_I(isat) = .false.
  end select

end subroutine satellite_trajectory_formula

!=============================================================================

subroutine open_satellite_output_files
  use ModIoUnit, ONLY : io_unit_new
  use ModIO, ONLY : nSatellite,Satellite_name,filename,&
       NamePlotDir,iUnitSat_I
  implicit none

  logical :: from_end = .true.
  integer :: isat, l1, l2
  character (len=4), Parameter :: IO_ext=".sat"
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  call set_oktest('open_satellite_output_files',oktest, oktest_me)

  do isat=1,nsatellite
     l1 = index(Satellite_name(isat),'/',from_end) + 1
     l2 = index(Satellite_name(isat),'.') - 1
     if (l1-1<=0) l1=1
     if (l2+1<=0) l2=len_trim(Satellite_name(isat))

     write(filename,'(a,i2.2,a)')trim(NamePlotDir)//&
          'satellite_',isat, &
          '_'//Satellite_name(isat)(l1:l2)//IO_ext
     if(oktest) write(*,*) 'isat,l1,l2: ',isat,l1,l2
     if(oktest) write(*,*) 'open_satellite_output_files: satellitename:', &
          Satellite_name(isat)
     if(oktest) write(*,*) 'open_satellite_output_files: filename:', &
          filename

     iUnitSat_I(iSat)=io_unit_new()
     open(iUnitSat_I(iSat),file=filename,status='unknown')
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

