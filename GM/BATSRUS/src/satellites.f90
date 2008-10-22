!^CFG COPYRIGHT UM
!==============================================================================
subroutine read_satellite_input_files
  use ModProcMH
  use ModMain, ONLY : nDim, lVerbose, TypeCoordSystem, StartTime
  use ModIO
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

subroutine set_satellite_flags(iSat)
  use ModProcMH
  use ModMain, ONLY : nDim,nI,nJ,nK,nBlockMax,PROCtest,unusedBLK
  use ModGeometry, ONLY : XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModGeometry, ONLY : UseCovariant               
  use ModIO, ONLY : iBLKsatellite,iPEsatellite,Xsatellite,&
       SatelliteInBLK,DoTrackSatellite_I,nSatellite
  use ModNumConst
  use ModMpi
  implicit none

  integer, intent(in) :: iSat
  integer :: iPE,iPEtmp, iBLK, iBLKtemp
  real    :: xSat,ySat,zSat
  integer :: i,j,k, iError
  real,dimension(nDim)::GenOut_D
  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  call set_oktest('set_satellite_flags',oktest, oktest_me)

  if (oktest_me) &
       write(*,*)'Starting set_satellite_flags',&
       nSatellite,', call set_satellite_positions'

  call set_satellite_positions(iSat)
     if(.not.DoTrackSatellite_I(iSat)) return !Position is not defined

     if(UseCovariant)then                  
        call xyz_to_gen(XSatellite(iSat,:),GenOut_D)
        xSat=GenOut_D(1)
        ySat=GenOut_D(2)
        zSat=GenOut_D(3)
     else                                   
        xSat=XSatellite(iSat,1)
        ySat=XSatellite(iSat,2)
        zSat=XSatellite(iSat,3)
     end if                                              

     iPE = -1
     iBLKtemp = -1

     do iBLK = 1, nBlockMax

        SatelliteInBLK(isat,iBLK) =.not.unusedBLK(iBLK).and.&
             xSat >  XyzStart_BLK(1,iBLK) - cHalf*dx_BLK(iBLK) .and. &
             xSat <= XyzStart_BLK(1,iBLK) + (nI-cHalf)*dx_BLK(iBLK) .and. &
             ySat >  XyzStart_BLK(2,iBLK) - cHalf*dy_BLK(iBLK) .and. &
             ySat <= XyzStart_BLK(2,iBLK) + (nJ-cHalf)*dy_BLK(iBLK) .and. &
             zSat >  XyzStart_BLK(3,iBLK) - cHalf*dz_BLK(iBLK) .and. &
             zSat <= XyzStart_BLK(3,iBLK) + (nK-cHalf)*dz_BLK(iBLK)

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
end subroutine set_satellite_flags


!=============================================================================

subroutine set_satellite_positions(iSat)
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,StartTime,time_simulation, &
       time_accurate
  use ModGeometry, ONLY : x1,x2,y1,y2,z1,z2,XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModGeometry, ONLY : TypeGeometry               
  use ModNumConst
  use ModIO
  implicit none

  integer, intent(in) :: iSat
  integer :: i, iBLK
  real    :: xSat,ySat,zSat
  real    :: dtime

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------
  if (iProc==0) &
       call set_oktest('set_satellite_positions',oktest, oktest_me)
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
!^CFG IF RAYTRACE BEGIN
subroutine sat_get_ray(iSatIn, SatRayVar_I)

  use ModRaytrace
  use ModMpi
  use ModMain,       ONLY: nI, nJ, nK, nCells, nBlock
  use ModGeometry,   ONLY: XyzStart_BLK, dx_BLK, dy_BLK, dz_BLK, UseCovariant
  use ModIO,         ONLY: iBLKsatellite, iPEsatellite, Xsatellite

  integer, intent(in) :: iSatIn
  real,    intent(out):: SatRayVar_I(5)
  
  character(len=*), parameter :: NameSub = 'get_sat_ray'
  integer  :: iDir, iBLK, iDim
  real     :: Xyz_D(3),  Dxyz_D(3), RayVars(3,2,nI,nJ,nK)
  real     :: Dx1, Dx2, Dy1, Dy2, Dz1, Dz2
  integer  :: i1, i2, j1, j2, k1, k2, iNear, jNear, kNear
  integer  :: i, j, k  ! Used in the FORALL statement

  !---------------------------------------------------------------------------

  ! Only use this if we're on the correct node.
  if (iProc /= iPEsatellite(iSatIn)) then
     do iDim=1,5
        SatRayVar_I(iDim) = 0.0
     enddo
     RETURN
  endif

  iBLK = iBLKSatellite(iSatIn)
  Dxyz_D(1) = dx_BLK(iBLK)
  Dxyz_D(2) = dy_BLK(iBLK)
  Dxyz_D(3) = dz_BLK(iBLK)

  if (UseCovariant) then 
     call xyz_to_gen(Xsatellite(iSatIn,:),Xyz_D)
  else
     Xyz_D = Xsatellite(iSatIn,:)
  endif

  ! Normalize coordinates
  do iDim=1,3
     Xyz_D(iDim) = (Xyz_D(iDim) - XyzStart_BLK(iDim,iBLK)) / Dxyz_D(iDim)
  end do

  ! Set location assuming point is inside block.
  i1 = floor(Xyz_D(1))
  j1 = floor(Xyz_D(2))  
  k1 = floor(Xyz_D(3))
  i2 = ceiling(Xyz_D(1))
  j2 = ceiling(Xyz_D(2))
  k2 = ceiling(Xyz_D(3))

  ! If Xyz_D is outside of block, change i,j,k in order to extrapolate.
  if(any( Xyz_D < 1) .or. any(Xyz_D > (/nI, nJ, nK/))) then
     i1 = min(nI-1, max(1, i1));   i2 = i1 + 1
     j1 = min(nJ-1, max(1, j1));   j2 = j1 + 1
     k1 = min(nK-1, max(1, k1));   k2 = k1 + 1
  endif

  ! Set interpolation weights
  Dx1 = Xyz_D(1) - i1; Dx2 = 1.0 - Dx1
  Dy1 = Xyz_D(2) - j1; Dy2 = 1.0 - Dy1
  Dz1 = Xyz_D(3) - k1; Dz2 = 1.0 - Dz1

  ! Calculate the nearest point.
  iNear = min( nI, max(nint(Xyz_D(1)),1) )
  jNear = min( nJ, max(nint(Xyz_D(2)),1) )
  kNear = min( nK, max(nint(Xyz_D(3)),1) )

  ! Copy ray tracing values to new array so allow changing of values.
  RayVars = ray(1:3,1:2,1:nI,1:nJ,1:nK,iBLK)

  ! Use the ray status of the nearest point to the satellite.
  SatRayVar_I(3) = RayVars(3, 1, iNear,jNear,kNear)

  ! For each direction along the ray, determine if all lines surrounding point
  ! are open or closed, then set SatRayVar_I accordingly.
  do iDir=1,2

     if ( any(RayVars(3,1,i1:i2,j1:j2,k1:k2) < 1) .or. &
          any(RayVars(3,1,i1:i2,j1:j2,k1:k2) == iDir) ) then
        ! One or more lines is open in direction iDir, must use nearest point.
        do iDim=1,2
           SatRayVar_I(iDim + 3*(iDir-1))=RayVars(iDim,iDir,iNear,jNear,kNear)
        end do

     else   ! All lines closed in direction iDir, interpolate.

        ! If the satellite is near the 0/360 degree boundary in longitude,
        ! the result of the interpolation will be incorrect.  Adjust
        ! longitudes accordingly.
        if (any(RayVars(2,iDir,i1:i2,j1:j2,k1:k2)>330.0) .AND. &
            any(RayVars(2,iDir,i1:i2,j1:j2,k1:k2)<30.0)) then

           do i=i1,i2
              do j=j1,j2
                 do k=k1,k2
                    if (RayVars(2,iDir,i,j,k) < 30.0) &
                       RayVars(2,iDir,i,j,k) = RayVars(2,iDir,i,j,k) + 360.0
                 enddo
              enddo
           enddo
           !forall(i=i1:i2,j=j1:j2,k=k1:k2,RayVars(2,iDir,i,j,k)<30.0)
           !   RayVars(2,iDir,i,j,k) = RayVars(2,iDir,i,j,k) + 360.0
           !end forall
        endif

        do iDim=1,2
           SatRayVar_I(iDim + 3*(iDir-1)) =                         &
                 Dz2*(   Dy2*(   Dx2*RayVars(iDim,iDir,i1,j1,k1)   &
                +                Dx1*RayVars(iDim,iDir,i2,j1,k1))  &
                +        Dy1*(   Dx2*RayVars(iDim,iDir,i1,j2,k1)   &
                +                Dx1*RayVars(iDim,iDir,i2,j2,k1))) &
                +Dz1*(   Dy2*(   Dx2*RayVars(iDim,iDir,i1,j1,k2)   &
                +                Dx1*RayVars(iDim,iDir,i2,j1,k2))  &
                +        Dy1*(   Dx2*RayVars(iDim,iDir,i1,j2,k2)   &
                +                Dx1*RayVars(iDim,iDir,i2,j2,k2)))
        end do
     endif

  end do

  ! Ensure that longitude wraps around 0/360 degree boundary correctly.
  if(SatRayVar_I(2)<000.0) SatRayVar_I(2) = SatRayVar_I(2) + 360.0
  if(SatRayVar_I(5)<000.0) SatRayVar_I(5) = SatRayVar_I(5) + 360.0
  if(SatRayVar_I(2)>360.0) SatRayVar_I(2) = SatRayVar_I(2) - 360.0
  if(SatRayVar_I(5)>360.0) SatRayVar_I(5) = SatRayVar_I(5) - 360.0

end subroutine sat_get_ray
!=============================================================================
!^CFG END RAYTRACE
