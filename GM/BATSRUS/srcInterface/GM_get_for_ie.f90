!^CFG COPYRIGHT UM
!^CMP FILE IE
!============================================
!                                           |
!     Magnetosphere/Ionosphere Coupling     |
!     Subroutines here are used with MPI    |
!============================================
module ModFieldAlignedCurrent

  use ModIeGrid
  
  implicit none
  save
  real, allocatable :: &
       FieldAlignedCurrent_II(:,:), bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)

contains

  !============================================================================

  subroutine init_mod_field_aligned_current(iSize,jSize)
    use ModProcMH,         ONLY: iProc, iComm
    integer, intent(in) :: iSize, jSize

    if(allocated(FieldAlignedCurrent_II)) RETURN

    call init_mod_ie_grid(iSize, jSize)

    allocate( bCurrent_VII(0:6, nThetaIono, nPhiIono), &
         bCurrentLocal_VII(0:6, nThetaIono, nPhiIono), &
         FieldAlignedCurrent_II(nThetaIono, nPhiIono))

  end subroutine init_mod_field_aligned_current

  !============================================================================
  subroutine clean_mod_field_aligned_current

    if(allocated(FieldAlignedCurrent_II)) &
         deallocate( bCurrent_VII, bCurrentLocal_VII, FieldAlignedCurrent_II)
  end subroutine clean_mod_field_aligned_current
 
end module ModFieldAlignedCurrent

!==============================================================================

subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

  use CON_axes, ONLY: transform_matrix
  use ModMain, ONLY: Time_Simulation, TypeCoordSystem
  use ModIeGrid
  use ModCurrent,            ONLY: calc_field_aligned_current
  use ModFieldAlignedCurrent,ONLY: FieldAlignedCurrent_II, &
       init_mod_field_aligned_current!, calc_field_aligned_current
  use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
       InvB_, RhoInvB_, pInvB_, xEnd_, yEnd_, zEnd_, CLOSEDRAY, GmSm_DD
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics, ONLY: No2Si_V, UnitP_, UnitRho_, UnitB_,rCurrents, UnitJ_
  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  integer, intent(in) :: iSize, jSize, nVar
  real, intent(out), dimension(iSize, jSize, nVar) :: Buffer_IIV
  integer :: i, j
  real :: Radius,iError, Phi, Theta, LatBoundary
  real, allocatable, dimension(:), save :: IE_lat, IE_lon
  real, allocatable, dimension(:,:,:)   :: bSm_DII 
  real :: XyzIono_D(3), RtpIono_D(3), Xyz_D(3), Lat,Lon, dLat,dLon
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  if(.not.allocated(FieldAlignedCurrent_II)) &
       call init_mod_field_aligned_current(iSize, jSize)

  if(.not.allocated(IE_lat)) &
       allocate(IE_lat(iSize), IE_lon(jSize))

  if(.not. allocated(bSm_DII))&
       allocate(bSm_DII(3,iSize,jSize))

  ! initialize all elements to zero on proc 0, others should not use it
  if(iProc == 0) Buffer_IIV = 0.0

  ! Put field aligned currents into the first "slice" of the buffer
  call calc_field_aligned_current(nThetaIono, nPhiIono, rIonosphere, &
       FieldAlignedCurrent_II, bSm_DII, LatBoundary, ThetaIono_I, PhiIono_I)

  ! Take radial component of FAC and put it into the buffer sent to IE
  ! The resulting FAC_r will be positive for radially outgoing current
  ! and negative for radially inward going currents.
  if(iProc==0)then
     do j=1, nPhiIono
        Phi = PhiIono_I(j)
        
        do i=1, nThetaIono
           Theta = ThetaIono_I(i)
           call sph_to_xyz(rIonosphere, Theta, Phi, XyzIono_D)
           FieldAlignedCurrent_II(i,j) = FieldAlignedCurrent_II(i,j) &
                * sum(bSm_DII(:,i,j)*XyzIono_D) &
                / (sqrt(sum((bSm_DII(:,i,j))**2)) * rIonosphere)
        end do
     end do

     ! Save the latitude boundary information to the equator 
     Buffer_IIV(:,:,1) = FieldAlignedCurrent_II(:,:)*No2Si_V(UnitJ_)
     Buffer_IIV(nThetaIono/2:nThetaIono/2+1,1,1) = LatBoundary                          

  end if
  
  deallocate(bSm_DII)

  if(DoTraceIE) then
     ! Load grid and convert to lat-lon in degrees
     IE_lat = 90.0 - cRadToDeg * ThetaIono_I(1:iSize)
     IE_lon =        cRadToDeg * PhiIono_I
     Radius = (6378.+100.)/6378.
     call integrate_ray_accurate(iSize, jSize, IE_lat, IE_lon, Radius, &
          'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')

     ! Only processor 0 has the resulting integrals,
     !   the others do not participate in the coupling
     if(iProc == 0)then
        where(RayResult_VII(InvB_,:,:)>0.)
           RayResult_VII(RhoInvB_,:,:)=RayResult_VII(RhoInvB_,:,:)/RayResult_VII(InvB_,:,:)
           RayResult_VII(pInvB_,:,:)=RayResult_VII(pInvB_,:,:)/RayResult_VII(InvB_,:,:)
        end where
        where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
           RayResult_VII(   InvB_,:,:) = -1.
           RayResult_VII(rhoInvB_,:,:) = 0.
           RayResult_VII(  pInvB_,:,:) = 0.
        end where
        Buffer_IIV(:,:,2) = RayResult_VII(   InvB_,:,:) / No2Si_V(UnitB_)
        Buffer_IIV(:,:,3) = RayResult_VII(rhoInvB_,:,:) * No2Si_V(UnitRho_)
        Buffer_IIV(:,:,4) = RayResult_VII(  pInvB_,:,:) * No2Si_V(UnitP_)

        ! Transformation matrix from default (GM) to SM coordinates
        GmSm_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')

        ! Loop to compute deltas
        do i=1,iSize
           do j=1,jSize
              Lat = -IE_Lat(i)
              Lon =  IE_Lon(j)
              if(RayResult_VII(InvB_,i,j)>1.e-10)then
                 XyzIono_D(1)=RayResult_VII(xEnd_,i,j)
                 XyzIono_D(2)=RayResult_VII(yEnd_,i,j)
                 XyzIono_D(3)=RayResult_VII(zEnd_,i,j)
                 XyzIono_D = matmul(GmSm_DD,XyzIono_D)
                 call xyz_to_sph(XyzIono_D, RtpIono_D)
                 Lat=90.-RtpIono_D(2)*cRadToDeg
                 Lon=    RtpIono_D(3)*cRadToDeg
              end if

              ! Shift to closer to local value, even if outside 0-360
              if( (IE_Lon(j)-Lon)<-180. ) Lon = Lon-360.
              if( (IE_Lon(j)-Lon)> 180. ) Lon = Lon+360.

              ! Compute deltas
              dLat = abs(Lat) - abs(IE_Lat(i))
              dLon =     Lon  -     IE_Lon(j)

              ! Put into exchange buffer
              Buffer_IIV(i,j,5) = dLat
              Buffer_IIV(i,j,6) = dLon
           end do
        end do
     end if
     deallocate(RayIntegral_VII, RayResult_VII)
  end if

  if(DoTest)write(*,*)NameSub,': finished'
end subroutine GM_get_for_ie


!**************************************************************
subroutine ieb_plot(iFile)
  use CON_line_extract,  ONLY: line_get, line_clean
  use CON_planet_field,  ONLY: map_planet_field
  use CON_axes,          ONLY: transform_matrix
  use ModIoUnit,         ONLY: UnitTmp_
  use ModAdvance,        ONLY: nVar, Ux_, Uz_, Bx_, Bz_
  use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, time_accurate, n_step
  use ModNumConst,       ONLY: cRadToDeg, cDegToRad
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitX_, rBody
  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
  use ModIO,             ONLY: StringDateOrTime, NamePlotDir
  use ModRaytrace,       ONLY: RayResult_VII, RayIntegral_VII
  use ModIeGrid
  implicit none

  integer, intent(in) :: iFile

  character (len=*), parameter :: NameSub='ieb_plot'
  character (len=80) :: FileName,stmp
  character (len=2) :: Coord
  character (len=1) :: NS
  integer :: i,j,k, nLat,nLon, nLine, nTP, iStart,iEnd, iLat,iLon, OC
  integer :: iPoint, nPoint, nVarOut, iHemisphere, nFile
  real :: PlotVar_V(0:4+nVar)
  real :: Radius, Lat,Lon, Theta,Phi, LonOC
  real :: XyzIono_D(3), RtpIono_D(3), Xyz_D(3)
  real :: Gsm2Smg_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real :: Smg2Gsm_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real, allocatable :: PlotVar_VI(:,:), IE_lat(:), IE_lon(:)
  logical :: MapDown
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  nLat=nThetaIono
  nLon=int((nPhiIono-1)/5)
  if(.not.allocated(IE_lat)) allocate(IE_lat(nLat), IE_lon(nLon))

  ! Load grid and convert to lat-lon in degrees
  IE_lat = 90.0 - cRadToDeg * ThetaIono_I(1:nThetaIono)
  do i=1,nLon
     j=1+(i-1)*5
     IE_lon(i) =        cRadToDeg * PhiIono_I(j)
  end do
  Radius = (6378.+100.)/6378.
  nTP=int( (rBody-Radius)/.1 )

  call integrate_ray_accurate(nLat, nLon, IE_lat, IE_lon, Radius, 'extract_I')

  if(iProc == 0)then

     ! Transformation matrix from default (GM) to SM coordinates
     Gsm2Smg_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')
     Smg2Gsm_DD = transform_matrix(time_simulation,'SMG','GSM')

     call line_get(nVarOut, nPoint)
     if(nPoint>0)then
        !PlotVar_VI variables = 'iLine l x y z rho ux uy uz bx by bz p'
        allocate(PlotVar_VI(0:nVarOut, nPoint))
        call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)

        do nFile=1,4

           if(nFile==1)then
              Coord = 'SM';  NS = 'N'
           elseif(nFile==2)then
              Coord = 'GM';  NS = 'N'
           elseif(nFile==3)then
              Coord = 'SM';  NS = 'S'
           elseif(nFile==4)then
              Coord = 'GM';  NS = 'S'
           end if
           FileName=trim(NamePlotDir)//'IEB-'//trim(Coord)//'-'//trim(NS)
           if(time_accurate)then
              call get_time_string
              FileName = trim(FileName) // "_t" // StringDateOrTime
           end if
           write(FileName,'(a,i7.7,a)') trim(FileName)//"_n", n_step,".dat"

           open( UnitTmp_, FILE=trim(FileName), STATUS="replace")
           if(Coord == 'GM')then
              write(UnitTmp_,'(a)')'TITLE="IE B traces (GM Coordinates)"'
           else
              write(UnitTmp_,'(a)')'TITLE="IE B traces (SM Coordinates)"'
           end if
           write(UnitTmp_,'(a)')'VARIABLES="X [R]", "Y [R]", "Z [R]", "Lat", "Lon", "OC"'

           k=0
           LonOC=-1.
           do iPoint = 1, nPoint
              nLine=PlotVar_VI(0,iPoint)
              if(k /= nLine)then
                 !\\
                 ! finish previous line
                 if(k/=0)then
                    iEnd = iPoint-1
                    MapDown = .false.
                    Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
                    if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)MapDown = .true.
                    j=(1+iEnd-iStart)+(nTP+1)
                    if(MapDown)j=j+(nTP+1)
                    OC=-1; if(MapDown)OC=2
                    if(MapDown .and. LonOC/=Lon)OC=1

!\
!                     write(UnitTmp_,'(a,2f7.2,a,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
!                          ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!-
                    write(UnitTmp_,'(a,2f7.2,a,a,f7.2,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
                         ', STRANDID=1, SOLUTIONTIME=',Lon, &
                         ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!/
                    write(stmp,'(f8.2)')Lat
                    write(UnitTmp_,'(a,a,a)') 'AUXDATA LAT="',trim(adjustl(stmp)),'"'
                    write(stmp,'(f8.2)')Lon
                    write(UnitTmp_,'(a,a,a)') 'AUXDATA LON="',trim(adjustl(stmp)),'"'

                    ! Convert to SMG Cartesian coordinates on the surface of the ionosphere
                    Theta = cDegToRad*(90.0 - Lat)     
                    Phi = cDegToRad*Lon
                    call sph_to_xyz(Radius, Theta, Phi, XyzIono_D)
                    Xyz_D=XyzIono_D
                    if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                    write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    do i=1,nTP
                       ! Map from the ionosphere to rBody
                       call map_planet_field(time_simulation, XyzIono_D, 'SMG NORM', &
                            Radius+i*.1, Xyz_D, iHemisphere)
                       if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                       write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    end do
                    do i=iStart,iEnd
                       ! Convert vectors to SM coordinates
                       PlotVar_V = PlotVar_VI(:, i)
                       PlotVar_V(2:4) = matmul(Gsm2Smg_DD,PlotVar_V(2:4))
                       PlotVar_V(2:4) = PlotVar_V(2:4) * Si2No_V(UnitX_)
                       Xyz_D = PlotVar_V(2:4)
                       if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                       write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    end do
                    if(MapDown)then
                       Xyz_D=PlotVar_V(2:4)
                       do i=nTP,0,-1
                          ! Map from rBody to the ionosphere
                          call map_planet_field(time_simulation, Xyz_D, 'SMG NORM', &
                               Radius+i*.1, XyzIono_D, iHemisphere)
                          if(Coord == 'GM') XyzIono_D = matmul(Smg2Gsm_DD,XyzIono_D)
                          write(UnitTmp_, *) XyzIono_D,Lat,Lon,OC
                       end do
                    end if
                 end if

                 !\\
                 ! start new line counters
                 k=nLine
                 iStart = iPoint
                 iLon=1+((nLine-1)/nLat)
                 iLat=nLine-(iLon-1)*nLat
                 Lon=IE_lon(iLon)
                 Lat=IE_lat(iLat)
                 if(NS == 'N')then
                    if(Lat<0.)k=0
                 else
                    if(Lat>0.)k=0
                 end if
              end if
           end do

           !\\
           ! finish last line
           if(k/=0)then
              iEnd = nPoint
              MapDown = .false.
              Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
              if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)MapDown = .true.
              j=(1+iEnd-iStart)+(nTP+1)
              if(MapDown)j=j+(nTP+1)
              OC=-1; if(MapDown)OC=2
              if(MapDown .and. LonOC/=Lon)OC=1

!\
!              write(UnitTmp_,'(a,2f7.2,a,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
!                   ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!-
              write(UnitTmp_,'(a,2f7.2,a,a,f7.2,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
                   ', STRANDID=1, SOLUTIONTIME=',Lon, &
                   ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!/
              write(stmp,'(f8.2)')Lat
              write(UnitTmp_,'(a,a,a)') 'AUXDATA LAT="',trim(adjustl(stmp)),'"'
              write(stmp,'(f8.2)')Lon
              write(UnitTmp_,'(a,a,a)') 'AUXDATA LON="',trim(adjustl(stmp)),'"'

              ! Convert to SMG Cartesian coordinates on the surface of the ionosphere
              Theta = cDegToRad*(90.0 - Lat)     
              Phi = cDegToRad*Lon
              call sph_to_xyz(Radius, Theta, Phi, XyzIono_D)
              Xyz_D=XyzIono_D
              if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
              write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              do i=1,nTP
                 ! Map from the ionosphere to rBody
                 call map_planet_field(time_simulation, XyzIono_D, 'SMG NORM', &
                      Radius+i*.1, Xyz_D, iHemisphere)
                 if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                 write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              end do
              do i=iStart,iEnd
                 ! Convert vectors to SM coordinates
                 PlotVar_V = PlotVar_VI(:, i)
                 PlotVar_V(2:4) = matmul(Gsm2Smg_DD,PlotVar_V(2:4))
                 PlotVar_V(2:4) = PlotVar_V(2:4) * Si2No_V(UnitX_)
                 Xyz_D = PlotVar_V(2:4)
                 if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                 write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              end do
              if(MapDown)then
                 Xyz_D=PlotVar_V(2:4)
                 do i=nTP,0,-1
                    ! Map from the ionosphere to rBody
                    call map_planet_field(time_simulation, Xyz_D, 'SMG NORM', &
                         Radius+i*.1, XyzIono_D, iHemisphere)
                    if(Coord == 'GM') XyzIono_D = matmul(Smg2Gsm_DD,XyzIono_D)
                    write(UnitTmp_, *) XyzIono_D,Lat,Lon,OC
                 end do
              end if
           end if

           close(UnitTmp_)
        end do

        deallocate(PlotVar_VI)
     end if
     call line_clean
  end if

  if(allocated(RayIntegral_VII)) deallocate(RayIntegral_VII)
  if(allocated(RayResult_VII))   deallocate(RayResult_VII)

  if(DoTest)write(*,*)NameSub,': finished'
end subroutine ieb_plot
