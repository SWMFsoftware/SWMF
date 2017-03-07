!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE

!============================================================================

module GM_couple_ie

  use ModProcMH
  use ModIeCoupling, ONLY: rIonosphere, &
       nThetaIono, nPhiIono, ThetaIono_I, PhiIono_I,&
       IonoPotential_II, IonoJouleHeating_II, &
       SigmaHall_II, SigmaPedersen_II, jHall_DII, jPedersen_DII, &
       init_ie_grid, calc_grad_ie_potential, map_jouleheating_to_inner_bc

  implicit none
  save

  private ! except

  public:: GM_get_for_ie
  public:: GM_get_info_for_ie
  public:: GM_put_from_ie

contains

  !============================================================================
  subroutine GM_get_info_for_ie(nVar, NameVar_I)

    use ModMain,             ONLY: TypeFaceBc_I, body1_
    use ModGroundMagPerturb, ONLY: nMagTotal

    ! Pass number and names of variables requested from IE

    integer, intent(out) :: nVar
    character(len=*), intent(out), optional:: NameVar_I(:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='GM_get_info_for_ie'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Always pass electric potential
    nVar = 1
    if(present(NameVar_I)) NameVar_I(nVar) = 'potential'

    ! Joule heating is needed for empirical inner boundary
    if(TypeFaceBc_I(body1_) == 'ionosphereoutflow')then
       nVar = nVar+1
       if(present(NameVar_I)) NameVar_I(nVar) = 'jouleheat'
    end if

    ! Conductances are needed for magnetometer calculations
    if(nMagTotal > 0)then
       nVar = nVar + 1
       if(present(NameVar_I)) NameVar_I(nVar) = 'sigmahall'
       nVar = nVar + 1
       if(present(NameVar_I)) NameVar_I(nVar) = 'sigmapedersen'
    end if

    if(DoTestMe)then
       write(*,*) NameSub//' set nVar=', nVar
       if(present(NameVar_I)) write(*,*) NameSub//' set NameVar_I=', &
            NameVar_I
    end if

  end subroutine GM_get_info_for_ie

  !============================================================================
  subroutine GM_get_for_ie(Buffer_IIV, iSize, jSize, nVar)

    ! Send the following information from GM to IE on the IE grid:
    !  1. radial component of the field-aligned-currents (FACs)
    !  2. latitude boundary (LatBoundary) scalar
    !  3. field line tracing information if DoTraceIE is true

    use CON_axes, ONLY: transform_matrix
    use ModMain, ONLY: Time_Simulation, TypeCoordSystem
    use ModCurrent,            ONLY: calc_field_aligned_current
    use ModRaytrace, ONLY: DoTraceIE, RayResult_VII, RayIntegral_VII, &
         InvB_, RhoInvB_, pInvB_, xEnd_, yEnd_, zEnd_, CLOSEDRAY, GmSm_DD
    use ModNumConst, ONLY: cRadToDeg
    use ModProcMH,         ONLY: iProc
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitP_, UnitRho_, UnitB_, UnitJ_
    use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
    use CON_coupler, ONLY: Grid_C, IE_

    integer, intent(in) :: iSize, jSize, nVar
    real,    intent(out):: Buffer_IIV(iSize,jSize,nVar)

    integer :: i, j
    real :: Radius, Phi, Theta, LatBoundary
    real, allocatable:: FieldAlignedCurrent_II(:,:), bSm_DII(:,:,:)
    real, allocatable:: IE_lat(:), IE_lon(:)
    real :: XyzIono_D(3), RtpIono_D(3), Lat,Lon, dLat,dLon
    logical :: DoTest, DoTestMe

    character (len=*), parameter :: NameSub='GM_get_for_ie'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,': starting'

    if(nThetaIono < 1) call init_ie_grid( &
         Grid_C(IE_) % Coord1_I, &
         Grid_C(IE_) % Coord2_I, &
         Grid_C(IE_) % Coord3_I(1))

    allocate(FieldAlignedCurrent_II(iSize,jSize), bSm_DII(3,iSize,jSize))

    ! Put field aligned currents into the first variable of the buffer
    call calc_field_aligned_current(nThetaIono, nPhiIono, rIonosphere, &
         FieldAlignedCurrent_II, bSm_DII, LatBoundary, ThetaIono_I, PhiIono_I)

    ! Take radial component of FAC and put it into the buffer sent to IE
    ! The resulting FAC_r will be positive for radially outgoing current
    ! and negative for radially inward going currents.
    if(iProc==0)then
       ! initialize all elements to zero on proc 0, others should not use it
       Buffer_IIV = 0.0

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

    deallocate(FieldAlignedCurrent_II, bSm_DII)

    if(DoTraceIE) then
       allocate(IE_lat(iSize), IE_lon(jSize))

       ! Load grid and convert to lat-lon in degrees
       IE_lat = 90.0 - cRadToDeg * ThetaIono_I(1:iSize)
       IE_lon =        cRadToDeg * PhiIono_I
       Radius = (6378.+100.)/6378.
       call integrate_ray_accurate(iSize, jSize, IE_lat, IE_lon, Radius, &
            'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')

       ! Only processor 0 has the resulting integrals,
       !   the others do not participate in the coupling
       if(iProc == 0)then
          where(RayResult_VII(InvB_,:,:) > 0.)
             RayResult_VII(RhoInvB_,:,:) = RayResult_VII(RhoInvB_,:,:) &
                  /RayResult_VII(InvB_,:,:)
             RayResult_VII(pInvB_,:,:)   = RayResult_VII(pInvB_,:,:) &
                  /RayResult_VII(InvB_,:,:)
          end where
          where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
             RayResult_VII(   InvB_,:,:) = -1.
             RayResult_VII(rhoInvB_,:,:) = 0.
             RayResult_VII(  pInvB_,:,:) = 0.
          end where
          Buffer_IIV(:,:,2) = RayResult_VII(   InvB_,:,:) &
               * No2Si_V(UnitX_)/No2Si_V(UnitB_)
          Buffer_IIV(:,:,3) = RayResult_VII(RhoInvB_,:,:) * No2Si_V(UnitRho_)
          Buffer_IIV(:,:,4) = RayResult_VII(  pInvB_,:,:) * No2Si_V(UnitP_)

          ! Transformation matrix from default (GM) to SM coordinates
          GmSm_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')

          ! Loop to compute deltas
          do j=1,jSize
             do i=1,iSize
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
       deallocate(IE_lat, IE_lon, RayIntegral_VII, RayResult_VII)
    end if

    if(DoTest)write(*,*)NameSub,': finished'

  end subroutine GM_get_for_ie

  !============================================================================
  subroutine GM_put_from_ie(Buffer_IIV, iSize, jSize, nVar, NameVar_I)

    ! Receive nVar variables listed in NameVar_I from IE on the IE grid
    ! and store them into the allocatable arrays in ModIeCoupling

    use ModPhysics, ONLY: &
         Si2No_V, UnitX_, UnitElectric_, UnitPoynting_, UnitJ_
    use CON_coupler, ONLY: Grid_C, IE_

    integer, intent(in) :: iSize, jSize, nVar
    real,    intent(in) :: Buffer_IIV(iSize,jSize,nVar)
    character(len=*), intent(in):: NameVar_I(nVar)

    integer:: iVar

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='GM_put_from_ie'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)NameSub,': iSize,jSiz,nVare=', iSize, jSize, nVar

    if(nThetaIono < 1) call init_ie_grid( &
         Grid_C(IE_) % Coord1_I, &
         Grid_C(IE_) % Coord2_I, &
         Grid_C(IE_) % Coord3_I(1))

    do iVar = 1, nVAr
       select case(NameVar_I(iVar))
       case('potential')
          if(.not. allocated(IonoPotential_II)) &
               allocate(IonoPotential_II(nThetaIono, nPhiIono))
          ! Electric potential has units of [E]*[x] (in SI: V/m*m=V)
          IonoPotential_II = Buffer_IIV(:,:,iVar) &
               *Si2No_V(UnitElectric_)*Si2No_V(UnitX_)
          call calc_grad_ie_potential
       case('jouleheat')
          if (.not. allocated(IonoJouleHeating_II)) &
               allocate(IonoJouleHeating_II(nThetaIono, nPhiIono))
          ! Joule heating variable has units of Poynting flux
          IonoJouleHeating_II = Buffer_IIV(:,:,2) &
               *Si2No_V(UnitPoynting_)

          ! scale Joule heating from ionosphere to GM inner BC
          call map_jouleheating_to_inner_bc
       case('sigmahall')
          if(.not. allocated(SigmaHall_II)) &
               allocate(SigmaHall_II(nThetaIono,nPhiIono))

          ! Height integrated conductivty has units of [j]*[x]/[E] (in SI:A/V)
          SigmaHall_II = Buffer_IIV(:,:,iVar) &
               *Si2No_V(UnitJ_)*Si2No_V(UnitX_)/Si2No_V(UnitElectric_)

          ! The ionosphere currents will need recalculation, so deallocate them
          if(allocated(jHall_DII)) deallocate(jHall_DII)

       case('sigmapedersen')
          if(.not. allocated(SigmaPedersen_II)) &
               allocate(SigmaPedersen_II(nThetaIono,nPhiIono))

          ! Height integrated conductivty has units of [j]*[x]/[E] (in SI:A/V)
          SigmaPedersen_II = Buffer_IIV(:,:,iVar) &
               *Si2No_V(UnitJ_)*Si2No_V(UnitX_)/Si2No_V(UnitElectric_)

          ! The ionosphere currents will need recalculation, so deallocate them
          if(allocated(jPedersen_DII)) deallocate(jPedersen_DII)
       case default
          call stop_mpi(NameSub//' unknown variable='//NameVar_I(iVar))
       end select
    end do

    if(DoTest)write(*,*)NameSub,': done'

  end subroutine GM_put_from_ie

end module GM_couple_ie
