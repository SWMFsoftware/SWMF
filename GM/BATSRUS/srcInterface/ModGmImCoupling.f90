!^CFG COPYRIGHT UM
!^CMP FILE IM
!^CFG FILE RAYTRACE
module ModGmImCoupling

  use ModMpi
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use CON_coupler, ONLY: Grid_C, IM_, ncells_decomposition_d

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace, ONLY : ray,rayface
  use ModPhysics, ONLY: No2Si_V, Si2No_V, &
       UnitP_, UnitRho_, UnitTemperature_, UnitB_, &
       Bdp, DipoleStrengthSi, rCurrents, rBody
  implicit none

  character (len=*), parameter :: NameMod='ModGmImCoupling'

  ! IM Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the IM grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RCM_lat, RCM_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_PE_vol,  MHD_SUM_vol, MHD_tmp, &
       MHD_PE_rho,  MHD_SUM_rho, &
       MHD_PE_p,    MHD_SUM_p, &
       MHD_PE_Beq,  MHD_Beq, &
       MHD_PE_Xeq,  MHD_Xeq, &
       MHD_PE_Yeq,  MHD_Yeq, &
       MHD_Fluxerror
  real, parameter :: noValue=-99999.
  real :: eqB,xL,yL,zL
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s6,s8,factor1,factor2

  integer :: iError

  logical :: dbg0=.false.

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

  logical :: DoTestTec, DoTestIdl

contains

  subroutine allocate_gm_im(iSizeIn,jSizeIn)
    use CON_comp_param, ONLY: IM_

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_im'

    if(allocated(MHD_lat_boundary)) deallocate(MHD_lat_boundary)
    if(allocated(MHD_PE_vol))       deallocate(MHD_PE_vol)
    if(allocated(MHD_SUM_vol))      deallocate(MHD_SUM_vol)
    if(allocated(MHD_PE_rho))       deallocate(MHD_PE_rho)
    if(allocated(MHD_SUM_rho))      deallocate(MHD_SUM_rho)
    if(allocated(MHD_PE_p))         deallocate(MHD_PE_p)
    if(allocated(MHD_SUM_p))        deallocate(MHD_SUM_p)
    if(allocated(MHD_tmp))          deallocate(MHD_tmp)
    if(allocated(MHD_PE_Beq))       deallocate(MHD_PE_Beq)
    if(allocated(MHD_PE_Xeq))       deallocate(MHD_PE_Xeq)
    if(allocated(MHD_PE_Yeq))       deallocate(MHD_PE_Yeq)
    if(allocated(MHD_Beq))          deallocate(MHD_Beq)
    if(allocated(MHD_Xeq))          deallocate(MHD_Xeq)
    if(allocated(MHD_Yeq))          deallocate(MHD_Yeq)
    if(allocated(MHD_Fluxerror))    deallocate(MHD_Fluxerror)

    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(RCM_lat))then
       nCells_D=ncells_decomposition_d(IM_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RCM_lat(iSize), RCM_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       RCM_lat = 90.0 - Grid_C(IM_) % Coord1_I * cRadToDeg
       RCM_lon =        Grid_C(IM_) % Coord2_I * cRadToDeg
    end if

    ! Arrays needed for the field line integrals
    allocate( MHD_SUM_vol(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_vol")
    MHD_SUM_vol = 0.

    allocate( MHD_SUM_rho(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_rho")
    MHD_SUM_rho = 0.

    allocate( MHD_SUM_p(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_SUM_p")
    MHD_SUM_p = 0.

    allocate( MHD_Beq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Beq")
    MHD_Beq = 0.

    allocate( MHD_Xeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Xeq")
    MHD_Xeq = 0.

    allocate( MHD_Yeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Yeq")
    MHD_Yeq = 0.

    allocate( MHD_PE_vol(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_vol")
    MHD_PE_vol = 0.

    allocate( MHD_PE_rho(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_rho")
    MHD_PE_rho = 0.

    allocate( MHD_PE_p(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_p")
    MHD_PE_p = 0.

    allocate( MHD_PE_Beq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Beq")
    MHD_PE_Beq = noValue

    allocate( MHD_PE_Xeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Xeq")
    MHD_PE_Xeq = noValue

    allocate( MHD_PE_Yeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_PE_Yeq")
    MHD_PE_Yeq = noValue

    allocate( MHD_tmp(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_tmp")
    MHD_tmp = 0.

    allocate( MHD_Fluxerror(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Fluxerror")
    MHD_Fluxerror = 0.

    allocate( MHD_lat_boundary(jsize), stat=iError )
    call alloc_check(iError,"MHD_lat_boundary")
    MHD_lat_boundary = 0

  end subroutine allocate_gm_im

  !============================================================================
  subroutine write_integrated_data_tec
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: j2, nCall=0
    real :: tmpT, tmpV1,tmpV2, lonShift
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
         ', "Xeq", "Yeq"', &
         ', "Volume", "Volume**(-2/3)"', &
         ', "MHD `r", "MHD p", "MHD T", "Beq"', &
         ', "FluxError"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)*Si2No_V(UnitP_))/(MHD_SUM_rho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          write(UNITTMP_,'(2i4,12G14.6)') j2,i,RCM_lon(j)+lonShift,RCM_lat(i),&
               MHD_lat_boundary(j), &
               MHD_Xeq(i,j),MHD_Yeq(i,j), &
               tmpV1,tmpV2, &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j), &
               MHD_Fluxerror(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_tec

  !============================================================================
  subroutine write_integrated_data_idl

    use ModIoUnit, ONLY: UNITTMP_
    use ModMain,   ONLY: time_simulation
    CHARACTER (LEN=100) :: filename
    integer :: nCall = 0
    !-------------------------------------------------------------------------

    !write values to plot file
    nCall = nCall+1
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCall,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,7
    write(UNITTMP_,'(3i4)')            jSize+1,iSize
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)') 'Lon Lat Xeq Yeq vol rho p Beq FluxError nothing'
    do i=isize,1,-1
       do j=1,jsize
          write(UNITTMP_,'(100(1pe18.10))') &
               RCM_lon(j),       &
               RCM_lat(i),       &
               MHD_Xeq(i,j),     &
               MHD_Yeq(i,j),     &
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j), &
               MHD_SUM_p(i,j),   &
               MHD_Beq(i,j),     &
               MHD_FluxError(i,j)
       end do
       write(UNITTMP_,'(100(1pe18.10))') &
               RCM_lon(1)+360.0, &
               RCM_lat(i),       &
               MHD_Xeq(i,1),     &
               MHD_Yeq(i,1),     &
               MHD_SUM_vol(i,1), &
               MHD_SUM_rho(i,1), &
               MHD_SUM_p(i,1),   &
               MHD_Beq(i,1),     &
               MHD_FluxError(i,1)               
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_idl

  !===========================================================================
  subroutine process_integrated_data

    integer :: iLoc_I(1), iEquator, iNorthPole
    !-----------------------------------------------------------------------
    if(dbg0)then
       if(iProc==0)write(*,*)' -D-'
    end if

    where(MHD_SUM_vol>0.) &
         MHD_SUM_p=MHD_SUM_p/MHD_SUM_vol

    where(MHD_SUM_vol>0.) &
         MHD_SUM_rho=MHD_SUM_rho/MHD_SUM_vol

    !Set volume floor
    MHD_SUM_vol = max(1.E-8,MHD_SUM_vol)

    if(dbg0)then
       if(iProc==0)write(*,*)' -E-'
    end if

    ! If the field-line tracer returned a good value, we may need to
    ! add contribution to V from the part of the field line that
    ! goes inside the body. If the field-line tracer did not
    ! return a good value, we will compute total V assuming dipole.
    Ri=(6378.+100.)/6378.
    Factor = 2. * (Ri**4) / abs(Bdp)
    do i=1,isize
       Colat = (90.0 - RCM_lat(i))*cDegToRad
       s2    = sin(colat)**2
       s6    = s2**3
       s8    = s6*s2

       Ci=abs(cos(colat))

       if( s2 < Ri/Rbody )then
          !Fieldline goes beyond Rbody, add piece of fieldline volume
          Cs=sqrt(1.-(Rbody/Ri)*s2)
          FCiCs = (Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) &
               - (1./7.)*(Ci**7-Cs**7)
          !!!CHANGE!!!
          if (s8 /= 0.0) then
             Vol = Factor*FCiCs/s8
          else
             Vol = 0.0
          endif
          where(MHD_SUM_vol(i,:)>1.1E-8)
             MHD_SUM_vol(i,:)=MHD_SUM_vol(i,:) + Vol
          end where
       end if

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values

          !Compute the full analytic volume
          FCiCs = Ci - Ci**3 + (3./5.)*Ci**5 - (1./7.)*Ci**7
          !!!CHANGE!!!
          if (s8 /= 0.0) then
             Vol = factor*FCiCs/s8
          else
             Vol = 0.0
          endif

          !Compute equatorial B value for dipole at this latitude
          eqB = abs(Bdp)*s6/Ri**3

          if( s2 > Ri/Rbody )then
             ! Fieldline stays inside of Rbody

             ! Recompute exact volume
             MHD_SUM_vol(i,:)=Vol

             ! Fix the grid inside Rbody
             MHD_Xeq(i,:) = (Ri/s2)*cos(RCM_lon(:)*cDegToRad)
             MHD_Yeq(i,:) = (Ri/s2)*sin(RCM_lon(:)*cDegToRad)

             ! Fix the equatorial B value
             MHD_Beq(i,:) = eqB

          else
             ! Fieldline stays inside of Rcurrents but outside Rbody
             ! Weight analytic formula proportional to the distance 
             ! of the field line from rBody within the rBody-rCurrents range

             Factor1= (Ri/Rbody - s2)/ (Ri/Rbody - Ri/Rcurrents)
             Factor2= 1.0 - Factor1

             ! Check if numerical volume exists
             where(MHD_SUM_vol(i,:)>1.1E-8)
                ! Blend numerical volume with exact volume
                MHD_SUM_vol(i,:) = Factor1*MHD_SUM_vol(i,:) + Factor2*Vol
             end where

          end if
       end if
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -F-'
    end if

    ! find index for latitude closest to equator
    iLoc_I   = minloc(abs(RCM_lat))
    iEquator = iLoc_I(1)

    ! find index for latitude closest to north pole
    iLoc_I = maxloc(RCM_lat)
    iNorthPole = iLoc_I(1)

    !set open fieldline values
    do j = 1, jSize

       ! Initialize the index of the last closed field line to be at
       ! the highest latitude in the RCM grid
       i0 = iNorthPole

       do i = iEquator, iNorthPole, sign(1, iNorthPole - iEquator)

          ! Find first open or very stretched field line 
          if(MHD_SUM_vol(i,j) < 1.1E-8 .or. &
               abs(MHD_Xeq(i,j)) > 200.0 .or. &
               abs(MHD_Yeq(i,j)) > 200.0 ) then
             i0 = i+1
             EXIT
          end if
       end do

       ! Save the index of the "last" closed field line into MHD_lat_boundary
       MHD_lat_boundary(j) = i0
!       write(*,*) "finding closed field-line ",j,i0,MHD_lat_boundary(j)

    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -H-'
    end if

    ! Set impossible values for open fieldlines
    ! except for Xeq and Yeq where the last closed values are used
    ! which is useful when the equatorial grid is plotted
    do j=1,jsize
       i = int(MHD_lat_boundary(j))
       MHD_SUM_vol(1:i-1,j) = -1.
       MHD_SUM_rho(1:i-1,j) = -1.
       MHD_SUM_p  (1:i-1,j) = -1.
       MHD_Beq    (1:i-1,j) = -1.
       MHD_Xeq    (1:i-1,j) = MHD_Xeq(i,j)
       MHD_Yeq    (1:i-1,j) = MHD_Yeq(i,j)
    end do

    !dimensionalize values
    where(MHD_SUM_vol > 0.)
       MHD_SUM_vol = MHD_SUM_vol / No2Si_V(UnitB_)
       MHD_SUM_rho = MHD_SUM_rho * No2Si_V(UnitRho_)
       MHD_SUM_p   = MHD_SUM_p   * No2Si_V(UnitP_)
    elsewhere
       MHD_SUM_vol = -1.
       MHD_SUM_rho = -1.
       MHD_SUM_p   = -1.
    end where

    where(MHD_Beq > 0.)
       MHD_Beq = MHD_Beq * No2Si_V(UnitB_)
    elsewhere
       MHD_Beq = -1.
    end where

    if(dbg0)then
       if(iProc==0)write(*,*)' -I-'
    end if

    if(DoTestTec .or. DoTestIdl)then
       do j=1,jsize-1; do i=1,isize-1
          if ( MHD_SUM_vol(i  ,j) > 0.0 .AND. &
               MHD_SUM_vol(i+1,j) > 0.0 .AND. &
               MHD_SUM_vol(i,j+1) > 0.0 .AND. &
               MHD_SUM_vol(i+1,j+1)> 0.0) THEN
             MHD_Fluxerror(i,j) = &
                  0.25E+9*(MHD_Beq(i,j  ) + MHD_Beq(i+1,j  ) + &
                  MHD_Beq(i,j+1) + MHD_Beq(i+1,j+1)) * &
                  0.5*(ABS((MHD_Xeq(i,j+1)-MHD_Xeq(i,j))* &
                  (MHD_Yeq(i+1,j)-MHD_Yeq(i,j))  &
                  - (MHD_Yeq(i,j+1)-MHD_Yeq(i,j)) &
                  * (MHD_Xeq(i+1,j)-MHD_Xeq(i,j))) &
                  + ABS((MHD_Xeq(i+1,j)-MHD_Xeq(i+1,j+1)) &
                  *(MHD_Yeq(i,j+1)-MHD_Yeq(i+1,j+1)) &
                  -(MHD_Yeq(i+1,j)-MHD_Yeq(i+1,j+1)) &
                  *(MHD_Xeq(i,j+1)-MHD_Xeq(i+1,j+1))))/&
                  (ABS(DipoleStrengthSi)*(SIN(Rcm_lat(i)*cDegToRad)**2 &
                  -SIN(Rcm_lat(i+1)*cDegToRad)**2)* &
                  (RCM_lon(j+1)-RCM_lon(j))*cDegToRad )- 1.0
          ELSE
             MHD_Fluxerror (i,j) = 0.0
          END IF
       end do; end do
       MHD_fluxerror(:,jsize) = MHD_Fluxerror (:,1)
       MHD_fluxerror(isize,:) = MHD_FLuxerror(isize-1,:)
    end if
    if(dbg0)then
       if(iProc==0)write(*,*)' -J-'
    end if

  end subroutine process_integrated_data

end module ModGmImCoupling
