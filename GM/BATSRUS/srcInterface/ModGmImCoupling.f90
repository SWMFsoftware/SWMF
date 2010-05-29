!^CFG COPYRIGHT UM
!^CMP FILE IM
!^CFG FILE RAYTRACE
module ModGmImCoupling

  use ModMpi
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use CON_coupler, ONLY: Grid_C, IM_, ncells_decomposition_d

  use ModProcMH
  use ModMain, ONLY: nI,nJ,nK,n_step,nBlockMax,unusedBLK,DoMultiFluidIMCoupling
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace, ONLY : ray,rayface
  use ModPhysics, ONLY: No2Si_V, Si2No_V, &
       UnitP_, UnitRho_, UnitTemperature_, UnitB_, &
       Bdp, DipoleStrengthSi, rCurrents, rBody
  implicit none

  character(len=*), parameter :: NameMod='ModGmImCoupling'

  ! IM Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the IM grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RCM_lat, RCM_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_SUM_vol, MHD_tmp, &
       MHD_SUM_rho, MHD_HpRho, MHD_OpRho, &
       MHD_SUM_p, MHD_HpP, MHD_OpP, &
       MHD_Beq, &
       MHD_Xeq, &
       MHD_Yeq, &
       MHD_Fluxerror
  real, parameter :: noValue=-99999.
  real :: eqB,xL,yL,zL
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s6,s8,factor1,factor2

  integer :: iError

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

  logical :: DoTestTec, DoTestIdl

  ! This is for the GM-IM/RAM coupling                                      
  real, allocatable :: StateLine_VI(:,:)

contains

  subroutine allocate_gm_im(iSizeIn,jSizeIn)
    use CON_comp_param, ONLY: IM_

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_im'

    if(allocated(MHD_lat_boundary)) deallocate(MHD_lat_boundary)
    if(allocated(MHD_SUM_vol))      deallocate(MHD_SUM_vol)
    if(allocated(MHD_SUM_rho))      deallocate(MHD_SUM_rho)
    if(allocated(MHD_SUM_p))        deallocate(MHD_SUM_p)
    if(allocated(MHD_tmp))          deallocate(MHD_tmp)
    if(allocated(MHD_Beq))          deallocate(MHD_Beq)
    if(allocated(MHD_Xeq))          deallocate(MHD_Xeq)
    if(allocated(MHD_Yeq))          deallocate(MHD_Yeq)
    if(allocated(MHD_Fluxerror))    deallocate(MHD_Fluxerror)
    
    if(DoMultiFluidIMCoupling)then
       if(allocated(MHD_HpRho))        deallocate(MHD_HpRho)
       if(allocated(MHD_HpP))          deallocate(MHD_HpP)
       if(allocated(MHD_OpRho))        deallocate(MHD_OpRho)
       if(allocated(MHD_OpP))          deallocate(MHD_OpP)
    end if
    
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

    if(DoMultiFluidIMCoupling)then
       allocate( MHD_Hprho(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Hprho")
       MHD_Hprho = 0.
       
       allocate( MHD_Oprho(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Oprho")
       MHD_Oprho = 0.
       
       allocate( MHD_Hpp(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Hpp")
       MHD_Hpp = 0.
       
       allocate( MHD_Opp(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Opp")
       MHD_Opp = 0.
    end if
    
    allocate( MHD_Beq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Beq")
    MHD_Beq = 0.

    allocate( MHD_Xeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Xeq")
    MHD_Xeq = 0.

    allocate( MHD_Yeq(isize,jsize), stat=iError )
    call alloc_check(iError,"MHD_Yeq")
    MHD_Yeq = 0.

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
    real :: tmpT, tmpV1,tmpV2, lonShift,tmpHpT,tmpOpT
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    if(.not.DoMultiFluidIMCoupling)then
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
            ', "Xeq", "Yeq"', &
            ', "Volume", "Volume**(-2/3)"', &
            ', "MHD `r", "MHD p", "MHD T", "Beq"', &
            ', "FluxError"'
    else
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
            ', "Xeq", "Yeq"', &
            ', "Volume", "Volume**(-2/3)"', &
            ', "MHD `r", "MHD p", "MHD T"', &
            ', "MHD Hp`r", "MHD Hp p", "MHD Hp T"', &
            ', "MHD Op`r", "MHD Op p", "MHD Op T", "Beq"', &
            ', "FluxError"'
    end if
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)*Si2No_V(UnitP_))/(MHD_SUM_rho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
          if(DoMultiFluidIMCoupling)then
             tmpHpT=-1.; if(MHD_Hprho(i,j)>0.) &
                  tmpHpT = ((MHD_Hpp(i,j)*Si2No_V(UnitP_))/(MHD_Hprho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
             tmpOpT=-1.; if(MHD_Oprho(i,j)>0.) &
                  tmpOpT = ((MHD_Opp(i,j)*Si2No_V(UnitP_))/(MHD_Oprho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
          end if
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          if(.not.DoMultiFluidIMCoupling)then
             write(UNITTMP_,'(2i4,12G14.6)') j2,i,RCM_lon(j)+lonShift,RCM_lat(i),&
                  MHD_lat_boundary(j), &
                  MHD_Xeq(i,j),MHD_Yeq(i,j), &
                  tmpV1,tmpV2, &
                  MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j), &
                  MHD_Fluxerror(i,j)
          else
             write(UNITTMP_,'(2i4,18G14.6)') j2,i,RCM_lon(j)+lonShift,RCM_lat(i),&
                  MHD_lat_boundary(j), &
                  MHD_Xeq(i,j),MHD_Yeq(i,j), &
                  tmpV1,tmpV2, &
                  MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,&
                  MHD_Hprho(i,j), MHD_Hpp(i,j), tmpHpT, &
                  MHD_Oprho(i,j), MHD_Opp(i,j), tmpOpT,MHD_Beq(i,j), &
                  MHD_Fluxerror(i,j)
          end if
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

    integer :: iLoc_I(1),iEquator,iNorthPole

    if(.not.DoMultiFluidIMCoupling)then
       where(MHD_SUM_vol>0.)
          MHD_SUM_p   = MHD_SUM_p/MHD_SUM_vol
          MHD_SUM_rho = MHD_SUM_rho/MHD_SUM_vol
       end where
    else
       where(MHD_SUM_vol>0.)
          MHD_SUM_p   = MHD_SUM_p/MHD_SUM_vol
          MHD_SUM_rho = MHD_SUM_rho/MHD_SUM_vol
          MHD_Hprho = MHD_Hprho/MHD_SUM_vol
          MHD_Oprho = MHD_Oprho/MHD_SUM_vol
          MHD_Hpp   = MHD_Hpp/MHD_SUM_vol
          MHD_Opp   = MHD_Opp/MHD_SUM_vol
       end where
    end if

    !Set volume floor
    MHD_SUM_vol = max(1.E-8,MHD_SUM_vol)

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

    ! find index for latitude closest to equator
    iLoc_I   = minloc(abs(RCM_lat))
    iEquator = iLoc_I(1)

    ! find index for latitude closest to north pole
    iLoc_I = maxloc(RCM_lat)
    iNorthPole = iLoc_I(1)

    iEquator = iEquator + sign(1, iNorthPole - iEquator) 
    !set open fieldline values
    do j=1,jsize

      ! Initialize the index of the last closed field line to be at
       ! the highest latitude in the RCM grid
       i0 = iNorthPole

       do i = iEquator, iNorthPole, sign(1, iNorthPole - iEquator)

          if(MHD_SUM_vol(i,j) < 1.1E-8 .or. &
               abs(MHD_Xeq(i,j)) > 200.0 .or. &
               abs(MHD_Yeq(i,j)) > 200.0 ) then
             i0=i+1
             exit
          end if
       end do

       ! Save the index of the "last" closed field line into MHD_lat_boundary
       MHD_lat_boundary(j)=i0
      ! write(*,*) "finding closed field-line ",j,i0,MHD_lat_boundary(j)

    end do

    ! Set impossible values for open fieldlines
    ! except for Xeq and Yeq where the last closed values are used
    ! which is useful when the equatorial grid is plotted
    do j=1,jsize
       i = int(MHD_lat_boundary(j))
       MHD_Beq    (1:i-1,j) = -1.
       MHD_Xeq    (1:i-1,j) = MHD_Xeq(i,j)
       MHD_Yeq    (1:i-1,j) = MHD_Yeq(i,j)
       MHD_SUM_vol(1:i-1,j) = -1.
    
       if(.not.DoMultiFluidIMCoupling)then
          MHD_SUM_rho(1:i-1,j) = -1.
          MHD_SUM_p  (1:i-1,j) = -1.
       else
          MHD_SUM_rho(1:i-1,j) = -1.
          MHD_SUM_p  (1:i-1,j) = -1.
          MHD_HpRho(1:i-1,j) = -1.
          MHD_OpRho(1:i-1,j) = -1.
          MHD_HpP  (1:i-1,j) = -1.
          MHD_OpP  (1:i-1,j) = -1.
       endif
    end do

    ! Dimensionalize values
    ! Note: dimensions of "MHD_sum_vol" is Distance/Magnetic field.
    ! The distance unit is planetary radius in GM, and the same is
    ! used in RCM, so only the magnetic unit is converted to SI units.
    ! Similarly Xeq and Yeq (equatorial crossing coords) remain in 
    ! normalized units.
    if(.not. DoMultiFluidIMCoupling)then
       where(MHD_SUM_vol > 0.)
          MHD_SUM_vol = MHD_SUM_vol / No2Si_V(UnitB_)
          MHD_SUM_rho = MHD_SUM_rho * No2Si_V(UnitRho_)
          MHD_SUM_p   = MHD_SUM_p   * No2Si_V(UnitP_)
       elsewhere
          MHD_SUM_vol = -1.
          MHD_SUM_rho = -1.
          MHD_SUM_p   = -1.
       end where
    else
       where(MHD_SUM_vol > 0.)
          MHD_SUM_vol = MHD_SUM_vol / No2Si_V(UnitB_)
          MHD_SUM_rho = MHD_SUM_rho * No2Si_V(UnitRho_)
          MHD_SUM_p   = MHD_SUM_p   * No2Si_V(UnitP_)
          MHD_HpRho = MHD_HpRho * No2Si_V(UnitRho_)
          MHD_OpRho = MHD_OpRho * No2Si_V(UnitRho_)
          MHD_HpP   = MHD_HpP * No2Si_V(UnitP_)
          MHD_OpP   = MHD_OpP * No2Si_V(UnitP_)
       elsewhere
          MHD_SUM_vol = -1.
          MHD_SUM_rho = -1.
          MHD_SUM_p   = -1.
          MHD_HpRho = -1.
          MHD_OpRho = -1.
          MHD_HpP   = -1.
          MHD_OpP   = -1.
       end where
    endif

    where(MHD_Beq > 0.)
       MHD_Beq = MHD_Beq * No2Si_V(UnitB_)
    elsewhere
       MHD_Beq = -1.
    end where

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

  end subroutine process_integrated_data

end module ModGmImCoupling
