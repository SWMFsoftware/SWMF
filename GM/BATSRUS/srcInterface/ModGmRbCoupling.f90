!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
module ModGmRbCoupling

  use ModMpi
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use CON_coupler, ONLY: Grid_C, RB_, ncells_decomposition_d

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,nBlockMax,unusedBLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace, ONLY : ray,rayface
  use ModPhysics, ONLY: unitSI_p, unitSI_rho, unitSI_temperature, unitSI_b, &
       Bdp, Bdp_dim, rCurrents, rBody
  implicit none

  character (len=*), parameter :: NameMod='ModGmRbCoupling'

  ! RB Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the RB grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RB_lat, RB_lon

  integer :: i,j,k, i0,i1,i2, j0,j1,j2, n, iBLK

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_SUM_vol, &
       MHD_SUM_rho, &
       MHD_SUM_p, &
       MHD_Beq, &
       MHD_Xeq, &
       MHD_Yeq
  real, parameter :: noValue=-99999.
  real :: qb(3),eqB,xL,yL,zL
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s6,s8,factor1,factor2

  integer, parameter :: maxMessages=10
  integer :: itag, NewMsg, iError
  integer :: nRECVrequests, RECVrequests(maxMessages), &
       nSENDrequests, SENDrequests(maxMessages), &
       MESGstatus(MPI_STATUS_SIZE, maxMessages)  

  logical :: dbg0=.false.

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

contains

  subroutine allocate_gm_rb(iSizeIn,jSizeIn)

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_rb'

    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(RB_lat))then
       nCells_D=ncells_decomposition_d(RB_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RB_lat(iSize), RB_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       RB_lat = 90.0 - Grid_C(RB_) % Coord1_I * cRadToDeg
       RB_lon =        Grid_C(RB_) % Coord2_I * cRadToDeg
    end if

    if(.not.allocated(MHD_SUM_vol))then
       allocate( MHD_SUM_vol(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_vol")
       MHD_SUM_vol = 0.
    end if

    if(.not.allocated(MHD_SUM_rho))then
       allocate( MHD_SUM_rho(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_rho")
       MHD_SUM_rho = 0.
    end if

    if(.not.allocated(MHD_SUM_p))then
       allocate( MHD_SUM_p(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_p")
       MHD_SUM_p = 0.
    end if

    if(.not.allocated(MHD_Beq))then
       allocate( MHD_Beq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Beq")
       MHD_Beq = 0.
    end if

    if(.not.allocated(MHD_Xeq))then
       allocate( MHD_Xeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Xeq")
       MHD_Xeq = 0.
    end if

    if(.not.allocated(MHD_Yeq))then
       allocate( MHD_Yeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Yeq")
       MHD_Yeq = 0.
    end if

    if(.not.allocated(MHD_lat_boundary))then
       allocate( MHD_lat_boundary(jsize), stat=iError )
       call alloc_check(iError,"MHD_lat_boundary")
       MHD_lat_boundary = 0
    end if

  end subroutine allocate_gm_rb

  !============================================================================
  subroutine write_integrated_data_tec
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: j2, nCall=0
    real :: tmpT, tmpV1,tmpV2, lonShift
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
         ', "Xeq", "Yeq"', &
         ', "Volume", "Volume**(-2/3)"', &
         ', "MHD `r", "MHD p", "MHD T", "Beq"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)/unitSI_p)/(MHD_SUM_rho(i,j)/unitSI_rho)) &
               * unitSI_temperature
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          write(UNITTMP_,'(2i4,12G14.6)') j2,i,RB_lon(j)+lonShift,RB_lat(i), &
               MHD_lat_boundary(j), &
               MHD_Xeq(i,j),MHD_Yeq(i,j), &
               tmpV1,tmpV2, &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j)
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
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,6
    write(UNITTMP_,'(3i4)')            jSize+1,iSize
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)')            'Lon Lat Xeq Yeq vol rho p Beq nothing'
    do i=isize,1,-1
       do j=1,jsize
          write(UNITTMP_,'(100(1pe18.10))') &
               RB_lon(j),RB_lat(i), &
               MHD_Xeq(i,j),MHD_Yeq(i,j),&
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Beq(i,j)
       end do
       write(UNITTMP_,'(100(1pe18.10))') &
               RB_lon(1)+360.0, RB_lat(i), &
               MHD_Xeq(i,1),MHD_Yeq(i,1),&
               MHD_SUM_vol(i,1), &
               MHD_SUM_rho(i,1),MHD_SUM_p(i,1),MHD_Beq(i,1)
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_idl

  !===========================================================================
  subroutine process_integrated_data

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
       Colat = (90.0 - RB_lat(i))*cDegToRad
       s2    = sin(colat)**2
       s6    = s2**3
       s8    = s6*s2

       Ci=cos(colat)

       if( s2 < Ri/Rbody )then
          !Fieldline goes beyond Rbody, add piece of fieldline volume
          Cs=sqrt(1.-(Rbody/Ri)*s2)
          FCiCs = (Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) &
               - (1./7.)*(Ci**7-Cs**7)
          Vol   = Factor*FCiCs/s8

          where(MHD_SUM_vol(i,:)>1.1E-8)
             MHD_SUM_vol(i,:)=MHD_SUM_vol(i,:) + Vol
          end where
       end if

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values

          !Compute the full analytic volume
          FCiCs = Ci - Ci**3 + (3./5.)*Ci**5 - (1./7.)*Ci**7
          Vol   = factor*FCiCs/s8

          !Compute equatorial B value for dipole at this latitude
          eqB = abs(Bdp)*s6/Ri**3

          if( s2 > Ri/Rbody )then
             ! Fieldline stays inside of Rbody

             ! Recompute exact volume
             MHD_SUM_vol(i,:)=Vol

             ! Fix the grid inside Rbody
             MHD_Xeq(i,:) = (Ri/s2)*cos(RB_lon(:)*cDegToRad)
             MHD_Yeq(i,:) = (Ri/s2)*sin(RB_lon(:)*cDegToRad)

             ! Fix the equatorial B value
             MHD_Beq(i,:) = eqB

          else
             ! Fieldline stays inside of Rcurrents but outside Rbody
             ! Weight analytic formula proportional to the distance 
             ! of the field line from rBody within the rBody-rCurrents range

             Factor1= (Ri/Rbody - s2)/ (Ri/Rbody - Ri/Rcurrents)
             Factor2= 1.0 - factor1

             ! Check if numerical volume exists
             where(MHD_SUM_vol(i,:)>1.1E-8)
                ! Blend numerical volume with exact volume
                MHD_SUM_vol(i,:) = Factor1*MHD_SUM_vol(i,:) + Factor2*Vol
             elsewhere
                ! Use analytic volume
                MHD_SUM_vol(i,:)=Vol

                ! Fix the grid
                MHD_Xeq(i,:) = (Ri/s2)*cos(RB_lon(:)*cDegToRad)
                MHD_Yeq(i,:) = (Ri/s2)*sin(RB_lon(:)*cDegToRad)

                ! Fix the equatorial B value
                MHD_Beq(i,:) = eqB
             end where

          end if
       end if
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -F-'
    end if

    !set open fieldline values
    do j=1,jsize
       i0=isize
       do i=isize,1,-1
          if(MHD_SUM_vol(i,j)<1.1E-8 .OR. &
               MHD_Xeq(i,j)<=noValue .OR. &
               MHD_Yeq(i,j)<=noValue ) then
             i0=i+1
             exit
          end if
       end do

       i=i0
       MHD_lat_boundary(j)=i

       MHD_SUM_vol(1:i-1,j)=MHD_SUM_vol(i,j)
       MHD_Beq    (1:i-1,j)=MHD_Beq    (i,j)
       MHD_Xeq    (1:i-1,j)=MHD_Xeq    (i,j)
       MHD_Yeq    (1:i-1,j)=MHD_Yeq    (i,j)
    end do

    if(dbg0)then
       if(iProc==0)write(*,*)' -G-'
    end if

    !error checking
    do j=1,jsize; do i=isize,1,-1
       IF (iProc ==0) then
          if (mhd_sum_vol(i,j) > 1.0E-8 .and. &
               (MHD_Xeq(i,j) < -999. .or. MHD_Yeq(i,j) < -999.)) then
             print*,'999999:',i,j,MHD_Xeq(i,j),MHD_Yeq(i,j), mhd_sum_vol(i,j)
             call stop_mpi('stanislav')
          end if
          if (mhd_lat_boundary(j) < i) then
             if (mhd_sum_vol(i,j) < 1.1E-8) then
                print*,'7777',i,j,mhd_lat_boundary(j),mhd_sum_vol(i,j), &
                     MHD_Xeq(i,j),MHD_Yeq(i,j)
                CALL STOP_MPI ('HOLE ')
             end if
          end if
       end if
    end do; enddo

    if(dbg0)then
       if(iProc==0)write(*,*)' -H-'
    end if

    !set values for open fieldlines
    do j=1,jsize
       MHD_SUM_vol(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_rho(1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_SUM_p  (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Beq    (1:int(MHD_lat_boundary(j))-1,j)=-1.
       MHD_Xeq    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Xeq(int(MHD_lat_boundary(j)),j)
       MHD_Yeq    (1:int(MHD_lat_boundary(j))-1,j)=MHD_Yeq(int(MHD_lat_boundary(j)),j)
    end do

    !dimensionalize values
    where(MHD_SUM_vol>0.)
       MHD_SUM_vol = MHD_SUM_vol  / unitSI_B
       MHD_SUM_rho = MHD_SUM_rho * unitSI_rho
       MHD_SUM_p   = MHD_SUM_p * unitSI_p
    elsewhere
       MHD_SUM_vol = -1.
       MHD_SUM_rho = -1.
       MHD_SUM_p   = -1.
    end where

    where(MHD_Beq>0.)
       MHD_Beq = MHD_Beq * unitSI_B
    elsewhere
       MHD_Beq = -1.
    end where

    if(dbg0)then
       if(iProc==0)write(*,*)' -I-'
    end if

    if(dbg0)then
       if(iProc==0)write(*,*)' -J-'
    end if

  end subroutine process_integrated_data

end module ModGmRbCoupling
