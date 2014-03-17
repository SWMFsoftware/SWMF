!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
!^CMP FILE IM
!==========================================================================
subroutine GM_put_from_im(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  !call stop_mpi('RCM is OFF')

  use CON_coupler
  use CON_world,      ONLY: get_comp_info
  use CON_comp_param, ONLY: lNameVersion
  use ModImPressure                              ! Storage for IM pressure
  use ModNumConst
  use ModMain, ONLY : n_step,time_simulation, DoMultiFluidIMCoupling, &
       DoAnisoPressureIMCoupling
  use ModIoUnit, ONLY: UNITTMP_
  use ModProcMH, ONLY: iProc
  implicit none
  CHARACTER (LEN=80) :: filename
  character(len=*), parameter :: NameSub='GM_put_from_im'

  integer, intent(in) :: iSizeIn,jSizeIn,nVar
  real, intent(in) :: Buffer_IIV(iSizeIn,jSizeIn,nVar)
  character(len=*), intent(in) :: NameVar
  character(len=lNameVersion) :: NameVersionIm
  integer :: nCells_D(2), iError, i,j
  integer, parameter :: pres_=1, dens_=2, parpres_=3, bmin_=4, &
       Hpres_=3,Opres_=4,Hdens_=5,Odens_=6
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoMultiFluidIMCoupling)then
     if(NameVar /= 'p:rho:Hpp:Opp:Hprho:Oprho') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  else if(DoAnisoPressureIMCoupling)then
     if(NameVar /= 'p:rho:ppar:bmin') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  else
     if(NameVar /= 'p:rho') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end if

  nCells_D=ncells_decomposition_d(IM_)
  if( iSizeIn /= nCells_D(1) .or. jSizeIn /= nCells_D(2) ) then

     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSizeIn,jSizeIn,nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  if(.not.allocated(IM_lat))then
     ! Allocate IM_lat, IM_lon, IM_p, IM_dens
     call im_pressure_init(iSizeIn, jSizeIn)
     ! Set up IM ionospheric grid and store.
     ! Latitude specification is module specific, we must set it up
     ! according to the IM module we have selected.
     ! Determine version of IM:
     call get_comp_info(IM_, NameVersion=NameVersionIm)
     if(NameVersionIm(1:3) == 'RAM')then
        ! HEIDI and RAM-SCB have similar equatorial grids.
        IM_lat = Grid_C(IM_) % Coord1_I
     else
        ! RCM uses a CoLat based grid, information is stored in 
        ! module grid information.
        IM_lat = (cHalfPi - Grid_C(IM_) % Coord1_I) * cRadToDeg
     end if
     IM_lon = Grid_C(IM_)% Coord2_I * cRadToDeg
  end if

  ! Store IM variable for internal use
  IM_p    = Buffer_IIV(:,:,pres_)
  IM_dens = Buffer_IIV(:,:,dens_)
  iNewPIm  = iNewPIm + 1

  ! for multifluid                                               
  if(DoMultiFluidIMCoupling)then
     IM_Hpp = Buffer_IIV(:,:,Hpres_)
     IM_Opp = Buffer_IIV(:,:,Opres_)
     IM_Hpdens = Buffer_IIV(:,:,Hdens_)
     IM_Opdens = Buffer_IIV(:,:,Odens_)
  endif

  ! for anisotropic pressure
  if(DoAnisoPressureIMCoupling)then
     IM_ppar = Buffer_IIV(:,:,parpres_)
     IM_bmin = Buffer_IIV(:,:,bmin_)
  end if

  if(DoTest)call write_IMvars_tec  ! TecPlot output
  if(DoTest)call write_IMvars_idl  ! IDL     output

contains

  !============================================================================
  subroutine write_IMvars_tec
    integer :: j2
    real :: lonShift
    !-------------------------------------------------------------------------
    if(iProc /= 0)RETURN

    !write values to plot file
    write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".dat"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    if(DoMultiFluidIMCoupling)then
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat",&                
            &"IM pressure", "IM density", &                                   
            &"IM Hp pressure", "IM Hp density", &                             
            &"IM Op pressure", "IM Op density"'
    else if(DoAnisoPressureIMCoupling)then
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat",&                
            &"IM pressure", "IM density", &
            &"IM parallel pressure", "IM minimum B"'
    else
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat",&                
            &"IM pressure", "IM density"'
    end if
    write(UNITTMP_,'(a,i4,a,i4,a)') &
         'ZONE T="IM Pressure", I=',jSizeIn+1,', J=',iSizeIn,', K=1, F=POINT'
    do i=1,iSizeIn
       do j2=1,jSizeIn+1
          j=j2; if(j2==jSizeIn+1) j=1
          lonShift=0.; if(j2==jSizeIn+1) lonShift=360.
          if(DoMultiFluidIMCoupling)then
             write(UNITTMP_,'(2i4,8G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i), &
                  IM_p(i,j),IM_dens(i,j), &
                  IM_Hpp(i,j),IM_Hpdens(i,j),IM_Opp(i,j), IM_Opdens(i,j)
          else if(DoAnisoPressureIMCoupling)then
             write(UNITTMP_,'(2i4,6G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i), &
                  IM_p(i,j),IM_dens(i,j),IM_ppar(i,j),IM_bmin(i,j)
          else
             write(UNITTMP_,'(2i4,4G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i), &
                  IM_p(i,j),IM_dens(i,j)
          endif
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_IMvars_tec

  !============================================================================
  subroutine write_IMvars_idl
    if(iProc /= 0)RETURN

    !write values to plot file
    write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".out"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open file "//filename)
    write(UNITTMP_,'(a79)')            'IM pressure_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,2
    write(UNITTMP_,'(3i4)')            jSizeIn,iSizeIn
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    if(DoMultiFluidIMCoupling)then
       write(UNITTMP_,'(a79)')'Lon Lat p rho Hpp Hprho Opp Oprho nothing'
    else if(DoAnisoPressureIMCoupling)then
       write(UNITTMP_,'(a79)')'Lon Lat p rho ppar bmin nothing'
    else
       write(UNITTMP_,'(a79)')'Lon Lat p rho nothing'
    endif
     do i=iSizeIn,1,-1
       do j=1,jSizeIn
          if(DoMultiFluidIMCoupling)then
             write(UNITTMP_,'(100(1pe18.10))') &
                  IM_lon(j),IM_lat(i),IM_p(i,j),IM_dens(i,j), &
                  IM_Hpp(i,j), IM_Hpdens(i,j), IM_Opp(i,j), IM_Opdens(i,j)
          else if(DoAnisoPressureIMCoupling)then
             write(UNITTMP_,'(100(1pe18.10))') &
                  IM_lon(j),IM_lat(i),IM_p(i,j),IM_dens(i,j), &
                  IM_ppar(i,j), IM_bmin(i,j)
          else
             write(UNITTMP_,'(100(1pe18.10))') &
                  IM_lon(j),IM_lat(i),IM_p(i,j),IM_dens(i,j)
          endif
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_IMvars_idl



end subroutine GM_put_from_im
