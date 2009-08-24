! Wrapper for Internal Magnetosphere (IM) component
!=============================================================================
subroutine IM_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use ModProcIM
  use RCM_variables, ONLY: NameRcmDir, iUnitOut, StringPrefix, STDOUT_, &
       DoRestart, iDtRcm, iDtPlot, asci_flag, nFilesPlot, iDnPlot, &
       plot_area, plot_var, plot_format, &
       x_h, x_o, L_dktime, sunspot_number, f107, doy, &
       ipot, ibnd_type, &
       precipitation_tau
  use ModReadParam
  use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case

  implicit none

  character (len=*), parameter :: NameSub='IM_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  !LOCAL VARIABLES:
  character (len=100) :: NameCommand, StringPlot
  character (len=20)  :: StringTypeIonosphere='', StringTypeOuterboundary='',&
                         LossFactorString=''
  logical             :: DoEcho=.false.
  logical             :: UseStrict=.true.
  integer :: iFile
  real :: FractionH,FractionO, SunspotNumber,F107MonthlyMean,DayOfYear,tau_in=0.0
  !-------------------------------------------------------------------------

  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                         &
          Use=.true.,                           &
          NameVersion='RCM (De Zeeuw-Sazykin)', &
          Version=2.0)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     if(nProc>1)call CON_stop(&
          'IM_init_mpi: IM_ERROR this version can run on 1 PE only!')
  case('CHECK')
     ! We should check and correct parameters here
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
     RETURN
  case('READ')
     if(iProc==0)write(*,*) NameSub,': READ iSession =',i_session_read(),&
          ' iLine=',i_line_read(),' nLine =',n_line_read()
     do
        if(.not.read_line() ) EXIT
        if(.not.read_command(NameCommand)) CYCLE

        select case(NameCommand)
        case("#STRICT")
           call read_var('UseStrict',UseStrict)
        case("#RCMDIR")
           call read_var('NameRcmDir',NameRcmDir)
           call fix_dir_name(NameRcmDir)
           if(iProc==0)call check_dir(NameRcmDir)
        case("#ASCII")
           call read_var('IsAscii',asci_flag)
        case("#RESTART")
           call read_var('DoRestart',DoRestart)
        case("#TIMESTEP")
           call read_var('iDtRcm',iDtRcm)
        case("#SAVEPLOT")
           call read_var('nPlotFile',nFilesPlot)
           do iFile=1,nFilesPlot
              call read_var('StringPlot',StringPlot)
              call lower_case(StringPlot)

              ! Plotting frequency
              call read_var('DnSavePlot',iDnPlot(iFile))
              call read_var('iDtSavePlot',iDtPlot(iFile))

              ! Extract geometry string value:
              if(index(StringPlot,'2d')>0)then
                 plot_area(ifile)='2d_'
              elseif(index(StringPlot,'3d')>0)then
                 plot_area(ifile)='3d_'
              else
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#SAVEPLOT geometry string (2d/3d) is missing')
              end if

              ! Extract variable string value:
              if(index(StringPlot,'min')>0)then
                 plot_var(ifile)='min'
              elseif(   index(StringPlot,'max')>0 &
                   .or. index(StringPlot,'rcm')>0 )then
                 plot_var(ifile)='max'
              else
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#SAVEPLOT variable string (min/max/rcm) is missing')
              end if

              ! Extract format string value:
              if(index(StringPlot,'idl')>0)then
                 plot_format(ifile)='idl'
              elseif(index(StringPlot,'tec')>0)then 
                 plot_format(ifile)='tec'
              else
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#SAVEPLOT file format (idl/tec) is missing')
              end if

           end do
        case("#COMPOSITION")
           call read_var('FractionH',FractionH)
           call read_var('FractionO',FractionO)
           x_h=FractionH; x_o=FractionO
        case("#CHARGEEXCHANGE")
           call read_var('UseChargeExchange',L_dktime)
           if(L_dktime)then
              call read_var('SunspotNumber',SunspotNumber)
              call read_var('F107MonthlyMean',F107MonthlyMean)
              call read_var('DayOfYear',DayOfYear)
              ! Convert from default real to real(r_prec)
              sunspot_number= SunspotNumber
              f107          = F107MonthlyMean
              doy           = DayOfYear
           end if
        case ("#IONOSPHERE")
              call read_var('TypeIonosphere',StringTypeIonosphere)
              select case (StringTypeIonosphere)
              case ('IE','Ie','ie','iE')
                 ipot = 6
              case ('RCM','rcm','Rcm')
                 ipot = 4
              case default
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#IONOSPHERE parameter has an illegal value')
              end select
        case ("#OUTERBOUNDARY")
              call read_var ('TypeOuterBoundary',StringTypeOuterboundary)
              select case (StringTypeOuterboundary)
                 case ('max','MAX','Max')
                    ibnd_type = 5
                 case ('ellipse','ELLIPSE','Ellipse')
                    ibnd_type = 4
                 case default
                    call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#OUTERBOUNDARY has an illegal value')
              end select
        case ("#PRECIPITATION")
              call read_var ('LossFactore',tau_in)
                   precipitation_tau(1) = tau_in
              call read_var ('LossFactorH',tau_in)
                   precipitation_tau(2) = tau_in
              call read_var ('LossFactorO',tau_in)
                   precipitation_tau(3) = tau_in
              if (precipitation_tau(1) < 0 .or. precipitation_tau(1) > 1.0)&
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#LossFactore has an illegal value')
              if (precipitation_tau(2) < 0 .or. precipitation_tau(2) > 1.0)&
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#LossFactorH has an illegal value')
              if (precipitation_tau(3) < 0 .or. precipitation_tau(3) > 1.0)&
                 call CON_stop('ERROR in IM/RCM2/src/IM_wrapper.f90: '// &
                      '#LossFactorO has an illegal value')
        case default
           if(iProc==0) then
              write(*,'(a,i4,a)')NameSub//' IM_ERROR at line ',i_line_read(),&
                   ' invalid command '//trim(NameCommand)
              if(UseStrict)call CON_stop('Correct PARAM.in!')
           end if
        end select
     end do
  case('STDOUT')
     iUnitOut=STDOUT_
     StringPrefix='IM:'
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''
  case('GRID')
     call IM_set_grid
  case default
     call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine IM_set_param

!============================================================================
subroutine IM_set_grid
  use ModProcIM
  use ModNumConst
  use CON_coupler, ONLY: set_grid_descriptor, is_proc, IM_
  use Rcm_variables, ONLY: iSize, jSize, colat, aloct, Ri
  use RCM_io, ONLY: read_grid
  implicit none
  character (len=*), parameter :: NameSub='IM_set_grid'
  real :: Radius_I(1)
  logical :: IsInitialized=.false.
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)'IM_set_grid_descriptor called, IsInitialized=',&
       IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  if(is_proc(IM_))call read_grid()

  if (nProc<0)then
     colat=0.; aloct=0.; Ri=1.
  end if

  Radius_I(1) = Ri*1000.0 ! radial size of the ionosphere in meters

  ! IM grid size in generalized coordinates
  call set_grid_descriptor( IM_,                 & ! component index
       nDim=2,                                   & ! dimensionality
       nRootBlock_D=(/1,1/),                     & ! single block
       nCell_D=(/iSize,jSize/),                  & ! size of cell based grid
       XyzMin_D=(/cHalf, cHalf/),                & ! min gen.coords for cells
       XyzMax_D=(/iSize+cHalf,jSize+cHalf/),     & ! max gen.coords for cells
       TypeCoord='SMG',                          & ! solar magnetic coord
       Coord1_I=real(colat(1:iSize,1)),          & ! colatitudes
       Coord2_I=real(aloct(1,1:jSize)),          & ! longitudes
       Coord3_I=Radius_I,                        & ! radial size in meters
       IsPeriodic_D=(/.false.,.true./))            ! periodic in longitude

end subroutine IM_set_grid
!==============================================================================
subroutine IM_print_variables(NameSource)

  use rcm_variables
  use ModNumConst
  use ModIoUnit, ONLY: UNITTMP_
  implicit none
  character(len=*), parameter :: NameSub='IM_print_variables'

  character(len=*),intent(in) :: NameSource
  integer            :: nFile=0
  character(len=100) :: NameFile
  character(len=100) :: NameVar
  integer            :: i,j
  real               :: Lat,Lon
  !--------------------------------------------------------------------------
  select case(NameSource)
  case('IE')
     NameVar='j i lon lat jr pot sigmaH sigmaP'
  case('GM')
     if(.not.DoMultiFluidGMCoupling)then
        NameVar='j i lon lat density pressure vm xmin ymin bmin temperature'
     else
        NameVar='j i lon lat Hpdensity Opdensity Hpressure Opressure vm xmin,&
             &ymin bmin Hptemperature Optemperature'
     end if
  case default
     write(*,*)NameSub,': incorrect NameSource=',NameSource
     RETURN
  end select

  nFile=nFile+1
  write(NameFile,'(a,i1,a)')'IM_from_'//NameSource//'_',nFile,'.dat'
  open(UNITTMP_,file=NameFile)
  write(UNITTMP_,'(a)')trim(NameVar)

  do i=1,iSize
     do j=1,jSize
        Lon = (        aloct(i,j))*(180./cPi)
        Lat = (cHalfPi-colat(i,j))*(180./cPi)
        select case(NameSource)
        case('IE')
           write(UNITTMP_,'(2i4,6G14.6)')j,i,Lon,Lat,v(i,j),birk_mhd(i,j),&
                sigmaH_mhd(i,j),sigmaP_mhd(i,j)
        case('GM')
           if(.not.DoMultiFluidGMCoupling)then
              write(UNITTMP_,'(2i4,9G14.6)')j,i,Lon,Lat,density(i,j),pressure(i,j),&
                   vm(i,j),xmin(i,j),ymin(i,j),bmin(i,j),temperature(i,j)
           else
              write(UNITTMP_,'(2i4,12G14.6)')j,i,Lon,Lat,densityHp(i,j),densityOp(i,j),&
                   pressureHp(i,j),pressureOp(i,j),vm(i,j),xmin(i,j),ymin(i,j), &
                   bmin(i,j),temperatureHp(i,j),temperatureOp(i,j)
           end if
        end select
     end do
  end do
  close(UNITTMP_)

end subroutine IM_print_variables

!==============================================================================
subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

  ! Provide current for IE
  ! The value should be interpolated from nPoints with
  ! indexes stored in Index and weights stored in Weight
  ! The variables should be put into Buff_V

  use RCM_variables, ONLY: birk, eavg, eflux, iSize, jSize, colat, aloct
  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  implicit none
  character(len=*), parameter :: NameSub='IM_get_for_ie'

  integer,intent(in)            :: nPoint, iPointStart, nVar
  real,intent(out)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight

  integer :: iLat, iLon, iBlock, iPoint
  real    :: w

  Buff_V = 0.0

  do iPoint = iPointStart, iPointStart + nPoint - 1

     iLat   = Index % iCB_II(1,iPoint)
     iLon   = Index % iCB_II(2,iPoint)
     iBlock = Index % iCB_II(3,iPoint)
     w      = Weight % Weight_I(iPoint)

     if(iBlock/=1)then
        write(*,*)NameSub,': iPoint,Index % iCB_II=',&
             iPoint,Index%iCB_II(:,iPoint)
        call CON_stop(NameSub//&
             ' SWMF_ERROR iBlock should be 1=North in IM-IE coupling')
     end if

     if (iLon > jSize) iLon = mod(iLon,jSize)
     if(iLat<1 .or. iLat>iSize*2 .or. iLon<0 .or. iLon>jSize)then
        write(*,*)'iLat,iLon=',iLat, iSize, iLon, jSize
        call CON_stop(NameSub//' SWMF_ERROR index out of range')
     end if

     ! Only worry about the northern hemisphere....  IE can fix the southern hemisphere.

     if (iLat <= iSize .and. iLon <= jSize) then
        Buff_V(1) = Buff_V(1) - w * birk(iLat,iLon)/2 / 1.0e6
        Buff_V(2) = Buff_V(2) + w * eflux(iLat,iLon,1)
        Buff_V(3) = Buff_V(3) + w * eavg(iLat,iLon,1) / 1000.0
     endif

!     if (iLat > iSize .and. iLon <= jSize) &
!          Buff_V(1) = Buff_V(1) + w * birk(2*iSize-iLat+1,iLon)/2 / 1.0e6

!     if (iLat > IONO_nTheta .and. iLon <= IONO_nPsi) &
!          Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_RCM_JR(2*IONO_nTheta-iLat+1,iLon)

  end do

!!! aloct = angle from noon on radians
!!! colat 
!!! birk - with ghost cells (positive down, summed over two hems microamps/m2)
!!! eflux - lat, lon, species (1 - electrons) (total between 2 hems ergs/cm2/s)
!!! eavg - lat, lon, species (1 - electrons) (keV)

end subroutine IM_get_for_ie

!============================================================================
subroutine IM_put_from_ie_mpi

  call CON_stop('IM_put_from_ie_mpi cannot be used by RCM2!')

end subroutine IM_put_from_ie_mpi

!==============================================================================
subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  use CON_router,   ONLY: IndexPtrType, WeightPtrType
  use RCM_variables, ONLY: v, birk_mhd, iSize, jSize, sigmaH_mhd,sigmaP_mhd

  implicit none
  character(len=*), parameter   :: NameSub='IM_put_from_ie'
  integer,intent(in)            :: nPoint, iPointStart, nVar
  real, intent(in)              :: Buff_V(nVar)
  type(IndexPtrType),intent(in) :: Index
  type(WeightPtrType),intent(in):: Weight
  logical,intent(in)            :: DoAdd
  integer :: iBlock,i,j
  !---------------------------------------------------------------------------
  if(nPoint>1)then
     write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
          nPoint,iPointStart,Weight % Weight_I
     call CON_stop(NameSub//': should be called with 1 point')
  end if
  if(DoAdd)then
     write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
          nPoint,iPointStart,Weight % Weight_I
     write(*,*)NameSub,': WARNING DoAdd is true'
  end if

  i = Index % iCB_II(1,iPointStart)
  j = Index % iCB_II(2,iPointStart)

  if(i<1.or.i>isize.or.j<1.or.j>jsize)then
     write(*,*)'i,j,DoAdd=',i,j,DoAdd
     call CON_stop('IM_put_from_ie: index out of range')
  end if

  if(DoAdd)then
     v(i,j)        = v(i,j)        + Buff_V(1)
     birk_mhd(i,j) = birk_mhd(i,j) + Buff_V(2)
     sigmaH_mhd(i,j) = sigmaH_mhd(i,j) + Buff_V(3)
     sigmaP_mhd(i,j) = sigmaP_mhd(i,j) + Buff_V(4)
  else
     v(i,j)        = Buff_V(1)
     birk_mhd(i,j) = Buff_V(2)
     sigmaH_mhd(i,j) = Buff_V(3)
     sigmaP_mhd(i,j) = Buff_V(4)
  end if

end subroutine IM_put_from_ie
!==============================================================================
subroutine IM_put_from_ie_complete

  use RCM_variables, ONLY: v, birk_mhd, iSize, jSize, sigmaH_mhd,sigmaP_mhd, n_gc
  implicit none

  call wrap_around_ghostcells(v,isize,jsize,n_gc)
  call wrap_around_ghostcells(birk_mhd,isize,jsize,n_gc)
  call wrap_around_ghostcells(sigmaH_mhd,isize,jsize,n_gc)
  call wrap_around_ghostcells(sigmaP_mhd,isize,jsize,n_gc)

end subroutine IM_put_from_ie_complete
!==============================================================================
subroutine IM_put_from_gm_line

  call CON_stop('IM_put_from_gm_line should not be called for IM/RCM2')

end subroutine IM_put_from_gm_line
!==============================================================================
subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

  use RCM_variables
  implicit none
  character (len=*),parameter :: NameSub='IM_put_from_gm'

  integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
  real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
  character (len=*),intent(in)       :: NameVar

  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6, &
       HpRho_=7, OpRho_=8,HpP_=9,OpP_=10
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(nVarIn > 6)DoMultiFluidGMCoupling = .true.
  if(DoTest)write(*,*)NameSub,' starting with NameVar=',NameVar
  if(DoTest) call write_data

  if(.not. DoMultiFluidGMCoupling)then
     if(NameVar /= 'vol:z0x:z0y:bmin:rho:p') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  else
     if(NameVar /= 'vol:z0x:z0y:bmin:rho:p:Hprho:Oprho:Hpp:Opp') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end if
 
  DoneGmCoupling = .true.
  if(.not. DoMultiFluidGMCoupling)then
     if(iSizeIn /= iSize .or. jSizeIn /= jSize .or. nVarIn /= p_)then
        write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn,nVarIn
        call CON_stop(NameSub//' SWMF_ERROR')
     end if
  else
     if(iSizeIn /= iSize .or. jSizeIn /= jSize .or. nVarIn /= OpP_)then
        write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn,nVarIn
        call CON_stop(NameSub//' SWMF_ERROR')
     end if
  end if
  vm(1:isize,1:jsize)   = Buffer_IIV(:,:,vol_)
  ! Convert GM volume into IM volume variable. Change units (m -> km) first
  vm(1:isize,1:jsize) = vm(1:isize,1:jsize) / 1.0e+9
  where(vm(1:isize,1:jsize)>0.) &
       vm(1:isize,1:jsize)=vm(1:isize,1:jsize)**(-2./3.) 
  xmin(1:isize,1:jsize) = Buffer_IIV(:,:,z0x_)
  ymin(1:isize,1:jsize) = Buffer_IIV(:,:,z0y_)
  bmin(1:isize,1:jsize) = Buffer_IIV(:,:,bmin_)
  density(1:isize,1:jsize) = Buffer_IIV(:,:,rho_)/xmass(2)/1.0E+6 ! in cm-3
  pressure(1:isize,1:jsize) = Buffer_IIV(:,:,p_)
  where(Buffer_IIV(:,:,rho_) /= 0.0) 
     temperature (1:iSize,1:jSize) = &
          Buffer_IIV(:,:,p_)/(Buffer_IIV(:,:,rho_)/xmass(2))/1.6E-19 ! in K
  elsewhere
     temperature (1:iSize,1:jSize) =5000.0
  end where

  call wrap_around_ghostcells(vm,isize,jsize,n_gc)
  call wrap_around_ghostcells(xmin,isize,jsize,n_gc)
  call wrap_around_ghostcells(ymin,isize,jsize,n_gc)
  call wrap_around_ghostcells(bmin,isize,jsize,n_gc)
  call wrap_around_ghostcells(density, isize, jsize, n_gc)
  call wrap_around_ghostcells(temperature, isize, jsize, n_gc)

  if(DoMultiFluidGMCoupling)then
     ! MultiFluid                                                                       
     densityHp(1:isize,1:jsize) = Buffer_IIV(:,:,HpRho_)/xmass(2)/1.0E+6 ! in cm-3       
     densityOp(1:isize,1:jsize) = Buffer_IIV(:,:,OpRho_)/xmass(3)/1.0E+6 ! in cm-3       
     pressureHp(1:isize,1:jsize)= Buffer_IIV(:,:,HpP_)
     pressureOp(1:isize,1:jsize)= Buffer_IIV(:,:,OpP_)

     where(Buffer_IIV(:,:,Hprho_) /= 0.0)
        temperatureHp (1:iSize,1:jSize) = &
             Buffer_IIV(:,:,Hpp_)/(Buffer_IIV(:,:,Hprho_)/xmass(2))/1.6E-19
        !in eV,assuming p=nkT, [p/n]=J,J/e-->eV                                       
     elsewhere
        temperatureHp (1:iSize,1:jSize) =5000.0
     end where
     where(Buffer_IIV(:,:,Oprho_) /= 0.0)
        temperatureOp (1:iSize,1:jSize) = &
             Buffer_IIV(:,:,Opp_)/(Buffer_IIV(:,:,Oprho_)/xmass(3))/1.6E-19 ! in ev     
     elsewhere
        temperatureOp (1:iSize,1:jSize) =5000.0
     end where

     call wrap_around_ghostcells(densityHp, isize, jsize, n_gc)
     call wrap_around_ghostcells(densityOp, isize, jsize, n_gc)
     call wrap_around_ghostcells(temperatureHp, isize, jsize, n_gc)
     call wrap_around_ghostcells(temperatureOp, isize, jsize, n_gc)

  endif

contains

  !============================================================================
  !write values sent to IM from GM
  subroutine write_data
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: i,j
    integer, save :: nCall=0
    !-------------------------------------------------------------------------

    nCall=nCall+1
    write(filename,'(a,i5.5,a)')"gm2im_debug_",nCall,".dat"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="gm2im debug values"'
    if(.not. DoMultiFluidGMCoupling)then
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "vol", "z0x", "z0y", "bmin", "rho", "p"'
    else
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "vol", "z0x", "z0y", "bmin",&           
            &"rho","p","Hprho","Oprho","Hpp","Opp"'
    end if
    write(UNITTMP_,'(a,i4,a,i4,a)') &
         'ZONE T="SAVE", I=',jsize,', J=',isize,', K=1, F=POINT'
    do i=1,iSizeIn; do j=1,jSizeIn
       if(.not. DoMultiFluidGMCoupling)then
          write(UNITTMP_,'(2i4,6G14.6)') j,i, &
               Buffer_IIV(i,j,vol_),Buffer_IIV(i,j,z0x_),Buffer_IIV(i,j,z0y_), &
               Buffer_IIV(i,j,bmin_),Buffer_IIV(i,j,rho_),Buffer_IIV(i,j,p_)
       else
          !multi-fluid                                                                
          write(UNITTMP_,'(2i4,8G14.6)') j,i, &
               Buffer_IIV(i,j,vol_),Buffer_IIV(i,j,z0x_),Buffer_IIV(i,j,z0y_), &
               Buffer_IIV(i,j,bmin_),Buffer_IIV(i,j,rho_),Buffer_IIV(i,j,p_), &
               Buffer_IIV(i,j,Hprho_),Buffer_IIV(i,j,Oprho_), &
               Buffer_IIV(i,j,Hpp_),Buffer_IIV(i,j,Opp_)
       end if
       end do; end do
    CLOSE(UNITTMP_)
  end subroutine write_data

end subroutine IM_put_from_gm

!==============================================================================
subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
  ! Puts satellite locations and names from GM into IM variables.
  !!!DTW 2007
  use RCM_variables, ONLY: nImSats, DoWriteSats, NameSat_I, SatLoc_3I
  use ModNumConst,   ONLY: cDegToRad
  
  implicit none
  character (len=*),parameter :: NameSub='IM_put_sat_from_gm'

  ! Arguments
  integer, intent(in)            :: nSats
  real, intent(in)               :: Buffer_III(3,2,nSats)
  character(len=100), intent(in) :: Buffer_I(nSats)

  ! Internal variables
  integer :: iError, iSat, l1, l2
  !--------------------------------------------------------------------------- 
  ! Activate satellite writing in RCM
  DoWriteSats = .true.
  nImSats = nSats

  ! Check allocation of sat tracing variables
  if(allocated(SatLoc_3I)) deallocate(SatLoc_3I)
  if(allocated(NameSat_I)) deallocate(NameSat_I)

  allocate(SatLoc_3I(3,2,nImSats), stat=iError)
  allocate(NameSat_I(nImSats),     stat=iError)

  ! Assign incoming values, remove path and extension from name.
  SatLoc_3I = Buffer_III
  do iSat=1, nSats
     l1 = index(Buffer_I(iSat), '/', back=.true.) + 1
     l2 = index(Buffer_I(iSat), '.') - 1
     if (l1-1<=0) l1=1
     if (l2+1<=0) l2=len_trim(Buffer_I(iSat))
     NameSat_I(iSat) = Buffer_I(iSat)(l1:l2)
  end do

  ! Change to correct units (degrees to radians)
  SatLoc_3I(1,2,:) = (90. - SatLoc_3I(1,2,:)) * cDegToRad
  SatLoc_3I(2,2,:) =        SatLoc_3I(2,2,:)  * cDegToRad

end subroutine IM_put_sat_from_gm

!==============================================================================
subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  use CON_time, ONLY : get_time
  use RCM_variables
  use ModNumConst, ONLY: cRadToDeg
  implicit none
  character (len=*),parameter :: NameSub='IM_get_for_gm'

  integer, intent(in)                                :: iSizeIn,jSizeIn,nVar
  real, dimension(iSizeIn,jSizeIn,nVar), intent(out) :: Buffer_IIV
  character (len=*),intent(in)                       :: NameVar

  !LOCAL VARIABLES:
  real :: tSimulation
  integer :: iTimeStart
  integer, parameter :: pres_=1, dens_=2, Hpres_=3,Opres_=4,Hdens_=5,Odens_=6

  integer :: i,j,k
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  if(nVar > 2)DoMultiFluidGMCoupling = .true.

  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if (DoTestMe) &
       write(*,*)NameSub,' starting with iSizeIn,jSizeIn,nVar,NameVar=',&
       iSizeIn,jSizeIn,nVar,NameVar

  if(DoMultiFluidGMCoupling)then
     if(NameVar /= 'p:rho:Hpp:Opp:Hprho:Oprho') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  else
     if(NameVar /= 'p:rho') &
          call CON_stop(NameSub//' invalid NameVar='//NameVar)
  endif

  if(IsUninitialized)then
     if(DoTestMe)write(*,*) NameSub,' call RCM_advec(1...)'
     if(.not.DoneGmCoupling)call CON_stop(NameSub//&
          ' SWMF_ERROR: IM/RCM has not been coupled with GM')
     call get_time(tSimulationOut = tSimulation)
     iTimeStart=nint(tSimulation)

     call RCM_advec (1, iTimeStart, 99999, 0)
     IsUninitialized = .false.
     if(DoTestMe)write(*,*) NameSub,' done RCM_advec(1...)'
  end if

  if(iSizeIn /= iSize .or. jSizeIn /= jSize)then
     write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  Buffer_IIV = 0.

  !Fill pressure and density
  do i=1,iSize; do j=1,jSize
     if( i<imin_j(j) .or. vm(i,j) <= 0.0 ) then
        Buffer_IIV(i,j,pres_) = -1.
        Buffer_IIV(i,j,dens_) = -1.
     else
        do k=1,kcsize
           Buffer_IIV(i,j,pres_) = Buffer_IIV(i,j,pres_) + &
                vm(i,j)**2.5*eeta(i,j,k)*ABS(alamc(k))
           Buffer_IIV(i,j,dens_) = Buffer_IIV(i,j,dens_) + &
                eeta(i,j,k)*vm(i,j)**1.5 * xmass(ikflavc(k))
        end do
     end if
     if(DoMultiFluidGMCoupling)then
        !  Multifluid case                                                         
        if( i<imin_j(j) .or. vm(i,j) <= 0.0 ) then
           Buffer_IIV(i,j,Hpres_) = -1.
           Buffer_IIV(i,j,Opres_) = -1.
           Buffer_IIV(i,j,Hdens_) = -1.
           Buffer_IIV(i,j,Odens_) = -1.
        else
           do k=kmin(2),kmax(2)
              Buffer_IIV(i,j,Hpres_) = Buffer_IIV(i,j,Hpres_) + &
                   vm(i,j)**2.5*eeta(i,j,k)*ABS(alamc(k))
              Buffer_IIV(i,j,Hdens_) = Buffer_IIV(i,j,Hdens_) + &
                   eeta(i,j,k)*vm(i,j)**1.5 * xmass(ikflavc(k))
           end do
           do k=kmin(3),kmax(3)
              Buffer_IIV(i,j,Opres_) = Buffer_IIV(i,j,Opres_) + &
                   vm(i,j)**2.5*eeta(i,j,k)*ABS(alamc(k))
              Buffer_IIV(i,j,Odens_) = Buffer_IIV(i,j,Odens_) + &
                   eeta(i,j,k)*vm(i,j)**1.5 * xmass(ikflavc(k))
           end do
        end if
     end if

     ! Only a not-a-number can be less than zero and larger than one
     if(  .not. Buffer_IIV(i,j,pres_) > 0 .and. &
          .not. Buffer_IIV(i,j,pres_) < 1) then
        write(*,*)NameSub,': ERROR IN PRESSURE'
        write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,pres_)
        write(*,*)NameSub,': Lon,Lat[dg]=', &
             aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
        write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
        write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
        write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
        write(*,*)NameSub,': alamc      =',alamc
        call CON_stop(NameSub // ' ERROR: Not a number found in IM pressure !')
     end if
     if(  .not. Buffer_IIV(i,j,dens_) > 0 .and. &
          .not. Buffer_IIV(i,j,dens_) < 1) then
        write(*,*)NameSub,': ERROR IN DENSITY'
        write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,dens_)
        write(*,*)NameSub,': Lon,Lat[dg]=', &
             aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
        write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
        write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
        write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
        write(*,*)NameSub,': ikflavc(k) =',ikflavc(k)
        write(*,*)NameSub,': xmass(ikflavc(k))=',xmass(ikflavc(k))
        call CON_stop(NameSub // ' ERROR: Not a number found in IM density !')
     end if

     !multi-fluid
     if(DoMultiFluidGMCoupling)then
        if(  .not. Buffer_IIV(i,j,Hpres_) > 0 .and. &
             .not. Buffer_IIV(i,j,Hpres_) < 1) then
           write(*,*)NameSub,': ERROR IN PRESSURE'
           write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Hpres_)
           write(*,*)NameSub,': Lon,Lat[dg]=', &
                aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
           write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
           write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
           write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
           write(*,*)NameSub,': alamc      =',alamc
           call CON_stop(NameSub // ' ERROR: Not a number found in IM Hp pressure !')
        end if
        if(  .not. Buffer_IIV(i,j,Hdens_) > 0 .and. &
             .not. Buffer_IIV(i,j,Hdens_) < 1) then
           write(*,*)NameSub,': ERROR IN DENSITY'
           write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Hdens_)
           write(*,*)NameSub,': Lon,Lat[dg]=', &
                aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
           write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
           write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
           write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
           write(*,*)NameSub,': ikflavc(k) =',ikflavc(k)
           write(*,*)NameSub,': xmass(ikflavc(k))=',xmass(ikflavc(k))
           call CON_stop(NameSub // ' ERROR: Not a number found in IM Hp density !')
        end if
        if(  .not. Buffer_IIV(i,j,Opres_) > 0 .and. &
             .not. Buffer_IIV(i,j,Opres_) < 1) then
           write(*,*)NameSub,': ERROR IN PRESSURE'
           write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Opres_)
           write(*,*)NameSub,': Lon,Lat[dg]=', &
                aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
           write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
           write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
           write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
           write(*,*)NameSub,': alamc      =',alamc
           call CON_stop(NameSub // ' ERROR: Not a number found in IM Op pressure !')
        end if
        if(  .not. Buffer_IIV(i,j,Odens_) > 0 .and. &
             .not. Buffer_IIV(i,j,Odens_) < 1) then
           write(*,*)NameSub,': ERROR IN DENSITY'
           write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Odens_)
           write(*,*)NameSub,': Lon,Lat[dg]=', &
                aloct(i,j)*cRadToDeg,90.0-colat(i,j)*cRadToDeg
           write(*,*)NameSub,': imin_j(j)  =',imin_j(j)
           write(*,*)NameSub,': vm(i,j)    =',vm(i,j)
           write(*,*)NameSub,': eeta(i,j,:)=',eeta(i,j,:)
           write(*,*)NameSub,': ikflavc(k) =',ikflavc(k)
           write(*,*)NameSub,': xmass(ikflavc(k))=',xmass(ikflavc(k))
           call CON_stop(NameSub // ' ERROR: Not a number found in IM Op density !')
        end if
     end if
  end do; end do

  where(Buffer_IIV(:,:,pres_) > 0.0) &
       Buffer_IIV(:,:,pres_) = Buffer_IIV(:,:,pres_) * 1.67E-35

  ! Units of rcm_mass_density are kg/m3
  where(Buffer_IIV(:,:,dens_) > 0.0) &
       Buffer_IIV(:,:,dens_) = Buffer_IIV(:,:,dens_) / 6.37E+15

  if(DoMultiFluidGMCoupling)then
     ! MultiFluid                                                                  
     where(Buffer_IIV(:,:,Hpres_) > 0.0) &
          Buffer_IIV(:,:,Hpres_) = Buffer_IIV(:,:,Hpres_) * 1.67E-35
     where(Buffer_IIV(:,:,Opres_) > 0.0) &
          Buffer_IIV(:,:,Opres_) = Buffer_IIV(:,:,Opres_) * 1.67E-35
     ! Units of rcm_mass_density are kg/m3                                                
     where(Buffer_IIV(:,:,Hdens_) > 0.0) &
          Buffer_IIV(:,:,Hdens_) = Buffer_IIV(:,:,Hdens_) / 6.37E+15
     where(Buffer_IIV(:,:,Odens_) > 0.0) &
          Buffer_IIV(:,:,Odens_) = Buffer_IIV(:,:,Odens_) / 6.37E+15
  endif

  if(DoTestMe)write(*,*) NameSub,' finished'

end subroutine IM_get_for_gm

!==============================================================================

subroutine IM_init_session(iSession, TimeSimulation)

  use RCM_variables, ONLY: iUnitOut

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_init_session'

  logical :: IsUninitialized = .true.
  !---------------------------------------------------------------------------

  ! IM will be initialized after being coupled to GM
  RETURN

end subroutine IM_init_session

!==============================================================================

subroutine IM_finalize(TimeSimulation)

  use RCM_variables, ONLY: iUnitOut

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  !LOCAL VARIABLES:
  character(len=*), parameter :: NameSub='IM_finalize'
  !---------------------------------------------------------------------------

  call IM_write_prefix; write(iUnitOut,*) &
       NameSub,' at TimeSimulation=',TimeSimulation

  call RCM_advec (3, 0, 0, 0)

end subroutine IM_finalize

!==============================================================================

subroutine IM_save_restart(TimeSimulation)

  use RCM_variables, ONLY: iUnitOut

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_save_restart'

  call IM_write_prefix; write(iUnitOut,*) &
       NameSub,' at TimeSimulation=',TimeSimulation

  call RCM_advec (4, nint(TimeSimulation), nint(TimeSimulation), 0)

end subroutine IM_save_restart

!BOP ==========================================================================
!ROUTINE: IM_run - run IM
!INTERFACE:
subroutine IM_run(TimeSimulation,TimeSimulationLimit)

  !USES:
  use RCM_variables, ONLY: iDtRcm, IsUninitialized, DoneGmCoupling
  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  !LOCAL VARIABLES:
  integer :: iDtNow, iTimeEnd, iTimeStart
  character(len=*), parameter :: NameSub='IM_run'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)

  iDtNow = ceiling(min(real(iDtRcm), &
       max(1.0,TimeSimulationLimit - TimeSimulation)))

  if(DoTest)write(*,*)NameSub,' TimeSim,Limit,iDtNow=',&
       TimeSimulation,TimeSimulationLimit,iDtNow

  iTimeStart = nint(TimeSimulation)
  iTimeEnd   = iTimeStart + iDtNow ! we only want to make a single step !!!
  call IM_write_prefix; write(*,*)'iTimeStart, iTimeEnd, iDtNow', &
       iTimeStart, iTimeEnd, iDtNow

  if(IsUninitialized)then
     if(.not.DoneGmCoupling)call CON_stop(NameSub//&
          ' SWMF_ERROR: IM/RCM has not been coupled with GM')
     call RCM_advec (1, iTimeStart, 99999, 0)
     IsUninitialized = .false.
  end if

  call RCM_advec (2, iTimeStart, iTimeEnd, iDtNow)

  ! return time at the end of the time step to CON
  TimeSimulation   = TimeSimulation + iDtNow 


end subroutine IM_run

!===========================================================================

subroutine IM_write_prefix

  use RCM_variables, ONLY: iUnitOut, STDOUT_, StringPrefix

  implicit none

  if(iUnitOut==STDOUT_)write(*,'(a)',ADVANCE='NO')trim(StringPrefix)

end subroutine IM_write_prefix

!===========================================================================
