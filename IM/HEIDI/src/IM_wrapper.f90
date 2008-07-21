! Wrapper for Internal Magnetosphere (IM) component
!=============================================================================

subroutine IM_set_param(CompInfo,TypeAction)
  
  use CON_comp_info
  use ModProcIM
  use ModHeidiMain
  use ModReadParam, only: i_session_read
  use ModUtilities, ONLY: fix_dir_name, check_dir, lower_case
   
  implicit none
  character (len=*), parameter :: NameSub='IM_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  
  !LOCAL VARIABLES:
  character (len=100) :: NameCommand, StringPlot
  logical             :: DoEcho=.false.
  logical             :: UseStrict=.true.  
  
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                         &
          Use=.true.,                           &
          NameVersion='HEIDI (Liemohn)', &
          Version=1.1)
  case('MPI')
!     call MPI_INIT(iError)
!     iComm= MPI_COMM_WORLD
!     call MPI_COMM_RANK(iComm, iProc, iError)
!     call MPI_COMM_SIZE(iComm, nProc, iError)   
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     if(nProc>4)call CON_stop(&
          'IM_init_mpi: IM_ERROR this version can run on 4 PE !')
  case('READ')
     call heidi_read
  case('CHECK')
     !We should check and correct parameters here
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',i_session_read()
     RETURN
     call heidi_check
  case('STDOUT')
  case('FILEOUT')
  case('GRID')
     call IM_set_grid

  case default
     call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)

  end select
end subroutine IM_set_param

!============================================================================
subroutine IM_set_grid
  use ModIonosphere
  use ModHeidiSize
  use ModHeidiMain
  use ModProcIM
  use ModNumConst
  use CON_coupler, ONLY: set_grid_descriptor, is_proc, IM_
!  use heidi_read
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

  Radius_I(1) = 110*1000.0 ! radial size of the ionosphere in meters

  ! Total hack, because I don't actually know the real grid we will use in the 
  ! actual code.

  call set_grid_descriptor(                        &
       IM_,                          &! component index
       nDim=2,                       &! dimensionality
       nRootBlock_D=(/1,1/),         &! north+south hemispheres
       nCell_D =(/IONO_nTheta - 1,IONO_nPsi - 1/), &! size of node based grid
       XyzMin_D=(/cOne, cOne/),      &! min colat and longitude indexes
       XyzMax_D=(/real(IONO_nTheta-1),&
       real(IONO_nPsi)/),            &! max colat and longitude indexes
       TypeCoord='SMG',                            &! solar magnetic coord.
       Coord1_I=IONO_NORTH_Theta(:,1),             &! colatitudes
       Coord2_I=IONO_NORTH_Psi(1,:),               &! longitudes
       Coord3_I=(/Radius_I/),                          &! radial size in meters
       IsPeriodic_D=(/.false.,.true./))

!!!  ! IM grid size in generalized coordinates
!!!  call set_grid_descriptor( IM_,                         & ! component index
!!!       nDim=4,                                           & ! dimensionality
!!!       nRootBlock_D=(/1,1/),                             & ! single block
!!!       nCell_D=(/nT,nR/),                          & ! size of cell based grid
!!!       XyzMin_D=(/cHalf, cHalf/),                        & ! min gen.coords for cells
!!!       XyzMax_D=(/nT+cHalf,nR+cHalf/),             & ! max gen.coords for cells
!!!       TypeCoord='SMG',                                   & ! solar magnetic coord
!!!       Coord1_I=phi,  & ! magnetic local times
!!!       Coord2_I=z,    & ! l-shell?
!!!!       Coord1_I=real(colat(1:NT,1)),                     & ! colatitudes
!!!!       Coord2_I=real(aloct(1,1:NR)),                     & ! longitudes
!!!!       Coord3_I=Radius_I,                                & ! radial size in meters
!!!!       Coord4_I=NPA(0:90)                                & ! Grid in pitch angle
!!!!       COORD5_I=NE                                       & ! Grid in energy
!!!       IsPeriodic_D=(/.true.,.false./))    ! periodic in longitude
!!!!       IsPeriodic_D=(/.true.,.false.,.false.,.false./))    ! periodic in longitude

end subroutine IM_set_grid
!==============================================================================
subroutine IM_print_variables(NameSource)

  implicit none

  character(len=*),intent(in) :: NameSource


  select case(NameSource)
  case('IE')
     write(*,*) 'im_print_variables called from IE'
  case('GM')
     write(*,*) 'im_print_variables called from GM'
  case default
     write(*,*) 'im_print_variables: incorrect NameSource=',NameSource
     RETURN
  end select

  write(*,*) "this routine doesn't do anything!"

end subroutine IM_print_variables

!!!  integer            :: nFile=0
!!!  character(len=100) :: NameFile
!!!  character(len=100) :: NameVar
!!!  integer            :: i,j
!!!  real               :: Lat,Lon
!!!  !--------------------------------------------------------------------------
!!!  select case(NameSource)
!!!  case('IE')
!!!     NameVar='j i lon lat jr pot sigmaH sigmaP'
!!!  case('GM')
!!!     NameVar='j i lon lat density pressure vm xmin ymin bmin temperature'
!!!  case default
!!!     write(*,*)NameSub,': incorrect NameSource=',NameSource
!!!     RETURN
!!!  end select
!!!
!!!  nFile=nFile+1
!!!  write(NameFile,'(a,i1,a)')'IM_from_'//NameSource//'_',nFile,'.dat'
!!!  open(UNITTMP_,file=NameFile)
!!!  write(UNITTMP_,'(a)')trim(NameVar)
!!!
!!!  do i=1,iSize
!!!     do j=1,jSize
!!!        Lon = (        aloct(i,j))*(180./cPi)
!!!        Lat = (cHalfPi-colat(i,j))*(180./cPi)
!!!        select case(NameSource)
!!!        case('IE')
!!!           write(UNITTMP_,'(2i4,6G14.6)')j,i,Lon,Lat,v(i,j),birk_mhd(i,j),&
!!!                sigmaH_mhd(i,j),sigmaP_mhd(i,j)
!!!        case('GM')
!!!           write(UNITTMP_,'(2i4,9G14.6)')j,i,Lon,Lat,density(i,j),pressure(i,j),&
!!!                vm(i,j),xmin(i,j),ymin(i,j),bmin(i,j),temperature(i,j)
!!!        end select
!!!     end do
!!!  end do
!!!  close(UNITTMP_)
!!!
!!!end subroutine IM_print_variables

!==============================================================================
subroutine IM_get_for_ie(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  write(*,*) 'This is not working'  

end subroutine IM_get_for_ie

!==============================================================================
subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

  use CON_coupler,   ONLY: IndexPtrType, WeightPtrType
  use ModIonosphere, ONLY: IONO_NORTH_PHI, IONO_nTheta, IONO_nPsi

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

  if(i<1.or.i>IONO_nTheta.or.j<1.or.j>IONO_nPsi)then
     write(*,*)'i,j,DoAdd=',i,j,DoAdd
     call CON_stop('IM_put_from_ie: index out of range')
  end if

  if(DoAdd)then
     IONO_NORTH_PHI(i,j)        = IONO_NORTH_PHI(i,j)        + Buff_V(1)
  else
     IONO_NORTH_PHI(i,j)        = Buff_V(1)
  end if

end subroutine IM_put_from_ie
!==============================================================================
subroutine IM_put_from_ie_complete

  write(*,*) "Don't know what this is really supposed to do.  I think that it is"
  write(*,*) "Supposed to be applying periodic boundaries...?"

end subroutine IM_put_from_ie_complete
!==============================================================================
subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

write(*,*) 'This is not working' 

end subroutine IM_put_from_gm

!==============================================================================
subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
  ! Puts satellite locations and names from GM into IM variables.
  !!!DTW 2007

write(*,*) 'This is not working' 

end subroutine IM_put_sat_from_gm

!==============================================================================
subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

write(*,*) 'This is not working'  

end subroutine IM_get_for_gm

!==============================================================================

 subroutine IM_init_session(iSession, TimeSimulation)
    implicit none
    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession       ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    logical :: IsUninitialized = .true.
    if(IsUninitialized)then
       call heidi_init

       IsUninitialized = .false.
    end if
  end subroutine IM_init_session

!==============================================================================
 subroutine IM_finalize(TimeSimulation)
    use ModProcIM
    use ModInit,ONLY:nS
  
    implicit none

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time
    integer:: iSpecies
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do iSpecies=1,NS
       CLOSE(15+iSpecies)          ! Closes continuous output file
    end do

    CLOSE(13)	            ! Closes sw1 input file
    CLOSE(15)		    ! Closes sw2 input file
    CLOSE(14)               ! Closes MPA input file
    CLOSE(16)               ! Closes SOPA input file
    CLOSE(18)               ! Closes FPOT input file
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
    call MPI_finalize(iError)
  end subroutine IM_finalize

! =============================================================================
 subroutine IM_run(TimeSimulation,TimeSimulationLimit)
    implicit none

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout) :: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
    call heidi_run 
  end subroutine IM_run

!===========================================================================

  subroutine IM_save_restart(TimeSimulation)
    implicit none

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IM_save_restart'

  end subroutine IM_save_restart


!===========================================================================


