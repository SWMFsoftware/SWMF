!BOP
!MODULE: GM_put_from_ie_swmf - transform IE potential to GM boundary velocities
!IROUTINE:GMIE_set_grid - set 2D spherical grid in GM_ to get IE_ potential
!INTERFACE:
subroutine GMIE_set_grid
  !USES:
  use CON_coupler
  use ModInnerBC
  use ModMain,ONLY:NameThisComp
  use ModProcMH,ONLY:iProc
!EOP
  implicit none
  integer,dimension(2)::nCell_D
  integer,parameter::GmIeGrid_=MaxComp+1
  integer::iError
  real::rBoundary=cOne
  real,allocatable,dimension(:)::Colat_I,Longit_I
  logical::DoTest,DoTestMe

  call CON_set_do_test('test_grids',DoTest,DoTestMe)

  if(.not.(done_dd_init(GM_).and.done_dd_init(IE_)))&
       call CON_stop('GMIE_grid should be intialized after GM_ and IE_')

  nCell_D=ncells_decomposition_D(IE_)
  if(NameThisComp=='GM') then
     call set_mapping_param
     call init_inner_bc_arrays(nCell_D(1)+1,nCell_D(2)+1)
  end if

  allocate(Colat_I(1:2*nCell_D(1)+1),&
       Longit_I(1:nCell_D(2)+1),stat=iError)
  call check_allocate(iError,'GMIE_set_grid')
  Colat_I=cZero;Longit_I=0



  call init_decomposition(GmIeGrid_,GM_,2)
  if(NameThisComp=='GM')then
     Colat_I(1:IONO_nTheta)=IONO_NORTH_Theta(:,1)
     Colat_I(IONO_nTheta:2*IONO_nTheta-1)=IONO_SOUTH_Theta(:,1)
     Longit_I=IONO_NORTH_Psi(1,:)
     rBoundary=IONO_Radius_Mag_Boundary/IONO_Radius
  end if
  call set_coord_system(GmIeGrid_, &
       Grid_C(GM_)%TypeCoord,&
       Colat_I,&
       Longit_I,&
       (/rBoundary/))
  if(NameThisComp=='GM'.and.iProc==0)call get_root_decomposition(&
       GridID_=GmIeGrid_,&
       iRootMapDim_D=(/2,1/),&
       XyzMin_D=(/cOne,cOne/),&
       XyzMax_D=(/real(2*IONO_nTheta-1),real(IONO_nPsi)/),&
       nCells_D=nCell_D,&
       iBlock_I=(/1,2/))
  call bcast_decomposition(GmIeGrid_)
  if(DoTest)call test_global_message_pass(GmIeGrid_)
end subroutine GMIE_set_grid

!IROUTINE: GM_put_from_ie_swmf - put the interpolated IE potential to GmIe grid
!INTERFACE:
subroutine GM_put_from_ie_swmf(nPartial,iPut,Put,W,DoAdd,State_V,nVar)
  !USES:
  use CON_router
  use ModInnerBC,ONLY:IONO_NORTH_Phi_BC,IONO_SOUTH_Phi_BC
  implicit none
!INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iPut,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::State_V
  character(LEN=*),parameter::SubName='GM_put_from_ie_swmf'
!REVISION HISTORY:
!03SEPT03     I.Sokolov <igorsok@umich.edu> - prototype/code/prolog
!EOP

  if(DoAdd)call CON_stop(SubName//' DoAdd can not be .true.')
  select case(Put%iCB_II(3,iPut))
  case(1)
     IONO_NORTH_Phi_BC(Put%iCB_II(1,iPut),Put%iCB_II(2,iPut))=State_V(1)
  case(2)
     IONO_SOUTH_Phi_BC(Put%iCB_II(1,iPut),Put%iCB_II(2,iPut))=State_V(1)
  case default
     call CON_stop(SubName//'  - wrong block number:',Put%iCB_II(3,iPut))
  end select
end subroutine GM_put_from_ie_swmf
!BOP
!IROUTINE:transform_phi_bc_to_u_bc - transform  GmIe potential  to GM boundary  velocities
!INTERFACE: 
subroutine transform_phi_bc_to_u_bc
  !USES:
  use ModProcMH
  use ModInnerBC
  use CON_coupler
!EOP
  implicit none

  integer,parameter::GmIeGrid_=MaxComp+1
  integer::nCell_D(2),nSize,iError

  if(iProc<0)return
  nCell_D=ncells_decomposition_d(GmIeGrid_)
  nSize=product(nCell_D+1)
  if(nProc>1)then
     call MPI_BCAST(IONO_NORTH_Phi_BC(1,1),nSize,MPI_REAL,0,iComm,iError)
     call MPI_BCAST(IONO_SOUTH_Phi_BC(1,1),nSize,MPI_REAL,0,iComm,iError)
  end if
  call ionosphere_magBCs(&
     IONO_NORTH_PHI_BC,&
     IONO_NORTH_ETh_BC,&
     IONO_NORTH_EPs_BC,&
     IONO_NORTH_UR_BC, &
     IONO_NORTH_UTh_BC,&
     IONO_NORTH_UPs_BC,&
     IONO_Radius_Mag_Boundary,  &
     IONO_NORTH_Theta,      &
     IONO_NORTH_Psi,        &
     IONO_nTheta,              &
     IONO_nPsi)
  call ionosphere_magBCs(&
     IONO_SOUTH_PHI_BC,&
     IONO_SOUTH_ETh_BC,&
     IONO_SOUTH_EPs_BC,&
     IONO_SOUTH_UR_BC, &
     IONO_SOUTH_UTh_BC,&
     IONO_SOUTH_UPs_BC,&
     IONO_Radius_Mag_Boundary,  &
     IONO_SOUTH_Theta,      &
     IONO_SOUTH_Psi,        &
     IONO_nTheta,              &
     IONO_nPsi)
end subroutine transform_phi_bc_to_u_bc
