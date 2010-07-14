!^CFG COPYRIGHT UM
!^CMP FILE SP
!\
! This coupler employs the following global arrays.
! SP_Xyz_DI - is the array of the Lagrangian points.
! A part of this array, as well as the mask 'SP_IsInIH' is only availbale at the
! PE set, at which either IH or SP  run. This coordinates are expressed in terms of
! the length units and with respect to the frame of reference defined in IH.
! Another part of this array, as well as the mask 'SP_IsInSC' is only availbale at the
! PE set, at which either SC or SP  run. This coordinates are expressed in terms of
! the length units and with respect to the frame of reference defined in SC.
! SP_XyzSP - in the array of all Lagrangian points, in units and in the frame of 
! reference defined at SP, is available only at the PE set, at which SP runs.
!/
Module CON_couple_mh_sp
  use CON_coupler
  use CON_global_message_pass
  use CON_axes
  implicit none
  private !Except
  public::couple_mh_sp_init
  public::couple_ih_sp              !^CMP IF IH
  public::couple_sc_sp              !^CMP IF SC

  type(GridDescriptorType),save::SP_GridDescriptor !Target

  type(GridDescriptorType),save::IH_GridDescriptor !Source  !^CMP IF IH
  type(RouterType),save,private::RouterIhSp                 !^CMP IF IH

  type(GridDescriptorType),save::SC_GridDescriptor !Source  !^CMP IF SC
  type(RouterType),save,private::RouterScSp                 !^CMP IF SC

  logical,save::DoInit=.true.
  real,allocatable,dimension(:,:)::XyzTemp_DI

  real,dimension(:,:),pointer ::Xyz_DI
  logical,dimension(:),pointer :: Is_I
  integer,parameter::nPointMax=5000
  integer::nPoint=0
  integer::iPoint
  integer::iError
  real::bDxyz_I(1:6)!The interpolated values of full B and DXyz
  real::DsResolution,XyzLine_D(3)
  real,save::rBoundIh                !^CMP IF IH
  real,save::rBoundSc                !^CMP IF SC
  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real,dimension(3,3)::ScToIh_DD,ScToSp_DD,IhToSp_DD
  real :: tNow

contains
  !========================================================================
  subroutine xyz_point_mapping(&
       SP_nDim,SP_Xyz_D,MH_nDim,MH_Xyz_D,IsInterfacePoint)
    use ModConst
    integer,intent(in)::SP_nDim,MH_nDim
    real,dimension(SP_nDim),intent(in)::SP_Xyz_D
    real,dimension(MH_nDim),intent(out)::MH_Xyz_D
    logical,intent(out)::IsInterfacePoint
    MH_Xyz_D=XyzTemp_DI(:,iPoint)
    IsInterfacePoint=.true.
  end subroutine xyz_point_mapping
  !==================================================================
  subroutine SP_put_a_line_point(nPartial,&
       iPutStart,&
       Put,&
       w,&
       DoAdd,&
       Buff_I,nVar)
    implicit none
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::w
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer::iCell
    real:: Weight
    !-------------------------------------------------------------------------
    Weight=w%Weight_I(iPutStart)
    if(DoAdd)then
       bDxyz_I = bDxyz_I + Buff_I*Weight
    else
       bDxyz_I = Buff_I*Weight
    end if
  end subroutine SP_put_a_line_point
  !==================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst
    interface
       subroutine IH_get_a_line_point(nPartial,&         !^CMP IF IH BEGIN
            iGetStart,&
            Get,&
            Weight,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::Weight
         real,dimension(nVar),intent(out)::Buff_I
       end subroutine IH_get_a_line_point                !^CMP END IH
       subroutine SC_get_a_line_point(nPartial,&         !^CMP IF SC BEGIN
            iGetStart,&
            Get,&
            Weight,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::Weight
         real,dimension(nVar),intent(out)::Buff_I
       end subroutine SC_get_a_line_point               !^CMP END SC
    end interface

    logical::DoneRestart
    integer::nU_I(2)
    !======================================================================
    if(.not.DoInit)return
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    DoInit=.false.
    !The initialization can be done only once
    

    call get_time(tSimulationOut=tNow)


    !Check, if the lagrangian meshes are saved in restart files


    if(is_proc(SP_))then
       call read_global_vector('SP_Xyz_DI')
       DoneRestart=used_vector('SP_Xyz_DI')
       call SP_get_line_param(&
            DsResolution,XyzLine_D&
            ,rBoundSc &              !^CMP IF SC
            ,rBoundIh &              !^CMP IF IH
            ) 
    end if
    call MPI_BCAST(DoneRestart,1,MPI_LOGICAL,i_proc0(SP_),i_comm(),iError)
    if(use_comp(IH_))call MPI_bcast(&                     !^CMP IF IH
         rBoundIh,1,MPI_REAL,i_proc0(SP_),i_comm(),iError)!^CMP IF IH

    if(use_comp(SC_))call MPI_bcast(&                     !^CMP IF SC
         rBoundSc,1,MPI_REAL,i_proc0(SP_),i_comm(),iError)!^CMP IF SC

    if(.not.DoneRestart)then
       !Construct auxiliary SP_grid with ONE point per line
       if(is_proc0(SP_))&
            call get_root_decomposition(&
            SP_,&                    !GridDescriptor to be constructed
            iRootMapDim_D=(/1/),&    !The block amount, along each direction(D)
            XyzMin_D=(/cHalf/),&     !Minimal gen. coordinates, along each D 
            XyzMax_D=(/cHalf+cOne/),&!Maximal gen. coordinates, along each D
            nCells_D=(/1/),&
            PE_I=(/0/))              !PEs layout, throughout all the blocks
       call bcast_decomposition(SP_)
       

       call set_standard_grid_descriptor(SP_,GridDescriptor=&
            SP_GridDescriptor)

       !Set the initial point for the line

       if(is_proc0(SP_))then
          if(sum(XyzLine_D**2)<cOne)then
             !Calculate the Earth position in SGI
             XyzLine_D = XyzPlanetHgi_D   !In HGI, in meters
             if(use_comp(IH_))then             !^CMP IF IH BEGIN
                XyzLine_D=matmul(&
                     transform_matrix(tNow,'HGI',&
                     Grid_C(IH_) % TypeCoord),XyzLine_D)/Grid_C(IH_)%UnitX
             else                              !^CMP END IH
                XyzLine_D=matmul(&
                     transform_matrix(tNow,'HGI',&
                     Grid_C(SC_) % TypeCoord),XyzLine_D)/Grid_C(SC_)%UnitX
             end if                            !^CMP IF IH
             write(*,*)'Actual Earth position is at', XyzLine_D
          end if
       end if

       allocate(XyzTemp_DI(3,nPointMax),stat=iError)
       call check_allocate(iError,NameSub//': XyzTemp_DI')

   
       if(is_proc(SP_))then
          iPoint=1
          XyzTemp_DI(:,iPoint)=XyzLine_D(:)
       end if
       
       if(use_comp(IH_))then       !^CMP IF IH BEGIN
          call IH_synchronize_refinement(i_proc0(IH_),i_comm())
  
          call trace_line(IH_,&
               IH_GridDescriptor,RouterIhSp,rBoundIh**2,IH_get_a_line_point)
          call MPI_bcast(nPoint,1,MPI_INTEGER,&
               i_proc0(SP_),i_comm(),iError)

          if(DoTest.and.is_proc0())write(*,*)'nPointIH=',nPoint
       end if                      !^CMP END IH

       if(use_comp(SC_))then       !^CMP IF SC BEGIN
          call SC_synchronize_refinement(i_proc0(SC_),i_comm())
          call MPI_bcast(rBoundSc,1,MPI_REAL,&
               i_proc0(SP_),i_comm(),iError)
          if(Grid_C(SC_)%TypeCoord/=Grid_C(IH_)%TypeCoord& !^CMP IF IH BEGIN
               .and.is_proc0(SP_).and.use_comp(IH_))&
               XyzTemp_DI(:,iPoint)=matmul(transform_matrix(&
               tNow,Grid_C(IH_)%TypeCoord,Grid_C(SC_)%TypeCoord),&
               XyzTemp_DI(:,iPoint))                       !^CMP END IH

          call trace_line(SC_,&
               SC_GridDescriptor,RouterScSp,rBoundSc**2,SC_get_a_line_point)
          call MPI_bcast(nPoint,1,MPI_INTEGER,&
               i_proc0(SP_),i_comm(),iError)
          if(DoTest.and.is_proc0())write(*,*)'nPoint',nPoint
       end if                      !^CMP END SC
    
       !\
       !Reset SP_domain_decomposition
       !/
       if(is_proc0(SP_))&
            call get_root_decomposition(&
            SP_,&                            !GridDescroptor to be constructed
            iRootMapDim_D=(/1/),&            !The block amount, along each direction(D)
            XyzMin_D=(/cHalf/),&             !Minimal gen. coordinates, along each D 
            XyzMax_D=(/cHalf+real(nPoint)/),&!Maximal gen. coordinates, along each D
            nCells_D=(/nPoint/))
       call bcast_decomposition(SP_)
       call clean_grid_descriptor(SP_GridDescriptor)
       call set_standard_grid_descriptor(SP_,GridDescriptor=&
            SP_GridDescriptor)


       !\
       !Save initial positions of all the grid points
       !/
       call allocate_vector('SP_Xyz_DI',3,SP_GridDescriptor)
       if(is_proc0())&
            write(*,*)'Allocated global vector SP_Xyz_DI, dims=',&
            ubound_vector('SP_Xyz_DI')
       if(is_proc0(SP_))then
          call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
          do iPoint=1,nPoint
             Xyz_DI(:,iPoint)= XyzTemp_DI(:,nPoint-iPoint+1)
          end do
          nullify(Xyz_DI)
          if(DoTest)call save_global_vector('SP_Xyz_DI',iFileIn=0)
       end if
       deallocate(XyzTemp_DI)
    else
       !In this case all the point coordinates are
       !available on SP PEs
       if(is_proc0(SP_))then
          nU_I=ubound_vector('SP_Xyz_DI')
          nPoint=nU_I(2)
          call get_root_decomposition(&
               SP_,&
               iRootMapDim_D=(/1/),&
               XyzMin_D=(/cHalf/),&
               XyzMax_D=(/cHalf+real(nPoint)/),&
               nCells_D=(/nPoint/))
       end if
       call bcast_decomposition(SP_)
       call set_standard_grid_descriptor(SP_,GridDescriptor=&
            SP_GridDescriptor)
       nU_I(2:2)=ncells_decomposition_d(SP_)
       nPoint=nU_I(2)
       if(.not.is_proc0(SP_))call allocate_vector('SP_Xyz_DI',3,nPoint)
       if(use_comp(SC_))then  !^CMP IF SC BEGIN
          call set_standard_grid_descriptor(SC_,GridDescriptor=&
               SC_GridDescriptor)
          call init_router(SC_GridDescriptor,SP_GridDescriptor,&
               RouterScSp)    !^CMP END SC
       end if
       if(use_comp(IH_))then  !^CMP IF IH BEGIN
          call set_standard_grid_descriptor(IH_,GridDescriptor=&
               IH_GridDescriptor)
          call init_router(IH_GridDescriptor,SP_GridDescriptor,&
               RouterIhSp)
       end if                 !^CMP END IH
    end if

    call bcast_global_vector('SP_Xyz_DI',i_proc0(SP_),i_comm())
       
    
    !\
    !Set masks
    !/
    if(use_comp(IH_))then !^CMP IF IH BEGIN
       if(RouterIhSp%IsProc)then
          call allocate_mask('SP_IsInIH','SP_Xyz_DI')
          call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
          if(is_proc0().and.DoTest)&
               write(*,*)'Mask SP_IsInIH count=',count_mask('SP_IsInIH')
       end if
    end if                !^CMP END IH
    if(use_comp(SC_))then !^CMP IF SC BEGIN
       if(RouterScSp%IsProc)then
          call allocate_mask('SP_IsInSC','SP_Xyz_DI')
          call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
          if(is_proc0().and.DoTest)&
               write(*,*)'Mask SP_IsInSC count=',count_mask('SP_IsInSC')
       end if
    end if                !^CMP END SC
    if(is_proc(SP_))then
       call allocate_vector('SP_XyzSP',3,nPoint)
       if(use_comp(SC_))call transform_to_sp_from(SC_)!^CMP IF SC
       if(use_comp(IH_))call transform_to_sp_from(IH_)!^CMP IF IH
       if(DoTest.and.is_proc0(SP_))&
            call save_global_vector('SP_XyzSP',iFileIn=0)
    end if

  contains
    !==================================================================
    subroutine trace_line(iComp,Gd,Router,R2,MH_get_a_line_point)
      use CON_global_message_pass
      interface
         subroutine MH_get_a_line_point(nPartial,&
              iGetStart,&
              Get,&
              Weight,&
              Buff_I,nVar)
           use CON_router
           implicit none
           integer,intent(in)::nPartial,iGetStart,nVar
           type(IndexPtrType),intent(in)::Get
           type(WeightPtrType),intent(in)::Weight
           real,dimension(nVar),intent(out)::Buff_I
         end subroutine MH_get_a_line_point
      end interface
      
      integer,intent(in)::iComp
      type(GridDescriptorType),intent(out)::Gd
      type(RouterType),intent(out)::Router
      real,intent(in)::R2
      
      integer::iPointStart
      real,save::SignBDotR
      !--------------------------------------------------------!
      call set_standard_grid_descriptor(iComp,GridDescriptor=Gd)
      call init_router(Gd,SP_GridDescriptor,Router)
      if(.not.Router%IsProc)return
      call MPI_bcast(iPoint,1,MPI_INTEGER,&
           Router%iProc0Target,Router%iCommUnion,iError)
      if(iPoint>1)then
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint-1),6,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
      else
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint),3,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
      end if
      iPointStart=iPoint
      if(DoTest.and.is_proc0(SP_))write(*,*)'Start line from xyz=',&
           XyzTemp_DI(:,iPointStart),' at PE=',i_proc()
      do while (dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint))>=R2)
         if(iPoint>1)then
            if(dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint))>=&
                 dot_product(XyzTemp_DI(:,iPoint-1),XyzTemp_DI(:,iPoint-1)))then
               if(is_proc0(SP_))nPoint=iPoint
               EXIT
            end if
         end if
!         if(DoTest.and.is_proc0(SP_))write(*,*)NameSub//': done ',iPoint,&
!              ' Line Point, Xyz=',XyzTemp_DI(:,iPoint),' at PE=',i_proc()
         if(iPoint==nPointMax)call CON_stop(&
              'Insufficient number of points in couple_ih_sp =',iPoint)
         
         call set_router(& 
              GridDescriptorSource=Gd,&
              GridDescriptorTarget=SP_GridDescriptor,&
              Router=Router,&
              mapping=xyz_point_mapping,&
              interpolate=interpolation_fix_reschange)     
         call global_message_pass(Router,&
              nVar=6,&
              fill_buffer=MH_get_a_line_point,&
              apply_buffer=SP_put_a_line_point)
         if(is_proc0(SP_))then
            if(iPoint==iPointStart)then
               SignBDotR=cOne
               XyzTemp_DI(:,iPoint+1)=XyzTemp_DI(:,iPoint)-DsResolution*&
                    cHalf*bDxyz_I(1:3)/sqrt(dot_product(&
                    bDxyz_I(1:3),bDxyz_I(1:3)))*&
                    sqrt(dot_product(&
                    bDxyz_I(4:6),bDxyz_I(4:6)))          
               if(dot_product(XyzTemp_DI(:,iPoint+1),XyzTemp_DI(:,iPoint+1))>&
                    dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint)))&
                    SignBDotR=-cOne
            end if
            XyzTemp_DI(:,iPoint+1)=XyzTemp_DI(:,iPoint)-DsResolution*&
                 cHalf*bDxyz_I(1:3)/sqrt(dot_product(&
                 bDxyz_I(1:3),bDxyz_I(1:3)))*&
                 sqrt(dot_product(&
                 bDxyz_I(4:6),bDxyz_I(4:6)))*&
                 SignBDotR
            nPoint=iPoint
         end if
         iPoint=iPoint+1
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint),3,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
         call set_router(& 
              GridDescriptorSource=Gd,&
              GridDescriptorTarget=SP_GridDescriptor,&
              Router=Router,&
              mapping=xyz_point_mapping,&
              interpolate=interpolation_fix_reschange)     
         call global_message_pass(Router,&
              nVar=6,&
              fill_buffer=MH_get_a_line_point,&
              apply_buffer=SP_put_a_line_point)
         if(is_proc0(SP_))then     
            XyzTemp_DI(:,iPoint)=XyzTemp_DI(:,iPoint-1)-DsResolution*&
                 bDxyz_I(1:3)/sqrt(dot_product(&
                 bDxyz_I(1:3),bDxyz_I(1:3)))*&
                 sqrt(dot_product(&
                 bDxyz_I(4:6),bDxyz_I(4:6)))*&
                 SignBDotR
         end if
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint),3,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
      end do
    end subroutine trace_line
  end subroutine couple_mh_sp_init
  !==================================================================!
  subroutine transform_to_sp_from(iComp)
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::SP_LocalXyz_DI
    logical,pointer,dimension(:)::Is_I
    integer::nU_I(2),i
    real,dimension(3,3)::MhToSp_DD
    character(LEN=2)::NameComp
    real::LengthRatio
    call get_comp_info(iComp,Name=NameComp)
    MhToSp_DD=transform_matrix(tNow,&
            Grid_C(iComp)%TypeCoord,&
            Grid_C(SP_)%TypeCoord)
    if(DoTest)write(*,*)'Transform SP coordinates from '//NameComp
    call associate_with_global_mask(Is_I,'SP_IsIn'//NameComp)
    call associate_with_global_vector(SP_LocalXyz_DI,'SP_XyzSP')
    nU_I=ubound(SP_LocalXyz_DI)
    if(DoTest)write(*,*)nU_I
    LengthRatio=Grid_C(iComp)%UnitX/Grid_C(SP_)%UnitX
    do i=1,nU_I(2)
       if(.not.Is_I(i))CYCLE
       SP_LocalXyz_DI(:,i)=matmul(MhToSp_DD,&
            point_state_v('SP_Xyz_DI',3,i))*LengthRatio
    end do
    nullify(SP_LocalXyz_DI,Is_I)
  end subroutine transform_to_sp_from
    
  !==================================================================
  !^CMP IF IH BEGIN
  subroutine couple_ih_sp(DataInputTime)     
    use CON_global_message_pass
    interface
       subroutine SP_put_from_mh(nPartial,&
            iPutStart,&
            Put,&
            Weight,&
            DoAdd,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iPutStart,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::Weight
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::Buff_I
       end subroutine SP_put_from_mh
    end interface

    real,intent(in)::DataInputTime
    real,dimension(3)::Xyz_D
    !-------------------------------------------------------------------------

    if(.not.RouterIhSp%IsProc)return

    tNow=DataInputTime
    IhToSp_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    ScToIh_DD=transform_matrix(tNow,&                   !^CMP IF SC
         Grid_C(SC_)%TypeCoord, Grid_C(IH_)%TypeCoord)  !^CMP IF SC


    call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion)
    call bcast_global_vector('SP_Xyz_DI',&
         RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion,&
         'SP_IsInIH')
    call set_router(& 
         GridDescriptorSource=IH_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterIhSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInIH',&
         interpolate=interpolation_fix_reschange)

    if(is_proc(SP_))then
       call SP_put_input_time(DataInputTime)
       call transform_to_sp_from(IH_)
    end if

    call global_message_pass(RouterIhSp,&
         nVar=8,&
         fill_buffer=IH_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    !^CMP IF SC BEGIN
    !This coupler is performed after SC-SP coupling, so that 
    !on SP the updated coordinates are available for those
    !points which passed from SC to IH

    if(use_comp(SC_))then              
       if(is_proc0(SP_))then
          !Check the points which passed from SC to IH:
          call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
          call associate_with_global_mask(Is_I,'SP_IsInIH')

          do iPoint=1,nPoint
             if(Is_I(iPoint))CYCLE
             !Consider a point which before was   not in IH
             Xyz_D=matmul(ScToIh_DD,Xyz_DI(:,iPoint))*&
                  (Grid_C(SC_)%UnitX/Grid_C(IH_)%UnitX) 
             if(is_in_ih(Xyz_D))& ! Now the point is in IH
                  Xyz_DI(:,iPoint)=Xyz_D
             !..that is why we convert it to IH coordinates
          end do
          nullify(Xyz_DI)
          nullify(Is_I)
       end if
       call bcast_global_vector('SP_Xyz_DI',&
            RouterIhSp%iProc0Target,&
            RouterIhSp%iCommUnion)           !^CMP END SC
    end if
    call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
  end subroutine couple_ih_sp
  !==================================================================
  logical function is_in_ih(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    is_in_ih=dot_product(Xyz_D,Xyz_D)>=rBoundIh**2.and.&
         all(Xyz_D<=xyz_max_d(IH_)).and.all(Xyz_D>=xyz_min_d(IH_))
  end function is_in_ih
  !==================================================================!        
   subroutine IH_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
       BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call IH_get_for_sp(&
         nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)
    
    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State3_V(BuffX_:BuffZ_),&
         Grid_C(IH_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(IhToSp_DD,State_V(Bx_:Bz_))
  end subroutine IH_get_for_sp_and_transform
  !^CMP END IH
  !=========================================================================
  !^CMP IF SC BEGIN
  subroutine couple_sc_sp(DataInputTime)
    use CON_global_message_pass
    interface
       subroutine SP_put_from_mh(nPartial,&
            iPutStart,&
            Put,&
            Weight,&
            DoAdd,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iPutStart,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::Weight
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::Buff_I
       end subroutine SP_put_from_mh
    end interface
    real,intent(in)::DataInputTime
 
    if(.not.RouterScSp%IsProc)return

    tNow=DataInputTime
    ScToSP_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)

    call SC_synchronize_refinement(RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion)
    call bcast_global_vector('SP_Xyz_DI',&
         RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion,&
         'SP_IsInSC')
    call set_router(& 
         GridDescriptorSource=SC_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterScSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInSC',&
         interpolate=interpolation_fix_reschange)
    if(is_proc(SP_))then
       call SP_put_input_time(DataInputTime)  
       call transform_to_sp_from(SC_)
    end if
    call global_message_pass(RouterScSp,&
         nVar=8,&
         fill_buffer=SC_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
  end subroutine couple_sc_sp
  !-------------------------------------------------------------------------
  logical function is_in_sc(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    real::R2
    R2=dot_product(Xyz_D,Xyz_D)
    if(use_comp(IH_))then            !^CMP IF IH BEGIN
       is_in_sc=R2>=rBoundSc**2.and.R2<rBoundIh**2
    else                             !^CMP END IH
       is_in_sc=R2>=rBoundSc**2.and.&
            all(Xyz_D<=xyz_max_d(SC_)).and.all(Xyz_D>=xyz_min_d(SC_))
    end if                           !^CMP IF IH
  end function is_in_sc             
  !--------------------------------------------------------------------------
   subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)
    
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
       BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call SC_get_for_sp(&
         nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)
    
    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State3_V(BuffX_:BuffZ_),&
         Grid_C(SC_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(ScToSp_DD,State_V(Bx_:Bz_))
  end subroutine SC_get_for_sp_and_transform
  !=========================================================================
  !^CMP END SC
end Module CON_couple_mh_sp
