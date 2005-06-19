!^CFG COPYRIGHT UM
!^CMP FILE SP
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
  real::BAndDXyz_I(1:6)!The interpolated values of full B and DXyz
  real::DsResolution,XyzLine_D(3)
  real,save::RBoundIH                !^CMP IF IH
  real,save::RBoundSC                !^CMP IF SC
  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real,dimension(3,3)::ScToIh_DD

  real :: tNow

!====================================================================!
!The coupler is generalized for the following assumptions:
!The IH coordinate system is assumed to be inertial (non-rotatating).
!The SP coordinate system is assumed to be the same as in IH.
!The SC coordinate system can be arbitrary, if IH is used, otherwise 
!it is inertial and coincides with SP
!The unit of length in IH,SC is rSun
contains
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
    !The initialization can be done only twice
    
    !Initialize grid
    call init_decomposition(SP_,SP_,1)

    !Check, if the lagrangian meshes are saved in restart files


    if(is_proc0(SP_))then
       DoneRestart=.false.
       if(used_mask('SP_IsInSC'))DoneRestart=.true. !^CMP IF SC
       if(used_mask('SP_IsInIH'))DoneRestart=.true. !^CMP IF IH
    end if
    call MPI_BCAST(DoneRestart,1,MPI_LOGICAL,i_proc0(SP_),i_comm(),iError)
    if(.not.DoneRestart)then
       !Construct auxiliary SP_grid with ONE point per line
       if(is_proc0(SP_))&
            call get_root_decomposition(&
            SP_,&                             !GridDescriptor to be constructed
            iRootMapDim_D=(/1/),&             !The block amount, along each direction(D)
            XyzMin_D=(/cHalf/),&              !Minimal gen. coordinates, along each D 
            XyzMax_D=(/cHalf+cOne/),&        !Maximal gen. coordinates, along each D
            nCells_D=(/1/),&
            PE_I=(/0/))                     !PEs layout, throughout all the blocks
       call bcast_decomposition(SP_)
       

       call set_standard_grid_descriptor(SP_,GridDescriptor=&
            SP_GridDescriptor)

       call get_time(tSimulationOut=tNow)
       !Set the initial point for the line

       if(is_proc0(SP_))then
          call SP_get_line_param(&
            DsResolution,XyzLine_D&
            ,RBoundSC &              !^CMP IF SC
            ,RBoundIH &              !^CMP IF IH
            )
          if(sum(XyzLine_D**2)<cOne)then
             !Calculate the Earth position in SGI
             XyzLine_D = XyzPlanetHgi_D/rSun   !In HGI
             if(use_comp(IH_))then             !^CMP IF IH BEGIN
                XyzLine_D=matmul(&
                     transform_matrix(tNow,'HGI',&
                     Grid_C(IH_) % TypeCoord),XyzLine_D)
             else                              !^CMP END IH
                XyzLine_D=matmul(&
                     transform_matrix(tNow,'HGI',&
                     Grid_C(SC_) % TypeCoord),XyzLine_D)
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
  
          call MPI_bcast(RBoundIH,1,MPI_REAL,&
               i_proc0(SP_),i_comm(),iError)
          call trace_line(IH_,&
               IH_GridDescriptor,RouterIhSp,RBoundIH**2,IH_get_a_line_point)
          call MPI_bcast(nPoint,1,MPI_INTEGER,&
               i_proc0(SP_),i_comm(),iError)

          if(DoTest.and.is_proc0())write(*,*)'nPointIH=',nPoint
       end if                      !^CMP END IH

       if(use_comp(SC_))then       !^CMP IF SC BEGIN
          call SC_synchronize_refinement(i_proc0(SC_),i_comm())
          call MPI_bcast(RBoundSC,1,MPI_REAL,&
               i_proc0(SP_),i_comm(),iError)
          if(Grid_C(SC_)%TypeCoord/=Grid_C(IH_)%TypeCoord& !^CMP IF IH BEGIN
               .and.is_proc0(SP_))&
               XyzTemp_DI(:,nPoint)=matmul(transform_matrix(&
               tNow,Grid_C(IH_)%TypeCoord,Grid_C(SC_)%TypeCoord),&
               XyzTemp_DI(:,nPoint))                       !^CMP END IH


          call trace_line(SC_,&
               SC_GridDescriptor,RouterScSp,RBoundSC**2,SC_get_a_line_point)
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
    end if

    call bcast_global_vector('SP_Xyz_DI',i_proc0(SP_),i_comm())
       
    
    !\
    !Set masks
    !/
    if(use_comp(IH_))then !^CMP IF IH BEGIN
       if(RouterIhSp%IsProc)then
          if(.not.DoneRestart.or..not.is_proc0(SP_))&
               call allocate_mask('SP_IsInIH','SP_Xyz_DI')
          call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
          if(is_proc0().and.DoTest)&
               write(*,*)'Mask SP_IsInIH count=',count_mask('SP_IsInIH')
       end if
    end if                !^CMP END IH
    if(use_comp(SC_))then !^CMP IF SC BEGIN
       if(RouterScSp%IsProc)then
          if(.not.DoneRestart.or..not.is_proc0(SP_))&
               call allocate_mask('SP_IsInSC','SP_Xyz_DI')
          call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
          if(is_proc0().and.DoTest)&
               write(*,*)'Mask SP_IsInSC count=',count_mask('SP_IsInSC')
       end if
    end if                !^CMP END SC
  contains
    !==================================================================
    subroutine trace_line(CompID_,GD,Router,R2,MH_get_a_line_point)
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
      
      integer,intent(in)::CompID_
      type(GridDescriptorType),intent(out)::GD
      type(RouterType),intent(out)::Router
      real,intent(in)::R2
      
      integer::iPointStart
      real,save::SignBDotR
      !--------------------------------------------------------!
      call set_standard_grid_descriptor(CompID_,GridDescriptor=GD)
      call init_router(GD,SP_GridDescriptor,Router)
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
         if(DoTest.and.is_proc0(SP_))write(*,*)NameSub//': done ',iPoint,&
              ' Line Point, Xyz=',XyzTemp_DI(:,iPoint),' at PE=',i_proc()
         if(iPoint==nPointMax)call CON_stop(&
              'Insufficient number of points in couple_ih_sp =',iPoint)
         
         call set_router(& 
              GridDescriptorSource=GD,&
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
                    cHalf*BAndDXyz_I(1:3)/sqrt(dot_product(&
                    BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
                    sqrt(dot_product(&
                    BAndDXyz_I(4:6),BAndDXyz_I(4:6)))          
               if(dot_product(XyzTemp_DI(:,iPoint+1),XyzTemp_DI(:,iPoint+1))>&
                    dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint)))&
                    SignBDotR=-cOne
            end if
            XyzTemp_DI(:,iPoint+1)=XyzTemp_DI(:,iPoint)-DsResolution*&
                 cHalf*BAndDXyz_I(1:3)/sqrt(dot_product(&
                 BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
                 sqrt(dot_product(&
                 BAndDXyz_I(4:6),BAndDXyz_I(4:6)))*&
                 SignBDotR
            nPoint=iPoint
         end if
         iPoint=iPoint+1
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint),3,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
         call set_router(& 
              GridDescriptorSource=GD,&
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
                 BAndDXyz_I(1:3)/sqrt(dot_product(&
                 BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
                 sqrt(dot_product(&
                 BAndDXyz_I(4:6),BAndDXyz_I(4:6)))*&
                 SignBDotR
         end if
         call MPI_bcast(&
              XyzTemp_DI(1,iPoint),3,MPI_REAL,&
              Router%iProc0Target,Router%iCommUnion,iError)
      end do
    end subroutine trace_line
  end subroutine couple_mh_sp_init
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
       W,&
       DoAdd,&
       Buff_I,nVar)
    implicit none
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer::iCell
    real:: Weight
    !-------------------------------------------------------------------------
    Weight=W%Weight_I(iPutStart)
    if(DoAdd)then
       BAndDXyz_I=BAndDXyz_I+Buff_I(:)*Weight
    else
       BAndDXyz_I=Buff_I(:)*Weight
    end if
  end subroutine SP_put_a_line_point
  !==================================================================
  subroutine couple_ih_sp(DataInputTime)     !^CMP IF IH BEGIN
    use CON_global_message_pass
    interface
       subroutine IH_get_for_sp(nPartial,&
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
       end subroutine IH_get_for_sp
    end interface
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
    !-------------------------------------------------------------------------

    if(.not.RouterIhSp%IsProc)return
    tNow=DataInputTime
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
    if(is_proc(SP_))call SP_put_input_time(DataInputTime)     
    call global_message_pass(RouterIhSp,&
         nVar=8,&
         fill_buffer=IH_get_for_sp,&
         apply_buffer=SP_put_from_mh)
    !^CMP IF SC BEGIN
    !This coupler is performed after SC-SP coupling, so that 
    !on SP the updated coordinates are available for those
    !points which passed from SC to IH

    if(use_comp(SC_))then
       !Reset the transformation matrix 
       !Check the points which passed from SC to IH:
       ScToIh_DD=transform_matrix(DataInputTime,&
            Grid_C(SC_)%TypeCoord, Grid_C(IH_)%TypeCoord)

       call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
       call associate_with_global_mask(Is_I,'SP_IsInIH')

       do iPoint=1,nPoint
          if(.not.Is_I(iPoint) &      !This point before was   not in IH
               .and.is_in_ih(Xyz_DI(:,iPoint)))& !... and now it is in IH
               Xyz_DI(:,iPoint)=matmul(ScToIh_DD,Xyz_DI(:,iPoint)) 
                !..that is why we convert it to IH coordinates
       end do
       nullify(Xyz_DI)
       nullify(Is_I)
       call bcast_global_vector('SP_Xyz_DI',&
            RouterIhSp%iProc0Target,&
            RouterIhSp%iCommUnion)           !^CMP END SC
    end if
    call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
  end subroutine couple_ih_sp
  !==================================================================
  logical function is_in_ih(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    is_in_ih=dot_product(Xyz_D,Xyz_D)>=RBoundIH**2.and.&
         all(Xyz_D<=xyz_max_d(IH_)).and.all(Xyz_D>=xyz_min_d(IH_))
  end function is_in_ih                      !^CMP END IH
  !=========================================================================
  subroutine couple_sc_sp(DataInputTime)     !^CMP IF SC BEGIN
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
    ScToIh_DD=transform_matrix(DataInputTime,&
         Grid_C(SC_)%TypeCoord, Grid_C(IH_)%TypeCoord)

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
    if(is_proc(SP_))call SP_put_input_time(DataInputTime)     
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
       is_in_sc=R2>=RBoundSC**2.and.R2<RBoundIH**2
    else                             !^CMP END IH
       is_in_sc=R2>=RBoundSC**2.and.&
            all(Xyz_D<=xyz_max_d(SC_)).and.all(Xyz_D>=xyz_min_d(SC_))
    end if                           !^CMP IF IH
  end function is_in_sc             
  !--------------------------------------------------------------------------
   subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,W,State_V,nVar)
    
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
       BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call SC_get_for_sp(&
         nPartial,iGetStart,Get,W,State3_V,nVar+3)
    State_V=State3_V(1:nVar)
    
    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State3_V(BuffX_:BuffZ_),&
         Grid_C(SC_)%TypeCoord,Grid_C(IH_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(ScToIh_DD,State_V(Bx_:Bz_))
  end subroutine SC_get_for_sp_and_transform
  !=========================================================================
end Module CON_couple_mh_sp
